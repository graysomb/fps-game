#!/usr/bin/env python3
"""Paired end-to-end adaptive gate. Never enables adaptive mode or edits defaults.

Includes cold construction and topology events in each trial's amortized mean.
Warm-up is a separate complete trial, not omitted slow frames from measured runs.
Only the harness's timed physics/render interval is called frame time; validation
and PNG capture are deliberately outside it. Raw reports remain alongside results.
"""
import argparse
import json
import hashlib
import shutil
import math
import os
from pathlib import Path
import random
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ["adaptive-solid", "adaptive-hollow", "adaptive-thin", "adaptive-irregular", "adaptive-coarse-irregular",
             "adaptive-fracture", "dynamic-freefall", "overhang-impact",
             "pillar-impact-drift", "tether-thin-wall-ccd", "sleep-wake-floating",
             "voxel-strain-shear-fracture"]

def percentile(values, q):
    values = sorted(values)
    if not values:
        return None
    i = (len(values) - 1) * q
    lo = int(i)
    return values[lo] + (values[min(lo + 1, len(values)-1)]-values[lo])*(i-lo)

def paired_ci(differences, seed=731, samples=10000):
    if len(differences) < 10 or not all(math.isfinite(x) for x in differences):
        return None
    rng = random.Random(seed)
    n = len(differences)
    draws = [statistics.mean(differences[rng.randrange(n)] for _ in range(n))
             for _ in range(samples)]
    return [percentile(draws, .025), percentile(draws, .975)]

def gate(pairs):
    if len(pairs) < 10 or any(not p[m]["valid"] for p in pairs for m in ("off", "on")):
        return {"eligible": False, "status": "INVALID_OR_INSUFFICIENT_TRIALS"}
    if not all(p["on"].get("coarsened", False) for p in pairs):
        return {"eligible": False, "status": "INELIGIBLE_NO_COARSENING"}
    metrics = ("physics_ms", "frame_ms", "frame_p95_ms", "event_p95_ms")
    intervals = {}
    for metric in metrics:
        deltas = [p["on"][metric]-p["off"][metric] for p in pairs
                  if p["on"][metric] is not None and p["off"][metric] is not None]
        intervals[metric] = paired_ci(deltas)
    improvement = intervals["physics_ms"] is not None and intervals["physics_ms"][1] < 0
    regressions = [m for m in metrics[1:] if intervals[m] is not None and intervals[m][0] > 0]
    same_syncs = all(p["on"]["syncs"] <= p["off"]["syncs"] for p in pairs)
    eligible = improvement and not regressions and same_syncs
    return {"eligible": eligible,
            "status": "PASS" if eligible else "REGRESSION" if regressions or not same_syncs else "INCONCLUSIVE_OR_SLOWER",
            "paired_mean_difference_ci95_ms": intervals,
            "regressions": regressions, "no_added_syncs": same_syncs}

def run(binary, scenario, mode, output, steps, runtime=ROOT):
    output.mkdir(parents=True, exist_ok=True)
    command = [str(binary), "--physics=gpu-metal", f"--adaptive-physics={mode}",
               f"--debug-scenario={scenario}", f"--debug-steps={steps}",
               "--debug-capture-steps=0", f"--debug-output={output}"]
    env = dict(os.environ, FPS_ADAPTIVE_BENCHMARK="1")
    env.pop("FPS_ADAPTIVE_PROFILE", None)
    with (output / "run.log").open("w") as log:
        result = subprocess.run(command, cwd=runtime, env=env, stdout=log,
                                stderr=subprocess.STDOUT, timeout=300)
    report = output / "report.json"
    timings = output / "report.json.timings.jsonl"
    rows = [json.loads(line) for line in timings.read_text().splitlines()] if timings.exists() else []
    data = json.loads(report.read_text()) if report.exists() else {}
    valid = result.returncode == 0 and data.get("passed") is True and not data.get("fallback")
    valid = valid and data.get("activeBackend") == "gpu-metal" and len(rows) == steps
    valid = valid and all(not row["error"] and all(math.isfinite(row[k]) and row[k] >= 0
                        for k in ("physicsMs", "frameMs", "gpuMs", "waitMs")) for row in rows)
    summary = {"valid": bool(valid), "report": str(report), "rows": rows,
               "coarsened": any(r["scratchBytes"] > 0 and r["leaves"] < r["fineCells"] for r in rows)}
    if rows:
        summary.update(physics_ms=statistics.mean(r["physicsMs"] for r in rows),
                       frame_ms=statistics.mean(r["frameMs"] for r in rows),
                       frame_p95_ms=percentile([r["frameMs"] for r in rows], .95),
                       syncs=sum(r["syncs"] for r in rows),
                       scratch_bytes=max(r["scratchBytes"] for r in rows))
    else:
        summary.update(physics_ms=None, frame_ms=None, frame_p95_ms=None, syncs=0, scratch_bytes=0)
    summary["event_p95_ms"] = None
    return summary

def pair_events(pair):
    # Use the union of event steps so both modes measure the same time windows.
    events = set()
    for mode in ("off", "on"):
        last_generation = 0
        for row in pair[mode]["rows"]:
            if row["repacked"] or row["generation"] != last_generation:
                events.add(row["step"])
            last_generation = row["generation"]
    for mode in ("off", "on"):
        pair[mode]["event_p95_ms"] = percentile(
            [r["frameMs"] for r in pair[mode]["rows"] if r["step"] in events], .95)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / ".build/adaptive/fps_ray_gpu")
    parser.add_argument("--output", type=Path, default=ROOT / ".build/adaptive/benchmarks")
    parser.add_argument("--scenarios", nargs="+", default=SCENARIOS)
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--steps", type=int, default=120)
    args = parser.parse_args()
    if args.runs < 10 or args.steps < 60:
        parser.error("the release gate requires at least 10 pairs and 60 steps")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    runtime = args.output / "runtime"
    if runtime.exists():
        parser.error("output already contains a benchmark snapshot; choose a new output directory")
    runtime.mkdir()
    binary = runtime / "fps_ray_gpu"
    shutil.copy2(args.binary.resolve(), binary)
    shutil.copytree(ROOT / "shaders", runtime / "shaders")
    shader_hash = hashlib.sha256()
    for path in sorted((runtime / "shaders").rglob("*")):
        if path.is_file():
            shader_hash.update(str(path.relative_to(runtime)).encode())
            shader_hash.update(path.read_bytes())
    all_results = {"schema": 1, "seed": 731, "runs": args.runs, "steps": args.steps,
                   "binarySha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                   "shaderSha256": shader_hash.hexdigest(),
                   "cubeSize": os.environ.get("FPS_CUBE_SIZE", "8"), "scenarios": {}}
    rng = random.Random(731)
    for scenario in args.scenarios:
        for mode in ("off", "on"):
            run(binary, scenario, mode, args.output/scenario/f"warmup-{mode}", args.steps, runtime)
        pairs = []
        for trial in range(args.runs):
            modes = ["off", "on"]
            rng.shuffle(modes)
            pair = {mode: run(binary, scenario, mode,
                             args.output/scenario/f"pair-{trial:02}-{mode}", args.steps, runtime) for mode in modes}
            pair_events(pair)
            pairs.append(pair)
            print(f"{scenario} pair {trial+1}/{args.runs}: fine={pair['off']['physics_ms']} adaptive={pair['on']['physics_ms']}", flush=True)
        all_results["scenarios"][scenario] = {"gate": gate(pairs), "pairs": pairs}
        (args.output/"results.json").write_text(json.dumps(all_results, indent=2, allow_nan=False)+"\n")
        print(scenario, all_results["scenarios"][scenario]["gate"], flush=True)
    return 0 if all(x["gate"]["eligible"] for x in all_results["scenarios"].values()) else 1

if __name__ == "__main__":
    raise SystemExit(main())
