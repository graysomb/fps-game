#!/usr/bin/env python3
"""Benchmark identical fine and greedy bodies on selected physics backends."""
import argparse
import json
import pathlib
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_BIN = ROOT / ".build/bin/fps_ray_sized_cpu"
OUT = ROOT / ".build/greedy-benchmark"


def run(binary, backend, mode, label, batch):
    scenario = "greedy-activation-irregular" + ("-fine" if mode == "fine" else "")
    directory = OUT / f"{backend}-{mode}-{label}-{time.time_ns()}"
    subprocess.run([str(binary), f"--physics={backend}", f"--debug-scenario={scenario}",
                    "--debug-steps=620", f"--debug-physics-batch={batch}",
                    f"--debug-output={directory}"],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    report = json.loads((directory / "report.json").read_text())
    if report["activeBackend"] != backend:
        raise RuntimeError(f"{backend} executed as {report['activeBackend']}: "
                           f"{report.get('fallbackReason', '')}")
    if report["greedy"].get("fixedStepBatch") != batch:
        raise RuntimeError("debug harness did not report the requested fixed-step batch")
    return report["greedy"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=pathlib.Path, default=DEFAULT_BIN)
    parser.add_argument("--backends", nargs="+", default=["cpu-st", "cpu-mt"])
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--fixed-step-batch", type=int, default=1)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    if not 1 <= args.fixed_step_batch <= 8:
        parser.error("--fixed-step-batch must be between 1 and 8")
    OUT.mkdir(parents=True, exist_ok=True)
    results = {}
    for backend in args.backends:
        for mode in ("fine", "greedy"):
            run(args.binary, backend, mode, "warmup", args.fixed_step_batch)
            runs = [run(args.binary, backend, mode, f"run-{i+1}", args.fixed_step_batch) for i in range(args.runs)]
            assert all(r["physicsSteps"] == 600 for r in runs)
            totals = [r["physicsTotalMs"] for r in runs]
            results[f"{backend}-{mode}"] = {
                "runsMs": totals, "medianMs": statistics.median(totals),
                "rangeMs": [min(totals), max(totals)],
                "stagesMedianMs": {key: statistics.median(r[key] for r in runs) for key in
                    ("integrationMs", "vgsMs", "dynamicContactMs", "staticContactMs",
                     "interfaceRefreshMs", "passiveRefreshMs")},
                "particles": runs[0]["independent"], "groups": runs[0]["groups"]}
    path = OUT / "summary.json"
    path.write_text(json.dumps(results, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    main()
