#!/usr/bin/env python3
"""Run debug harness scenarios and measure physics execution timing."""

import argparse
import json
import os
import pathlib
import statistics
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]

DEFAULT_SCENARIOS = [
    "dynamic-freefall",
    "falling-pillar-sleep",
    "overhang-impact",
    "pillar-impact-drift",
    "tether-thin-wall-ccd",
    "sleep-wake-floating",
    "sleep-recycle",
    "creative-reset-cache",
]

def run_scenario(binary_path: str, scenario: str, steps: int = 60) -> dict:
    cmd = [binary_path, f"--debug-scenario={scenario}", f"--debug-steps={steps}"]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    wall_ms = (time.perf_counter() - t0) * 1000.0

    passed = proc.returncode == 0
    report_path = None
    output_lines = proc.stdout.splitlines() + proc.stderr.splitlines()
    for line in output_lines:
        if "report=" in line:
            for part in line.split():
                if part.startswith("report="):
                    report_path = part.split("=", 1)[1]
                    break

    step_times = []
    if report_path and os.path.exists(report_path):
        try:
            with open(report_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                for s in data.get("samples", []):
                    ms = s.get("physicsStepMs", 0.0)
                    if ms > 0.0:
                        step_times.append(ms)
        except Exception as e:
            print(f"Warning: failed reading report {report_path}: {e}", file=sys.stderr)

    return {
        "scenario": scenario,
        "passed": passed,
        "returncode": proc.returncode,
        "wall_ms": wall_ms,
        "samples_count": len(step_times),
        "step_mean_ms": statistics.mean(step_times) if step_times else 0.0,
        "step_median_ms": statistics.median(step_times) if step_times else 0.0,
        "step_min_ms": min(step_times) if step_times else 0.0,
        "step_max_ms": max(step_times) if step_times else 0.0,
        "all_step_times": step_times,
    }


def benchmark(binary: str, scenarios: list[str], steps: int, runs: int) -> dict:
    results = {}
    for scn in scenarios:
        scenario_runs = []
        for r in range(runs):
            res = run_scenario(binary, scn, steps)
            scenario_runs.append(res)
        means = [r["step_mean_ms"] for r in scenario_runs if r["step_mean_ms"] > 0]
        medians = [r["step_median_ms"] for r in scenario_runs if r["step_median_ms"] > 0]
        walls = [r["wall_ms"] for r in scenario_runs]
        all_passed = all(r["passed"] for r in scenario_runs)
        results[scn] = {
            "passed": all_passed,
            "runs": runs,
            "mean_step_ms": statistics.mean(means) if means else 0.0,
            "median_step_ms": statistics.mean(medians) if medians else 0.0,
            "mean_wall_ms": statistics.mean(walls) if walls else 0.0,
            "runs_data": scenario_runs,
        }
    return results


def print_table(results: dict, baseline: dict = None):
    header = f"| {'Scenario':<24} | {'Status':<6} | {'Physics Step (ms)':<17} | {'Wall Time (ms)':<14} |"
    if baseline:
        header += f" {'Speedup':<10} |"
    sep = f"|{'-'*26}|{'-'*8}|{'-'*19}|{'-'*16}|"
    if baseline:
        sep += f"{'-'*12}|"
    print(header)
    print(sep)
    for scn, data in results.items():
        status = "PASS" if data["passed"] else "FAIL"
        step_ms = f"{data['mean_step_ms']:.3f} ms"
        wall_ms = f"{data['mean_wall_ms']:.1f} ms"
        row = f"| {scn:<24} | {status:<6} | {step_ms:<17} | {wall_ms:<14} |"
        if baseline and scn in baseline:
            base_ms = baseline[scn]["mean_step_ms"]
            cur_ms = data["mean_step_ms"]
            if base_ms > 0 and cur_ms > 0:
                speedup = (base_ms / cur_ms)
                pct = ((base_ms - cur_ms) / base_ms) * 100.0
                row += f" {speedup:.2f}x ({pct:+.1f}%) |"
            else:
                row += f" {'N/A':<10} |"
        elif baseline:
            row += f" {'N/A':<10} |"
        print(row)


def main():
    parser = argparse.ArgumentParser(description="FPS Game Physics Benchmark")
    parser.add_argument("--binary", default=str(ROOT / ".build/bin/fps_ray_gpu"), help="Path to binary")
    parser.add_argument("--scenarios", nargs="*", default=DEFAULT_SCENARIOS, help="Scenarios to run")
    parser.add_argument("--steps", type=int, default=60, help="Steps per scenario")
    parser.add_argument("--runs", type=int, default=3, help="Number of repetitions per scenario")
    parser.add_argument("--save-baseline", help="Save benchmark results to JSON file")
    parser.add_argument("--compare", help="Compare results against baseline JSON file")
    args = parser.parse_args()

    print(f"Running benchmark on {args.binary} ({args.runs} runs, {args.steps} steps)...")
    res = benchmark(args.binary, args.scenarios, args.steps, args.runs)

    baseline_data = None
    if args.compare and os.path.exists(args.compare):
        with open(args.compare, "r", encoding="utf-8") as f:
            baseline_data = json.load(f)

    print()
    print_table(res, baseline_data)
    print()

    if args.save_baseline:
        with open(args.save_baseline, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2)
        print(f"Saved baseline to {args.save_baseline}")


if __name__ == "__main__":
    main()
