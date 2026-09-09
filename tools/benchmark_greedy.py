#!/usr/bin/env python3
"""Run the identical irregular body through fine and greedy CPU ST/MT physics."""
import json
import pathlib
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN = ROOT / ".build/bin/fps_ray_sized_cpu"
OUT = ROOT / ".build/greedy-benchmark"
MODES = (("cpu-st", "fine"), ("cpu-st", "greedy"),
         ("cpu-mt", "fine"), ("cpu-mt", "greedy"))


def run(backend, mode, label):
    scenario = "greedy-activation-irregular" + ("-fine" if mode == "fine" else "")
    directory = OUT / f"{backend}-{mode}-{label}-{time.time_ns()}"
    subprocess.run([str(BIN), f"--physics={backend}", f"--debug-scenario={scenario}",
                    "--debug-steps=620", f"--debug-output={directory}"],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    return json.loads((directory / "report.json").read_text())["greedy"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = {}
    for backend, mode in MODES:
        run(backend, mode, "warmup")
        runs = [run(backend, mode, f"run-{i+1}") for i in range(3)]
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
