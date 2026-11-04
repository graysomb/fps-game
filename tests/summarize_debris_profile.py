"""Summarize completed baseline/optimized wave-seven logs; no timing data is invented."""
import csv
import json
import re
from pathlib import Path
result = {
    "baseline": "de12b07",
    "platform": "Windows, OpenGL 4.3 GPU, NVIDIA Quadro RTX 4000, 1280x720",
    "method": "Sequential runs, 60 FPS cap, actual elapsed dt, same seed and scripted workload; 100 seconds released before wave reset.",
    "limitations": [
        "Live physics and AI trajectories can diverge between runs.",
        "Frame work excludes intentional cap sleep. Phase means average approximately two-second reporting blocks.",
        "Metal is unavailable on this host; total GPU timing is not portable.",
    ],
    "runs": {},
    "comparisons": {},
}
for scenario in ("shield", "tower", "control"):
    for variant in ("baseline", "optimized"):
        name = f"debris-{scenario}-{variant}"
        path = Path("artifacts") / (name + ".log")
        if not path.exists():
            continue
        lines = path.read_text().splitlines()
        header = next(x for x in lines if x.startswith("frame,")).split(",")
        rows = [dict(zip(header, x.split(","))) for x in lines
                if x[:1].isdigit() and len(x.split(",")) == len(header)]
        summaries = {}
        for line in lines:
            m = re.match(r"SUMMARY (\S+) median_ms=([\d.]+) p95_ms=([\d.]+)", line)
            if m:
                summaries[m[1]] = {"median_ms": float(m[2]), "p95_ms": float(m[3])}
        if len(summaries) != 4:
            continue
        with path.with_suffix(".csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, header)
            writer.writeheader()
            writer.writerows(rows)
        release = 12 if scenario == "control" else 40
        windows = {}
        for label, low, high in (("tether", 10, release),
                                 ("early_released", release, release+30),
                                 ("settled", release+60, release+100)):
            selected = [r for r in rows if low < float(r["frame"])/60 <= high]
            means = {k: sum(float(r[k]) for r in selected)/len(selected)
                     for k in header if k.endswith("_ms")} if selected else {}
            means["activation_sleep_wake_ms"] = sum(means.get(k, 0) for k in
                                                    ("activation_ms", "wake_ms", "sleep_ms"))
            windows[label] = means
        result["runs"][name] = {"summary": summaries, "windows": windows}
    before = result["runs"].get(f"debris-{scenario}-baseline")
    after = result["runs"].get(f"debris-{scenario}-optimized")
    if before and after:
        result["comparisons"][scenario] = {
            window: 100 * (1-after["windows"][window]["activation_sleep_wake_ms"] /
                           before["windows"][window]["activation_sleep_wake_ms"])
            for window in ("tether", "early_released", "settled")
        }
Path("artifacts/debris-optimization-results.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result["comparisons"], indent=2))
