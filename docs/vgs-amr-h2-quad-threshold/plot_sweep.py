#!/usr/bin/env python3
"""Plot rendered terminal cells for 4x-per-level curvature thresholds."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
by_threshold = {}
for name in ("pilot.csv", "detail.csv"):
    with (OUT / name).open(newline="") as file:
        for row in csv.DictReader(file):
            parsed = {"root_threshold": float(row["root_threshold"])}
            parsed.update({k: int(v) for k, v in row.items() if k != "root_threshold"})
            threshold = parsed["root_threshold"]
            if threshold in by_threshold and by_threshold[threshold] != parsed:
                raise ValueError(f"Conflicting results for root threshold {threshold}")
            by_threshold[threshold] = parsed
rows = [by_threshold[key] for key in sorted(by_threshold)]
with (OUT / "threshold-sweep-all.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

fig = plt.figure(figsize=(12.8, 8.0), dpi=160)
grid = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.18, top=0.85,
                        hspace=0.47, wspace=0.22)
axes = (fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]))
fig.suptitle("Adaptive cube drop · 4× curvature threshold per finer level",
             fontsize=18, fontweight="bold", y=0.96)

series = (
    ("leaves120", "Frame 120", "#a85c50", "o"),
    ("leaves140", "Frame 140", "#625b8c", "s"),
    ("leaves240", "Frame 240", "#16796d", "^"),
)
subsets = (rows,
           [row for row in rows if row["root_threshold"] <= 1],
           [row for row in rows if 1 <= row["root_threshold"] <= 10])
titles = ("Full root-threshold sweep", "Detail: root thresholds 0.05–1",
          "Detail: root thresholds 1–10")
for index, (ax, subset) in enumerate(zip(axes, subsets)):
    x = [row["root_threshold"] for row in subset]
    for key, label, color, marker in series:
        ax.plot(x, [row[key] for row in subset], label=label, color=color,
                linewidth=2.3, marker=marker, markersize=5)
    ax.set_title(titles[index], fontweight="bold")
    ax.set_yscale("log")
    ax.set_ylim(0.8, 50000)
    ax.set_xlabel("Root curvature threshold (s⁻²)")
    ax.set_ylabel("Rendered terminal cells")
    ax.grid(True, which="major", color="#cfcfcf", linewidth=0.6, alpha=0.8)
    ax.grid(True, which="minor", color="#d8d8d8", linewidth=0.4, alpha=0.5)
    if index == 0:
        ax.set_xscale("log")
        ax.set_xticks((0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                      labels=("0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50"))
        ax.set_xlim(0.045, 55)
        ax.legend(loc="upper right", frameon=False)
    elif index == 1:
        ax.set_xlim(0.03, 1.02)
    else:
        ax.set_xlim(0.95, 10.25)

fig.text(0.5, 0.105,
         "T(level) = 4^level × T(root) · h²-scaled VGS material · breaking disabled",
         ha="center", fontsize=10, color="#555b68")
fig.text(0.5, 0.072,
         f"0.15625 m effective finest cells · 240 frames · flat drop · Metal · {len(rows)} runs passed without fallback",
         ha="center", fontsize=10, color="#555b68")
fig.savefig(OUT / "active-resolution.png", dpi=160)
fig.savefig(OUT / "active-resolution.svg")
