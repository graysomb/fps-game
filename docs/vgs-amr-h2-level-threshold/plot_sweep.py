#!/usr/bin/env python3
"""Plot the 2^level curvature-threshold sweep and its low-threshold detail."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
by_threshold = {}
for name in ("threshold-sweep.csv", "low-threshold-sweep.csv", "mid-threshold-sweep.csv"):
    with (OUT / name).open(newline="") as file:
        for row in csv.DictReader(file):
            parsed = {k: int(v) for k, v in row.items()}
            threshold = parsed["root_threshold"]
            if threshold in by_threshold and by_threshold[threshold] != parsed:
                raise ValueError(f"Conflicting results for root threshold {threshold}")
            by_threshold[threshold] = parsed
rows = [by_threshold[key] for key in sorted(by_threshold)]
with (OUT / "threshold-sweep-all.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

fig, axes = plt.subplots(3, 2, figsize=(12.96, 11.5), dpi=160)
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.88,
                    wspace=0.25, hspace=0.55)
fig.suptitle("Adaptive cube drop · threshold doubles at each finer level",
             fontsize=18, fontweight="bold", y=0.96)

active_series = (
    ("leaves120", "Frame 120", "#a85c50", "o"),
    ("leaves140", "Frame 140", "#625b8c", "s"),
    ("leaves240", "Frame 240", "#16796d", "^"),
)
tree_series = (
    ("refinements", "Refinements", "#a85c50", "o"),
    ("coarsenings", "Coarsenings", "#625b8c", "s"),
)

subsets = (rows,
           [r for r in rows if r["root_threshold"] <= 10],
           [r for r in rows if 10 <= r["root_threshold"] <= 30])
for row_index, subset in enumerate(subsets):
    for column_index, series in enumerate((active_series, tree_series)):
        ax = axes[row_index, column_index]
        x = [r["root_threshold"] for r in subset]
        for key, label, color, marker in series:
            ax.plot(x, [r[key] for r in subset], label=label, color=color,
                    marker=marker, linewidth=2.2, markersize=5)
        ax.set_yscale("log")
        ax.grid(True, which="major", color="#d2d2d2", linewidth=0.6, alpha=0.75)
        ax.grid(True, which="minor", color="#d2d2d2", linewidth=0.4, alpha=0.45)
        ax.set_xlabel("Root curvature threshold (s⁻²)")
        ax.set_xlim(((0, 51), (0.7, 10.3), (9.5, 30.5))[row_index])
        if row_index == 1:
            ax.set_xticks(range(1, 11))
        elif row_index == 2:
            ax.set_xticks(range(10, 31, 2))
        if column_index == 0:
            ax.set_ylabel("Rendered terminal cells")
            ax.set_ylim(0.8, 50000)
        else:
            ax.set_ylabel("Transitions over 240 frames")
            ax.set_ylim(((40, 250000), (20000, 250000), (50, 100000))[row_index])
        if row_index == 0:
            ax.legend(loc="upper right", frameon=False, fontsize=9)

axes[0, 0].set_title("Active resolution · full sweep", fontweight="bold")
axes[0, 1].set_title("Tree activity · full sweep", fontweight="bold")
axes[1, 0].set_title("Active resolution · root thresholds 1–10", fontweight="bold")
axes[1, 1].set_title("Tree activity · root thresholds 1–10", fontweight="bold")
axes[2, 0].set_title("Active resolution · root thresholds 10–30", fontweight="bold")
axes[2, 1].set_title("Tree activity · root thresholds 10–30", fontweight="bold")

fig.text(0.5, 0.08,
         "T(level) = 2^level × T(root) · scaled VGS shear/stretch/volume · breaking disabled",
         ha="center", fontsize=10, color="#555b68")
fig.text(0.5, 0.055,
         f"0.15625 m effective finest cells · 240 frames · flat drop · Metal · all {len(rows)} runs passed without fallback",
         ha="center", fontsize=10, color="#555b68")
fig.savefig(OUT / "threshold-sweep.png", dpi=160)
fig.savefig(OUT / "threshold-sweep.svg")
