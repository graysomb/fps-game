# Adaptive VGS threshold after size-squared strength scaling

The attempt to make VGS material response more consistent across voxel sizes
produces a **phase-change-like instability in the adaptive tree** for the cube
drop. Here, “phase change” describes the sudden switch in active tree resolution,
not a material phase transition. A one-unit change in the curvature threshold
separates a mostly coarse tree from one that stays mostly refined.

The test is the 240-frame `adaptive-cube-drop` scenario on Metal, with a 0.25 m
requested PBD size and 0.15625 m effective finest cell size. The threshold is
applied after every PBD substep. All tested runs passed without GPU fallback.
The solver still uses **50 s⁻²**; the other thresholds were temporary diagnostic
builds.

| Curvature threshold (s⁻²) | Terminal cells at frame 120 | Terminal cells at frame 240 | Refinements | Coarsenings |
| ---: | ---: | ---: | ---: | ---: |
| 50 | 8 | 1 | 96 | 96 |
| 35 | 8 | 1 | 113 | 113 |
| 34 | 8 | 1 | 113 | 113 |
| 33 | 2,528 | 27,434 | 71,923 | 68,004 |
| 32 | 22,212 | 27,385 | 71,269 | 67,357 |
| 25 | 31,585 | 30,402 | 64,986 | 60,643 |

The root curvature peaks at 68.6 s⁻² near frame 109. The
[curvature-versus-scale figure](curvature-vs-scale.png) shows the 99th percentile
at each level over time, alongside the root response before and after scaling.
The [threshold-sweep figure](threshold-sweep.png) shows the abrupt tree response.
Full sweep counts are in [threshold-sweep.csv](threshold-sweep.csv).

The per-level CSVs include inactive descendants whose positions follow their
active ancestors through trilinear interpolation. Their measured curvature is
therefore not an independent fine-scale signal. The CSV samples once per frame,
while AMR makes decisions every substep. The [scaled](curvature-by-level-h2.csv)
and [previous](curvature-by-level-before-h2.csv) distributions should be read
with those limits in mind.

**Implication:** the current single refine/coarsen threshold cannot produce a
stable intermediate resolution for this drop. Threshold 34 barely passes the
first refinement level; threshold 33 triggers persistent near-full refinement.
Before selecting a lower default, the tree needs a way to control the cascade
and retention of active descendants.
