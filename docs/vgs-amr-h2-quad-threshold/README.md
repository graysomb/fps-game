# Adaptive cube drop with a threshold that quadruples at each finer level

This follow-up uses the h²-scaled VGS material from `coarse3-inv-scale` and
changes the curvature threshold factor from 2 to 4. VGS breaking is disabled.
The plotted root threshold `T₀` becomes `T(level) = 4^level T₀`. A terminal
parent refines when its own curvature exceeds its level threshold. A parent
coarsens only when its own curvature is below its threshold, the mean child
curvature is below the child level threshold, and no grandchildren are active.

The flat 10×10×10 static cube is activated as 32×32×32 PBD cells at 0.15625 m
effective finest size. AMR starts with the root VGS and runs after every PBD
substep. All 37 measured 240-frame Metal runs passed without fallback. This
historical solver has no active-cell cap.

The [active-resolution chart](active-resolution.png) records rendered terminal
cells at frames 120, 140, and 240. It includes a logarithmic full sweep and
detail views around the fractional root thresholds where the response changes.
At frame 240, root threshold 0.1 retains 29,849 cells, 0.5 retains 7,099,
0.9 retains 64, and 1 retains eight. The response has small nonmonotonic
regions, so lines connect measurements without smoothing.

The [combined CSV](threshold-sweep-all.csv) contains every plotted measurement
and refine/coarsen totals. `pilot.csv` and `detail.csv` hold the two run batches.
Run `python3 plot_sweep.py` here to regenerate the PNG and SVG.

The code on this branch supports factors 2 and 4 through
`VGS_ADAPTIVE_THRESHOLD_LEVEL_FACTOR`, with 4 as the current default.
`VGS_ADAPTIVE_CURVATURE_THRESHOLD` remains the root threshold and defaults to
50 s⁻². The sweep overrides it for each run; at the default root threshold,
the cube stays at one terminal cell in this scenario.
