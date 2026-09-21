# Adaptive cube drop with a threshold that doubles at each finer level

This experiment starts from the h²-scaled VGS material in commit `407a10d`.
VGS shear strength `alpha`, stretch strength `1-beta`, and volume restoration
strength are multiplied by `(finest PBD edge / constraint edge)^2`. Constraint
breaking is disabled.

The adaptive drop activates the flat 10×10×10 static cube as a 32×32×32 PBD
cube with 0.15625 m effective finest cells. Only the root VGS is active at
the start. Adaptive refinement and coarsening run after every PBD substep.
All 38 threshold runs used Metal for 240 frames, passed, and had no fallback.

For a root threshold `T₀`, the threshold at level `l` is `2^l T₀`.
A terminal node refines when its own curvature exceeds its level threshold.
A parent coarsens when its curvature is below its own threshold, the mean of
its eight children's curvatures is below the child level's threshold, and it
has no active grandchildren. The default root threshold remains 50 s⁻²;
the other values were diagnostic builds.

The [chart](threshold-sweep.png) plots rendered terminal cells at frames 120,
140, and 240, plus total refinements and coarsenings. It includes all root
thresholds from 1 through 35, plus 40, 45, and 50. The dense 10–30 panel
shows that threshold 11 retains 15,149 cells at frame 140 and 428 at frame
240. At threshold 12, it has eight cells at frame 140 and one at frame 240.
The earlier late-time transition between thresholds 6 and 7 also remains:
12,314 versus one cell at frame 240.

The plot uses logarithmic vertical axes. The [combined CSV](threshold-sweep-all.csv)
contains every plotted measurement. The three source CSVs separate the
original reference points, low-threshold detail, and added 11–29 points.
Run `python3 plot_sweep.py` in this directory to regenerate the PNG and SVG.
