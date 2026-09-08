# Greedy cube activation

`--activation-cubes=greedy` enables CPU-only multiscale activation. Add
`--physics=cpu-st` or `--physics=cpu-mt`; GPU backends reject this mode.

The activation path collects the selected six-connected static component,
covers each compatible voxel type with deterministic integer cubes, creates the
ordinary unit child voxels, and binds their particles to the cube controls.
Children remain responsible for rendering and collision geometry. The current
implementation therefore reduces integrated degrees of freedom and VGS work,
but does not reduce child voxel allocation.

Run the large validation scene with:

```sh
tools/build_sized_tests.sh
.build/bin/fps_ray_sized_cpu --physics=cpu-mt \
  --debug-scenario=greedy-activation-irregular --debug-steps=600 \
  --debug-gif-fps=12 --debug-output=.build/greedy-demo
```

The harness holds the structure static for 20 fixed steps, then activates it.
Child colors identify their controlling cube size and white wireframes show
coarse cube bounds. `greedy-activation-irregular-fine` runs the same structure
with ordinary fine degrees of freedom for comparison.

Deletion, fracture, and other topology changes materialize the coarse island
back to fine particles first. The handoff corrects velocity to preserve total
linear and angular momentum. A later activation may therefore materialize an
existing coarse island before binding the new one.
