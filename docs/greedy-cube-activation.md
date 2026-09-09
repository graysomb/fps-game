# Greedy cube activation

`--activation-cubes=greedy` enables CPU-only multiscale activation. Add
`--physics=cpu-st` or `--physics=cpu-mt`; GPU backends reject this mode.

The activation path collects the selected six-connected static component,
covers each compatible voxel type with deterministic integer cubes, creates the
ordinary unit child voxels, and binds their particles to the cube controls.
Children remain render geometry only. Coarse parent corners provide collision
samples with radius scaled by parent edge length; a child point participates
only when another parent uses it as an interface corner. The implementation
reduces integrated degrees of freedom and VGS work, but does not reduce child
voxel allocation.

Run the large validation scene with:

```sh
tools/build_sized_tests.sh
.build/bin/fps_ray_sized_cpu --physics=cpu-mt \
  --debug-scenario=greedy-activation-irregular --debug-steps=620 \
  --debug-gif-fps=12 --debug-output=.build/greedy-demo
```

The harness holds the structure static for 20 fixed steps, then activates it and
executes 600 measured physics steps.
Child colors identify their controlling cube size and white wireframes show
coarse cube bounds. `greedy-activation-irregular-fine` runs the same structure
with ordinary fine degrees of freedom for comparison.

## Current benchmark

Measured on an Apple M2 Pro with the `-O2` CPU build using
`python3 tools/benchmark_greedy.py`. Each case runs once as warm-up and three
times for measurement. The table reports the median and range for 600 fixed
steps after the 20-step preview; setup, rendering, validation, and activation
are excluded.

| Backend | Fine median (range) | Greedy median (range) | Greedy speedup |
|---|---:|---:|---:|
| CPU ST | 2,236.55 ms (2,232.59–2,239.53) | 1,129.47 ms (1,123.88–1,149.37) | 1.98x |
| CPU MT | 2,060.28 ms (2,057.63–2,066.97) | 1,441.18 ms (1,436.53–1,456.21) | 1.43x |

The cover itself takes about 190 ms once at activation. Fine mode integrates
2,335 particles. Greedy mode retains 2,324 collision/render samples but
integrates only 177 controls, a 92.4% reduction, and evaluates 113 parent shapes
instead of 1,692 unit shapes, a 93.3% reduction.

The greedy path now uses pool-index lookup, persistent contact scratch, the
original parallel VGS accumulation, and one static-contact pass per substep.
Its remaining largest cost is dynamic contact generation: the CPU ST median
spends about 563 ms there, compared with 184 ms in VGS. The full per-stage
medians are written to `.build/greedy-benchmark/summary.json`.

Deletion, fracture, and other topology changes materialize the coarse island
back to fine particles first. The handoff corrects velocity to preserve total
linear and angular momentum. A later activation may therefore materialize an
existing coarse island before binding the new one.
