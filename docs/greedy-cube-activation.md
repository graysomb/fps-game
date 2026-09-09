# Greedy cube activation

`--activation-cubes=greedy` enables multiscale activation in every physics
backend. Select `--physics=cpu-st`, `--physics=cpu-mt`, `--physics=gpu-gl43`, or
`--physics=gpu-metal`. A platform build still exposes only its native GPU API.

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
.build/bin/fps_ray_sized_cpu --activation-cubes=greedy --physics=cpu-mt \
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
| CPU ST | 2,242.56 ms (2,236.27–2,244.56) | 979.47 ms (975.15–980.83) | 2.29x |
| CPU MT | 2,019.33 ms (1,982.02–2,108.31) | 1,293.57 ms (1,290.81–1,294.73) | 1.56x |
| Metal | 2,700.35 ms (2,687.83–2,700.98) | 3,697.08 ms (3,693.52–3,712.95) | 0.73x |

The cover itself takes about 190 ms once at activation. Fine mode integrates
2,335 particles. Greedy mode retains 2,324 collision/render samples but
integrates only 177 controls, a 92.4% reduction, and evaluates 113 parent shapes
instead of 1,692 unit shapes, a 93.3% reduction.

The greedy path now uses pool-index lookup, persistent contact scratch, the
original parallel VGS accumulation, and one static-contact pass per substep.
Its remaining largest CPU ST cost is dynamic contact generation: the median
spends about 573 ms there, compared with 185 ms in VGS. The full per-stage
medians are written to `.build/greedy-benchmark/summary.json`.
Pass GPU backends and the matching executable explicitly, for example
`python3 tools/benchmark_greedy.py --binary .build/bin/fps_ray_gpu
--backends gpu-metal`. The script rejects failed fixtures rather than recording
a fallback as a GPU result.

Use `--fixed-step-batch=1..8` to measure multiple fixed steps in one GPU
command buffer. The debug harness exposes the same control as
`--debug-physics-batch=1..8`. Normal gameplay already passes every fixed step
accumulated for the current rendered frame to one `gpu_physics_steps` call.

The OpenGL 4.3 and Metal pipelines consume the same flattened dependency rows
built by `greedy_coarse_bind`. Direct particles keep the original kernels;
mapped parent corners expand corrections into their independent sources and
refresh derived positions between collision and VGS stages. Passive children
are refreshed only after finalizing the controls. Collision hash reach is
computed from the largest participating parent radius. In resident mode the
dependency buffers remain uploaded until topology changes.

Greedy floor contact is separate from static-grid contact. At topology upload,
the GPU packer builds constraint islands and assigns each floor sample a compact
conflict batch. An island with at most 512 controls is solved by one 128-thread
workgroup: controls are cached in 16 KiB of threadgroup memory, interface points
are derived directly from their source rows, and barriers preserve batch order.
The controls are written back once after the island finishes. Larger islands use
the global-memory batch path. The same analytic pass handles world bounds and
player boxes; terrain and static blocks retain the general static-contact kernel.

On the irregular fixture, the floor pack contains one island, 177 controls, 258
contact samples, and 38 batches. With eight fixed steps per command buffer, an
A/B run reduced Metal dispatch encoding from 97.1 ms to 30.6 ms and total debug
physics time from 3,419.8 ms to 2,924.5 ms. The final positions differed from
the retained fallback by at most 0.00045 voxel units.

The 600-step irregular fixture passes on native Metal with 177 controls, 2,147
derived points, zero measured dependency error and floor penetration, and less
than `0.003` final parent strain. OpenGL and Metal share the validated buffer-slot
and pipeline-mode contract; native OpenGL execution must be run on a GL 4.3
platform.

The Metal path now emits child instance matrices from the final GPU particle
state. On Apple silicon, `DrawMeshInstanced` consumes the shared Metal buffer
directly; the CPU readback contains only the 177 independent controls. OpenGL
4.3 executes the same matrix kernel and currently stages its result for the
OpenGL renderer.

Batch size materially changes the crossover on this fixture. Three measured
runs on the same M2 Pro produced:

| Fixed steps per command buffer | Fine median (range) | Greedy median (range) | Greedy speedup |
|---:|---:|---:|---:|
| 1 | 2,700.35 ms (2,687.83–2,700.98) | 3,697.08 ms (3,693.52–3,712.95) | 0.73x |
| 8 | 3,672.77 ms (3,655.55–4,315.08) | 2,959.69 ms (2,956.29–2,991.17) | 1.24x |

The eight-step batch reduces greedy time by about 20%. Fine physics becomes
slower with such a large command buffer on this workload, so batch size should
remain a measured backend setting rather than a universal default.

Deletion, fracture, and other topology changes materialize the coarse island
back to fine particles first. The handoff corrects velocity to preserve total
linear and angular momentum. A later activation may therefore materialize an
existing coarse island before binding the new one.
