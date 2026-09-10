# Greedy cube activation

`--activation-cubes=greedy` enables multiscale activation in every physics
backend. Select `--physics=cpu-st`, `--physics=cpu-mt`, `--physics=gpu-gl43`, or
`--physics=gpu-metal`. A platform build still exposes only its native GPU API.

The activation path collects the selected six-connected static component,
covers each compatible voxel type with deterministic integer cubes, creates the
ordinary unit child voxels, and binds their particles to the cube controls.
Passive children remain render geometry only. When a child point is used as a
neighboring parent corner, topology construction promotes that shared point to
an independent particle and adds one interpolation attachment to the owning
coarse cube. Parent VGS and contacts therefore use the original direct particle
path. Coarse parent corners provide collision samples with radius scaled by
parent edge length. The implementation reduces integrated degrees of freedom
and VGS work, but does not reduce child voxel allocation.

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

## Pre-promotion benchmark baseline

These measurements describe the former mapped-correction implementation and are
kept as a comparison point. They must be rerun before quoting performance for
the promoted-interface implementation. Measured on an Apple M2 Pro with the `-O2` CPU build using
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

The former greedy path used pool-index lookup, persistent contact scratch, the
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

The OpenGL 4.3 and Metal pipelines consume the same compact rows built by
`greedy_coarse_bind`. Every parent corner is now a direct solver particle. An
interface corner keeps a short row only for its interpolation attachment to its
owner; passive child rows are used only by final rendering refresh. VGS and
collision kernels no longer traverse dependency rows or atomically scatter a
corner correction. Collision hash reach is computed from the largest
participating parent radius. In resident mode the rows remain uploaded until
topology changes.

Greedy floor contact is separate from static-grid contact. At topology upload,
the GPU packer builds constraint islands and assigns each floor sample a compact
conflict batch. An island with at most 512 controls is solved by one 128-thread
workgroup: controls are cached in 16 KiB of threadgroup memory, promoted interface points
are handled as direct controls, and barriers preserve batch order.
The controls are written back once after the island finishes. Larger islands use
the global-memory batch path. The same analytic pass handles world bounds and
player boxes; terrain and static blocks retain the general static-contact kernel.

In the promoted topology, attachment edges are included when building floor
islands, so floor corrections and their owning coarse controls stay in the same
island. Interface contacts are direct contacts. GPU attachments are ordinary
Jacobi constraints: all promoted rows evaluate in parallel, atomically emit
mass-weighted endpoint corrections into the original correction buffer, and use
the original correction-apply kernel. A topology-time attachment-incidence
diagonal preconditioner prevents high-valence coarse controls from being
under-relaxed by the apply kernel's averaging. The GPU runs two attachment
iterations after each main VGS iteration and 24 after the final contact/VGS
pass. Apply clears the correction buffer, so redundant reset dispatches between
these passes are elided. Attachment conflict batches and attachment-specific
threadgroup solvers are no longer built or executed.

For native Metal profiling, `FPS_METAL_GPU_TIME=1` reports command-buffer GPU
time without inserting counter barriers. `FPS_METAL_STAGE_PROFILE=1` samples
each compute encoder and reports time by solver mode; this deliberately changes
scheduling and should be used for relative stage shares rather than headline
timing. `FPS_GPU_DEPENDENCY_STATS=1` reports direct VGS fan-in and compact attachment
work at topology upload.

The former mapped topology emitted 10,316 atomic additions per VGS pass and its
hottest control received 197 contributions. The promoted irregular topology has
904 direct parent-corner references across 113 shapes and emits 3,616 VGS atomic
additions; its hottest VGS particle receives six contributions. It adds 163
interpolation attachments with 474 compact source entries. The headless 600-step
CPU invariant keeps attachment error below `1e-5` unit edges and preserves the
full material mass. The same invariant executes 600 actual Metal steps in
eight-step command buffers, checks each batch, and keeps attachment error below
`1e-5` unit edges with finite state. A direct GPU-timestamp run took 1,786.62 ms
on the development M2 Pro; this is an invariant-run measurement rather than a
fine-versus-greedy gameplay benchmark. OpenGL and Metal share the validated
buffer-slot and pipeline-mode contract; native OpenGL execution must be run on
a GL 4.3 platform. Visual validation still requires an unlocked desktop session.

The corresponding stage-profile run reported 3,057.28 ms of sampled GPU time.
Parallel attachment evaluation used 1,466.82 ms (47.98%), and its shared apply
passes used part of the 535.62 ms (17.52%) correction-application total. Pair
contacts used 444.36 ms (14.53%), while direct VGS used 167.26 ms (5.47%). Stage
profiling inserts counter sampling and is intended for proportions rather than
headline timing. Compared with the final ordered Gauss–Seidel implementation's
2,641.54 ms direct run, ordinary Jacobi constraints reduce this invariant run by
32.4%.

The Metal path now emits child instance matrices from the final GPU particle
state. On Apple silicon, `DrawMeshInstanced` consumes the shared Metal buffer
directly; the CPU readback contains only the 340 independent controls on the irregular fixture. OpenGL
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

Deletion and topology changes temporarily materialize the coarse controls into
their real unit children. The handoff corrects velocity to preserve total
linear and angular momentum. Activation then binds the combined old and newly
activated children. Child membership is stored explicitly and recovered through
stable voxel identities, so voxel-pool compaction no longer invalidates a group.

Wake events rebuild the current cover at the safe point between physics steps.
The rebuilt interpolation rows use immutable rest-grid coordinates, which
preserves a group's current translated or rotated pose. Fracture is evaluated
on the refreshed unit children. When a child crosses the existing
strain/shear/hinge thresholds, the two unit children adjoining each broken face
leave the cover and remain awake. The implementation temporarily materializes
the controls, runs the original pointer-based `break_face_link` particle-clone
path, and immediately applies the greedy cover again to every unaffected child.
The local seam remains ordinary fine physics while the rest of the structure
stays coarse. Later cracks repeat the same local refinement. CPU fracture
handoff occurs at the end of the current substep. GPU fracture handoff occurs
after the current fixed-step command-buffer batch, invalidates resident
topology, and uploads the mixed fine/coarse topology on the next batch.

Run the deterministic visual fracture fixture with:

```sh
.build/bin/fps_ray_cpu --activation-cubes=greedy --physics=cpu-mt \
  --debug-scenario=greedy-fracture-demo --debug-steps=300 \
  --debug-gif-fps=20 --debug-output=.build/greedy-fracture-demo
```

Red children are the locally refined fracture seam. The other colors identify
the sizes produced when the unaffected region is covered again.
