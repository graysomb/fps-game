# Experimental adaptive Metal physics

Adaptive physics is **off by default**. Enable it explicitly with
`--physics=gpu-metal --adaptive-physics=on`; use `--adaptive-physics=off` for the
matched fine path. No workload is automatically promoted to adaptive mode.
Fewer leaves are not evidence of lower total physics cost.

## Implementation

`adaptive_mesh.c` is a CPU test oracle, not the runtime mesher. It verifies
coverage and uses exhaustive positive-area face adjacency for balancing.
`adaptive_metal.m` and `shaders/adaptive/mesh.metal` implement the resident mesh.
Leaves carry a fine-grid Morton origin, level, material, activation domain, and
stable fine-cell source. The supported coordinate interval is `[-2^20, 2^20)`.
Coordinates are rest coordinates; current deformed positions never determine
sibling identity. Duplicate coordinates within a domain invalidate the build.

Initial occupancy is stably radix sorted by domain and Morton origin. Eight
aligned, compatible, unprotected siblings merge through hierarchical exclusive
scan and scatter. Replacement and ordered child expansion preserve Morton
order, so subsequent leaf sorts are unnecessary. Surface cells and all existing
collision samples stay fine. Material, tether, fixed, sleeping, broken, and
already substantially deformed cells are protected by the integration adapter.

Initial balancing uses six privately owned compact leaf-ID requests per leaf,
compaction, radix sorting, uniqueness, and owner-computed splits. Fine-side
queries detect a coarse neighbor without sampling its large face. Fracture
balancing uses an exact owner-computed search of every fine tile on a coarse
face, not four quadrant samples. Schedules are bounded by the supported depth;
inactive stages skip data work, and final GPU validation rejects remaining
violations. A split propagation chain travels toward successively coarser
levels, so its length cannot exceed that depth. The independent oracle tests
include an off-center tiny neighbor requiring multiple rounds.

`adaptive_physics.m` and `shaders/adaptive/physics.metal` integrate independent
controls and evaluate one sized VGS shape per leaf. Existing integration,
collision radii, contact kernels, fine fracture evaluation, and VGS projection
equations are reused. Fine interior shapes are not also solved. Trilinear
dependent points carry no duplicate inertia; their mass is distributed once
to independent controls. Two Jacobi attachment passes follow each shape pass.
Constraints write private contribution slots; controls gather sorted incident
contributions. New meshing/interface/topology kernels contain no atomics.
Existing solver contact atomics remain.

Stable shape slots and per-domain contribution ranges keep unchanged control,
mass, and interface records resident when another domain fractures. GPU
break detection determines whether a cached Metal indirect command buffer runs.
Unchanged topology executes a zero-length rebuild range. Refinement, balancing,
state interpolation, and per-domain linear/angular momentum corrections finish
before the next substep consumes the replacement. There are no new command
submissions, CPU waits, count downloads, or CPU convergence polls. Diagnostics
are consumed after the existing frame-end completion boundary.

The current branch uses VGS-as-glue: fracture disables a fine voxel's shape
constraint. This implementation preserves existing particle/connectivity IDs;
it never welds coincident positions by coordinate alone. It does not import
`coarse`'s global pointer topology or CPU fracture/recoarsening implementation.

## Build and validation

The isolated build requires the existing raylib macOS build produced by the
project's normal build setup. It does not overwrite the tracked app bundle.

```sh
sh tools/build_adaptive_tests.sh
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 .build/adaptive/mesh_test
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 .build/adaptive/physics_test
python3 -m unittest discover -s tests -p test_adaptive_gate.py
python3 tools/validate_gpu_layout.py
```

Mesh tests compare GPU leaves to the independent CPU oracle and check exact
occupied volume, material, alignment, source identity, and non-overlap. Cases
include holes, negative coordinates, protection, shuffled sparse inputs,
separate domains, hierarchical scans, and invalid inputs. Physics tests cover
shared particle IDs, distinct coincident IDs, repeated GPU refinement,
unchanged-domain state, rest-state constraints, relative mass error `1e-5`, and
relative linear/angular momentum error `1e-4`.

Game debug runs retain existing scenario assertions and additionally check
contact sample retention, independent simulation membership, mass conservation,
and interface separation bounded by `0.05` fine-cell widths. These checks do not
establish equivalence of interior deformation.

## Measurements

```sh
python3 tools/benchmark_adaptive.py \
  --output .build/adaptive/my-measurement --runs 10 --steps 120
```

Each output directory must be new. The runner snapshots the executable and
shaders, records their hashes, performs a separate warm-up run for each mode,
then runs at least ten randomized-order pairs with identical settings. Reports
and every frame's JSONL timings remain alongside `results.json`. Construction
and topology events are included in amortized physics cost; slow initial frames
are not discarded. Event p95 uses the union of event steps from both modes.
The gate requires a negative upper endpoint of the paired bootstrap 95% CI for
physics time, no statistically supported frame/event regression, no extra
synchronizations, successful Metal correctness runs, and actual coarsening.
Exit status 1 means at least one workload failed the gate, not necessarily that
the executable crashed. Inconclusive measurements require more pairs.

Physics timing includes CPU packing, reconstruction/commit, and GPU waiting.
`frameMs` also includes the harness's voxel draw submission; it is not a
display/presentation latency measurement. Correctness validation and PNG
capture are outside those intervals. Both modes use the same capture settings.
The JSONL includes GPU duration, CPU wait/pack/upload/readback cost, transfer
bytes, dispatches, independent controls, leaf count, generation, and scratch
allocation size.

For separate stage attribution:

```sh
FPS_ADAPTIVE_BENCHMARK=1 FPS_ADAPTIVE_PROFILE=1 \
  .build/adaptive/fps_ray_gpu --physics=gpu-metal --adaptive-physics=on \
  --debug-scenario=adaptive-solid --debug-steps=120 --debug-capture-steps=0 \
  --debug-output=.build/adaptive/profile
```

`stageGpuTicks` records construction, shape, interface, contact, and adaptive
fracture/rebuild timestamp intervals, in that order. Availability is explicit.
Counter resolution stays in the existing command buffer. On devices lacking
dispatch-boundary sampling, profiling inserts encoder boundaries, so the
paired speed runner always disables this instrumentation. These stage intervals
are diagnostic, not an exhaustive partition of end-to-end cost.

## Release limitations

This is an experimental implementation, not a completed release qualification:

- CPU world repacking currently recreates the adaptive context. Per-domain
  caching applies within a resident generation; it does not yet survive a
  global CPU repack. The old fine particle backing storage and gameplay
  reconstruction are retained. GPU integration is reduced; backing allocation
  is not.
- Rest-space domains currently use the existing activation/collision group at
  upload. A fully independent persistent arena lifecycle and fine-bond
  connectivity/corner reconstruction across global repacks remain future work.
- Broken fine cells remain protected, and resident fracture only refines.
  There is no resident recoarsening. Whole-world repack can rebuild unaffected
  regions; a permanent per-active-lifetime refinement history is not yet stored.
- Failed GPU generations are rejected before canonical CPU commit and use the
  existing backend failure handoff. Per-arena allocation-failure fallback and
  retaining a separate last-complete GPU topology are not implemented.
- Existing fine collision samples/radii are retained, but comprehensive paired
  contact-penetration and interior-deformation qualification remains necessary.
  `adaptive-irregular` exercises cavities/steps; `adaptive-coarse-irregular`
  reuses the historical `coarse` tunnel-and-tail occupancy with this branch's
  normal activation and solver path.
- No OpenGL adaptive implementation or measured eligibility classifier is
  enabled. A failing performance gate must keep the fine default.

The existing `tether-active-then-static` 120-step test also fails its pull-distance
threshold on the untouched starting binary. It must not be reported as a newly
introduced adaptive regression or silently treated as a passing release gate.
