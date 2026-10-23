# Physical Foundry golem laboratory

Run `.build/bin/fps_ray_gpu.exe --enemy-lab=golem --physics=gpu` from the
repository root, or use the CPU executable with `--physics=cpu-st`. The golem
also follows the mixed-enemy case when cycling the ordinary Enemy Lab.
This is an offline laboratory fixture. The live Foundry boss and network
protocol are unchanged.

## Controls

| Keyboard | Controller | Action |
| --- | --- | --- |
| F3 | D-pad Right | Idle / walk / slam / sweep / hand launch |
| F4 | Y / Triangle | Recreate the complete fixture |
| F5 | X / Square | Pause / resume |
| F6 | RB / R1 | One paused simulation tick |
| F7 | Back / Select | Player invulnerability |
| F8 | D-pad Down | Toggle regeneration |
| F9 | — | Skeleton overlay |
| F10 | LB / L1 | Release the right arm |
| F11 | — | Release the left shin and foot |
| F12 | — | Release the chest cover |
| Escape | — | Return to the menu |

Normal movement, weapon, mining and tether controls remain available. The
console's existing equipment commands can supply the cutter, sword or upgraded
tether. The HUD reports attached cells, total body/debris cells, constraints,
mount error and physics time.

## Construction and simulation

`golem_body.h` is shared by the original procedural preview and the physical
fixture: **283 cells in 18 segments**, including the rectangular knee plates,
hands, neck, reactor and thin visor. There is no golem-specific voxel limit or
reduced-detail physics body. Allocation follows the authored definition and is
subject to the engine's existing world capacity.

All cells retain mass and participate in the ordinary particle solver. The
renderer draws their simulated corners. Four distributed, compliant mounting
points attach each connected segment to the invisible animated skeleton.
Soft distributed pose motors supply rotational acceleration; they avoid forcing
one cube's welds to accelerate an entire arm. Animation transforms are smoothed
over 150 ms so a visual pose transition cannot teleport a physical mounting
target. These motors disappear when their cells lose connectivity.

The 483 authored face connections preserve the original gaps and dimensions.
They use affine corner constraints and the engine's ordinary face-link graph,
so cuts, fracture, cluster selection and tether acquisition use real voxels.
Sustained estimated mount loads above eight times segment weight release the
mount after 0.12 seconds. Released parts keep their velocity and internal links;
they can collide, break further and be grabbed or thrown. Lower-arm and lower-leg
mounts depend on their upper segments. Leg loss drives the existing kneeling
animation. Adjacent articulated segments exclude self-collision; unrelated
segments, debris, actors and the world still collide.

The hand attack launches the twelve-cell right fist toward the nearest human
player, leaving the forearm attached. Its initial velocity targets 22 m/s
horizontal travel with gravity compensation. With no player target it fires
sixteen meters forward. Release and velocity assignment happen once before
physics upload, so CPU and GPU receive the same launch.

Regeneration starts after 1.25 seconds and restores a complete missing limb
assembly in one transaction. The arm's upper arm, forearm and hand remain
separate articulated groups; both legs share a repair transaction. Surviving
cells keep their identity. Every new cell and weld is ready before physics
advances. Allocation failure rolls back newly created cells and retries.

Placement uses oriented cell geometry against actors, terrain and debris.
Only overlapping, unheld debris belonging to this golem is reclaimed, with a
quarter-second cyan shrink effect. Other fallen pieces remain available to
throw. Actors, held debris and unrelated structures block repair until cleared.
The lab HUD reports waiting, blocked, rebuilding and recovering; pending cells
have correctly oriented cyan outlines (orange when blocked). Leg recovery only
starts after the assembly exists and takes 0.8 seconds. Missing geometry is
included in floor clearance so the replacement cannot spawn under the floor.

Mount loads are evaluated once per fixed tick from committed physics state,
not stale CPU particle positions midway through a GPU dispatch. The recovery
interval excludes motor-driven pose correction from overload detachment;
weapon destruction and connection cuts still work throughout recovery.
Detached cells use normal debris expiry and never enter terrain restoration
queues.

The golem uses twelve constraint iterations; additional shape passes visit only
golem cells. Ordinary gameplay retains its existing iteration count. CPU, GLSL
and Metal share the cuboid rest dimensions and constraint representation, with
the existing sixteen GPU buffer bindings. The Metal entry-point name is versioned
to reject an incompatible old precompiled library. GPU failure restores the
animation checkpoint before the existing CPU fallback runs.

## Reproduce validation and captures

```powershell
./tools/validate_golem_physics.ps1 -Backend 33
./tools/validate_golem_physics.ps1 -Backend 43 -Capture
python tools/golem_physics_media.py
./tools/validate_golem_physics.ps1 -Backend 43 -SkipBuild -StressSeconds 600
./tools/validate_golem_physics.ps1 -Backend 43 -Accuracy
```

The fixture checks exact initial cell geometry, idle stability, locomotion,
attacks, hand launch velocity and travel, kneeling, regeneration, face cuts, upgraded
tether forces, obstructed regeneration, identity tracking after compaction,
finite particle positions and teardown. The optional stress run repeatedly
detaches and regenerates limbs while exercising the normal debris recycler.

Captures are written to `artifacts/golem/physics`: initial and settled PNGs,
action and damage stills, and seven six-second GIFs generated from real physics
steps. They are separate from the earlier scripted visual-only GIFs.

The test prints full physics-step median and p95 times for the isolated 283-cell
body. This is a physics prototype with a measurable simulation cost, not a
claim about four-player Foundry performance. Metal execution requires a Mac and
cannot be validated by the Windows runner.

### Windows results, 2026-10-22

Both OpenGL 3.3/CPU and OpenGL 4.3/GPU builds passed the complete behavior
fixture. Firefight integration (including four-view rendering), default tether
grapple, static support/collision cache, packet-format and shader-layout
regressions passed. A forced GPU validation failure successfully fell back to
CPU, recovered to GPU and preserved the animation clock and complete body.

The ten-minute GPU detachment/regeneration run peaked at **341 cells** including
debris, with no terrain-restoration entries. Both Windows game binaries were
rebuilt. Metal source/layout was updated but Metal compilation and execution
were not available on this Windows machine.

Isolated release-build measurements, 360 fixed steps at 60 Hz, with no other
test/build processes running (full physics step, excluding drawing):

| Backend | Median | p95 |
| --- | ---: | ---: |
| OpenGL 4.3 compute | 5.734 ms | 8.698 ms |
| Single-thread CPU | 28.294 ms | 29.276 ms |

The CPU cost is substantial. This fixture establishes physical behavior and
visual fidelity; further solver optimization is needed before live-boss rollout.

## Hand accuracy diagnostic

`-Accuracy` fires the actual animated hand at stationary torso-height targets
6, 12, 20 and 30 meters from the golem, at -45, 0 and +45 degrees. It samples
the twelve physical hand cells at 60 Hz until their first floor contact.
After launch the target actor is disabled so its collision cannot obscure the
unimpeded trajectory measurement. This tests aiming, not player damage.
Output includes nearest fist-center distance, nearest enclosing-box distance
to the target point, and distance short of the target along the firing direction
at floor contact (negative means overshoot). Box distance is a conservative
envelope measurement, not a hit-confirmation test. Moving targets are not tested.

The current gravity-only launch is inaccurate at long range. OpenGL 4.3 center
miss distances across the three angles are 0.12–0.18 m at 6 m, 0.50–1.07 m
at 12 m, 4.43–6.02 m at 20 m, and 12.77–14.51 m at 30 m. The engine damps
velocity each physics substep; the launch calculation does not compensate for
that damping. These are diagnostic results, not passing accuracy thresholds.
Single-thread CPU reproduces the same pattern: center misses of 0.13–0.19 m,
0.52–1.09 m, 4.62–6.09 m and 13.07–14.61 m respectively. Raw results are in
`artifacts/golem/physics/accuracy-gpu.csv` and `accuracy-cpu.txt`.

### Experimental damping comparison

The validation script accepts `-Damping 0.999` (default `0.99`), builds a
separate executable named `test-43-damping-0.999.exe`, and prints its damping
value. Normal game builds retain `0.99`; this experiment does not retune them.

```powershell
./tools/validate_golem_physics.ps1 -Backend 43 -Damping 0.999
./tools/validate_golem_physics.ps1 -Backend 43 -Damping 0.999 -SkipBuild -Accuracy
./tools/validate_golem_physics.ps1 -Backend 43 -Damping 0.999 -SkipBuild -StressSeconds 120
```

GPU nearest fist-center miss, across the same three stationary target angles:

| Range | Damping 0.99 | Damping 0.999 |
| --- | ---: | ---: |
| 6 m | 0.12–0.18 m | 0.05–0.22 m |
| 12 m | 0.50–1.07 m | 0.17–0.21 m |
| 20 m | 4.43–6.02 m | 0.34–0.54 m |
| 30 m | 12.77–14.51 m | 1.13–1.42 m |

CPU at 0.999 shows the same improvement, with center misses of 0.05–0.23 m,
0.16–0.25 m, 0.36–0.72 m and 1.31–1.82 m respectively. The CPU and GPU
trajectories are similar, not identical.

Fist forward travel over the first third of a second increases from 5.83 m to
7.11 m. All nine target points at 6–20 m enter the fist's bounding box;
the three 30 m points still miss that envelope by 0.38–0.71 m. Bounding-box
overlap is not an exact actor collision test. The long-range shot still needs
damping-aware aim; merely reducing damping does not make it exact.

At 0.999 the full stability suite passes on GPU and single-thread CPU (using
the same GL 4.3 test binary with `--physics=cpu-st`): idle, walking, slam,
sweep, launched hand, arm/leg loss, kneeling, regeneration, cut links, tether
detachment, obstructed regrowth, compaction and cleanup. The 120-second GPU
stress run peaks at 340 cells including debris, with finite positions and zero
terrain-restoration queue entries. It ends with 258 attached cells and 283 total;
obstruction-aware regrowth is not asserted to finish during that stress run.
This is a golem-fixture test, not a full-map destruction or multiplayer soak.
Logs use `damping-0999-*`, `accuracy-0999-*`, and `stress-0999-gpu.txt` in
`artifacts/golem/physics`.

## Lab startup and atomic repair regression

The startup failure came from the normal debris-hit path testing inactive
combatants. The lab's cleared enemy slot retained a position inside the new
golem, and contact removed 14 body cells on the first frame. Inactive slots now
skip debris hits. Mounted golem cells use physical collisions and are not
consumed as disposable projectiles when a player touches them.

Run the actual lab initialization plus gameplay physics/maintenance loop with:

```powershell
./tools/validate_golem_physics.ps1 -Backend 43 -Lab
./tools/validate_golem_physics.ps1 -Backend 33 -Lab
./tools/validate_golem_physics.ps1 -Backend 43 -Repair -Capture
./tools/validate_golem_physics.ps1 -Backend 33 -Repair
```

`-Lab` checks 60 seconds each at 30, 60 and 144 FPS plus catch-up frames.
`-Repair` retains fallen debris and checks complete segment connectivity,
both-leg regrowth, arm/hand/chest/head replacement, partial allocation rollback,
individual-cell loss, held debris and actor obstruction. The regular suite also
checks terrain obstruction, weapons/tether detachment and launched-hand travel.
Leg-regrowth captures are in `artifacts/golem/physics/leg-regrowth`.

Validated on Windows after the repair fix: OpenGL 3.3/CPU and OpenGL 4.3/GPU
passed all three 60-second lab runs and catch-up frames, the original action
suite, and the atomic repair suite. Damage during standing recovery and
canceling a debris dissolve passed on both backends. A forced GPU validation
failure fell back to CPU, completed regeneration and recovered to GPU. The
four-player Firefight/Crosscut integration regression passed. Both packaged
Windows game executables were rebuilt; Metal runtime validation remains
unavailable on this machine.

Run `python tools/golem_physics_media.py` after the capture commands to produce
`physical-leg-regrowth.gif` alongside the existing action GIFs.
