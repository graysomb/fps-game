# Enemy Lab golem combat

The physical golem now starts an autonomous local encounter in Enemy Lab. The
existing Firefight Foundry and network protocol are unchanged. This is not yet
an online boss. All 283 authored cells and existing articulated connections
are retained.

## Playing

F3 / D-pad Right cycles Combat, Idle, Walk, Slam, Sweep, Hand Launch, then Combat.
F4 / top face button resets the encounter. F5 / left face pauses, F6 / RB steps,
F7 / Back toggles invulnerability (off initially), F8 / D-pad Down toggles limb
repair, F9 shows the skeleton, F10 / LB removes an arm, F11 removes a shin, and
F12 removes chest armor. Escape exits. The existing console `give` commands can
supply weapons. The compact HUD displays separate shield and reactor bars.

Destroy limbs normally or pull them off with the tether. Leg loss kneels the
boss and opens its chest. Shoot the exposed reactor: 150 Matter shields 200
permanent reactor health. Damage beyond remaining shield carries into health.
The shield regenerates at 15 Matter/second after four seconds without a reactor
hit. Limb hits do not postpone recharge. Glob and tether drain only the shield;
glob suppresses recharge for five seconds. Railgun penetration and rocket
occlusion retain their normal rules. Health never regenerates. Reactor death
stops all repairs and releases the whole body as physical debris.

The boss tracks living humans and holograms, navigates around intact lab
structures at 2 m/s, and turns at 90 degrees/second. Slam, sweep, and launched
hand deal 40, 25, and 50 Matter respectively, once per attack per victim. Windups
are 1, 0.8, and 1.2 seconds; cooldowns are 4, 3, and 7 seconds. Melee attacks step
closer during the early windup and lock aim for the final quarter second.
Existing exposure/execution and invulnerability rules apply. The hand launch
solves the actual damped substep trajectory instead of changing game damping.

Unclaimed golem fragments collide physically but cannot multiply an attack's
damage through the generic debris handler. Player-thrown pieces retain normal
debris damage. Detached pieces expire after 45 unheld seconds; holding a piece
resets this timer. This bounds long encounters without reducing the body design.

## Implementation boundaries

`golem_combat.inc` owns decisions, shield/health, telegraphs, deduplication, and
physical hand tracking. The controller supplies pose inputs to the existing
renderer and articulated solver. The reactor reserves one non-rendered enemy
slot for shared weapon targeting; it does not run generic bot AI, movement, or
actor collision projection. Oriented cuboid ray tests use simulated body corners
and a small seam tolerance. Chest armor remains genuine physical cover.

Combat preparation and damage resolution bracket a committed fixed step.
GPU attempts only sample poses, so fallback cannot replay attacks, launches,
shield regeneration, or damage. Manual fixtures keep the original animations;
combat uses an optional reach target and lower sweep stance to reach rat height.
Commanded attack motion receives mount-load grace, as repairs already do;
held segments still evaluate overload, and structural cuts always apply.
Shield effects use ordinary depth-tested geometry, work without bloom, and
allocate no runtime meshes or materials.

## Reproduction

- `./tools/validate_golem_combat.ps1 -Backend 43`
- `./tools/validate_golem_combat.ps1 -Backend 33`
- `./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -Fallback`
- `./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -Preview`
- `python tools/golem_combat_media.py` (requires Pillow)
- `./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -StressSeconds 600`
- Existing `validate_golem_physics.ps1` covers manual poses and atomic repairs.

Captures go to `artifacts/golem/combat`. The stress test reports whole fixture
frame median/p95, early/late means, peak cells, and action counts; its cell-count
assertion is a test guard, not a production voxel cap.

## Validation results (Windows, 2026-10-25)

- OpenGL 3.3 CPU and OpenGL 4.3 GPU combat suites passed, including actual
  single-hit slam/sweep damage and physical hand hits at 6, 12, 20 and 25 m.
- Shield overflow, nonlethal glob/tether damage, suppression, reactor persistence,
  rocket cover protection, invulnerability, friendly fire, limb cancellation,
  chest opening, hologram targeting, reset and exit checks passed.
- Forced GPU validation failure around hand launch produced one launch and one
  hit with one clock advance. CPU fallback and subsequent GPU recovery passed.
- Existing GPU manual-pose/physics regression and CPU atomic repair tests passed,
  including allocation rollback, held debris, actor obstruction and recovery damage.
- A 600-second GPU simulation run completed 60 slams, 45 sweeps and 59 launches.
  Peak voxel count was 606; final count was 534. Every sampled transform remained
  finite. No body-detail cap was introduced.
- Whole fixture simulation/update time (rendering excluded, O1 test build):
  **8.383 ms median, 9.241 ms p95**. Mean of the first and last two minutes was
  **8.667 ms and 8.729 ms**, respectively; no progressive slowdown was observed.
  This is a lab measurement, not a full Firefight or split-screen benchmark.
- Updated optimized Windows CPU/GPU executables and visually reviewed attack,
  reactor/shield and defeat captures are available. LAN and live Firefight boss
  integration remain outside this local-lab phase.
