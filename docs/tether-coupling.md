# Weak tether coupling

The weak human tether uses one shared force between the existing forward hold
point and the grabbed face. Its player mass is 128 corner-mass units: sixteen
isolated standard voxels. Connected-group mass is cached by target identity,
voxel count and collision-topology generation. Sleeping groups wake before
measurement, and fixed particles make the group an anchored target.

Unique corner particles are counted with their material ownership weights.
A welded corner represents the material of every cell referencing it. Omitting
this weight incorrectly makes welded structures lighter and introduces apparent
momentum errors after the shape solver projects them. The CPU uses `glue_count`
and the GPU's solid-particle reference count carries the same weight.

The reduced-mass spring uses angular frequency 10 rad/s and critical damping.
Force caps limit player acceleration to 30 m/s² and free-group acceleration to
60 m/s²; relative closing speed is capped at 12 m/s. The player receives impulse
`J`, and the selected movable face particles share `-J`. CPU and GPU apply the
same prepared force, including its ownership weighting. The upgraded tether
and bots retain their existing integrator behavior.

## Scheduling and movement

Coupling is prepared once per simulated fixed tick, outside the GPU attempt and
its CPU fallback. The force remains immutable across the tick's two substeps;
fallback therefore reapplies voxel integration without replaying player recoil.
Catch-up ticks commit separately while weak tethers are active. Fixed anchors
also run coupling when the world has no dynamic particles.

GPU particles use a weak-player bit mask alongside their existing legacy tether
owner, allowing simultaneous weak pulls on one corner. This reuses the existing
tether buffer and uniform array; there are no additional buffer bindings or
per-frame heap allocations. The Metal entry point is versioned to reject an
old precompiled kernel with the previous tether encoding.

Normal walking brakes are 400 m/s². They remain active for light targets. Heavy
groups (at least player mass) and fixed anchors use nonzero 12 m/s² ground drag
while held, allowing a 30 m/s² grapple to pull the player along the ground.
Airborne weak-grapple carry bypasses no-input air braking and uses the existing
17 m/s horizontal grapple ceiling. Release stops force without canceling carry;
landing clears carry and restores normal brakes. Reset, true death and revival
clear the carry state.

Keyboard, gamepad and network movement share a walking acceleration limiter.
Holding a tether no longer lets movement input accelerate to the 17 m/s grapple
ceiling: the walk motor stops adding speed at the normal 5 m/s. If physical
pulling has already produced a higher speed, input can steer or brake that
momentum but cannot increase its magnitude. Actual tether impulses and airborne
release momentum retain their existing behavior. This fixes the reproduced
5-to-17 m/s walking boost without changing the paired physics impulse.

Protocol 28 adds target identity, point velocity and inverse mass to weak-tether
snapshots, plus an airborne-carry flag in player state. Client prediction uses
the same impulse function, never modifies authoritative voxels, and restores
carry/ground state before replaying unacknowledged inputs. Target velocity is
authoritative rather than estimated from visual-position differences.

## Validation

```powershell
./tools/validate_tether_coupling.ps1 -Backend 43
./tools/validate_tether_coupling.ps1 -Backend 33
```

Tests cover shared impulse balance, one-cube 16:1 velocity response, connected
and shared particles, two players pulling the same group, anchored motion,
upward grappling, zero-Matter eligibility, loose-cube ground braking, airborne
release and landing, collisions, periodic distance, anchor deletion, prediction
and lifecycle resets. The isolated physical test disables gravity and global
velocity damping in its own executable, retaining the full shape solver.
Walking tests compare untethered, tethered and predicted input, then exercise
forward, sideways and diagonal movement with a settled loose cube at 30/60/144
FPS. They also check steering and braking while carrying airborne momentum.

Measured momentum error, relative to the exchanged impulse:

| Group | Players | OpenGL compute | Single-thread CPU |
| --- | ---: | ---: | ---: |
| 1 cube | 1 | 0.00371% | 0.00369% |
| 4 welded cubes | 1 | 0.01035% | 0.01197% |
| 16 welded cubes | 1 | 0.00291% | 0.00553% |
| 4 welded cubes | 2 | 0.00676% | 0.00650% |

The pass threshold is 1%. A forced GPU validation failure repeats these checks
through CPU fallback and verifies that recoil is not applied twice.

With normal gravity and damping, the tested fixed anchor pulls about 4.06 m in
one second. Settled loose-cube hold drift is below 0.0001 m over the final second.
The 30-vs-144-FPS grounded-grapple position difference is about 0.035 m at one
second. Real localhost host/client tests cover both anchored and loose targets,
explicit mass replication, release and prediction. Existing drain and recent
playtest regressions also pass.

Windows CPU/OpenGL 3.3 and GPU/OpenGL 4.3 are supported by the harness. Metal
source is updated; Metal compilation and execution are unavailable on this
Windows machine. Gravity, global damping, contact forces and input are external
to the tether; conservation measurements intentionally isolate them.
