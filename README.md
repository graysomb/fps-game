


This is a vibe coded physics engine I made. Its kinda like Minecraft but fully simulated. I wanted it to be an FPS game but i'm not sure if thats going to work out. Writing the physics was fun though! its currently a jacobi style PBD simulation with [voxel gram-schmidt shape corrections](https://dl.acm.org/doi/10.1145/3677388.3696322) and glue. The voxels have a life cycle where they are greedily meshed (just like minecraft) until physically active, after which they enter the PBD loop.

Currently the game motives interaction with the physics by requiring environment based kills or melee

runs in linux and windows so far

# Warning SLOP


# Raylib Split-Screen FPS with Voxel PBD Physics

This project is a split-screen First Person Shooter (FPS) prototype built using [Raylib](https://www.raylib.com/). It features a custom Voxel-based Position Based Dynamics (PBD) physics engine that supports destructible environments, multiscale voxels, and structural integrity.

## Physics Engine Overview

The core of the game is a custom physics engine implemented in `fps_ray.c`. It uses a **Position Based Dynamics (PBD)** approach, which modifies particle positions directly to satisfy constraints, resulting in stable and controllable simulations.

### The Simulation Loop
The physics simulation follows a standard PBD loop:
1.  **Prediction:** Particles move based on their current velocity and external forces (gravity).
2.  **Constraint Solving:** The solver iteratively corrects particle positions to satisfy physical constraints (Shape, Glue, Collision).
3.  **Integration:** Velocities are updated based on the corrected positions.

### Physics backend architecture

The solver has four platform implementations behind three runtime tiers: OpenGL
4.3 compute on Windows/Linux, Metal compute on macOS, a persistent CPU worker pool,
and single-threaded CPU. They implement the same ordered PBD stages. The CPU-owned
voxel and particle arrays remain canonical, so a failed GPU
dispatch can be validated, read back, and continued by the CPU without restarting
the frame.

The GPU path packs only particles referenced by active voxels into SSBOs. Stable
particle-pool IDs are mapped to compact indices, preserving shared corners created
by glue. Static voxels are uploaded only when the static-grid generation changes.
Each substep dispatches prediction, dynamic hash construction, scene and particle
collisions, Jacobi apply passes, VGS iterations, static collisions, glue-break
evaluation, GPU particle splitting, topology rebuilding, wake propagation, and
velocity finalization. The compact result is committed back to the canonical CPU
arrays once per rendered-frame batch.

The CPU PBD loop uses a lightweight, persistent thread pool to avoid per-frame
thread creation.
Key ideas:
*   **Parallel-for dispatcher:** A small worker pool runs chunked ranges over arrays (particles/voxels).
*   **Jacobi-style accumulators:** Constraints write into per-particle correction accumulators (`corr_sum`, `corr_weight`) using atomic float adds. A later “apply” pass updates positions, preserving determinism within a single iteration.
*   **Stable particle snapshots:** Each constraint iteration snapshots `active_particles` into a local array before building the spatial hash and gathering collisions, avoiding concurrent mutation while workers read.
*   **Thread-safe spatial hash:** Particle hash buckets are cleared and built in parallel; inserts use atomic exchange for lock-free binning.
*   **Execution order:** Integrate → Hash Clear/Build → Collisions → Apply → VGS/Apply iterations → Static collisions → Glue/lifecycle commit → Velocity update.

This layout mirrors a GPU compute pipeline and keeps contention low while remaining safe on CPU.

### Key Global Definitions & Tuning
The physics behavior is controlled by several preprocessor definitions in `fps_ray.c`. Tweaking these values allows you to trade off performance for simulation quality.

*   **`PBD_SUBSTEPS` (Default: 2)**
    *   Number of physics sub-steps per frame.
    *   **Effect:** Higher values increase simulation stability and collision accuracy but cost more CPU time.
*   **`PBD_CONSTRAINT_ITERS` (Default: 6)**
    *   Number of solver iterations per sub-step.
    *   **Effect:** Higher values make materials "stiffer" (less rubbery) but increase computational cost.
*   **`VOXEL_SIZE` (Default: 0.5f)**
    *   The world-space size of a single voxel unit.
*   **`VGS_ALPHA` (0.75) & `VGS_BETA` (0.35)**
    *   Parameters for the Voxel Gram-Schmidt (VGS) shape matching.
    *   **Effect:** Control how strictly voxels maintain their cubic shape and volume.
*   **`GLUE_BREAK_STRAIN` (Default: 0.4f)**
    *   The amount of stretching allowed before a glue bond breaks.
    *   **Effect:** Lower values make structures more fragile; higher values make them tougher.
*   **`GLUE_BREAK_HINGE_ANGLE_DEG` (Default: 20.0f)**
    *   The angle of bending allowed before a glue bond snaps.
    *   **Effect:** Controls the rotational structural integrity.
*   **`COLLISION_RELAXATION` (Default: 0.99f)**
    *   Collision response factor.
    *   **Effect:** Controls how "bouncy" collisions are (1.0 is perfectly elastic).

## Constraints & Solvers

The engine implements three primary types of constraints:

### 1. Shape Matching (VGS)
*   **Function:** `solve_voxel_shape`
*   **Description:** Ensures that the 8 particles forming a voxel maintain a cubic shape.
*   **Method:** Uses **Voxel Gram-Schmidt (VGS)** orthogonalization. It computes the best-fit rotation and volume for the deformed particles and pulls them back towards a valid cube configuration.
*   **Simplifications:** Assumes uniform mass distribution for stability.

### 2. Voxel Glue
*   **Function** `glue_neighbor_faces_for_voxel`
*   **Description** forces voxels sharing a glue constraint to share corner particles

### 3. Collision Constraints
*   **Function:** `solve_particle_collisions`
*   **Description:** Handles interactions between dynamic voxels, the static environment, and players.
*   **Method:**
    *   **Particle-vs-Static:** Checks dynamic particles against the static voxel grid map.
    *   **Particle-vs-Player:** Pushes particles out of player bounding boxes (and deals damage).
    *   **Particle-vs-Particle:** Simple sphere-based collision between dynamic voxels.

## Optimizations & Assumptions

*   **Static "Sleeping" System:** Voxels that come to rest on the ground or against static walls are "frozen" and removed from the active simulation list to save performance. They wake up if their support is destroyed.
*   **Early Outs:**
    *   Shape matching skips if the voxel is not deformed.
*   **Approximations:** Collisions use predicted positions for stability, effectively "looking ahead" to prevent tunneling.

## Controls

### Player 1 (Keyboard)
*   **Move:** WASD
*   **Look:** F/H (yaw), T/G (pitch)
*   **Jump:** Space
*   **Shoot:** Left Ctrl (Cost: 1 Matter)
*   **Melee:** Z (Harvests 10 Matter from voxels)
*   **Build:** E (Cost: 10 Matter)
*   **Gravity Tether:** R (Hold to grab, release to throw)

### Player 2 (Keyboard)
*   **Move:** IJKL
*   **Look:** Arrow Keys
*   **Jump:** Right Shift
*   **Shoot:** Right Ctrl (Cost: 1 Matter)
*   **Melee:** M (Harvests 10 Matter from voxels)
*   **Build:** O (Cost: 10 Matter)
*   **Gravity Tether:** P (Hold to grab, release to throw)

### Gamepad
*   **Move:** Left Stick
*   **Look:** Right Stick
*   **Jump:** A / Down Face Button
*   **Shoot:** Right Trigger
*   **Melee:** B / Right Face Button
*   **Build:** X / Left Face Button
*   **Gravity Tether:** Y / Top Face Button

### Global
*   **Reset Game:** R (Main Menu / Game Over)
*   **Settings:** +/- (Active Players), 1-4 (Input Type), [ / ] (Winning Score), D (Drone Intro)
*   **Creative Mode:** C (from Main Menu)

## Creative Mode

Creative mode lets you fly around an empty world and build voxel maps that can be saved and loaded for multiplayer.

### Creative Controls
*   **Move:** WASD
*   **Look:** F/H (yaw), T/G (pitch)
*   **Up/Down:** E / Q
*   **Faster Fly:** Left Shift
*   **Brush Size:** [ / ] (1–20)
*   **Place Voxels:** Left Ctrl (grid-aligned)
*   **Remove Voxels:** Left Alt
*   **Cycle Block Color:** V (next) / B (prev)
*   **Cycle Pickup Type:** Tab
*   **Place Pickup:** P
*   **Remove Nearest Pickup:** Backspace
*   **Activate Structure (Make Dynamic):** K
*   **Save Map Slot:** Ctrl + 1/2/3
*   **Load Map Slot:** Alt + 1/2/3
*   **Exit Creative:** M or Esc

### Creative Controls (Player 2 Keyboard)
*   **Move:** IJKL
*   **Look:** Arrow Keys
*   **Up/Down:** O / U
*   **Brush Size:** Numpad - / +
*   **Cycle Block Color:** Numpad 1/2
*   **Cycle Pickup Type:** Numpad *
*   **Place Voxels:** Right Ctrl
*   **Remove Voxels:** Right Shift
*   **Place Pickup:** Numpad 0
*   **Remove Nearest Pickup:** Numpad .
*   **Activate Structure:** Numpad Enter

### Creative Controls (Gamepad)
*   **Move:** Left Stick
*   **Look:** Right Stick
*   **Up/Down:** RB / LB
*   **Faster Fly:** RT (hold)
*   **Brush Size:** D-pad Left/Right
*   **Cycle Block Color:** D-pad Up/Down
*   **Cycle Pickup Type:** A
*   **Place Voxels:** RT (press)
*   **Remove Voxels:** LT (press)
*   **Place Pickup:** X
*   **Remove Nearest Pickup:** B
*   **Activate Structure:** Y

### Custom Map in Multiplayer
*   **Toggle Custom Map:** U (Main Menu)
*   **Cycle Slot:** L (Main Menu)

## Constructor-Brawler Gameplay Loop

This section outlines the planned "Constructor-Brawler" gameplay loop.

### Phase 1: Matter Core & Melee Harvest
*   **Player Stats:** Health and Ammo are replaced by a single **Matter** resource (0-100).
*   **Exposed State:** When Matter reaches 0, the player becomes **Exposed**.
*   **Melee:** Short-range raycast. 
    *   Hitting a voxel destroys it and grants +10 Matter.
    *   Hitting an enemy applies knockback. If the enemy is Exposed, they are killed.
*   **HUD:** Large Matter bar. Flashes transparent blue on damage to shields (Matter > 0) and transparent red when Exposed.

### Phase 2: Gun & Build Mechanics
*   **Shooting:** Consumes Matter. Prevented if Matter is 0.
    *   Bullets deal damage to enemy Matter.
    *   Bullets deal 0 damage to Exposed players (must use Melee or Physics).
*   **Building:** Consumes 10 Matter to place a static voxel. Edits regenerate through the recycle pipeline but last longer than dynamic debris.

### Phase 3: Physics Finisher (Gravity Tether)
*   **Gravity Tether:** 
    *   **Grab:** Raycast to find a dynamic voxel or pull a chunk off a static structure. Applies spring force towards the player's hand.
    *   **Throw:** Release to apply a massive impulse in the look direction.
*   **Physics Damage:** High-velocity objects hitting an Exposed player will kill them.

## Building

The release layout contains a launcher plus GPU-capable and CPU-only game binaries.
The launcher tries the platform-native GPU backend first, then multithreaded CPU,
then single-threaded CPU. Windows/Linux use OpenGL 4.3 compute. macOS uses Metal
compute alongside raylib's OpenGL 3.3 renderer.

### Windows (MSYS2/UCRT64)

Point `RAYLIB_SOURCE_DIR` at a raylib 5.5 source tree or pass `-RaylibSource`:

```powershell
$env:RAYLIB_SOURCE_DIR = "D:\raylib-master\raylib-master"
.\build.ps1 -Configuration release
```

### Linux

```bash
RAYLIB_SOURCE_DIR=/path/to/raylib ./build.sh release
```

### macOS

Install the Xcode command-line tools and point `RAYLIB_SOURCE_DIR` at raylib 5.5:

```bash
export RAYLIB_SOURCE_DIR=/path/to/raylib
sh ./build-macos.sh release native
open ".build/FPS Game.app"
```

Use `universal` instead of `native` to build arm64 and x86_64 slices together:

```bash
sh ./build-macos.sh release universal
```

The Mac GPU binary links Metal physics with the same CPU MT and CPU ST fallbacks.
The app bundle is ad-hoc signed for local play; distribution outside the local Mac
still requires a Developer ID signature and notarization.

The build scripts place `fps_ray`, `fps_ray_gpu`, `fps_ray_cpu`, and the shader assets in
`.build/bin` (with `.exe` suffixes on Windows).

### Physics backend controls

```bash
fps_ray --physics=auto --physics-report
fps_ray --physics=gpu
fps_ray --physics=gpu-gl43
fps_ray --physics=gpu-metal
fps_ray --physics=cpu-mt
fps_ray --physics=cpu-st
```

`gpu` selects the native GPU implementation (`gpu-gl43` or `gpu-metal`). `auto` is
the normal fallback chain. An explicitly forced GPU mode fails instead of silently
selecting another backend. CPU MT still degrades to CPU ST when no worker can be
created. The active backend and its last-step time are displayed beside the FPS
counter.

The GPU backend keeps compact dynamic particle and voxel state resident between
frames. It uploads static colliders only when the static hash changes, executes every
fixed step due for the rendered frame as one batch, and performs one synchronized
readback for that batch. Use the legacy compact-upload path for A/B comparisons and
print cumulative transfer timings at shutdown with:

```bash
fps_ray --physics=gpu --gpu-transfer-mode=legacy --gpu-transfer-stats
fps_ray --physics=gpu --gpu-transfer-mode=resident --gpu-transfer-stats
```

Smoke runs still create a short-lived window because the GPU backend needs a
graphics context:

```bash
fps_ray --physics=gpu --physics-smoke=60 --physics-smoke-voxels=256
fps_ray --physics=cpu-mt --physics-smoke=60 --physics-smoke-voxels=256
fps_ray --physics=gpu --physics-smoke=80 --physics-smoke-voxels=4096 \
  --physics-smoke-batch=8 --gpu-transfer-stats
```

`FPS_FORCE_GPU_INIT_FAILURE`, `FPS_FORCE_GPU_STEP_FAILURE`, and
`FPS_FORCE_CPU_ST` are available for testing the fallback paths.

### Automated visual physics debugging

Named debug scenarios build isolated voxel structures, run a deterministic fixed-step
simulation, capture diagnostic PNGs, write a JSON report, and exit without requiring
gameplay input. The window is hidden by default but a graphics context is still
created for off-screen rendering.

```bash
fps_ray --physics=gpu --debug-scenario=sleep-wake-floating
fps_ray --physics=cpu-st --debug-scenario=activation-floating --debug-show-window
fps_ray --physics=cpu-mt --debug-scenario=dynamic-freefall --debug-steps=120 \
  --debug-capture-steps=0,1,30,60,120 --debug-output=debug-artifacts/custom-run
```

Available scenarios are:

* `dynamic-freefall`: creates an already-dynamic cluster to isolate basic integration.
* `activation-floating`: activates an unsupported static cluster through the tether path.
* `sleep-wake-floating`: sleeps a supported dynamic cluster, removes its support, then
  activates it through the tether path.
* `tether-throw-floating`: sleeps a supported cluster, removes its support, then releases
  a tether-held voxel toward it through the same proximity-activation path used by gameplay.
* `overhang-impact`: drops an ordinary active voxel onto the unsupported end of a static
  cantilever and checks that the overhang activates and deflects while its root stays supported.
* `large-fall-restoration`: releases a 360-voxel pillar above a static pad, waits for
  the dynamic-to-static restoration attempt, and records partial conversion, post-restore
  speed, spatial spread, and any subsequent reactivation. It automatically runs at least
  650 fixed steps because the failure happens after the structure appears to settle.

Without `--debug-output`, artifacts are written beneath
`.build/bin/debug-artifacts/<scenario>/<active-backend>`. Each backend directory contains
`report.json`, named setup-phase images when applicable, and images for simulation steps
`0,1,5,15,30,60` by default. Debug runs require at least 60 steps so the report can make
a meaningful fall assertion.

Run all backends sequentially through the launcher with:

```bash
fps_ray --debug-matrix --debug-scenario=sleep-wake-floating
```

Matrix mode creates `gpu-gl43` (Windows/Linux) or `gpu-metal` (macOS), `cpu-mt`, and `cpu-st` subdirectories plus a combined
`matrix.json`. An unavailable GPU is reported as `UNAVAILABLE`; assertion, capture, or
execution failures from an available backend fail the matrix. The per-backend reports
include activation and sleep-transition results, centroid motion, particle mass and
simulation membership checks, static glue counts, hash generations, backend fallback
details, timing samples, failures, and capture paths. GPU scenario reports also read back
the first-step SSBOs and record static-hash differences or stale cluster cells, uploaded
`simId` mismatches or missing corners, and zero-mass tagged corners. The same counters
appear in the capture overlay. Overhang reports additionally compare root and tip motion
at step 30, require a retained static root anchor, and assert measurable relative bending.
