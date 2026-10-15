


This is a vibe coded physics engine I made. Its kinda like Minecraft but fully simulated. I wanted it to be an FPS game but i'm not sure if thats going to work out. Writing the physics was fun though! its currently a jacobi style PBD simulation with [voxel gram-schmidt shape corrections](https://dl.acm.org/doi/10.1145/3677388.3696322) and glue. The voxels have a life cycle where they are greedily meshed (just like minecraft) until physically active, after which they enter the PBD loop.

Currently the game motives interaction with the physics by requiring environment based kills or melee

runs in linux and windows so far

# Warning SLOP


# Raylib Split-Screen FPS with Voxel PBD Physics

## Couch sessions and Firefight

### Combat clarity and grenades

The main menu scans a slowly rotating rat and enemy roster. Left/right browses the
entries; each explains its ability and counter. Scans advance automatically every ten
seconds without spawning gameplay actors.

Throw a grenade with **C** (keyboard P1), **period** (keyboard P2), or **D-pad Up**.
Carry two, start with one, and receive at least one at Firefight intermissions.
Reviving does not grant ammunition. Grenades detonate half a second after their first bounce, but never before
1.5 seconds after release. Later bounces do not reset the countdown; grenades still
in flight wait for their first impact. Orange supply stations add one grenade and refill after 30 seconds. Grenades
work in offline combat and LAN PvP; creative, paused, and downed players cannot throw.

**Explosions hurt everyone**, including the thrower, teammates, and other enemies.
Sapper charges and grenades share 60 damage falling to zero at four units, cover checks
before terrain destruction, a two-unit demolition radius, and a four-unit debris impulse.
Blast-driven debris retains friendly fire for ten seconds or until it sleeps.
Enemies abandon work and retreat from nearby live bombs; the Walking Foundry holds its
position. Sapper charges retain their four-second fuse and remain shootable to disarm.

Melee removes one of a Shieldbearer's four vertical shield sections per swing. The
shield has translucent blue panels and opaque outlines. Lobbers compute a damped arc
from the held structure's center, lock their target during wind-up, and reposition if
cover blocks the arc. Swept collision checks ignore the projectile's own glued group.
Standard 60 Hz physics damping is preserved.

Deaths produce a brief flash and collapsing angular fragments, distance-based sound,
and kill confirmation. Enemy activity uses a 32-voice spatial sound pool; split-screen
plays each cue once for the nearest living local listener. Swarmers have a bright top
marking. When one to three enemies remain without reinforcements, radar hints start
after eight seconds and refresh last known positions every four seconds.
Debug particle spheres stay off unless explicitly enabled.

Foundry supplies occupy health stations at (0, +/-10), Matter at (+/-12, 0), charged
shots at (-14, +/-20), gold tethers at (14, +/-20), and grenades at (0, +/-28).
Crosscut retains four outer Matter stations and central gold tether, with four grenade
stations based on (18, 6) and four charged-shot stations based on (9, 6), rotated by
quarter turns. Supported clear positions are chosen within two units; supplies drop
when support is destroyed. Custom layouts and weapon loadouts remain unchanged.

In Enemy Lab, **F8 / D-pad Down** cycles grenade range, bomb avoidance, death feedback,
close shield melee, and elevated Lobber fixtures. F4 resets the selected fixture.
`tools/validate_firefight.ps1` (and `-Gpu`) checks real structure flights, damage,
supplies, inventory, events, and two-process LAN late join with dropped snapshots and
duplicate inputs, alongside the existing Firefight regressions.
LAN peers require protocol 8; use the same build on every machine.

### Specialist enemies and Enemy Lab

Specialists use destructible equipment and telegraphed abilities rather than the
normal enemy's automatic fire. Goliaths harvest material to restore armor;
Shieldbearers turn their shields toward threats; Sappers plant disarmable
four-second demolition charges that launch surrounding debris; Lobbers seek and throw
real 8–24-voxel groups using gold tethers; Builders construct
harvestable barricades and ramps; Anchors brace nearby allies against knockback.
Break two of the Walking Foundry's four orange support assemblies to open its core
permanently. Its HUD health bar can be depleted with bullets, melee, or thrown debris;
it continues attacking with group throws, a telegraphed stomp, and bounded reinforcements.
Constructors search and approach clear building sites. Scrapper encounters now use Goliaths;
the legacy `scrapper` debug name aliases `goliath`.
Lobbers never generate replacement ammunition: control the arena's scrap supplies.

Open **Enemy Lab** from the couch menu or advanced setup (`L`). It has isolated
fixtures and a mixed encounter, a fixed seed, and invulnerability enabled by
default. `F3` / LB selects the next case; `F4` / top face resets it; `F5` / left
face pauses; `F6` / RB advances one simulation step; `F7` / Back toggles
invulnerability; Escape leaves. The overlay shows the ability state, target,
timers, equipment count, and intended counterplay. Normal movement, shooting,
melee, and tether controls remain available.

Launch directly with `--enemy-lab=sapper`, `--enemy-lab=mixed`, or
`--enemy-lab=all` (cycles isolated cases every 30 simulation seconds). Valid
enemy names are `goliath`, `shieldbearer`, `sapper`, `lobber`, `builder`,
`anchor`, and `foundry`. The offline Firefight console also accepts `spawn NAME`
and uses reserved enemy slots rather than replacing human seats.

`--enemy-smoke=all` runs deterministic behavior checks and a mixed physics/render
stress scenario, returns nonzero on failure, writes enemy previews under
`artifacts/`, and reports four-view timing. A single name selects its behavior
checks before the mixed scenario. Run from the repository root so shaders are
available; select `--physics=cpu-mt` or `--physics=gpu` explicitly. The normal
`tools/validate_firefight.ps1` build workflow runs these checks too, plus autonomous
Constructor/Lobber scenarios inside the full arena, controller hot-plug tests,
power-up replenishment, boss lifecycle, and demolition impulse checks.

Specialists and the lab are offline-only. Ordinary enemies retain player-count
scaling; specialist introductions stay at one, with one specialist active in
waves 5–10 and at most two in wave 11. Endless repeats the twelve-wave cycle.

Firefight ramps from swarmers to normals to Goliaths. Solo opening wave totals are
3 swarmers; 6 swarmers; 6 swarmers plus 1 normal; then 6 swarmers plus 2 normals and
1 Goliath. Wave 5 repeats a Goliath patrol. Waves 6–10 introduce Shieldbearer, Sapper,
Lobber, Builder, then Anchor alongside six swarmers and a normal. Wave 11 mixes
those five specialists with nine swarmers, two normals, and two Goliaths. Wave 12 pits
the team against the Walking Foundry and at most six swarmer reinforcements.
The first three waves
cap simultaneous solo enemies at 2, 3, and 4 respectively, and the opening two
waves space reinforcements four seconds apart. Totals scale with player count.

Local Versus defaults to **Crosscut**, a symmetric quarry for 2–4 players.
Crosscut occupies a 64×64-unit arena; Foundry Holdout occupies 80×80 units, with
reinforced factory bays, roofed side routes, and staggered full-height cover. Other worlds retain their
original bounds. LAN peers must use the same build (protocol 8 includes human collision bounds).
Use left/right in the lobby to cycle arenas; Tab or the top face button toggles
the custom map slot selected in advanced setup. Crosscut is also world 9 there.
Four raised terraces surround short bridges on orange breakable supports.
The ground circuit remains available after the crossings collapse. Contest the
central gold tether or collect Matter at four equally spaced outer supplies.
Screened perimeter respawns favor cover and distance from living opponents.
Rematch restores the selected arena and its destructible geometry.

The offline menu offers Firefight, Local Versus, Learn Combat, controller settings,
and the existing advanced world/LAN tools. Local menus automatically select a connected
controller (D-pad/left stick, bottom face to confirm, right face to go back), or use
keyboard arrows/Enter/Escape when none is connected. Menu hot-plugging does not change
existing gameplay seats; entering a new lobby assigns its initiating device to player one.
Additional controllers join the lobby with Start. `K` adds a second keyboard seat;
keyboard seats ready with Enter and Right Shift. Up to four humans can play
alongside twelve active firefight enemies. Controller disconnects pause offline
play; reconnect or press Start on an unused controller to explicitly reassign it.

Firefight defaults to twelve authored waves on **Foundry Holdout**, with a victory
screen and rematch. Left/right in the lobby selects Endless Survival; in Endless,
Tab or the top face button switches between Foundry and the current advanced map.
The arena has a raised holdout, two stair approaches, a lower circulation route,
and two heavy slabs on orange breakable supports. Match destruction persists
between waves; rematching rebuilds the arena. Opposing flank routes offer charged-shot
and gold-tether power-ups with 30-second respawns and intermission replenishment.
Eight real scrap bundles restock between waves only where their supply sites are clear.
Ordinary walls/roofs are thicker and supports are redundant; marked demolition points
remain vulnerable without changing global material strength.

In cooperative Firefight, lethal damage downs a player for 15 seconds. Hold LB
(keyboard P1 `X`, P2 comma) within two units and clear sight for three seconds to
revive them. Completing a revive costs 25 Matter and restores 35 Matter plus two
seconds of invulnerability. You must have more than 25 Matter to revive. Taking
damage, moving out of range, or releasing cancels progress without charging.
Bleeding out spends one of three shared respawn reserves; exhausted players wait
for the next wave. Solo players use reserves directly. Clearing a wave restores
the team, and every second cleared wave awards a reserve.

Preparation lasts 20 seconds. Vote for cheaper building (5 Matter), improved
harvesting (15 instead of 10 Matter), or 25% more tether throw impulse. The team
gets one bonus for the next wave; ties default to cheaper building. Everyone can
ready to skip the remaining break after five seconds. Top face (keyboard `V`/`N`)
toggles the vote panel so the normal combat bindings remain available for building
and harvesting. The vote panel uses D-pad left/right and bottom face to ready;
keyboard players use their look keys and jump.

Learn Combat teaches harvesting, exposure, melee finishing, and debris finishing
through completed actions. Leave through the pause menu at any time. Per-seat
sensitivity, deadzone, and invert-Y settings are saved in `couch_settings.cfg`.
The new Firefight modes are local-only; existing LAN PvP and creative protocols
remain four-player modes.

For a quick check without raylib or a display, run
`python tools/validate_headless.py` with Python 3.9+ and GCC on your PATH.
This compiles and runs the Firefight rules and network protocol tests, then
checks the shared GPU shader layouts. Use `--cc clang` to select another C
compiler, or `--sanitize` where AddressSanitizer and UndefinedBehaviorSanitizer
are available. CI runs these checks with sanitizers on Linux for pushes and
pull requests.

After building raylib with `build.ps1`, run `tools/validate_firefight.ps1` and
`tools/validate_firefight.ps1 -Gpu` to build and check CPU/OpenGL variants. These
also produce `.build/couch/fps_ray_cpu.exe` and
`.build/couch/fps_ray_gpu.exe`, alongside the matching shaders and runtime. Run
these from `.build/couch` (or run the normal game from the repository root).
Tests cover wave scaling, votes, combatant
capacity, revives/reserves, solo play, tutorial progression, rendering layouts,
and returning to PvP. `tools/validate_gpu_layout.py` checks shared shader layouts;
Metal execution and real-controller couch play still require platform/device QA.

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
by glue. Particle and fluid state stays resident across both substeps and unchanged
frames; static voxels are uploaded only when the static-grid generation changes.
The dynamic hash is sized to the active population instead of clearing the maximum
262,144 buckets. Each substep dispatches prediction, dynamic hash construction,
scene and particle collisions, Jacobi apply passes, VGS iterations, static collisions,
glue-break evaluation, GPU particle splitting, topology rebuilding, wake propagation,
and velocity finalization. PBF builds one 64-entry neighbor cache per particle and
reuses it for density, position correction, and viscosity, with a full hash scan only
for dense particles that overflow the cache. Fluid-only scenes omit the solid pair,
VGS, and glue-break passes. The compact result is committed back to the canonical CPU
arrays once per rendered-frame batch. Topology mutations such as actually breaking
glue or converting a sleeping cluster back to static remain ordered CPU commits.

The CPU PBD loop uses a lightweight, persistent thread pool to avoid per-frame
thread creation.
Key ideas:
*   **Parallel-for dispatcher:** A small worker pool runs chunked ranges over arrays (particles/voxels).
*   **Jacobi-style accumulators:** Constraints write into per-particle correction accumulators (`corr_sum`, `corr_weight`) using atomic float adds. A later “apply” pass updates positions, preserving determinism within a single iteration.
*   **Stable particle snapshots:** Each constraint iteration snapshots `active_particles` into a local array before building the spatial hash and gathering collisions, avoiding concurrent mutation while workers read.
*   **Thread-safe spatial hash:** Only buckets touched by the previous CPU step are cleared; inserts use atomic exchange for lock-free binning.
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
*   **Jump:** B / Right Face Button
*   **Shoot:** Right Trigger
*   **Melee:** A / Down Face Button
*   **Build:** RB / Right Bumper
*   **Gravity Tether:** LT / Left Trigger

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

## Bullet and tether gun

The shared bullet/tether gun is generated from beveled low-poly meshes and uses
the GLSL 330 shaders in `shaders/weapon*`. It appears in first-person views and on
standard players/bots; claw opening and cyan energy intensify while tethering.
Successful shots drive muzzle flash and recoil. Gold tether uses amber energy.
While Gold Tether is charged or held, two pulsing gold ribbons orbit the gun through
its vertex/fragment shaders. They follow recoil and melee on either gun, remain
visible with bloom disabled, and disappear when the charge is spent. The charged
state also travels in the existing LAN player-state flags byte.
Dynamic Shot immediately equips a larger charcoal-and-silver launcher with a
hollow octagonal barrel, transparent chamber and three glowing yellow cubes that
drift and rotate inside it. These cubes are visual particles; the power-up still
fires the existing gold physics voxels. The launcher also supports tethering and
grip-butt melee, and returns to the standard gun when the power-up clears. Its
appearance uses the existing LAN power-up flag without a protocol change.
Projectile collision, damage, and tether forces retain
their existing behavior.

Weapon bloom is isolated per viewport and respects world occlusion. Pass
`--no-weapon-bloom` to disable its extra render passes. Melee rotates the gun around
the grip and strikes with the grip heel. Attacks, mining and misses use the
same 0.32-second swing; the heel follows the existing hit sweep during its active
window. Rigid pose transforms preserve the flat colors, edge lines and glow, and
the tether follows the animated muzzle. Hit detection, damage, mining and cooldown
are unchanged. Goliaths retain their arm attack. The gun is hidden in creative mode
and respawn. LAN peers must all use protocol
version 5, which adds successful-shot sequence/age and powered-shot visual state.

After building raylib with `build.ps1`, run the Windows visual checks from the
repository root:

```powershell
.\tools\validate_weapon.ps1
```

This compiles and runs OpenGL 3.3/4.3 previews, shot-state and network serialization
tests, melee pose/hit checks, and the actual gameplay renderer with one through four
viewports, swing phases and remote-player swings. PNGs are
saved beneath `artifacts/weapon`; logs and executables are in `.build/weapon`.
The preview includes wall occlusion, extreme pitch, gold tether, firing, resize,
and a rendering-only four-view bloom comparison. Launcher checks cover transparent
glass, cube containment over two minutes, pickup/reversion and replicated state.
Gameplay captures also include launcher floating, firing, tethering and melee
animation frames. Timing excludes gameplay physics.
On other platforms, compile `tests/render_weapon_preview.c` with
`weapon_renderer.c` and raylib, then run it from the repository root with an
existing output directory as its first argument (optional second: `--no-bloom`).

## LAN multiplayer

LAN play uses a host-authoritative ENet session with up to four players total;
each machine may contribute one or more local split-screen players. The host
alone runs voxel physics; clients receive the
static world and compact dynamic-voxel transforms. Player input and snapshots run
at 60 Hz. Clients predict only their own movement/look, then reconcile and replay
unacknowledged inputs. Shooting, melee, building, tethering, creative edits,
pickups, scores, and voxel physics remain authoritative on the host.

Choose the number of LAN-local players with `-`/`+` on the main menu, then press
`H` to create a lobby on port 27015 or `J` to join a host on the same machine. Clients
can toggle ready with `R`; the host starts the match
with `Enter` after everyone has joined. Direct-IP and alternate-port startup is
available from the command line:

```text
fps_ray --lan-host
fps_ray --lan-host=27125
fps_ray --lan-connect=192.168.1.25:27125
fps_ray --lan-host --lan-local-players=2
fps_ray --lan-connect=192.168.1.25:27125 --lan-local-players=2
fps_ray --lan-host --lan-creative
fps_ray --lan-connect=192.168.1.25:27015 --lan-creative
```

Add `--net-stats` to print per-second tick, peer, RTT, and byte counters. The wire
format uses explicit endian-safe field encoding, so Windows, Linux, Intel macOS,
and Apple Silicon macOS can share a session. ENet 1.3.18 is vendored under
`third_party/enet` and built by all three platform scripts.

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

### Position-based fluids

Creative mode can place runtime PBF cells that reuse the voxel corner-particle
storage. Press `C` for player one, keypad `3` for player two, or the left stick
button on a controller to switch the creative brush between `SOLID` and `FLUID`.
Fluid cells are intentionally excluded from saved map occupancy, glue, solid-voxel
sleep, activation, projectile, and recycle paths. They collide with immutable static
geometry without waking it and exchange collision corrections with solids that
are already active. A settled fluid-only scene runs at one quarter rate after a
conservative stability window and immediately returns to full rate when a player
approaches or any solid physics particle becomes active.

The implementation follows the density constraint, artificial-pressure surface
correction, and XSPH viscosity stages from Macklin and Mueller's *Position Based
Fluids*. A 0.5 m fluid cell contains eight independent particles on a 0.25 m
lattice. Each physics substep runs two PBF Jacobi iterations with rest density
1000 kg/m3, `k=0.2`, `n=4`, `dq=0.2h`, and XSPH viscosity 0.15. The same stages
run through the normal GPU -> CPU MT -> CPU ST fallback chain.

Press `F5` to switch between low-poly particle spheres and a faceted approximate
hull. Surface mode bins particles into 0.25 m cells, builds a coarse scalar density
field, averages its edge crossings into one surface-net vertex per boundary cell,
and connects those vertices with flat-shaded triangles. With GPU physics the four
surface-net passes consume the resident particle SSBO and write renderable vertex
buffers directly; CPU backends retain the cached CPU mesh builder. The startup mode
can also be selected with:

```bash
fps_ray --fluid-render=particles
fps_ray --fluid-render=surface
```
Smoke runs still create a short-lived window because the GPU backend needs a
graphics context:

```bash
fps_ray --physics=gpu --physics-smoke=60 --physics-smoke-voxels=256
fps_ray --physics=cpu-mt --physics-smoke=60 --physics-smoke-voxels=256
fps_ray --physics=gpu --physics-smoke=80 --physics-smoke-voxels=4096 \
  --physics-smoke-batch=8 --gpu-transfer-stats
```

Runtime GPU failures immediately continue on CPU, then retry the GPU after 5,
15, and 60 seconds. A successful retry repacks state from the CPU and must pass
one GPU step before fallback status is cleared. Validation failures receive one
retry; transient dispatch/upload failures receive three. Permanent startup
failures are not retried.

`FPS_FORCE_GPU_INIT_FAILURE`, `FPS_FORCE_GPU_STEP_FAILURE`,
`FPS_FORCE_GPU_STEP_FAILURE_ONCE`, `FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE`,
`FPS_GPU_RETRY_IMMEDIATE`, and `FPS_FORCE_CPU_ST` are available for testing the
fallback and recovery paths.

Collision filtering can be measured independently on every backend. Set
`FPS_COLLISION_DYNAMIC_LEGACY=1` to put every simulated particle back into the
dynamic pair/static-collision passes, or `FPS_COLLISION_STATIC_LEGACY=1` to use
per-voxel static AABBs instead of the shared greedy surface cache. Setting both
provides the pre-filter baseline. Debug JSON samples report simulated versus
collision particles and static occupied-cell, surface-cell, patch, and patch-reference
counts.

High-speed voxels released by the tether use a short-lived swept-particle CCD
pass against the static and active voxel grids. It is enabled by default and shared by the
GPU, CPU MT, and CPU ST paths. Use `--no-tether-throw-ccd` for an endpoint-only
A/B run, or `--tether-throw-ccd` to enable it explicitly. Tracked tether impactors
remain outside the structural glue graph when their target wakes and rebound with
15% normal restitution. The swept response queues the projectile's equal-and-opposite
impulse for a static hit cell, then distributes it by mass across the activated wall
island as soon as that island has dynamic particles.

### Automated visual physics debugging

Named debug scenarios build isolated voxel structures, run a deterministic fixed-step
simulation, capture diagnostic PNGs, write a JSON report, and exit without requiring
gameplay input. The window is hidden by default but a graphics context is still
created for off-screen rendering.

When a calm dynamic island reaches its sleep cutoff, each voxel center is compared
with its original grid cell. Islands entirely within one quarter voxel and whose
original cells remain available are restored to static voxels; displaced or blocked
islands retain their exact particle pose and sleep in place. The normal equilibrium
strain and shear limits still prevent significantly deformed islands from snapping.

```bash
fps_ray --physics=gpu --debug-scenario=sleep-wake-floating
fps_ray --physics=cpu-st --debug-scenario=activation-floating --debug-show-window
fps_ray --physics=cpu-mt --debug-scenario=dynamic-freefall --debug-steps=120 \
  --debug-capture-steps=0,1,30,60,120 --debug-output=debug-artifacts/custom-run
fps_ray --physics=gpu --debug-scenario=falling-pillar-sleep --debug-steps=1250 \
  --debug-gif --debug-gif-fps=20 --debug-output=debug-artifacts/pillar-video
```

Available scenarios are:

* `dynamic-freefall`: creates an already-dynamic cluster to isolate basic integration.
* `creative-reset-cache`: populates the shared static surface cache, enters the empty
  creative world, and verifies that no render or collision patches survive the reset.
* `activation-floating`: activates an unsupported static cluster through the tether path.
* `sleep-wake-floating`: sleeps a supported dynamic cluster, removes its support, then
  activates it through the tether path.
* `sleep-snap-original`: activates a grounded 4×4×6 static pillar, settles it within
  the near-origin tolerance, and verifies that all 96 voxels snap back to their original
  static cells.
* `sleep-snap-no-rebuild`: continues the tether-activated pillar test beyond the delayed
  restoration deadline and verifies that its already-correct static cells are not removed
  and rebuilt from obsolete activation snapshots.
* `sleep-recycle`: ages a sleeping dynamic island into the delayed restoration queue and
  verifies that the complete island is removed atomically.
* `falling-pillar-sleep`: tips a 90-voxel pillar, lets it bend and settle, then runs past
  the legacy 10-second dynamic lifetime boundary to catch pose-snapping conversions.
* `tether-throw-floating`: sleeps a supported cluster, removes its support, then releases
  a tether-held voxel toward it through the same proximity-activation path used by gameplay.
* `overhang-impact`: drops an ordinary active voxel onto the unsupported end of a static
  cantilever and checks that the overhang activates and deflects while its root stays supported.
* `pillar-impact-drift`: drops a 500-voxel undercut pillar and bounds its late sliding speed
  and displacement.
* `single-corner-hinge`: verifies that a one-corner glue remnant breaks at the configured
  hinge-angle limit.
* `tether-active-then-static`: tethers an active voxel and then a static cluster voxel,
  verifying that activation/array compaction does not lose or redirect the tether target.
* `tether-active-glued-static`: tethers an active voxel glued to a static anchor and verifies
  whether stress/strain fracture checks trigger under tension to break the static glue bond.
* `tether-thin-wall-ccd`: throws a tether voxel fast enough to cross a one-voxel-thick
  wall between endpoint samples and verifies that swept particle rays stop it, rebound
  it, transfer forward momentum to the wall, and do not join their constraint graphs.
* `cube-stress`: stress test with scalable dynamic cube (size set by `FPS_CUBE_SIZE`).
* `pbf-container`: drops 512 fluid particles into a tall immutable voxel container, captures
  both render modes, and checks particle membership, density, settling, containment, and that
  fluid never activates the container.
* `pbf-dynamic-container`: drops the same fluid volume onto an elevated active PBD floor and
  checks that the deforming solid supports it on GPU, multi-threaded CPU, and single-threaded CPU.
* `pbf-scaling`: benchmarks a configurable cubic PBF volume with convex-hull rendering on every
  step. Use `--debug-pbf-side=N`; the particle count is `8 * N^3`. Reports include average
  physics, hull-render, and total frame times plus derived unthrottled FPS after warm-up.
  GPU reports also include timer-query compute/hull time, packing/upload/readback wall time,
  bytes transferred per frame, dispatch count, hull backend, and maximum neighbor overflow.
* `pbf-dense-neighbors`: compresses 512 particles into one smoothing neighborhood and verifies
  that GPU cache overflow is detected and routed through the complete hash-scan fallback.

For example, benchmark 4,096 particles across all three backends with:

```bash
fps_ray --debug-matrix --debug-scenario=pbf-scaling --debug-pbf-side=8 \
  --debug-steps=120 --debug-capture-steps=0
```

Without `--debug-output`, artifacts are written beneath
`.build/bin/debug-artifacts/<scenario>/<active-backend>`. Each backend directory contains
`report.json`, named setup-phase images when applicable, and images for simulation steps
`0,1,5,15,30,60` by default. Debug runs require at least 60 steps so the report can make
a meaningful fall assertion.

Pass `--debug-gif` to additionally stream a 640x360 looping `simulation.gif` into the
backend output directory. `--debug-gif-fps=N` selects 1-60 playback frames per second
(20 by default) and automatically enables GIF recording. The encoder is built into the
harness, so it does not require FFmpeg, Python, or ImageMagick.

Run all backends sequentially through the launcher with:

```bash
fps_ray --debug-matrix --debug-scenario=sleep-wake-floating
```

Matrix mode creates `gpu-gl43` (Windows/Linux) or `gpu-metal` (macOS), `cpu-mt`, and `cpu-st` subdirectories plus a combined
`matrix.json`. An unavailable GPU is reported as `UNAVAILABLE`; assertion, capture, or
execution failures from an available backend fail the matrix. The per-backend reports
include activation and sleep-transition results, centroid motion, particle mass and
simulation membership checks, static glue counts, hash generations, backend fallback
details, collision-shell and static-surface counts, timing samples, failures, and capture
paths. GPU scenario reports also read back
the first-step SSBOs and record static-hash differences or stale cluster cells, uploaded
`simId` mismatches or missing corners, and zero-mass tagged corners. The same counters
appear in the capture overlay. Overhang reports additionally compare root and tip motion
at step 30, require a retained static root anchor, and assert measurable relative bending.

### Rat player characters

Human-controlled characters use procedural low-poly rats with small black eyes, a transparent one-eye scouter,
and player-colored gear. Bots retain their original appearances. Rat bodies, status flashes,
walking/airborne/downed poses, and first-person weapon paws share the same renderer.
Human feet sit one unit below the unchanged eye anchor. Their collision and damage box
is 0.78 units wide and 1.12 units tall; bot bounds are unchanged. Walking follows distance
traveled with planted feet, short terrain adjustments, and a landing blend. The scouter
rim and gear share weapon bloom and remain luminous with `--no-weapon-bloom`.
LAN peers must use protocol 10, which also synchronizes chunk regeneration previews.

Run `tools/validate_character.ps1` on Windows after building raylib. It checks both
OpenGL versions, shader fallback, mesh data, emission occlusion, split-screen play,
ground contact, swept walls/ceilings, platform landings, and a loopback host/client snapshot.
Reference views and weapon-paw captures are
written to `artifacts/rat`; logs and preview executables are in `.build/rat`.

### Combat feedback tuning

Shield sections absorb six ordinary projectile hits; one melee strike breaks a section.
Anchors restore 8 Matter/second to up to three visible allies within eight units, without
stacking. Exposed allies get a two-second finishing window before protection returns.
Healing never rebuilds lost armor or shield sections.

Enemy deaths leave one neutral, tetherable body voxel. Up to 32 remain for 45 seconds;
held bodies are protected from cleanup. Footsteps use quiet, grounded movement-driven
noise taps. Grenade and Sapper explosions show a brief four-unit blast envelope and
light smoke; cover still blocks damage. Tether throws converge on the reticle target
without carrying residual sideways velocity, while gravity still affects their flight.

Enemy Lab F8 / D-pad Down cycles through eight fixtures, including Anchor recovery,
footstep audition, and an offset tether throw. Regression runs export blast previews and
three footstep WAV samples under `artifacts/`.


### Chunk regeneration and destruction profiling

Terrain retains its existing restoration deadlines. Eligible missing cells are grouped
by the 32-cell (16-metre) rendering chunk, previewed with a 0.6-second reconstruction
scan, then restored together. Newer damage waits for its own deadline. Intact terrain
and construction are preserved. Actors and dynamic debris can delay an overlapping
chunk without preventing free chunks from regenerating. Collision, support and mesh
updates are coalesced per batch; the preview is non-solid.

On Windows, after building raylib, run:

```powershell
.\tools\profile_destruction.ps1 -Map foundry -Backend gpu
.\tools\profile_destruction.ps1 -Map versus -Backend gpu
.\tools\profile_destruction.ps1 -ActivateStructure -Frames 1800
.\tools\profile_destruction.ps1 -ActivateStructure -Activation tether -Frames 360
.\tools\profile_destruction.ps1 -ActivateStructure -Backend cpu-st -Frames 360 -StateChecksums
```

The hidden-window harness disables AI and injects 24 grenades by default. Reports in
`artifacts/profile-*` include phase averages, individual regeneration commits and
frame-time percentiles. `-RestorationOnly` isolates restoration and rendering without
physics; its numbers are not whole-game frame times. The full workload separates the
quiet tail (frames 2401 onward) from active destruction. Large active-physics costs
and occasional global collision/support refresh spikes can still remain.

Activation removes static voxels in a support-invalidation batch, then retains the
existing full support refresh. Tether activation invalidates only affected mesh chunks.
Activation reports include removal/support timings, column-scan counts and individual
activation-window frames. `-StateChecksums` compares captured voxel, particle, support,
glue and restoration state outside measured intervals; use `cpu-st` for deterministic
follow-on physics comparisons. Profiling counters compile only into the test harness.
`tests/activation_batch_test.c` checks nesting, empty/failed operations, index compaction
and preserved restoration deadlines; the mesh regression checks tether chunk seams.

Static support now caches face connections by voxel identity and recomputes beliefs
for affected connected components. Collision patches and GPU collider ranges persist
per chunk; normal edits update affected cells and buffer ranges before physics runs.
Interior contacts resolve before edge contacts so chunk seams do not add sideways
impulses. Cache allocations are reused and cleared on world/backend reset.

`cache-frame` lines report support nodes, neighbor rebuilds, collision chunks, patch
counts, static upload bytes, and cache timing. GPU static timing includes CPU cache
work when it is first requested there; do not add these overlapping timings.
Set `FPS_SUPPORT_FULL_REBUILD=1` or `FPS_COLLISION_CACHE_FULL=1` to use the original
full support or surface-cache builders for comparison. Unexpected bulk edits and
incremental-cache capacity failures use synchronous full rebuilds.
`tests/static_cache_test.c` compares support against the full calculation, checks
exposed-face coverage, probes chunk seams and hash deletion, verifies GPU readbacks,
and exercises allocation reuse, restoration, and full-cache recovery. It runs in
both CPU and GPU variants of the validation script.

`tools/validate_firefight.ps1` includes chunk geometry and regeneration regressions.
`tests/chunk_regeneration_test.c` covers deadlines, boundaries, occupancy, cancellation,
capacity, preview depth and multipart network messages. Set `FPS_REGEN_SHADER_FALLBACK=1`
to exercise the translucent fallback. `tests/combat_clarity_lan.c` additionally checks
preview resynchronization and the authoritative terrain commit over loopback ENet.

### Cellular voxel surfaces

Solid voxel faces share four GPU-simulated second-order binary patterns, updated
at 20 Hz in a 256x256 atlas. The rule uses each cell's previous state as well as
its current neighborhood, with short fading activity trails. It follows the default
rule in https://xcont.com/cell/single_2nd.js. Each cube gets an atlas offset and
rotation from its original grid position; remote dynamic proxies use their stable
network identity. Local face coordinates keep the pattern attached when
voxels activate or rotate; static greedy meshes are not rebuilt for animation.
Material colors and black outlines remain intact, and fine patterns fade with
screen-space distance. Patterns are cosmetic and not synchronized over LAN.
Use `--no-voxel-ca` for the original plain surfaces. Shader/resource failures fall
back automatically; `FPS_VOXEL_CA_FALLBACK=1` exercises that fallback explicitly.
`tests/voxel_ca_test.c` compares GPU steps against a CPU reference, checks update
limits, sustained activity and resource reinitialization, exports eight seconds of
animation frames, and measures one/four-view rendering. The chunk-mesh regression
checks static/dynamic visual continuity and rotated sleeping cubes with the
simulation held fixed.
