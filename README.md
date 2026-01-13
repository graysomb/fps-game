# Raylib Split-Screen FPS with Voxel PBD Physics

This project is a split-screen First Person Shooter (FPS) prototype built using [Raylib](https://www.raylib.com/). It features a custom Voxel-based Position Based Dynamics (PBD) physics engine that supports destructible environments, multiscale voxels, and structural integrity.

## Physics Engine Overview

The core of the game is a custom physics engine implemented in `fps_ray.c`. It uses a **Position Based Dynamics (PBD)** approach, which modifies particle positions directly to satisfy constraints, resulting in stable and controllable simulations.

### The Simulation Loop
The physics simulation follows a standard PBD loop:
1.  **Prediction:** Particles move based on their current velocity and external forces (gravity).
2.  **Constraint Solving:** The solver iteratively corrects particle positions to satisfy physical constraints (Shape, Glue, Collision).
3.  **Integration:** Velocities are updated based on the corrected positions.

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

### 2. Glue Constraints
*   **Function:** `solve_voxel_glue`
*   **Description:** Connects adjacent voxels together to form larger structures.
*   **Method:**
    *   Uses **Barycentric Coordinates** to attach a "fine" voxel corner to a specific point on a "coarse" voxel's face.
    *   Includes a **Hinge Limit** to prevent unnatural bending.
    *   Supports **Multiscale Connections** (e.g., attaching a small voxel to a large 2x2x2 voxel).
*   **Breaking:** Bonds break if they are stretched beyond `GLUE_BREAK_STRAIN` or bent beyond `GLUE_BREAK_HINGE_ANGLE_DEG`.

### 3. Collision Constraints
*   **Function:** `solve_particle_collisions`
*   **Description:** Handles interactions between dynamic voxels, the static environment, and players.
*   **Method:**
    *   **Particle-vs-Static:** Checks dynamic particles against the static voxel grid map.
    *   **Particle-vs-Player:** Pushes particles out of player bounding boxes (and deals damage).
    *   **Particle-vs-Particle:** Simple sphere-based collision between dynamic voxels.
*   **Optimizations:**
    *   **Skip Internal:** Collisions between glued voxels are skipped to improve performance.
    *   **Span Awareness:** Collision detection adjusts for the size (span) of the voxel.

## Optimizations & Assumptions

*   **Static "Sleeping" System:** Voxels that come to rest on the ground or against static walls are "frozen" and removed from the active simulation list to save performance. They wake up if their support is destroyed.
*   **Multiscale Voxels:** The engine supports voxels of different sizes (Spans). A span-2 voxel represents a 2x2x2 block but is simulated as a single rigid body, significantly reducing the particle count for large debris.
*   **Early Outs:**
    *   Shape matching skips if the voxel is not deformed.
    *   Glue solving skips if the connected voxels are both static.
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

### Linux
```bash
gcc -fopenmp fps_ray.c -o fps_ray $(pkg-config --cflags --libs raylib) -lGL -lm -lpthread -ldl -lrt -lX11
```
### Windows (MSYS2)
```bash
gcc -fopenmp fps_ray.c -o fps_ray.exe -I/mingw64/include -L/mingw64/lib -lraylib -lopengl32 -lgdi32 -lwinmm
```
