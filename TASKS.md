# Project Tasks

This file tracks the current and future development tasks for the voxel FPS game.

## Gameplay Features (Constructor-Brawler Loop)

### Phase 1: Matter Core & Melee Harvest
- [ ] **Refactor Player Struct** (fps_ray.c)
    - Remove: `health`, `shield`, `ammo`.
    - Add: `float matter` (0-100), `bool isExposed`, `float matterMax`.
    - Initialize `matter` to 100 on spawn.
- [ ] **Implement Melee()**
    - Input: `Z` (Player 0), `M` (Player 1).
    - Logic: Raycast forward (short range, e.g., 2.0 units).
    - Hit Voxel: `DestroyVoxel()`, `matter += 10` (clamp to max). matter edits regenerate through recycle pipline
    - Hit Enemy: Apply Knockback force. If `enemy.isExposed`, trigger `KillPlayer()`.
- [ ] **Update HUD**
    - Replace Health/Shield/Ammo bars with a single large "Matter" bar.
    - Visual feedback for "Exposed" state (e.g., bar turns red/cracked).
    - damage to sheilds cuases players screen to flash transparent blue
    - exposed cuases players screen to flash transparent red

### Phase 2: Gun & Build Mechanics
- [ ] **Update Shoot()**
    - Input: `Left Ctrl` (Player 0), `Right Ctrl` (Player 1).
    - Cost: `matter -= 1` per shot. Prevent shooting if `matter <= 0`.
    - Damage: Deal `damage` to `enemy.matter`. feedback that player was hit
    - **Logic:** If `enemy.matter <= 0`, set `isExposed = true`. Bullets deal 0 damage to `isExposed` players.
- [ ] **Implement Build()**
    - Input: `E` (Player 0), `O` (Player 1).
    - Logic: Raycast to find placement target.
    - Cost: `matter -= 10`.
    - Effect: Place a static voxel 
    - keep the current recycle loop but make players edits last longer

### Phase 3: Physics Finisher (Gravity Tether)
- [ ] **Implement Gravity Tether**
    - Input: `R` (Player 0), `P` (Player 1).
    - State: `Holding` vs `Idle`.
    - **Grab:** Raycast to find voxel (pull dynamic voxel/cluster, or activate and pull off chunk of static voxels). Apply spring force to pull towards player's "hand" position.
    - **Throw:** On release, apply massive impulse to object in look direction.
- [ ] **Physics Damage Logic** (should largely exist in smush pipline)
    - In `UpdateParticles()` or collision resolution:
    - If Object Velocity > Threshold AND hits Player: 
    - If `Player.isExposed`: `KillPlayer()`.
    - keep respawn timer and screen



## Graphics and Rendering

- [ ] Implement a skybox with a proper texture instead of a solid color.
- [ ] Improve voxel texturing and lighting.
- [ ] Add particle effects for explosions, bullet impacts, and destruction.
- [ ] Add post-processing effects (e.g., bloom, motion blur).

## Audio

- [ ] Add sound effects for shooting, jumping, walking, and voxel collisions.
- [ ] Add background music and a simple audio mixer.

## UI/UX

- [ ] Create a main menu with options to start the game and exit.
- [ ] Implement a pause menu.
- [ ] Add support for gamepads.
- [ ] Improve the in-game HUD to be more visually appealing.

## Code and Architecture

- [ ] Refactor the code to improve organization (e.g., separate files for player, voxel physics, rendering).
- [ ] Add more comments to the code to explain complex parts.
- [ ] Set up a build system (e.g., Makefile, CMake) to simplify compilation.

## Networking

- [ ] Implement basic networking for online multiplayer.
- [ ] ENet or make my own?
- [ ] Add a server browser and lobby system.

## Fractal Mode
- [ ] Implement modular fractal world where each node is essentially this code
- [ ] dirty chunking the greedymeshing with modular greedy meshing
- [ ] Implement basic networking for online multiplayer.

## Level building
- [ ] whats the vibe
- [ ] biomes?
- [ ] fractal generation
- [ ] super position

## physics
- [ ] add ray casting to player collisions
- [ ] Implement rigid body physics (floating bodies fall)
- [ ] Implement torque
- [ ] Implement intervoxel sticking

# voxel types
- [ ] goop you move slowly throught
- [ ] bouncy
- [ ] slippy
- [ ] portal