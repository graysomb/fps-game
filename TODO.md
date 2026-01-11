# TODO: "Constructor-Brawler" Gameplay Loop Implementation

## Phase 1: Matter Core & Melee Harvest
- [x] **Refactor Player Struct** (fps_ray.c)
    - Remove: `health`, `shield`, `ammo`.
    - Add: `float matter` (0-100), `bool isExposed`, `float matterMax`.
    - Initialize `matter` to 100 on spawn.
- [x] **Implement Melee()**
    - Input: `Z` (Player 0), `M` (Player 1).
    - Logic: Raycast forward (short range, e.g., 2.0 units).
    - Hit Voxel: `DestroyVoxel()`, `matter += 10` (clamp to max). matter edits regenerate through recycle pipline
    - Hit Enemy: Apply Knockback force. If `enemy.isExposed`, trigger `KillPlayer()`.
- [x] **Update HUD**
    - Replace Health/Shield/Ammo bars with a single large "Matter" bar.
    - Visual feedback for "Exposed" state (e.g., bar turns red/cracked).
    - damage to sheilds cuases players screen to flash transparent blue
    - exposed cuases players screen to flash transparent red

## Phase 2: Gun & Build Mechanics
- [x] **Update Shoot()**
    - Input: `Left Ctrl` (Player 0), `Right Ctrl` (Player 1).
    - Cost: `matter -= 1` per shot. Prevent shooting if `matter <= 0`.
    - Damage: Deal `damage` to `enemy.matter`. feedback that player was hit
    - **Logic:** If `enemy.matter <= 0`, set `isExposed = true`. Bullets deal 0 damage to `isExposed` players.
- [x] **Implement Build()**
    - Input: `E` (Player 0), `O` (Player 1).
    - Logic: Raycast to find placement target.
    - Cost: `matter -= 10`.
    - Effect: Place a static voxel 
    - keep the current recycle loop but make players edits last longer

## Phase 3: Physics Finisher (Gravity Tether)
- [x] **Implement Gravity Tether**
    - Input: `R` (Player 0), `P` (Player 1).
    - State: `Holding` vs `Idle`.
    - **Grab:** Raycast to find voxel (pull dynamic voxel/cluster, or activate and pull off chunk of static voxels). Apply spring force to pull towards player's "hand" position.
    - **Throw:** On release, apply massive impulse to object in look direction.
- [x] **Physics Damage Logic** (should largely exist in smush pipline)
    - In `UpdateParticles()` or collision resolution:
    - If Object Velocity > Threshold AND hits Player: 
    - If `Player.isExposed`: `KillPlayer()`.
    - keep respawn timer and screen

## Phase 4: Controls & Cleanup
- [x] **Map Controls**
    - P0: Shoot(L-Ctrl), Melee(Z), Build(E), Tether(R).
    - P1: Shoot(R-Ctrl), Melee(M), Build(O), Tether(P).
    - Gamepad support updates.
- [x] **Balancing**
    - Tune `matter` costs and harvest rates.
    - Tune physics damage thresholds.

## Current Progress
- [x] Plan documented.
