# Optimization Plan: Hybrid CPU/GPU Physics Pipeline

## Overview
We aim to maximize performance by leveraging both the CPU (multithreading) and GPU (Compute Shaders) in parallel.
**Key Insight:** Voxels in this engine **do not share corner particles** (each voxel owns its 8 corners).
- **Consequence:** `solve_voxel_shape` is "embarrassingly parallel." No race conditions between voxels.
- **Goal:** Move heavy math to GPU, keep complex logic on CPU.

## 1. Immediate CPU Multithreading (The "Brain")
Since voxels are independent, we can immediately parallelize the shape matching loop on the CPU without complex synchronization.

- [ ] **Parallelize `solve_voxel_shape`:**
    - Use OpenMP (`#pragma omp parallel for`) or a thread pool for the main loop in `simulate_voxel_pbd`.
    - **Target:** `for (int i = 0; i < voxel_count; ++i) { solve_voxel_shape(&voxels[i]); }`
- [ ] **Parallelize Glue Logic:**
    - `solve_voxel_glue` updates shared particle positions between different voxels.
    - **Strategy:** Use "Graph Coloring" or Atomic operations if moved to GPU. On CPU, use a "Gather-Scatter" approach or simply thread the independent clusters.

## 2. GPU Compute Shaders (The "Muscle")
Move the dense floating-point operations to Compute Shaders to relieve the CPU and bypass memory bottlenecks.

### A. Data Structure Setup (SSBOs)
- [ ] **Particle Buffer:** Linear array of `struct Particle { vec4 pos; vec4 prev_pos; ... }`.
- [ ] **Voxel Buffer:** Linear array of Voxel indices/properties.
- [ ] **Grid/Table Buffer:** Flatten the `table_get` 3D lookup into a 1D SSBO or 3D Texture to accelerate neighbor lookups.

### B. Shader Porting
- [ ] **`solve_voxel_shape.comp`:**
    - **Input:** Voxel SSBO.
    - **Logic:** Compute centroid, principal axes, and target positions for the 8 owned particles.
    - **Output:** Write new `predicted_pos` to Particle SSBO.
    - **Benefit:** Massive parallelism (1 thread per voxel).
- [ ] **`solve_dynamic_collisions.comp`:**
    - **Status:** Partially implemented in `particle_collision.comp`.
    - **Task:** Finalize the "Grid Constraint" mode to replace the slow CPU `table_get` loop.
    - **Logic:** Each particle checks its cell neighbors in the SSBO/Texture grid.
- [ ] **`solve_voxel_glue.comp`:**
    - **Input:** Glue Constraint SSBO (pairs of particle indices).
    - **Logic:** Apply distance constraints.
    - **Safety:** Use `atomicAdd` for position deltas to handle race conditions where multiple constraints pull the same particle.

## 3. The Hybrid Pipeline
- [ ] **Execution Flow:**
    1.  **CPU:** Dispatch Compute Shaders (`Integrate` -> `Shape` -> `Collisions` -> `Glue`).
    2.  **CPU (Parallel):** While GPU crunches physics, CPU performs:
        - AI / Bot Logic.
        - Input processing.
        - **Glue Cluster Analysis:** Check for broken constraints (readback from GPU or calculated lazily) to determine if chunks need to split.
    3.  **Sync:** Barrier wait for Physics Compute to finish before Rendering.

## 4. Specific Refactoring Tasks
- [ ] **`fps_ray.c`**: Isolate `solve_voxel_shape` into a clean function that can be wrapped in a compute dispatch.
- [ ] **`particle_collision.comp`**: Ensure the grid structure matches the CPU's `table_get` logic for seamless transition.

---

# CPU AI Refactor (Constructor-Brawler Update)

The current CPU AI logic in `fps_ray.c` (`UpdateBot`) is designed for a traditional Health/Ammo pickup loop. It needs to be overhauled to support the new Matter-based economy and "Exposed" death mechanics.

## 1. Resource Management (Harvesting)
- [ ] **Deprecate Pickup Seeking:** Remove the loop that searches for `pickups[]` array (Ammo/Health boxes).
- [ ] **Implement "Harvest" State:**
    - **Trigger:** When `bot->matter` is low (e.g., < 25% or < `MATTER_SHOT_COST * 5`).
    - **Search:** Scan for the nearest active **Static Voxel** (Environment).
    - **Pathing:** Move within `MELEE_RANGE` (2.0 units) of the target voxel.
    - **Action:** Synthesize Melee Input (`KEY_Z` equivalent) to destroy the voxel and gain +10 Matter.

## 2. Combat Logic Upgrade
The bot must distinguish between "Shielded" enemies (Matter > 0) and "Exposed" enemies (Matter == 0).

### Phase A: Shield Breaking (Target Matter > 0)
- [ ] **Behavior:** Maintain medium range.
- [ ] **Action:** Use `FireVoxel` (Shoot) to deplete enemy Matter.
- [ ] **Constraint:** Stop shooting if `bot->matter` is critical (save reserve for Panic Builds or Shield buffer).

### Phase B: Finisher (Target `isExposed` == true)
- [ ] **Behavior:** Aggressive closing of distance. Shooting deals 0 damage now, so behavior must change.
- [ ] **Melee Kill:** Move to close range (< 2.0 units) and trigger Melee to knockback/kill.
- [ ] **Physics Kill (Hard Bots):**
    - Identify loose debris or rip a chunk from the wall using Gravity Tether (`KEY_R`).
    - Aim and throw object at the exposed player.

## 3. Self-Preservation (Exposed State)
- [ ] **Trigger:** High priority override when `bot->isExposed` is true.
- [ ] **Flee:** Move directly away from nearest enemy.
- [ ] **Panic Build:** Occasionally trigger `perform_build` while fleeing to place blocks behind (breaking line of sight).
- [ ] **Emergency Harvest:** Aggressively seek any nearby voxel to Melee harvest, restoring Matter to > 0 to regain "Shields".

## 4. Architecture Refactor
- [ ] **State Machine / Utility System:** Refactor the monolithic `UpdateBot` function into weighted utility functions:
    - `CalculateUtility_Harvest()`
    - `CalculateUtility_Combat()`
    - `CalculateUtility_Flee()`
- [ ] **Input Synthesis:** Ensure the bot function can simulate the new key inputs:
    - Melee (`KEY_Z` / `KEY_M`)
    - Build (`KEY_E` / `KEY_O`)
    - Tether (`KEY_R` / `KEY_P`)