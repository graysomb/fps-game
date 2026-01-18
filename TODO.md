# Plan for Integrating Physics Engine & Cleaning Codebase

## Overview
The goal is to fully integrate the new PBD physics engine from `new_physics_engine.c` into `fps_ray.c` while restoring critical gameplay systems (activation/deactivation, recycling) and removing obsolete multiscale voxel logic.

## 1. Cleanup: Remove Multiscale Voxel Logic (`fps_ray.c`)
- [ ] **Remove `span` from `Voxel` struct**: The `int span` field is obsolete. All voxels should be uniform size.
- [ ] **Remove `addVoxelSized`**: Replace usages with `addVoxel` (assuming span=1). however still needs to work with `glue_neighbor_faces_for_voxel`
- [ ] **Remove `emit_multiscale_voxels_from_units`**: Simplify to just emitting single unit voxels.
- [ ] **Remove Multiscale/Span Debug Logging**: Remove all code blocks guarded by `if (v->span > 1)` or logging "span".
- [ ] **Simplify Voxel Initialization**: Update `init_voxel_struct` to remove span calculation and bounding box logic related to variable sizes.

## 1.5. Remove old glue logic
- [ ] remove logic realted to `solve_voxel_glue` and functions that call it
- [ ] replace its usage with `glue_neighbor_faces_for_voxel` in functions related to activation/deactivation, recycling

## 2. Integration: Activation & Deactivation (Sleep/Wake)
- [ ] **Verify `Voxel` struct**: Ensure `fps_ray.c`'s `Voxel` struct has necessary fields for the belief system (`activationBelief`, `freezeBelief`, `sleepFrames`, etc.) consistent with the desired logic.
- [ ] **Port/Fix Belief System**: Ensure `update_static_voxel_belief` and `update_dynamic_activation_beliefs` are correctly implemented and called in the main loop.
- [ ] **Integrate with New Solver**:
    - Ensure the new PBD solver (`simulate_voxel_pbd`) respects the `simulate` flag managed by the activation system.
    - Ensure voxels transitioning from static to dynamic (and vice-versa) properly initialize/reset their particle states (positions, velocities) in the new solver.
- [ ] **Fix `check_voxel_sleeping`**: Ensure the condition for putting a voxel to sleep (velocity checks, support checks) works with the new particle velocities.
- [ ] **Fix `deactivate_sleeping_voxels`** Update to work with new pbd physics engine. ensure that entire clusters of glued voxels are correctly converted to static
- [ ] **Fix `activate_static_voxels_near_dynamic`** ensure it works properly with the new type of glue utilizing shared particles
- [ ] **Fix `batch_glued_dynamic_voxels`** ensure it works properly with the new type of glue utilizing shared particles

## 3. Integration: Recycling Pipeline
- [ ] **Fix `RecycleQueue`**: Ensure the `recycleQueue` is properly managed (push/pop).
- [ ] **Implement `recycle_dead_voxels`**: Ensure voxels that fall out of bounds or are "dead" are recycled or reset correctly.
- [ ] **Restore Debris Logic**: Ensure `STATIC_DEBRIS_OWNER` logic works for temporary static voxels (debris) to be cleaned up or recycled.

## 4. Physics: Static Collision
- [ ] **Implement Static Collision**: The new engine currently lacks collision between dynamic particles and static voxels.
    - port `solve_static_collisions`into to the collision system of the PBD loop at `gather_particle_collisions` and remove depricated solve_static_collisions
    - This should check dynamic particles against the `static_table` (spatial hash) and apply position corrections.

## 5. Final Verification
- [ ] **Verify Build**: Ensure `fps_ray.c` compiles without errors.
- [ ] **Verify Gameplay**:
    - Voxels should fall and stack.
    - Voxels should go to sleep (turn static) when stable.
    - Voxels should wake up (turn dynamic) when disturbed.
    - Static voxels should collide with dynamic ones (no falling through floor).
    - Multiscale artifacts should be gone.
