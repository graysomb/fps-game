# Bugs TODO
- [x] Smush always shows for any damage, should only show for kills by smush
- [x] need a slight indication in the coloring if a static voxel cannot be activated
- [x] smush should give 4x points not just 2x, total points should appear at the top right of players hud in points/total points in gold. They should glow/jiggle when they are updated
- [x] Multi span voxels don't do multi span damage, multispan voxels should do n*voxel damage to players
- [x] The shield recharge delay is based on the first box of hot hit, not any voxel hit. every voxel hit should cause shield recharge delay to start over.
- [x] The walls and floor need some kind of generated texture so that you can tell how close you are.
- [x] higher k/d should result in higher bullet damage for players.
- [x] lower k/d should result in larger bullets with bigger hit boxes.

- [] lower k/d results in bigger boxes but it cuases the player to hurt themselves. prevent type 0,1,2 bullets from hitting oneself.


# Feature TODO
- [] add a grappling hook. shoots out from players as a gray bullet voxel attached by a grey line to the player, pulls the player to bullet collision location. IT does not activated voxels and works on static, dynamic, and dynamic multiscale. it lets them stick there for 5 seconds, on keyboard make it e and o, on game pad make it a bumper

# Physics Performance TODO

## Phase 0: Rendering Optimization (Critical Bottleneck)
- [x] Implement `DrawMeshInstanced` for dynamic voxels to replace immediate-mode `drawCubeMan`.
  - [x] Create a static `Mesh` for a unit cube (reuse `GenMeshCube`).
  - [x] Implement a buffer management system to collect `Matrix` transforms and `Color` data for all dynamic voxels each frame.
  - [x] Batch draw calls by material/color to minimize state changes, or use a shader that accepts per-instance color.

## Phase 1: Spatial Hash Architecture (The Foundation)
- [x] Split `table` into `static_table` and `dynamic_table`.
- [x] Modify `physics_step` to only rebuild `dynamic_table`.
- [x] Create `init_static_hash` to populate `static_table` once.
- [x] Implement incremental updates: Call `static_table_insert` / `static_table_remove` immediately when static voxels are created or destroyed, avoiding full rebuilds.
- [x] Update `table_get` to check both tables (or specific one via flag).
- [ ] Optimize `voxel_table_register` for Span-N: Register only the surface shell (faces) to reduce hash collisions, assuming internal volume is implied.

## Phase 2: Neighbor Discovery (Broadphase Optimization)
- [ ] Refactor `gather_neighbor_voxels` with a fast-path for `span == 1` (direct 26-neighbor offsets).
- [ ] Refactor `gather_neighbor_voxels` for `span > 1` to iterate only the surface shell of the voxel's grid bounds, avoiding the internal volume (`O(Span^2)` vs `O(Span^3)`).
- [ ] Replace `gather_static_voxels_near_point` (used for particle collisions) with a direct 3x3x3 grid lookup around the particle's position, eliminating the large-radius search for Span-N corners.

## Phase 3: Collision Resolution (Solver Simplification)
- [ ] Simplify `solve_static_collisions` for Span-N: Gather unique static neighbors first, then compute a single aggregate "floor plane" or "wall plane" constraint instead of iterating `resolve_span_static_overlap` for every contacting static voxel.
- [ ] Optimize `nudge_voxel_bottom_above_static`: Switch from full-bottom-face scanning to a sparse sampling approach (check 4 corners + center) to determine ground height.

## Phase 4: Structural Logic (Glue & Sleep)
- [ ] Optimize `rebuild_glue_constraints`: Rewrite neighbor finding to check specific adjacent grid cells directly, skipping generic `gather_face_neighbors`.
- [ ] Optimize `voxel_connected_to_static_world`: Replace bounds-based surface loop with direct grid lookups at the voxel's boundary cells.

## Phase 5: Parallelization (CPU Multithreading)
- [ ] **Thread Pool & Job System**: Implement a simple fixed-size thread pool (worker threads = core count) to manage task dispatch.
- [ ] **Parallel Independent Stages**:
  - [ ] Parallelize `integrate_particles`: Independent integration of velocity/gravity per voxel.
  - [ ] Parallelize `solve_static_collisions`: Read-only access to static geometry allows safe parallel execution.
  - [ ] Parallelize `solve_voxel_shape`: Shape matching constraints are internal to each voxel (no neighbor dependencies).
  - [ ] Parallelize `update_particle_velocities`: Final velocity and position commit is independent.
- [ ] **Parallel Inter-Voxel Stages (Jacobi Solver)**:
  - [ ] Refactor `solve_dynamic_collisions` and `solve_voxel_glue` to use a **Jacobi** approach (Accumulate & Apply) instead of Gauss-Seidel (Immediate Apply):
    - [ ] **Accumulate Phase**: Threads calculate constraint deltas (collision/glue) and add them to a thread-local or atomic accumulator (`delta_pos`, `constraint_count`) per particle.
    - [ ] **Apply Phase**: Parallel loop to apply `predicted_pos += delta_pos / constraint_count` after all constraints are evaluated.
  - [ ] *Note:* This avoids complex locking/coloring schemes but may require tuning constraint iteration counts as Jacobi convergence is different.
- [ ] **Structural Logic Isolation**: Keep `split_strained_voxels` and `rebuild_voxel_hash` single-threaded (or strictly serialized) at the start/end of the frame to manage array resizing/reindexing safely.
