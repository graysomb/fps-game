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