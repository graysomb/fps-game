# `fps_ray.c` Removal Candidates

This document tracks dead, uncalled, or obsolete functions, structures, and variables in [`fps_ray.c`](fps_ray.c) identified via Clang AST analysis and call-graph reachability.

### 📊 Summary of Potential Savings
- **Memory Footprint:** **~240 MB** of dead static arrays in the BSS segment.
- **Lines of Code:** **2,697 lines** (~15.6% of the 17,310-line file across 99 functions).
- **Compiler Cleanliness:** Resolves all 89 Clang `-Wunused` warnings.

---

## 💾 1. Unused Static Global Arrays (~240 MB Memory)

These arrays are allocated at file-scope in `fps_ray.c` but are never accessed in active gameplay. Removing them immediately frees this memory:

- [x] `table[HASH_SIZE]` ([fps_ray.c:L3186](fps_ray.c#L3186))
  - **Type:** `Bucket[524288]` (16 bytes per bucket)
  - **Size:** **8.39 MB**
  - **Reason:** Commented in code as `// Deprecated, will be replaced by dynamic_table`. Zero references in the codebase.
- [x] `gluedNeighborHashKeys[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE]` ([fps_ray.c:L2820](fps_ray.c#L2820))
  - **Type:** `int[131072][128]`
  - **Size:** **67.11 MB**
  - **Reason:** Part of obsolete glue neighbor hash. Only accessed by uncalled `glue_neighbor_hash_*` functions.
- [x] `gluedNeighborHashStamp[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE]` ([fps_ray.c:L2822](fps_ray.c#L2822))
  - **Type:** `uint32_t[131072][128]`
  - **Size:** **67.11 MB**
  - **Reason:** Part of obsolete glue neighbor hash. Only accessed by uncalled `glue_neighbor_hash_*` functions.
- [x] `gluedNeighborList[MAX_VOXELS][MAX_FACE_NEIGHBORS]` ([fps_ray.c:L2818](fps_ray.c#L2818))
  - **Type:** `int[131072][64]`
  - **Size:** **33.55 MB**
  - **Reason:** Part of obsolete glue adjacency system. Only accessed by uncalled `glue_adjacency_*` functions.
- [x] `glueConstraintPeakViolation[MAX_VOXELS * 48]` ([fps_ray.c:L2815](fps_ray.c#L2815))
  - **Type:** `float[6291456]`
  - **Size:** **25.17 MB**
  - **Reason:** Only accessed by uncalled `reset_glue_constraint_peaks()`.
- [x] `gluedNeighborHashIndex[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE]` ([fps_ray.c:L2821](fps_ray.c#L2821))
  - **Type:** `uint8_t[131072][128]`
  - **Size:** **16.78 MB**
  - **Reason:** Part of obsolete glue neighbor hash. Only accessed by uncalled `glue_neighbor_hash_*` functions.
- [x] `gluedNeighborRefCounts[MAX_VOXELS][MAX_FACE_NEIGHBORS]` ([fps_ray.c:L2819](fps_ray.c#L2819))
  - **Type:** `uint16_t[131072][64]`
  - **Size:** **16.78 MB**
  - **Reason:** Part of obsolete glue adjacency system. Only accessed by uncalled `glue_adjacency_*` functions.
- [x] `yzPosList[MAX_VOXELS]`, `yzNegList[MAX_VOXELS]` ([fps_ray.c:L14385](fps_ray.c#L14385))
  - **Type:** `int[131072]` × 2
  - **Size:** **1.05 MB**
  - **Reason:** Scratch arrays for uncalled `merge_rects_on_plane()`. Never referenced outside declaration.
- [x] `xzPosList[MAX_VOXELS]`, `xzNegList[MAX_VOXELS]` ([fps_ray.c:L14386](fps_ray.c#L14386))
  - **Type:** `int[131072]` × 2
  - **Size:** **1.05 MB**
  - **Reason:** Scratch arrays for uncalled `merge_rects_on_plane()`. Never referenced outside declaration.
- [x] `xyPosList[MAX_VOXELS]`, `xyNegList[MAX_VOXELS]` ([fps_ray.c:L14387](fps_ray.c#L14387))
  - **Type:** `int[131072]` × 2
  - **Size:** **1.05 MB**
  - **Reason:** Scratch arrays for uncalled `merge_rects_on_plane()`. Never referenced outside declaration.
- [x] `staticBeliefQueue[MAX_VOXELS]` ([fps_ray.c:L972](fps_ray.c#L972))
  - **Type:** `int[131072]`
  - **Size:** **524 KB**
  - **Reason:** Declared but never referenced anywhere in the codebase.
- [x] `gluedNeighborEpoch[MAX_VOXELS]` ([fps_ray.c:L2823](fps_ray.c#L2823))
  - **Type:** `uint32_t[131072]`
  - **Size:** **524 KB**
  - **Reason:** Only accessed by uncalled `glue_neighbor_hash_*` functions.
- [x] `glueAdjacencyDirtyList[MAX_VOXELS]` ([fps_ray.c:L2825](fps_ray.c#L2825))
  - **Type:** `int[131072]`
  - **Size:** **524 KB**
  - **Reason:** Only accessed by uncalled `glue_adjacency_*` functions.
- [x] `glueClusterVisitedTemp[MAX_VOXELS]` ([fps_ray.c:L2827](fps_ray.c#L2827))
  - **Type:** `unsigned char[131072]`
  - **Size:** **131 KB**
  - **Reason:** Only accessed by uncalled `batch_glued_dynamic_voxels()`.
- [x] `glueAdjacencyDirtyFlags[MAX_VOXELS]` ([fps_ray.c:L2824](fps_ray.c#L2824))
  - **Type:** `uint8_t[131072]`
  - **Size:** **131 KB**
  - **Reason:** Only accessed by uncalled `glue_adjacency_*` functions.
- [x] `gluedNeighborCounts[MAX_VOXELS]` ([fps_ray.c:L2817](fps_ray.c#L2817))
  - **Type:** `uint8_t[131072]`
  - **Size:** **131 KB**
  - **Reason:** Only accessed by uncalled `glue_adjacency_*` functions.

---

## ⚙️ 2. Unused Global Variables (Scalar & Config)

### A. Zero Usages Anywhere (Dead Declarations)
- [x] `playerSpawnYaw[MAX_PLAYERS]` ([fps_ray.c:L783](fps_ray.c#L783)) - Spawns compute yaw directly.
- [x] `sfxVolumes` ([fps_ray.c:L1420](fps_ray.c#L1420)) - Unused audio array.
- [x] `debugSpanCollisionLogBudget` ([fps_ray.c:L1370](fps_ray.c#L1370))
- [x] `debugLogGlue` ([fps_ray.c:L1371](fps_ray.c#L1371))
- [x] `skipGlueClusterCollisions` ([fps_ray.c:L1373](fps_ray.c#L1373))
- [x] `debugLogDynamicVoxels` ([fps_ray.c:L1435](fps_ray.c#L1435))
- [x] `debugLogSmushSpawns` ([fps_ray.c:L1441](fps_ray.c#L1441))
- [x] `debugGlueBuildLogBudget` ([fps_ray.c:L1459](fps_ray.c#L1459))
- [x] `debugGlueSolveLogBudget` ([fps_ray.c:L1460](fps_ray.c#L1460))
- [x] `debugGlueBreakLogBudget` ([fps_ray.c:L1461](fps_ray.c#L1461))
- [x] `DEBUG_GLUE_BUILD_LOG_INIT` ([fps_ray.c:L1462](fps_ray.c#L1462))
- [x] `DEBUG_GLUE_SOLVE_LOG_INIT` ([fps_ray.c:L1463](fps_ray.c#L1463))
- [x] `DEBUG_GLUE_BREAK_LOG_INIT` ([fps_ray.c:L1464](fps_ray.c#L1464))
- [x] `DEBUG_ACTIVATION_LOG_INIT` ([fps_ray.c:L1466](fps_ray.c#L1466))

### B. Used Exclusively Inside Dead Functions
- [x] `debugTagOffset[16][3]` ([fps_ray.c:L963](fps_ray.c#L963)) - Only in dead `init_debug_tag_offsets` / `apply_debug_tag_offset`.
- [x] `renderAllDynamicFaces` ([fps_ray.c:L1374](fps_ray.c#L1374)) - Only in dead `compute_dynamic_face_visibility`.
- [x] `glueDirections[3]` ([fps_ray.c:L2213](fps_ray.c#L2213)) - Only in dead `voxels_face_direction`.
- [x] `glueAdjacencyDirtyCount` ([fps_ray.c:L2828](fps_ray.c#L2828)) - Only in dead `glue_adjacency_*`.
- [x] `glueAdjacencyDirtyAll` ([fps_ray.c:L2829](fps_ray.c#L2829)) - Only in dead `glue_adjacency_*`.

---

## 🔤 3. Unused Local Variables (Compiler Warnings)

- [x] `z1` ([fps_ray.c:L6916](fps_ray.c#L6916)) in `carveFortGateOnX` - Calculated but never read.
- [x] `minx`, `maxx`, `miny`, `maxy`, `minz`, `maxz` ([fps_ray.c:L9091-L9096](fps_ray.c#L9091-L9096)) in `update_projectiles` - Extents computed but unused.
- [x] `scale` ([fps_ray.c:L14600](fps_ray.c#L14600)) in `render_world_cubes_instanced` - Assigned `VOXEL_SIZE` but never read.
- [x] `countFrame` ([fps_ray.c:L16492](fps_ray.c#L16492)) in `main` - Incremented on debug frames but never used in any logic.

---

## 🧱 4. Unused Structures & Types

- [x] `GlueDirection` ([fps_ray.c:L2207-L2211](fps_ray.c#L2207-L2211))
  ```c
  typedef struct {
      int dx, dy, dz;
      int faceA[4];
      int faceB[4];
  } GlueDirection;
  ```
  *Reason:* Only used by `glueDirections` and `voxels_face_direction`, both of which are uncalled.
- [x] `enum SfxId` ([fps_ray.c:L1377](fps_ray.c#L1377))
  *Reason:* Sound functions take raw integer literals; the enum tags are not used symbolically.

---

## ✂️ 5. Uncalled / Unreachable Functions (99 Functions, 2,697 Lines)

### A. Obsolete World Builders (19 functions, 715 lines)
Active world generation uses `buildTestWorld` and `buildDemo`. These older prototypes are never called:
- [x] `buildProceduralWorld` ([fps_ray.c:L7065-L7105](fps_ray.c#L7065-L7105)) (41 lines)
- [x] `buildBloodWorld` ([fps_ray.c:L7107-L7207](fps_ray.c#L7107-L7207)) (101 lines)
- [x] `buildDebugWorld` ([fps_ray.c:L7209-L7438](fps_ray.c#L7209-L7438)) (230 lines)
- [x] `buildStackedPillar` ([fps_ray.c:L6805-L6825](fps_ray.c#L6805-L6825)) (21 lines)
- [x] `buildSplitPlatform` ([fps_ray.c:L6829-L6872](fps_ray.c#L6829-L6872)) (44 lines)
- [x] `buildLeg2x2Notched` ([fps_ray.c:L6874-L6883](fps_ray.c#L6874-L6883)) (10 lines)
- [x] `buildFortWalls` ([fps_ray.c:L6898-L6911](fps_ray.c#L6898-L6911)) (14 lines)
- [x] `carveFortGateOnX` ([fps_ray.c:L6914-L6930](fps_ray.c#L6914-L6930)) (17 lines)
- [x] `carveFortWindows` ([fps_ray.c:L6933-L6953](fps_ray.c#L6933-L6953)) (21 lines)
- [x] `buildFortRoof` ([fps_ray.c:L6955-L6963](fps_ray.c#L6955-L6963)) (9 lines)
- [x] `buildGateFrame` ([fps_ray.c:L6965-L6981](fps_ray.c#L6965-L6981)) (17 lines)
- [x] `add_static_box_at_grid` ([fps_ray.c:L6662-L6670](fps_ray.c#L6662-L6670)) (9 lines)
- [x] `add_static_disc_at_grid` ([fps_ray.c:L6673-L6683](fps_ray.c#L6673-L6683)) (11 lines)
- [x] `add_static_cylinder_at_grid` ([fps_ray.c:L6688-L6692](fps_ray.c#L6688-L6692)) (5 lines)
- [x] `add_static_ring_at_grid` ([fps_ray.c:L6698-L6712](fps_ray.c#L6698-L6712)) (15 lines)
- [x] `add_static_boulder_at_grid` ([fps_ray.c:L6716-L6731](fps_ray.c#L6716-L6731)) (16 lines)
- [x] `add_static_ring_with_door_at_grid` ([fps_ray.c:L6739-L6760](fps_ray.c#L6739-L6760)) (22 lines)
- [x] `add_static_hollow_cylinder_at_grid` ([fps_ray.c:L6766-L6774](fps_ray.c#L6766-L6774)) (9 lines)
- [x] `add_static_arch_at_grid` ([fps_ray.c:L6779-L6801](fps_ray.c#L6779-L6801)) (23 lines)

### B. Obsolete Dynamic Glue & Adjacency System (22 functions, 574 lines)
- [x] `batch_glued_dynamic_voxels` ([fps_ray.c:L11117-L11279](fps_ray.c#L11117-L11279)) (163 lines)
- [x] `gather_neighbor_voxels` ([fps_ray.c:L11569-L11601](fps_ray.c#L11569-L11601)) (33 lines)
- [x] `particles_are_glued_pair` ([fps_ray.c:L11544-L11566](fps_ray.c#L11544-L11566)) (23 lines)
- [x] `voxels_face_direction` ([fps_ray.c:L11470-L11524](fps_ray.c#L11470-L11524)) (55 lines)
- [x] `voxels_are_glued` ([fps_ray.c:L11910-L11935](fps_ray.c#L11910-L11935)) (26 lines)
- [x] `voxels_share_edge_or_corner` ([fps_ray.c:L11940-L11950](fps_ray.c#L11940-L11950)) (11 lines)
- [x] `compute_voxel_center_and_mass` ([fps_ray.c:L11953-L11969](fps_ray.c#L11953-L11969)) (17 lines)
- [x] `deactivate_glue_constraints_between` ([fps_ray.c:L10398-L10427](fps_ray.c#L10398-L10427)) (30 lines)
- [x] `rebuild_glue_adjacency_if_dirty` ([fps_ray.c:L3032-L3068](fps_ray.c#L3032-L3068)) (37 lines)
- [x] `mark_glue_adjacency_dirty_for_voxel` ([fps_ray.c:L2835-L2853](fps_ray.c#L2835-L2853)) (19 lines)
- [x] `reset_glue_constraint_peaks` ([fps_ray.c:L2831-L2833](fps_ray.c#L2831-L2833)) (3 lines)
- [x] `glue_adjacency_add_ref_oneway` ([fps_ray.c:L2948-L2980](fps_ray.c#L2948-L2980)) (33 lines)
- [x] `glue_adjacency_add_ref_pair` ([fps_ray.c:L2982-L2985](fps_ray.c#L2982-L2985)) (4 lines)
- [x] `glue_adjacency_remove_ref_oneway` ([fps_ray.c:L2987-L3025](fps_ray.c#L2987-L3025)) (39 lines)
- [x] `glue_adjacency_remove_ref_pair` ([fps_ray.c:L3027-L3030](fps_ray.c#L3027-L3030)) (4 lines)
- [x] `glue_neighbor_hash_seed` ([fps_ray.c:L2863-L2865](fps_ray.c#L2863-L2865)) (3 lines)
- [x] `glue_neighbor_hash_find` ([fps_ray.c:L2868-L2904](fps_ray.c#L2868-L2904)) (37 lines)
- [x] `glue_neighbor_hash_insert` ([fps_ray.c:L2906-L2919](fps_ray.c#L2906-L2919)) (14 lines)
- [x] `glue_neighbor_hash_rebuild` ([fps_ray.c:L2921-L2928](fps_ray.c#L2921-L2928)) (8 lines)
- [x] `glue_direction_label_from_delta` ([fps_ray.c:L2220-L2228](fps_ray.c#L2220-L2228)) (9 lines)

### C. Obsolete Physics & Particle Routines (16 functions, 638 lines)
- [x] `physics_step` ([fps_ray.c:L8788-L8915](fps_ray.c#L8788-L8915)) (128 lines: old physics solver superseded by `pbd_step_simulation`)
- [x] `solve_voxel_shape` ([fps_ray.c:L9785-L9952](fps_ray.c#L9785-L9952)) (168 lines)
- [x] `reset_voxel_shape_to_rest` ([fps_ray.c:L9763-L9782](fps_ray.c#L9763-L9782)) (20 lines)
- [x] `solve_dynamic_collisions` ([fps_ray.c:L12062-L12065](fps_ray.c#L12062-L12065)) (4 lines)
- [x] `apply_uniform_velocity` ([fps_ray.c:L12180-L12190](fps_ray.c#L12180-L12190)) (11 lines)
- [x] `split_voxel_at` ([fps_ray.c:L12192-L12198](fps_ray.c#L12192-L12198)) (7 lines)
- [x] `split_strained_voxels` ([fps_ray.c:L12200-L12203](fps_ray.c#L12200-L12203)) (4 lines)
- [x] `cull_dust_voxels` ([fps_ray.c:L12205-L12245](fps_ray.c#L12205-L12245)) (41 lines)
- [x] `accumulate_simulated_corner_deltas` ([fps_ray.c:L9449-L9468](fps_ray.c#L9449-L9468)) (20 lines)
- [x] `mark_simulated_particles` ([fps_ray.c:L9470-L9480](fps_ray.c#L9470-L9480)) (11 lines)
- [x] `apply_interior_sync` ([fps_ray.c:L9482-L9506](fps_ray.c#L9482-L9506)) (25 lines)
- [x] `reset_particle_mass_and_flags` ([fps_ray.c:L9330-L9341](fps_ray.c#L9330-L9341)) (12 lines)
- [x] `apply_shell_effective_mass` ([fps_ray.c:L9360-L9382](fps_ray.c#L9360-L9382)) (23 lines)
- [x] `update_voxel_coarsening_state` ([fps_ray.c:L9261-L9328](fps_ray.c#L9261-L9328)) (68 lines)
- [x] `copy_particle_snapshot_range` ([fps_ray.c:L12393-L12397](fps_ray.c#L12393-L12397)) (5 lines)

### D. Obsolete Face Rendering & Greedy Meshing (5 functions, 287 lines)
- [x] `merge_rects_on_plane` ([fps_ray.c:L14214-L14382](fps_ray.c#L14214-L14382)) (169 lines)
- [x] `compute_dynamic_face_visibility` ([fps_ray.c:L13901-L14002](fps_ray.c#L13901-L14002)) (102 lines)
- [x] `drawCubeMan` ([fps_ray.c:L14024-L14090](fps_ray.c#L14024-L14090)) (67 lines)
- [x] `compute_voxel_face_visibility` ([fps_ray.c:L14005-L14021](fps_ray.c#L14005-L14021)) (17 lines)
- [x] `drawCubeEdges` ([fps_ray.c:L13882-L13898](fps_ray.c#L13882-L13898)) (17 lines)

### E. Dead Spatial Hash Wrappers (5 functions, 30 lines)
- [x] `table_remove` ([fps_ray.c:L3403-L3409](fps_ray.c#L3403-L3409)) (7 lines)
- [x] `table_set` ([fps_ray.c:L3296-L3298](fps_ray.c#L3296-L3298)) (3 lines)
- [x] `hashVoxel` ([fps_ray.c:L3291-L3293](fps_ray.c#L3291-L3293)) (3 lines)
- [x] `unit_voxel_grid_index` ([fps_ray.c:L4280-L4283](fps_ray.c#L4280-L4283)) (4 lines)
- [x] `emit_static_voxels_from_units` ([fps_ray.c:L4322-L4333](fps_ray.c#L4322-L4333)) (12 lines)

### F. Dead Contact & Geometry Math Utilities (17 functions, 194 lines)
- [x] `voxel_predicted_bounds` ([fps_ray.c:L1701-L1724](fps_ray.c#L1701-L1724)) (24 lines)
- [x] `voxel_visibility_bounds` ([fps_ray.c:L1727-L1736](fps_ray.c#L1727-L1736)) (10 lines)
- [x] `axis_contact_state` ([fps_ray.c:L1741-L1755](fps_ray.c#L1741-L1755)) (15 lines)
- [x] `bounds_overlap` ([fps_ray.c:L1759-L1761](fps_ray.c#L1759-L1761)) (3 lines)
- [x] `bounds_overlap_length` ([fps_ray.c:L1765-L1770](fps_ray.c#L1765-L1770)) (6 lines)
- [x] `face_blocked_by_voxel` ([fps_ray.c:L2016-L2065](fps_ray.c#L2016-L2065)) (50 lines)
- [x] `voxel_touching_axes` ([fps_ray.c:L2069-L2113](fps_ray.c#L2069-L2113)) (45 lines)
- [x] `ranges_touch_int` ([fps_ray.c:L2115-L2117](fps_ray.c#L2115-L2117)) (3 lines)
- [x] `ranges_overlap_int` ([fps_ray.c:L2119-L2121](fps_ray.c#L2119-L2121)) (3 lines)
- [x] `voxels_share_edge_or_corner_rest` ([fps_ray.c:L2123-L2147](fps_ray.c#L2123-L2147)) (25 lines)
- [x] `voxels_share_face_rest` ([fps_ray.c:L2149-L2168](fps_ray.c#L2149-L2168)) (20 lines)
- [x] `voxel_rest_axis_min` ([fps_ray.c:L2170-L2176](fps_ray.c#L2170-L2176)) (7 lines)
- [x] `voxel_rest_axis_max` ([fps_ray.c:L2178-L2184](fps_ray.c#L2178-L2184)) (7 lines)
- [x] `voxel_rest_corner_axis_coord` ([fps_ray.c:L2186-L2191](fps_ray.c#L2186-L2191)) (6 lines)
- [x] `voxel_rest_corner_world` ([fps_ray.c:L2193-L2199](fps_ray.c#L2193-L2199)) (7 lines)
- [x] `face_normal_predicted` ([fps_ray.c:L3071-L3087](fps_ray.c#L3071-L3087)) (17 lines)
- [x] `face_normal_rest` ([fps_ray.c:L3090-L3103](fps_ray.c#L3090-L3103)) (14 lines)
- [x] `get_face_corners_for_direction` ([fps_ray.c:L3108-L3113](fps_ray.c#L3108-L3113)) (6 lines)
- [x] `order_coarse_fine_pair` ([fps_ray.c:L3118-L3123](fps_ray.c#L3118-L3123)) (6 lines)
- [x] `mix_color_channel` ([fps_ray.c:L6529-L6531](fps_ray.c#L6529-L6531)) (3 lines)

### G. Dead Debug & Test Stub Functions (15 functions, 159 lines)
- [x] `debug_cluster_tag_label` ([fps_ray.c:L1490-L1508](fps_ray.c#L1490-L1508)) (19 lines)
- [x] `debug_should_log_tag_break` ([fps_ray.c:L1511-L1519](fps_ray.c#L1511-L1519)) (9 lines)
- [x] `init_debug_tag_offsets` ([fps_ray.c:L6589-L6602](fps_ray.c#L6589-L6602)) (14 lines)
- [x] `apply_debug_tag_offset` ([fps_ray.c:L6579-L6586](fps_ray.c#L6579-L6586)) (8 lines)
- [x] `add_dynamic_span_voxel_at_grid_tag` ([fps_ray.c:L6606-L6618](fps_ray.c#L6606-L6618)) (13 lines)
- [x] `add_dynamic_span_voxel_at_grid` ([fps_ray.c:L6622-L6624](fps_ray.c#L6622-L6624)) (3 lines)
- [x] `add_dynamic_unit_block_tag` ([fps_ray.c:L6629-L6639](fps_ray.c#L6629-L6639)) (11 lines)
- [x] `add_dynamic_unit_block` ([fps_ray.c:L6643-L6645](fps_ray.c#L6643-L6645)) (3 lines)
- [x] `build_oblique_voxel_pyramid` ([fps_ray.c:L6533-L6576](fps_ray.c#L6533-L6576)) (44 lines)
- [x] `grid_to_world_g` ([fps_ray.c:L6654-L6656](fps_ray.c#L6654-L6656)) (3 lines)
- [x] `update_static_voxel_belief` ([fps_ray.c:L5767-L5807](fps_ray.c#L5767-L5807)) (41 lines)
- [x] `activate_static_neighbors_of_region` ([fps_ray.c:L8704-L8769](fps_ray.c#L8704-L8769)) (66 lines)
- [x] `particle_hash_clear_range` ([fps_ray.c:L2547-L2553](fps_ray.c#L2547-L2553)) (7 lines)
- [x] `face_local_coords` ([fps_ray.c:L2665-L2679](fps_ray.c#L2665-L2679)) (15 lines)
- [x] `voxel_center_near_grid` ([fps_ray.c:L11102-L11114](fps_ray.c#L11102-L11114)) (13 lines)

---

## 🎨 6. Unused Shader Files (18 Shaders, 3,085 Lines)

The active game engine ([`fps_ray.c`](fps_ray.c) and [`physics_gpu_common.inc`](physics_gpu_common.inc)) only loads:
- `shaders/instanced_voxel_hack.vert` & `.frag` (dynamic voxels)
- `shaders/orb.vert` & `.frag` (orb rendering)
- `shaders/voxel_simple.vert` & `.frag` (greedy meshed static voxels)
- `shaders/pbd/pbd_pipeline.comp` & `pbd_pipeline.metal` (GPU physics compute pipeline)

All other 18 shader files in `shaders/` are unused leftovers:

### A. Unconnected Placeholder Shaders (46 lines)
- [x] [`shaders/skybox_shader.vert`](shaders/skybox_shader.vert) (35 lines) - Created for planned skybox feature, but never loaded.
- [x] [`shaders/skybox_shader.frag`](shaders/skybox_shader.frag) (11 lines) - Created for planned skybox feature, but never loaded.

### B. Legacy OpenGL Compute Shaders (1,765 lines)
*(Leftovers from an earlier compute-shader prototype referenced only in orphaned `DestructiveCSSceneObject.h`. The active engine uses `shaders/pbd/pbd_pipeline.*`)*
- [x] [`shaders/particle_collision.comp`](shaders/particle_collision.comp) (789 lines)
- [x] [`shaders/particle_vgs_face.comp`](shaders/particle_vgs_face.comp) (404 lines)
- [x] [`shaders/particle_vgs_voxel.comp`](shaders/particle_vgs_voxel.comp) (277 lines)
- [x] [`shaders/particle_uniform_grid.comp`](shaders/particle_uniform_grid.comp) (149 lines)
- [x] [`shaders/particle_make_static.comp`](shaders/particle_make_static.comp) (97 lines)
- [x] [`shaders/particle_prefix_sum.comp`](shaders/particle_prefix_sum.comp) (49 lines)

### C. Legacy Rendering Shaders (1,274 lines)
*(Leftovers from an earlier custom OpenGL rendering pipeline, superseded by Raylib and `voxel_simple`)*
- [x] [`shaders/voxel.vert`](shaders/voxel.vert) (139 lines) & [`shaders/voxel.frag`](shaders/voxel.frag) (208 lines) - Replaced by `voxel_simple`.
- [x] [`shaders/voxel_skin.vert`](shaders/voxel_skin.vert) (140 lines) & [`shaders/voxel_skin.frag`](shaders/voxel_skin.frag) (200 lines)
- [x] [`shaders/particle.vert`](shaders/particle.vert) (72 lines) & [`shaders/particle.frag`](shaders/particle.frag) (203 lines)
- [x] [`shaders/object_shader.vert`](shaders/object_shader.vert) (42 lines) & [`shaders/object_shader.frag`](shaders/object_shader.frag) (206 lines)
- [x] [`shaders/debug_AABB.vert`](shaders/debug_AABB.vert) (55 lines) & [`shaders/debug_AABB.frag`](shaders/debug_AABB.frag) (9 lines)

### D. Orphaned Prototype Header (1,361 lines)
- [x] [`DestructiveCSSceneObject.h`](DestructiveCSSceneObject.h) (1,361 lines) - Old C++ header that referenced the legacy shaders above. Never included or compiled anywhere in the project.

---

## ⚡ 7. Solver Architecture Rework & VGS Glue Simplification (~384 MB Memory)

As part of the solver overhaul (inverting particle/voxel coupling so **Particle = Persistent Node** and **Voxel = Breakable VGS Glue** without runtime particle cloning):
- [x] `glueConstraints[MAX_VOXELS * 48]` (~384 MB BSS memory) and `GlueConstraint` struct.
- [x] `compact_glue_constraints()`, `deactivate_glue_constraints_for_voxel()`, and `rebuild_glue_constraints()` unneeded resets.
- [x] `detach_face_particles()`, `break_face_link()`, and `detach_all_glue_faces()` (runtime particle-duplication splitting routines).
- [x] `process_break_masks()` and `gather_voxel_break_masks()` replaced with parallel `evaluate_voxel_fracture()`.
- [x] Atomic float accumulator contention in Jacobi shape matching eliminated via non-atomic 8-octant `jacobi_slots[8]` scatter-gather.
- [x] Unglued free particles (`glue_count <= 0`) rendered as 3D triangles (`DrawTriangle3D`).

