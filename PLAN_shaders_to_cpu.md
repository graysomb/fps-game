# Shader-to-CPU Port Plan

1. **Shared Types Header**
   - Mirror `XpbdParticle`, `FrictionConstraint`, `FaceConstraint`, `VoxelConstraint`, and `DESTRUCTIVE_SYSTEM_DATA` from the GLSL compute shaders in `shaders/particle_collision.comp`, `shaders/particle_vgs_face.comp`, and `shaders/particle_vgs_voxel.comp`.
   - Provide host-friendly layout (packed structs, C-friendly vector types) that matches the shader memory layout. Reuse existing math types where available.

2. **Math Helpers**
   - Create a lightweight inline math helper set (vectors, dot/cross, length, normalize, clamp, mix, project) mirroring the GLSL utilities leveraged by the shaders (`proj`, `div_free`, etc.).
   - Prefer wrapping existing project math functionality; otherwise implement equivalents in a dedicated header the solvers can include.

3. **Collision Solver Header**
   - Translate `SolveConstraintsGrid`, `WallCollision`, `ProjectileCollision`, `UpdateVandX`, and supporting functions from `shaders/particle_collision.comp` into a CPU header (`pbd_collisions.h`).
   - Adapt uniform-grid iteration over `mStart`, `mCount`, and `mContent` buffers to the CPU-side data the PBD loop already builds.

4. **Voxel Shape Solver Header**
   - Port the voxel Gram–Schmidt solver (`VGS`) from `shaders/particle_vgs_voxel.comp` into `pbd_voxel_shape.h`, preserving orthogonalisation, volume correction, and static-voxel early-outs.
   - Expose an interface that accepts voxel particle pointers plus the system settings so the PBD loop can call it directly.

5. **Face Constraint Solver Header**
   - Convert `ProjFaceConstraint` (and supporting helpers such as constraint breaking and strain tests) from `shaders/particle_vgs_face.comp` into `pbd_face_constraints.h`.
   - Surface functions for evaluating strain limits, toggling constraint activity, and applying the projective corrections used in the shader.

6. **Integrate with PBD Loop**
   - Replace ad-hoc CPU implementations with calls into the new headers inside the existing PBD Gauss–Seidel loop (predict → collisions → voxel → face → velocity update).
   - Ensure partition ordering (collision partitions, voxel iterations, axis-wise face batches) mirrors the shader pipeline.

7. **Validation and Instrumentation**
   - Build small regression helpers comparing GPU and CPU corrections for reference scenarios (single voxel, face pair, collision pair) to guard against regressions.
   - Add temporary debug hooks or assertions (e.g., volume sign, penetration thresholds) while integrating, then remove or gate them before final submission.

