# Updated Break Flow (No `process_pending_breaks()`)

The break flow is now handled directly in `solve_voxel_shape()` on a per-face basis. There is no longer a `pending_full_break` flag or a `process_pending_breaks()` pass.

### How Breaking Works Now

1. **Strain + Shear Checks:** During `solve_voxel_shape()`, the voxel’s strain and shear are evaluated.
2. **Per-Face Decisions:** Each affected face is checked independently. If a threshold is exceeded for a given axis or shear plane, only the corresponding faces are broken.
3. **Breaking a Face:** `break_face_link(voxel, face)` performs the actual detach:
   - Clones shared particles on that face via `detach_face_particles()`.
   - Clears `glued_faces` on both the voxel and its neighbor.

### In Summary

Instead of breaking all faces at once, the simulation now breaks only the faces whose particles exceed stress/strain thresholds. This yields more localized fractures while still using the same particle-cloning logic.
