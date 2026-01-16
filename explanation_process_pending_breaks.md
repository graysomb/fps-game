# How `process_pending_breaks()` Currently Works

The `process_pending_breaks()` function is a crucial part of the voxel destruction simulation. It is responsible for carrying out the "breaking" of voxel connections when a voxel has been flagged for destruction due to excessive stress.

Here is a breakdown of how it works:

1.  **Iterates Through All Voxels:** The function begins by looping through every voxel currently in the simulation.

2.  **Checks for Pending Break Flag:** For each voxel, it checks the boolean flag `voxel->pending_full_break`. This flag is the key. It's set to `true` within the `solve_voxel_shape()` function if the voxel's deformation (strain or shear) exceeds predefined thresholds (`STRAIN_BREAK_THRESHOLD` or `SHEAR_BREAK_THRESHOLD`).

3.  **Conditional Execution:** If `pending_full_break` is `false`, the function does nothing for that voxel and continues to the next one. This is an efficient way to skip voxels that are not under stress.

4.  **Breaks All Face Links:** If `pending_full_break` is `true`, the function then enters a loop that iterates through all 6 faces of the voxel (from `face = 0` to `5`). In each iteration, it calls `break_face_link(voxel, face)`.

5.  **`break_face_link()` Action:** The `break_face_link()` function does the actual work of severing the connection. It:
    *   Finds the neighboring voxel on that specific face.
    *   Calls `detach_face_particles()` for the current voxel's face and the neighbor's corresponding face. This function replaces the shared particles on that face with new, un-shared particles, effectively breaking the "glue".
    *   Sets the `glued_faces` flag for both voxels to `false` for that face.

6.  **Resets the Flag:** After iterating through all 6 faces and breaking all links, the function resets `voxel->pending_full_break` to `false`. This is important to ensure that the voxel doesn't get broken again in the next simulation step unless it is flagged again due to new stresses.

### In Summary

The current implementation of `process_pending_breaks()` is a simple but effective "all or nothing" approach. When a voxel is determined to be unstable, it is completely severed from all of its neighbors at once. It doesn't allow for partial breaks (e.g., only one face breaking). The decision to break is made elsewhere (`solve_voxel_shape`), and `process_pending_breaks` is the function that executes that decision.
