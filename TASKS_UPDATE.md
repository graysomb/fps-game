# Update Plan: Matter Bullets & Tether Mechanics

## 1. Matter Bullets
**Goal:** Visually replace the "weird boxel thing" with a glowing blue orb, prevent it from activating static voxels, and implement damage/cost logic.

### 1.1 Assets & Shaders
- [x] **Create Shaders:** Created `shaders/orb.vert` and `shaders/orb.frag` for the glowing blue orb effect.
    - `orb.vert`: Standard transform + normal/view direction calculation.
    - `orb.frag`: Blue color with fresnel-based rim lighting and center glow.

### 1.2 `fps_ray.c` Implementation
- [ ] **Load Shader:**
    - Initialize `Shader orbShader` in `main()` or `InitGame()`.
    - Load `shaders/orb.vert` and `shaders/orb.frag`.
- [ ] **Update `FireVoxel`:**
    - Ensure bullet cost (`MATTER_SHOT_COST`) is deducted (verify existing logic).
    - Set `voxel->isBullet = true`.
    - Set `voxel->type = 0` (assuming type 0 is the blue bullet).
    - **Crucial:** Ensure `voxel->glueEligible = false` so it doesn't stick to things.
- [ ] **Update Physics/Collision (`handle_pbd_projectile_hits` or `physics_step`):**
    - **Collision:** Implement Raycast vs Player AABB/Capsule.
        - Iterate through all players.
        - Check collision with bullet ray.
        - If hit:
            - Deduct health/matter from target.
            - If `matter <= 0`, set `isExposed = true`.
            - Destroy bullet.
    - **No Activation:** In `activate_static_voxels_near_dynamic`, add a check:
        - `if (dynamic->isBullet) continue;`
        - This prevents bullets from waking up static terrain on flyby or impact.
    - **Static Collision:** Ensure bullet destroys/impacts static voxels without waking neighbors (unless intended to destroy).
- [ ] **Rendering (`DrawVoxels`):**
    - Iterate through all active voxels.
    - If `v->isBullet && v->type == 0`:
        - **Skip** the standard `DrawMeshInstanced` logic for this voxel.
        - **Draw** using `orbShader` with a Sphere mesh (use `GenMeshSphere` or `DrawSphere` helper if compatible with custom shader, otherwise use `DrawMesh`).

## 2. Tether Mechanics ("Rip")
**Goal:** Tether should rip a *single* voxel from the static structure without activating the entire cluster.

### 2.1 `fps_ray.c` Logic
- [ ] **Create `rip_single_static_voxel(int voxel_idx, int activator)`:**
    - **Input:** Index of the static voxel to rip.
    - **Logic:**
        1. `remove_voxel_index(voxel_idx)`: Remove it from the static arrays.
        2. `addVoxel(...)`: Spawn a new *dynamic* voxel at the same position.
        3. **Do NOT** call `glue_dynamic_voxel_to_static_neighbors` or `collect_static_activation_cluster`.
        4. Return the index of the new dynamic voxel.
- [ ] **Update `start_tether`:**
    - Use `first_voxel_hit` to find the target.
    - If target is static:
        - Call `rip_single_static_voxel`.
        - Set `player->tetherVoxel` to the returned index.
        - Set `player->tetherHolding = true`.
- [ ] **Prevent Neighbor Wake-up:**
    - In `activate_static_voxels_near_dynamic`:
        - Add check: `if (tetherTag[i] > 0) continue;` (Assuming `tetherTag` marks held voxels).
        - This ensures the single ripped voxel doesn't immediately wake up the wall it was ripped from.

## 3. General Cleanup
- [ ] Remove old debug code or unused voxel types if they conflict with the new bullet/tether logic.
- [ ] Verify `MATTER_SHOT_COST` definition (ensure it exists).
