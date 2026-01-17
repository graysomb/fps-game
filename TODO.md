# Multithreading Plan for PBD Loop (CPU & Future GPU)

## Goal
Safe multithreading of the PBD physics loop with minimal overhead, preparing for a future Compute Shader migration.
The loop uses a Jacobi solver approach (Accumulate -> Apply), which is inherently parallel-friendly compared to Gauss-Seidel.

## 1. Infrastructure & Helpers
- [ ] **Thread Pool**: Implement a lightweight, persistent thread pool (worker threads spin/wait on condition variable). Avoid spawning/joining threads per frame.
- [ ] **Atomic Float Add**: Implement `atomic_add_float(volatile float* target, float value)` using a CAS (Compare-And-Swap) loop.
    - *Why*: Standard `stdatomic.h` does not support `atomic_fetch_add` for floats.
    - *Shader Parity*: Maps to `atomicAdd` (GLSL) or `InterlockedAdd` (HLSL) for floats (often needing a CAS loop implementation in shaders too, or extensions).

## 2. Kernel B: Particle Operations (Thread over Active Particles)
*Scope: Operations iterating `active_particles` array.*

### 2.1 Integration & Reset
- [ ] **Parallelize `integrate_particles`**: Simple parallel-for. No dependencies.
- [ ] **Parallelize `reset_particle_accumulators`**: Simple parallel-for. Clears `corr_sum` and `corr_weight`.

### 2.2 Spatial Hashing (`build_particle_hash`)
- [ ] **Parallel Clear**: Parallel-for to set `particle_hash_head[]` to -1.
- [ ] **Parallel Build**: Parallel-for over particles.
    - Calculate hash `h`.
    - **Atomic Insert**: Use `atomic_exchange` (or `atomic_compare_exchange`) on `particle_hash_head[h]` to safely insert the particle into the linked list.
    - *Note*: Order of particles in the bin will become non-deterministic, which is acceptable.

### 2.3 Collision Gathering (`gather_particle_collisions`)
- [ ] **Parallel Loop**: Iterate particles in chunks.
- [ ] **Read Phase**: Read neighbors via `particle_hash` (safe, read-only during this phase).
- [ ] **Write Phase**:
    - When a collision is found, calculate correction `delta`.
    - **Atomic Accumulate**: Use `atomic_add_float` to update `p->corr_sum` (x, y, z) and `p->corr_weight` for *both* particles involved.
    - *Optimization*: Compute `delta` locally, only lock/atomic-add once per interaction.

### 2.4 Apply Accumulators (`apply_particle_accumulators`)
- [ ] **Parallelize**: Simple parallel-for.
    - Read `corr_sum` / `corr_weight`.
    - Update `predicted_pos`.
    - Reset accumulators here (optimization to merge with Reset step, reducing memory passes).

## 3. Kernel A: Voxel Operations (Thread over Active Voxels)
*Scope: Operations iterating `voxels` array.*

### 3.1 Shape Constraints (`gather_voxel_shape_constraints`)
- [ ] **Parallel Loop**: Iterate active voxels (skip inactive/air).
- [ ] **Read Phase**: Read positions of the 8 corner particles.
- [ ] **Compute**: Run VGS (Volume/Gradient/Strain) constraint logic logic locally.
- [ ] **Write Phase**:
    - Calculate corrections for the 8 particles.
    - **Atomic Accumulate**: Use `atomic_add_float` to update `corr_sum` and `corr_weight` on the shared corner particles.

## 4. Execution Flow (Per Frame)
1. `Dispatch(Kernel_B_Integrate)`
2. `Dispatch(Kernel_Hash_Clear)`
3. `Dispatch(Kernel_Hash_Build)`
4. **Loop (Constraint Iterations)**:
    a. `Dispatch(Kernel_B_Reset_Accumulators)` (Or merged into Apply)
    b. `Dispatch(Kernel_B_Collisions)`  <-- *Heavy contention potential*
    c. `Dispatch(Kernel_A_Constraints)` <-- *Heavy contention potential*
    d. `Barrier` (Wait for all accumulators to be finalized)
    e. `Dispatch(Kernel_B_Apply)`
5. `Dispatch(Kernel_B_Update_Velocities)`

## 5. Future Shader Migration Strategy
- **Data Layout**: Transition `Particle` structs to SoA (Structure of Arrays) in usage (e.g., `pos_x[]`, `pos_y[]`, `pos_z[]`). This is cache-friendly for CPU and required for coalesced access on GPU.
- **Compute Shaders**:
    - `Kernel A` becomes a Compute Shader dispatching `(VOXEL_COUNT / 64, 1, 1)`.
    - `Kernel B` becomes a Compute Shader dispatching `(PARTICLE_COUNT / 64, 1, 1)`.
    - **SSBOs**: Particles and Voxels reside in Shader Storage Buffer Objects.
    - **Grid**: Use a standard Grid/Binning compute shader approach (e.g., Counting Sort or Atomic binning).
