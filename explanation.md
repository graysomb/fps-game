# Explanation of `simulate_voxel_pbd()`

The `simulate_voxel_pbd` function implements a Position-Based Dynamics (PBD) simulation for voxels. PBD is a method for simulating physical phenomena such as cloth, fluids, and deformable objects. Instead of directly computing forces and integrating them to update velocities and positions, PBD directly manipulates the positions of particles to satisfy a set of constraints. This approach is known for its stability and control over the simulation's behavior.

The simulation is broken down into a series of steps, each of which is further broken down into sub-steps. Within each sub-step, the following operations are performed:

1.  **`integrate_particles(sub_dt)`**: This function predicts the next position of each particle based on its current velocity and external forces (gravity). This is the "dynamics" part of PBD.
2.  **`solve_particle_collisions(sub_dt)`**: This function resolves collisions between particles and with the environment. It adjusts the predicted positions of particles to prevent them from interpenetrating.
3.  **`solve_voxel_shape(voxel)`**: This is the core of the voxel simulation. It applies shape-matching constraints to each voxel. The goal is to preserve the original shape of the voxel, even as it deforms. This is where the "glue" and "shape constraints" come into play. The `solve_voxel_shape` function uses a technique called Voxel Gram-Schmidt shape matching to keep each cell near its rest state. It also checks for strain and shear, and if they exceed a certain threshold, it breaks only the affected faces by detaching their shared particles.
4.  **`update_particle_velocities(sub_dt)`**: After the particle positions have been adjusted to satisfy all constraints, this function updates the particle velocities based on the change in position.

Now, let's look at each of the sub-functions in more detail:

### `integrate_particles(sub_dt)`

This function implements a simple semi-implicit Euler integration scheme. For each particle, it updates its predicted position based on its current velocity and the force of gravity. The `VELOCITY_DAMPING` constant is used to gradually reduce the velocity of the particles over time, which helps to stabilize the simulation.

### `solve_particle_collisions(sub_dt)`

This function resolves collisions between particles and between particles and the environment. It does this by iterating through all pairs of particles and checking if they are interpenetrating. If they are, it adjusts their predicted positions to move them apart. It also handles collisions with the floor and player bounding boxes.

### `solve_voxel_shape(voxel)`

This function is the heart of the voxel simulation. It implements a shape-matching constraint that tries to preserve the original shape of each voxel. It does this by computing the current principal axes of the voxel and then comparing them to the rest shape. If the voxel has deformed, the function computes a correction that will move the particles back towards their rest positions.

The `solve_voxel_shape` function also implements the "glue" constraint. It checks the strain and shear on each face of the voxel. If the strain or shear exceeds a certain threshold, it immediately breaks the affected face by detaching the shared particles with `break_face_link`.

### `update_particle_velocities(sub_dt)`

This function updates the velocity of each particle based on the change in its position. This is necessary because the PBD solver directly manipulates the positions of the particles, so the velocities need to be updated to reflect these changes.

### How it all works together

The `simulate_voxel_pbd` function works by first predicting the next position of each particle, then resolving collisions and applying shape-matching constraints to adjust the predicted positions. Finally, it updates the velocities of the particles based on the adjusted positions. This process is repeated for a fixed number of sub-steps and constraint iterations to ensure that all constraints are satisfied.

The "glue" between voxels is implemented by sharing particles between adjacent voxels. When the strain or shear on a face exceeds a certain threshold, the "glue" is broken by detaching the shared particles. This allows the voxels to break apart and form clusters.

The shape constraints are what give the voxel clusters their rigidity. The `solve_voxel_shape` function tries to preserve the original shape of each voxel, which in turn preserves the shape of the cluster.

In summary, the `simulate_voxel_pbd` function uses a combination of PBD, shape matching, and a "glue" constraint to simulate the behavior of voxel clusters. The PBD solver predicts the next position of each particle, the collision detection and response system prevents particles from interpenetrating, and the shape-matching and "glue" constraints preserve the shape and connectivity of the voxel clusters.
