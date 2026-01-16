# Voxel and Particle Conservation on Break

The question of whether the number of voxels and particles is conserved when a voxel breaks is a good one, and the answer reveals a key aspect of the simulation's design.

### Voxel Conservation: Yes

The number of **voxels is conserved** during the breaking process.

When a voxel face exceeds the strain or shear thresholds, `break_face_link()` is called directly to detach that face. These functions only deal with the particles at the corners of the voxel's faces. They do not create or destroy any `Voxel` structs. The `voxel_count` remains unchanged throughout this process. A voxel is only ever added with `addVoxel()` and there is no corresponding `removeVoxel()` function called during a break (destruction of voxels happens in other parts of the code, for example when a bullet hits a voxel).

### Particle Conservation: No

The number of **particles is not conserved** when a voxel breaks. In fact, the number of particles increases.

Here's why:

1.  **Shared Particles:** Initially, adjacent voxels are "glued" together by sharing the same `Particle` pointers at their joining faces. A single particle can be a corner for up to 8 different voxels, and its `refcount` (reference count) will reflect this.

2.  **Detaching Faces:** When a voxel face breaks, the simulation needs to allow the two previously joined faces to move apart. To do this, the shared particles along that face must be duplicated.

3.  **`detach_face_particles()` and `particle_clone()`:** This is handled by the `detach_face_particles()` function. It iterates through the four particles of the face being broken. If a particle's `refcount` is greater than 1, it means the particle is shared. The function then calls `particle_clone()` to create a brand new, identical particle.

4.  **New Becomes Old:** The breaking voxel then updates its `particles` array to point to this new particle. The old, shared particle has its `refcount` decremented via `particle_release()`.

The end result is that where there was once one shared particle, there are now two independent particles, one for each of the separated voxels. This is what allows the voxels to fly apart as separate objects.

### Conclusion

In short, the simulation conserves the number of voxels but not the number of particles when a break occurs. The creation of new particles is the fundamental mechanism that allows for the "un-gluing" and separation of voxel clusters.
