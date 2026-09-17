#ifndef FPS_ADAPTIVE_PHYSICS_H
#define FPS_ADAPTIVE_PHYSICS_H
#include "adaptive_metal.h"
#include "physics_gpu_metal.h"
typedef struct {
    int32_t xyz[3]; uint32_t arena;
    uint32_t incident[8];
    float inverse_mass, base_inverse_mass;
    uint32_t count, simulated;
} AdaptiveParticleInput;
_Static_assert(sizeof(AdaptiveParticleInput)==64,"adaptive particle input ABI");
typedef struct AdaptivePhysics AdaptivePhysics;
AdaptivePhysics *adaptive_physics_create(void *device,const AdaptiveInput *,uint32_t cells,
    const AdaptiveParticleInput *,uint32_t particles,char *error,size_t error_size);
void adaptive_physics_destroy(AdaptivePhysics *);
bool adaptive_physics_build(AdaptivePhysics *,void *encoder,void *particle_buffer,void *voxel_buffer,
                           void *sim_buffer,void *control_buffer,const FpsGpuUniforms *);
bool adaptive_physics_stage(AdaptivePhysics *,void *encoder,void *particle_buffer,void *voxel_buffer,
                           void *sim_buffer,void *control_buffer,const FpsGpuUniforms *,uint32_t mode);
enum { ADAPTIVE_STAGE_SHAPE=30, ADAPTIVE_STAGE_ATTACH=31, ADAPTIVE_STAGE_REFRESH=32,
       ADAPTIVE_STAGE_FRACTURE=33 };
bool adaptive_physics_independent(const AdaptivePhysics *,uint32_t particle);
double adaptive_physics_interface_error(const AdaptivePhysics *,void *particle_buffer);
double adaptive_physics_mass_error(const AdaptivePhysics *,void *particle_buffer);
uint64_t adaptive_physics_dispatches(const AdaptivePhysics *,bool include_gpu_indirect);
size_t adaptive_physics_bytes(const AdaptivePhysics *);
void adaptive_physics_diagnostics(const AdaptivePhysics *,uint32_t *leaves,uint32_t *generation,uint32_t *error);
#endif
