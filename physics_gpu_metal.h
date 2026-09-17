#ifndef FPS_PHYSICS_GPU_METAL_H
#define FPS_PHYSICS_GPU_METAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "physics_gpu_layout.h"

typedef struct FpsGpuUniforms {
    int32_t particle_count;
    int32_t sim_count;
    int32_t voxel_count;
    int32_t static_collider_count;
    int32_t hash_size;
    int32_t static_hash_size;
    int32_t active_players;
    int32_t break_damp_frames;
    int32_t mode;
    int32_t vgs_color;
    int32_t integer_padding[2];
    float dt;
    float voxel_size;
    float floor_size;
    float gravity;
    float velocity_damping;
    float sor;
    float collision_relaxation;
    float vgs_alpha;
    float vgs_beta;
    float vgs_epsilon;
    float strain_threshold;
    float shear_threshold;
    float tether_spring;
    float tether_damping;
    float rest_grid_step;
    float padding;
    float players[4][4];
    float tether_targets[4][4];
} FpsGpuUniforms;

_Static_assert(sizeof(FpsGpuUniforms) == 240, "Metal uniform layout mismatch");

enum { FPS_METAL_BUILD, FPS_METAL_SHAPE, FPS_METAL_INTERFACE, FPS_METAL_CONTACT, FPS_METAL_FRACTURE, FPS_METAL_STAGES };
typedef struct FpsMetalProfileInfo {
    double last_gpu_exec_ms;
    double last_wait_ms;
    double last_encode_ms;
    int last_dispatch_count;
    bool stage_counters_available;
    uint64_t stage_ticks[FPS_METAL_STAGES];
} FpsMetalProfileInfo;

#if defined(FPS_GPU_METAL)
bool fps_metal_initialize(const char *library_path, long long *max_buffer_size,
                          char *error, size_t error_capacity);
void fps_metal_shutdown(void);
void *fps_metal_buffer_create(size_t size);
void fps_metal_buffer_destroy(void *buffer);
bool fps_metal_buffer_update(void *buffer, const void *data, size_t size, size_t offset);
bool fps_metal_buffer_read(void *buffer, void *data, size_t size, size_t offset);
void fps_metal_bind_buffer(int slot, void *buffer);
void fps_metal_get_uniforms(FpsGpuUniforms *uniforms);
void fps_metal_set_uniforms(const FpsGpuUniforms *uniforms);
void fps_metal_set_vgs_color(int color);
void fps_metal_profile_mark(unsigned stage,bool end);
void *fps_metal_device(void);
void *fps_metal_encoder(void);
void fps_metal_restore_bindings(void);
bool fps_metal_begin_batch(void);
bool fps_metal_dispatch(int mode, int count);
bool fps_metal_dispatch_indirect(int mode, void *buffer, size_t offset);
bool fps_metal_end_batch(void);
bool fps_metal_commit_batch(bool wait);
bool fps_metal_wait_batch(void);
void fps_metal_get_profile_info(FpsMetalProfileInfo *out);
#endif

#endif
