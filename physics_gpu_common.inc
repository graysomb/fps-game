/* Included by fps_ray.c after the voxel/PBD helpers are defined. */
#include "physics_gpu_metal.h"

typedef struct { float pos_radius[4], prev_inv_mass[4], predicted_base_inv_mass[4], velocity[4]; } GpuParticleState;
typedef struct { uint32_t bits[4]; } GpuCorrectionState;
typedef struct { int32_t value[4]; } GpuCellState;
typedef struct { int32_t value[4]; } GpuCollisionMeta;
typedef int32_t GpuTetherState;
typedef struct { int32_t value[4]; } GpuStaticCell;
typedef struct {
    float center[4];
    float bounds_min[4];
    float bounds_max[4];
} GpuStaticCollider;
typedef struct { int32_t particle_count, sim_count, particle_capacity, overflow; } GpuControlState;
typedef struct { int32_t value[4]; } GpuTopologyState;
typedef struct {
    uint32_t particle[8];
    float pos_rest_edge[4];
    float velocity_rest_volume[4];
    int32_t flags[4];
    float bounds_min[4];
    float bounds_max[4];
    uint32_t lifecycle[4];
} GpuVoxelState;

_Static_assert(sizeof(GpuParticleState) == 64, "GpuParticleState/std430 mismatch");
_Static_assert(sizeof(GpuCorrectionState) == 16, "GpuCorrectionState/std430 mismatch");
_Static_assert(sizeof(GpuCellState) == 16, "GpuCellState/std430 mismatch");
_Static_assert(sizeof(GpuCollisionMeta) == 16, "GpuCollisionMeta/std430 mismatch");
_Static_assert(sizeof(GpuVoxelState) == 128, "GpuVoxelState/std430 mismatch");
_Static_assert(sizeof(GpuStaticCell) == 16, "GpuStaticCell/std430 mismatch");
_Static_assert(sizeof(GpuStaticCollider) == 48, "GpuStaticCollider/std430 mismatch");
_Static_assert(sizeof(GpuControlState) == 16, "GpuControlState/std430 mismatch");

typedef struct {
    uintptr_t id;
    size_t capacity;
} GpuBuffer;

enum {
    GPU_LOC_MODE = 0,
    GPU_LOC_PARTICLE_COUNT,
    GPU_LOC_SIM_COUNT,
    GPU_LOC_VOXEL_COUNT,
    GPU_LOC_STATIC_COLLIDER_COUNT,
    GPU_LOC_HASH_SIZE,
    GPU_LOC_STATIC_HASH_SIZE,
    GPU_LOC_ACTIVE_PLAYERS,
    GPU_LOC_DT,
    GPU_LOC_VOXEL_SIZE,
    GPU_LOC_FLOOR_SIZE,
    GPU_LOC_GRAVITY,
    GPU_LOC_VELOCITY_DAMPING,
    GPU_LOC_SOR,
    GPU_LOC_COLLISION_RELAXATION,
    GPU_LOC_VGS_ALPHA,
    GPU_LOC_VGS_BETA,
    GPU_LOC_VGS_EPSILON,
    GPU_LOC_STRAIN_THRESHOLD,
    GPU_LOC_SHEAR_THRESHOLD,
    GPU_LOC_TETHER_SPRING,
    GPU_LOC_TETHER_DAMPING,
    GPU_LOC_BREAK_DAMP_FRAMES,
    GPU_LOC_REST_GRID_STEP,
    GPU_LOC_PLAYERS,
    GPU_LOC_TETHER_TARGETS,
    GPU_LOC_COUNT
};

typedef struct {
    bool ready;
    bool resident_valid;
    unsigned int program;
    int loc[GPU_LOC_COUNT];
    long long max_ssbo_size;
    GpuBuffer particle_buffer;
    GpuBuffer correction_buffer;
    GpuBuffer cell_buffer;
    GpuBuffer sim_id_buffer;
    GpuBuffer voxel_buffer;
    GpuBuffer hash_head_buffer;
    GpuBuffer hash_next_buffer;
    GpuBuffer static_cell_buffer;
    GpuBuffer static_collider_buffer;
    GpuBuffer refcount_buffer;
    GpuBuffer control_buffer;
    GpuBuffer clone_parent_buffer;
    GpuBuffer dispatch_buffer;
    GpuBuffer topology_buffer;
    GpuBuffer tether_buffer;
    GpuBuffer collision_meta_buffer;
    GpuParticleState *particle_stage;
    GpuCellState *cell_stage;
    GpuTetherState *tether_stage;
    GpuTetherState *tether_compare_stage;
    GpuCollisionMeta *collision_meta_stage;
    Particle **compact_to_particle;
    uint32_t *sim_ids;
    uint32_t *collision_ids;
    GpuVoxelState *voxel_stage;
    GpuStaticCell *static_stage;
    GpuStaticCollider *static_collider_stage;
    int *compact_to_voxel;
    int32_t *refcount_stage;
    int32_t *clone_parent_stage;
    GpuControlState control_stage;
    GpuTopologyState *topology_stage;
    int *world_to_dynamic;
    int *pool_to_compact;
    size_t particle_capacity;
    size_t sim_capacity;
    size_t voxel_capacity;
    size_t static_collider_capacity;
    uint64_t static_uploaded_generation;
    uint64_t collision_uploaded_generation;
    int compact_count;
    int uploaded_particle_count;
    int particle_slot_capacity;
    int sim_slot_capacity;
    int dynamic_voxel_count;
    int collision_count;
    int static_collider_count;
    uint64_t upload_bytes;
    uint64_t readback_bytes;
    uint64_t sync_count;
    uint64_t full_uploads;
    uint64_t resident_reuses;
    uint64_t static_uploads;
    uint64_t clones_created;
    double pack_ms;
    double upload_ms;
    double dispatch_ms;
    double readback_ms;
    double commit_ms;
} GpuPhysicsState;

static GpuPhysicsState gpuPhysics = { 0 };

static bool gpu_resize_array(void **ptr, size_t count, size_t item_size) {
    if (count == 0) count = 1;
    if (item_size != 0 && count > SIZE_MAX / item_size) return false;
    void *next = realloc(*ptr, count * item_size);
    if (!next) return false;
    *ptr = next;
    return true;
}

static void gpu_buffer_destroy(GpuBuffer *buffer) {
    if (buffer->id != 0) {
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
        fps_metal_buffer_destroy((void *)buffer->id);
#else
        rlUnloadShaderBuffer((unsigned int)buffer->id);
#endif
    }
    buffer->id = 0;
    buffer->capacity = 0;
}

static bool gpu_buffer_ensure(GpuBuffer *buffer, size_t size) {
    if (size == 0) size = 16;
    if (size > UINT_MAX || (gpuPhysics.max_ssbo_size > 0 && (long long)size > gpuPhysics.max_ssbo_size)) {
        return false;
    }
    if (buffer->id != 0 && buffer->capacity >= size) return true;
    size_t capacity = 256;
    while (capacity < size) {
        if (capacity > SIZE_MAX / 2) return false;
        capacity *= 2;
    }
    if (gpuPhysics.max_ssbo_size > 0 && (long long)capacity > gpuPhysics.max_ssbo_size) return false;
    gpu_buffer_destroy(buffer);
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    buffer->id = (uintptr_t)fps_metal_buffer_create(capacity);
#else
    buffer->id = (uintptr_t)rlLoadShaderBuffer((unsigned int)capacity, NULL, RL_DYNAMIC_COPY);
#endif
    if (buffer->id == 0) return false;
    buffer->capacity = capacity;
    return true;
}

static void gpu_buffer_update(GpuBuffer *buffer, const void *data, size_t size, size_t offset) {
    if (!buffer || !buffer->id || !data || size == 0) return;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    (void)fps_metal_buffer_update((void *)buffer->id, data, size, offset);
#else
    rlUpdateShaderBuffer((unsigned int)buffer->id, data, (unsigned int)size, (unsigned int)offset);
#endif
}

static void gpu_buffer_read(GpuBuffer *buffer, void *data, size_t size, size_t offset) {
    if (!buffer || !buffer->id || !data || size == 0) return;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    (void)fps_metal_buffer_read((void *)buffer->id, data, size, offset);
#else
    rlReadShaderBuffer((unsigned int)buffer->id, data, (unsigned int)size, (unsigned int)offset);
#endif
}

static void gpu_physics_release_resources(void) {
    gpu_buffer_destroy(&gpuPhysics.particle_buffer);
    gpu_buffer_destroy(&gpuPhysics.correction_buffer);
    gpu_buffer_destroy(&gpuPhysics.cell_buffer);
    gpu_buffer_destroy(&gpuPhysics.sim_id_buffer);
    gpu_buffer_destroy(&gpuPhysics.voxel_buffer);
    gpu_buffer_destroy(&gpuPhysics.hash_head_buffer);
    gpu_buffer_destroy(&gpuPhysics.hash_next_buffer);
    gpu_buffer_destroy(&gpuPhysics.static_cell_buffer);
    gpu_buffer_destroy(&gpuPhysics.static_collider_buffer);
    gpu_buffer_destroy(&gpuPhysics.refcount_buffer);
    gpu_buffer_destroy(&gpuPhysics.control_buffer);
    gpu_buffer_destroy(&gpuPhysics.clone_parent_buffer);
    gpu_buffer_destroy(&gpuPhysics.dispatch_buffer);
    gpu_buffer_destroy(&gpuPhysics.topology_buffer);
    gpu_buffer_destroy(&gpuPhysics.tether_buffer);
    gpu_buffer_destroy(&gpuPhysics.collision_meta_buffer);
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    fps_metal_shutdown();
#else
    if (gpuPhysics.program != 0) rlUnloadShaderProgram(gpuPhysics.program);
#endif
    gpuPhysics.program = 0;
    gpuPhysics.ready = false;
    gpuPhysics.resident_valid = false;
}

static void gpu_physics_shutdown(void) {
    if (gpuTransferStatsEnabled && (gpuPhysics.full_uploads || gpuPhysics.sync_count)) {
        fprintf(stderr,
                "gpu-transfer mode=%s uploads=%llu reuses=%llu staticUploads=%llu syncs=%llu clones=%llu uploadBytes=%llu readbackBytes=%llu packMs=%.3f uploadMs=%.3f dispatchMs=%.3f readbackMs=%.3f commitMs=%.3f\n",
                gpuTransferMode == GPU_TRANSFER_RESIDENT ? "resident" : "legacy",
                (unsigned long long)gpuPhysics.full_uploads,
                (unsigned long long)gpuPhysics.resident_reuses,
                (unsigned long long)gpuPhysics.static_uploads,
                (unsigned long long)gpuPhysics.sync_count,
                (unsigned long long)gpuPhysics.clones_created,
                (unsigned long long)gpuPhysics.upload_bytes,
                (unsigned long long)gpuPhysics.readback_bytes,
                gpuPhysics.pack_ms, gpuPhysics.upload_ms, gpuPhysics.dispatch_ms,
                gpuPhysics.readback_ms, gpuPhysics.commit_ms);
    }
    gpu_physics_release_resources();
    free(gpuPhysics.particle_stage);
    free(gpuPhysics.cell_stage);
    free(gpuPhysics.tether_stage);
    free(gpuPhysics.tether_compare_stage);
    free(gpuPhysics.collision_meta_stage);
    free(gpuPhysics.compact_to_particle);
    free(gpuPhysics.sim_ids);
    free(gpuPhysics.collision_ids);
    free(gpuPhysics.voxel_stage);
    free(gpuPhysics.static_stage);
    free(gpuPhysics.static_collider_stage);
    free(gpuPhysics.compact_to_voxel);
    free(gpuPhysics.refcount_stage);
    free(gpuPhysics.clone_parent_stage);
    free(gpuPhysics.topology_stage);
    free(gpuPhysics.world_to_dynamic);
    free(gpuPhysics.pool_to_compact);
    memset(&gpuPhysics, 0, sizeof(gpuPhysics));
}

static void gpu_fail(const char *reason) {
    snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "%s", reason ? reason : "GPU failure");
    physicsBackend.sticky_fallback = true;
    physicsBackend.gpu_available = false;
    gpu_physics_release_resources();
}

static bool gpu_stage_reserve(size_t particles, size_t sims, size_t voxels_count) {
    if (particles > gpuPhysics.particle_capacity) {
        size_t cap = gpuPhysics.particle_capacity ? gpuPhysics.particle_capacity : 256;
        while (cap < particles) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.particle_stage, cap, sizeof(*gpuPhysics.particle_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.cell_stage, cap, sizeof(*gpuPhysics.cell_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.tether_stage, cap, sizeof(*gpuPhysics.tether_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.tether_compare_stage, cap, sizeof(*gpuPhysics.tether_compare_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.collision_meta_stage, cap, sizeof(*gpuPhysics.collision_meta_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.compact_to_particle, cap, sizeof(*gpuPhysics.compact_to_particle)) ||
            !gpu_resize_array((void **)&gpuPhysics.refcount_stage, cap, sizeof(*gpuPhysics.refcount_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.clone_parent_stage, cap, sizeof(*gpuPhysics.clone_parent_stage))) return false;
        gpuPhysics.particle_capacity = cap;
    }
    if (sims > gpuPhysics.sim_capacity) {
        size_t cap = gpuPhysics.sim_capacity ? gpuPhysics.sim_capacity : 256;
        while (cap < sims) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.sim_ids, cap, sizeof(*gpuPhysics.sim_ids)) ||
            !gpu_resize_array((void **)&gpuPhysics.collision_ids, cap, sizeof(*gpuPhysics.collision_ids))) return false;
        gpuPhysics.sim_capacity = cap;
    }
    if (voxels_count > gpuPhysics.voxel_capacity) {
        size_t cap = gpuPhysics.voxel_capacity ? gpuPhysics.voxel_capacity : 256;
        while (cap < voxels_count) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.voxel_stage, cap, sizeof(*gpuPhysics.voxel_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.compact_to_voxel, cap, sizeof(*gpuPhysics.compact_to_voxel)) ||
            !gpu_resize_array((void **)&gpuPhysics.topology_stage, cap * 2u, sizeof(*gpuPhysics.topology_stage))) return false;
        gpuPhysics.voxel_capacity = cap;
    }
    return true;
}

static int gpu_add_particle(Particle *p) {
    if (!p) return -1;
    ptrdiff_t pool_id = p - particles_pool;
    if (pool_id < 0 || pool_id >= particle_pool_count) return -1;
    int mapped = gpuPhysics.pool_to_compact[pool_id];
    if (mapped >= 0) return mapped;
    int id = gpuPhysics.compact_count++;
    gpuPhysics.pool_to_compact[pool_id] = id;
    gpuPhysics.compact_to_particle[id] = p;
    GpuParticleState *g = &gpuPhysics.particle_stage[id];
    g->pos_radius[0] = p->pos.x; g->pos_radius[1] = p->pos.y; g->pos_radius[2] = p->pos.z; g->pos_radius[3] = p->radius;
    g->prev_inv_mass[0] = p->prev_pos.x; g->prev_inv_mass[1] = p->prev_pos.y; g->prev_inv_mass[2] = p->prev_pos.z; g->prev_inv_mass[3] = p->inv_mass;
    g->predicted_base_inv_mass[0] = p->predicted_pos.x; g->predicted_base_inv_mass[1] = p->predicted_pos.y;
    g->predicted_base_inv_mass[2] = p->predicted_pos.z; g->predicted_base_inv_mass[3] = p->base_inv_mass;
    g->velocity[0] = p->vel.x; g->velocity[1] = p->vel.y; g->velocity[2] = p->vel.z; g->velocity[3] = 0.0f;
    gpuPhysics.cell_stage[id].value[0] = p->cell_x;
    gpuPhysics.cell_stage[id].value[1] = p->cell_y;
    gpuPhysics.cell_stage[id].value[2] = p->cell_z;
    gpuPhysics.cell_stage[id].value[3] = p->break_timer;
    gpuPhysics.collision_meta_stage[id].value[0] = p->rest_cell_x;
    gpuPhysics.collision_meta_stage[id].value[1] = p->rest_cell_y;
    gpuPhysics.collision_meta_stage[id].value[2] = p->rest_cell_z;
    gpuPhysics.collision_meta_stage[id].value[3] = p->collision_group;
    gpuPhysics.tether_stage[id] = -1;
    return id;
}

static uint32_t gpu_glued_face_mask(const Voxel *v) {
    uint32_t mask = 0;
    for (int face = 0; face < 6; ++face) if (v->glued_faces[face]) mask |= 1u << face;
    return mask;
}

static bool gpu_pack_world(void) {
    int dynamic_count = 0;
    for (int i = 0; i < voxel_count; ++i) if (voxels[i].simulate) dynamic_count++;
    size_t maximum_references = (size_t)sim_particle_count + (size_t)dynamic_count * 8u;
    if (maximum_references > (size_t)particle_pool_count) maximum_references = (size_t)particle_pool_count;
    size_t cpu_free_slots = (size_t)free_particle_count + (size_t)(MAX_PARTICLES - particle_pool_count);
    size_t possible_clones = (size_t)dynamic_count * 8u;
    if (possible_clones > cpu_free_slots) possible_clones = cpu_free_slots;
    size_t slot_capacity = maximum_references + possible_clones;
    if (slot_capacity > (size_t)MAX_PARTICLES) slot_capacity = MAX_PARTICLES;
    if (!gpu_stage_reserve(slot_capacity, slot_capacity, (size_t)dynamic_count)) return false;
    if (!gpuPhysics.world_to_dynamic) {
        gpuPhysics.world_to_dynamic = malloc((size_t)MAX_VOXELS * sizeof(*gpuPhysics.world_to_dynamic));
        if (!gpuPhysics.world_to_dynamic) return false;
    }
    if (voxel_count > 0) memset(gpuPhysics.world_to_dynamic, 0xff,
                                (size_t)voxel_count * sizeof(*gpuPhysics.world_to_dynamic));
    if (!gpuPhysics.pool_to_compact) {
        gpuPhysics.pool_to_compact = malloc((size_t)MAX_PARTICLES * sizeof(*gpuPhysics.pool_to_compact));
        if (!gpuPhysics.pool_to_compact) return false;
    }
    if (particle_pool_count > 0) memset(gpuPhysics.pool_to_compact, 0xff, (size_t)particle_pool_count * sizeof(*gpuPhysics.pool_to_compact));
    gpuPhysics.compact_count = 0;

    for (int i = 0; i < sim_particle_count; ++i) {
        int id = gpu_add_particle(sim_particles[i]);
        if (id < 0) return false;
        gpuPhysics.sim_ids[i] = (uint32_t)id;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        for (int c = 0; c < 8; ++c) if (gpu_add_particle(v->particles[c]) < 0) return false;
    }

    gpuPhysics.uploaded_particle_count = gpuPhysics.compact_count;
    gpuPhysics.particle_slot_capacity = (int)slot_capacity;
    gpuPhysics.sim_slot_capacity = (int)slot_capacity;
    gpuPhysics.collision_count = 0;
    for (int i = 0; i < collision_particle_count; ++i) {
        Particle *particle = collision_particles[i];
        if (!particle) continue;
        ptrdiff_t pool_id = particle - particles_pool;
        if (pool_id < 0 || pool_id >= particle_pool_count) continue;
        int compact_id = gpuPhysics.pool_to_compact[pool_id];
        if (compact_id < 0 || compact_id >= gpuPhysics.compact_count) continue;
        gpuPhysics.collision_ids[gpuPhysics.collision_count++] = (uint32_t)compact_id;
    }
    for (int i = 0; i < gpuPhysics.compact_count; ++i) {
        Particle *p = gpuPhysics.compact_to_particle[i];
        gpuPhysics.refcount_stage[i] = p ? p->refcount : 0;
        gpuPhysics.clone_parent_stage[i] = INT_MAX;
    }
    for (int i = gpuPhysics.compact_count; i < gpuPhysics.particle_slot_capacity; ++i) {
        gpuPhysics.compact_to_particle[i] = NULL;
        gpuPhysics.refcount_stage[i] = 0;
        gpuPhysics.clone_parent_stage[i] = INT_MAX;
    }
    for (int i = 0; i < gpuPhysics.collision_count; ++i) {
        uint32_t id = gpuPhysics.collision_ids[i];
        if (id < (uint32_t)gpuPhysics.particle_slot_capacity)
            gpuPhysics.clone_parent_stage[id] |= INT_MIN;
    }
    gpuPhysics.control_stage = (GpuControlState){
        gpuPhysics.compact_count, sim_particle_count, gpuPhysics.particle_slot_capacity, 0
    };
    if (getenv("FPS_COLLISION_DIAGNOSTICS")) {
        fprintf(stderr, "collision-pack sim=%d collision=%d particles=%d capacity=%d\n",
                sim_particle_count, gpuPhysics.collision_count,
                gpuPhysics.compact_count, gpuPhysics.sim_slot_capacity);
    }

    gpuPhysics.dynamic_voxel_count = 0;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        int dynamic_id = gpuPhysics.dynamic_voxel_count++;
        gpuPhysics.compact_to_voxel[dynamic_id] = i;
        gpuPhysics.world_to_dynamic[i] = dynamic_id;
        GpuVoxelState *g = &gpuPhysics.voxel_stage[dynamic_id];
        memset(g, 0, sizeof(*g));
        for (int c = 0; c < 8; ++c) {
            g->particle[c] = UINT32_MAX;
            if (v->particles[c]) {
                ptrdiff_t pool_id = v->particles[c] - particles_pool;
                if (pool_id >= 0 && pool_id < particle_pool_count) {
                    int mapped = gpuPhysics.pool_to_compact[pool_id];
                    if (mapped >= 0) g->particle[c] = (uint32_t)mapped;
                }
            }
        }
        g->pos_rest_edge[0] = v->pos.x; g->pos_rest_edge[1] = v->pos.y; g->pos_rest_edge[2] = v->pos.z; g->pos_rest_edge[3] = v->rest_edge;
        g->velocity_rest_volume[0] = v->vel.x; g->velocity_rest_volume[1] = v->vel.y;
        g->velocity_rest_volume[2] = v->vel.z; g->velocity_rest_volume[3] = v->rest_volume;
        // flags.x means awake dynamic. Sleeping dynamics remain packed as
        // zero-mass collision anchors but skip voxel-level simulation kernels.
        g->flags[0] = v->sleeping ? 0 : 1; g->flags[1] = v->isBullet ? 1 : 0;
        g->flags[2] = v->type; g->flags[3] = v->simulate_dofs ? 1 : 0;
        VoxelWorldBounds bounds;
        voxel_world_bounds(v, &bounds);
        g->bounds_min[0] = bounds.minx; g->bounds_min[1] = bounds.miny; g->bounds_min[2] = bounds.minz;
        g->bounds_max[0] = bounds.maxx; g->bounds_max[1] = bounds.maxy; g->bounds_max[2] = bounds.maxz;
        g->bounds_min[3] = (float)v->wake_timer;
        g->bounds_max[3] = (float)v->wake_timer;
        g->lifecycle[0] = (uint32_t)((v->skipCollisionVelocityFrames > 0) ? v->skipCollisionVelocityFrames : 0);
        g->lifecycle[1] = v->wake_source ? 1u : 0u;
        g->lifecycle[2] = gpu_glued_face_mask(v) | (v->glueEligible ? 0x80000000u : 0u);
        g->lifecycle[3] = v->break_mask;
    }

    for (int dynamic_id = 0; dynamic_id < gpuPhysics.dynamic_voxel_count; ++dynamic_id) {
        int world_id = gpuPhysics.compact_to_voxel[dynamic_id];
        Voxel *v = &voxels[world_id];
        GpuTopologyState *a = &gpuPhysics.topology_stage[dynamic_id * 2];
        GpuTopologyState *b = &gpuPhysics.topology_stage[dynamic_id * 2 + 1];
        for (int face = 0; face < 6; ++face) {
            int neighbor_world = table_get(v->gx + face_offsets[face][0],
                                           v->gy + face_offsets[face][1],
                                           v->gz + face_offsets[face][2]);
            int neighbor_dynamic = (neighbor_world >= 0 && neighbor_world < voxel_count)
                ? gpuPhysics.world_to_dynamic[neighbor_world] : -1;
            if (face < 4) a->value[face] = neighbor_dynamic;
            else b->value[face - 4] = neighbor_dynamic;
        }
        b->value[2] = dynamic_id;
        b->value[3] = world_id;
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate || v->type != 0 || v->isBullet || tetherTag[i] <= 0) continue;
        int player = tetherTag[i] - 1;
        if (player < 0 || player >= MAX_PLAYERS) continue;
        for (int c = 0; c < 8; ++c) {
            ptrdiff_t pool_id = v->particles[c] - particles_pool;
            if (pool_id < 0 || pool_id >= particle_pool_count) continue;
            int id = gpuPhysics.pool_to_compact[pool_id];
            if (id < 0 || gpuPhysics.tether_stage[id] >= 0) continue;
            gpuPhysics.tether_stage[id] = player;
        }
    }
    return true;
}

static uint32_t gpu_static_hash_coord(int x, int y, int z) {
    uint32_t h = (uint32_t)x * 73856093u;
    h ^= (uint32_t)y * 19349663u;
    h ^= (uint32_t)z * 83492791u;
    return h & (HASH_SIZE - 1u);
}

static float gpu_int_bits_as_float(int value) {
    uint32_t bits = (uint32_t)value;
    float result;
    memcpy(&result, &bits, sizeof(result));
    return result;
}

static void gpu_fill_static_patch(GpuStaticCollider *collider, const Patch *patch,
                                  int patch_index, int next) {
    memset(collider, 0, sizeof(*collider));
    float plane = (float)(patch->layer + (patch->positive ? 1 : 0)) * VOXEL_SIZE;
    float min_i = (float)patch->i0 * VOXEL_SIZE;
    float max_i = (float)(patch->i0 + patch->di) * VOXEL_SIZE;
    float min_j = (float)patch->j0 * VOXEL_SIZE;
    float max_j = (float)(patch->j0 + patch->dj) * VOXEL_SIZE;
    if (patch->plane == 0) {
        collider->center[2] = patch->positive ? 1.0f : -1.0f;
        collider->bounds_min[0] = min_i; collider->bounds_max[0] = max_i;
        collider->bounds_min[1] = min_j; collider->bounds_max[1] = max_j;
        collider->bounds_min[2] = collider->bounds_max[2] = plane;
    } else if (patch->plane == 1) {
        collider->center[1] = patch->positive ? 1.0f : -1.0f;
        collider->bounds_min[0] = min_i; collider->bounds_max[0] = max_i;
        collider->bounds_min[2] = min_j; collider->bounds_max[2] = max_j;
        collider->bounds_min[1] = collider->bounds_max[1] = plane;
    } else {
        collider->center[0] = patch->positive ? 1.0f : -1.0f;
        collider->bounds_min[1] = min_i; collider->bounds_max[1] = max_i;
        collider->bounds_min[2] = min_j; collider->bounds_max[2] = max_j;
        collider->bounds_min[0] = collider->bounds_max[0] = plane;
    }
    collider->center[3] = gpu_int_bits_as_float(next);
    collider->bounds_min[3] = gpu_int_bits_as_float(patch_index);
}

static bool gpu_rebuild_static_cells(void) {
    bool use_surface_cache = getenv("FPS_COLLISION_STATIC_LEGACY") == NULL;
    if (use_surface_cache && !ensure_static_surface_cache()) return false;
    if (!gpuPhysics.static_stage) {
        gpuPhysics.static_stage = malloc((size_t)HASH_SIZE * sizeof(*gpuPhysics.static_stage));
        if (!gpuPhysics.static_stage) return false;
    }
    for (int i = 0; i < HASH_SIZE; ++i) {
        gpuPhysics.static_stage[i].value[0] = 0;
        gpuPhysics.static_stage[i].value[1] = 0;
        gpuPhysics.static_stage[i].value[2] = 0;
        gpuPhysics.static_stage[i].value[3] = -1;
    }
    int static_voxel_count = 0;
    for (int i = 0; i < voxel_count; ++i) if (!voxels[i].simulate) static_voxel_count++;
    size_t required = (use_surface_cache ? (size_t)staticSurfacePatchRefCount : 0u) +
                      (size_t)static_voxel_count;
    if (required > gpuPhysics.static_collider_capacity) {
        size_t cap = gpuPhysics.static_collider_capacity ? gpuPhysics.static_collider_capacity : 256;
        while (cap < required) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.static_collider_stage, cap,
                              sizeof(*gpuPhysics.static_collider_stage))) return false;
        gpuPhysics.static_collider_capacity = cap;
    }
    gpuPhysics.static_collider_count = 0;

    for (int cell_index = 0; use_surface_cache && cell_index < staticSurfaceCellCount; ++cell_index) {
        const StaticSurfaceCell *surface_cell = &staticSurfaceCells[cell_index];
        int head = -1;
        for (int item = 0; item < surface_cell->count; ++item) {
            int patch_index = staticSurfacePatchRefs[surface_cell->start + item];
            if (patch_index < 0 || patch_index >= patchCount) continue;
            int collider_id = gpuPhysics.static_collider_count++;
            gpu_fill_static_patch(&gpuPhysics.static_collider_stage[collider_id],
                                  &patches[patch_index], patch_index, head);
            head = collider_id;
        }
        uint32_t h = gpu_static_hash_coord(surface_cell->x, surface_cell->y, surface_cell->z);
        while (gpuPhysics.static_stage[h].value[3] != -1) h = (h + 1u) & (HASH_SIZE - 1u);
        gpuPhysics.static_stage[h].value[0] = surface_cell->x;
        gpuPhysics.static_stage[h].value[1] = surface_cell->y;
        gpuPhysics.static_stage[h].value[2] = surface_cell->z;
        gpuPhysics.static_stage[h].value[3] = head;
    }

    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *v = &voxels[i];
        if (v->simulate) continue;
        int collider_id = gpuPhysics.static_collider_count++;
        GpuStaticCollider *collider = &gpuPhysics.static_collider_stage[collider_id];
        VoxelWorldBounds bounds;
        voxel_world_bounds(v, &bounds);
        collider->center[0] = v->pos.x; collider->center[1] = v->pos.y;
        collider->center[2] = v->pos.z; collider->center[3] = 0.0f;
        collider->bounds_min[0] = bounds.minx; collider->bounds_min[1] = bounds.miny;
        collider->bounds_min[2] = bounds.minz; collider->bounds_min[3] = 0.0f;
        collider->bounds_max[0] = bounds.maxx; collider->bounds_max[1] = bounds.maxy;
        collider->bounds_max[2] = bounds.maxz; collider->bounds_max[3] = 0.0f;
        for (int z = v->min_gz; z <= v->max_gz; ++z) {
            for (int y = v->min_gy; y <= v->max_gy; ++y) {
                for (int x = v->min_gx; x <= v->max_gx; ++x) {
                    uint32_t h = gpu_static_hash_coord(x, y, z);
                    for (int probe = 0; probe < HASH_SIZE; ++probe) {
                        GpuStaticCell *entry = &gpuPhysics.static_stage[h];
                        if (entry->value[3] == -1 ||
                            (entry->value[0] == x && entry->value[1] == y && entry->value[2] == z)) {
                            entry->value[0] = x; entry->value[1] = y; entry->value[2] = z;
                            entry->value[3] = -(collider_id + 2);
                            break;
                        }
                        h = (h + 1u) & (HASH_SIZE - 1u);
                    }
                }
            }
        }
    }
    return true;
}

static void gpu_bind_world_buffers(void) {
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    GpuBuffer *buffers[FPS_GPU_BUFFER_COUNT] = {
        &gpuPhysics.particle_buffer, &gpuPhysics.correction_buffer, &gpuPhysics.cell_buffer,
        &gpuPhysics.sim_id_buffer, &gpuPhysics.voxel_buffer, &gpuPhysics.hash_head_buffer,
        &gpuPhysics.hash_next_buffer, &gpuPhysics.static_cell_buffer, &gpuPhysics.tether_buffer,
        &gpuPhysics.collision_meta_buffer, &gpuPhysics.static_collider_buffer,
        &gpuPhysics.refcount_buffer, &gpuPhysics.control_buffer, &gpuPhysics.clone_parent_buffer,
        &gpuPhysics.dispatch_buffer, &gpuPhysics.topology_buffer
    };
    for (int i = 0; i < FPS_GPU_BUFFER_COUNT; ++i)
        fps_metal_bind_buffer(i, (void *)buffers[i]->id);
#else
    rlBindShaderBuffer((unsigned int)gpuPhysics.particle_buffer.id, 0);
    rlBindShaderBuffer((unsigned int)gpuPhysics.correction_buffer.id, 1);
    rlBindShaderBuffer((unsigned int)gpuPhysics.cell_buffer.id, 2);
    rlBindShaderBuffer((unsigned int)gpuPhysics.sim_id_buffer.id, 3);
    rlBindShaderBuffer((unsigned int)gpuPhysics.voxel_buffer.id, 4);
    rlBindShaderBuffer((unsigned int)gpuPhysics.hash_head_buffer.id, 5);
    rlBindShaderBuffer((unsigned int)gpuPhysics.hash_next_buffer.id, 6);
    rlBindShaderBuffer((unsigned int)gpuPhysics.static_cell_buffer.id, 7);
    rlBindShaderBuffer((unsigned int)gpuPhysics.tether_buffer.id, 8);
    rlBindShaderBuffer((unsigned int)gpuPhysics.collision_meta_buffer.id, 9);
    rlBindShaderBuffer((unsigned int)gpuPhysics.static_collider_buffer.id, 10);
    rlBindShaderBuffer((unsigned int)gpuPhysics.refcount_buffer.id, 11);
    rlBindShaderBuffer((unsigned int)gpuPhysics.control_buffer.id, 12);
    rlBindShaderBuffer((unsigned int)gpuPhysics.clone_parent_buffer.id, 13);
    rlBindShaderBuffer((unsigned int)gpuPhysics.dispatch_buffer.id, 14);
    rlBindShaderBuffer((unsigned int)gpuPhysics.topology_buffer.id, 15);
#endif
}

static bool gpu_upload_world(void) {
    double upload_started = GetTime();
    bool static_changed = gpuPhysics.static_uploaded_generation != staticHashGeneration;
    if (static_changed && !gpu_rebuild_static_cells()) return false;
    size_t particle_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuParticleState);
    size_t cell_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuCellState);
    size_t sim_bytes = (size_t)sim_particle_count * sizeof(uint32_t);
    size_t voxel_bytes = (size_t)gpuPhysics.dynamic_voxel_count * sizeof(GpuVoxelState);
    size_t static_bytes = (size_t)HASH_SIZE * sizeof(GpuStaticCell);
    size_t static_collider_bytes = (size_t)gpuPhysics.static_collider_count * sizeof(GpuStaticCollider);
    size_t collision_meta_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuCollisionMeta);
    size_t particle_capacity_bytes = (size_t)gpuPhysics.particle_slot_capacity * sizeof(GpuParticleState);
    size_t correction_capacity_bytes = (size_t)gpuPhysics.particle_slot_capacity * sizeof(GpuCorrectionState);
    size_t cell_capacity_bytes = (size_t)gpuPhysics.particle_slot_capacity * sizeof(GpuCellState);
    size_t sim_capacity_bytes = (size_t)gpuPhysics.sim_slot_capacity * sizeof(uint32_t);
    size_t collision_id_bytes = (size_t)gpuPhysics.collision_count * sizeof(uint32_t);
    size_t collision_meta_capacity_bytes = (size_t)gpuPhysics.particle_slot_capacity * sizeof(GpuCollisionMeta);
    if (!gpu_buffer_ensure(&gpuPhysics.particle_buffer, particle_capacity_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.correction_buffer, correction_capacity_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.cell_buffer, cell_capacity_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.sim_id_buffer, sim_capacity_bytes * 2u) ||
        !gpu_buffer_ensure(&gpuPhysics.voxel_buffer, voxel_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.hash_head_buffer, (size_t)PARTICLE_HASH_SIZE * sizeof(int32_t)) ||
        !gpu_buffer_ensure(&gpuPhysics.hash_next_buffer, sim_capacity_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.static_cell_buffer, static_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.static_collider_buffer, static_collider_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.tether_buffer, (size_t)gpuPhysics.particle_slot_capacity * sizeof(GpuTetherState)) ||
        !gpu_buffer_ensure(&gpuPhysics.collision_meta_buffer, collision_meta_capacity_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.refcount_buffer, (size_t)gpuPhysics.particle_slot_capacity * sizeof(int32_t)) ||
        !gpu_buffer_ensure(&gpuPhysics.control_buffer, sizeof(GpuControlState)) ||
        !gpu_buffer_ensure(&gpuPhysics.clone_parent_buffer, (size_t)gpuPhysics.particle_slot_capacity * sizeof(int32_t)) ||
        !gpu_buffer_ensure(&gpuPhysics.dispatch_buffer, 64) ||
        !gpu_buffer_ensure(&gpuPhysics.topology_buffer,
                           (size_t)gpuPhysics.dynamic_voxel_count * 2u * sizeof(GpuTopologyState))) return false;
    if (particle_bytes) gpu_buffer_update(&gpuPhysics.particle_buffer, gpuPhysics.particle_stage, particle_bytes, 0);
    if (cell_bytes) gpu_buffer_update(&gpuPhysics.cell_buffer, gpuPhysics.cell_stage, cell_bytes, 0);
    if (sim_bytes) gpu_buffer_update(&gpuPhysics.sim_id_buffer, gpuPhysics.sim_ids, sim_bytes, 0);
    if (collision_id_bytes) gpu_buffer_update(&gpuPhysics.sim_id_buffer,
                                              gpuPhysics.collision_ids, collision_id_bytes,
                                              sim_capacity_bytes);
    uint32_t dispatch_stage[16] = { 0 };
    dispatch_stage[8] = (uint32_t)gpuPhysics.collision_count;
    dispatch_stage[9] = (uint32_t)gpuPhysics.sim_slot_capacity;
    dispatch_stage[12] = getenv("FPS_COLLISION_STATIC_LEGACY") ? 0u : 1u;
    gpu_buffer_update(&gpuPhysics.dispatch_buffer, dispatch_stage, sizeof(dispatch_stage), 0);
    if (voxel_bytes) gpu_buffer_update(&gpuPhysics.voxel_buffer, gpuPhysics.voxel_stage, voxel_bytes, 0);
    if (gpuPhysics.compact_count) gpu_buffer_update(&gpuPhysics.tether_buffer, gpuPhysics.tether_stage,
                                                    (size_t)gpuPhysics.compact_count * sizeof(GpuTetherState), 0);
    if (collision_meta_bytes) gpu_buffer_update(&gpuPhysics.collision_meta_buffer,
                                                gpuPhysics.collision_meta_stage, collision_meta_bytes, 0);
    if (gpuPhysics.compact_count) {
        gpu_buffer_update(&gpuPhysics.refcount_buffer, gpuPhysics.refcount_stage,
                          (size_t)gpuPhysics.compact_count * sizeof(int32_t), 0);
        gpu_buffer_update(&gpuPhysics.clone_parent_buffer, gpuPhysics.clone_parent_stage,
                          (size_t)gpuPhysics.compact_count * sizeof(int32_t), 0);
    }
    gpu_buffer_update(&gpuPhysics.control_buffer, &gpuPhysics.control_stage,
                      sizeof(gpuPhysics.control_stage), 0);
    if (gpuPhysics.dynamic_voxel_count) {
        gpu_buffer_update(&gpuPhysics.topology_buffer, gpuPhysics.topology_stage,
                          (size_t)gpuPhysics.dynamic_voxel_count * 2u * sizeof(GpuTopologyState), 0);
    }
    if (static_changed) {
        gpu_buffer_update(&gpuPhysics.static_cell_buffer, gpuPhysics.static_stage, static_bytes, 0);
        if (static_collider_bytes) gpu_buffer_update(&gpuPhysics.static_collider_buffer,
                                                     gpuPhysics.static_collider_stage,
                                                     static_collider_bytes, 0);
        gpuPhysics.static_uploaded_generation = staticHashGeneration;
        gpuPhysics.static_uploads++;
    }
    gpu_bind_world_buffers();
    gpuPhysics.collision_uploaded_generation = collisionMetadataGeneration;
    gpuPhysics.upload_bytes += (uint64_t)(particle_bytes + cell_bytes + sim_bytes + voxel_bytes +
        (size_t)gpuPhysics.compact_count * sizeof(GpuTetherState) + collision_meta_bytes +
        (size_t)gpuPhysics.compact_count * sizeof(int32_t) * 2u + sizeof(GpuControlState) +
        (size_t)gpuPhysics.dynamic_voxel_count * 2u * sizeof(GpuTopologyState) +
        collision_id_bytes + sizeof(dispatch_stage));
    if (static_changed) gpuPhysics.upload_bytes += (uint64_t)(static_bytes + static_collider_bytes);
    gpuPhysics.full_uploads++;
    gpuPhysics.upload_ms += (GetTime() - upload_started) * 1000.0;
    return true;
}

#if defined(__APPLE__) && defined(FPS_GPU_METAL)
static bool gpu_dispatch(int mode, int count) { return fps_metal_dispatch(mode, count); }
static bool gpu_dispatch_indirect_offset(int mode, size_t offset) {
    return fps_metal_dispatch_indirect(mode, (void *)gpuPhysics.dispatch_buffer.id, offset);
}
static bool gpu_dispatch_indirect(int mode) { return gpu_dispatch_indirect_offset(mode, 0); }
#elif defined(GRAPHICS_API_OPENGL_43)
static bool gpu_dispatch(int mode, int count) {
    if (count <= 0) return true;
    while (glGetError() != GL_NO_ERROR) { }
    rlSetUniform(gpuPhysics.loc[GPU_LOC_MODE], &mode, SHADER_UNIFORM_INT, 1);
    unsigned int groups = ((unsigned int)count + 127u) / 128u;
    rlComputeShaderDispatch(groups, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
    return glGetError() == GL_NO_ERROR;
}

static bool gpu_dispatch_indirect_offset(int mode, size_t offset) {
    while (glGetError() != GL_NO_ERROR) { }
    rlSetUniform(gpuPhysics.loc[GPU_LOC_MODE], &mode, SHADER_UNIFORM_INT, 1);
    glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER, (unsigned int)gpuPhysics.dispatch_buffer.id);
    glDispatchComputeIndirect((GLintptr)offset);
    glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER, 0);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
    return glGetError() == GL_NO_ERROR;
}
static bool gpu_dispatch_indirect(int mode) { return gpu_dispatch_indirect_offset(mode, 0); }
#else
static bool gpu_dispatch(int mode, int count) { (void)mode; (void)count; return false; }
static bool gpu_dispatch_indirect(int mode) { (void)mode; return false; }
static bool gpu_dispatch_indirect_offset(int mode, size_t offset) { (void)mode; (void)offset; return false; }
#endif

static void gpu_set_uniforms(float dt) {
    int particle_count = gpuPhysics.compact_count;
    int dynamic_voxel_count = gpuPhysics.dynamic_voxel_count;
    int static_collider_count = gpuPhysics.static_collider_count;
    // This must exactly match the number of hash heads cleared below.  Using
    // slot capacity here while clearing by collision count leaves stale heads
    // in the upper buckets and can create self-referential chains.
    int hash_size = particle_hash_size_for_count(gpuPhysics.collision_count);
    int static_hash_size = HASH_SIZE;
    float voxel_size = VOXEL_SIZE, floor_size = FLOOR_SIZE, gravity = GRAVITY;
    float damping = VELOCITY_DAMPING, sor = PBD_SOR_FACTOR, collision = COLLISION_RELAXATION;
    float alpha = VGS_ALPHA, beta = VGS_BETA, eps = VGS_EPS;
    float strain = STRAIN_BREAK_THRESHOLD, shear = SHEAR_BREAK_THRESHOLD;
    float tether_spring = TETHER_SPRING, tether_damping = TETHER_DAMPING;
    int break_damp_frames = BREAK_DAMP_FRAMES;
    float rest_grid_step = VOXEL_SIZE * 0.5f;
    float player_data[MAX_PLAYERS][4];
    float tether_data[MAX_PLAYERS][4];
    for (int i = 0; i < MAX_PLAYERS; ++i) {
        player_data[i][0] = players[i].pos.x; player_data[i][1] = players[i].pos.y; player_data[i][2] = players[i].pos.z;
        player_data[i][3] = (i < activePlayers && players[i].respawn_timer <= 0.0f) ? PLAYER_SIZE * 0.5f : -1.0f;
        tether_data[i][0] = tetherTargetByPlayer[i].x;
        tether_data[i][1] = tetherTargetByPlayer[i].y;
        tether_data[i][2] = tetherTargetByPlayer[i].z;
        tether_data[i][3] = 1.0f;
    }
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    FpsGpuUniforms uniforms = { 0 };
    uniforms.particle_count = particle_count; uniforms.sim_count = sim_particle_count;
    uniforms.voxel_count = dynamic_voxel_count; uniforms.static_collider_count = static_collider_count;
    uniforms.hash_size = hash_size; uniforms.static_hash_size = static_hash_size;
    uniforms.active_players = activePlayers; uniforms.break_damp_frames = break_damp_frames;
    uniforms.dt = dt; uniforms.voxel_size = voxel_size; uniforms.floor_size = floor_size;
    uniforms.gravity = gravity; uniforms.velocity_damping = damping; uniforms.sor = sor;
    uniforms.collision_relaxation = collision; uniforms.vgs_alpha = alpha; uniforms.vgs_beta = beta;
    uniforms.vgs_epsilon = eps; uniforms.strain_threshold = strain; uniforms.shear_threshold = shear;
    uniforms.tether_spring = tether_spring; uniforms.tether_damping = tether_damping;
    uniforms.rest_grid_step = rest_grid_step;
    memcpy(uniforms.players, player_data, sizeof(player_data));
    memcpy(uniforms.tether_targets, tether_data, sizeof(tether_data));
    fps_metal_set_uniforms(&uniforms);
#else
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PARTICLE_COUNT], &particle_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_SIM_COUNT], &sim_particle_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VOXEL_COUNT], &dynamic_voxel_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_STATIC_COLLIDER_COUNT], &static_collider_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_HASH_SIZE], &hash_size, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_STATIC_HASH_SIZE], &static_hash_size, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_ACTIVE_PLAYERS], &activePlayers, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_DT], &dt, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VOXEL_SIZE], &voxel_size, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_FLOOR_SIZE], &floor_size, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_GRAVITY], &gravity, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VELOCITY_DAMPING], &damping, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_SOR], &sor, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_COLLISION_RELAXATION], &collision, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VGS_ALPHA], &alpha, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VGS_BETA], &beta, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VGS_EPSILON], &eps, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_STRAIN_THRESHOLD], &strain, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_SHEAR_THRESHOLD], &shear, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_TETHER_SPRING], &tether_spring, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_TETHER_DAMPING], &tether_damping, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_BREAK_DAMP_FRAMES], &break_damp_frames, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_REST_GRID_STEP], &rest_grid_step, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PLAYERS], player_data, SHADER_UNIFORM_VEC4, MAX_PLAYERS);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_TETHER_TARGETS], tether_data, SHADER_UNIFORM_VEC4, MAX_PLAYERS);
#endif
}

static bool gpu_validate_and_commit(void) {
    double read_started = GetTime();
    bool collision_diagnostics = getenv("FPS_COLLISION_DIAGNOSTICS") != NULL;
    gpu_buffer_read(&gpuPhysics.control_buffer, &gpuPhysics.control_stage,
                    sizeof(gpuPhysics.control_stage), 0);
    int used_particles = gpuPhysics.control_stage.particle_count;
    if (collision_diagnostics) {
        uint32_t dispatch_stage[12] = { 0 };
        gpu_buffer_read(&gpuPhysics.dispatch_buffer, dispatch_stage, sizeof(dispatch_stage), 0);
        fprintf(stderr,
                "collision-readback particles=%d/%d uploaded=%d sim=%d/%d flags=0x%x "
                "dispatch-full=%u collision=%u/%u overflow=%u chainGuards=%u\n",
                used_particles, gpuPhysics.particle_slot_capacity, gpuPhysics.uploaded_particle_count,
                gpuPhysics.control_stage.sim_count, gpuPhysics.sim_slot_capacity,
                gpuPhysics.control_stage.overflow, dispatch_stage[0], dispatch_stage[8],
                dispatch_stage[9], dispatch_stage[10], dispatch_stage[11]);
    }
    if (used_particles < gpuPhysics.uploaded_particle_count ||
        used_particles > gpuPhysics.particle_slot_capacity ||
        gpuPhysics.control_stage.sim_count < 0 ||
        gpuPhysics.control_stage.sim_count > gpuPhysics.sim_slot_capacity ||
        (gpuPhysics.control_stage.overflow & FPS_GPU_CONTROL_OVERFLOW) != 0) {
        if (collision_diagnostics) fprintf(stderr, "collision-validation: invalid GPU control state\n");
        return false;
    }
    size_t particle_bytes = (size_t)used_particles * sizeof(GpuParticleState);
    size_t cell_bytes = (size_t)used_particles * sizeof(GpuCellState);
    size_t voxel_bytes = (size_t)gpuPhysics.dynamic_voxel_count * sizeof(GpuVoxelState);
    if (particle_bytes) gpu_buffer_read(&gpuPhysics.particle_buffer, gpuPhysics.particle_stage, particle_bytes, 0);
    if (cell_bytes) gpu_buffer_read(&gpuPhysics.cell_buffer, gpuPhysics.cell_stage, cell_bytes, 0);
    if (voxel_bytes) gpu_buffer_read(&gpuPhysics.voxel_buffer, gpuPhysics.voxel_stage, voxel_bytes, 0);
    if (used_particles > gpuPhysics.uploaded_particle_count) {
        size_t offset = (size_t)gpuPhysics.uploaded_particle_count * sizeof(int32_t);
        size_t bytes = (size_t)(used_particles - gpuPhysics.uploaded_particle_count) * sizeof(int32_t);
        gpu_buffer_read(&gpuPhysics.clone_parent_buffer,
                        gpuPhysics.clone_parent_stage + gpuPhysics.uploaded_particle_count,
                        bytes, offset);
    }
    gpuPhysics.sync_count++;
    gpuPhysics.readback_bytes += (uint64_t)(sizeof(GpuControlState) + particle_bytes + cell_bytes + voxel_bytes +
        (size_t)(used_particles - gpuPhysics.uploaded_particle_count) * sizeof(int32_t));
    gpuPhysics.readback_ms += (GetTime() - read_started) * 1000.0;
    double commit_started = GetTime();
    for (int i = 0; i < used_particles; ++i) {
        GpuParticleState *g = &gpuPhysics.particle_stage[i];
        for (int c = 0; c < 3; ++c) {
            if (!isfinite(g->pos_radius[c]) || !isfinite(g->predicted_base_inv_mass[c]) ||
                !isfinite(g->velocity[c]) || fabsf(g->pos_radius[c]) > 100000.0f) {
                if (collision_diagnostics) fprintf(stderr, "collision-validation: invalid particle %d component %d\n", i, c);
                return false;
            }
        }
    }
    for (int i = 0; i < gpuPhysics.dynamic_voxel_count; ++i) {
        GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        if (!isfinite(g->pos_rest_edge[0]) || !isfinite(g->pos_rest_edge[1]) || !isfinite(g->pos_rest_edge[2])) {
            if (collision_diagnostics) fprintf(stderr, "collision-validation: invalid voxel %d\n", i);
            return false;
        }
    }
    for (int i = gpuPhysics.uploaded_particle_count; i < used_particles; ++i) {
        int parent_id = gpuPhysics.clone_parent_stage[i] & INT_MAX;
        if (parent_id < 0 || parent_id >= i || !gpuPhysics.compact_to_particle[parent_id]) {
            if (collision_diagnostics) fprintf(stderr, "collision-validation: clone %d has parent %d raw=0x%x\n",
                                               i, parent_id, (unsigned)gpuPhysics.clone_parent_stage[i]);
            return false;
        }
        Particle *clone = particle_clone(gpuPhysics.compact_to_particle[parent_id]);
        if (!clone) {
            if (collision_diagnostics) fprintf(stderr, "collision-validation: CPU clone allocation failed for %d\n", i);
            return false;
        }
        gpuPhysics.compact_to_particle[i] = clone;
        ptrdiff_t pool_id = clone - particles_pool;
        if (pool_id >= 0 && pool_id < MAX_PARTICLES) gpuPhysics.pool_to_compact[pool_id] = i;
    }
    for (int i = 0; i < used_particles; ++i) {
        Particle *p = gpuPhysics.compact_to_particle[i];
        if (!p) return false;
        GpuParticleState *g = &gpuPhysics.particle_stage[i];
        p->pos = (Vector3){g->pos_radius[0], g->pos_radius[1], g->pos_radius[2]};
        p->radius = g->pos_radius[3];
        p->prev_pos = (Vector3){g->prev_inv_mass[0], g->prev_inv_mass[1], g->prev_inv_mass[2]};
        p->inv_mass = g->prev_inv_mass[3];
        p->predicted_pos = (Vector3){g->predicted_base_inv_mass[0], g->predicted_base_inv_mass[1], g->predicted_base_inv_mass[2]};
        p->base_inv_mass = g->predicted_base_inv_mass[3];
        p->vel = (Vector3){g->velocity[0], g->velocity[1], g->velocity[2]};
        p->cell_x = gpuPhysics.cell_stage[i].value[0]; p->cell_y = gpuPhysics.cell_stage[i].value[1];
        p->cell_z = gpuPhysics.cell_stage[i].value[2]; p->break_timer = gpuPhysics.cell_stage[i].value[3];
    }
    bool topology_changed = used_particles != gpuPhysics.uploaded_particle_count;
    gpuPhysics.clones_created += (uint64_t)(used_particles - gpuPhysics.uploaded_particle_count);
    for (int i = 0; i < gpuPhysics.dynamic_voxel_count; ++i) {
        int voxel_id = gpuPhysics.compact_to_voxel[i];
        if (voxel_id < 0 || voxel_id >= voxel_count) return false;
        Voxel *v = &voxels[voxel_id]; GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        for (int corner = 0; corner < 8; ++corner) {
            uint32_t particle_id = g->particle[corner];
            if (particle_id >= (uint32_t)used_particles || !gpuPhysics.compact_to_particle[particle_id]) return false;
            Particle *desired = gpuPhysics.compact_to_particle[particle_id];
            Particle *old = v->particles[corner];
            if (old != desired) {
                particle_retain(desired);
                if (old) particle_release(old);
                v->particles[corner] = desired;
                topology_changed = true;
            }
        }
        if (v->type == 0 && !v->isBullet) {
            v->pos = (Vector3){g->pos_rest_edge[0], g->pos_rest_edge[1], g->pos_rest_edge[2]};
            v->vel = (Vector3){g->velocity_rest_volume[0], g->velocity_rest_volume[1], g->velocity_rest_volume[2]};
            v->skipCollisionVelocityFrames = (int)g->lifecycle[0];
        }
        v->wake_source = g->lifecycle[1] != 0;
        v->wake_timer = (int)g->bounds_min[3];
        v->break_mask = 0;
        for (int face = 0; face < 6; ++face) {
            bool glued = (g->lifecycle[2] & (1u << face)) != 0;
            if (v->glued_faces[face] != glued) topology_changed = true;
            v->glued_faces[face] = glued;
        }
    }
    for (int i = gpuPhysics.uploaded_particle_count; i < used_particles; ++i) {
        Particle *clone = gpuPhysics.compact_to_particle[i];
        if (clone) particle_release(clone);
    }
    if (topology_changed) collisionTopologyDirty = true;
    gpuPhysics.compact_count = used_particles;
    gpuPhysics.uploaded_particle_count = used_particles;
    gpuPhysics.resident_valid = true;
    gpuPhysics.commit_ms += (GetTime() - commit_started) * 1000.0;
    return true;
}

static bool gpu_resident_matches_cpu(void) {
    if (!gpuPhysics.resident_valid || gpuPhysics.static_uploaded_generation != staticHashGeneration ||
        gpuPhysics.collision_uploaded_generation != collisionMetadataGeneration) return false;
    int dynamic_count = 0;
    for (int i = 0; i < voxel_count; ++i) if (voxels[i].simulate) dynamic_count++;
    if (dynamic_count != gpuPhysics.dynamic_voxel_count) return false;
    for (int i = 0; i < gpuPhysics.compact_count; ++i) {
        Particle *p = gpuPhysics.compact_to_particle[i];
        if (!p || !p->active) continue;
        GpuParticleState *g = &gpuPhysics.particle_stage[i];
        if (p->pos.x != g->pos_radius[0] || p->pos.y != g->pos_radius[1] || p->pos.z != g->pos_radius[2] ||
            p->vel.x != g->velocity[0] || p->vel.y != g->velocity[1] || p->vel.z != g->velocity[2] ||
            p->inv_mass != g->prev_inv_mass[3] || p->base_inv_mass != g->predicted_base_inv_mass[3] ||
            p->break_timer != gpuPhysics.cell_stage[i].value[3]) return false;
    }
    for (int i = 0; i < gpuPhysics.dynamic_voxel_count; ++i) {
        int world_id = gpuPhysics.compact_to_voxel[i];
        if (world_id < 0 || world_id >= voxel_count || !voxels[world_id].simulate) return false;
        Voxel *v = &voxels[world_id];
        GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        if ((v->sleeping ? 0 : 1) != g->flags[0] ||
            v->type != g->flags[2] || (v->isBullet ? 1 : 0) != g->flags[1] ||
            (v->simulate_dofs ? 1 : 0) != g->flags[3]) return false;
        uint32_t glue_mask = gpu_glued_face_mask(v);
        if ((g->lifecycle[2] & 0x3fu) != glue_mask) return false;
        for (int corner = 0; corner < 8; ++corner) {
            Particle *p = v->particles[corner];
            if (!p) return false;
            ptrdiff_t pool_id = p - particles_pool;
            if (pool_id < 0 || pool_id >= particle_pool_count) return false;
            int compact_id = gpuPhysics.pool_to_compact[pool_id];
            if (compact_id < 0 || g->particle[corner] != (uint32_t)compact_id) return false;
        }
    }
    return true;
}

static bool gpu_refresh_tethers(void) {
    if (gpuPhysics.compact_count <= 0) return true;
    memcpy(gpuPhysics.tether_compare_stage, gpuPhysics.tether_stage,
           (size_t)gpuPhysics.compact_count * sizeof(*gpuPhysics.tether_compare_stage));
    for (int i = 0; i < gpuPhysics.compact_count; ++i) gpuPhysics.tether_stage[i] = -1;
    for (int i = 0; i < gpuPhysics.dynamic_voxel_count; ++i) {
        int world_id = gpuPhysics.compact_to_voxel[i];
        if (world_id < 0 || world_id >= voxel_count) return false;
        Voxel *v = &voxels[world_id];
        if (v->type != 0 || v->isBullet || tetherTag[world_id] <= 0) continue;
        int player = tetherTag[world_id] - 1;
        if (player < 0 || player >= MAX_PLAYERS) continue;
        for (int corner = 0; corner < 8; ++corner) {
            ptrdiff_t pool_id = v->particles[corner] - particles_pool;
            if (pool_id < 0 || pool_id >= particle_pool_count) continue;
            int id = gpuPhysics.pool_to_compact[pool_id];
            if (id < 0 || id >= gpuPhysics.compact_count || gpuPhysics.tether_stage[id] >= 0) continue;
            gpuPhysics.tether_stage[id] = player;
        }
    }
    bool changed = memcmp(gpuPhysics.tether_compare_stage, gpuPhysics.tether_stage,
                          (size_t)gpuPhysics.compact_count * sizeof(*gpuPhysics.tether_compare_stage)) != 0;
    if (changed) {
        gpu_buffer_update(&gpuPhysics.tether_buffer, gpuPhysics.tether_stage,
                          (size_t)gpuPhysics.compact_count * sizeof(GpuTetherState), 0);
        gpuPhysics.upload_bytes += (uint64_t)((size_t)gpuPhysics.compact_count * sizeof(GpuTetherState));
    }
    gpu_bind_world_buffers();
    gpuPhysics.resident_reuses++;
    return true;
}

static bool gpu_dispatch_substep(float sub_dt) {
    gpu_set_uniforms(sub_dt);
    int hash_size = particle_hash_size_for_count(gpuPhysics.collision_count);
    bool ok = gpu_dispatch(GPU_MODE_PREPARE_INDIRECT, 1) &&
              gpu_dispatch_indirect(GPU_MODE_INTEGRATE) &&
              gpu_dispatch_indirect(GPU_MODE_RESET) &&
              gpu_dispatch(GPU_MODE_HASH_CLEAR, hash_size) &&
              gpu_dispatch_indirect_offset(GPU_MODE_HASH_BUILD, 16) &&
              gpu_dispatch_indirect_offset(GPU_MODE_PAIR_COLLISIONS, 16) &&
              gpu_dispatch_indirect(GPU_MODE_APPLY);
    for (int it = 0; ok && it < PBD_CONSTRAINT_ITERS; ++it) {
        ok = gpu_dispatch_indirect(GPU_MODE_RESET) &&
             gpu_dispatch(GPU_MODE_VGS, gpuPhysics.dynamic_voxel_count) &&
             gpu_dispatch_indirect(GPU_MODE_APPLY);
    }
    ok = ok && gpu_dispatch_indirect_offset(GPU_MODE_STATIC_COLLISIONS, 16) &&
         gpu_dispatch(GPU_MODE_BREAK_MASK, gpuPhysics.dynamic_voxel_count) &&
         gpu_dispatch(GPU_MODE_SPLIT_BREAKS, gpuPhysics.dynamic_voxel_count) &&
         gpu_dispatch(GPU_MODE_WAKE_GATHER, gpuPhysics.dynamic_voxel_count) &&
         gpu_dispatch(GPU_MODE_WAKE_APPLY, gpuPhysics.dynamic_voxel_count) &&
         gpu_dispatch(GPU_MODE_PREPARE_INDIRECT, 1) &&
         gpu_dispatch_indirect(GPU_MODE_FINALIZE_PARTICLES) &&
         gpu_dispatch(GPU_MODE_FINALIZE_VOXELS, gpuPhysics.dynamic_voxel_count);
    return ok;
}

static bool gpu_rebuild_topology_if_needed(void) {
    return gpu_dispatch(GPU_MODE_TOPOLOGY_REBUILD_SERIAL, 1);
}

static int gpu_physics_steps(float dt, int fixed_steps) {
    if (fixed_steps <= 0 || getenv("FPS_FORCE_GPU_STEP_FAILURE") != NULL) return 0;
    if (gpuTransferMode == GPU_TRANSFER_RESIDENT && gpu_resident_matches_cpu()) {
        if (!gpu_refresh_tethers()) return 0;
    } else {
        gpuPhysics.resident_valid = false;
        double pack_started = GetTime();
        if (!gpu_pack_world()) return 0;
        gpuPhysics.pack_ms += (GetTime() - pack_started) * 1000.0;
        if (!gpu_upload_world()) return 0;
    }
    float sub_dt = dt / (float)PBD_SUBSTEPS;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    if (!fps_metal_begin_batch()) return 0;
#else
    rlEnableShader(gpuPhysics.program);
#endif
    double dispatch_started = GetTime();
    bool ok = true;
    for (int fixed = 0; ok && fixed < fixed_steps; ++fixed) {
        for (int step = 0; step < PBD_SUBSTEPS; ++step) {
            if (!gpu_dispatch_substep(sub_dt)) { ok = false; break; }
        }
        if (ok) ok = gpu_rebuild_topology_if_needed();
    }
    gpuPhysics.dispatch_ms += (GetTime() - dispatch_started) * 1000.0;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    bool batch_finished = fps_metal_end_batch();
#else
    rlDisableShader();
    bool batch_finished = true;
#endif
    if (!ok || !batch_finished || !gpu_validate_and_commit()) return 0;
    return fixed_steps;
}

static bool gpu_physics_init(void) {
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    physicsBackend.gpu_compiled = true;
    if (getenv("FPS_FORCE_GPU_INIT_FAILURE") != NULL) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "forced GPU initialization failure");
        return false;
    }
    if (!fps_metal_initialize("shaders/pbd/pbd_pipeline.metallib", &gpuPhysics.max_ssbo_size,
                              physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason))) {
        return false;
    }
    gpuPhysics.pool_to_compact = malloc((size_t)MAX_PARTICLES * sizeof(*gpuPhysics.pool_to_compact));
    if (!gpuPhysics.pool_to_compact) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "GPU particle map allocation failed");
        gpu_physics_release_resources();
        return false;
    }
    gpuPhysics.ready = true;
    gpuPhysics.static_uploaded_generation = 0;
    physicsBackend.gpu_available = true;
    return true;
#elif !defined(GRAPHICS_API_OPENGL_43)
    physicsBackend.gpu_compiled = false;
#if defined(__APPLE__)
    snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "binary was built without Metal compute support");
#else
    snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "binary was built without OpenGL 4.3 compute support");
#endif
    return false;
#else
    physicsBackend.gpu_compiled = true;
    if (getenv("FPS_FORCE_GPU_INIT_FAILURE") != NULL) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "forced GPU initialization failure");
        return false;
    }
    if (rlGetVersion() != RL_OPENGL_43) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "raylib did not create an OpenGL 4.3 context");
        return false;
    }
    GLint max_bindings = 0, max_invocations = 0;
    GLint64 max_block = 0;
    glGetIntegerv(GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS, &max_bindings);
    glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &max_invocations);
    glGetInteger64v(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &max_block);
    if (max_bindings < FPS_GPU_BUFFER_COUNT || max_invocations < 128 || max_block < 1024 * 1024) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason),
                 "insufficient compute limits (bindings=%d invocations=%d block=%lld)",
                 max_bindings, max_invocations, (long long)max_block);
        return false;
    }
    gpuPhysics.max_ssbo_size = (long long)max_block;
    char *source = LoadFileText("shaders/pbd/pbd_pipeline.comp");
    if (!source) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "missing shaders/pbd/pbd_pipeline.comp");
        return false;
    }
    unsigned int shader = glCreateShader(GL_COMPUTE_SHADER);
    const char *shader_source = source;
    glShaderSource(shader, 1, &shader_source, NULL);
    glCompileShader(shader);
    GLint shader_compiled = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_compiled);
    if (shader_compiled == GL_FALSE) {
        char compile_log[2048] = { 0 };
        glGetShaderInfoLog(shader, (GLsizei)sizeof(compile_log) - 1, NULL, compile_log);
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason),
                 "PBD compute shader compilation failed: %.1800s", compile_log);
        glDeleteShader(shader);
        shader = 0;
    }
    UnloadFileText(source);
    if (shader == 0) {
        if (!physicsBackend.fallback_reason[0])
            snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "PBD compute shader compilation failed");
        return false;
    }
    gpuPhysics.program = glCreateProgram();
    glAttachShader(gpuPhysics.program, shader);
    glLinkProgram(gpuPhysics.program);
    GLint program_linked = GL_FALSE;
    glGetProgramiv(gpuPhysics.program, GL_LINK_STATUS, &program_linked);
    if (program_linked == GL_FALSE) {
        char link_log[2048] = { 0 };
        glGetProgramInfoLog(gpuPhysics.program, (GLsizei)sizeof(link_log) - 1, NULL, link_log);
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason),
                 "PBD compute shader link failed: %.1800s", link_log);
        glDeleteProgram(gpuPhysics.program);
        gpuPhysics.program = 0;
    }
    glDeleteShader(shader);
    if (gpuPhysics.program == 0) {
        if (!physicsBackend.fallback_reason[0])
            snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "PBD compute shader link failed");
        return false;
    }
    static const char *names[GPU_LOC_COUNT] = {
        "uMode", "uParticleCount", "uSimCount", "uVoxelCount", "uStaticColliderCount", "uHashSize", "uStaticHashSize",
        "uActivePlayers", "uDt", "uVoxelSize", "uFloorSize", "uGravity", "uVelocityDamping",
        "uSor", "uCollisionRelaxation", "uVgsAlpha", "uVgsBeta", "uVgsEpsilon",
        "uStrainBreakThreshold", "uShearBreakThreshold", "uTetherSpring", "uTetherDamping",
        "uBreakDampFrames", "uRestGridStep", "uPlayerPosHalf", "uTetherTarget"
    };
    for (int i = 0; i < GPU_LOC_COUNT; ++i) {
        gpuPhysics.loc[i] = rlGetLocationUniform(gpuPhysics.program, names[i]);
        if (gpuPhysics.loc[i] < 0) {
            snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "compute shader uniform missing: %s", names[i]);
            gpu_physics_release_resources();
            return false;
        }
    }
    gpuPhysics.pool_to_compact = malloc((size_t)MAX_PARTICLES * sizeof(*gpuPhysics.pool_to_compact));
    if (!gpuPhysics.pool_to_compact) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "GPU particle map allocation failed");
        gpu_physics_release_resources();
        return false;
    }
    gpuPhysics.ready = true;
    gpuPhysics.static_uploaded_generation = 0;
    physicsBackend.gpu_available = true;
    return true;
#endif
}
