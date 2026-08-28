/* Included by fps_ray.c after the voxel/PBD helpers are defined. */

typedef struct { float pos_radius[4], prev_inv_mass[4], predicted_base_inv_mass[4], velocity[4]; } GpuParticleState;
typedef struct { uint32_t bits[4]; } GpuCorrectionState;
typedef struct { int32_t value[4]; } GpuCellState;
typedef struct { float value[4]; } GpuTetherState;
typedef struct { int32_t value[4]; } GpuStaticCell;
typedef struct { float density_lambda[4], delta[4]; } GpuFluidState;
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
_Static_assert(sizeof(GpuVoxelState) == 128, "GpuVoxelState/std430 mismatch");
_Static_assert(sizeof(GpuStaticCell) == 16, "GpuStaticCell/std430 mismatch");
_Static_assert(sizeof(GpuFluidState) == 32, "GpuFluidState/std430 mismatch");

typedef struct {
    unsigned int id;
    size_t capacity;
} GpuBuffer;

enum {
    GPU_MODE_RESET = 0,
    GPU_MODE_INTEGRATE,
    GPU_MODE_HASH_CLEAR,
    GPU_MODE_HASH_BUILD,
    GPU_MODE_SCENE_COLLISIONS,
    GPU_MODE_PAIR_COLLISIONS,
    GPU_MODE_APPLY,
    GPU_MODE_VGS,
    GPU_MODE_STATIC_COLLISIONS,
    GPU_MODE_BREAK_MASK,
    GPU_MODE_FINALIZE_PARTICLES,
    GPU_MODE_FINALIZE_VOXELS,
    GPU_MODE_PBF_LAMBDA,
    GPU_MODE_PBF_DELTA,
    GPU_MODE_PBF_APPLY,
    GPU_MODE_PBF_STATIC_COLLISIONS,
    GPU_MODE_PBF_VISCOSITY,
    GPU_MODE_PBF_APPLY_VISCOSITY
};

enum {
    GPU_LOC_MODE = 0,
    GPU_LOC_PARTICLE_COUNT,
    GPU_LOC_SIM_COUNT,
    GPU_LOC_FLUID_COUNT,
    GPU_LOC_VOXEL_COUNT,
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
    GPU_LOC_PBF_H,
    GPU_LOC_PBF_REST_DENSITY,
    GPU_LOC_PBF_PARTICLE_MASS,
    GPU_LOC_PBF_LAMBDA_EPSILON,
    GPU_LOC_PBF_SCORR_K,
    GPU_LOC_PBF_SCORR_N,
    GPU_LOC_PBF_SCORR_DELTA_Q,
    GPU_LOC_PBF_VISCOSITY,
    GPU_LOC_PLAYERS,
    GPU_LOC_COUNT
};

typedef struct {
    bool ready;
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
    GpuBuffer tether_buffer;
    GpuBuffer fluid_buffer;
    GpuParticleState *particle_stage;
    GpuCellState *cell_stage;
    GpuTetherState *tether_stage;
    Particle **compact_to_particle;
    uint32_t *sim_ids;
    GpuVoxelState *voxel_stage;
    GpuStaticCell *static_stage;
    GpuFluidState *fluid_stage;
    int *pool_to_compact;
    size_t particle_capacity;
    size_t sim_capacity;
    size_t voxel_capacity;
    uint64_t static_uploaded_generation;
    int compact_count;
    int committed_substeps;
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
    if (buffer->id != 0) rlUnloadShaderBuffer(buffer->id);
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
    buffer->id = rlLoadShaderBuffer((unsigned int)capacity, NULL, RL_DYNAMIC_COPY);
    if (buffer->id == 0) return false;
    buffer->capacity = capacity;
    return true;
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
    gpu_buffer_destroy(&gpuPhysics.tether_buffer);
    gpu_buffer_destroy(&gpuPhysics.fluid_buffer);
    if (gpuPhysics.program != 0) rlUnloadShaderProgram(gpuPhysics.program);
    gpuPhysics.program = 0;
    gpuPhysics.ready = false;
}

static void gpu_physics_shutdown(void) {
    gpu_physics_release_resources();
    free(gpuPhysics.particle_stage);
    free(gpuPhysics.cell_stage);
    free(gpuPhysics.tether_stage);
    free(gpuPhysics.compact_to_particle);
    free(gpuPhysics.sim_ids);
    free(gpuPhysics.voxel_stage);
    free(gpuPhysics.static_stage);
    free(gpuPhysics.fluid_stage);
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
            !gpu_resize_array((void **)&gpuPhysics.fluid_stage, cap, sizeof(*gpuPhysics.fluid_stage)) ||
            !gpu_resize_array((void **)&gpuPhysics.compact_to_particle, cap, sizeof(*gpuPhysics.compact_to_particle))) return false;
        gpuPhysics.particle_capacity = cap;
    }
    if (sims > gpuPhysics.sim_capacity) {
        size_t cap = gpuPhysics.sim_capacity ? gpuPhysics.sim_capacity : 256;
        while (cap < sims) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.sim_ids, cap, sizeof(*gpuPhysics.sim_ids))) return false;
        gpuPhysics.sim_capacity = cap;
    }
    if (voxels_count > gpuPhysics.voxel_capacity) {
        size_t cap = gpuPhysics.voxel_capacity ? gpuPhysics.voxel_capacity : 256;
        while (cap < voxels_count) cap *= 2;
        if (!gpu_resize_array((void **)&gpuPhysics.voxel_stage, cap, sizeof(*gpuPhysics.voxel_stage))) return false;
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
    g->velocity[0] = p->vel.x; g->velocity[1] = p->vel.y; g->velocity[2] = p->vel.z;
    g->velocity[3] = (float)p->material;
    memset(&gpuPhysics.fluid_stage[id], 0, sizeof(gpuPhysics.fluid_stage[id]));
    gpuPhysics.cell_stage[id].value[0] = p->cell_x;
    gpuPhysics.cell_stage[id].value[1] = p->cell_y;
    gpuPhysics.cell_stage[id].value[2] = p->cell_z;
    gpuPhysics.cell_stage[id].value[3] = p->break_timer;
    memset(&gpuPhysics.tether_stage[id], 0, sizeof(gpuPhysics.tether_stage[id]));
    return id;
}

static uint32_t gpu_glued_face_mask(const Voxel *v) {
    uint32_t mask = 0;
    for (int face = 0; face < 6; ++face) if (v->glued_faces[face]) mask |= 1u << face;
    return mask;
}

static bool gpu_pack_world(void) {
    size_t maximum_references = (size_t)sim_particle_count + (size_t)voxel_count * 8u;
    if (maximum_references > (size_t)particle_pool_count) maximum_references = (size_t)particle_pool_count;
    if (!gpu_stage_reserve(maximum_references,
                           (size_t)sim_particle_count + (size_t)fluid_particle_count,
                           (size_t)voxel_count)) return false;
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
    for (int i = 0; i < fluid_particle_count; ++i) {
        Particle *particle = fluid_particles[i];
        if (!particle) return false;
        int id = gpu_add_particle(particle);
        if (id < 0) return false;
        gpuPhysics.sim_ids[sim_particle_count + i] = (uint32_t)id;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        for (int c = 0; c < 8; ++c) if (gpu_add_particle(v->particles[c]) < 0) return false;
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        memset(g, 0, sizeof(*g));
        for (int c = 0; c < 8; ++c) {
            g->particle[c] = UINT32_MAX;
            if (v->simulate && v->particles[c]) {
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
        g->flags[0] = v->simulate ? 1 : 0; g->flags[1] = v->isBullet ? 1 : 0;
        g->flags[2] = v->type; g->flags[3] = v->simulate_dofs ? 1 : 0;
        VoxelWorldBounds bounds;
        voxel_world_bounds(v, &bounds);
        g->bounds_min[0] = bounds.minx; g->bounds_min[1] = bounds.miny; g->bounds_min[2] = bounds.minz;
        g->bounds_max[0] = bounds.maxx; g->bounds_max[1] = bounds.maxy; g->bounds_max[2] = bounds.maxz;
        g->lifecycle[0] = (uint32_t)((v->skipCollisionVelocityFrames > 0) ? v->skipCollisionVelocityFrames : 0);
        g->lifecycle[1] = v->wake_source ? 1u : 0u;
        g->lifecycle[2] = gpu_glued_face_mask(v);
        g->lifecycle[3] = v->break_mask;
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
            if (id < 0 || gpuPhysics.tether_stage[id].value[3] > 0.5f) continue;
            gpuPhysics.tether_stage[id].value[0] = tetherTargetByPlayer[player].x;
            gpuPhysics.tether_stage[id].value[1] = tetherTargetByPlayer[player].y;
            gpuPhysics.tether_stage[id].value[2] = tetherTargetByPlayer[player].z;
            gpuPhysics.tether_stage[id].value[3] = 1.0f;
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

static bool gpu_rebuild_static_cells(void) {
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
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *v = &voxels[i];
        if (v->simulate) continue;
        for (int z = v->min_gz; z <= v->max_gz; ++z) {
            for (int y = v->min_gy; y <= v->max_gy; ++y) {
                for (int x = v->min_gx; x <= v->max_gx; ++x) {
                    uint32_t h = gpu_static_hash_coord(x, y, z);
                    for (int probe = 0; probe < HASH_SIZE; ++probe) {
                        GpuStaticCell *entry = &gpuPhysics.static_stage[h];
                        if (entry->value[3] < 0 ||
                            (entry->value[0] == x && entry->value[1] == y && entry->value[2] == z)) {
                            entry->value[0] = x; entry->value[1] = y; entry->value[2] = z; entry->value[3] = i;
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

static bool gpu_upload_world(void) {
    size_t particle_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuParticleState);
    size_t correction_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuCorrectionState);
    size_t cell_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuCellState);
    size_t sim_bytes = ((size_t)sim_particle_count + (size_t)fluid_particle_count) * sizeof(uint32_t);
    size_t voxel_bytes = (size_t)voxel_count * sizeof(GpuVoxelState);
    size_t fluid_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuFluidState);
    size_t next_bytes = ((size_t)sim_particle_count + (size_t)fluid_particle_count) * sizeof(int32_t);
    size_t static_bytes = (size_t)HASH_SIZE * sizeof(GpuStaticCell);
    if (!gpu_buffer_ensure(&gpuPhysics.particle_buffer, particle_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.correction_buffer, correction_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.cell_buffer, cell_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.sim_id_buffer, sim_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.voxel_buffer, voxel_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.hash_head_buffer, (size_t)PARTICLE_HASH_SIZE * sizeof(int32_t)) ||
        !gpu_buffer_ensure(&gpuPhysics.hash_next_buffer, next_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.static_cell_buffer, static_bytes) ||
        !gpu_buffer_ensure(&gpuPhysics.tether_buffer, (size_t)gpuPhysics.compact_count * sizeof(GpuTetherState)) ||
        !gpu_buffer_ensure(&gpuPhysics.fluid_buffer, fluid_bytes)) return false;
    if (particle_bytes) rlUpdateShaderBuffer(gpuPhysics.particle_buffer.id, gpuPhysics.particle_stage, (unsigned int)particle_bytes, 0);
    if (cell_bytes) rlUpdateShaderBuffer(gpuPhysics.cell_buffer.id, gpuPhysics.cell_stage, (unsigned int)cell_bytes, 0);
    if (sim_bytes) rlUpdateShaderBuffer(gpuPhysics.sim_id_buffer.id, gpuPhysics.sim_ids, (unsigned int)sim_bytes, 0);
    if (voxel_bytes) rlUpdateShaderBuffer(gpuPhysics.voxel_buffer.id, gpuPhysics.voxel_stage, (unsigned int)voxel_bytes, 0);
    if (gpuPhysics.compact_count) rlUpdateShaderBuffer(gpuPhysics.tether_buffer.id, gpuPhysics.tether_stage,
                                                       (unsigned int)((size_t)gpuPhysics.compact_count * sizeof(GpuTetherState)), 0);
    if (fluid_bytes) rlUpdateShaderBuffer(gpuPhysics.fluid_buffer.id, gpuPhysics.fluid_stage,
                                          (unsigned int)fluid_bytes, 0);
    if (gpuPhysics.static_uploaded_generation != staticHashGeneration) {
        if (!gpu_rebuild_static_cells()) return false;
        rlUpdateShaderBuffer(gpuPhysics.static_cell_buffer.id, gpuPhysics.static_stage, (unsigned int)static_bytes, 0);
        gpuPhysics.static_uploaded_generation = staticHashGeneration;
    }
    rlBindShaderBuffer(gpuPhysics.particle_buffer.id, 0);
    rlBindShaderBuffer(gpuPhysics.correction_buffer.id, 1);
    rlBindShaderBuffer(gpuPhysics.cell_buffer.id, 2);
    rlBindShaderBuffer(gpuPhysics.sim_id_buffer.id, 3);
    rlBindShaderBuffer(gpuPhysics.voxel_buffer.id, 4);
    rlBindShaderBuffer(gpuPhysics.hash_head_buffer.id, 5);
    rlBindShaderBuffer(gpuPhysics.hash_next_buffer.id, 6);
    rlBindShaderBuffer(gpuPhysics.static_cell_buffer.id, 7);
    rlBindShaderBuffer(gpuPhysics.tether_buffer.id, 8);
    rlBindShaderBuffer(gpuPhysics.fluid_buffer.id, 9);
    return true;
}

#if defined(GRAPHICS_API_OPENGL_43)
static bool gpu_dispatch(int mode, int count) {
    if (count <= 0) return true;
    while (glGetError() != GL_NO_ERROR) { }
    rlSetUniform(gpuPhysics.loc[GPU_LOC_MODE], &mode, SHADER_UNIFORM_INT, 1);
    unsigned int groups = ((unsigned int)count + 127u) / 128u;
    rlComputeShaderDispatch(groups, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);
    return glGetError() == GL_NO_ERROR;
}
#else
static bool gpu_dispatch(int mode, int count) { (void)mode; (void)count; return false; }
#endif

static void gpu_set_uniforms(float dt) {
    int particle_count = gpuPhysics.compact_count;
    int fluid_count = fluid_particle_count;
    int hash_size = PARTICLE_HASH_SIZE;
    int static_hash_size = HASH_SIZE;
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PARTICLE_COUNT], &particle_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_SIM_COUNT], &sim_particle_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_FLUID_COUNT], &fluid_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_VOXEL_COUNT], &voxel_count, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_HASH_SIZE], &hash_size, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_STATIC_HASH_SIZE], &static_hash_size, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_ACTIVE_PLAYERS], &activePlayers, SHADER_UNIFORM_INT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_DT], &dt, SHADER_UNIFORM_FLOAT, 1);
    float voxel_size = VOXEL_SIZE, floor_size = FLOOR_SIZE, gravity = GRAVITY;
    float damping = VELOCITY_DAMPING, sor = PBD_SOR_FACTOR, collision = COLLISION_RELAXATION;
    float alpha = VGS_ALPHA, beta = VGS_BETA, eps = VGS_EPS;
    float strain = STRAIN_BREAK_THRESHOLD, shear = SHEAR_BREAK_THRESHOLD;
    float tether_spring = TETHER_SPRING, tether_damping = TETHER_DAMPING;
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
    float pbf_h = PBF_SMOOTHING_RADIUS;
    float pbf_rest_density = PBF_REST_DENSITY;
    float pbf_particle_mass = PBF_PARTICLE_MASS;
    float pbf_lambda_epsilon = PBF_LAMBDA_EPSILON;
    float pbf_scorr_k = PBF_SCORR_K;
    float pbf_scorr_n = PBF_SCORR_N;
    float pbf_scorr_delta_q = PBF_SCORR_DELTA_Q;
    float pbf_viscosity = PBF_XSPH_VISCOSITY;
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_H], &pbf_h, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_REST_DENSITY], &pbf_rest_density, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_PARTICLE_MASS], &pbf_particle_mass, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_LAMBDA_EPSILON], &pbf_lambda_epsilon, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_SCORR_K], &pbf_scorr_k, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_SCORR_N], &pbf_scorr_n, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_SCORR_DELTA_Q], &pbf_scorr_delta_q, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PBF_VISCOSITY], &pbf_viscosity, SHADER_UNIFORM_FLOAT, 1);
    float player_data[MAX_PLAYERS][4];
    for (int i = 0; i < MAX_PLAYERS; ++i) {
        player_data[i][0] = players[i].pos.x; player_data[i][1] = players[i].pos.y; player_data[i][2] = players[i].pos.z;
        player_data[i][3] = (i < activePlayers && players[i].respawn_timer <= 0.0f) ? PLAYER_SIZE * 0.5f : -1.0f;
    }
    rlSetUniform(gpuPhysics.loc[GPU_LOC_PLAYERS], player_data, SHADER_UNIFORM_VEC4, MAX_PLAYERS);
}

static bool gpu_validate_and_commit(void) {
    size_t particle_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuParticleState);
    size_t cell_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuCellState);
    size_t voxel_bytes = (size_t)voxel_count * sizeof(GpuVoxelState);
    size_t fluid_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuFluidState);
    if (particle_bytes) rlReadShaderBuffer(gpuPhysics.particle_buffer.id, gpuPhysics.particle_stage, (unsigned int)particle_bytes, 0);
    if (cell_bytes) rlReadShaderBuffer(gpuPhysics.cell_buffer.id, gpuPhysics.cell_stage, (unsigned int)cell_bytes, 0);
    if (voxel_bytes) rlReadShaderBuffer(gpuPhysics.voxel_buffer.id, gpuPhysics.voxel_stage, (unsigned int)voxel_bytes, 0);
    if (fluid_particle_count > 0 && fluid_bytes) {
        rlReadShaderBuffer(gpuPhysics.fluid_buffer.id, gpuPhysics.fluid_stage, (unsigned int)fluid_bytes, 0);
    }
    for (int i = 0; i < gpuPhysics.compact_count; ++i) {
        GpuParticleState *g = &gpuPhysics.particle_stage[i];
        for (int c = 0; c < 3; ++c) {
            if (!isfinite(g->pos_radius[c]) || !isfinite(g->predicted_base_inv_mass[c]) ||
                !isfinite(g->velocity[c]) || fabsf(g->pos_radius[c]) > 100000.0f) return false;
        }
    }
    for (int i = 0; i < voxel_count; ++i) {
        GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        if (!isfinite(g->pos_rest_edge[0]) || !isfinite(g->pos_rest_edge[1]) || !isfinite(g->pos_rest_edge[2])) return false;
    }
    for (int i = 0; i < gpuPhysics.compact_count; ++i) {
        Particle *p = gpuPhysics.compact_to_particle[i];
        GpuParticleState *g = &gpuPhysics.particle_stage[i];
        p->pos = (Vector3){g->pos_radius[0], g->pos_radius[1], g->pos_radius[2]};
        p->radius = g->pos_radius[3];
        p->prev_pos = (Vector3){g->prev_inv_mass[0], g->prev_inv_mass[1], g->prev_inv_mass[2]};
        p->inv_mass = g->prev_inv_mass[3];
        p->predicted_pos = (Vector3){g->predicted_base_inv_mass[0], g->predicted_base_inv_mass[1], g->predicted_base_inv_mass[2]};
        p->base_inv_mass = g->predicted_base_inv_mass[3];
        p->vel = (Vector3){g->velocity[0], g->velocity[1], g->velocity[2]};
        if (p->material == PARTICLE_MATERIAL_FLUID) {
            GpuFluidState *fluid = &gpuPhysics.fluid_stage[i];
            if (!isfinite(fluid->density_lambda[0]) || !isfinite(fluid->density_lambda[1])) return false;
            p->pbf_density = fluid->density_lambda[0];
            p->pbf_lambda = fluid->density_lambda[1];
        }
        p->cell_x = gpuPhysics.cell_stage[i].value[0]; p->cell_y = gpuPhysics.cell_stage[i].value[1];
        p->cell_z = gpuPhysics.cell_stage[i].value[2]; p->break_timer = gpuPhysics.cell_stage[i].value[3];
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i]; GpuVoxelState *g = &gpuPhysics.voxel_stage[i];
        if (v->simulate && !v->isBullet &&
            (v->type == VOXEL_TYPE_SOLID || voxel_is_fluid(v))) {
            v->pos = (Vector3){g->pos_rest_edge[0], g->pos_rest_edge[1], g->pos_rest_edge[2]};
            v->vel = (Vector3){g->velocity_rest_volume[0], g->velocity_rest_volume[1], g->velocity_rest_volume[2]};
            v->skipCollisionVelocityFrames = (int)g->lifecycle[0];
        }
        v->wake_source = g->lifecycle[1] != 0;
        v->break_mask = (uint8_t)(g->lifecycle[3] & 0x3fu);
    }
    return true;
}

static bool gpu_run_substep(float sub_dt) {
    if (!gpu_pack_world() || !gpu_upload_world()) return false;
    rlEnableShader(gpuPhysics.program);
    gpu_set_uniforms(sub_dt);
    int dynamic_count = sim_particle_count + fluid_particle_count;
    bool ok = gpu_dispatch(GPU_MODE_INTEGRATE, dynamic_count) &&
              gpu_dispatch(GPU_MODE_RESET, dynamic_count) &&
              gpu_dispatch(GPU_MODE_HASH_CLEAR, PARTICLE_HASH_SIZE) &&
              gpu_dispatch(GPU_MODE_HASH_BUILD, dynamic_count) &&
              gpu_dispatch(GPU_MODE_SCENE_COLLISIONS, dynamic_count) &&
              gpu_dispatch(GPU_MODE_PAIR_COLLISIONS, dynamic_count) &&
              gpu_dispatch(GPU_MODE_APPLY, dynamic_count);
    for (int it = 0; ok && it < PBF_DENSITY_ITERS; ++it) {
        ok = gpu_dispatch(GPU_MODE_PBF_LAMBDA, fluid_particle_count) &&
             gpu_dispatch(GPU_MODE_PBF_DELTA, fluid_particle_count) &&
             gpu_dispatch(GPU_MODE_PBF_APPLY, fluid_particle_count) &&
             gpu_dispatch(GPU_MODE_PBF_STATIC_COLLISIONS, fluid_particle_count);
    }
    for (int it = 0; ok && it < PBD_CONSTRAINT_ITERS; ++it) {
        ok = gpu_dispatch(GPU_MODE_RESET, dynamic_count) &&
             gpu_dispatch(GPU_MODE_VGS, voxel_count) &&
             gpu_dispatch(GPU_MODE_APPLY, dynamic_count);
    }
    ok = ok && gpu_dispatch(GPU_MODE_STATIC_COLLISIONS, dynamic_count) &&
         gpu_dispatch(GPU_MODE_BREAK_MASK, voxel_count) &&
         gpu_dispatch(GPU_MODE_FINALIZE_PARTICLES, dynamic_count) &&
         gpu_dispatch(GPU_MODE_PBF_VISCOSITY, fluid_particle_count) &&
         gpu_dispatch(GPU_MODE_PBF_APPLY_VISCOSITY, fluid_particle_count) &&
         gpu_dispatch(GPU_MODE_FINALIZE_VOXELS, voxel_count);
    rlDisableShader();
    if (!ok || !gpu_validate_and_commit()) return false;
    process_break_masks();
    update_wake_timers();
    return true;
}

static int gpu_physics_step(float dt) {
    gpuPhysics.committed_substeps = 0;
    if (getenv("FPS_FORCE_GPU_STEP_FAILURE") != NULL) return 0;
    float sub_dt = dt / (float)PBD_SUBSTEPS;
    for (int step = 0; step < PBD_SUBSTEPS; ++step) {
        if (!gpu_run_substep(sub_dt)) return gpuPhysics.committed_substeps;
        gpuPhysics.committed_substeps++;
    }
    return gpuPhysics.committed_substeps;
}

static bool gpu_physics_init(void) {
#if !defined(GRAPHICS_API_OPENGL_43)
    physicsBackend.gpu_compiled = false;
    snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "binary was built without OpenGL 4.3 compute support");
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
    if (max_bindings < 10 || max_invocations < 128 || max_block < 1024 * 1024) {
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
    unsigned int shader = rlCompileShader(source, RL_COMPUTE_SHADER);
    UnloadFileText(source);
    if (shader == 0) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "PBD compute shader compilation failed");
        return false;
    }
    gpuPhysics.program = rlLoadComputeShaderProgram(shader);
    glDeleteShader(shader);
    if (gpuPhysics.program == 0) {
        snprintf(physicsBackend.fallback_reason, sizeof(physicsBackend.fallback_reason), "PBD compute shader link failed");
        return false;
    }
    static const char *names[GPU_LOC_COUNT] = {
        "uMode", "uParticleCount", "uSimCount", "uFluidCount", "uVoxelCount", "uHashSize", "uStaticHashSize",
        "uActivePlayers", "uDt", "uVoxelSize", "uFloorSize", "uGravity", "uVelocityDamping",
        "uSor", "uCollisionRelaxation", "uVgsAlpha", "uVgsBeta", "uVgsEpsilon",
        "uStrainBreakThreshold", "uShearBreakThreshold", "uTetherSpring", "uTetherDamping",
        "uPbfH", "uPbfRestDensity", "uPbfParticleMass", "uPbfLambdaEpsilon",
        "uPbfScorrK", "uPbfScorrN", "uPbfScorrDeltaQ", "uPbfViscosity",
        "uPlayerPosHalf"
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
