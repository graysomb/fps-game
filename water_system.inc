/* Sparse 8^3 water chunks. Voxels remain a read-only obstacle source. */
typedef struct WaterChunk {
    uint16_t mass[WATER_CHUNK_CELLS];
    uint16_t scratch[WATER_CHUNK_CELLS];
    uint16_t batch_snapshot[WATER_CHUNK_CELLS];
    uint32_t velocity[WATER_CHUNK_CELLS];
    uint32_t velocity_scratch[WATER_CHUNK_CELLS];
    uint32_t velocity_batch_snapshot[WATER_CHUNK_CELLS];
    uint8_t static_solid[WATER_CHUNK_CELLS];
    uint8_t dynamic_solid[WATER_CHUNK_CELLS];
    uint8_t previous_solid[WATER_CHUNK_CELLS];
} WaterChunk;

typedef struct WaterSystem {
    WaterChunk *chunks[WATER_TILE_COUNT];
    uint16_t resident_tile_ids[WATER_TILE_COUNT];
    uint16_t resident_tile_count;
    uint8_t tile_active[WATER_TILE_COUNT];
    uint8_t batch_tile_active[WATER_TILE_COUNT];
    uint8_t tile_still_steps[WATER_TILE_COUNT];
    uint8_t tile_changed[WATER_TILE_COUNT];
    uint8_t tile_stats_dirty[WATER_TILE_COUNT];
    uint8_t tile_render_dirty[WATER_TILE_COUNT];
    uint16_t active_tile_ids[WATER_TILE_COUNT];
    uint16_t active_tile_count;
    uint16_t tile_wet_cells[WATER_TILE_COUNT];
    uint16_t tile_moving_cells[WATER_TILE_COUNT];
    uint8_t tile_max_component_speed[WATER_TILE_COUNT];
    uint64_t tile_visible_mass[WATER_TILE_COUNT];
    uint64_t tile_trapped_mass[WATER_TILE_COUNT];
    uint64_t static_generation;
    uint64_t mass_generation;
    uint64_t dynamic_generation;
    uint64_t dynamic_mask_hash;
    uint64_t resident_generation;
    bool initialized;
    bool render_dirty;
    bool cpu_scratch_valid;
    WaterDiagnostics diagnostics;
} WaterSystem;

static WaterSystem waterSystem;
static inline void water_tile_coords(int tile, int *tx, int *ty, int *tz);
static inline void water_coords_from_index(int index, int *gx, int *gy, int *gz);
static void water_wake_cell_and_neighbors(int gx, int gy, int gz);

typedef struct WaterVelocity {
    int32_t x, y, z;
} WaterVelocity;

static inline int32_t water_clamp_velocity(int32_t value) {
    if (value < -WATER_VELOCITY_MAX) return -WATER_VELOCITY_MAX;
    if (value > WATER_VELOCITY_MAX) return WATER_VELOCITY_MAX;
    return value;
}

static inline uint32_t water_pack_velocity(WaterVelocity velocity) {
    return ((uint32_t)(uint8_t)water_clamp_velocity(velocity.x)) |
           ((uint32_t)(uint8_t)water_clamp_velocity(velocity.y) << 8) |
           ((uint32_t)(uint8_t)water_clamp_velocity(velocity.z) << 16);
}

static inline int32_t water_unpack_component(uint32_t packed, int shift) {
    int32_t value = (int32_t)((packed >> shift) & 255u);
    return value > 127 ? value - 256 : value;
}

static inline WaterVelocity water_unpack_velocity(uint32_t packed) {
    return (WaterVelocity){ water_unpack_component(packed, 0),
                            water_unpack_component(packed, 8),
                            water_unpack_component(packed, 16) };
}

static inline WaterVelocity water_accelerated_velocity(uint32_t packed, uint32_t mass) {
    if (!mass) return (WaterVelocity){ 0, 0, 0 };
    WaterVelocity velocity = water_unpack_velocity(packed);
    velocity.y = water_clamp_velocity(velocity.y - WATER_VELOCITY_GRAVITY);
    return velocity;
}

static inline WaterVelocity water_mix_velocity(uint32_t mass_a, WaterVelocity a,
                                                uint32_t mass_b, WaterVelocity b) {
    uint32_t total = mass_a + mass_b;
    if (!total) return (WaterVelocity){ 0, 0, 0 };
    return (WaterVelocity){
        ((int32_t)mass_a * a.x + (int32_t)mass_b * b.x) / (int32_t)total,
        ((int32_t)mass_a * a.y + (int32_t)mass_b * b.y) / (int32_t)total,
        ((int32_t)mass_a * a.z + (int32_t)mass_b * b.z) / (int32_t)total
    };
}

static WaterChunk *water_chunk_get(int tile, bool create) {
    if (tile < 0 || tile >= WATER_TILE_COUNT) return NULL;
    WaterChunk *chunk = waterSystem.chunks[tile];
    if (chunk || !create) return chunk;
    chunk = (WaterChunk *)calloc(1, sizeof(*chunk));
    if (!chunk) return NULL;
    int tx, ty, tz;
    water_tile_coords(tile, &tx, &ty, &tz);
    for (int local = 0; local < WATER_CHUNK_CELLS; ++local) {
        int lx = local & 7, lz = (local >> 3) & 7, ly = local >> 6;
        int gx = WATER_MIN_X + tx * 8 + lx;
        int gy = WATER_MIN_Y + ty * 8 + ly;
        int gz = WATER_MIN_Z + tz * 8 + lz;
        int voxel = table_get_static_only(gx, gy, gz);
        chunk->static_solid[local] = (uint8_t)(voxel >= 0 && voxel < voxel_count &&
                                               !voxels[voxel].simulate);
        chunk->previous_solid[local] = chunk->static_solid[local];
    }
    waterSystem.chunks[tile] = chunk;
    waterSystem.resident_tile_ids[waterSystem.resident_tile_count++] = (uint16_t)tile;
    if (++waterSystem.resident_generation == 0) waterSystem.resident_generation = 1;
    waterSystem.tile_stats_dirty[tile] = 1;
    waterSystem.tile_render_dirty[tile] = 1;
    return chunk;
}

static inline int water_chunk_local_for_index(int index, int *tile_out) {
    int y = index / (WATER_SIZE_X * WATER_SIZE_Z);
    int rem = index - y * WATER_SIZE_X * WATER_SIZE_Z;
    int z = rem / WATER_SIZE_X;
    int x = rem - z * WATER_SIZE_X;
    int tx = x / WATER_TILE_SIZE, ty = y / WATER_TILE_SIZE, tz = z / WATER_TILE_SIZE;
    *tile_out = tx + WATER_TILE_COUNT_X * (tz + WATER_TILE_COUNT_Z * ty);
    return (x & 7) + WATER_TILE_SIZE * ((z & 7) + WATER_TILE_SIZE * (y & 7));
}

static inline WaterChunk *water_chunk_for_index(int index, bool create, int *local_out) {
    int tile;
    int local = water_chunk_local_for_index(index, &tile);
    if (local_out) *local_out = local;
    return water_chunk_get(tile, create);
}

static inline uint16_t water_mass_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk ? chunk->mass[local] : 0;
}

static inline uint16_t water_scratch_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk ? chunk->scratch[local] : 0;
}

static inline uint32_t water_velocity_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk ? chunk->velocity[local] : 0;
}

static inline uint32_t water_velocity_scratch_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk ? chunk->velocity_scratch[local] : 0;
}

static inline bool water_static_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk && chunk->static_solid[local];
}

static inline bool water_dynamic_index(int index) {
    int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    return chunk && chunk->dynamic_solid[local];
}

static void water_free_chunks(void) {
    for (int i = 0; i < waterSystem.resident_tile_count; ++i) {
        int tile = waterSystem.resident_tile_ids[i];
        free(waterSystem.chunks[tile]);
        waterSystem.chunks[tile] = NULL;
    }
    waterSystem.resident_tile_count = 0;
}

static void water_mark_tile_dirty(int tile) {
    if (tile < 0 || tile >= WATER_TILE_COUNT) return;
    waterSystem.tile_stats_dirty[tile] = 1;
    waterSystem.tile_render_dirty[tile] = 1;
    waterSystem.render_dirty = true;
}

static void water_note_cpu_mass_edit(void) {
    if (++waterSystem.mass_generation == 0) waterSystem.mass_generation = 1;
}

static WaterVelocity water_velocity_from_world(Vector3 world_velocity) {
    float scale = 127.0f * PBD_MAX_STEP_DT / VOXEL_SIZE;
    return (WaterVelocity){
        water_clamp_velocity((int32_t)lroundf(world_velocity.x * scale)),
        water_clamp_velocity((int32_t)lroundf(world_velocity.y * scale)),
        water_clamp_velocity((int32_t)lroundf(world_velocity.z * scale))
    };
}

static void water_add_impulse_index(int index, Vector3 world_velocity, int32_t upward_boost) {
    int local;
    WaterChunk *chunk = water_chunk_for_index(index, false, &local);
    if (!chunk || !chunk->mass[local]) return;
    WaterVelocity old = water_unpack_velocity(chunk->velocity[local]);
    WaterVelocity impulse = water_velocity_from_world(world_velocity);
    old.x = water_clamp_velocity(old.x + impulse.x);
    old.y = water_clamp_velocity(old.y + impulse.y + upward_boost);
    old.z = water_clamp_velocity(old.z + impulse.z);
    uint32_t packed = water_pack_velocity(old);
    if (packed == chunk->velocity[local]) return;
    chunk->velocity[local] = packed;
    chunk->velocity_scratch[local] = packed;
    water_note_cpu_mass_edit();
    int gx, gy, gz;
    water_coords_from_index(index, &gx, &gy, &gz);
    water_wake_cell_and_neighbors(gx, gy, gz);
}

static void water_mark_tile_and_neighbors_dirty(int tile) {
    int tx, ty, tz;
    water_tile_coords(tile, &tx, &ty, &tz);
    static const int dx[7] = { 0, 1, -1, 0, 0, 0, 0 };
    static const int dy[7] = { 0, 0, 0, 1, -1, 0, 0 };
    static const int dz[7] = { 0, 0, 0, 0, 0, 1, -1 };
    for (int d = 0; d < 7; ++d) {
        int nx = tx + dx[d], ny = ty + dy[d], nz = tz + dz[d];
        if (nx < 0 || nx >= WATER_TILE_COUNT_X || ny < 0 || ny >= WATER_TILE_COUNT_Y ||
            nz < 0 || nz >= WATER_TILE_COUNT_Z) continue;
        water_mark_tile_dirty(nx + WATER_TILE_COUNT_X * (nz + WATER_TILE_COUNT_Z * ny));
    }
}

static inline bool water_runtime_enabled(void) {
    return netTransport.role == NET_ROLE_OFFLINE && netRequestedRole == NET_ROLE_OFFLINE;
}

static inline bool water_in_domain(int gx, int gy, int gz) {
    return gx >= WATER_MIN_X && gx <= WATER_MAX_X &&
           gy >= WATER_MIN_Y && gy <= WATER_MAX_Y &&
           gz >= WATER_MIN_Z && gz <= WATER_MAX_Z;
}

static inline int water_index_unchecked(int gx, int gy, int gz) {
    int x = gx - WATER_MIN_X;
    int y = gy - WATER_MIN_Y;
    int z = gz - WATER_MIN_Z;
    return x + WATER_SIZE_X * (z + WATER_SIZE_Z * y);
}

static inline void water_coords_from_index(int index, int *gx, int *gy, int *gz) {
    int y = index / (WATER_SIZE_X * WATER_SIZE_Z);
    int rem = index - y * WATER_SIZE_X * WATER_SIZE_Z;
    int z = rem / WATER_SIZE_X;
    int x = rem - z * WATER_SIZE_X;
    *gx = x + WATER_MIN_X;
    *gy = y + WATER_MIN_Y;
    *gz = z + WATER_MIN_Z;
}

static inline int water_tile_index_for_cell(int gx, int gy, int gz) {
    int tx = (gx - WATER_MIN_X) / WATER_TILE_SIZE;
    int ty = (gy - WATER_MIN_Y) / WATER_TILE_SIZE;
    int tz = (gz - WATER_MIN_Z) / WATER_TILE_SIZE;
    return tx + WATER_TILE_COUNT_X * (tz + WATER_TILE_COUNT_Z * ty);
}

static inline void water_tile_coords(int tile, int *tx, int *ty, int *tz) {
    *ty = tile / (WATER_TILE_COUNT_X * WATER_TILE_COUNT_Z);
    int rem = tile - *ty * WATER_TILE_COUNT_X * WATER_TILE_COUNT_Z;
    *tz = rem / WATER_TILE_COUNT_X;
    *tx = rem - *tz * WATER_TILE_COUNT_X;
}

static void water_wake_tile_xyz(int tx, int ty, int tz) {
    if (tx < 0 || tx >= WATER_TILE_COUNT_X ||
        ty < 0 || ty >= WATER_TILE_COUNT_Y ||
        tz < 0 || tz >= WATER_TILE_COUNT_Z) return;
    int tile = tx + WATER_TILE_COUNT_X * (tz + WATER_TILE_COUNT_Z * ty);
    if (!water_chunk_get(tile, true)) return;
    waterSystem.tile_active[tile] = 1;
    waterSystem.tile_still_steps[tile] = 0;
    water_mark_tile_dirty(tile);
}

static void water_wake_cell_and_neighbors(int gx, int gy, int gz) {
    if (!water_in_domain(gx, gy, gz)) return;
    int tx = (gx - WATER_MIN_X) / WATER_TILE_SIZE;
    int ty = (gy - WATER_MIN_Y) / WATER_TILE_SIZE;
    int tz = (gz - WATER_MIN_Z) / WATER_TILE_SIZE;
    water_wake_tile_xyz(tx, ty, tz);
    water_wake_tile_xyz(tx + 1, ty, tz);
    water_wake_tile_xyz(tx - 1, ty, tz);
    water_wake_tile_xyz(tx, ty + 1, tz);
    water_wake_tile_xyz(tx, ty - 1, tz);
    water_wake_tile_xyz(tx, ty, tz + 1);
    water_wake_tile_xyz(tx, ty, tz - 1);
}

static inline bool water_cell_blocked_index(int index) {
    return water_static_index(index) || water_dynamic_index(index);
}

static void water_init(void) {
    if (waterSystem.initialized) return;
    memset(&waterSystem, 0, sizeof(waterSystem));
    waterSystem.initialized = true;
    waterSystem.static_generation = UINT64_MAX;
    waterSystem.render_dirty = true;
    waterSystem.cpu_scratch_valid = true;
    memset(waterSystem.tile_stats_dirty, 1, sizeof(waterSystem.tile_stats_dirty));
    memset(waterSystem.tile_render_dirty, 1, sizeof(waterSystem.tile_render_dirty));
    waterSystem.diagnostics.backend = "cpu-st";
}

static void water_reset(void) {
    water_gpu_shutdown_backend();
    water_free_chunks();
    memset(&waterSystem, 0, sizeof(waterSystem));
    waterSystem.initialized = true;
    waterSystem.static_generation = UINT64_MAX;
    waterSystem.render_dirty = true;
    memset(waterSystem.tile_stats_dirty, 1, sizeof(waterSystem.tile_stats_dirty));
    memset(waterSystem.tile_render_dirty, 1, sizeof(waterSystem.tile_render_dirty));
    waterSystem.diagnostics.backend = "cpu-st";
}

static void water_shutdown(void) {
    water_free_chunks();
    memset(&waterSystem, 0, sizeof(waterSystem));
}

static uint32_t water_get_mass(int gx, int gy, int gz) {
    if (!waterSystem.initialized || !water_runtime_enabled() ||
        !water_in_domain(gx, gy, gz)) return 0;
    return water_mass_index(water_index_unchecked(gx, gy, gz));
}

static bool water_set_mass(int gx, int gy, int gz, uint32_t mass) {
    if (!waterSystem.initialized) water_init();
    if (!water_runtime_enabled() || !water_in_domain(gx, gy, gz)) return false;
    int index = water_index_unchecked(gx, gy, gz);
    if (water_cell_blocked_index(index) || table_get(gx, gy, gz) >= 0) return false;
    if (mass > WATER_MAX_MASS) mass = WATER_MAX_MASS;
    int local; WaterChunk *chunk = water_chunk_for_index(index, mass != 0, &local);
    if (!chunk) return mass == 0;
    if (chunk->mass[local] == mass && chunk->scratch[local] == mass &&
        chunk->velocity[local] == 0 && chunk->velocity_scratch[local] == 0) return true;
    chunk->mass[local] = (uint16_t)mass;
    chunk->scratch[local] = (uint16_t)mass;
    chunk->velocity[local] = 0;
    chunk->velocity_scratch[local] = 0;
    water_note_cpu_mass_edit();
    water_wake_cell_and_neighbors(gx, gy, gz);
    water_mark_tile_and_neighbors_dirty(water_tile_index_for_cell(gx, gy, gz));
    waterSystem.render_dirty = true;
    return true;
}

static void water_add_region(int minx, int maxx, int miny, int maxy,
                             int minz, int maxz, uint32_t mass) {
    if (mass > WATER_MAX_MASS) mass = WATER_MAX_MASS;
    if (minx < WATER_MIN_X) minx = WATER_MIN_X;
    if (maxx > WATER_MAX_X) maxx = WATER_MAX_X;
    if (miny < WATER_MIN_Y) miny = WATER_MIN_Y;
    if (maxy > WATER_MAX_Y) maxy = WATER_MAX_Y;
    if (minz < WATER_MIN_Z) minz = WATER_MIN_Z;
    if (maxz > WATER_MAX_Z) maxz = WATER_MAX_Z;
    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            for (int x = minx; x <= maxx; ++x) {
                uint32_t old = water_get_mass(x, y, z);
                uint32_t next = old + mass;
                if (next < old || next > WATER_MAX_MASS) next = WATER_MAX_MASS;
                (void)water_set_mass(x, y, z, next);
            }
        }
    }
}

static void water_remove_region(int minx, int maxx, int miny, int maxy,
                                int minz, int maxz) {
    if (!waterSystem.initialized || !water_runtime_enabled()) return;
    if (minx < WATER_MIN_X) minx = WATER_MIN_X;
    if (maxx > WATER_MAX_X) maxx = WATER_MAX_X;
    if (miny < WATER_MIN_Y) miny = WATER_MIN_Y;
    if (maxy > WATER_MAX_Y) maxy = WATER_MAX_Y;
    if (minz < WATER_MIN_Z) minz = WATER_MIN_Z;
    if (maxz > WATER_MAX_Z) maxz = WATER_MAX_Z;
    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            for (int x = minx; x <= maxx; ++x) {
                int index = water_index_unchecked(x, y, z);
                int local; WaterChunk *chunk = water_chunk_for_index(index, false, &local);
                if (!chunk) continue;
                chunk->mass[local] = 0;
                chunk->scratch[local] = 0;
                chunk->velocity[local] = 0;
                chunk->velocity_scratch[local] = 0;
                water_note_cpu_mass_edit();
                water_wake_cell_and_neighbors(x, y, z);
                water_mark_tile_and_neighbors_dirty(water_tile_index_for_cell(x, y, z));
            }
        }
    }
    waterSystem.render_dirty = true;
}

static float water_sample_fill(Vector3 world_position) {
    int gx = (int)floorf(world_position.x / VOXEL_SIZE);
    int gy = (int)floorf(world_position.y / VOXEL_SIZE);
    int gz = (int)floorf(world_position.z / VOXEL_SIZE);
    return (float)water_get_mass(gx, gy, gz) / (float)WATER_MAX_MASS;
}

static WorldCell world_cell_get(int gx, int gy, int gz) {
    int voxel = table_get(gx, gy, gz);
    if (voxel >= 0) return (WorldCell){ WORLD_CELL_SOLID, voxel, 0 };
    uint32_t mass = water_get_mass(gx, gy, gz);
    if (mass > 0) return (WorldCell){ WORLD_CELL_WATER, -1, mass };
    return (WorldCell){ WORLD_CELL_EMPTY, -1, 0 };
}

static void water_rebuild_static_mask(void) {
    if (waterSystem.static_generation == staticHashGeneration) return;
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r];
        WaterChunk *chunk = waterSystem.chunks[tile];
        int tx, ty, tz; water_tile_coords(tile, &tx, &ty, &tz);
        for (int ly = 0; ly < WATER_TILE_SIZE; ++ly)
            for (int lz = 0; lz < WATER_TILE_SIZE; ++lz)
                for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                int gx = WATER_MIN_X + tx * 8 + lx;
                int gy = WATER_MIN_Y + ty * 8 + ly;
                int gz = WATER_MIN_Z + tz * 8 + lz;
                int local = lx + 8 * (lz + 8 * ly);
                int voxel = table_get_static_only(gx, gy, gz);
                uint8_t blocked = (uint8_t)(voxel >= 0 && voxel < voxel_count && !voxels[voxel].simulate);
                chunk->static_solid[local] = blocked;
            }
    }
    waterSystem.static_generation = staticHashGeneration;
}

static void water_rebuild_dynamic_mask(void) {
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        WaterChunk *chunk = waterSystem.chunks[waterSystem.resident_tile_ids[r]];
        memset(chunk->dynamic_solid, 0, sizeof(chunk->dynamic_solid));
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
        float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
        for (int c = 0; c < 8; ++c) {
            Particle *p = v->particles[c];
            if (!p) continue;
            if (p->pos.x < minx) minx = p->pos.x;
            if (p->pos.y < miny) miny = p->pos.y;
            if (p->pos.z < minz) minz = p->pos.z;
            if (p->pos.x > maxx) maxx = p->pos.x;
            if (p->pos.y > maxy) maxy = p->pos.y;
            if (p->pos.z > maxz) maxz = p->pos.z;
        }
        if (minx == FLT_MAX) {
            minx = maxx = v->pos.x;
            miny = maxy = v->pos.y;
            minz = maxz = v->pos.z;
        }
        int gx0 = (int)floorf(minx / VOXEL_SIZE);
        int gy0 = (int)floorf(miny / VOXEL_SIZE);
        int gz0 = (int)floorf(minz / VOXEL_SIZE);
        int gx1 = (int)floorf((maxx - 1e-5f) / VOXEL_SIZE);
        int gy1 = (int)floorf((maxy - 1e-5f) / VOXEL_SIZE);
        int gz1 = (int)floorf((maxz - 1e-5f) / VOXEL_SIZE);
        if (gx0 < WATER_MIN_X) gx0 = WATER_MIN_X;
        if (gy0 < WATER_MIN_Y) gy0 = WATER_MIN_Y;
        if (gz0 < WATER_MIN_Z) gz0 = WATER_MIN_Z;
        if (gx1 > WATER_MAX_X) gx1 = WATER_MAX_X;
        if (gy1 > WATER_MAX_Y) gy1 = WATER_MAX_Y;
        if (gz1 > WATER_MAX_Z) gz1 = WATER_MAX_Z;
        for (int gy = gy0; gy <= gy1; ++gy)
            for (int gz = gz0; gz <= gz1; ++gz)
                for (int gx = gx0; gx <= gx1; ++gx)
                {
                    int index = water_index_unchecked(gx, gy, gz), local;
                    WaterChunk *chunk = water_chunk_for_index(index, false, &local);
                    if (chunk) {
                        if (!chunk->previous_solid[local] && chunk->mass[local] > 0)
                            water_add_impulse_index(index, v->vel, 0);
                        chunk->dynamic_solid[local] = 1;
                    }
                }
    }
    uint64_t hash = UINT64_C(1469598103934665603);
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r];
        WaterChunk *chunk = waterSystem.chunks[tile];
        for (int i = 0; i < 512; ++i) if (chunk->dynamic_solid[i]) {
            hash ^= ((uint64_t)(uint32_t)tile << 9) | (uint32_t)i;
            hash *= UINT64_C(1099511628211);
        }
    }
    if (hash != waterSystem.dynamic_mask_hash) {
        waterSystem.dynamic_mask_hash = hash;
        if (++waterSystem.dynamic_generation == 0) waterSystem.dynamic_generation = 1;
    }
}

static uint32_t water_displace_from_cell(int source_index) {
    uint32_t remaining = water_mass_index(source_index);
    if (remaining == 0) return 0;
    WaterVelocity source_velocity = water_unpack_velocity(water_velocity_index(source_index));
    enum { BFS_SIDE = 17, BFS_CAPACITY = 17 * 17 * 17 };
    int queue[BFS_CAPACITY];
    uint8_t depth[BFS_CAPACITY], seen[BFS_CAPACITY];
    memset(seen, 0, sizeof(seen));
    uint32_t head = 0, tail = 0;
    queue[tail] = source_index;
    depth[tail++] = 0;
    int sgx, sgy, sgz; water_coords_from_index(source_index, &sgx, &sgy, &sgz);
    seen[8 + BFS_SIDE * (8 + BFS_SIDE * 8)] = 1;
    uint32_t moved = 0;
    static const int dx[6] = { 0, 1, -1, 0, 0, 0 };
    static const int dy[6] = { 1, 0, 0, 0, 0, -1 };
    static const int dz[6] = { 0, 0, 0, 1, -1, 0 };
    while (head < tail && remaining > 0) {
        int index = queue[head];
        int cell_depth = depth[head++];
        int gx, gy, gz;
        water_coords_from_index(index, &gx, &gy, &gz);
        if (index != source_index && !water_cell_blocked_index(index)) {
            uint32_t current = water_mass_index(index);
            uint32_t capacity = WATER_MAX_MASS - current;
            uint32_t amount = remaining < capacity ? remaining : capacity;
            if (amount > 0) {
                int local; WaterChunk *chunk = water_chunk_for_index(index, true, &local);
                if (!chunk) break;
                WaterVelocity target_velocity = water_unpack_velocity(chunk->velocity[local]);
                WaterVelocity mixed = water_mix_velocity(current, target_velocity,
                                                          amount, source_velocity);
                chunk->mass[local] = (uint16_t)(current + amount);
                chunk->scratch[local] = chunk->mass[local];
                chunk->velocity[local] = water_pack_velocity(mixed);
                chunk->velocity_scratch[local] = chunk->velocity[local];
                water_note_cpu_mass_edit();
                remaining -= amount;
                moved += amount;
                water_wake_cell_and_neighbors(gx, gy, gz);
            }
        }
        if (cell_depth >= 8) continue;
        for (int d = 0; d < 6; ++d) {
            int nx = gx + dx[d], ny = gy + dy[d], nz = gz + dz[d];
            if (!water_in_domain(nx, ny, nz)) continue;
            int rx = nx - sgx + 8, ry = ny - sgy + 8, rz = nz - sgz + 8;
            int seen_index = rx + BFS_SIDE * (rz + BFS_SIDE * ry);
            if (seen[seen_index]) continue;
            seen[seen_index] = 1;
            if (tail >= BFS_CAPACITY) continue;
            queue[tail] = water_index_unchecked(nx, ny, nz);
            depth[tail++] = (uint8_t)(cell_depth + 1);
        }
    }
    int source_local; WaterChunk *source_chunk = water_chunk_for_index(source_index, true, &source_local);
    source_chunk->mass[source_local] = (uint16_t)remaining;
    source_chunk->scratch[source_local] = (uint16_t)remaining;
    source_chunk->velocity[source_local] = 0;
    source_chunk->velocity_scratch[source_local] = 0;
    water_note_cpu_mass_edit();
    return moved;
}

static void water_rebuild_obstacles_and_displace(void) {
    if (!waterSystem.initialized || !water_runtime_enabled()) return;
    water_rebuild_static_mask();
    water_rebuild_dynamic_mask();
    int resident_count = waterSystem.resident_tile_count;
    for (int r = 0; r < resident_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r];
        WaterChunk *chunk = waterSystem.chunks[tile];
        int tx, ty, tz; water_tile_coords(tile, &tx, &ty, &tz);
        for (int local = 0; local < 512; ++local) {
        int lx = local & 7, lz = (local >> 3) & 7, ly = local >> 6;
        int gx = WATER_MIN_X + tx * 8 + lx;
        int gy = WATER_MIN_Y + ty * 8 + ly;
        int gz = WATER_MIN_Z + tz * 8 + lz;
        int i = water_index_unchecked(gx, gy, gz);
        uint8_t blocked = chunk->static_solid[local] || chunk->dynamic_solid[local];
        if (blocked && !chunk->previous_solid[local] && chunk->mass[local] > 0) {
            waterSystem.diagnostics.displaced_mass += water_displace_from_cell(i);
        } else if (!blocked && chunk->previous_solid[local]) {
            int gx, gy, gz;
            water_coords_from_index(i, &gx, &gy, &gz);
            water_wake_cell_and_neighbors(gx, gy, gz);
        }
        if (blocked != chunk->previous_solid[local]) {
            water_mark_tile_and_neighbors_dirty(water_tile_index_for_cell(gx, gy, gz));
        }
        chunk->previous_solid[local] = blocked;
        }
    }
}

static void water_build_active_tile_list(void) {
    waterSystem.active_tile_count = 0;
    for (int tile = 0; tile < WATER_TILE_COUNT; ++tile) {
        uint8_t active = waterSystem.tile_active[tile];
        waterSystem.batch_tile_active[tile] = active;
        waterSystem.tile_changed[tile] = 0;
        if (active) waterSystem.active_tile_ids[waterSystem.active_tile_count++] = (uint16_t)tile;
    }
}

static inline bool water_neighbor_tile_active(int tx, int ty, int tz) {
    if (tx < 0 || tx >= WATER_TILE_COUNT_X || ty < 0 || ty >= WATER_TILE_COUNT_Y ||
        tz < 0 || tz >= WATER_TILE_COUNT_Z) return false;
    return waterSystem.batch_tile_active[tx + WATER_TILE_COUNT_X *
           (tz + WATER_TILE_COUNT_Z * ty)] != 0;
}

enum {
    WATER_DIRECTION_XP = 0, WATER_DIRECTION_XM = 1,
    WATER_DIRECTION_ZP = 2, WATER_DIRECTION_ZM = 3,
    WATER_DIRECTION_UP = 4
};

static inline bool water_index_active(int index) {
    int tile;
    (void)water_chunk_local_for_index(index, &tile);
    return waterSystem.batch_tile_active[tile] != 0;
}

static inline bool water_laterally_supported_index(int index) {
    int y = index / (WATER_SIZE_X * WATER_SIZE_Z);
    if (y == 0) return true;
    int below = index - WATER_SIZE_X * WATER_SIZE_Z;
    return water_cell_blocked_index(below) ||
           water_scratch_index(below) >= WATER_LATERAL_SUPPORT_MASS;
}

static inline WaterVelocity water_effective_velocity_index(int index) {
    uint32_t mass = water_scratch_index(index);
    if (!mass || water_cell_blocked_index(index)) return (WaterVelocity){ 0, 0, 0 };
    WaterVelocity velocity = water_unpack_velocity(water_velocity_scratch_index(index));
    int y = index / (WATER_SIZE_X * WATER_SIZE_Z);
    bool impact = y == 0;
    if (y > 0) {
        int below = index - WATER_SIZE_X * WATER_SIZE_Z;
        if (water_cell_blocked_index(below)) {
            impact = true;
        } else if (water_scratch_index(below) >= WATER_LATERAL_SUPPORT_MASS) {
            WaterVelocity below_velocity =
                water_unpack_velocity(water_velocity_scratch_index(below));
            impact = below_velocity.y >= -WATER_VELOCITY_IMPACT_THRESHOLD;
        }
    }
    if (impact && velocity.y < 0) {
        velocity.y = velocity.y < -WATER_VELOCITY_IMPACT_THRESHOLD
            ? ((-velocity.y) * 3) / 4 : 0;
    }
    return velocity;
}

static inline int32_t water_velocity_for_direction(WaterVelocity velocity, int direction) {
    if (direction == WATER_DIRECTION_XP) return velocity.x;
    if (direction == WATER_DIRECTION_XM) return -velocity.x;
    if (direction == WATER_DIRECTION_ZP) return velocity.z;
    if (direction == WATER_DIRECTION_ZM) return -velocity.z;
    return velocity.y;
}

static inline int32_t water_surface_splash_speed(int index) {
    int32_t fall_speed = 0;
    bool impact = false;
    for (int depth = 0; depth < 8; ++depth) {
        uint32_t mass = water_scratch_index(index);
        if (!mass || water_cell_blocked_index(index)) return 0;
        WaterVelocity velocity = water_unpack_velocity(water_velocity_scratch_index(index));
        if (-velocity.y > fall_speed) fall_speed = -velocity.y;
        int y = index / (WATER_SIZE_X * WATER_SIZE_Z);
        if (y == 0) { impact = true; break; }
        int below = index - WATER_SIZE_X * WATER_SIZE_Z;
        if (water_cell_blocked_index(below)) { impact = true; break; }
        if (water_scratch_index(below) < WATER_LATERAL_SUPPORT_MASS) return 0;
        WaterVelocity below_velocity =
            water_unpack_velocity(water_velocity_scratch_index(below));
        if (below_velocity.y >= -WATER_VELOCITY_IMPACT_THRESHOLD) {
            impact = true;
            break;
        }
        index = below;
    }
    if (!impact || fall_speed <= WATER_VELOCITY_IMPACT_THRESHOLD) return 0;
    uint32_t hash = (uint32_t)index * UINT32_C(747796405) + UINT32_C(2891336453);
    if ((hash >> 30) != 0) return 0;
    int32_t splash = fall_speed * 2;
    return splash > WATER_VELOCITY_MAX ? WATER_VELOCITY_MAX : splash;
}

static inline uint32_t water_directional_flow(int from, int to, int direction) {
    if (from < 0 || from >= WATER_CELL_COUNT || to < 0 || to >= WATER_CELL_COUNT ||
        !water_index_active(from) || !water_index_active(to) ||
        water_cell_blocked_index(from) || water_cell_blocked_index(to)) return 0;
    uint32_t a = water_scratch_index(from), b = water_scratch_index(to);
    if (!a || b >= WATER_MAX_MASS) return 0;
    if (direction != WATER_DIRECTION_UP) {
        int to_y = to / (WATER_SIZE_X * WATER_SIZE_Z);
        if (to_y > 0) {
            int below_to = to - WATER_SIZE_X * WATER_SIZE_Z;
            uint32_t below_mass = water_scratch_index(below_to);
            if (water_surface_splash_speed(below_to) > 0 ||
                (below_mass < WATER_LATERAL_SUPPORT_MASS &&
                 water_effective_velocity_index(below_to).y > 0)) return 0;
        }
    }
    uint32_t base = 0;
    if (direction != WATER_DIRECTION_UP && water_laterally_supported_index(from) && a > b)
        base = (a - b) >> 3;
    WaterVelocity effective_velocity = water_effective_velocity_index(from);
    int32_t speed = water_velocity_for_direction(effective_velocity, direction);
    int32_t surface_splash = water_surface_splash_speed(from);
    bool rising_droplet = direction == WATER_DIRECTION_UP && surface_splash == 0 &&
                          a < WATER_LATERAL_SUPPORT_MASS && effective_velocity.y > 0;
    if (direction != WATER_DIRECTION_UP &&
        (surface_splash > 0 ||
         (a < WATER_LATERAL_SUPPORT_MASS && effective_velocity.y > 0))) return 0;
    if (direction == WATER_DIRECTION_UP && surface_splash == 0 &&
        effective_velocity.y > 0 && !rising_droplet) return 0;
    if (direction == WATER_DIRECTION_UP) {
        if (surface_splash > speed) speed = surface_splash;
    } else {
        int32_t outward_splash = effective_velocity.y > surface_splash
            ? effective_velocity.y : surface_splash;
        if (outward_splash > 0) speed += outward_splash / 2;
    }
    if (speed > WATER_VELOCITY_MAX) speed = WATER_VELOCITY_MAX;
    bool splash_launch = direction == WATER_DIRECTION_UP && surface_splash > 0;
    uint32_t momentum = rising_droplet ? a : (speed > 0
        ? (uint32_t)(((uint64_t)a * (uint32_t)speed) / (127u * 4u)) : 0);
    uint32_t flow = base + momentum;
    uint32_t source_cap = rising_droplet ? a : (splash_launch ? a / 4u : a / 5u);
    uint32_t destination_cap = (rising_droplet || splash_launch)
        ? WATER_MAX_MASS - b : (WATER_MAX_MASS - b) / 5u;
    if (flow > source_cap) flow = source_cap;
    if (flow > destination_cap) flow = destination_cap;
    return flow;
}

static inline WaterVelocity water_carried_velocity(int from, int direction) {
    WaterVelocity velocity = water_effective_velocity_index(from);
    if (direction == WATER_DIRECTION_UP) {
        int32_t surface_splash = water_surface_splash_speed(from);
        if (surface_splash > 0) {
            if (surface_splash > velocity.y) velocity.y = surface_splash;
            uint32_t hash = (uint32_t)from * UINT32_C(747796405) + UINT32_C(2891336453);
            velocity.x = water_clamp_velocity(velocity.x + ((hash & 1u) ? 64 : -64));
            velocity.z = water_clamp_velocity(velocity.z + ((hash & 2u) ? 64 : -64));
        }
    }
    if (direction == WATER_DIRECTION_XP)
        velocity.x = water_clamp_velocity(velocity.x + WATER_VELOCITY_HORIZONTAL_KICK);
    else if (direction == WATER_DIRECTION_XM)
        velocity.x = water_clamp_velocity(velocity.x - WATER_VELOCITY_HORIZONTAL_KICK);
    else if (direction == WATER_DIRECTION_ZP)
        velocity.z = water_clamp_velocity(velocity.z + WATER_VELOCITY_HORIZONTAL_KICK);
    else if (direction == WATER_DIRECTION_ZM)
        velocity.z = water_clamp_velocity(velocity.z - WATER_VELOCITY_HORIZONTAL_KICK);
    return velocity;
}

static inline WaterChunk *water_neighbor_chunk(int tx, int ty, int tz) {
    if (tx < 0 || tx >= WATER_TILE_COUNT_X || ty < 0 || ty >= WATER_TILE_COUNT_Y ||
        tz < 0 || tz >= WATER_TILE_COUNT_Z) return NULL;
    return waterSystem.chunks[tx + WATER_TILE_COUNT_X * (tz + WATER_TILE_COUNT_Z * ty)];
}

static void water_gravity_tiles(int start, int end, int worker_id, void *user) {
    (void)worker_id; (void)user;
    for (int job = start; job < end; ++job) {
        int tile = waterSystem.active_tile_ids[job];
        int tx, ty, tz;
        water_tile_coords(tile, &tx, &ty, &tz);
        WaterChunk *chunk = waterSystem.chunks[tile];
        WaterChunk *below_chunk = water_neighbor_chunk(tx, ty - 1, tz);
        WaterChunk *above_chunk = water_neighbor_chunk(tx, ty + 1, tz);
        int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
        for (int ly = 0; ly < WATER_TILE_SIZE; ++ly) {
            int y = y0 + ly;
            bool below_active = ly > 0 || water_neighbor_tile_active(tx, ty - 1, tz);
            bool above_active = ly + 1 < WATER_TILE_SIZE || water_neighbor_tile_active(tx, ty + 1, tz);
            for (int lz = 0; lz < WATER_TILE_SIZE; ++lz) {
                for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                    int local = lx + 8 * (lz + 8 * ly);
                    uint32_t source_mass = chunk->mass[local];
                    WaterVelocity source_velocity =
                        water_accelerated_velocity(chunk->velocity[local], source_mass);
                    if (chunk->static_solid[local] || chunk->dynamic_solid[local]) {
                        chunk->scratch[local] = (uint16_t)source_mass;
                        chunk->velocity_scratch[local] = 0;
                        continue;
                    }
                    uint32_t down = 0, incoming = 0;
                    WaterChunk *down_owner = ly > 0 ? chunk : below_chunk;
                    int down_local = ly > 0 ? local - 64 : lx + 8 * (lz + 8 * 7);
                    if (source_velocity.y <= 0 && y > 0 && below_active && down_owner &&
                        !down_owner->static_solid[down_local] && !down_owner->dynamic_solid[down_local]) {
                        uint32_t capacity = WATER_MAX_MASS - down_owner->mass[down_local];
                        down = source_mass < capacity ? source_mass : capacity;
                    }
                    WaterChunk *up_owner = ly < 7 ? chunk : above_chunk;
                    int up_local = ly < 7 ? local + 64 : lx + 8 * lz;
                    WaterVelocity incoming_velocity = (WaterVelocity){ 0, 0, 0 };
                    if (y + 1 < WATER_SIZE_Y && above_active && up_owner &&
                        !up_owner->static_solid[up_local] && !up_owner->dynamic_solid[up_local]) {
                        uint32_t capacity = WATER_MAX_MASS - source_mass;
                        uint32_t above = up_owner->mass[up_local];
                        incoming_velocity = water_accelerated_velocity(
                            up_owner->velocity[up_local], above);
                        if (incoming_velocity.y <= 0)
                            incoming = above < capacity ? above : capacity;
                    }
                    uint32_t retained = source_mass - down;
                    uint32_t next_mass = retained + incoming;
                    WaterVelocity next_velocity = water_mix_velocity(
                        retained, source_velocity, incoming, incoming_velocity);
                    chunk->scratch[local] = (uint16_t)next_mass;
                    chunk->velocity_scratch[local] = next_mass
                        ? water_pack_velocity(next_velocity) : 0;
                }
            }
        }
    }
}

static void water_horizontal_tiles(int start, int end, int worker_id, void *user) {
    (void)worker_id; (void)user;
    for (int job = start; job < end; ++job) {
        int tile = waterSystem.active_tile_ids[job];
        int tx, ty, tz;
        water_tile_coords(tile, &tx, &ty, &tz);
        WaterChunk *chunk = waterSystem.chunks[tile];
        int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
        bool tile_changed = false;
        for (int ly = 0; ly < WATER_TILE_SIZE; ++ly) {
            int y = y0 + ly;
            for (int lz = 0; lz < WATER_TILE_SIZE; ++lz) {
                int z = z0 + lz;
                for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                    int local = lx + 8 * (lz + 8 * ly);
                    int x = x0 + lx;
                    int index = x + WATER_SIZE_X * (z + WATER_SIZE_Z * y);
                    uint32_t next = chunk->scratch[local];
                    uint32_t next_velocity = 0;
                    if (!chunk->static_solid[local] && !chunk->dynamic_solid[local]) {
                        int neighbors[5] = {
                            x + 1 < WATER_SIZE_X ? index + 1 : -1,
                            x > 0 ? index - 1 : -1,
                            z + 1 < WATER_SIZE_Z ? index + WATER_SIZE_X : -1,
                            z > 0 ? index - WATER_SIZE_X : -1,
                            y + 1 < WATER_SIZE_Y ? index + WATER_SIZE_X * WATER_SIZE_Z : -1
                        };
                        int below = y > 0 ? index - WATER_SIZE_X * WATER_SIZE_Z : -1;
                        uint32_t outgoing[5] = { 0 }, incoming[5] = { 0 };
                        uint32_t outgoing_sum = 0, incoming_sum = 0;
                        for (int n = 0; n < 5; ++n) if (neighbors[n] >= 0) {
                            outgoing[n] = water_directional_flow(index, neighbors[n], n);
                            outgoing_sum += outgoing[n];
                        }
                        static const int reverse[4] = {
                            WATER_DIRECTION_XM, WATER_DIRECTION_XP,
                            WATER_DIRECTION_ZM, WATER_DIRECTION_ZP
                        };
                        for (int n = 0; n < 4; ++n) if (neighbors[n] >= 0) {
                            incoming[n] = water_directional_flow(neighbors[n], index, reverse[n]);
                            incoming_sum += incoming[n];
                        }
                        if (below >= 0) {
                            incoming[4] = water_directional_flow(below, index, WATER_DIRECTION_UP);
                            incoming_sum += incoming[4];
                        }
                        uint32_t retained = chunk->scratch[local] - outgoing_sum;
                        next = retained + incoming_sum;
                        if (next) {
                            WaterVelocity source_velocity = water_effective_velocity_index(index);
                            int32_t momentum_x = (int32_t)retained * source_velocity.x;
                            int32_t momentum_y = (int32_t)retained * source_velocity.y;
                            int32_t momentum_z = (int32_t)retained * source_velocity.z;
                            for (int n = 0; n < 4; ++n) if (incoming[n]) {
                                WaterVelocity carried = water_carried_velocity(neighbors[n], reverse[n]);
                                momentum_x += (int32_t)incoming[n] * carried.x;
                                momentum_y += (int32_t)incoming[n] * carried.y;
                                momentum_z += (int32_t)incoming[n] * carried.z;
                            }
                            if (incoming[4]) {
                                WaterVelocity carried = water_carried_velocity(below, WATER_DIRECTION_UP);
                                momentum_x += (int32_t)incoming[4] * carried.x;
                                momentum_y += (int32_t)incoming[4] * carried.y;
                                momentum_z += (int32_t)incoming[4] * carried.z;
                            }
                            WaterVelocity velocity = {
                                (momentum_x / (int32_t)next) * 7 / 8,
                                (momentum_y / (int32_t)next) * 7 / 8,
                                (momentum_z / (int32_t)next) * 7 / 8
                            };
                            next_velocity = water_pack_velocity(velocity);
                        }
                    }
                    if (chunk->mass[local] != next ||
                        chunk->velocity[local] != next_velocity) tile_changed = true;
                    chunk->mass[local] = (uint16_t)next;
                    chunk->velocity[local] = next_velocity;
                }
            }
        }
        if (tile_changed) {
            waterSystem.tile_changed[tile] = 1;
        }
    }
}

static void water_sync_tile_buffers(int tile) {
    WaterChunk *chunk = waterSystem.chunks[tile];
    if (chunk) {
        memcpy(chunk->scratch, chunk->mass, sizeof(chunk->mass));
        memcpy(chunk->velocity_scratch, chunk->velocity, sizeof(chunk->velocity));
    }
}

static void water_update_tile_sleep(int still_increment) {
    for (int job = 0; job < waterSystem.active_tile_count; ++job) {
        int tile = waterSystem.active_tile_ids[job];
        if (waterSystem.tile_changed[tile]) {
            waterSystem.tile_still_steps[tile] = 0;
            int tx, ty, tz;
            water_tile_coords(tile, &tx, &ty, &tz);
            water_wake_tile_xyz(tx + 1, ty, tz);
            water_wake_tile_xyz(tx - 1, ty, tz);
            water_wake_tile_xyz(tx, ty + 1, tz);
            water_wake_tile_xyz(tx, ty - 1, tz);
            water_wake_tile_xyz(tx, ty, tz + 1);
            water_wake_tile_xyz(tx, ty, tz - 1);
        } else {
            int still = (int)waterSystem.tile_still_steps[tile] + still_increment;
            if (still > WATER_SLEEP_STEPS) still = WATER_SLEEP_STEPS;
            waterSystem.tile_still_steps[tile] = (uint8_t)still;
            if (still < WATER_SLEEP_STEPS) continue;
            waterSystem.tile_active[tile] = 0;
            waterSystem.tile_still_steps[tile] = WATER_SLEEP_STEPS;
            water_sync_tile_buffers(tile);
        }
    }
}

static void water_reclaim_empty_chunks(void) {
    for (int r = waterSystem.resident_tile_count - 1; r >= 0; --r) {
        int tile = waterSystem.resident_tile_ids[r];
        WaterChunk *chunk = waterSystem.chunks[tile];
        if (!chunk || waterSystem.tile_active[tile]) continue;
        bool wet = false;
        for (int i = 0; i < WATER_CHUNK_CELLS; ++i) wet |= chunk->mass[i] != 0;
        if (wet) continue;
        free(chunk);
        waterSystem.chunks[tile] = NULL;
        waterSystem.resident_tile_ids[r] =
            waterSystem.resident_tile_ids[--waterSystem.resident_tile_count];
        waterSystem.tile_still_steps[tile] = 0;
        waterSystem.tile_changed[tile] = 0;
        waterSystem.tile_wet_cells[tile] = 0;
        waterSystem.tile_moving_cells[tile] = 0;
        waterSystem.tile_max_component_speed[tile] = 0;
        waterSystem.tile_visible_mass[tile] = 0;
        waterSystem.tile_trapped_mass[tile] = 0;
        waterSystem.tile_stats_dirty[tile] = 0;
        waterSystem.tile_render_dirty[tile] = 1;
        waterSystem.render_dirty = true;
        if (++waterSystem.resident_generation == 0) waterSystem.resident_generation = 1;
    }
}

static void water_recompute_tile_stats(int tile) {
    WaterChunk *chunk = waterSystem.chunks[tile];
    if (!chunk) {
        waterSystem.tile_wet_cells[tile] = 0;
        waterSystem.tile_moving_cells[tile] = 0;
        waterSystem.tile_max_component_speed[tile] = 0;
        waterSystem.tile_visible_mass[tile] = 0;
        waterSystem.tile_trapped_mass[tile] = 0;
        waterSystem.tile_stats_dirty[tile] = 0;
        return;
    }
    int tx, ty, tz;
    water_tile_coords(tile, &tx, &ty, &tz);
    int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
    uint16_t wet = 0, moving = 0;
    uint8_t max_speed = 0;
    uint64_t visible = 0, trapped = 0;
    for (int ly = 0; ly < WATER_TILE_SIZE; ++ly)
        for (int lz = 0; lz < WATER_TILE_SIZE; ++lz) {
            int row = x0 + WATER_SIZE_X * ((z0 + lz) + WATER_SIZE_Z * (y0 + ly));
            for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                int i = row + lx;
                int local = lx + 8 * (lz + 8 * ly);
                uint32_t mass = chunk->mass[local];
                if (!mass) continue;
                WaterVelocity velocity = water_unpack_velocity(chunk->velocity[local]);
                int speed_x = velocity.x < 0 ? -velocity.x : velocity.x;
                int speed_y = velocity.y < 0 ? -velocity.y : velocity.y;
                int speed_z = velocity.z < 0 ? -velocity.z : velocity.z;
                int speed = speed_x > speed_y ? speed_x : speed_y;
                if (speed_z > speed) speed = speed_z;
                if (speed > 0) moving++;
                if (speed > max_speed) max_speed = (uint8_t)speed;
                if (water_cell_blocked_index(i)) trapped += mass;
                else { wet++; visible += mass; }
            }
        }
    waterSystem.tile_wet_cells[tile] = wet;
    waterSystem.tile_moving_cells[tile] = moving;
    waterSystem.tile_max_component_speed[tile] = max_speed;
    waterSystem.tile_visible_mass[tile] = visible;
    waterSystem.tile_trapped_mass[tile] = trapped;
    waterSystem.tile_stats_dirty[tile] = 0;
}

static void water_refresh_diagnostics(void) {
    uint32_t wet = 0, active = 0, moving = 0, max_speed = 0;
    uint64_t visible = 0, trapped = 0;
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r];
        if (waterSystem.tile_stats_dirty[tile]) water_recompute_tile_stats(tile);
        active += waterSystem.tile_active[tile] != 0;
        wet += waterSystem.tile_wet_cells[tile];
        moving += waterSystem.tile_moving_cells[tile];
        if (waterSystem.tile_max_component_speed[tile] > max_speed)
            max_speed = waterSystem.tile_max_component_speed[tile];
        visible += waterSystem.tile_visible_mass[tile];
        trapped += waterSystem.tile_trapped_mass[tile];
    }
    waterSystem.diagnostics.active_tiles = active;
    waterSystem.diagnostics.active_cells = active * (WATER_TILE_SIZE * WATER_TILE_SIZE * WATER_TILE_SIZE);
    waterSystem.diagnostics.wet_cells = wet;
    waterSystem.diagnostics.moving_cells = moving;
    waterSystem.diagnostics.max_component_speed = max_speed;
    waterSystem.diagnostics.visible_mass = visible;
    waterSystem.diagnostics.trapped_mass = trapped;
    waterSystem.diagnostics.resident_chunks = waterSystem.resident_tile_count;
    waterSystem.diagnostics.resident_bytes =
        (uint64_t)waterSystem.resident_tile_count * (uint64_t)sizeof(WaterChunk);
    waterSystem.diagnostics.velocity_bytes =
        (uint64_t)waterSystem.resident_tile_count * WATER_CHUNK_CELLS *
        (sizeof(uint32_t) * 3u);
}

static void water_step_batch(float fixed_dt, int fixed_steps) {
    (void)fixed_dt;
    if (!waterSystem.initialized || !water_runtime_enabled() || fixed_steps <= 0) return;
    if (fixed_steps > PBD_MAX_ACCUM_STEPS) fixed_steps = PBD_MAX_ACCUM_STEPS;
    waterSystem.diagnostics.last_gpu_upload_ms = 0.0;
    waterSystem.diagnostics.last_gpu_dispatch_ms = 0.0;
    waterSystem.diagnostics.last_gpu_readback_ms = 0.0;
    double started = pbd_time_now_ms();
    water_rebuild_obstacles_and_displace();
    double after_obstacles = pbd_time_now_ms();
    water_build_active_tile_list();
    if (waterSystem.active_tile_count == 0) {
        water_reclaim_empty_chunks();
        water_refresh_diagnostics();
        double finished = pbd_time_now_ms();
        waterSystem.diagnostics.processed_cells = 0;
        waterSystem.diagnostics.last_obstacle_ms = after_obstacles - started;
        waterSystem.diagnostics.last_solver_ms = 0.0;
        waterSystem.diagnostics.last_maintenance_ms = finished - after_obstacles;
        waterSystem.diagnostics.last_step_ms = finished - started;
        return;
    }
    bool gpu_requested = physics_backend_is_gpu(physicsBackend.active);
    if (gpu_requested) for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        WaterChunk *chunk = waterSystem.chunks[waterSystem.resident_tile_ids[r]];
        memcpy(chunk->batch_snapshot, chunk->mass, sizeof(chunk->mass));
        memcpy(chunk->velocity_batch_snapshot, chunk->velocity, sizeof(chunk->velocity));
    }
    if (gpu_requested && water_gpu_run_batch(fixed_steps)) {
        for (int job = 0; job < waterSystem.active_tile_count; ++job) {
            int tile = waterSystem.active_tile_ids[job];
            int tx, ty, tz;
            water_tile_coords(tile, &tx, &ty, &tz);
            WaterChunk *chunk = waterSystem.chunks[tile];
            int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
            for (int ly = 0; ly < WATER_TILE_SIZE && !waterSystem.tile_changed[tile]; ++ly)
                for (int lz = 0; lz < WATER_TILE_SIZE && !waterSystem.tile_changed[tile]; ++lz) {
                    int row = x0 + WATER_SIZE_X * ((z0 + lz) + WATER_SIZE_Z * (y0 + ly));
                    for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                        int local = lx + 8 * (lz + 8 * ly);
                        if (chunk->mass[local] != chunk->batch_snapshot[local] ||
                            chunk->velocity[local] != chunk->velocity_batch_snapshot[local]) {
                            waterSystem.tile_changed[tile] = 1;
                            break;
                        }
                    }
                }
            if (waterSystem.tile_changed[tile]) water_mark_tile_and_neighbors_dirty(tile);
        }
        double after_solver = pbd_time_now_ms();
        water_update_tile_sleep(fixed_steps);
        water_reclaim_empty_chunks();
        water_refresh_diagnostics();
        double finished = pbd_time_now_ms();
        waterSystem.diagnostics.processed_cells =
            (uint32_t)waterSystem.active_tile_count * 512u * (uint32_t)fixed_steps;
        waterSystem.diagnostics.last_obstacle_ms = after_obstacles - started;
        waterSystem.diagnostics.last_solver_ms = after_solver - after_obstacles;
        waterSystem.diagnostics.last_maintenance_ms = finished - after_solver;
        waterSystem.diagnostics.last_step_ms = finished - started;
        return;
    }
    if (gpu_requested) {
        for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
            WaterChunk *chunk = waterSystem.chunks[waterSystem.resident_tile_ids[r]];
            memcpy(chunk->mass, chunk->batch_snapshot, sizeof(chunk->mass));
            memcpy(chunk->scratch, chunk->batch_snapshot, sizeof(chunk->mass));
            memcpy(chunk->velocity, chunk->velocity_batch_snapshot, sizeof(chunk->velocity));
            memcpy(chunk->velocity_scratch, chunk->velocity_batch_snapshot,
                   sizeof(chunk->velocity));
        }
        waterSystem.cpu_scratch_valid = true;
    }
    bool old_force_st = physics_force_single_cpu;
    physics_force_single_cpu = physicsBackend.active == PHYSICS_BACKEND_CPU_ST;
    waterSystem.diagnostics.backend = physics_force_single_cpu ? "cpu-st" : "cpu-mt";
    if (!waterSystem.cpu_scratch_valid) {
        for (int job = 0; job < waterSystem.active_tile_count; ++job)
            water_sync_tile_buffers(waterSystem.active_tile_ids[job]);
        waterSystem.cpu_scratch_valid = true;
    }
    for (int step = 0; step < fixed_steps; ++step) {
        pbd_parallel_for(0, waterSystem.active_tile_count, water_gravity_tiles, NULL);
        pbd_parallel_for(0, waterSystem.active_tile_count, water_horizontal_tiles, NULL);
    }
    double after_solver = pbd_time_now_ms();
    for (int job = 0; job < waterSystem.active_tile_count; ++job) {
        int tile = waterSystem.active_tile_ids[job];
        if (waterSystem.tile_changed[tile]) water_mark_tile_and_neighbors_dirty(tile);
    }
    bool any_mass_changed = false;
    for (int job = 0; job < waterSystem.active_tile_count; ++job)
        any_mass_changed |= waterSystem.tile_changed[waterSystem.active_tile_ids[job]] != 0;
    if (any_mass_changed) water_note_cpu_mass_edit();
    water_update_tile_sleep(fixed_steps);
    physics_force_single_cpu = old_force_st;
    water_reclaim_empty_chunks();
    water_refresh_diagnostics();
    double finished = pbd_time_now_ms();
    waterSystem.diagnostics.processed_cells =
        (uint32_t)waterSystem.active_tile_count * 512u * (uint32_t)fixed_steps;
    waterSystem.diagnostics.last_obstacle_ms = after_obstacles - started;
    waterSystem.diagnostics.last_solver_ms = after_solver - after_obstacles;
    waterSystem.diagnostics.last_maintenance_ms = finished - after_solver;
    waterSystem.diagnostics.last_step_ms = finished - started;
}

static void water_disturb(Vector3 world_position, Vector3 direction, uint32_t limit) {
    if (!waterSystem.initialized || !water_runtime_enabled() || limit == 0) return;
    int gx = (int)floorf(world_position.x / VOXEL_SIZE);
    int gy = (int)floorf(world_position.y / VOXEL_SIZE);
    int gz = (int)floorf(world_position.z / VOXEL_SIZE);
    if (!water_in_domain(gx, gy, gz)) return;
    int source = water_index_unchecked(gx, gy, gz);
    if (water_cell_blocked_index(source)) return;
    uint32_t amount = water_mass_index(source);
    if (amount > limit) amount = limit;
    if (!amount) return;
    int tx = gx, ty = gy + 1, tz = gz;
    if (fabsf(direction.x) >= fabsf(direction.z) && fabsf(direction.x) > 0.05f) tx += direction.x > 0 ? 1 : -1;
    else if (fabsf(direction.z) > 0.05f) tz += direction.z > 0 ? 1 : -1;
    if (!water_in_domain(tx, ty, tz)) return;
    int target = water_index_unchecked(tx, ty, tz);
    if (water_cell_blocked_index(target)) return;
    uint32_t capacity = WATER_MAX_MASS - water_mass_index(target);
    if (amount > capacity) amount = capacity;
    int source_local, target_local;
    WaterChunk *source_chunk = water_chunk_for_index(source, false, &source_local);
    WaterChunk *target_chunk = water_chunk_for_index(target, true, &target_local);
    if (!source_chunk || !target_chunk) return;
    WaterVelocity source_velocity = water_unpack_velocity(source_chunk->velocity[source_local]);
    WaterVelocity target_velocity = water_unpack_velocity(target_chunk->velocity[target_local]);
    WaterVelocity mixed_velocity = water_mix_velocity(target_chunk->mass[target_local],
                                                       target_velocity, amount,
                                                       source_velocity);
    source_chunk->mass[source_local] -= (uint16_t)amount;
    source_chunk->scratch[source_local] = source_chunk->mass[source_local];
    if (!source_chunk->mass[source_local]) source_chunk->velocity[source_local] = 0;
    source_chunk->velocity_scratch[source_local] = source_chunk->velocity[source_local];
    target_chunk->mass[target_local] += (uint16_t)amount;
    target_chunk->scratch[target_local] = target_chunk->mass[target_local];
    target_chunk->velocity[target_local] = water_pack_velocity(mixed_velocity);
    target_chunk->velocity_scratch[target_local] = target_chunk->velocity[target_local];
    water_add_impulse_index(source, direction, 0);
    water_add_impulse_index(target, direction, 0);
    water_note_cpu_mass_edit();
    water_wake_cell_and_neighbors(gx, gy, gz);
    water_wake_cell_and_neighbors(tx, ty, tz);
    waterSystem.render_dirty = true;
}

static int water_write_map(FILE *fp) {
    if (!fp) return 0;
    if (!water_runtime_enabled()) {
        fprintf(fp, "WATER 0\n");
        return 0;
    }
    int count = 0;
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        WaterChunk *chunk = waterSystem.chunks[waterSystem.resident_tile_ids[r]];
        for (int local = 0; local < 512; ++local) count += chunk->mass[local] > 0;
    }
    fprintf(fp, "WATER %d\n", count);
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r], tx, ty, tz;
        WaterChunk *chunk = waterSystem.chunks[tile];
        water_tile_coords(tile, &tx, &ty, &tz);
        for (int local = 0; local < 512; ++local) {
            if (!chunk->mass[local]) continue;
            int lx = local & 7, lz = (local >> 3) & 7, ly = local >> 6;
            int gx = WATER_MIN_X + tx * 8 + lx;
            int gy = WATER_MIN_Y + ty * 8 + ly;
            int gz = WATER_MIN_Z + tz * 8 + lz;
            fprintf(fp, "%d %d %d %u\n", gx, gy, gz, (unsigned)chunk->mass[local]);
        }
    }
    return count;
}

static bool water_read_map_entries(FILE *fp, int count) {
    if (!fp || count < 0 || count > WATER_CELL_COUNT) return false;
    water_reset();
    bool enabled = water_runtime_enabled();
    for (int i = 0; i < count; ++i) {
        int gx, gy, gz;
        unsigned int mass;
        if (fscanf(fp, "%d %d %d %u\n", &gx, &gy, &gz, &mass) != 4) return false;
        if (!water_in_domain(gx, gy, gz) || mass > WATER_MAX_MASS) return false;
        if (enabled) {
            (void)water_set_mass(gx, gy, gz, mass);
        }
    }
    waterSystem.render_dirty = true;
    water_refresh_diagnostics();
    return true;
}

static const WaterDiagnostics *water_diagnostics(void) {
    water_refresh_diagnostics();
    return &waterSystem.diagnostics;
}
