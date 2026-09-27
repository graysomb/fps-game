/* Sparse 8^3 water chunks. Voxels remain a read-only obstacle source. */
typedef struct WaterChunk {
    uint16_t mass[WATER_CHUNK_CELLS];
    uint16_t scratch[WATER_CHUNK_CELLS];
    uint16_t batch_snapshot[WATER_CHUNK_CELLS];
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
    if (chunk->mass[local] == mass && chunk->scratch[local] == mass) return true;
    chunk->mass[local] = (uint16_t)mass;
    chunk->scratch[local] = (uint16_t)mass;
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
                    if (chunk) chunk->dynamic_solid[local] = 1;
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
                chunk->mass[local] = (uint16_t)(current + amount);
                chunk->scratch[local] = chunk->mass[local];
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

static inline uint32_t water_flow_chunks(const WaterChunk *from_chunk, int from,
                                         const WaterChunk *to_chunk, int to,
                                         bool edge_active, bool from_supported) {
    if (!edge_active || !from_chunk || !to_chunk ||
        !from_supported ||
        from_chunk->static_solid[from] || from_chunk->dynamic_solid[from] ||
        to_chunk->static_solid[to] || to_chunk->dynamic_solid[to]) return 0;
    uint32_t a = from_chunk->scratch[from], b = to_chunk->scratch[to];
    return a > b ? (a - b) >> 3 : 0;
}

static inline bool water_laterally_supported(const WaterChunk *below_chunk,
                                              int below_local, int y) {
    if (y == 0) return true;
    if (!below_chunk) return false;
    if (below_chunk->static_solid[below_local] || below_chunk->dynamic_solid[below_local])
        return true;
    return below_chunk->scratch[below_local] >= WATER_LATERAL_SUPPORT_MASS;
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
                    if (chunk->static_solid[local] || chunk->dynamic_solid[local]) {
                        chunk->scratch[local] = (uint16_t)source_mass;
                        continue;
                    }
                    uint32_t down = 0, incoming = 0;
                    WaterChunk *down_owner = ly > 0 ? chunk : below_chunk;
                    int down_local = ly > 0 ? local - 64 : lx + 8 * (lz + 8 * 7);
                    if (y > 0 && below_active && down_owner &&
                        !down_owner->static_solid[down_local] && !down_owner->dynamic_solid[down_local]) {
                        uint32_t capacity = WATER_MAX_MASS - down_owner->mass[down_local];
                        down = source_mass < capacity ? source_mass : capacity;
                    }
                    WaterChunk *up_owner = ly < 7 ? chunk : above_chunk;
                    int up_local = ly < 7 ? local + 64 : lx + 8 * lz;
                    if (y + 1 < WATER_SIZE_Y && above_active && up_owner &&
                        !up_owner->static_solid[up_local] && !up_owner->dynamic_solid[up_local]) {
                        uint32_t capacity = WATER_MAX_MASS - source_mass;
                        uint32_t above = up_owner->mass[up_local];
                        incoming = above < capacity ? above : capacity;
                    }
                    chunk->scratch[local] = (uint16_t)(source_mass - down + incoming);
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
        WaterChunk *xp = water_neighbor_chunk(tx + 1, ty, tz);
        WaterChunk *xm = water_neighbor_chunk(tx - 1, ty, tz);
        WaterChunk *zp = water_neighbor_chunk(tx, ty, tz + 1);
        WaterChunk *zm = water_neighbor_chunk(tx, ty, tz - 1);
        WaterChunk *below = water_neighbor_chunk(tx, ty - 1, tz);
        WaterChunk *xp_below = water_neighbor_chunk(tx + 1, ty - 1, tz);
        WaterChunk *xm_below = water_neighbor_chunk(tx - 1, ty - 1, tz);
        WaterChunk *zp_below = water_neighbor_chunk(tx, ty - 1, tz + 1);
        WaterChunk *zm_below = water_neighbor_chunk(tx, ty - 1, tz - 1);
        int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
        bool tile_changed = false;
        for (int ly = 0; ly < WATER_TILE_SIZE; ++ly) {
            int y = y0 + ly;
            for (int lz = 0; lz < WATER_TILE_SIZE; ++lz) {
                int z = z0 + lz;
                for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                    int local = lx + 8 * (lz + 8 * ly);
                    uint32_t next = chunk->scratch[local];
                    if (!chunk->static_solid[local] && !chunk->dynamic_solid[local]) {
                        WaterChunk *owners[4] = {
                            lx < 7 ? chunk : xp, lx > 0 ? chunk : xm,
                            lz < 7 ? chunk : zp, lz > 0 ? chunk : zm
                        };
                        int neighbors[4] = {
                            lx < 7 ? local + 1 : 8 * (lz + 8 * ly),
                            lx > 0 ? local - 1 : 7 + 8 * (lz + 8 * ly),
                            lz < 7 ? local + 8 : lx + 8 * (8 * ly),
                            lz > 0 ? local - 8 : lx + 8 * (7 + 8 * ly)
                        };
                        bool enabled[4] = {
                            lx + 1 < WATER_TILE_SIZE || water_neighbor_tile_active(tx + 1, ty, tz),
                            lx > 0 || water_neighbor_tile_active(tx - 1, ty, tz),
                            lz + 1 < WATER_TILE_SIZE || water_neighbor_tile_active(tx, ty, tz + 1),
                            lz > 0 || water_neighbor_tile_active(tx, ty, tz - 1)
                        };
                        WaterChunk *below_owners[4] = {
                            ly > 0 ? owners[0] : (lx < 7 ? below : xp_below),
                            ly > 0 ? owners[1] : (lx > 0 ? below : xm_below),
                            ly > 0 ? owners[2] : (lz < 7 ? below : zp_below),
                            ly > 0 ? owners[3] : (lz > 0 ? below : zm_below)
                        };
                        int below_local = ly > 0 ? local - 64 : local + 448;
                        int neighbor_below[4] = {
                            ly > 0 ? neighbors[0] - 64 : neighbors[0] + 448,
                            ly > 0 ? neighbors[1] - 64 : neighbors[1] + 448,
                            ly > 0 ? neighbors[2] - 64 : neighbors[2] + 448,
                            ly > 0 ? neighbors[3] - 64 : neighbors[3] + 448
                        };
                        WaterChunk *source_below = ly > 0 ? chunk : below;
                        bool source_supported = water_laterally_supported(source_below,
                                                                          below_local, y);
                        if (x0 + lx + 1 >= WATER_SIZE_X) enabled[0] = false;
                        if (x0 + lx <= 0) enabled[1] = false;
                        if (z + 1 >= WATER_SIZE_Z) enabled[2] = false;
                        if (z <= 0) enabled[3] = false;
                        uint32_t outgoing = 0, incoming = 0;
                        for (int n = 0; n < 4; ++n) {
                            if (!enabled[n]) continue;
                            bool neighbor_supported = water_laterally_supported(
                                below_owners[n], neighbor_below[n], y);
                            outgoing += water_flow_chunks(chunk, local, owners[n], neighbors[n],
                                                          true, source_supported);
                            incoming += water_flow_chunks(owners[n], neighbors[n], chunk, local,
                                                          true, neighbor_supported);
                        }
                        next = chunk->scratch[local] - outgoing + incoming;
                    }
                    if (chunk->mass[local] != next) tile_changed = true;
                    chunk->mass[local] = (uint16_t)next;
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
    if (chunk) memcpy(chunk->scratch, chunk->mass, sizeof(chunk->mass));
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
        waterSystem.tile_visible_mass[tile] = 0;
        waterSystem.tile_trapped_mass[tile] = 0;
        waterSystem.tile_stats_dirty[tile] = 0;
        return;
    }
    int tx, ty, tz;
    water_tile_coords(tile, &tx, &ty, &tz);
    int x0 = tx * WATER_TILE_SIZE, y0 = ty * WATER_TILE_SIZE, z0 = tz * WATER_TILE_SIZE;
    uint16_t wet = 0;
    uint64_t visible = 0, trapped = 0;
    for (int ly = 0; ly < WATER_TILE_SIZE; ++ly)
        for (int lz = 0; lz < WATER_TILE_SIZE; ++lz) {
            int row = x0 + WATER_SIZE_X * ((z0 + lz) + WATER_SIZE_Z * (y0 + ly));
            for (int lx = 0; lx < WATER_TILE_SIZE; ++lx) {
                int i = row + lx;
                int local = lx + 8 * (lz + 8 * ly);
                uint32_t mass = chunk->mass[local];
                if (!mass) continue;
                if (water_cell_blocked_index(i)) trapped += mass;
                else { wet++; visible += mass; }
            }
        }
    waterSystem.tile_wet_cells[tile] = wet;
    waterSystem.tile_visible_mass[tile] = visible;
    waterSystem.tile_trapped_mass[tile] = trapped;
    waterSystem.tile_stats_dirty[tile] = 0;
}

static void water_refresh_diagnostics(void) {
    uint32_t wet = 0, active = 0;
    uint64_t visible = 0, trapped = 0;
    for (int r = 0; r < waterSystem.resident_tile_count; ++r) {
        int tile = waterSystem.resident_tile_ids[r];
        if (waterSystem.tile_stats_dirty[tile]) water_recompute_tile_stats(tile);
        active += waterSystem.tile_active[tile] != 0;
        wet += waterSystem.tile_wet_cells[tile];
        visible += waterSystem.tile_visible_mass[tile];
        trapped += waterSystem.tile_trapped_mass[tile];
    }
    waterSystem.diagnostics.active_tiles = active;
    waterSystem.diagnostics.active_cells = active * (WATER_TILE_SIZE * WATER_TILE_SIZE * WATER_TILE_SIZE);
    waterSystem.diagnostics.wet_cells = wet;
    waterSystem.diagnostics.visible_mass = visible;
    waterSystem.diagnostics.trapped_mass = trapped;
    waterSystem.diagnostics.resident_chunks = waterSystem.resident_tile_count;
    waterSystem.diagnostics.resident_bytes =
        (uint64_t)waterSystem.resident_tile_count * (uint64_t)sizeof(WaterChunk);
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
                        if (chunk->mass[local] != chunk->batch_snapshot[local]) {
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
    source_chunk->mass[source_local] -= (uint16_t)amount;
    source_chunk->scratch[source_local] = source_chunk->mass[source_local];
    target_chunk->mass[target_local] += (uint16_t)amount;
    target_chunk->scratch[target_local] = target_chunk->mass[target_local];
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
