/*
 * Dense, finite-mass cellular-automata water.  This deliberately lives beside
 * the voxel/PBD storage: voxels are only consulted to build read-only masks.
 */
typedef struct WaterSystem {
    uint32_t mass[WATER_CELL_COUNT];
    uint32_t scratch[WATER_CELL_COUNT];
    uint32_t pre_batch[WATER_CELL_COUNT];
    uint32_t batch_snapshot[WATER_CELL_COUNT];
    uint8_t static_solid[WATER_CELL_COUNT];
    uint8_t dynamic_solid[WATER_CELL_COUNT];
    uint8_t previous_solid[WATER_CELL_COUNT];
    uint8_t tile_active[WATER_TILE_COUNT];
    uint8_t tile_still_steps[WATER_TILE_COUNT];
    uint32_t bfs_queue[WATER_CELL_COUNT];
    uint8_t bfs_depth[WATER_CELL_COUNT];
    uint32_t bfs_seen[WATER_CELL_COUNT];
    uint32_t bfs_stamp;
    uint64_t static_generation;
    bool initialized;
    bool render_dirty;
    WaterDiagnostics diagnostics;
} WaterSystem;

static WaterSystem waterSystem;

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
    waterSystem.tile_active[tile] = 1;
    waterSystem.tile_still_steps[tile] = 0;
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
    return waterSystem.static_solid[index] || waterSystem.dynamic_solid[index];
}

static void water_init(void) {
    if (waterSystem.initialized) return;
    memset(&waterSystem, 0, sizeof(waterSystem));
    waterSystem.initialized = true;
    waterSystem.static_generation = UINT64_MAX;
    waterSystem.render_dirty = true;
    waterSystem.diagnostics.backend = "cpu-st";
}

static void water_reset(void) {
    water_gpu_shutdown_backend();
    memset(&waterSystem, 0, sizeof(waterSystem));
    waterSystem.initialized = true;
    waterSystem.static_generation = UINT64_MAX;
    waterSystem.render_dirty = true;
    waterSystem.diagnostics.backend = "cpu-st";
}

static void water_shutdown(void) {
    memset(&waterSystem, 0, sizeof(waterSystem));
}

static uint32_t water_get_mass(int gx, int gy, int gz) {
    if (!waterSystem.initialized || !water_runtime_enabled() ||
        !water_in_domain(gx, gy, gz)) return 0;
    return waterSystem.mass[water_index_unchecked(gx, gy, gz)];
}

static bool water_set_mass(int gx, int gy, int gz, uint32_t mass) {
    if (!waterSystem.initialized) water_init();
    if (!water_runtime_enabled() || !water_in_domain(gx, gy, gz)) return false;
    int index = water_index_unchecked(gx, gy, gz);
    if (water_cell_blocked_index(index) || table_get(gx, gy, gz) >= 0) return false;
    if (mass > WATER_MAX_MASS) mass = WATER_MAX_MASS;
    if (waterSystem.mass[index] == mass && waterSystem.scratch[index] == mass) return true;
    waterSystem.mass[index] = mass;
    waterSystem.scratch[index] = mass;
    water_wake_cell_and_neighbors(gx, gy, gz);
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
                waterSystem.mass[index] = 0;
                waterSystem.scratch[index] = 0;
                water_wake_cell_and_neighbors(x, y, z);
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
    for (int gy = WATER_MIN_Y; gy <= WATER_MAX_Y; ++gy) {
        for (int gz = WATER_MIN_Z; gz <= WATER_MAX_Z; ++gz) {
            for (int gx = WATER_MIN_X; gx <= WATER_MAX_X; ++gx) {
                int index = water_index_unchecked(gx, gy, gz);
                int voxel = table_get_static_only(gx, gy, gz);
                waterSystem.static_solid[index] = (uint8_t)(voxel >= 0 && voxel < voxel_count && !voxels[voxel].simulate);
            }
        }
    }
    waterSystem.static_generation = staticHashGeneration;
}

static void water_rebuild_dynamic_mask(void) {
    memset(waterSystem.dynamic_solid, 0, sizeof(waterSystem.dynamic_solid));
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
                    waterSystem.dynamic_solid[water_index_unchecked(gx, gy, gz)] = 1;
    }
}

static uint32_t water_displace_from_cell(int source_index) {
    uint32_t remaining = waterSystem.mass[source_index];
    if (remaining == 0) return 0;
    if (++waterSystem.bfs_stamp == 0) {
        memset(waterSystem.bfs_seen, 0, sizeof(waterSystem.bfs_seen));
        waterSystem.bfs_stamp = 1;
    }
    uint32_t head = 0, tail = 0;
    waterSystem.bfs_queue[tail] = (uint32_t)source_index;
    waterSystem.bfs_depth[tail++] = 0;
    waterSystem.bfs_seen[source_index] = waterSystem.bfs_stamp;
    uint32_t moved = 0;
    static const int dx[6] = { 0, 1, -1, 0, 0, 0 };
    static const int dy[6] = { 1, 0, 0, 0, 0, -1 };
    static const int dz[6] = { 0, 0, 0, 1, -1, 0 };
    while (head < tail && remaining > 0) {
        int index = (int)waterSystem.bfs_queue[head];
        int depth = waterSystem.bfs_depth[head++];
        int gx, gy, gz;
        water_coords_from_index(index, &gx, &gy, &gz);
        if (index != source_index && !water_cell_blocked_index(index)) {
            uint32_t current = waterSystem.mass[index];
            uint32_t capacity = WATER_MAX_MASS - current;
            uint32_t amount = remaining < capacity ? remaining : capacity;
            if (amount > 0) {
                waterSystem.mass[index] = current + amount;
                waterSystem.scratch[index] = current + amount;
                remaining -= amount;
                moved += amount;
                water_wake_cell_and_neighbors(gx, gy, gz);
            }
        }
        if (depth >= 8) continue;
        for (int d = 0; d < 6; ++d) {
            int nx = gx + dx[d], ny = gy + dy[d], nz = gz + dz[d];
            if (!water_in_domain(nx, ny, nz)) continue;
            int next = water_index_unchecked(nx, ny, nz);
            if (waterSystem.bfs_seen[next] == waterSystem.bfs_stamp) continue;
            waterSystem.bfs_seen[next] = waterSystem.bfs_stamp;
            waterSystem.bfs_queue[tail] = (uint32_t)next;
            waterSystem.bfs_depth[tail++] = (uint8_t)(depth + 1);
        }
    }
    waterSystem.mass[source_index] = remaining;
    waterSystem.scratch[source_index] = remaining;
    return moved;
}

static void water_rebuild_obstacles_and_displace(void) {
    if (!waterSystem.initialized || !water_runtime_enabled()) return;
    water_rebuild_static_mask();
    water_rebuild_dynamic_mask();
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        uint8_t blocked = (uint8_t)water_cell_blocked_index(i);
        if (blocked && !waterSystem.previous_solid[i] && waterSystem.mass[i] > 0) {
            waterSystem.diagnostics.displaced_mass += water_displace_from_cell(i);
        } else if (!blocked && waterSystem.previous_solid[i]) {
            int gx, gy, gz;
            water_coords_from_index(i, &gx, &gy, &gz);
            water_wake_cell_and_neighbors(gx, gy, gz);
        }
        waterSystem.previous_solid[i] = blocked;
    }
    waterSystem.render_dirty = true;
}

static inline uint32_t water_down_flow(const uint32_t *source, int index) {
    if (water_cell_blocked_index(index)) return 0;
    int gx, gy, gz;
    water_coords_from_index(index, &gx, &gy, &gz);
    if (gy <= WATER_MIN_Y) return 0;
    int below = index - WATER_SIZE_X * WATER_SIZE_Z;
    if (water_cell_blocked_index(below)) return 0;
    uint32_t available = WATER_MAX_MASS - source[below];
    return source[index] < available ? source[index] : available;
}

static void water_gravity_range(int start, int end, int worker_id, void *user) {
    (void)worker_id;
    (void)user;
    const uint32_t *source = waterSystem.mass;
    uint32_t *dest = waterSystem.scratch;
    for (int index = start; index < end; ++index) {
        int gx, gy, gz;
        water_coords_from_index(index, &gx, &gy, &gz);
        int tile = water_tile_index_for_cell(gx, gy, gz);
        if (!waterSystem.tile_active[tile] || water_cell_blocked_index(index)) {
            dest[index] = source[index];
            continue;
        }
        uint32_t down = water_down_flow(source, index);
        uint32_t incoming = 0;
        if (gy < WATER_MAX_Y) {
            int above = index + WATER_SIZE_X * WATER_SIZE_Z;
            incoming = water_down_flow(source, above);
        }
        dest[index] = source[index] - down + incoming;
    }
}

static inline uint32_t water_horizontal_flow(const uint32_t *source,
                                             int from, int to) {
    if (to < 0 || to >= WATER_CELL_COUNT ||
        water_cell_blocked_index(from) || water_cell_blocked_index(to)) return 0;
    uint32_t a = source[from], b = source[to];
    return a > b ? (a - b) >> 3 : 0;
}

static void water_horizontal_range(int start, int end, int worker_id, void *user) {
    (void)worker_id;
    (void)user;
    const uint32_t *source = waterSystem.scratch;
    uint32_t *dest = waterSystem.mass;
    for (int index = start; index < end; ++index) {
        int gx, gy, gz;
        water_coords_from_index(index, &gx, &gy, &gz);
        int tile = water_tile_index_for_cell(gx, gy, gz);
        if (!waterSystem.tile_active[tile] || water_cell_blocked_index(index)) {
            dest[index] = source[index];
            continue;
        }
        int neighbor[4] = {
            gx < WATER_MAX_X ? index + 1 : -1,
            gx > WATER_MIN_X ? index - 1 : -1,
            gz < WATER_MAX_Z ? index + WATER_SIZE_X : -1,
            gz > WATER_MIN_Z ? index - WATER_SIZE_X : -1
        };
        uint32_t outgoing = 0, incoming = 0;
        for (int n = 0; n < 4; ++n) {
            int j = neighbor[n];
            if (j < 0) continue;
            outgoing += water_horizontal_flow(source, index, j);
            incoming += water_horizontal_flow(source, j, index);
        }
        dest[index] = source[index] - outgoing + incoming;
    }
}

static void water_update_tile_sleep(int still_increment) {
    uint8_t changed_tiles[WATER_TILE_COUNT] = { 0 };
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        int gx, gy, gz;
        water_coords_from_index(i, &gx, &gy, &gz);
        int tile = water_tile_index_for_cell(gx, gy, gz);
        if (waterSystem.tile_active[tile] && waterSystem.mass[i] != waterSystem.pre_batch[i]) {
            changed_tiles[tile] = 1;
        }
    }
    for (int tile = 0; tile < WATER_TILE_COUNT; ++tile) {
        if (!waterSystem.tile_active[tile]) continue;
        if (changed_tiles[tile]) {
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
        }
    }
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        int gx, gy, gz;
        water_coords_from_index(i, &gx, &gy, &gz);
        int tile = water_tile_index_for_cell(gx, gy, gz);
        if (!waterSystem.tile_active[tile]) waterSystem.scratch[i] = waterSystem.mass[i];
    }
}

static void water_refresh_diagnostics(void) {
    uint32_t wet = 0, active = 0;
    uint64_t visible = 0, trapped = 0;
    for (int tile = 0; tile < WATER_TILE_COUNT; ++tile) active += waterSystem.tile_active[tile] != 0;
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        uint32_t mass = waterSystem.mass[i];
        if (!mass) continue;
        if (water_cell_blocked_index(i)) trapped += mass;
        else { wet++; visible += mass; }
    }
    waterSystem.diagnostics.active_tiles = active;
    waterSystem.diagnostics.active_cells = active * (WATER_TILE_SIZE * WATER_TILE_SIZE * WATER_TILE_SIZE);
    waterSystem.diagnostics.wet_cells = wet;
    waterSystem.diagnostics.visible_mass = visible;
    waterSystem.diagnostics.trapped_mass = trapped;
}

static void water_step_batch(float fixed_dt, int fixed_steps) {
    (void)fixed_dt;
    if (!waterSystem.initialized || !water_runtime_enabled() || fixed_steps <= 0) return;
    if (fixed_steps > PBD_MAX_ACCUM_STEPS) fixed_steps = PBD_MAX_ACCUM_STEPS;
    double started = pbd_time_now_ms();
    water_rebuild_obstacles_and_displace();
    memcpy(waterSystem.batch_snapshot, waterSystem.mass, sizeof(waterSystem.mass));
    if (physics_backend_is_gpu(physicsBackend.active) && water_gpu_run_batch(fixed_steps)) {
        memcpy(waterSystem.pre_batch, waterSystem.batch_snapshot, sizeof(waterSystem.mass));
        water_update_tile_sleep(fixed_steps);
        waterSystem.render_dirty = true;
        water_refresh_diagnostics();
        waterSystem.diagnostics.last_step_ms = pbd_time_now_ms() - started;
        return;
    }
    if (physics_backend_is_gpu(physicsBackend.active)) {
        memcpy(waterSystem.mass, waterSystem.batch_snapshot, sizeof(waterSystem.mass));
        memcpy(waterSystem.scratch, waterSystem.batch_snapshot, sizeof(waterSystem.mass));
    }
    bool old_force_st = physics_force_single_cpu;
    physics_force_single_cpu = physicsBackend.active == PHYSICS_BACKEND_CPU_ST;
    waterSystem.diagnostics.backend = physics_force_single_cpu ? "cpu-st" : "cpu-mt";
    for (int step = 0; step < fixed_steps; ++step) {
        memcpy(waterSystem.pre_batch, waterSystem.mass, sizeof(waterSystem.mass));
        pbd_parallel_for(0, WATER_CELL_COUNT, water_gravity_range, NULL);
        pbd_parallel_for(0, WATER_CELL_COUNT, water_horizontal_range, NULL);
        water_update_tile_sleep(1);
    }
    physics_force_single_cpu = old_force_st;
    waterSystem.render_dirty = true;
    water_refresh_diagnostics();
    waterSystem.diagnostics.last_step_ms = pbd_time_now_ms() - started;
}

static void water_disturb(Vector3 world_position, Vector3 direction, uint32_t limit) {
    if (!waterSystem.initialized || !water_runtime_enabled() || limit == 0) return;
    int gx = (int)floorf(world_position.x / VOXEL_SIZE);
    int gy = (int)floorf(world_position.y / VOXEL_SIZE);
    int gz = (int)floorf(world_position.z / VOXEL_SIZE);
    if (!water_in_domain(gx, gy, gz)) return;
    int source = water_index_unchecked(gx, gy, gz);
    if (water_cell_blocked_index(source)) return;
    uint32_t amount = waterSystem.mass[source];
    if (amount > limit) amount = limit;
    if (!amount) return;
    int tx = gx, ty = gy + 1, tz = gz;
    if (fabsf(direction.x) >= fabsf(direction.z) && fabsf(direction.x) > 0.05f) tx += direction.x > 0 ? 1 : -1;
    else if (fabsf(direction.z) > 0.05f) tz += direction.z > 0 ? 1 : -1;
    if (!water_in_domain(tx, ty, tz)) return;
    int target = water_index_unchecked(tx, ty, tz);
    if (water_cell_blocked_index(target)) return;
    uint32_t capacity = WATER_MAX_MASS - waterSystem.mass[target];
    if (amount > capacity) amount = capacity;
    waterSystem.mass[source] -= amount;
    waterSystem.scratch[source] = waterSystem.mass[source];
    waterSystem.mass[target] += amount;
    waterSystem.scratch[target] = waterSystem.mass[target];
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
    for (int i = 0; i < WATER_CELL_COUNT; ++i) if (waterSystem.mass[i] > 0) count++;
    fprintf(fp, "WATER %d\n", count);
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        if (waterSystem.mass[i] == 0) continue;
        int gx, gy, gz;
        water_coords_from_index(i, &gx, &gy, &gz);
        fprintf(fp, "%d %d %d %u\n", gx, gy, gz, waterSystem.mass[i]);
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
            int index = water_index_unchecked(gx, gy, gz);
            waterSystem.mass[index] = (uint32_t)mass;
            waterSystem.scratch[index] = (uint32_t)mass;
            water_wake_cell_and_neighbors(gx, gy, gz);
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
