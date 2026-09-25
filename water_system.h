#ifndef FPS_WATER_SYSTEM_H
#define FPS_WATER_SYSTEM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#define WATER_MIN_X (-40)
#define WATER_MIN_Y 0
#define WATER_MIN_Z (-40)
#define WATER_SIZE_X 80
#define WATER_SIZE_Y 80
#define WATER_SIZE_Z 80
#define WATER_MAX_X (WATER_MIN_X + WATER_SIZE_X - 1)
#define WATER_MAX_Y (WATER_MIN_Y + WATER_SIZE_Y - 1)
#define WATER_MAX_Z (WATER_MIN_Z + WATER_SIZE_Z - 1)
#define WATER_CELL_COUNT (WATER_SIZE_X * WATER_SIZE_Y * WATER_SIZE_Z)
#define WATER_MAX_MASS UINT32_C(65535)
#define WATER_RENDER_MIN_MASS UINT32_C(256)
#define WATER_TILE_SIZE 8
#define WATER_TILE_COUNT_X (WATER_SIZE_X / WATER_TILE_SIZE)
#define WATER_TILE_COUNT_Y (WATER_SIZE_Y / WATER_TILE_SIZE)
#define WATER_TILE_COUNT_Z (WATER_SIZE_Z / WATER_TILE_SIZE)
#define WATER_TILE_COUNT (WATER_TILE_COUNT_X * WATER_TILE_COUNT_Y * WATER_TILE_COUNT_Z)
#define WATER_SLEEP_STEPS 30

typedef enum WaterCellKind {
    WORLD_CELL_EMPTY = 0,
    WORLD_CELL_SOLID = 1,
    WORLD_CELL_WATER = 2
} WaterCellKind;

typedef struct WorldCell {
    WaterCellKind kind;
    int voxel_index;
    uint32_t water_mass;
} WorldCell;

typedef struct WaterDiagnostics {
    const char *backend;
    uint32_t active_tiles;
    uint32_t active_cells;
    uint32_t wet_cells;
    uint64_t visible_mass;
    uint64_t trapped_mass;
    uint64_t displaced_mass;
    double last_step_ms;
} WaterDiagnostics;

static void water_init(void);
static void water_reset(void);
static void water_shutdown(void);
static bool water_in_domain(int gx, int gy, int gz);
static uint32_t water_get_mass(int gx, int gy, int gz);
static bool water_set_mass(int gx, int gy, int gz, uint32_t mass);
static void water_add_region(int minx, int maxx, int miny, int maxy,
                             int minz, int maxz, uint32_t mass);
static void water_remove_region(int minx, int maxx, int miny, int maxy,
                                int minz, int maxz);
static float water_sample_fill(Vector3 world_position);
static WorldCell world_cell_get(int gx, int gy, int gz);
static void water_step_batch(float fixed_dt, int fixed_steps);
static bool water_gpu_run_batch(int fixed_steps);
static void water_gpu_shutdown_backend(void);
static void water_rebuild_obstacles_and_displace(void);
static void water_disturb(Vector3 world_position, Vector3 direction, uint32_t limit);
static int water_write_map(FILE *fp);
static bool water_read_map_entries(FILE *fp, int count);
static const WaterDiagnostics *water_diagnostics(void);

#endif
