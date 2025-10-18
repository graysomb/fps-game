#ifndef PBD_UNIFORM_GRID_H
#define PBD_UNIFORM_GRID_H

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "pbd_types.h"
#include "pbd_math.h"
#include "pbd_collisions.h"

typedef struct {
    int *counts;
    int *starts;
    int *content;
    int capacity_cells;
    int capacity_entries;
} PbdUniformGridCPU;

typedef PbdParticle *(*PbdParticleLookup)(int index, void *ctx);
typedef bool (*PbdSkipPairFn)(int particleA, int particleB, void *ctx);

static inline void pbd_system_configure(PbdSystemData *system,
                                        float dt,
                                        float omega_collision,
                                        Vector3 gravity,
                                        const Vector4 walls[6],
                                        int numParticles,
                                        int maxCellsPerElement)
{
    if (!system) return;
    system->dt = dt;
    system->omega_collision = omega_collision;
    system->numParticles = numParticles;
    system->maxCellsPerElement = maxCellsPerElement;
    system->gravity = (Vector4){ gravity.x, gravity.y, gravity.z, 0.0f };
    for (int i = 0; i < 6; ++i) {
        system->wallConstraints[i] = walls[i];
    }
    system->time += dt;
    system->compliance = 0.0f;
}

static inline bool pbd_uniform_grid_ensure_capacity(PbdUniformGridCPU *grid,
                                                    int numCells,
                                                    int numEntries)
{
    if (!grid) return false;

    if (numCells > grid->capacity_cells) {
        int *newCounts = (int *)realloc(grid->counts, (size_t)numCells * sizeof(int));
        if (!newCounts) return false;
        grid->counts = newCounts;

        int *newStarts = (int *)realloc(grid->starts, (size_t)numCells * sizeof(int));
        if (!newStarts) return false;
        grid->starts = newStarts;

        grid->capacity_cells = numCells;
    }

    if (numEntries > grid->capacity_entries) {
        int *newContent = (int *)realloc(grid->content, (size_t)numEntries * sizeof(int));
        if (!newContent) return false;
        grid->content = newContent;
        grid->capacity_entries = numEntries;
    }

    return true;
}

static inline bool pbd_uniform_grid_rebuild(PbdUniformGridCPU *grid,
                                            PbdSystemData *system,
                                            int particleCount,
                                            PbdParticleLookup lookup,
                                            void *lookupCtx,
                                            Vector3 gridMin,
                                            Vector3 gridMax,
                                            Vector3 cellSize)
{
    if (!grid || !system || !lookup || particleCount <= 0) {
        system->numCells = 0;
        return false;
    }

    int rx = (int)ceilf((gridMax.x - gridMin.x) / cellSize.x);
    int ry = (int)ceilf((gridMax.y - gridMin.y) / cellSize.y);
    int rz = (int)ceilf((gridMax.z - gridMin.z) / cellSize.z);
    if (rx < 1) rx = 1;
    if (ry < 1) ry = 1;
    if (rz < 1) rz = 1;

    const int numCells = rx * ry * rz;
    if (!pbd_uniform_grid_ensure_capacity(grid, numCells, particleCount)) {
        system->numCells = 0;
        return false;
    }

    memset(grid->counts, 0, (size_t)numCells * sizeof(int));
    memset(grid->starts, 0, (size_t)numCells * sizeof(int));

    system->gridMin = gridMin;
    system->gridMax = gridMax;
    system->cellSize = cellSize;
    system->gridRes[0] = rx;
    system->gridRes[1] = ry;
    system->gridRes[2] = rz;
    system->numCells = numCells;

    const float invCellX = 1.0f / cellSize.x;
    const float invCellY = 1.0f / cellSize.y;
    const float invCellZ = 1.0f / cellSize.z;

    for (int idx = 0; idx < particleCount; ++idx) {
        PbdParticle *p = lookup(idx, lookupCtx);
        if (!p) continue;
        Vector3 rel = pbd_v3_sub(p->predicted_pos, gridMin);
        int cx = (int)floorf(rel.x * invCellX);
        int cy = (int)floorf(rel.y * invCellY);
        int cz = (int)floorf(rel.z * invCellZ);
        cx = cx < 0 ? 0 : (cx >= rx ? rx - 1 : cx);
        cy = cy < 0 ? 0 : (cy >= ry ? ry - 1 : cy);
        cz = cz < 0 ? 0 : (cz >= rz ? rz - 1 : cz);
        int cell = cx + rx * (cy + ry * cz);
        grid->counts[cell]++;
    }

    int running = 0;
    for (int cell = 0; cell < numCells; ++cell) {
        grid->starts[cell] = running;
        running += grid->counts[cell];
        grid->counts[cell] = 0;
    }

    if (!pbd_uniform_grid_ensure_capacity(grid, numCells, running)) {
        system->numCells = 0;
        return false;
    }

    if (running > 0) {
        memset(grid->counts, 0, (size_t)numCells * sizeof(int));
        for (int idx = 0; idx < particleCount; ++idx) {
            PbdParticle *p = lookup(idx, lookupCtx);
            if (!p) continue;
            Vector3 rel = pbd_v3_sub(p->predicted_pos, gridMin);
            int cx = (int)floorf(rel.x * invCellX);
            int cy = (int)floorf(rel.y * invCellY);
            int cz = (int)floorf(rel.z * invCellZ);
            cx = cx < 0 ? 0 : (cx >= rx ? rx - 1 : cx);
            cy = cy < 0 ? 0 : (cy >= ry ? ry - 1 : cy);
            cz = cz < 0 ? 0 : (cz >= rz ? rz - 1 : cz);
            int cell = cx + rx * (cy + ry * cz);
            int start = grid->starts[cell];
            int offset = grid->counts[cell]++;
            grid->content[start + offset] = idx;
        }
    }

    return true;
}

static inline void pbd_project_collisions_with_grid(const PbdUniformGridCPU *grid,
                                                    const PbdSystemData *system,
                                                    int particleCount,
                                                    PbdParticleLookup lookup,
                                                    void *lookupCtx,
                                                    PbdSkipPairFn skipFn,
                                                    void *skipCtx)
{
    if (!grid || !system || !lookup || particleCount <= 0 || system->numCells <= 0) {
        return;
    }

    const int rx = system->gridRes[0];
    const int ry = system->gridRes[1];
    const int rz = system->gridRes[2];
    const Vector3 cellSize = system->cellSize;
    const Vector3 gridMin = system->gridMin;
    const float invCellX = 1.0f / cellSize.x;
    const float invCellY = 1.0f / cellSize.y;
    const float invCellZ = 1.0f / cellSize.z;
    const float eps = 1e-6f;

    for (int ix = 0; ix < particleCount; ++ix) {
        PbdParticle *pi = lookup(ix, lookupCtx);
        if (!pi) continue;
        float wi = pi->inv_mass;
        if (wi == 0.0f) continue;

        float ri = pi->radius;
        Vector3 dp = { 2.0f * ri, 2.0f * ri, 2.0f * ri };
        Vector3 minPos = pbd_v3_sub(pi->predicted_pos, dp);
        Vector3 maxPos = pbd_v3_add(pi->predicted_pos, dp);

        int cx_min = (int)floorf((minPos.x - gridMin.x) * invCellX);
        int cy_min = (int)floorf((minPos.y - gridMin.y) * invCellY);
        int cz_min = (int)floorf((minPos.z - gridMin.z) * invCellZ);
        int cx_max = (int)floorf((maxPos.x - gridMin.x) * invCellX);
        int cy_max = (int)floorf((maxPos.y - gridMin.y) * invCellY);
        int cz_max = (int)floorf((maxPos.z - gridMin.z) * invCellZ);

        cx_min = cx_min < 0 ? 0 : (cx_min >= rx ? rx - 1 : cx_min);
        cy_min = cy_min < 0 ? 0 : (cy_min >= ry ? ry - 1 : cy_min);
        cz_min = cz_min < 0 ? 0 : (cz_min >= rz ? rz - 1 : cz_min);
        cx_max = cx_max < 0 ? 0 : (cx_max >= rx ? rx - 1 : cx_max);
        cy_max = cy_max < 0 ? 0 : (cy_max >= ry ? ry - 1 : cy_max);
        cz_max = cz_max < 0 ? 0 : (cz_max >= rz ? rz - 1 : cz_max);

        Vector3 dx_total = { 0.0f, 0.0f, 0.0f };

        for (int cx = cx_min; cx <= cx_max; ++cx) {
            for (int cy = cy_min; cy <= cy_max; ++cy) {
                for (int cz = cz_min; cz <= cz_max; ++cz) {
                    int cell = cx + rx * (cy + ry * cz);
                    int start = grid->starts[cell];
                    int count = grid->counts[cell];
                    for (int t = 0; t < count; ++t) {
                        int jx = grid->content[start + t];
                        if (jx == ix) continue;
                        if (skipFn && skipFn(ix, jx, skipCtx)) continue;

                        PbdParticle *pj = lookup(jx, lookupCtx);
                        if (!pj) continue;
                        float wj = pj->inv_mass;
                        if (wj == 0.0f) continue;

                        float rj = pj->radius;
                        Vector3 nij = pbd_v3_sub(pi->predicted_pos, pj->predicted_pos);
                        float d = pbd_v3_length(nij) + eps;
                        float target = ri + rj;
                        if (d >= target) continue;

                        Vector3 uij = pbd_v3_scale(nij, 1.0f / d);
                        float pen = target - d;
                        if (pen <= 0.0f) continue;

                        float denom = wi + wj;
                        if (denom <= 0.0f) continue;

                        float h = 0.5f * pen;
                        float scale = wi * system->omega_collision * h / denom;
                        dx_total = pbd_v3_add(dx_total, pbd_v3_scale(uij, scale));
                    }
                }
            }
        }

        pi->predicted_pos = pbd_v3_add(pi->predicted_pos, dx_total);
    }
}

#endif // PBD_UNIFORM_GRID_H
