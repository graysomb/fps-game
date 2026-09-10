#ifndef MEGALITH_RASTERIZER_H
#define MEGALITH_RASTERIZER_H

#include "megalith_grammar.h"
#include "raylib.h"
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer for placing a voxel at grid coordinate (gx, gy, gz)
typedef void (*MegaVoxelPlotFunc)(int gx, int gy, int gz, Color color);

// Fast 3D integer hash for weathered stone texture and procedural erosion
static inline uint32_t mega_hash3(int x, int y, int z) {
    uint32_t h = (uint32_t)(x * 73856093 ^ y * 19349663 ^ z * 83492791);
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

// Rasterizes an abstract MegalithPlan into discrete 3D voxels
static inline void rasterize_megalith_plan(const MegalithPlan *plan,
                                           int center_gx, int center_gz, int base_gy,
                                           MegaVoxelPlotFunc plot) {
    if (!plan || !plot) return;

    int mv = 2; // 2 voxels per modular meter

    // Prehistoric Geological Palette
    Color col_granite_light  = (Color){ 150, 145, 135, 255 }; // Weathered sarsen / granite
    Color col_granite_dark   = (Color){ 115, 110, 105, 255 }; // Crevice / cleft granite
    Color col_lichen_gold    = (Color){ 175, 160, 85, 255 };  // Ancient yellow lichen fleck
    Color col_moss_green     = (Color){ 90, 115, 75, 255 };   // Damp north-face moss
    Color col_turf_top       = (Color){ 65, 105, 50, 255 };   // Grassy mound mantle
    Color col_earth_subsoil  = (Color){ 105, 75, 50, 255 };   // Packed earth barrow fill
    Color col_floor_slab     = (Color){ 140, 135, 125, 255 }; // Passage flagstones
    Color col_ash            = (Color){ 85, 80, 80, 255 };    // Hearth ash bed
    Color col_ember          = (Color){ 255, 125, 20, 255 };  // Sacred offering fire

    // 1. First Pass: Rasterize Earth Tumulus Mound (Leaving Chamber Void Open)
    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        if (n->type != MEGALITH_PRIMITIVE_MOUND) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int rx = (gx1 - gx0) / 2;
        int rz = (gz1 - gz0) / 2;
        int max_h = (n->box.max_y - n->box.min_y) * mv;

        for (int z = gz0; z < gz1; ++z) {
            float nz = (float)(z - mid_z) / (float)rz;
            for (int x = gx0; x < gx1; ++x) {
                float nx = (float)(x - mid_x) / (float)rx;
                float d2 = nx * nx + nz * nz;
                if (d2 >= 1.0f) continue;

                // Parabolic dome height
                int dome_h = (int)((1.0f - d2) * max_h);
                if (dome_h < 1) dome_h = 1;

                // Keep chamber and passage void hollow!
                // Chamber radius is ~8 voxels around (center_gx, center_gz)
                int dx_c = x - center_gx;
                int dz_c = z - center_gz;
                bool inside_chamber = (dx_c * dx_c + dz_c * dz_c < 7 * 7);
                // Passage extends along +Z from center_gz to center_gz + 28, width ~4
                bool inside_passage = (abs(dx_c) <= 3 && dz_c >= 0 && dz_c <= 28);

                for (int y = base_gy; y <= base_gy + dome_h; ++y) {
                    if (inside_chamber && y <= base_gy + 10) continue; // Vault void
                    if (inside_passage && y <= base_gy + 7) continue;  // Passage void

                    bool is_top = (y == base_gy + dome_h);
                    Color c = is_top ? col_turf_top : col_earth_subsoil;
                    plot(x, y, z, c);
                }
            }
        }
    }

    // 2. Second Pass: Rasterize Megalithic Stones (Orthostats, Capstones, Walls, Corbels)
    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        if (n->type == MEGALITH_PRIMITIVE_MOUND) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int h_span = gy1 - gy0;

        for (int y = gy0; y < gy1; ++y) {
            float y_ratio = (h_span > 0) ? (float)(y - gy0) / (float)h_span : 0.0f;

            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    uint32_t h = mega_hash3(x, y, z);

                    // Weathered Erosion Noise: Chipped corners & natural organic profile
                    bool is_edge_x = (x == gx0 || x == gx1 - 1);
                    bool is_edge_z = (z == gz0 || z == gz1 - 1);
                    bool is_corner = is_edge_x && is_edge_z;

                    if (n->type == MEGALITH_PRIMITIVE_ORTHOSTAT || n->type == MEGALITH_PRIMITIVE_WALL_SLAB) {
                        // Tapering toward the sky (broader base, pointed apex)
                        if (y_ratio > 0.7f && is_corner) continue; // Taper top corners
                        if (y_ratio > 0.85f && (is_edge_x || is_edge_z) && (h % 3 == 0)) continue;
                        // Surface erosion
                        if (is_corner && (h % 2 == 0)) continue;
                    } else if (n->type == MEGALITH_PRIMITIVE_CAPSTONE) {
                        // Bulging heavy convex table profile: rounded top corners
                        if (y == gy1 - 1 && is_corner) continue;
                        if (y == gy1 - 1 && (is_edge_x || is_edge_z) && (h % 4 == 0)) continue;
                    }

                    // Prehistoric Mottled Color Assignment
                    Color c = col_granite_light;
                    int variant = h % 10;
                    if (variant == 0 || variant == 1) c = col_granite_dark;
                    else if (variant == 2) c = col_moss_green;
                    else if (variant == 3 && (y == gy1 - 1)) c = col_lichen_gold;

                    if (n->type == MEGALITH_PRIMITIVE_FLOOR_SLAB) {
                        c = ((x + z) % 2 == 0) ? col_floor_slab : col_granite_dark;
                    } else if (n->type == MEGALITH_PRIMITIVE_HEARTH) {
                        bool is_center = (x == mid_x && z == mid_z);
                        c = is_center ? col_ember : col_ash;
                    }

                    plot(x, y, z, c);
                }
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // MEGALITH_RASTERIZER_H
