#ifndef HYPERBOREAN_RASTERIZER_H
#define HYPERBOREAN_RASTERIZER_H

#include "hyperborean_grammar.h"
#include "raylib.h"
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer for placing a voxel at grid coordinate (gx, gy, gz)
typedef void (*HyperVoxelPlotFunc)(int gx, int gy, int gz, Color color);

// Fast 3D integer hash for procedural erosion and weathering
static inline uint32_t hyper_hash3(int x, int y, int z) {
    uint32_t h = (uint32_t)(x * 73856093 ^ y * 19349663 ^ z * 83492791);
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

// Rasterizes an abstract HyperPlan into discrete 3D voxels
static inline void rasterize_hyperborean_plan(const HyperPlan *plan,
                                              int center_gx, int center_gz, int base_gy,
                                              HyperVoxelPlotFunc plot) {
    if (!plan || !plot) return;

    int mv = 2; // 2 voxels per modular meter

    // Color Palette: The Apollonian & Chthonic Contrast
    Color col_sarsen_light   = (Color){ 145, 140, 130, 255 }; // Weathered sarsen grey
    Color col_sarsen_dark    = (Color){ 100, 95, 90, 255 };   // Cleft granite
    Color col_lichen_gold    = (Color){ 185, 165, 80, 255 };  // Ancient yellow lichen
    Color col_moss_green     = (Color){ 80, 110, 65, 255 };   // Damp green moss
    Color col_marble_white   = (Color){ 242, 240, 235, 255 }; // Gleaming Pentelic white
    Color col_marble_shadow  = (Color){ 212, 210, 204, 255 }; // Fluting shadow
    Color col_gold_accent    = (Color){ 240, 195, 55, 255 };  // Chased gold leaf
    Color col_pbf_water      = (Color){ 50, 160, 220, 255 };  // Sacred PBF fluid
    Color col_ember          = (Color){ 255, 125, 25, 255 };  // Sacred fire core
    Color col_ember_yellow   = (Color){ 255, 220, 60, 255 };  // Hot flame tip
    Color col_bronze         = (Color){ 110, 90, 60, 255 };   // Antique bronze tripod
    Color col_cypress_trunk  = (Color){ 85, 60, 45, 255 };    // Trunk wood
    Color col_cypress_leaf   = (Color){ 35, 75, 40, 255 };    // Tall cypress foliage

    for (int i = 0; i < plan->node_count; ++i) {
        const HyperNode *n = &plan->nodes[i];

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int span_y = gy1 - gy0;

        // -------------------------------------------------------------------
        // 1. Megalithic Primitives (Menhirs, Heel Stone, Rough Lintels)
        // -------------------------------------------------------------------
        if (n->type == HYPER_PRIM_MENHIR || n->type == HYPER_PRIM_HEEL_STONE ||
            n->type == HYPER_PRIM_ROUGH_LINTEL) {
            for (int y = gy0; y < gy1; ++y) {
                float y_ratio = span_y > 0 ? (float)(y - gy0) / (float)span_y : 0.0f;
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        uint32_t h = hyper_hash3(x, y, z);
                        bool is_edge_x = (x == gx0 || x == gx1 - 1);
                        bool is_edge_z = (z == gz0 || z == gz1 - 1);
                        bool is_corner = is_edge_x && is_edge_z;

                        // Weathered chipping and skyward tapering
                        if (y_ratio > 0.75f && is_corner) continue;
                        if (y_ratio > 0.85f && (is_edge_x || is_edge_z) && (h % 3 == 0)) continue;
                        if (is_corner && (h % 2 == 0)) continue;

                        Color c = col_sarsen_light;
                        int v = h % 10;
                        if (v == 0 || v == 1) c = col_sarsen_dark;
                        else if (v == 2) c = col_moss_green;
                        else if (v == 3 && (y == gy1 - 1)) c = col_lichen_gold;

                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 2. Classical Fluted Marble Columns & Architraves
        // -------------------------------------------------------------------
        else if (n->type == HYPER_PRIM_FLUTED_COLUMN) {
            for (int y = gy0; y < gy1; ++y) {
                bool is_capital = (y == gy1 - 1);
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        // Fluting texture: corners slightly shadowed
                        bool is_edge = (x == gx0 || x == gx1 - 1 || z == gz0 || z == gz1 - 1);
                        Color c = is_capital ? col_gold_accent : (is_edge ? col_marble_shadow : col_marble_white);
                        plot(x, y, z, c);
                    }
                }
            }
        }
        else if (n->type == HYPER_PRIM_ARCHITRAVE || n->type == HYPER_PRIM_CYCLOPEAN_ARCHITRAVE) {
            for (int y = gy0; y < gy1; ++y) {
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        // Bottom fillet moulding / guttae accents
                        bool is_bottom = (y == gy0);
                        Color c = is_bottom ? col_marble_shadow : col_marble_white;
                        if (n->type == HYPER_PRIM_CYCLOPEAN_ARCHITRAVE && y == gy1 - 1 && ((x + z) % 3 == 0)) {
                            c = col_gold_accent; // Gilded triglyphs
                        }
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 3. Classical Doric Pediment with Gold Tympanum Relief
        // -------------------------------------------------------------------
        else if (n->type == HYPER_PRIM_PEDIMENT) {
            int width = gx1 - gx0;
            for (int y = gy0; y < gy1; ++y) {
                int dy = y - gy0;
                int inset = dy; // Inward pitch
                int px0 = gx0 + inset;
                int px1 = gx1 - inset;
                if (px0 >= px1) px0 = px1 - 1;

                for (int z = gz0; z < gz1; ++z) {
                    for (int x = px0; x < px1; ++x) {
                        bool is_cornice = (x == px0 || x == px1 - 1 || y == gy1 - 1);
                        // Gold leaf tympanum relief in the center field
                        Color c = is_cornice ? col_marble_white : col_gold_accent;
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 4. Sacred Sunken PBF Fluid Moat & Tholos Podium
        // -------------------------------------------------------------------
        else if (n->type == HYPER_PRIM_FLUID_BASIN) {
            // Fluids disabled for now
            /*
            int rx = (gx1 - gx0) / 2;
            int rz = (gz1 - gz0) / 2;
            for (int z = gz0; z < gz1; ++z) {
                float nz = (float)(z - mid_z) / (float)rz;
                for (int x = gx0; x < gx1; ++x) {
                    float nx = (float)(x - mid_x) / (float)rx;
                    float d2 = nx * nx + nz * nz;
                    if (d2 <= 1.0f && d2 >= 0.45f) {
                        // Annular circular fluid canal
                        for (int y = gy0; y < gy1; ++y) {
                            plot(x, y, z, col_pbf_water);
                        }
                    }
                }
            }
            */
        }
        else if (n->type == HYPER_PRIM_THOLOS_PODIUM) {
            int rx = (gx1 - gx0) / 2;
            int rz = (gz1 - gz0) / 2;
            for (int z = gz0; z < gz1; ++z) {
                float nz = (float)(z - mid_z) / (float)rz;
                for (int x = gx0; x < gx1; ++x) {
                    float nx = (float)(x - mid_x) / (float)rx;
                    float d2 = nx * nx + nz * nz;
                    if (d2 <= 1.0f) {
                        for (int y = gy0; y < gy1; ++y) {
                            bool is_rim = (d2 > 0.8f);
                            plot(x, y, z, is_rim ? col_marble_shadow : col_marble_white);
                        }
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 5. Eternal Sacred Flame Brazier & Votive Statues
        // -------------------------------------------------------------------
        else if (n->type == HYPER_PRIM_BRAZIER) {
            for (int y = gy0; y < gy1; ++y) {
                bool is_fire = (y > gy0);
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        if (!is_fire) {
                            plot(x, y, z, col_bronze); // Bronze tripod pedestal
                        } else {
                            bool is_center = (x == mid_x && z == mid_z);
                            plot(x, y, z, is_center ? col_ember_yellow : col_ember);
                        }
                    }
                }
            }
        }
        else if (n->type == HYPER_PRIM_STATUE) {
            for (int y = gy0; y < gy1; ++y) {
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        Color c = (n->material == HYPER_MAT_BRONZE) ? col_bronze : col_marble_white;
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 6. Sacred Grove Italian Cypresses
        // -------------------------------------------------------------------
        else if (n->type == HYPER_PRIM_CYPRESS) {
            for (int y = gy0; y < gy1; ++y) {
                bool is_trunk = (y <= gy0 + 1);
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        if (is_trunk) {
                            plot(x, y, z, col_cypress_trunk);
                        } else {
                            plot(x, y, z, col_cypress_leaf);
                        }
                    }
                }
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // HYPERBOREAN_RASTERIZER_H
