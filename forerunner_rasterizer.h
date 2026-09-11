#ifndef FORERUNNER_RASTERIZER_H
#define FORERUNNER_RASTERIZER_H

#include "forerunner_grammar.h"
#include "raylib.h"
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer for placing a voxel at grid coordinate (gx, gy, gz)
typedef void (*ForerunnerVoxelPlotFunc)(int gx, int gy, int gz, Color color);

// Fast 3D integer hash for metallic panel seams
static inline uint32_t forerunner_hash3(int x, int y, int z) {
    uint32_t h = (uint32_t)(x * 73856093 ^ y * 19349663 ^ z * 83492791);
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

// Rasterizes a ForerunnerPlan into 3D discrete voxels
static inline void rasterize_forerunner_plan(const ForerunnerPlan *plan,
                                             int center_gx, int center_gz, int base_gy,
                                             ForerunnerVoxelPlotFunc plot) {
    if (!plan || !plot) return;

    int mv = 2; // 2 voxels per modular meter

    // Forerunner Color Palette
    Color col_pewter_hull      = (Color){ 110, 120, 132, 255 }; // Cold titanium alloy
    Color col_pewter_dark      = (Color){ 78, 85, 96, 255 };    // Panel seam shadow
    Color col_pewter_bright    = (Color){ 145, 155, 168, 255 }; // Chamfer highlight
    Color col_basalt           = (Color){ 38, 42, 48, 255 };    // Abyssal bedrock
    Color col_hardlight_cyan   = (Color){ 45, 225, 255, 255 };  // Emissive cyan hard-light
    Color col_hardlight_core   = (Color){ 160, 245, 255, 255 }; // Pure energy core
    Color col_solar_amber      = (Color){ 255, 165, 35, 255 };  // Telemetry beam emitter
    Color col_bronze_trim      = (Color){ 145, 115, 75, 255 };  // Warm metallic accent

    // Compute minimum non-void modular Y level to anchor foundation at base_gy
    int min_non_void_y = 0;
    bool found_non_void = false;
    for (int i = 0; i < plan->node_count; ++i) {
        if (plan->nodes[i].is_void) continue;
        if (!found_non_void || plan->nodes[i].box.min_y < min_non_void_y) {
            min_non_void_y = plan->nodes[i].box.min_y;
            found_non_void = true;
        }
    }
    int offset_y = base_gy - min_non_void_y * mv;

    for (int i = 0; i < plan->node_count; ++i) {
        const ForerunnerNode *n = &plan->nodes[i];
        if (n->is_void) continue; // Voids define empty space in the chasm

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = offset_y  + n->box.min_y * mv;
        int gy1 = offset_y  + n->box.max_y * mv;
        if (gy0 < base_gy) gy0 = base_gy;
        if (gy1 < gy0) continue;

        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int span_y = gy1 - gy0;

        // -------------------------------------------------------------------
        // 1. Canted Angular Pylons (15-25 degree rake)
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------
        // 1. Canted Angular Pylons (15-25 degree rake)
        // -------------------------------------------------------------------
        if (n->type == FORERUNNER_PRIM_PYLON) {
            // Extrude pylon base downward to base_gy bedrock
            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    for (int y = gy0 - 1; y >= base_gy; --y) {
                        plot(x, y, z, col_pewter_dark);
                    }
                }
            }

            float cant_x = n->cant_angle_deg * (3.14159265f / 180.0f);
            float cant_z = n->cant_angle_z * (3.14159265f / 180.0f);
            float slope_x = tanf(fabsf(cant_x));
            float slope_z = tanf(fabsf(cant_z));

            for (int y = gy0; y < gy1; ++y) {
                int dy = y - gy0;
                int shift_x = (int)(dy * slope_x);
                int shift_z = (int)(dy * slope_z);

                int px0 = gx0;
                int px1 = gx1;
                if (n->cant_angle_deg > 0.0f) {
                    px0 = gx0 + shift_x;
                    px1 = gx1 + shift_x;
                } else if (n->cant_angle_deg < 0.0f) {
                    px0 = gx0 - shift_x;
                    px1 = gx1 - shift_x;
                }

                int pz0 = gz0;
                int pz1 = gz1;
                if (n->cant_angle_z > 0.0f) {
                    pz0 = gz0 + shift_z;
                    pz1 = gz1 + shift_z;
                } else if (n->cant_angle_z < 0.0f) {
                    pz0 = gz0 - shift_z;
                    pz1 = gz1 - shift_z;
                }

                // Seal stair-step overhang into layer y - 1 to maintain vertical 6-connectivity
                if (dy > 0) {
                    for (int z = pz0; z < pz1; ++z) {
                        for (int x = px0; x < px1; ++x) {
                            plot(x, y - 1, z, col_pewter_dark);
                        }
                    }
                }

                for (int z = pz0; z < pz1; ++z) {
                    for (int x = px0; x < px1; ++x) {
                        bool is_chamfer = (x == px0 || x == px1 - 1 || z == pz0 || z == pz1 - 1);
                        bool is_seam = ((y % 4 == 0) || ((x + z) % 6 == 0));

                        Color c = col_pewter_hull;
                        if (is_chamfer) c = col_pewter_bright;
                        else if (is_seam) c = col_pewter_dark;

                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 2. Heavy Chevron Lintel Beams & Terraces
        // -------------------------------------------------------------------
        else if (n->type == FORERUNNER_PRIM_LINTEL || n->type == FORERUNNER_PRIM_TERRACE ||
                 n->type == FORERUNNER_PRIM_BRIDGE_SPAN) {
            Color base_col = (n->material == FORERUNNER_MAT_DARK_BASALT) ? col_basalt :
                             (n->material == FORERUNNER_MAT_BRONZE_ACCENT) ? col_bronze_trim : col_pewter_hull;

            for (int y = gy0; y < gy1; ++y) {
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        bool is_edge = (x == gx0 || x == gx1 - 1 || z == gz0 || z == gz1 - 1);
                        Color c = is_edge ? col_pewter_dark : base_col;
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 3. Emissive Hard-Light Energy Bridge
        // -------------------------------------------------------------------
        else if (n->type == FORERUNNER_PRIM_HARDLIGHT_BRIDGE) {
            for (int y = gy0; y < gy1; ++y) {
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        bool is_center = (x == mid_x);
                        Color c = is_center ? col_hardlight_core : col_hardlight_cyan;
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 4. Vertical Telemetry Conduit Spires
        // -------------------------------------------------------------------
        else if (n->type == FORERUNNER_PRIM_SPIRE) {
            // Extrude spire base downward to base_gy bedrock
            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    for (int y = gy0 - 1; y >= base_gy; --y) {
                        plot(x, y, z, col_pewter_dark);
                    }
                }
            }

            for (int y = gy0; y < gy1; ++y) {
                float y_ratio = span_y > 0 ? (float)(y - gy0) / (float)span_y : 0.0f;
                int inset = (int)(y_ratio * 1.5f); // Taper upward
                int sx0 = gx0 + inset;
                int sx1 = gx1 - inset;
                int sz0 = gz0 + inset;
                int sz1 = gz1 - inset;
                if (sx0 >= sx1 || sz0 >= sz1) {
                    sx0 = mid_x; sx1 = mid_x + 1;
                    sz0 = mid_z; sz1 = mid_z + 1;
                }

                for (int z = sz0; z < sz1; ++z) {
                    for (int x = sx0; x < sx1; ++x) {
                        bool is_corner = (x == sx0 || x == sx1 - 1) && (z == sz0 || z == sz1 - 1);
                        Color c = is_corner ? col_pewter_bright : col_pewter_hull;
                        plot(x, y, z, c);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 5. Hovering Central Gravity Matrix Core
        // -------------------------------------------------------------------
        else if (n->type == FORERUNNER_PRIM_GRAVITY_CORE) {
            int rx = (gx1 - gx0) / 2;
            int rz = (gz1 - gz0) / 2;
            int ry = (gy1 - gy0) / 2;
            int mid_y = (gy0 + gy1) / 2;

            for (int y = gy0; y < gy1; ++y) {
                float ny = ry > 0 ? (float)(y - mid_y) / (float)ry : 0.0f;
                for (int z = gz0; z < gz1; ++z) {
                    float nz = rz > 0 ? (float)(z - mid_z) / (float)rz : 0.0f;
                    for (int x = gx0; x < gx1; ++x) {
                        float nx = rx > 0 ? (float)(x - mid_x) / (float)rx : 0.0f;
                        float d = fabsf(nx) + fabsf(ny) + fabsf(nz); // Octahedral energy diamond
                        if (d <= 1.2f) {
                            bool is_core = (d < 0.6f);
                            Color c = is_core ? col_hardlight_core : col_hardlight_cyan;
                            plot(x, y, z, c);
                        }
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        // 6. Recessed Cyan & Solar Amber Light Channels
        // -------------------------------------------------------------------
        else if (n->type == FORERUNNER_PRIM_LIGHT_CHANNEL) {
            Color c = (n->material == FORERUNNER_MAT_SOLAR_AMBER) ? col_solar_amber : col_hardlight_cyan;
            for (int y = gy0; y < gy1; ++y) {
                for (int z = gz0; z < gz1; ++z) {
                    for (int x = gx0; x < gx1; ++x) {
                        plot(x, y, z, c);
                    }
                }
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // FORERUNNER_RASTERIZER_H
