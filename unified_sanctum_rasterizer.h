#ifndef UNIFIED_SANCTUM_RASTERIZER_H
#define UNIFIED_SANCTUM_RASTERIZER_H

#include "unified_sanctum_grammar.h"
#include <math.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*SanctumPlotFn)(int gx, int gy, int gz, Color c);
typedef void (*SanctumVoidFn)(int gx, int gy, int gz);

static inline Color sanctum_color_variance(Color base, int h) {
    int v = ((h % 23) - 11);
    int r = (int)base.r + v;
    int g = (int)base.g + v;
    int b = (int)base.b + v;
    if (r < 0) r = 0; if (r > 255) r = 255;
    if (g < 0) g = 0; if (g > 255) g = 255;
    if (b < 0) b = 0; if (b > 255) b = 255;
    return (Color){ (unsigned char)r, (unsigned char)g, (unsigned char)b, base.a };
}

static inline void rasterize_sanctum_plan(const SanctumCitadelPlan *plan,
                                          int center_gx, int center_gz, int base_gy,
                                          SanctumPlotFn plot_fn,
                                          SanctumVoidFn void_fn) {
    if (!plan || !plot_fn) return;

    // Compute minimum non-void modular Y level to anchor foundation at base_gy
    int min_non_void_y = 0;
    bool found_non_void = false;
    for (int i = 0; i < plan->node_count; ++i) {
        if (plan->nodes[i].is_void) continue;
        if (!found_non_void || plan->nodes[i].y < min_non_void_y) {
            min_non_void_y = plan->nodes[i].y;
            found_non_void = true;
        }
    }
    int offset_y = base_gy - min_non_void_y;

    // PASS 1: Negative Mass Void Excavation (Chasms, Pits, Crypt Voids)
    for (int i = 0; i < plan->node_count; ++i) {
        const SanctumNode *node = &plan->nodes[i];
        if (!node->is_void) continue;

        int gx0 = center_gx + node->x - node->w / 2;
        int gx1 = center_gx + node->x + (node->w - 1) / 2;
        int gy0 = offset_y + node->y;
        int gy1 = offset_y + node->y + node->h - 1;
        int gz0 = center_gz + node->z - node->d / 2;
        int gz1 = center_gz + node->z + (node->d - 1) / 2;

        for (int y = gy0; y <= gy1; ++y) {
            for (int z = gz0; z <= gz1; ++z) {
                for (int x = gx0; x <= gx1; ++x) {
                    if (void_fn) {
                        void_fn(x, y, z);
                    }
                }
            }
        }
    }

    // PASS 2: Positive Mass Deposition (Megaliths, Marble, Titanium, Hard-Light, Fluids)
    for (int i = 0; i < plan->node_count; ++i) {
        const SanctumNode *node = &plan->nodes[i];
        if (node->is_void) continue;

        int gx0 = center_gx + node->x - node->w / 2;
        int gx1 = center_gx + node->x + (node->w - 1) / 2;
        int gy0 = offset_y + node->y;
        int gy1 = offset_y + node->y + node->h - 1;
        int gz0 = center_gz + node->z - node->d / 2;
        int gz1 = center_gz + node->z + (node->d - 1) / 2;

        float slope_x = tanf(node->cant_angle_x * DEG2RAD);
        float slope_z = tanf(node->cant_angle_z * DEG2RAD);
        float half_w = node->w * 0.5f;
        float half_d = node->d * 0.5f;

        for (int y = gy0; y <= gy1; ++y) {
            int dy = y - gy0; // Height from base of node
            float allowed_dx = half_w;
            float allowed_dz = half_d;

            if (node->cant_angle_x != 0.0f) {
                allowed_dx = half_w - dy * slope_x;
                if (allowed_dx < 0.5f) allowed_dx = 0.5f;
            }
            if (node->cant_angle_z != 0.0f) {
                allowed_dz = half_d - dy * slope_z;
                if (allowed_dz < 0.5f) allowed_dz = 0.5f;
            }

            for (int z = gz0; z <= gz1; ++z) {
                float dz = fabsf((float)(z - (center_gz + node->z)));
                if (node->cant_angle_z != 0.0f && dz > allowed_dz) continue;

                for (int x = gx0; x <= gx1; ++x) {
                    float dx = fabsf((float)(x - (center_gx + node->x)));
                    if (node->cant_angle_x != 0.0f && dx > allowed_dx) continue;

                    // Compute procedural variation hash for weathering & fluting
                    int h = (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
                    h = (h ^ (h >> 13)) * 0x85ebca6b;
                    h = abs(h);

                    Color c = node->color;

                    if (node->prim_type == SANCTUM_PRIM_STONE_ORTHOSTAT ||
                        node->prim_type == SANCTUM_PRIM_STONE_LINTEL) {
                        // Rough chthonic weathering
                        c = sanctum_color_variance(c, h);
                    } else if (node->prim_type == SANCTUM_PRIM_MARBLE_COLUMN) {
                        // Doric vertical column fluting: subtle lighting modulation along circumference
                        int angle_mod = abs((x - (center_gx + node->x)) * 3 + (z - (center_gz + node->z)) * 5) % 4;
                        if (angle_mod == 0) c = ColorBrightness(c, -0.06f);
                        else if (angle_mod == 2) c = ColorBrightness(c, 0.04f);
                    } else if (node->prim_type == SANCTUM_PRIM_HARDLIGHT_BRIDGE) {
                        // Pulsing luminescent center track
                        if (fabsf(dx) < 0.75f || fabsf(dz) < 0.75f) {
                            c = (Color){ 210, 252, 255, 255 }; // Bright cyan core line
                        }
                    } else if (node->prim_type == SANCTUM_PRIM_LIGHT_CHANNEL) {
                        c = (Color){ 45, 235, 255, 255 };
                    } else if (node->prim_type == SANCTUM_PRIM_PBF_WATER) {
                        c = (Color){ 50, 160, 220, 255 };
                    } else if (node->prim_type == SANCTUM_PRIM_COFFERED_CEILING) {
                        int cx = (int)(dx * 2.0f);
                        int cz = (int)(dz * 2.0f);
                        if ((cx % 4 == 0) || (cz % 4 == 0)) {
                            c = ColorBrightness(c, -0.15f);
                        }
                    } else if (node->prim_type == SANCTUM_PRIM_BRAZIER) {
                        c = (Color){ 255, 175, 40, 255 };
                    } else if (node->prim_type == SANCTUM_PRIM_OBELISK) {
                        // Tapered needle with gold pyramidion capstone at top
                        if (dy >= node->h - 2) {
                            c = (Color){ 255, 215, 60, 255 }; // Electrum / Gold pyramidion
                        } else {
                            c = sanctum_color_variance(c, h);
                        }
                    } else if (node->prim_type == SANCTUM_PRIM_ALTAR) {
                        // Stepped sacrificial plinth with glowing ember core
                        if (dy == node->h - 1 && fabsf(dx) < 1.0f && fabsf(dz) < 1.0f) {
                            c = (Color){ 255, 150, 30, 255 }; // Glowing sacrificial ember
                        } else {
                            c = ColorBrightness(c, ((y % 2 == 0) ? -0.05f : 0.05f));
                        }
                    }

                    plot_fn(x, y, z, c);
                }
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // UNIFIED_SANCTUM_RASTERIZER_H
