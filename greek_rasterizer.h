#ifndef GREEK_RASTERIZER_H
#define GREEK_RASTERIZER_H

#include "greek_grammar.h"
#include "raylib.h"
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer or callback type for placing a voxel at grid coordinate (gx, gy, gz)
typedef void (*VoxelPlotFunc)(int gx, int gy, int gz, Color color);

// Rasterizes an abstract TemplePlan into discrete grid voxels
static inline void rasterize_temple_plan(const TemplePlan *plan,
                                         int center_gx, int center_gz, int base_gy,
                                         VoxelPlotFunc plot) {
    if (!plan || !plot) return;

    int mv = plan->module_voxels;
    if (mv < 1) mv = 2;

    // Palette: Classic Ancient Greek Polychromy
    Color col_step_light   = (Color){ 225, 220, 210, 255 }; // Stylobate top tier
    Color col_step_dark    = (Color){ 180, 175, 165, 255 }; // Stereobate lower tiers
    Color col_wall         = (Color){ 235, 228, 215, 255 }; // Cella marble wall
    Color col_floor        = (Color){ 160, 75, 60, 255 };   // Terracotta mosaic floor
    Color col_column       = (Color){ 245, 242, 235, 255 }; // White fluted marble
    Color col_drum_seam    = (Color){ 205, 200, 190, 255 }; // Drum joint line
    Color col_lintel       = (Color){ 220, 215, 205, 255 }; // Architrave stone
    Color col_triglyph     = (Color){ 40, 90, 165, 255 };   // Aegean blue frieze
    Color col_roof_tile    = (Color){ 195, 85, 55, 255 };   // Terracotta clay tiles
    Color col_altar        = (Color){ 235, 185, 50, 255 };   // Gold / Bronze sacred altar

    // 1. Rasterize Stylobate / Stereobate Steps
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_STEP) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        Color c = (n->box.min_y >= -1) ? col_step_light : col_step_dark;
        for (int y = gy0; y < gy1; ++y) {
            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    plot(x, y, z, c);
                }
            }
        }
    }

    // 2. Rasterize Cella Floor & Interior
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_CELL) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        // Paved floor at base level
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, base_gy, z, col_floor);
            }
        }

        // Central Altar in the middle of the Cella
        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        for (int dy = 1; dy <= 2; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                for (int dx = -1; dx <= 1; ++dx) {
                    plot(mid_x + dx, base_gy + dy, mid_z + dz, col_altar);
                }
            }
        }
    }

    // 3. Rasterize Cella Masonry Walls
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_WALL) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        for (int y = gy0; y < gy1; ++y) {
            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    plot(x, y, z, col_wall);
                }
            }
        }
    }

    // 4. Rasterize Columns (Segmented Cylindrical Drums)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_COLUMN) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        for (int y = gy0; y < gy1; ++y) {
            // Drum seam every 3 voxels
            bool is_seam = ((y - gy0) % 3 == 0) || (y == gy1 - 1);
            Color drum_color = is_seam ? col_drum_seam : col_column;

            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    // Slight corner softening for round column feel if width > 2
                    if ((gx1 - gx0 >= 3) && (gz1 - gz0 >= 3)) {
                        if ((x == gx0 || x == gx1 - 1) && (z == gz0 || z == gz1 - 1)) {
                            continue; // cut 4 corners
                        }
                    }
                    plot(x, y, z, drum_color);
                }
            }
        }

        // Doric Capital at the top (flare out by 1 voxel if possible)
        int cap_y = gy1;
        for (int z = gz0 - 1; z <= gz1; ++z) {
            for (int x = gx0 - 1; x <= gx1; ++x) {
                plot(x, cap_y, z, col_drum_seam);
            }
        }
    }

    // 5. Rasterize Entablature & Architrave Lintels
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_LINTEL) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        for (int y = gy0; y < gy1; ++y) {
            for (int z = gz0; z < gz1; ++z) {
                for (int x = gx0; x < gx1; ++x) {
                    // Alternate blue triglyphs every 3 voxels along facade
                    bool is_triglyph = ((x + z) % 4 == 0) && (y == gy0 + 1);
                    plot(x, y, z, is_triglyph ? col_triglyph : col_lintel);
                }
            }
        }
    }

    // 6. Rasterize Triangular Pediment & Pitched Roof
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_PEDIMENT) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy   + n->box.min_y * mv;
        int gy1 = base_gy   + n->box.max_y * mv;

        int half_w = (gx1 - gx0) / 2;
        int mid_x = (gx0 + gx1) / 2;
        int ped_h = gy1 - gy0;
        if (ped_h < 2) ped_h = 2;

        for (int y = gy0; y <= gy0 + ped_h; ++y) {
            int cur_dy = y - gy0;
            // Width narrows as y rises (sloped gable roof)
            int inset = (cur_dy * half_w) / ped_h;
            int cur_x0 = gx0 + inset;
            int cur_x1 = gx1 - inset;

            for (int z = gz0; z < gz1; ++z) {
                for (int x = cur_x0; x < cur_x1; ++x) {
                    // Front and back gable tympanum wall
                    bool is_tympanum = (z == gz0 || z == gz1 - 1);
                    // Outer roof shell is tile, inner is wood/open
                    bool is_outer_roof = (x == cur_x0 || x == cur_x1 - 1 || y == gy0 + ped_h);

                    if (is_tympanum) {
                        plot(x, y, z, col_wall);
                    } else if (is_outer_roof) {
                        plot(x, y, z, col_roof_tile);
                    }
                }
            }
        }
    }

    Color col_court_paving = (Color){ 205, 200, 190, 255 }; // Courtyard stone tiles
    Color col_pool_rim     = (Color){ 230, 225, 215, 255 }; // Marble basin coping
    Color col_water        = (Color){ 50, 160, 220, 220 };   // Sparkling azure liquid
    Color col_trunk        = (Color){ 105, 75, 45, 255 };    // Olive timber trunk
    Color col_foliage      = (Color){ 75, 115, 60, 255 };    // Mediterranean olive leaves

    // 7. Rasterize Courtyard Pavement
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_COURTYARD) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, base_gy, z, col_court_paving);
            }
        }
    }

    // 8. Rasterize Sunken Reflection Pool (Basin & Fluid Cells)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_POOL) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                bool is_rim = (x == gx0 || x == gx1 - 1 || z == gz0 || z == gz1 - 1);
                if (is_rim) {
                    plot(x, base_gy, z, col_pool_rim);
                } else {
                    // Sunken bed
                    plot(x, base_gy - 1, z, col_floor);
                    // Liquid water surface
                    plot(x, base_gy, z, col_water);
                }
            }
        }
    }

    // 9. Rasterize Sacred Olive Trees (Trunk + Clustered Foliage)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_TREE) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy;
        int gy1 = base_gy + (n->box.max_y - n->box.min_y) * mv;

        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int trunk_h = (gy1 - gy0) * 2 / 3;

        // Trunk
        for (int y = gy0; y < gy0 + trunk_h; ++y) {
            plot(mid_x, y, mid_z, col_trunk);
            plot(mid_x + 1, y, mid_z, col_trunk);
            plot(mid_x, y, mid_z + 1, col_trunk);
            plot(mid_x + 1, y, mid_z + 1, col_trunk);
        }

        // Stepped rounded canopy
        int can_y0 = gy0 + trunk_h - 1;
        int can_y1 = gy1 + 1;
        int tree_h_mod = n->box.max_y - n->box.min_y;
        bool is_cypress = (tree_h_mod >= 5); // Tall slender Italian cypress

        Color col_foliage_tree = is_cypress ? (Color){ 35, 75, 40, 255 } : col_foliage;

        if (is_cypress) {
            // Slender vertical spire canopy
            for (int y = gy0 + 2; y <= gy1 + 1; ++y) {
                int r_layer = (y > gy1 - 1) ? 1 : 2;
                for (int dz = -r_layer; dz <= r_layer; ++dz) {
                    for (int dx = -r_layer; dx <= r_layer; ++dx) {
                        if (dx * dx + dz * dz <= r_layer * r_layer + 1) {
                            plot(mid_x + dx, y, mid_z + dz, col_foliage_tree);
                        }
                    }
                }
            }
        } else {
            // Sprawling broad olive canopy
            int can_r = (gx1 - gx0) / 2 + 1;
            for (int y = can_y0; y <= can_y1; ++y) {
                int dy = y - (can_y0 + can_y1) / 2;
                int r_layer = can_r - abs(dy) / 2;
                if (r_layer < 1) r_layer = 1;

                for (int dz = -r_layer; dz <= r_layer; ++dz) {
                    for (int dx = -r_layer; dx <= r_layer; ++dx) {
                        if (dx * dx + dz * dz <= r_layer * r_layer + 1) {
                            plot(mid_x + dx, y, mid_z + dz, col_foliage_tree);
                        }
                    }
                }
            }
        }
    }

    // 10. Rasterize Processional Garden Pathways
    Color col_path_a = (Color){ 215, 210, 195, 255 }; // Light flagstone
    Color col_path_b = (Color){ 190, 185, 175, 255 }; // Stepping stone accent
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_PATH) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                Color c = ((x + z) % 3 == 0) ? col_path_b : col_path_a;
                plot(x, base_gy, z, c);
            }
        }
    }

    // 11. Rasterize The Tholos (Circular Monopteros Rotunda)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_THOLOS) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;
        int r_tholos = (gx1 - gx0) / 2;

        // Podium steps (circular)
        for (int dz = -r_tholos; dz <= r_tholos; ++dz) {
            for (int dx = -r_tholos; dx <= r_tholos; ++dx) {
                if (dx * dx + dz * dz <= r_tholos * r_tholos) {
                    plot(mid_x + dx, base_gy, mid_z + dz, col_step_light);
                }
            }
        }

        // 6 Columns in a circle
        int col_h = 6;
        for (int col_i = 0; col_i < 6; ++col_i) {
            float angle = col_i * (3.14159265f / 3.0f);
            int cx = mid_x + (int)((r_tholos - 1) * cosf(angle));
            int cz = mid_z + (int)((r_tholos - 1) * sinf(angle));
            for (int y = base_gy + 1; y <= base_gy + col_h; ++y) {
                plot(cx, y, cz, col_column);
            }
        }

        // Circular Architrave ring
        int arch_y = base_gy + col_h + 1;
        for (int dz = -r_tholos; dz <= r_tholos; ++dz) {
            for (int dx = -r_tholos; dx <= r_tholos; ++dx) {
                int d2 = dx * dx + dz * dz;
                if (d2 >= (r_tholos - 2) * (r_tholos - 2) && d2 <= r_tholos * r_tholos) {
                    plot(mid_x + dx, arch_y, mid_z + dz, col_lintel);
                }
            }
        }

        // Conical Roof
        for (int layer = 0; layer <= 3; ++layer) {
            int y = arch_y + 1 + layer;
            int r_cone = r_tholos - layer;
            if (r_cone < 1) r_cone = 1;
            for (int dz = -r_cone; dz <= r_cone; ++dz) {
                for (int dx = -r_cone; dx <= r_cone; ++dx) {
                    if (dx * dx + dz * dz <= r_cone * r_cone) {
                        plot(mid_x + dx, y, mid_z + dz, col_roof_tile);
                    }
                }
            }
        }
        // Gold apex finial
        plot(mid_x, arch_y + 5, mid_z, col_altar);

        // Center sacred bronze tripod
        for (int y = base_gy + 1; y <= base_gy + 2; ++y) {
            plot(mid_x, y, mid_z, col_altar);
        }
    }

    // 12. Rasterize Exedrae (Philosopher Semicircular Marble Benches)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_EXEDRA) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        // Base plinth
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, base_gy, z, col_step_light);
                plot(x, base_gy + 1, z, col_column); // seat
            }
        }
        // High curved backrest along outer edge
        bool on_left = (gx0 < center_gx);
        int back_x = on_left ? gx0 : gx1 - 1;
        for (int z = gz0; z < gz1; ++z) {
            plot(back_x, base_gy + 2, z, col_wall);
            plot(back_x, base_gy + 3, z, col_wall);
        }
        // Armrest piers at the ends
        plot(on_left ? gx1 - 1 : gx0, base_gy + 2, gz0, col_step_dark);
        plot(on_left ? gx1 - 1 : gx0, base_gy + 2, gz1 - 1, col_step_dark);
    }

    // 13. Rasterize Naiskoi (Miniature Temple Votive Shrines)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_NAISKOS) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy;
        int gy1 = base_gy + (n->box.max_y - n->box.min_y) * mv;

        // Stepped base
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, gy0, z, col_step_light);
            }
        }
        // Solid rear wall
        bool on_left = (gx0 < center_gx);
        int rear_x = on_left ? gx0 : gx1 - 1;
        for (int y = gy0 + 1; y < gy1; ++y) {
            for (int z = gz0; z < gz1; ++z) {
                plot(rear_x, y, z, col_wall);
            }
        }
        // 2 Front mini columns
        int front_x = on_left ? gx1 - 1 : gx0;
        for (int y = gy0 + 1; y < gy1; ++y) {
            plot(front_x, y, gz0, col_column);
            plot(front_x, y, gz1 - 1, col_column);
        }
        // Entablature and mini pediment
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, gy1, z, col_roof_tile);
                plot(x, gy1 + 1, z, col_roof_tile);
            }
        }
        // Votive statuette in niche
        int mid_z = (gz0 + gz1) / 2;
        plot(on_left ? gx0 + 1 : gx1 - 2, gy0 + 1, mid_z, col_altar);
        plot(on_left ? gx0 + 1 : gx1 - 2, gy0 + 2, mid_z, col_altar);
    }

    // 14. Rasterize Garden Fountains (Marble Basin & PBF Fluid)
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_FOUNTAIN) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;

        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                bool is_rim = (x == gx0 || x == gx1 - 1 || z == gz0 || z == gz1 - 1);
                if (is_rim) {
                    plot(x, base_gy, z, col_pool_rim);
                    plot(x, base_gy + 1, z, col_pool_rim);
                } else {
                    plot(x, base_gy - 1, z, col_court_paving);
                    plot(x, base_gy, z, col_water); // live fluid
                }
            }
        }
        // Center spout column
        plot(mid_x, base_gy + 1, mid_z, col_column);
        plot(mid_x, base_gy + 2, mid_z, col_column);
    }

    // 15. Rasterize Pergolas (Timber Post-and-Beam Arbor with Vines)
    Color col_vine = (Color){ 55, 130, 45, 255 };
    Color col_grape = (Color){ 125, 45, 110, 255 };
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_PERGOLA) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int gy0 = base_gy;
        int gy1 = base_gy + (n->box.max_y - n->box.min_y) * mv;

        // Timber uprights at 4 corners and intervals
        for (int z = gz0; z < gz1; z += 4) {
            for (int y = gy0; y < gy1; ++y) {
                plot(gx0, y, z, col_trunk);
                plot(gx1 - 1, y, z, col_trunk);
            }
        }
        // Top crossbeams
        for (int z = gz0; z < gz1; ++z) {
            plot(gx0, gy1, z, col_trunk);
            plot(gx1 - 1, gy1, z, col_trunk);
            if (z % 2 == 0) {
                for (int x = gx0; x < gx1; ++x) {
                    plot(x, gy1, z, col_trunk);
                }
            }
        }
        // Hanging green vines and grapes
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                if ((x + z) % 3 != 0) {
                    plot(x, gy1 + 1, z, col_vine);
                } else if ((x * z) % 5 == 0) {
                    plot(x, gy1 - 1, z, col_grape); // hanging grape cluster
                }
            }
        }
    }

    // 16. Rasterize Ceremonial Braziers (Stone Fire Altars)
    Color col_fire1 = (Color){ 255, 135, 20, 255 }; // Bright flame
    Color col_fire2 = (Color){ 240, 55, 15, 255 };  // Deep ember
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_BRAZIER) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;
        int mid_x = (gx0 + gx1) / 2;
        int mid_z = (gz0 + gz1) / 2;

        // Base pedestal
        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                plot(x, base_gy, z, col_step_dark);
                plot(x, base_gy + 1, z, col_wall);
            }
        }
        // Bronze brazier bowl
        plot(mid_x, base_gy + 2, mid_z, col_altar);
        plot(mid_x + 1, base_gy + 2, mid_z, col_altar);
        plot(mid_x, base_gy + 2, mid_z + 1, col_altar);
        plot(mid_x + 1, base_gy + 2, mid_z + 1, col_altar);
        // Flaming embers
        plot(mid_x, base_gy + 3, mid_z, col_fire1);
        plot(mid_x + 1, base_gy + 3, mid_z, col_fire2);
        plot(mid_x, base_gy + 3, mid_z + 1, col_fire2);
        plot(mid_x + 1, base_gy + 3, mid_z + 1, col_fire1);
    }

    // 17. Rasterize Flowering Shrubs & Boxwood Borders
    Color col_boxwood = (Color){ 45, 95, 40, 255 };
    Color col_flower_lavender = (Color){ 145, 115, 185, 255 };
    Color col_flower_poppy    = (Color){ 220, 50, 40, 255 };
    Color col_flower_gold     = (Color){ 245, 195, 45, 255 };
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_SHRUB) continue;

        int gx0 = center_gx + n->box.min_x * mv;
        int gx1 = center_gx + n->box.max_x * mv;
        int gz0 = center_gz + n->box.min_z * mv;
        int gz1 = center_gz + n->box.max_z * mv;

        for (int z = gz0; z < gz1; ++z) {
            for (int x = gx0; x < gx1; ++x) {
                // Low boxwood hedge
                plot(x, base_gy, z, col_boxwood);

                // Varied blooms on top
                Color bloom = ((x + z) % 3 == 0) ? col_flower_lavender :
                              ((x + z) % 3 == 1) ? col_flower_poppy : col_flower_gold;
                plot(x, base_gy + 1, z, bloom);
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // GREEK_RASTERIZER_H
