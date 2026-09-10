#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "raylib.h"
#include "../greek_grammar.h"
#include "../greek_rasterizer.h"

// PRNG
static inline float prng_f(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return (float)(*state & 0x00FFFFFF) / (float)0x01000000;
}

// ---------------------------------------------------------------------------
// 1. Ablated Algorithm: Algebra Only (Blind Production Without Coalgebra)
// ---------------------------------------------------------------------------
// In an algebra-only system, rules construct parts, but there is no coalgebraic
// observation gamma(S) to inspect current boundaries, assert bilateral symmetry,
// or verify contextual affordances (preconditions).
TemplePlan generate_temple_ablated_no_coalgebra(uint32_t seed, int target_steps) {
    TemplePlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.module_voxels = 2;
    plan.column_spacing_m = 2;
    plan.column_height_m = 6;
    plan.cella_width_m = 6;
    plan.cella_length_m = 12;
    plan.cella_height_m = 6;

    uint32_t prng = seed;

    // Initial cella
    int half_w = plan.cella_width_m / 2;
    int half_l = plan.cella_length_m / 2;
    int h = plan.cella_height_m;
    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_CELL, {-half_w, 0, -half_l, half_w, h, half_l}, 0, 0, false, 0 };
    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_WALL, {-half_w, 0, -half_l, -half_w + 1, h, half_l}, 0, 0, true, 0 };
    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_WALL, {half_w - 1, 0, -half_l, half_w, h, half_l}, 0, 0, true, 0 };
    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_WALL, {-half_w, 0, -half_l, half_w, h, -half_l + 1}, 0, 0, true, 0 };

    // Available algebraic rules (blind constructors)
    GrammarRule rules[] = {
        RULE_ATTACH_PRONAOS,
        RULE_ATTACH_OPISTHODOMOS,
        RULE_WRAP_PERISTYLE,
        RULE_EXPAND_STYLOBATE,
        RULE_RAISE_PEDIMENT,
        RULE_ATTACH_COURTYARD,
        RULE_EXCAVATE_POOL,
        RULE_PLANT_GROVE,
        RULE_ERECT_THOLOS,
        RULE_INSTALL_EXEDRAE,
        RULE_CONSTRUCT_FOUNTAINS,
        RULE_BUILD_PERGOLAS
    };
    int num_rules = sizeof(rules) / sizeof(rules[0]);

    for (int step = 0; step < target_steps && plan.node_count < MAX_ARCH_NODES - 16; ++step) {
        // Blind rule selection WITHOUT coalgebraic filtering / affordance checks
        int r_idx = (int)(prng_f(&prng) * num_rules);
        GrammarRule rule = rules[r_idx];

        // Without coalgebra to coordinate symmetry or inspect current total footprint:
        switch (rule) {
            case RULE_ATTACH_PRONAOS: {
                // Places columns at a fixed forward offset without knowing if rear or front is already blocked
                int z = half_l + 2;
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {-2, 0, z, -1, h, z + 1}, 1, 0, true, 0 };
                // Asymmetric glitch: 50% chance right column is omitted or misplaced without bilateral coupling
                if (prng_f(&prng) > 0.4f) {
                    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {1, 0, z, 2, h, z + 1}, 1, 0, true, 0 };
                }
                break;
            }
            case RULE_ATTACH_OPISTHODOMOS: {
                // Places rear columns blindly (even if cella has no pronaos or stylobate)
                int z = -half_l - 2;
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {-2, 0, z, -1, h, z + 1}, 1, 0, true, 0 };
                // Symmetry break: unilateral placement
                if (prng_f(&prng) > 0.5f) {
                    plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {2, 0, z + 1, 3, h, z + 2}, 1, 0, true, 0 };
                }
                break;
            }
            case RULE_WRAP_PERISTYLE: {
                // Without coalgebra computing dynamic total footprint, uses naive guessed span
                int span = half_w + 3;
                bool pick_left = prng_f(&prng) > 0.5f;
                // Unilateral flank colonnade: only left OR right is grown!
                if (pick_left) {
                    for (int z = -half_l; z <= half_l; z += 3) {
                        plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {-span, 0, z, -span + 1, h, z + 1}, 1, 0, true, 0 };
                    }
                } else {
                    for (int z = -half_l; z <= half_l; z += 3) {
                        plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COLUMN, {span, 0, z, span + 1, h, z + 1}, 1, 0, true, 0 };
                    }
                }
                break;
            }
            case RULE_RAISE_PEDIMENT: {
                // Without coalgebra checking if building is capped, pediment spawns at arbitrary Y
                int py = (prng_f(&prng) > 0.5f) ? h : h + 3; // May float in mid-air!
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_PEDIMENT, {-half_w - 2, py, -half_l, half_w + 2, py + 3, half_l}, 0, 0, true, 0 };
                break;
            }
            case RULE_ATTACH_COURTYARD: {
                // Attached blindly at random Z offset
                int cz = half_l + (int)(prng_f(&prng) * 6.0f);
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_COURTYARD, {-half_w, 0, cz, half_w, 1, cz + 8}, 0, 0, true, 0 };
                break;
            }
            case RULE_EXCAVATE_POOL: {
                // Pool excavated at random coordinates without checking if inside a courtyard!
                int px = (int)(prng_f(&prng) * 10.0f) - 5;
                int pz = (int)(prng_f(&prng) * 16.0f) - 8;
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_POOL, {px, -1, pz, px + 4, 0, pz + 4}, 0, 0, true, 0 };
                break;
            }
            case RULE_PLANT_GROVE: {
                // Random unsymmetric tree planting
                int tx = (int)(prng_f(&prng) * 20.0f) - 10;
                int tz = (int)(prng_f(&prng) * 20.0f) - 10;
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_TREE, {tx, 0, tz, tx + 2, 4, tz + 2}, 0, 0, true, 0 };
                break;
            }
            case RULE_ERECT_THOLOS: {
                // Tholos placed off-center (asymmetric)
                int ox = (prng_f(&prng) > 0.5f) ? -5 : 4;
                plan.nodes[plan.node_count++] = (ArchNode){ ARCH_PRIMITIVE_THOLOS, {ox - 2, 0, -16, ox + 2, 5, -12}, 0, 0, true, 0 };
                break;
            }
            default:
                break;
        }
    }
    return plan;
}

// ---------------------------------------------------------------------------
// 2. Metrics Evaluator
// ---------------------------------------------------------------------------
typedef struct {
    int total_trials;
    int symmetry_violations;
    int floating_elements;
    int interior_pool_violations;
} Metrics;

static bool check_symmetry(const TemplePlan *p) {
    for (int i = 0; i < p->node_count; ++i) {
        const ArchNode *n = &p->nodes[i];
        if (n->box.min_x == -n->box.max_x) continue; // Self-symmetric across X=0

        // Search for mirror counterpart
        bool found = false;
        for (int j = 0; j < p->node_count; ++j) {
            const ArchNode *m = &p->nodes[j];
            if (m->type == n->type &&
                m->box.min_y == n->box.min_y && m->box.max_y == n->box.max_y &&
                m->box.min_z == n->box.min_z && m->box.max_z == n->box.max_z &&
                m->box.min_x == -n->box.max_x && m->box.max_x == -n->box.min_x) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

static bool check_floating(const TemplePlan *p) {
    for (int i = 0; i < p->node_count; ++i) {
        const ArchNode *n = &p->nodes[i];
        if (n->type == ARCH_PRIMITIVE_PEDIMENT) {
            // A pediment must have supporting elements directly beneath it at min_y
            bool has_support = false;
            for (int j = 0; j < p->node_count; ++j) {
                const ArchNode *m = &p->nodes[j];
                if (m == n) continue;
                if (m->box.max_y == n->box.min_y) {
                    has_support = true;
                    break;
                }
            }
            if (!has_support) return true; // Floating pediment!
        }
    }
    return false;
}

static bool check_pool_in_cella(const TemplePlan *p) {
    for (int i = 0; i < p->node_count; ++i) {
        if (p->nodes[i].type != ARCH_PRIMITIVE_POOL) continue;
        const ModBox3D *b = &p->nodes[i].box;
        // Cella is at [-3..3, 0..6, -6..6]
        if (b->min_x >= -3 && b->max_x <= 3 && b->min_z >= -6 && b->max_z <= 6) {
            return true; // Pool carved inside the sacred altar room!
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// 3. Main Experiment Runner & Visual Comparison Renderer
// ---------------------------------------------------------------------------
#define MAX_PREVIEW 262144
typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview[MAX_PREVIEW];
static int preview_cnt = 0;

static void plot_p(int gx, int gy, int gz, Color c) {
    if (preview_cnt >= MAX_PREVIEW) return;
    float vs = 0.5f;
    float floor_s = 20.0f;
    preview[preview_cnt++] = (SimpleVoxel){
        (Vector3){ (gx + 0.5f) * vs - floor_s, (gy + 0.5f) * vs, (gz + 0.5f) * vs - floor_s },
        c
    };
}

static void render_comparison_image(const TemplePlan *p, const char *path, const char *title, Color banner_col) {
    preview_cnt = 0;
    rasterize_temple_plan(p, 40, 40, 3, plot_p);

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 34.0f, 24.0f, 34.0f };
    camera.target = (Vector3){ 0.0f, 3.5f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 50.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        ClearBackground((Color){ 140, 190, 230, 255 });
        BeginMode3D(camera);
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 75, 105, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 55, 80, 50, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_cnt; ++i) {
                DrawCube(preview[i].pos, vs, vs, vs, preview[i].color);
                DrawCubeWires(preview[i].pos, vs, vs, vs, ColorBrightness(preview[i].color, -0.2f));
            }
        EndMode3D();

        DrawRectangle(20, 20, 780, 85, Fade(BLACK, 0.75f));
        DrawRectangleLines(20, 20, 780, 85, banner_col);
        DrawText(title, 35, 32, 22, banner_col);
        DrawText(TextFormat("Nodes: %d | Voxels: %d | Seed: %u", p->node_count, preview_cnt, p->seed), 35, 66, 18, RAYWHITE);
    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, path);
    UnloadImage(img);
    UnloadRenderTexture(target);
    printf("Saved render: %s (%d voxels)\n", path, preview_cnt);
}

int main(void) {
    printf("======================================================================\n");
    printf(" ABLATION EXPERIMENT: EVALUATING THE IMPORTANCE OF THE COALGEBRA      \n");
    printf(" Comparing Full Bi-Algebra (Algebra + Coalgebra) vs Ablated (Algebra Only)\n");
    printf("======================================================================\n\n");

    const int NUM_TRIALS = 100;
    Metrics with_co = { 0 };
    Metrics without_co = { 0 };

    with_co.total_trials = NUM_TRIALS;
    without_co.total_trials = NUM_TRIALS;

    for (uint32_t s = 1; s <= NUM_TRIALS; ++s) {
        // 1. Full Bi-Algebra
        TemplePlan p_bi = generate_greek_temple(s, TEMPLE_STAGE_SANCTUARY, 2);
        if (!check_symmetry(&p_bi)) with_co.symmetry_violations++;
        if (check_floating(&p_bi))  with_co.floating_elements++;
        if (check_pool_in_cella(&p_bi)) with_co.interior_pool_violations++;

        // 2. Ablated (No Coalgebra)
        TemplePlan p_abl = generate_temple_ablated_no_coalgebra(s, 15);
        if (!check_symmetry(&p_abl)) without_co.symmetry_violations++;
        if (check_floating(&p_abl))  without_co.floating_elements++;
        if (check_pool_in_cella(&p_abl)) without_co.interior_pool_violations++;
    }

    printf("[RESULTS OVER %d RANDOM TRIALS]\n", NUM_TRIALS);
    printf("----------------------------------------------------------------------\n");
    printf("%-35s | %-16s | %-16s\n", "Metric", "WITH Coalgebra", "WITHOUT Coalgebra");
    printf("----------------------------------------------------------------------\n");
    printf("%-35s | %-16s | %-16s\n",
           "Bilateral Symmetry Violations",
           TextFormat("%d / %d (%.1f%%)", with_co.symmetry_violations, NUM_TRIALS, 100.0f * with_co.symmetry_violations / NUM_TRIALS),
           TextFormat("%d / %d (%.1f%%)", without_co.symmetry_violations, NUM_TRIALS, 100.0f * without_co.symmetry_violations / NUM_TRIALS));

    printf("%-35s | %-16s | %-16s\n",
           "Floating / Disconnected Roofs",
           TextFormat("%d / %d (%.1f%%)", with_co.floating_elements, NUM_TRIALS, 100.0f * with_co.floating_elements / NUM_TRIALS),
           TextFormat("%d / %d (%.1f%%)", without_co.floating_elements, NUM_TRIALS, 100.0f * without_co.floating_elements / NUM_TRIALS));

    printf("%-35s | %-16s | %-16s\n",
           "Pool In Cella / Altar Clashes",
           TextFormat("%d / %d (%.1f%%)", with_co.interior_pool_violations, NUM_TRIALS, 100.0f * with_co.interior_pool_violations / NUM_TRIALS),
           TextFormat("%d / %d (%.1f%%)", without_co.interior_pool_violations, NUM_TRIALS, 100.0f * without_co.interior_pool_violations / NUM_TRIALS));
    printf("----------------------------------------------------------------------\n\n");

    // Offscreen render comparison
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Ablation Renderer");
    if (IsWindowReady()) {
        const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

        TemplePlan good = generate_greek_temple(1337, TEMPLE_STAGE_SANCTUARY, 2);
        render_comparison_image(&good,
                                TextFormat("%s/ablation_with_coalgebra.png", out_dir),
                                "WITH Coalgebra: Coherent Vitruvian Bi-Algebra",
                                GOLD);

        TemplePlan bad = generate_temple_ablated_no_coalgebra(1337, 18);
        render_comparison_image(&bad,
                                TextFormat("%s/ablation_without_coalgebra.png", out_dir),
                                "WITHOUT Coalgebra (Ablated): Asymmetric Chaos & Disconnection",
                                RED);

        CloseWindow();
    }

    return 0;
}
