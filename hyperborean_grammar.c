#include "hyperborean_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// PRNG helper
static inline float hyper_prng(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return (float)(*state & 0x00FFFFFF) / (float)0x01000000;
}

static inline HyperBox3D make_hyper_box(int x0, int y0, int z0, int x1, int y1, int z1) {
    return (HyperBox3D){ x0, y0, z0, x1, y1, z1 };
}

static int hyper_add_node(HyperPlan *plan, HyperPrimitiveType type, HyperMaterial mat,
                          HyperBox3D box, const int *parents, int parent_count, float weight) {
    if (plan->node_count >= MAX_HYPER_NODES) return -1;
    int id = plan->node_count++;
    HyperNode *n = &plan->nodes[id];
    n->type = type;
    n->material = mat;
    n->box = box;
    n->stability_weight = weight;
    n->parent_count = 0;
    if (parents && parent_count > 0) {
        int count = parent_count < MAX_HYPER_PARENTS ? parent_count : MAX_HYPER_PARENTS;
        for (int i = 0; i < count; ++i) {
            n->parent_ids[i] = parents[i];
        }
        n->parent_count = count;
    }
    return id;
}

// ---------------------------------------------------------------------------
// Syncretic Coalgebra γ(S)
// ---------------------------------------------------------------------------
int hyperborean_coalgebra_inspect(const HyperPlan *plan, HyperOpportunity *out_opps, int max_opps) {
    if (!plan || !out_opps || max_opps <= 0) return 0;
    int count = 0;

    // 1. Solstice Processional Avenue & Heel Stone Marker
    if (!plan->has_avenue && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_ALIGN_SOLSTICE_AVENUE;
        opp->primitive_type = HYPER_PRIM_HEEL_STONE;
        opp->material = HYPER_MAT_WEATHERED_SARSEN;
        opp->target_box = make_hyper_box(-1, 0, 24, 1, 7, 26);
        opp->parent_count = 0;
        opp->score = 1.0f;
    }

    // 2. Outer Cromlech Megalithic Ring (16 Sarsen Menhirs + Curved Lintels)
    if (plan->outer_stones_placed < 16 && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_ERECT_OUTER_HENGE;
        opp->primitive_type = HYPER_PRIM_MENHIR;
        opp->material = HYPER_MAT_WEATHERED_SARSEN;
        opp->target_box = make_hyper_box(-plan->r_outer_henge, 0, -plan->r_outer_henge,
                                         plan->r_outer_henge, 6, plan->r_outer_henge);
        opp->parent_count = 0;
        opp->score = 0.95f;
    }

    // 3. Middle Classical Greek Peristyle Ring (12 Fluted Marble Columns & Entablature)
    if (plan->outer_stones_placed >= 16 && plan->peristyle_cols_placed < 12 && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_ERECT_MARBLE_PERISTYLE;
        opp->primitive_type = HYPER_PRIM_FLUTED_COLUMN;
        opp->material = HYPER_MAT_PENTELIC_MARBLE;
        opp->target_box = make_hyper_box(-plan->r_mid_peristyle, 0, -plan->r_mid_peristyle,
                                         plan->r_mid_peristyle, 6, plan->r_mid_peristyle);
        opp->parent_count = 0;
        opp->score = 0.90f;
    }

    // 4. Inner Colossal Hybrid Trilithons with Classical Pediment
    if (plan->peristyle_cols_placed >= 12 && plan->trilithons_placed < 3 && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_RAISE_HYBRID_TRILITHONS;
        opp->primitive_type = HYPER_PRIM_HYBRID_TRILITHON;
        opp->material = HYPER_MAT_PENTELIC_MARBLE;
        opp->target_box = make_hyper_box(-plan->r_inner_trilithon, 0, -plan->r_inner_trilithon,
                                         plan->r_inner_trilithon, 9, plan->r_inner_trilithon);
        opp->parent_count = 0;
        opp->score = 0.85f;
    }

    // 5. Epicenter: Sunken Marble Tholos, PBF Moat & Bronze Fire Brazier
    if (plan->trilithons_placed >= 3 && !plan->has_tholos && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_EXCAVATE_THOLOS_BASIN;
        opp->primitive_type = HYPER_PRIM_THOLOS_PODIUM;
        opp->material = HYPER_MAT_PENTELIC_MARBLE;
        opp->target_box = make_hyper_box(-plan->r_tholos_core, 0, -plan->r_tholos_core,
                                         plan->r_tholos_core, 2, plan->r_tholos_core);
        opp->parent_count = 0;
        opp->score = 0.80f;
    }

    // 6. Sacred Grove & Statuary
    if (plan->has_tholos && !plan->has_moat && count < max_opps) {
        HyperOpportunity *opp = &out_opps[count++];
        opp->rule = HYPER_RULE_POPULATE_SACRED_GROVE;
        opp->primitive_type = HYPER_PRIM_CYPRESS;
        opp->material = HYPER_MAT_CYPRESS;
        opp->target_box = make_hyper_box(-18, 0, -18, 18, 8, 28);
        opp->parent_count = 0;
        opp->score = 0.75f;
    }

    return count;
}

// ---------------------------------------------------------------------------
// Syncretic Algebra α(S, r_t)
// ---------------------------------------------------------------------------
bool hyperborean_algebra_apply(HyperPlan *plan, const HyperOpportunity *opp) {
    if (!plan || !opp) return false;

    switch (opp->rule) {
        case HYPER_RULE_ALIGN_SOLSTICE_AVENUE: {
            // 1. Colossal Heel Stone Monolith marking midsummer sunrise
            hyper_add_node(plan, HYPER_PRIM_HEEL_STONE, HYPER_MAT_WEATHERED_SARSEN,
                           make_hyper_box(-1, 0, 24, 1, 7, 26), NULL, 0, 0.5f);

            // 2. Flanking Avenue Megaliths along Z in [18..24] at X = +-5
            for (int z = 18; z <= 24; z += 3) {
                hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                               make_hyper_box(-5, 0, z, -4, 4, z + 1), NULL, 0, 0.4f);
                hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                               make_hyper_box(4, 0, z, 5, 4, z + 1), NULL, 0, 0.4f);

                // Flanking Classical Votive Statues / Herms
                hyper_add_node(plan, HYPER_PRIM_STATUE, HYPER_MAT_PENTELIC_MARBLE,
                               make_hyper_box(-3, 0, z, -2, 2, z + 1), NULL, 0, 0.2f);
                hyper_add_node(plan, HYPER_PRIM_STATUE, HYPER_MAT_PENTELIC_MARBLE,
                               make_hyper_box(2, 0, z, 3, 2, z + 1), NULL, 0, 0.2f);
            }

            plan->has_avenue = true;
            if (plan->stage < HYPER_STAGE_AVENUE) plan->stage = HYPER_STAGE_AVENUE;
            return true;
        }

        case HYPER_RULE_ERECT_OUTER_HENGE: {
            // Outer Cromlech Ring: 16 massive sarsen menhirs with bridging curved lintels
            int r = plan->r_outer_henge; // 17
            int num_stones = 16;
            int menhir_ids[16];

            for (int i = 0; i < num_stones; ++i) {
                // Skip opening towards the solstice avenue at +Z
                float angle = (2.0f * (float)M_PI * (float)i) / (float)num_stones;
                int mx = (int)roundf(r * cosf(angle));
                int mz = (int)roundf(r * sinf(angle));

                if (mz > r - 3 && abs(mx) <= 4) {
                    menhir_ids[i] = -1;
                    continue; // Gateway portal threshold
                }

                menhir_ids[i] = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                               make_hyper_box(mx - 1, 0, mz - 1, mx + 1, 5, mz + 1),
                                               NULL, 0, 0.5f);
            }

            // Bridging curved lintels spanning adjacent menhirs
            for (int i = 0; i < num_stones; ++i) {
                int next = (i + 1) % num_stones;
                if (menhir_ids[i] >= 0 && menhir_ids[next] >= 0) {
                    float a1 = (2.0f * (float)M_PI * (float)i) / (float)num_stones;
                    float a2 = (2.0f * (float)M_PI * (float)next) / (float)num_stones;
                    int x1 = (int)roundf(r * cosf(a1));
                    int z1 = (int)roundf(r * sinf(a1));
                    int x2 = (int)roundf(r * cosf(a2));
                    int z2 = (int)roundf(r * sinf(a2));

                    int lx0 = (x1 < x2 ? x1 : x2) - 1;
                    int lx1 = (x1 > x2 ? x1 : x2) + 1;
                    int lz0 = (z1 < z2 ? z1 : z2) - 1;
                    int lz1 = (z1 > z2 ? z1 : z2) + 1;

                    int parents[2] = { menhir_ids[i], menhir_ids[next] };
                    hyper_add_node(plan, HYPER_PRIM_ROUGH_LINTEL, HYPER_MAT_WEATHERED_SARSEN,
                                   make_hyper_box(lx0, 5, lz0, lx1, 6, lz1),
                                   parents, 2, 0.4f);
                }
            }

            plan->outer_stones_placed = 16;
            if (plan->stage < HYPER_STAGE_OUTER_HENGE) plan->stage = HYPER_STAGE_OUTER_HENGE;
            return true;
        }

        case HYPER_RULE_ERECT_MARBLE_PERISTYLE: {
            // Middle Concentric Ring: 12 Fluted White Pentelic Marble Doric Columns
            int r = plan->r_mid_peristyle; // 11
            int num_cols = 12;
            int col_ids[12];

            for (int i = 0; i < num_cols; ++i) {
                // Keep solstice entrance opening at +Z
                float angle = (2.0f * (float)M_PI * (float)i) / (float)num_cols;
                int cx = (int)roundf(r * cosf(angle));
                int cz = (int)roundf(r * sinf(angle));

                if (cz > r - 2 && abs(cx) <= 3) {
                    col_ids[i] = -1;
                    continue; // Avenue portal
                }

                col_ids[i] = hyper_add_node(plan, HYPER_PRIM_FLUTED_COLUMN, HYPER_MAT_PENTELIC_MARBLE,
                                            make_hyper_box(cx, 0, cz, cx + 1, 5, cz + 1),
                                            NULL, 0, 0.45f);
            }

            // Polished White Marble Architrave & Frieze Lintel Beams
            for (int i = 0; i < num_cols; ++i) {
                int next = (i + 1) % num_cols;
                if (col_ids[i] >= 0 && col_ids[next] >= 0) {
                    float a1 = (2.0f * (float)M_PI * (float)i) / (float)num_cols;
                    float a2 = (2.0f * (float)M_PI * (float)next) / (float)num_cols;
                    int x1 = (int)roundf(r * cosf(a1));
                    int z1 = (int)roundf(r * sinf(a1));
                    int x2 = (int)roundf(r * cosf(a2));
                    int z2 = (int)roundf(r * sinf(a2));

                    int lx0 = (x1 < x2 ? x1 : x2);
                    int lx1 = (x1 > x2 ? x1 : x2) + 1;
                    int lz0 = (z1 < z2 ? z1 : z2);
                    int lz1 = (z1 > z2 ? z1 : z2) + 1;

                    int parents[2] = { col_ids[i], col_ids[next] };
                    hyper_add_node(plan, HYPER_PRIM_ARCHITRAVE, HYPER_MAT_PENTELIC_MARBLE,
                                   make_hyper_box(lx0, 5, lz0, lx1, 6, lz1),
                                   parents, 2, 0.4f);
                }
            }

            plan->peristyle_cols_placed = 12;
            if (plan->stage < HYPER_STAGE_CLASSICAL_PERISTYLE) plan->stage = HYPER_STAGE_CLASSICAL_PERISTYLE;
            return true;
        }

        case HYPER_RULE_RAISE_HYBRID_TRILITHONS: {
            // Inner Horseshoe of Colossal Hybrid Trilithons
            // Trilithon 1: Central Monumental Pedimented Trilithon (At Winter Solstice Sunset Z = -6)
            int left1  = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(-3, 0, -6, -1, 7, -4), NULL, 0, 0.5f);
            int right1 = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(1, 0, -6, 3, 7, -4), NULL, 0, 0.5f);
            int p1[2] = { left1, right1 };

            // Classical Doric Architrave
            int arch1 = hyper_add_node(plan, HYPER_PRIM_CYCLOPEAN_ARCHITRAVE, HYPER_MAT_PENTELIC_MARBLE,
                                       make_hyper_box(-4, 7, -6, 4, 8, -4), p1, 2, 0.5f);
            int p_ped[3] = { left1, right1, arch1 };
            // Classical Carved Doric Triangular Pediment with Gold-leaf Tympanum
            hyper_add_node(plan, HYPER_PRIM_PEDIMENT, HYPER_MAT_GOLD_ACCENT,
                           make_hyper_box(-4, 8, -6, 4, 10, -4), p_ped, 3, 0.6f);

            // Trilithon 2: Left Flank Hybrid Trilithon (X = -6, Z in [-2..2])
            int left2  = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(-6, 0, -3, -4, 6, -1), NULL, 0, 0.45f);
            int right2 = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(-6, 0, 1, -4, 6, 3), NULL, 0, 0.45f);
            int p2[2] = { left2, right2 };
            hyper_add_node(plan, HYPER_PRIM_ARCHITRAVE, HYPER_MAT_PENTELIC_MARBLE,
                           make_hyper_box(-6, 6, -3, -4, 7, 3), p2, 2, 0.4f);

            // Trilithon 3: Right Flank Hybrid Trilithon (X = +6, Z in [-2..2])
            int left3  = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(4, 0, -3, 6, 6, -1), NULL, 0, 0.45f);
            int right3 = hyper_add_node(plan, HYPER_PRIM_MENHIR, HYPER_MAT_WEATHERED_SARSEN,
                                        make_hyper_box(4, 0, 1, 6, 6, 3), NULL, 0, 0.45f);
            int p3[2] = { left3, right3 };
            hyper_add_node(plan, HYPER_PRIM_ARCHITRAVE, HYPER_MAT_PENTELIC_MARBLE,
                           make_hyper_box(4, 6, -3, 6, 7, 3), p3, 2, 0.4f);

            plan->trilithons_placed = 3;
            if (plan->stage < HYPER_STAGE_GREAT_TRILITHONS) plan->stage = HYPER_STAGE_GREAT_TRILITHONS;
            return true;
        }

        case HYPER_RULE_EXCAVATE_THOLOS_BASIN: {
            // Epicenter: Sunken Circular Marble Tholos Podium surrounded by Sacred PBF Moat
            // int r_moat = plan->r_tholos_core; // 3
            // 1. Concentric Sacred PBF Reflection Fluid Moat (disabled for now)
            // hyper_add_node(plan, HYPER_PRIM_FLUID_BASIN, HYPER_MAT_PBF_WATER,
            //                make_hyper_box(-r_moat, 0, -r_moat, r_moat + 1, 1, r_moat + 1),
            //                NULL, 0, 0.1f);

            // 2. Central Raised Tiered Marble Tholos Altar Podium (R = 2)
            hyper_add_node(plan, HYPER_PRIM_THOLOS_PODIUM, HYPER_MAT_PENTELIC_MARBLE,
                           make_hyper_box(-2, 0, -2, 2, 2, 2), NULL, 0, 0.3f);

            // 3. Central Sacrificial Bronze Tripod with Eternal Sacred Fire Embers
            hyper_add_node(plan, HYPER_PRIM_BRAZIER, HYPER_MAT_EMBER,
                           make_hyper_box(-1, 2, -1, 1, 3, 1), NULL, 0, 0.2f);

            plan->has_tholos = true;
            plan->has_moat = true;
            plan->has_brazier = true;
            return true;
        }

        case HYPER_RULE_POPULATE_SACRED_GROVE: {
            // Sacred Grove: Symmetrical Italian Cypresses & Votive Statuary in inter-ring clearings
            int angles_deg[] = { 45, 135, 225, 315 };
            for (int i = 0; i < 4; ++i) {
                float rad = (float)angles_deg[i] * (float)M_PI / 180.0f;
                // Outer ring cypresses (R = 14)
                int ox = (int)roundf(14.0f * cosf(rad));
                int oz = (int)roundf(14.0f * sinf(rad));
                hyper_add_node(plan, HYPER_PRIM_CYPRESS, HYPER_MAT_CYPRESS,
                               make_hyper_box(ox, 0, oz, ox + 1, 8, oz + 1), NULL, 0, 0.1f);

                // Inner ring bronze votives (R = 9)
                int ix = (int)roundf(9.0f * cosf(rad));
                int iz = (int)roundf(9.0f * sinf(rad));
                hyper_add_node(plan, HYPER_PRIM_STATUE, HYPER_MAT_BRONZE,
                               make_hyper_box(ix, 0, iz, ix + 1, 2, iz + 1), NULL, 0, 0.2f);
            }

            plan->stage = HYPER_STAGE_FULL_SANCTUM;
            return true;
        }

        default:
            return false;
    }
}

// ---------------------------------------------------------------------------
// Developmental Growth Loop
// ---------------------------------------------------------------------------
HyperPlan generate_hyperborean_structure(uint32_t seed, HyperStage target_stage) {
    HyperPlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.stage = HYPER_STAGE_AVENUE;

    // Astronomical axis: Midsummer sunrise ~ +Z
    plan.axis_dx = 0.0f;
    plan.axis_dz = 1.0f;

    plan.r_outer_henge     = 17;
    plan.r_mid_peristyle   = 11;
    plan.r_inner_trilithon = 6;
    plan.r_tholos_core     = 3;

    uint32_t prng = seed;
    int max_steps = 25;

    while (max_steps-- > 0) {
        if (plan.stage >= target_stage) {
            if (target_stage == HYPER_STAGE_AVENUE && plan.has_avenue) break;
            if (target_stage == HYPER_STAGE_OUTER_HENGE && plan.outer_stones_placed >= 16) break;
            if (target_stage == HYPER_STAGE_CLASSICAL_PERISTYLE && plan.peristyle_cols_placed >= 12) break;
            if (target_stage == HYPER_STAGE_GREAT_TRILITHONS && plan.trilithons_placed >= 3) break;
            if (target_stage == HYPER_STAGE_FULL_SANCTUM && plan.has_tholos) break;
        }

        HyperOpportunity opps[MAX_HYPER_OPPORTUNITIES];
        int opp_count = hyperborean_coalgebra_inspect(&plan, opps, MAX_HYPER_OPPORTUNITIES);
        if (opp_count == 0) break;

        // Weighted probabilistic selection
        float total_weight = 0.0f;
        for (int i = 0; i < opp_count; ++i) {
            total_weight += opps[i].score;
        }

        float r = hyper_prng(&prng) * total_weight;
        float accum = 0.0f;
        int chosen_idx = 0;
        for (int i = 0; i < opp_count; ++i) {
            accum += opps[i].score;
            if (r <= accum) {
                chosen_idx = i;
                break;
            }
        }

        hyperborean_algebra_apply(&plan, &opps[chosen_idx]);
    }

    return plan;
}

// ---------------------------------------------------------------------------
// ASCII Blueprint Floorplan Renderer
// ---------------------------------------------------------------------------
void print_hyperborean_blueprint(const HyperPlan *plan) {
    if (!plan) return;

    printf("\n=======================================================\n");
    printf(" HYPERBOREAN SUN-HENGE BLUEPRINT (Nodes: %d, Seed: %u)\n", plan->node_count, plan->seed);
    printf(" Key: [O] Sarsen Menhir, [#] Marble Column, [=] Architrave / Lintel,\n");
    printf("      [^] Doric Pediment, [~] Sacred PBF Moat, [@] Tholos Altar,\n");
    printf("      [*] Cypress Tree, [!] Heel Stone Solstice Marker\n");
    printf("=======================================================\n");

    int min_coord = -20;
    int max_coord = 26;
    int step = 2;

    for (int z = max_coord; z >= min_coord; z -= step) {
        printf("%3d | ", z);
        for (int x = min_coord; x <= max_coord; x += step) {
            char ch = ' ';

            // Find top-most node at (x, z)
            int top_y = -1;
            for (int i = 0; i < plan->node_count; ++i) {
                const HyperNode *n = &plan->nodes[i];
                if (x >= n->box.min_x && x < n->box.max_x &&
                    z >= n->box.min_z && z < n->box.max_z) {
                    if (n->box.max_y > top_y) {
                        top_y = n->box.max_y;
                        if (n->type == HYPER_PRIM_HEEL_STONE) ch = '!';
                        else if (n->type == HYPER_PRIM_PEDIMENT) ch = '^';
                        else if (n->type == HYPER_PRIM_BRAZIER) ch = '@';
                        else if (n->type == HYPER_PRIM_FLUID_BASIN) ch = '~';
                        else if (n->type == HYPER_PRIM_FLUTED_COLUMN) ch = '#';
                        else if (n->type == HYPER_PRIM_ARCHITRAVE || n->type == HYPER_PRIM_CYCLOPEAN_ARCHITRAVE || n->type == HYPER_PRIM_ROUGH_LINTEL) ch = '=';
                        else if (n->type == HYPER_PRIM_MENHIR) ch = 'O';
                        else if (n->type == HYPER_PRIM_CYPRESS) ch = '*';
                        else if (n->type == HYPER_PRIM_STATUE) ch = 's';
                        else ch = '.';
                    }
                }
            }
            printf("%c ", ch);
        }
        printf("\n");
    }

    printf("      ");
    for (int x = min_coord; x <= max_coord; x += step * 2) {
        printf("----");
    }
    printf("\n      ");
    for (int x = min_coord; x <= max_coord; x += step * 4) {
        printf("%-8d", x);
    }
    printf("\n\n");
}
