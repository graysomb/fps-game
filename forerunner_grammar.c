#include "forerunner_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// PRNG helper
static inline float forerunner_prng(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return (float)(*state & 0x00FFFFFF) / (float)0x01000000;
}

static inline ForerunnerBox3D make_forerunner_box(int x0, int y0, int z0, int x1, int y1, int z1) {
    return (ForerunnerBox3D){ x0, y0, z0, x1, y1, z1 };
}

static int forerunner_add_node(ForerunnerPlan *plan, ForerunnerPrimitiveType type, ForerunnerMaterial mat,
                               ForerunnerBox3D box, const int *parents, int parent_count,
                               float cant_angle, bool is_void, float mass_weight) {
    if (plan->node_count >= MAX_FORERUNNER_NODES) return -1;
    int id = plan->node_count++;
    ForerunnerNode *n = &plan->nodes[id];
    n->type = type;
    n->material = mat;
    n->box = box;
    n->cant_angle_deg = cant_angle;
    n->is_void = is_void;
    n->mass_weight = mass_weight;
    n->parent_count = 0;
    if (parents && parent_count > 0) {
        int count = parent_count < MAX_FORERUNNER_PARENTS ? parent_count : MAX_FORERUNNER_PARENTS;
        for (int i = 0; i < count; ++i) {
            n->parent_ids[i] = parents[i];
        }
        n->parent_count = count;
    }
    return id;
}

// ---------------------------------------------------------------------------
// Forerunner Coalgebra γ(S): Exposing Mass/Void Polarity Continuations
// ---------------------------------------------------------------------------
int forerunner_coalgebra_inspect(const ForerunnerPlan *plan, ForerunnerOpportunity *out_opps, int max_opps) {
    if (!plan || !out_opps || max_opps <= 0) return 0;
    int count = 0;

    // 1. Polarity Genesis (-): The Abyssal Negative Void
    if (!plan->has_chasm && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_EXCAVATE_CHASM;
        opp->primitive_type = FORERUNNER_PRIM_CHASM;
        opp->material = FORERUNNER_MAT_DARK_BASALT;
        opp->target_box = make_forerunner_box(-plan->chasm_width / 2, -plan->chasm_depth, -plan->chasm_length / 2,
                                              plan->chasm_width / 2, 0, plan->chasm_length / 2);
        opp->parent_count = 0;
        opp->score = 1.0f;
    }

    // 2. Polarity Response (+): Canted Chevron Monumental Gateway
    if (plan->has_chasm && !plan->has_gateway && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY;
        opp->primitive_type = FORERUNNER_PRIM_PYLON;
        opp->material = FORERUNNER_MAT_PEWTER;
        opp->target_box = make_forerunner_box(-10, 0, 14, 10, 15, 18);
        opp->parent_count = 0;
        opp->score = 0.95f;
    }

    // 3. Void Response (- to +): Suspended Skybridge with Hard-Light Span
    if (plan->has_gateway && !plan->has_bridge && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_SPAN_SKYBRIDGE;
        opp->primitive_type = FORERUNNER_PRIM_HARDLIGHT_BRIDGE;
        opp->material = FORERUNNER_MAT_HARDLIGHT_CYAN;
        opp->target_box = make_forerunner_box(-2, 2, -10, 2, 3, 14);
        opp->parent_count = 0;
        opp->score = 0.90f;
    }

    // 4. Mass Enclosure (+): Subterranean Terraced Vault Galleries
    if (plan->has_bridge && !plan->has_vault && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_TERRACE_VAULT;
        opp->primitive_type = FORERUNNER_PRIM_TERRACE;
        opp->material = FORERUNNER_MAT_PEWTER;
        opp->target_box = make_forerunner_box(-14, 0, -8, 14, 8, 8);
        opp->parent_count = 0;
        opp->score = 0.85f;
    }

    // 5. Apex Core (+ & -): Hovering Gravity Matrix & Telemetry Spire Pair
    if (plan->has_vault && !plan->has_core && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_RAISE_GRAVITY_CORE;
        opp->primitive_type = FORERUNNER_PRIM_GRAVITY_CORE;
        opp->material = FORERUNNER_MAT_HARDLIGHT_CYAN;
        opp->target_box = make_forerunner_box(-5, 0, -12, 5, 18, -8);
        opp->parent_count = 0;
        opp->score = 0.80f;
    }

    // 6. Detail Scale: Recessed Emissive Cyan Hard-Light Conduits
    if (plan->has_core && plan->light_channel_count < 4 && count < max_opps) {
        ForerunnerOpportunity *opp = &out_opps[count++];
        opp->rule = FORERUNNER_RULE_CARVE_LIGHT_CHANNELS;
        opp->primitive_type = FORERUNNER_PRIM_LIGHT_CHANNEL;
        opp->material = FORERUNNER_MAT_HARDLIGHT_CYAN;
        opp->target_box = make_forerunner_box(-10, 0, -12, 10, 16, 16);
        opp->parent_count = 0;
        opp->score = 0.75f;
    }

    return count;
}

// ---------------------------------------------------------------------------
// Forerunner Algebra α(S, o): Applying Motifs & Coupling Structures
// ---------------------------------------------------------------------------
bool forerunner_algebra_apply(ForerunnerPlan *plan, const ForerunnerOpportunity *opp) {
    if (!plan || !opp) return false;

    switch (opp->rule) {
        case FORERUNNER_RULE_EXCAVATE_CHASM: {
            // Negative Void Trench: X in [-6..6], Z in [-16..16], Y in [-8..0]
            int w = plan->chasm_width / 2; // 6
            int l = plan->chasm_length / 2; // 16
            int d = plan->chasm_depth; // 8

            // 1. Excavated Abyssal Void box (clears terrain voxels in rasterizer)
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d, -l, w, 0, l),
                                NULL, 0, 0.0f, true, 0.0f);

            // 2. Basalt bedrock floor trench at Y = -d
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d - 1, -l, w, -d, l),
                                NULL, 0, 0.0f, false, 1.0f);

            // 3. Primordial Floating Sentinel Beacon Spire at rear abyss
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-1, -d + 2, -l + 3, 1, 4, -l + 5),
                                NULL, 0, 0.0f, false, 0.5f);

            plan->has_chasm = true;
            if (plan->stage < FORERUNNER_STAGE_CHASM) plan->stage = FORERUNNER_STAGE_CHASM;
            return true;
        }

        case FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY: {
            // Monumental Entrance Gateway: Paired Canted Pylons (25 deg rake) & Chevron Lintel
            // Left Canted Pylon (X in [-9.. -4], Z in [14..18], Y in [0..14])
            int pylon_l = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(-9, 0, 14, -4, 14, 18),
                                              NULL, 0, 25.0f, false, 0.8f);

            // Right Canted Pylon (Reflected across X = 0: X in [4..9])
            int pylon_r = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(4, 0, 14, 9, 14, 18),
                                              NULL, 0, -25.0f, false, 0.8f);

            // Heavy Chevron Lintel bridging pylons at Y in [12..15]
            int parents[2] = { pylon_l, pylon_r };
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-8, 12, 14, 8, 15, 18),
                                parents, 2, 0.0f, false, 0.6f);

            // Glowing hard-light aperture border
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-4, 0, 16, -3, 12, 17),
                                parents, 1, 0.0f, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(3, 0, 16, 4, 12, 17),
                                parents, 1, 0.0f, false, 0.1f);

            plan->pylon_count += 2;
            plan->has_gateway = true;
            if (plan->stage < FORERUNNER_STAGE_GATEWAY) plan->stage = FORERUNNER_STAGE_GATEWAY;
            return true;
        }

        case FORERUNNER_RULE_SPAN_SKYBRIDGE: {
            // Suspended Skybridge across the chasm at Y = 2
            // 1. South Metallic Cantilever Ramp (Z in [8..14])
            int ramp_s = forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                             make_forerunner_box(-2, 2, 8, 2, 3, 14),
                                             NULL, 0, 0.0f, false, 0.7f);

            // 2. North Metallic Cantilever Ramp (Z in [-10.. -4])
            int ramp_n = forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                             make_forerunner_box(-2, 2, -10, 2, 3, -4),
                                             NULL, 0, 0.0f, false, 0.7f);

            // 3. Central Pure Hard-Light Energy Bridge (Z in [-4..8])
            int p_bridge[2] = { ramp_s, ramp_n };
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 2, -4, 2, 3, 8),
                                p_bridge, 2, 0.0f, false, 0.3f);

            // 4. Glowing Cyan Guide-Rails along the flanks (X = -2 and X = 2)
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-3, 3, -10, -2, 4, 14),
                                p_bridge, 2, 0.0f, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(2, 3, -10, 3, 4, 14),
                                p_bridge, 2, 0.0f, false, 0.1f);

            plan->has_bridge = true;
            if (plan->stage < FORERUNNER_STAGE_SKYBRIDGE) plan->stage = FORERUNNER_STAGE_SKYBRIDGE;
            return true;
        }

        case FORERUNNER_RULE_TERRACE_VAULT: {
            // Subterranean Terraced Vault Galleries flanking the chasm
            // Left Wall Galleries (X in [-14.. -6], Z in [-8..8])
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-14, 0, -8, -6, 2, 8), NULL, 0, 0.0f, false, 0.6f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-14, 2, -6, -8, 5, 6), NULL, 0, 0.0f, false, 0.6f);

            // Right Wall Galleries (Bilateral Reflection: X in [6..14], Z in [-8..8])
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(6, 0, -8, 14, 2, 8), NULL, 0, 0.0f, false, 0.6f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(8, 2, -6, 14, 5, 6), NULL, 0, 0.0f, false, 0.6f);

            // Canted Stabilizing Wall Buttresses
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-12, 0, -2, -6, 8, 2), NULL, 0, 20.0f, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(6, 0, -2, 12, 8, 2), NULL, 0, -20.0f, false, 0.8f);

            plan->pylon_count += 2;
            plan->has_vault = true;
            if (plan->stage < FORERUNNER_STAGE_VAULT) plan->stage = FORERUNNER_STAGE_VAULT;
            return true;
        }

        case FORERUNNER_RULE_RAISE_GRAVITY_CORE: {
            // Apex Monument: Hovering Installation Core & Twin Telemetry Spires
            // 1. Twin Towering Telemetry Spires at Z = -12, X = +-5 (Y in [0..18])
            int spire_l = forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(-6, 0, -13, -4, 18, -11),
                                              NULL, 0, 15.0f, false, 0.9f);
            int spire_r = forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(4, 0, -13, 6, 18, -11),
                                              NULL, 0, -15.0f, false, 0.9f);

            // Solar Amber Telemetry Emitters at Spire Tips
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(-6, 17, -13, -4, 19, -11),
                                NULL, 0, 0.0f, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(4, 17, -13, 6, 19, -11),
                                NULL, 0, 0.0f, false, 0.1f);

            // 2. Central Hovering Gravity Matrix Core (Suspended in mid-air at X in [-2..2], Y in [6..10], Z in [-13.. -9])
            int p_spires[2] = { spire_l, spire_r };
            forerunner_add_node(plan, FORERUNNER_PRIM_GRAVITY_CORE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 6, -13, 2, 10, -9),
                                p_spires, 2, 0.0f, false, 0.4f);

            // Core energy containment ring
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-3, 7, -14, 3, 9, -8),
                                p_spires, 2, 0.0f, false, 0.3f);

            plan->has_core = true;
            if (plan->stage < FORERUNNER_STAGE_CARTOGRAPHER) plan->stage = FORERUNNER_STAGE_CARTOGRAPHER;
            return true;
        }

        case FORERUNNER_RULE_CARVE_LIGHT_CHANNELS: {
            // Emissive Cyan Grooves carved along the length of the vault walls
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-7, 1, -8, -6, 2, 8), NULL, 0, 0.0f, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(6, 1, -8, 7, 2, 8), NULL, 0, 0.0f, false, 0.1f);

            // Bronze accent panels around portal base
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_BRONZE_ACCENT,
                                make_forerunner_box(-10, 0, 13, -9, 3, 15), NULL, 0, 0.0f, false, 0.2f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_BRONZE_ACCENT,
                                make_forerunner_box(9, 0, 13, 10, 3, 15), NULL, 0, 0.0f, false, 0.2f);

            plan->light_channel_count += 4;
            plan->stage = FORERUNNER_STAGE_CARTOGRAPHER;
            return true;
        }

        default:
            return false;
    }
}

// ---------------------------------------------------------------------------
// Developmental Growth Loop
// ---------------------------------------------------------------------------
ForerunnerPlan generate_forerunner_structure(uint32_t seed, ForerunnerStage target_stage) {
    ForerunnerPlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.stage = FORERUNNER_STAGE_CHASM;

    plan.chasm_width  = 12; // 6 meters half-width
    plan.chasm_depth  = 8;  // 4 meters deep drop
    plan.chasm_length = 32; // 16 meters half-length

    uint32_t prng = seed;
    int max_steps = 20;

    while (max_steps-- > 0) {
        if (plan.stage >= target_stage) {
            if (target_stage == FORERUNNER_STAGE_CHASM && plan.has_chasm) break;
            if (target_stage == FORERUNNER_STAGE_GATEWAY && plan.has_gateway) break;
            if (target_stage == FORERUNNER_STAGE_SKYBRIDGE && plan.has_bridge) break;
            if (target_stage == FORERUNNER_STAGE_VAULT && plan.has_vault) break;
            if (target_stage == FORERUNNER_STAGE_CARTOGRAPHER && plan.has_core) break;
        }

        ForerunnerOpportunity opps[MAX_FORERUNNER_OPPORTUNITIES];
        int opp_count = forerunner_coalgebra_inspect(&plan, opps, MAX_FORERUNNER_OPPORTUNITIES);
        if (opp_count == 0) break;

        // Probabilistic selection
        float total_weight = 0.0f;
        for (int i = 0; i < opp_count; ++i) {
            total_weight += opps[i].score;
        }

        float r = forerunner_prng(&prng) * total_weight;
        float accum = 0.0f;
        int chosen_idx = 0;
        for (int i = 0; i < opp_count; ++i) {
            accum += opps[i].score;
            if (r <= accum) {
                chosen_idx = i;
                break;
            }
        }

        forerunner_algebra_apply(&plan, &opps[chosen_idx]);
    }

    return plan;
}

// ---------------------------------------------------------------------------
// ASCII Blueprint Renderer
// ---------------------------------------------------------------------------
void print_forerunner_blueprint(const ForerunnerPlan *plan) {
    if (!plan) return;

    printf("\n=======================================================\n");
    printf(" FORERUNNER INSTALLATION BLUEPRINT (Nodes: %d, Seed: %u)\n", plan->node_count, plan->seed);
    printf(" Key: [X] Canted Pylon, [=] Lintel / Ramp, [#] Hard-Light Bridge,\n");
    printf("      [:] Light Channel, [!] Telemetry Spire, [O] Gravity Core, [ ] Chasm Void\n");
    printf("=======================================================\n");

    int min_coord = -16;
    int max_coord = 20;
    int step = 2;

    for (int z = max_coord; z >= min_coord; z -= step) {
        printf("%3d | ", z);
        for (int x = min_coord; x <= max_coord; x += step) {
            char ch = '.';

            int top_y = -99;
            for (int i = 0; i < plan->node_count; ++i) {
                const ForerunnerNode *n = &plan->nodes[i];
                if (x >= n->box.min_x && x < n->box.max_x &&
                    z >= n->box.min_z && z < n->box.max_z) {
                    if (n->is_void) {
                        if (top_y == -99) ch = ' '; // Chasm void
                    } else if (n->box.max_y > top_y) {
                        top_y = n->box.max_y;
                        if (n->type == FORERUNNER_PRIM_PYLON) ch = 'X';
                        else if (n->type == FORERUNNER_PRIM_LINTEL) ch = '=';
                        else if (n->type == FORERUNNER_PRIM_HARDLIGHT_BRIDGE) ch = '#';
                        else if (n->type == FORERUNNER_PRIM_BRIDGE_SPAN) ch = '=';
                        else if (n->type == FORERUNNER_PRIM_SPIRE) ch = '!';
                        else if (n->type == FORERUNNER_PRIM_GRAVITY_CORE) ch = 'O';
                        else if (n->type == FORERUNNER_PRIM_LIGHT_CHANNEL) ch = ':';
                        else if (n->type == FORERUNNER_PRIM_TERRACE) ch = '-';
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
