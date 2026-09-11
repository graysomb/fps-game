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
                               float cant_x, float cant_z, bool is_void, float mass_weight) {
    if (plan->node_count >= MAX_FORERUNNER_NODES) return -1;
    int id = plan->node_count++;
    ForerunnerNode *n = &plan->nodes[id];
    n->type = type;
    n->material = mat;
    n->box = box;
    n->cant_angle_deg = cant_x;
    n->cant_angle_z = cant_z;
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

    switch (plan->archetype) {
        case FORERUNNER_ARCHETYPE_CARTOGRAPHER: {
            if (!plan->has_chasm && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_EXCAVATE_CHASM;
                opp->primitive_type = FORERUNNER_PRIM_CHASM;
                opp->score = 1.0f;
            }
            if (plan->has_chasm && !plan->has_gateway && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY;
                opp->primitive_type = FORERUNNER_PRIM_PYLON;
                opp->score = 0.95f;
            }
            if (plan->has_gateway && !plan->has_bridge && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_SPAN_SKYBRIDGE;
                opp->primitive_type = FORERUNNER_PRIM_HARDLIGHT_BRIDGE;
                opp->score = 0.90f;
            }
            if (plan->has_bridge && !plan->has_vault && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_TERRACE_VAULT;
                opp->primitive_type = FORERUNNER_PRIM_TERRACE;
                opp->score = 0.85f;
            }
            if (plan->has_vault && !plan->has_core && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_RAISE_GRAVITY_CORE;
                opp->primitive_type = FORERUNNER_PRIM_GRAVITY_CORE;
                opp->score = 0.80f;
            }
            if (plan->has_core && plan->light_channel_count < 4 && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_CARVE_LIGHT_CHANNELS;
                opp->primitive_type = FORERUNNER_PRIM_LIGHT_CHANNEL;
                opp->score = 0.75f;
            }
            break;
        }

        case FORERUNNER_ARCHETYPE_CROSSROADS: {
            if (!plan->has_chasm && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_EXCAVATE_CRUCIFORM_CHASM;
                opp->primitive_type = FORERUNNER_PRIM_CHASM;
                opp->score = 1.0f;
            }
            if (plan->has_chasm && !plan->has_gateway && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ERECT_CORNER_PYLONS;
                opp->primitive_type = FORERUNNER_PRIM_PYLON;
                opp->score = 0.95f;
            }
            if (plan->has_gateway && !plan->has_bridge && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_SPAN_CROSS_BRIDGES;
                opp->primitive_type = FORERUNNER_PRIM_HARDLIGHT_BRIDGE;
                opp->score = 0.90f;
            }
            if (plan->has_bridge && !plan->has_vault && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_TERRACE_VAULT;
                opp->primitive_type = FORERUNNER_PRIM_TERRACE;
                opp->score = 0.85f;
            }
            if (plan->has_vault && !plan->has_core && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ADD_INTERSECTION_NEXUS;
                opp->primitive_type = FORERUNNER_PRIM_GRAVITY_CORE;
                opp->score = 0.80f;
            }
            if (plan->has_core && plan->light_channel_count < 4 && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_CARVE_LIGHT_CHANNELS;
                opp->primitive_type = FORERUNNER_PRIM_LIGHT_CHANNEL;
                opp->score = 0.75f;
            }
            break;
        }

        case FORERUNNER_ARCHETYPE_CRUCIBLE: {
            if (!plan->has_chasm && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_EXCAVATE_GRAVITY_PIT;
                opp->primitive_type = FORERUNNER_PRIM_CHASM;
                opp->score = 1.0f;
            }
            if (plan->has_chasm && !plan->has_gateway && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ERECT_RADIAL_PYLON_RING;
                opp->primitive_type = FORERUNNER_PRIM_PYLON;
                opp->score = 0.95f;
            }
            if (plan->has_gateway && !plan->has_bridge && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_SPAN_RADIAL_CATWALK;
                opp->primitive_type = FORERUNNER_PRIM_HARDLIGHT_BRIDGE;
                opp->score = 0.90f;
            }
            if (plan->has_bridge && !plan->has_vault && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_TERRACE_VAULT;
                opp->primitive_type = FORERUNNER_PRIM_TERRACE;
                opp->score = 0.85f;
            }
            if (plan->has_vault && !plan->has_core && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_SUSPEND_CENTRIFUGE_CORE;
                opp->primitive_type = FORERUNNER_PRIM_GRAVITY_CORE;
                opp->score = 0.80f;
            }
            if (plan->has_core && plan->light_channel_count < 4 && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_CARVE_LIGHT_CHANNELS;
                opp->primitive_type = FORERUNNER_PRIM_LIGHT_CHANNEL;
                opp->score = 0.75f;
            }
            break;
        }

        case FORERUNNER_ARCHETYPE_SPIRE: {
            if (!plan->has_chasm && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_EXCAVATE_ANNULAR_MOAT;
                opp->primitive_type = FORERUNNER_PRIM_CHASM;
                opp->score = 1.0f;
            }
            if (plan->has_chasm && !plan->has_gateway && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ERECT_CENTRAL_PINNACLE;
                opp->primitive_type = FORERUNNER_PRIM_SPIRE;
                opp->score = 0.95f;
            }
            if (plan->has_gateway && !plan->has_bridge && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ADD_FLYING_BUTTRESSES;
                opp->primitive_type = FORERUNNER_PRIM_PYLON;
                opp->score = 0.90f;
            }
            if (plan->has_bridge && !plan->has_vault && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_ATTACH_OBSERVATION_DECKS;
                opp->primitive_type = FORERUNNER_PRIM_TERRACE;
                opp->score = 0.85f;
            }
            if (plan->has_vault && !plan->has_core && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_RAISE_GRAVITY_CORE;
                opp->primitive_type = FORERUNNER_PRIM_GRAVITY_CORE;
                opp->score = 0.80f;
            }
            if (plan->has_core && plan->light_channel_count < 4 && count < max_opps) {
                ForerunnerOpportunity *opp = &out_opps[count++];
                opp->rule = FORERUNNER_RULE_CARVE_LIGHT_CHANNELS;
                opp->primitive_type = FORERUNNER_PRIM_LIGHT_CHANNEL;
                opp->score = 0.75f;
            }
            break;
        }

        default:
            break;
    }

    return count;
}

// ---------------------------------------------------------------------------
// Forerunner Algebra α(S, o): Applying Motifs & Coupling Structures
// ---------------------------------------------------------------------------
bool forerunner_algebra_apply(ForerunnerPlan *plan, const ForerunnerOpportunity *opp) {
    if (!plan || !opp) return false;

    // High-entropy hash derived from seed for procedural variation
    uint32_t h = (plan->seed * 2654435761u) ^ (plan->seed >> 16);
    h ^= (h << 13); h ^= (h >> 17); h ^= (h << 5);

    int w = plan->chasm_width / 2;
    int l = plan->chasm_length / 2;
    int d = plan->chasm_depth;

    switch (opp->rule) {
        // ===================================================================
        // 1. Cartographer Rules
        // ===================================================================
        case FORERUNNER_RULE_EXCAVATE_CHASM: {
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d, -l, w, 0, l),
                                NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d - 1, -l, w, -d, l),
                                NULL, 0, 0, 0, false, 1.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-1, -d + 2, -l + 3, 1, 4, -l + 5),
                                NULL, 0, 0, 0, false, 0.5f);
            plan->has_chasm = true;
            if (plan->stage < FORERUNNER_STAGE_FOUNDATION) plan->stage = FORERUNNER_STAGE_FOUNDATION;
            return true;
        }

        case FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY: {
            float cant = 18.0f + ((h >> 8) % 4) * 4.0f;
            int h_pylon = 13 + ((h >> 10) % 3) * 2;
            int p_inner = w - 1;
            int p_outer = p_inner + 4 + ((h >> 12) % 2);

            int pylon_l = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(-p_outer, 0, 14, -p_inner, h_pylon, 18),
                                              NULL, 0, cant, 0, false, 0.8f);
            int pylon_r = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(p_inner, 0, 14, p_outer, h_pylon, 18),
                                              NULL, 0, -cant, 0, false, 0.8f);

            int y_lintel_bot = h_pylon - 2;
            int y_lintel_top = h_pylon + 1;
            int parents[2] = { pylon_l, pylon_r };
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-p_outer + 1, y_lintel_bot, 14, p_outer - 1, y_lintel_top, 18),
                                parents, 2, 0, 0, false, 0.6f);

            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-p_inner, 0, 16, -p_inner + 1, y_lintel_bot, 17),
                                parents, 1, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(p_inner - 1, 0, 16, p_inner, y_lintel_bot, 17),
                                parents, 1, 0, 0, false, 0.1f);

            plan->pylon_count += 2;
            plan->has_gateway = true;
            if (plan->stage < FORERUNNER_STAGE_SUPERSTRUCTURE) plan->stage = FORERUNNER_STAGE_SUPERSTRUCTURE;
            return true;
        }

        case FORERUNNER_RULE_SPAN_SKYBRIDGE: {
            int y_br = 2 + ((h >> 14) % 2) * 2;
            int bw = 2;
            int ramp_style = (h >> 15) % 2;
            int z_s_start = (ramp_style == 0) ? 8 : 10;
            int z_n_end   = (ramp_style == 0) ? -4 : -8;

            int ramp_s = forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                             make_forerunner_box(-bw, y_br, z_s_start, bw, y_br + 1, 14),
                                             NULL, 0, 0, 0, false, 0.7f);
            int ramp_n = forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                             make_forerunner_box(-bw, y_br, -10, bw, y_br + 1, z_n_end),
                                             NULL, 0, 0, 0, false, 0.7f);

            int p_bridge[2] = { ramp_s, ramp_n };
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-bw, y_br, z_n_end, bw, y_br + 1, z_s_start),
                                p_bridge, 2, 0, 0, false, 0.3f);

            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-bw, y_br + 1, -10, -bw + 1, y_br + 2, 14),
                                p_bridge, 2, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(bw - 1, y_br + 1, -10, bw, y_br + 2, 14),
                                p_bridge, 2, 0, 0, false, 0.1f);

            // Anchor bridge approaches firmly to bedrock with structural support pylons
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-bw, 0, -10, bw, y_br, -8),
                                NULL, 0, 0, 0, false, 0.9f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-bw, 0, 12, bw, y_br, 14),
                                NULL, 0, 0, 0, false, 0.9f);

            plan->has_bridge = true;
            if (plan->stage < FORERUNNER_STAGE_TRANSIT) plan->stage = FORERUNNER_STAGE_TRANSIT;
            return true;
        }

        case FORERUNNER_RULE_TERRACE_VAULT: {
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-w - 7, 0, -8, -w + 1, 2, 8), NULL, 0, 0, 0, false, 0.6f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-w - 7, 2, -6, -w - 1, 5, 6), NULL, 0, 0, 0, false, 0.6f);

            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(w - 1, 0, -8, w + 7, 2, 8), NULL, 0, 0, 0, false, 0.6f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(w + 1, 2, -6, w + 7, 5, 6), NULL, 0, 0, 0, false, 0.6f);

            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-w - 5, 0, -2, -w + 1, 8, 2), NULL, 0, 20.0f, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(w - 1, 0, -2, w + 5, 8, 2), NULL, 0, -20.0f, 0, false, 0.8f);

            plan->pylon_count += 2;
            plan->has_vault = true;
            if (plan->stage < FORERUNNER_STAGE_TERRACES) plan->stage = FORERUNNER_STAGE_TERRACES;
            return true;
        }

        case FORERUNNER_RULE_RAISE_GRAVITY_CORE: {
            int h_spire = 17 + ((h >> 16) % 4) * 2;
            int sx = 4 + ((h >> 18) % 2);

            int spire_l = forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(-sx - 2, 0, -13, -sx, h_spire, -11),
                                              NULL, 0, 15.0f, 0, false, 0.9f);
            int spire_r = forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                              make_forerunner_box(sx, 0, -13, sx + 2, h_spire, -11),
                                              NULL, 0, -15.0f, 0, false, 0.9f);

            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(-sx - 2, h_spire - 1, -13, -sx, h_spire + 1, -11),
                                NULL, 0, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(sx, h_spire - 1, -13, sx + 2, h_spire + 1, -11),
                                NULL, 0, 0, 0, false, 0.1f);

            int y_core = 6 + ((h >> 20) % 3) * 2;
            int p_spires[2] = { spire_l, spire_r };
            forerunner_add_node(plan, FORERUNNER_PRIM_GRAVITY_CORE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, y_core, -13, 2, y_core + 4, -9),
                                p_spires, 2, 0, 0, false, 0.4f);
            // Lintel beam spans across both spires to eliminate structural air gap
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-sx - 1, y_core + 1, -13, sx + 1, y_core + 3, -9),
                                p_spires, 2, 0, 0, false, 0.6f);

            plan->has_core = true;
            if (plan->stage < FORERUNNER_STAGE_APEX) plan->stage = FORERUNNER_STAGE_APEX;
            return true;
        }

        case FORERUNNER_RULE_CARVE_LIGHT_CHANNELS: {
            int p_inner = w - 1;
            int p_outer = p_inner + 4 + ((h >> 12) % 2);

            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-w, 1, -8, -w + 1, 2, 8), NULL, 0, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(w - 1, 1, -8, w, 2, 8), NULL, 0, 0, 0, false, 0.1f);

            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_BRONZE_ACCENT,
                                make_forerunner_box(-p_outer - 1, 0, 13, -p_outer, 3, 15), NULL, 0, 0, 0, false, 0.2f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_BRONZE_ACCENT,
                                make_forerunner_box(p_outer, 0, 13, p_outer + 1, 3, 15), NULL, 0, 0, 0, false, 0.2f);

            plan->light_channel_count += 4;
            plan->stage = FORERUNNER_STAGE_APEX;
            return true;
        }

        // ===================================================================
        // 2. Crossroads Rules (Cruciform Abyss & Bi-Level Bridges)
        // ===================================================================
        case FORERUNNER_RULE_EXCAVATE_CRUCIFORM_CHASM: {
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d, -l, w, 0, l), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-w, -d - 1, -l, w, -d, l), NULL, 0, 0, 0, false, 1.0f);

            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d, -w, l, 0, w), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d - 1, -w, l, -d, w), NULL, 0, 0, 0, false, 1.0f);

            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-2, -d, -2, 2, -d + 3, 2), NULL, 0, 0, 0, false, 1.0f);

            plan->has_chasm = true;
            if (plan->stage < FORERUNNER_STAGE_FOUNDATION) plan->stage = FORERUNNER_STAGE_FOUNDATION;
            return true;
        }

        case FORERUNNER_RULE_ERECT_CORNER_PYLONS: {
            float cant = 18.0f + ((h >> 8) % 3) * 4.0f;
            int hp = 14 + ((h >> 10) % 3) * 2;
            int c1 = w + 1;
            int c2 = c1 + 5;

            int p_sw = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                           make_forerunner_box(-c2, 0, c1, -c1, hp, c2), NULL, 0, cant, -cant, false, 0.8f);
            int p_se = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                           make_forerunner_box(c1, 0, c1, c2, hp, c2), NULL, 0, -cant, -cant, false, 0.8f);
            int p_nw = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                           make_forerunner_box(-c2, 0, -c2, -c1, hp, -c1), NULL, 0, cant, cant, false, 0.8f);
            int p_ne = forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                           make_forerunner_box(c1, 0, -c2, c2, hp, -c1), NULL, 0, -cant, cant, false, 0.8f);

            int p_s[2] = { p_sw, p_se };
            int p_n[2] = { p_nw, p_ne };
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-c2, hp - 2, c1 - 4, c2, hp, c2), p_s, 2, 0, 0, false, 0.6f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-c2, hp - 2, -c2, c2, hp, -c1 + 4), p_n, 2, 0, 0, false, 0.6f);

            plan->pylon_count += 4;
            plan->has_gateway = true;
            if (plan->stage < FORERUNNER_STAGE_SUPERSTRUCTURE) plan->stage = FORERUNNER_STAGE_SUPERSTRUCTURE;
            return true;
        }

        case FORERUNNER_RULE_SPAN_CROSS_BRIDGES: {
            // Lower Bridge (North-South) at Y = 2
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 2, -l + 4, 2, 3, l - 4), NULL, 0, 0, 0, false, 0.4f);
            forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 2, l - 4, 2, 3, l), NULL, 0, 0, 0, false, 0.7f);
            forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 2, -l, 2, 3, -l + 4), NULL, 0, 0, 0, false, 0.7f);

            // Upper Bridge (East-West) at Y = 6 crossing OVER lower bridge
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-l + 4, 6, -2, l - 4, 7, 2), NULL, 0, 0, 0, false, 0.4f);
            forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(l - 4, 6, -2, l, 7, 2), NULL, 0, 0, 0, false, 0.7f);
            forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-l, 6, -2, -l + 4, 7, 2), NULL, 0, 0, 0, false, 0.7f);

            // Support piers anchoring cross bridges to bedrock
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, -2, 2, 6, 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, -l, 2, 2, -l + 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, l - 2, 2, 2, l), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-l, 0, -2, -l + 2, 6, 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(l - 2, 0, -2, l, 6, 2), NULL, 0, 0, 0, false, 0.8f);

            plan->has_bridge = true;
            if (plan->stage < FORERUNNER_STAGE_TRANSIT) plan->stage = FORERUNNER_STAGE_TRANSIT;
            return true;
        }

        case FORERUNNER_RULE_ADD_INTERSECTION_NEXUS: {
            // Central support pedestal under core
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 6, -2, 2, 9, 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_GRAVITY_CORE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 9, -2, 2, 13, 2), NULL, 0, 0, 0, false, 0.5f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-4, 10, -4, 4, 12, 4), NULL, 0, 0, 0, false, 0.4f);

            // 4 Corner vertical energy telemetry needles
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-w - 2, 0, w, -w, 18, w + 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(w, 0, w, w + 2, 18, w + 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-w - 2, 0, -w - 2, -w, 18, -w), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(w, 0, -w - 2, w + 2, 18, -w), NULL, 0, 0, 0, false, 0.8f);

            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(-w - 2, 17, w, -w, 19, w + 2), NULL, 0, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(w, 17, w, w + 2, 19, w + 2), NULL, 0, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(-w - 2, 17, -w - 2, -w, 19, -w), NULL, 0, 0, 0, false, 0.1f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(w, 17, -w - 2, w + 2, 19, -w), NULL, 0, 0, 0, false, 0.1f);

            plan->has_core = true;
            if (plan->stage < FORERUNNER_STAGE_APEX) plan->stage = FORERUNNER_STAGE_APEX;
            return true;
        }

        // ===================================================================
        // 3. Crucible Rules (Radial Pylon Ring & Gravity Pit)
        // ===================================================================
        case FORERUNNER_RULE_EXCAVATE_GRAVITY_PIT: {
            int r = w + 3;
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-r, -d, -r, r, 0, r), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-r, -d - 1, -r, r, -d, r), NULL, 0, 0, 0, false, 1.0f);

            plan->has_chasm = true;
            if (plan->stage < FORERUNNER_STAGE_FOUNDATION) plan->stage = FORERUNNER_STAGE_FOUNDATION;
            return true;
        }

        case FORERUNNER_RULE_ERECT_RADIAL_PYLON_RING: {
            int r = w + 2;
            int hp = 15 + ((h >> 10) % 3) * 2;
            float cant = 22.0f + ((h >> 8) % 3) * 4.0f;

            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-3, 0, r, 3, hp, r + 5), NULL, 0, 0, -cant, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-3, 0, -r - 5, 3, hp, -r), NULL, 0, 0, cant, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-r - 5, 0, -3, -r, hp, 3), NULL, 0, cant, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(r, 0, -3, r + 5, hp, 3), NULL, 0, -cant, 0, false, 0.8f);

            plan->pylon_count += 4;
            plan->has_gateway = true;
            if (plan->stage < FORERUNNER_STAGE_SUPERSTRUCTURE) plan->stage = FORERUNNER_STAGE_SUPERSTRUCTURE;
            return true;
        }

        case FORERUNNER_RULE_SPAN_RADIAL_CATWALK: {
            int r = w + 2;
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 3, 3, 2, 4, r), NULL, 0, 0, 0, false, 0.3f);
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 3, -r, 2, 4, -3), NULL, 0, 0, 0, false, 0.3f);
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-r, 3, -2, -3, 4, 2), NULL, 0, 0, 0, false, 0.3f);
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(3, 3, -2, r, 4, 2), NULL, 0, 0, 0, false, 0.3f);

            // Central support pillar grounding catwalk hub to bedrock
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, -2, 2, 3, 2), NULL, 0, 0, 0, false, 0.8f);

            forerunner_add_node(plan, FORERUNNER_PRIM_BRIDGE_SPAN, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-3, 3, -3, 3, 4, 3), NULL, 0, 0, 0, false, 0.6f);

            plan->has_bridge = true;
            if (plan->stage < FORERUNNER_STAGE_TRANSIT) plan->stage = FORERUNNER_STAGE_TRANSIT;
            return true;
        }

        case FORERUNNER_RULE_SUSPEND_CENTRIFUGE_CORE: {
            // Support pedestal resting on catwalk hub
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 4, -2, 2, 7, 2), NULL, 0, 0, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_GRAVITY_CORE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-2, 7, -2, 2, 11, 2), NULL, 0, 0, 0, false, 0.4f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-4, 8, -4, 4, 10, 4), NULL, 0, 0, 0, false, 0.3f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-6, 9, -6, 6, 11, 6), NULL, 0, 0, 0, false, 0.3f);

            plan->has_core = true;
            if (plan->stage < FORERUNNER_STAGE_APEX) plan->stage = FORERUNNER_STAGE_APEX;
            return true;
        }

        // ===================================================================
        // 4. Spire Citadel Rules (Towering Pinnacle & Flying Buttresses)
        // ===================================================================
        case FORERUNNER_RULE_EXCAVATE_ANNULAR_MOAT: {
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d, -l, l, 0, -w - 3), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d, w + 3, l, 0, l), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d, -w - 3, -w - 3, 0, w + 3), NULL, 0, 0, 0, true, 0.0f);
            forerunner_add_node(plan, FORERUNNER_PRIM_CHASM, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(w + 3, -d, -w - 3, l, 0, w + 3), NULL, 0, 0, 0, true, 0.0f);

            forerunner_add_node(plan, FORERUNNER_PRIM_TERRACE, FORERUNNER_MAT_DARK_BASALT,
                                make_forerunner_box(-l, -d - 1, -l, l, -d, l), NULL, 0, 0, 0, false, 1.0f);

            plan->has_chasm = true;
            if (plan->stage < FORERUNNER_STAGE_FOUNDATION) plan->stage = FORERUNNER_STAGE_FOUNDATION;
            return true;
        }

        case FORERUNNER_RULE_ERECT_CENTRAL_PINNACLE: {
            int hp = 22 + ((h >> 10) % 2) * 2;
            forerunner_add_node(plan, FORERUNNER_PRIM_SPIRE, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-3, 0, -3, 3, hp, 3), NULL, 0, 0, 0, false, 0.9f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LIGHT_CHANNEL, FORERUNNER_MAT_SOLAR_AMBER,
                                make_forerunner_box(-3, hp - 1, -3, 3, hp + 3, 3), NULL, 0, 0, 0, false, 0.2f);

            plan->pylon_count += 1;
            plan->has_gateway = true;
            if (plan->stage < FORERUNNER_STAGE_SUPERSTRUCTURE) plan->stage = FORERUNNER_STAGE_SUPERSTRUCTURE;
            return true;
        }

        case FORERUNNER_RULE_ADD_FLYING_BUTTRESSES: {
            float cant = 28.0f;
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, 4, 2, 14, 11), NULL, 0, 0, -cant, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-2, 0, -11, 2, 14, -4), NULL, 0, 0, cant, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-11, 0, -2, -4, 14, 2), NULL, 0, cant, 0, false, 0.8f);
            forerunner_add_node(plan, FORERUNNER_PRIM_PYLON, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(4, 0, -2, 11, 14, 2), NULL, 0, -cant, 0, false, 0.8f);

            plan->pylon_count += 4;
            plan->has_bridge = true;
            if (plan->stage < FORERUNNER_STAGE_TRANSIT) plan->stage = FORERUNNER_STAGE_TRANSIT;
            return true;
        }

        case FORERUNNER_RULE_ATTACH_OBSERVATION_DECKS: {
            forerunner_add_node(plan, FORERUNNER_PRIM_HARDLIGHT_BRIDGE, FORERUNNER_MAT_HARDLIGHT_CYAN,
                                make_forerunner_box(-6, 8, -6, 6, 9, 6), NULL, 0, 0, 0, false, 0.3f);
            forerunner_add_node(plan, FORERUNNER_PRIM_LINTEL, FORERUNNER_MAT_PEWTER,
                                make_forerunner_box(-7, 8, -7, 7, 10, 7), NULL, 0, 0, 0, false, 0.4f);

            plan->has_vault = true;
            if (plan->stage < FORERUNNER_STAGE_TERRACES) plan->stage = FORERUNNER_STAGE_TERRACES;
            return true;
        }

        default:
            return false;
    }
}

// ---------------------------------------------------------------------------
// Developmental Growth Loop
// ---------------------------------------------------------------------------
ForerunnerPlan generate_forerunner_structure(uint32_t seed, ForerunnerArchetype archetype, ForerunnerStage target_stage) {
    ForerunnerPlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.archetype = archetype;
    plan.stage = FORERUNNER_STAGE_FOUNDATION;

    // Fast hash to derive high-entropy architectural dimensions
    uint32_t h = (seed * 2654435761u) ^ (seed >> 16);
    h ^= (h << 13); h ^= (h >> 17); h ^= (h << 5);

    int half_w = 5 + (h % 3);                // Half-width: 5, 6, or 7 modules
    int depth  = 6 + ((h >> 3) % 3) * 2;     // Depth: 6, 8, or 10 modules
    int half_l = 14 + ((h >> 6) % 3) * 2;    // Half-length: 14, 16, or 18 modules

    plan.chasm_width  = half_w * 2;
    plan.chasm_depth  = depth;
    plan.chasm_length = half_l * 2;

    uint32_t prng = seed;
    int max_steps = 25;

    while (max_steps-- > 0) {
        if (plan.stage >= target_stage) {
            if (target_stage == FORERUNNER_STAGE_FOUNDATION && plan.has_chasm) break;
            if (target_stage == FORERUNNER_STAGE_SUPERSTRUCTURE && plan.has_gateway) break;
            if (target_stage == FORERUNNER_STAGE_TRANSIT && plan.has_bridge) break;
            if (target_stage == FORERUNNER_STAGE_TERRACES && plan.has_vault) break;
            if (target_stage == FORERUNNER_STAGE_APEX && plan.has_core) break;
        }

        ForerunnerOpportunity opps[MAX_FORERUNNER_OPPORTUNITIES];
        int opp_count = forerunner_coalgebra_inspect(&plan, opps, MAX_FORERUNNER_OPPORTUNITIES);
        if (opp_count == 0) break;

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
