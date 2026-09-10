#include "megalith_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Simple deterministic PRNG
static inline float mega_prng(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return (float)(*state & 0x00FFFFFF) / (float)0x01000000;
}

static inline MegaBox3D make_box3d(int min_x, int min_y, int min_z, int max_x, int max_y, int max_z) {
    return (MegaBox3D){ min_x, min_y, min_z, max_x, max_y, max_z };
}

static int add_stone(MegalithPlan *plan, MegalithPrimitiveType type, MegaBox3D box,
                     const int *parents, int parent_cnt, float roughness) {
    if (plan->node_count >= MAX_MEGALITH_NODES) return -1;
    int id = plan->node_count++;
    MegalithNode *node = &plan->nodes[id];
    node->node_id = id;
    node->type = type;
    node->box = box;
    node->support_count = (parent_cnt > 4) ? 4 : parent_cnt;
    for (int i = 0; i < node->support_count; ++i) {
        node->support_parents[i] = parents[i];
    }
    // Subtle organic tilt
    uint32_t seed = plan->seed + (uint32_t)id * 7919u;
    node->tilt_x = (mega_prng(&seed) - 0.5f) * 0.15f;
    node->tilt_z = (mega_prng(&seed) - 0.5f) * 0.15f;
    node->roughness = roughness;
    node->flags = 0;
    return id;
}

// ---------------------------------------------------------------------------
// Coalgebra gamma(S): Inspects support, void, boundary, axis, and terrain
// ---------------------------------------------------------------------------
int megalith_coalgebra_inspect(const MegalithPlan *plan, MegalithOpportunity *out_opps, int max_opps) {
    int count = 0;

    if (plan->archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE) {
        // Stage 0 -> 1: Menhir to Dolmen / Portal Tomb
        if (!plan->has_capstone && count < max_opps) {
            // Find existing orthostats
            int orthos[4];
            int ortho_cnt = 0;
            for (int i = 0; i < plan->node_count && ortho_cnt < 4; ++i) {
                if (plan->nodes[i].type == MEGALITH_PRIMITIVE_ORTHOSTAT) {
                    orthos[ortho_cnt++] = i;
                }
            }

            if (ortho_cnt < 2) {
                // gamma_void: Exposes pairing opportunity to erect opposing orthostat
                out_opps[count++] = (MegalithOpportunity){
                    .rule = MEGALITH_RULE_PAIR_ORTHOSTATS,
                    .target_box = make_box3d(2, 0, 0, 3, 5, 2),
                    .support_count = 0,
                    .stability_score = 1.0f,
                    .alignment_score = 0.8f,
                    .enclosure_score = 0.5f,
                    .total_score = 30.0f
                };
            } else {
                // gamma_support: Orthostats exist. Exposes capstone bridging opportunity
                int min_x = 999, max_x = -999, max_y = 0, min_z = 999, max_z = -999;
                for (int i = 0; i < ortho_cnt; ++i) {
                    const MegaBox3D *b = &plan->nodes[orthos[i]].box;
                    if (b->min_x < min_x) min_x = b->min_x;
                    if (b->max_x > max_x) max_x = b->max_x;
                    if (b->max_y > max_y) max_y = b->max_y;
                    if (b->min_z < min_z) min_z = b->min_z;
                    if (b->max_z > max_z) max_z = b->max_z;
                }
                // Capstone spans across all supports with slight overhang
                MegaBox3D cap_box = make_box3d(min_x - 1, max_y, min_z - 1, max_x + 1, max_y + 2, max_z + 1);
                MegalithOpportunity opp = {
                    .rule = MEGALITH_RULE_BRIDGE_CAPSTONE,
                    .target_box = cap_box,
                    .support_count = ortho_cnt,
                    .stability_score = 1.0f,
                    .alignment_score = 1.0f,
                    .enclosure_score = 0.8f,
                    .total_score = 40.0f
                };
                for (int i = 0; i < ortho_cnt; ++i) opp.support_parents[i] = orthos[i];
                out_opps[count++] = opp;
            }
        }

        // Stage 2: Chamber Enclosure (Polygonal radial vault) - once portal dolmen is established
        if (plan->has_capstone && !plan->has_corbel && count < max_opps) {
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_ENCLOSE_CHAMBER,
                .target_box = make_box3d(-4, 0, -4, 4, 6, 4),
                .support_count = 0,
                .stability_score = 0.9f,
                .alignment_score = 0.9f,
                .enclosure_score = 1.0f,
                .total_score = 45.0f
            };
        }

        // Stage 3: Solstice Passage Extension - once chamber is vaulted
        if (plan->has_corbel && plan->passage_length < 8 && count < max_opps) {
            int p_z = plan->chamber_cz + 4 + plan->passage_length;
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_EXTEND_PASSAGE,
                .target_box = make_box3d(-2, 0, p_z, 2, 4, p_z + 3),
                .support_count = 0,
                .stability_score = 0.9f,
                .alignment_score = 1.0f, // Aligned along sunrise axis
                .enclosure_score = 0.8f,
                .total_score = 40.0f
            };
        }

        // Stage 4: Earthen Tumulus Barrow & Peristalith Kerb Ring - once passage is complete
        if (plan->has_passage && !plan->has_mound && count < max_opps) {
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_ACCUMULATE_MOUND,
                .target_box = make_box3d(-14, 0, -12, 14, 8, 16),
                .support_count = 0,
                .stability_score = 1.0f,
                .alignment_score = 0.7f,
                .enclosure_score = 1.0f,
                .total_score = 50.0f
            };
        }
    } else {
        // Archetype B: Megalithic Stone Circle / Henge
        // gamma_boundary: Check circle completion
        if (plan->circle_stones_placed < plan->circle_stones_total && count < max_opps) {
            float angle = (2.0f * 3.14159265f * (float)plan->circle_stones_placed) / (float)plan->circle_stones_total;
            int r = plan->circle_radius;
            int sx = (int)roundf(plan->chamber_cx + r * cosf(angle));
            int sz = (int)roundf(plan->chamber_cz + r * sinf(angle));
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_EXPAND_CIRCLE,
                .target_box = make_box3d(sx - 1, 0, sz - 1, sx + 1, 5, sz + 1),
                .support_count = 0,
                .stability_score = 0.95f,
                .alignment_score = 0.8f,
                .enclosure_score = 0.9f,
                .total_score = 30.0f
            };
        }

        // Inner Monumental Trilithon Horseshoe - once outer ring is complete
        if (plan->circle_stones_placed >= plan->circle_stones_total && !plan->has_corbel && count < max_opps) {
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_RAISE_TRILITHON,
                .target_box = make_box3d(-5, 0, -5, 5, 7, 5),
                .support_count = 0,
                .stability_score = 1.0f,
                .alignment_score = 1.0f,
                .enclosure_score = 0.85f,
                .total_score = 45.0f
            };
        }

        // Processional Avenue - once inner trilithons are raised
        if (plan->has_corbel && !plan->has_avenue && count < max_opps) {
            out_opps[count++] = (MegalithOpportunity){
                .rule = MEGALITH_RULE_ALIGN_AVENUE,
                .target_box = make_box3d(-3, 0, plan->circle_radius, 3, 4, plan->circle_radius + 14),
                .support_count = 0,
                .stability_score = 0.9f,
                .alignment_score = 1.0f, // Solstice vector
                .enclosure_score = 0.6f,
                .total_score = 35.0f
            };
        }
    }

    return count;
}

// ---------------------------------------------------------------------------
// Algebra alpha(S, o): Synthesizes parts into composite megalithic state
// ---------------------------------------------------------------------------
bool megalith_algebra_apply(MegalithPlan *plan, const MegalithOpportunity *opp) {
    if (!plan || !opp) return false;

    switch (opp->rule) {
        case MEGALITH_RULE_ERECT_ORTHOSTAT: {
            add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, opp->target_box, NULL, 0, 0.4f);
            return true;
        }

        case MEGALITH_RULE_PAIR_ORTHOSTATS: {
            // Erect opposing stone to create entrance portal
            add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, opp->target_box, NULL, 0, 0.4f);
            return true;
        }

        case MEGALITH_RULE_BRIDGE_CAPSTONE: {
            // Place massive horizontal capstone resting on supports
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE, opp->target_box,
                      opp->support_parents, opp->support_count, 0.5f);
            plan->has_capstone = true;
            if (plan->stage < MEGALITH_STAGE_DOLMEN) plan->stage = MEGALITH_STAGE_DOLMEN;
            return true;
        }

        case MEGALITH_RULE_ENCLOSE_CHAMBER: {
            // Ring of inward-leaning orthostats forming burial vault
            int r = plan->chamber_radius;
            int h = 5;
            int num_wall_stones = 8;
            int wall_ids[8];
            for (int i = 0; i < num_wall_stones; ++i) {
                // Leave entrance gap along +Z
                if (i == 2) continue; // Entrance opening
                float angle = (2.0f * 3.14159265f * (float)i) / (float)num_wall_stones;
                int ox = (int)roundf(plan->chamber_cx + r * cosf(angle));
                int oz = (int)roundf(plan->chamber_cz + r * sinf(angle));
                wall_ids[i] = add_stone(plan, MEGALITH_PRIMITIVE_WALL_SLAB,
                                        make_box3d(ox - 1, 0, oz - 1, ox + 1, h, oz + 1),
                                        NULL, 0, 0.35f);
            }

            // Central ritual hearth / offering altar stone
            add_stone(plan, MEGALITH_PRIMITIVE_HEARTH,
                      make_box3d(plan->chamber_cx - 1, 0, plan->chamber_cz - 1,
                                 plan->chamber_cx + 1, 1, plan->chamber_cz + 1),
                      NULL, 0, 0.2f);
            plan->has_hearth = true;

            // Inward-stepping corbel roof slabs
            int corbel_parents[4] = { wall_ids[0], wall_ids[1], wall_ids[4], wall_ids[5] };
            add_stone(plan, MEGALITH_PRIMITIVE_CORBEL,
                      make_box3d(plan->chamber_cx - r + 1, h - 1, plan->chamber_cz - r + 1,
                                 plan->chamber_cx + r - 1, h + 1, plan->chamber_cz + r - 1),
                      corbel_parents, 4, 0.4f);

            // Great central capstone crowning the dome
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE,
                      make_box3d(plan->chamber_cx - 2, h + 1, plan->chamber_cz - 2,
                                 plan->chamber_cx + 2, h + 2, plan->chamber_cz + 2),
                      corbel_parents, 2, 0.6f);

            plan->has_corbel = true;
            plan->has_capstone = true;
            if (plan->stage < MEGALITH_STAGE_CHAMBER) plan->stage = MEGALITH_STAGE_CHAMBER;
            return true;
        }

        case MEGALITH_RULE_EXTEND_PASSAGE: {
            int p_z = plan->chamber_cz + 4 + plan->passage_length;
            int w = 2;
            int h = 4;
            // Left orthostat wall slab
            int left_id = add_stone(plan, MEGALITH_PRIMITIVE_WALL_SLAB,
                                    make_box3d(-w, 0, p_z, -w + 1, h, p_z + 3),
                                    NULL, 0, 0.3f);
            // Right orthostat wall slab
            int right_id = add_stone(plan, MEGALITH_PRIMITIVE_WALL_SLAB,
                                     make_box3d(w - 1, 0, p_z, w, h, p_z + 3),
                                     NULL, 0, 0.3f);
            // Paved flagstone passage floor
            add_stone(plan, MEGALITH_PRIMITIVE_FLOOR_SLAB,
                      make_box3d(-w + 1, 0, p_z, w - 1, 1, p_z + 3),
                      NULL, 0, 0.1f);

            // Transverse lintel / roof slab bridging the passage
            int parents[2] = { left_id, right_id };
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE,
                      make_box3d(-w, h, p_z, w, h + 1, p_z + 3),
                      parents, 2, 0.4f);

            plan->passage_length += 3;
            plan->has_passage = true;
            if (plan->stage < MEGALITH_STAGE_PASSAGE_GRAVE) plan->stage = MEGALITH_STAGE_PASSAGE_GRAVE;
            return true;
        }

        case MEGALITH_RULE_ACCUMULATE_MOUND: {
            // Earthen tumulus barrow covering the burial vault
            add_stone(plan, MEGALITH_PRIMITIVE_MOUND, opp->target_box, NULL, 0, 0.1f);

            // Peristalith: Ring of massive contiguous kerb stones delineating mound base
            int r_kerb = (opp->target_box.max_x - opp->target_box.min_x) / 2;
            int num_kerbs = 20;
            for (int k = 0; k < num_kerbs; ++k) {
                float angle = (2.0f * 3.14159265f * (float)k) / (float)num_kerbs;
                // Skip passage entrance gap at +Z
                int kz = (int)roundf(r_kerb * sinf(angle));
                if (kz > r_kerb - 3) continue; // Entrance portal threshold
                int kx = (int)roundf(r_kerb * cosf(angle));
                add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT,
                          make_box3d(kx - 1, 0, kz - 1, kx + 1, 2, kz + 1),
                          NULL, 0, 0.5f);
            }

            plan->has_mound = true;
            plan->stage = MEGALITH_STAGE_TUMULUS;
            return true;
        }

        case MEGALITH_RULE_EXPAND_CIRCLE: {
            // Add standing stone along circular arc
            add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, opp->target_box, NULL, 0, 0.45f);
            plan->circle_stones_placed++;
            return true;
        }

        case MEGALITH_RULE_RAISE_TRILITHON: {
            // Inner monumental horseshoe of 3 great trilithons (like Stonehenge)
            // Trilithon 1: Center rear (Winter Solstice Sunset axis)
            int left1  = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(-2, 0, -5, -1, 6, -4), NULL, 0, 0.4f);
            int right1 = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(1, 0, -5, 2, 6, -4), NULL, 0, 0.4f);
            int p1[2] = { left1, right1 };
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE, make_box3d(-3, 6, -5, 3, 7, -4), p1, 2, 0.5f);

            // Trilithon 2: Left flank
            int left2  = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(-5, 0, -1, -4, 6, 0), NULL, 0, 0.4f);
            int right2 = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(-5, 0, 2, -4, 6, 3), NULL, 0, 0.4f);
            int p2[2] = { left2, right2 };
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE, make_box3d(-5, 6, -1, -4, 7, 3), p2, 2, 0.5f);

            // Trilithon 3: Right flank
            int left3  = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(4, 0, -1, 5, 6, 0), NULL, 0, 0.4f);
            int right3 = add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(4, 0, 2, 5, 6, 3), NULL, 0, 0.4f);
            int p3[2] = { left3, right3 };
            add_stone(plan, MEGALITH_PRIMITIVE_CAPSTONE, make_box3d(4, 6, -1, 5, 7, 3), p3, 2, 0.5f);

            // Central Altar Stone / Offering slab
            add_stone(plan, MEGALITH_PRIMITIVE_HEARTH, make_box3d(-1, 0, -1, 1, 1, 1), NULL, 0, 0.2f);

            plan->has_corbel = true;
            plan->has_hearth = true;
            plan->stage = MEGALITH_STAGE_HENGE;
            return true;
        }

        case MEGALITH_RULE_ALIGN_AVENUE: {
            // Parallel megalithic rows guiding the approach along the solstice axis
            int z_start = plan->circle_radius + 1;
            int z_end   = z_start + 12;
            int avenue_w = 4;
            for (int z = z_start; z <= z_end; z += 3) {
                add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT,
                          make_box3d(-avenue_w, 0, z, -avenue_w + 1, 4, z + 1),
                          NULL, 0, 0.4f);
                add_stone(plan, MEGALITH_PRIMITIVE_ORTHOSTAT,
                          make_box3d(avenue_w - 1, 0, z, avenue_w, 4, z + 1),
                          NULL, 0, 0.4f);
            }
            plan->has_avenue = true;
            return true;
        }

        default:
            return false;
    }
}

// ---------------------------------------------------------------------------
// Prehistoric Developmental Growth Loop
// ---------------------------------------------------------------------------
MegalithPlan generate_megalith_structure(uint32_t seed, MegalithArchetype archetype, MegalithStage target_stage) {
    MegalithPlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.archetype = archetype;
    plan.stage = MEGALITH_STAGE_MENHIR;

    // Astronomical axis vector (Winter Solstice Sunrise ~ +Z)
    plan.axis_dx = 0.0f;
    plan.axis_dz = 1.0f;

    plan.chamber_cx = 0;
    plan.chamber_cz = 0;
    plan.chamber_radius = 4;
    plan.passage_length = 0;

    plan.circle_radius = 12;
    plan.circle_stones_placed = 0;
    plan.circle_stones_total = 16;

    uint32_t prng = seed;

    // 1. Initial Seed: Single massive standing stone (The Primordial Menhir)
    add_stone(&plan, MEGALITH_PRIMITIVE_ORTHOSTAT, make_box3d(-2, 0, 0, -1, 5, 2), NULL, 0, 0.4f);

    // 2. Developmental stepping loop: Weight and void decide which forms return
    int max_steps = 35;
    while (max_steps-- > 0) {
        if (plan.stage >= target_stage) {
            if (archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE && target_stage == MEGALITH_STAGE_TUMULUS && plan.has_mound) break;
            if (archetype == MEGALITH_ARCHETYPE_STONE_CIRCLE && target_stage == MEGALITH_STAGE_HENGE && plan.has_avenue) break;
            if (archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE && target_stage == MEGALITH_STAGE_DOLMEN && plan.has_capstone) break;
            if (archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE && target_stage == MEGALITH_STAGE_CHAMBER && plan.has_corbel) break;
            if (archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE && target_stage == MEGALITH_STAGE_PASSAGE_GRAVE && plan.passage_length >= 6) break;
            if (archetype == MEGALITH_ARCHETYPE_STONE_CIRCLE && target_stage <= MEGALITH_STAGE_DOLMEN && plan.circle_stones_placed >= 8) break;
        }

        MegalithOpportunity opps[MAX_MEGALITH_OPPORTUNITIES];
        int opp_count = megalith_coalgebra_inspect(&plan, opps, MAX_MEGALITH_OPPORTUNITIES);
        if (opp_count == 0) break;

        // Weighted probabilistic selection
        float total_weight = 0.0f;
        for (int i = 0; i < opp_count; ++i) {
            total_weight += opps[i].total_score;
        }

        float r = mega_prng(&prng) * total_weight;
        float accum = 0.0f;
        int chosen_idx = 0;
        for (int i = 0; i < opp_count; ++i) {
            accum += opps[i].total_score;
            if (r <= accum) {
                chosen_idx = i;
                break;
            }
        }

        megalith_algebra_apply(&plan, &opps[chosen_idx]);
    }

    return plan;
}

// ---------------------------------------------------------------------------
// ASCII Blueprint Diagnostic Renderer for Megaliths
// ---------------------------------------------------------------------------
void print_megalith_blueprint(const MegalithPlan *plan) {
    if (!plan) return;

    int min_x = -16, max_x = 16;
    int min_z = -16, max_z = 24;
    int w = max_x - min_x + 1;
    int h = max_z - min_z + 1;

    char *grid = (char *)malloc(w * h);
    memset(grid, ' ', w * h);

    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        char sym = '?';
        switch (n->type) {
            case MEGALITH_PRIMITIVE_ORTHOSTAT:  sym = 'O'; break; // Standing stone
            case MEGALITH_PRIMITIVE_CAPSTONE:   sym = '='; break; // Capstone / Lintel
            case MEGALITH_PRIMITIVE_WALL_SLAB:  sym = '#'; break; // Chamber / passage wall
            case MEGALITH_PRIMITIVE_CORBEL:     sym = '^'; break; // Corbel vault
            case MEGALITH_PRIMITIVE_MOUND:      sym = '~'; break; // Earth tumulus
            case MEGALITH_PRIMITIVE_FLOOR_SLAB: sym = '.'; break; // Paved floor
            case MEGALITH_PRIMITIVE_HEARTH:     sym = '@'; break; // Ritual hearth
            default: continue;
        }

        for (int z = n->box.min_z; z < n->box.max_z; ++z) {
            for (int x = n->box.min_x; x < n->box.max_x; ++x) {
                int gx = x - min_x;
                int gz = z - min_z;
                if (gx >= 0 && gx < w && gz >= 0 && gz < h) {
                    char prev = grid[gz * w + gx];
                    // Precedence: Hearth > Capstone > Orthostat > Wall > Corbel > Mound
                    if (prev == '@' || prev == '=') continue;
                    if (prev == 'O' && sym != '@' && sym != '=') continue;
                    grid[gz * w + gx] = sym;
                }
            }
        }
    }

    const char *arch_name = (plan->archetype == MEGALITH_ARCHETYPE_PASSAGE_GRAVE) ?
                            "Passage Grave & Barrow" : "Stone Circle & Henge";

    printf("\n=======================================================\n");
    printf(" MEGALITH BLUEPRINT: %s (Nodes: %d, Seed: %u)\n", arch_name, plan->node_count, plan->seed);
    printf(" Key: [O] Orthostat / Standing Stone, [=] Capstone, [#] Wall Slab,\n");
    printf("      [^] Corbel Vault, [@] Hearth / Altar, [.] Floor, [~] Tumulus Mound\n");
    printf("=======================================================\n");

    for (int z = h - 1; z >= 0; z -= 2) {
        printf("%3d | ", z + min_z);
        for (int x = 0; x < w; ++x) {
            putchar(grid[z * w + x]);
        }
        putchar('\n');
    }
    printf("      ");
    for (int x = 0; x < w; ++x) putchar('-');
    printf("\n      ");
    for (int x = 0; x < w; x += 4) printf("%-4d", x + min_x);
    printf("\n\n");

    free(grid);
}
