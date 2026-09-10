#include "greek_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Simple deterministic PRNG for grammar selection
static uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    if (x == 0) x = 0x12345678;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

static float prng_uniform(uint32_t *state) {
    return (float)(xorshift32(state) & 0xFFFFFF) / (float)0xFFFFFF;
}

static ModBox3D make_box(int min_x, int min_y, int min_z, int max_x, int max_y, int max_z) {
    ModBox3D b = { min_x, min_y, min_z, max_x, max_y, max_z };
    return b;
}

static ModBox3D compute_total_footprint(const TemplePlan *plan) {
    if (plan->node_count == 0) return make_box(0, 0, 0, 0, 0, 0);
    ModBox3D res = plan->nodes[0].box;
    for (int i = 1; i < plan->node_count; ++i) {
        if (plan->nodes[i].type == ARCH_PRIMITIVE_STEP) continue; // ignore basement steps
        const ModBox3D *b = &plan->nodes[i].box;
        if (b->min_x < res.min_x) res.min_x = b->min_x;
        if (b->max_x > res.max_x) res.max_x = b->max_x;
        if (b->min_y < res.min_y) res.min_y = b->min_y;
        if (b->max_y > res.max_y) res.max_y = b->max_y;
        if (b->min_z < res.min_z) res.min_z = b->min_z;
        if (b->max_z > res.max_z) res.max_z = b->max_z;
    }
    return res;
}

static void add_node(TemplePlan *plan, ArchPrimitiveType type, ModBox3D box, int col_count, int col_spacing, bool exterior) {
    if (plan->node_count >= MAX_ARCH_NODES) return;
    ArchNode *n = &plan->nodes[plan->node_count++];
    n->type = type;
    n->box = box;
    n->column_count = col_count;
    n->column_spacing = col_spacing;
    n->is_exterior = exterior;
    n->flags = 0;
}

// Initial seed S_0: Construct sacred rectangular cella
static void make_initial_cella(TemplePlan *plan) {
    plan->node_count = 0;
    plan->site_count = 0;
    plan->complexity = 1;
    plan->stage = TEMPLE_STAGE_CELLA;
    plan->has_pediment = false;
    plan->has_stylobate = false;

    int half_w = plan->cella_width_m / 2;
    int half_l = plan->cella_length_m / 2;
    int h = plan->cella_height_m;

    // 1. Cella interior chamber volume
    add_node(plan, ARCH_PRIMITIVE_CELL, make_box(-half_w + 1, 0, -half_l + 1, half_w - 1, h, half_l - 1), 0, 0, false);

    // 2. Cella Walls (Masonry with 1-module thickness)
    // Left wall (-X)
    add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-half_w, 0, -half_l, -half_w + 1, h, half_l), 0, 0, true);
    // Right wall (+X)
    add_node(plan, ARCH_PRIMITIVE_WALL, make_box(half_w - 1, 0, -half_l, half_w, h, half_l), 0, 0, true);
    // Back wall (-Z)
    add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-half_w, 0, -half_l, half_w, h, -half_l + 1), 0, 0, true);
    // Front wall (+Z) - split into left jamb, right jamb, leaving portal in center
    int door_half_w = 1;
    add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-half_w, 0, half_l - 1, -door_half_w, h, half_l), 0, 0, true);
    add_node(plan, ARCH_PRIMITIVE_WALL, make_box(door_half_w, 0, half_l - 1, half_w, h, half_l), 0, 0, true);
    // Lintel above doorway
    int door_height = (h * 2) / 3;
    add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(-door_half_w, door_height, half_l - 1, door_half_w, h, half_l), 0, 0, true);

    // 3. Exposed Boundary Sites for Coalgebra γ
    // Site 0: Front face (+Z)
    plan->sites[0] = (GrowthSite){
        .site_id = 0,
        .face = SITE_FACE_FRONT,
        .boundary = make_box(-half_w, 0, half_l, half_w, h, half_l),
        .is_on_symmetry_axis = true,
        .symmetry_partner_id = -1,
        .active = true
    };
    // Site 1: Back face (-Z)
    plan->sites[1] = (GrowthSite){
        .site_id = 1,
        .face = SITE_FACE_BACK,
        .boundary = make_box(-half_w, 0, -half_l, half_w, h, -half_l),
        .is_on_symmetry_axis = true,
        .symmetry_partner_id = -1,
        .active = true
    };
    // Site 2: Left flank (-X)
    plan->sites[2] = (GrowthSite){
        .site_id = 2,
        .face = SITE_FACE_LEFT,
        .boundary = make_box(-half_w, 0, -half_l, -half_w, h, half_l),
        .is_on_symmetry_axis = false,
        .symmetry_partner_id = 3,
        .active = true
    };
    // Site 3: Right flank (+X)
    plan->sites[3] = (GrowthSite){
        .site_id = 3,
        .face = SITE_FACE_RIGHT,
        .boundary = make_box(half_w, 0, -half_l, half_w, h, half_l),
        .is_on_symmetry_axis = false,
        .symmetry_partner_id = 2,
        .active = true
    };
    // Site 4: Top roofline (+Y)
    plan->sites[4] = (GrowthSite){
        .site_id = 4,
        .face = SITE_FACE_TOP,
        .boundary = make_box(-half_w, h, -half_l, half_w, h, half_l),
        .is_on_symmetry_axis = true,
        .symmetry_partner_id = -1,
        .active = true
    };
    plan->site_count = 5;
}

// Coalgebra: Observes S and returns candidate opportunities
int greek_coalgebra_inspect(const TemplePlan *plan, GrowthOpportunity *out_opps, int max_opps) {
    int count = 0;

    // Rule 1: Front porch (Pronaos)
    if (plan->stage == TEMPLE_STAGE_CELLA && plan->sites[0].active && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = 0,
            .partner_site_index = -1,
            .rule = RULE_ATTACH_PRONAOS,
            .score = 25.0f // High priority to crown entrance
        };
    }

    // Rule 2: Mirrored rear porch (Opisthodomos)
    if (plan->stage == TEMPLE_STAGE_PROSTYLE && plan->sites[1].active && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = 1,
            .partner_site_index = -1,
            .rule = RULE_ATTACH_OPISTHODOMOS,
            .score = 20.0f // Natural bilateral continuation
        };
    }

    // Rule 3: Peristyle Colonnade Wrap (Lateral expansion paired across symmetry axis)
    if ((plan->stage == TEMPLE_STAGE_PROSTYLE || plan->stage == TEMPLE_STAGE_AMPHIPROSTYLE) &&
        plan->sites[2].active && plan->sites[3].active && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = 2,
            .partner_site_index = 3, // Left and Right coupled simultaneously
            .rule = RULE_WRAP_PERISTYLE,
            .score = 15.0f
        };
    }

    // Rule 4: Stepped Stylobate Podium (Base)
    if (!plan->has_stylobate && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = -1,
            .partner_site_index = -1,
            .rule = RULE_EXPAND_STYLOBATE,
            .score = (plan->stage >= TEMPLE_STAGE_PROSTYLE) ? 18.0f : 5.0f
        };
    }

    // Rule 5: Entablature and Triangular Pediment (Roof)
    if (!plan->has_pediment && plan->sites[4].active && count < max_opps) {
        // Pediment caps building once columns/walls are present
        float score = (plan->stage >= TEMPLE_STAGE_PROSTYLE) ? 30.0f : 1.0f;
        out_opps[count++] = (GrowthOpportunity){
            .site_index = 4,
            .partner_site_index = -1,
            .rule = RULE_RAISE_PEDIMENT,
            .score = score
        };
    }
    // Rule 6: Peristyle Courtyard Forecourt (Temenos)
    if ((plan->stage == TEMPLE_STAGE_PERIPTERAL || plan->stage == TEMPLE_STAGE_AMPHIPROSTYLE) &&
        !plan->has_courtyard && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = 0,
            .partner_site_index = -1,
            .rule = RULE_ATTACH_COURTYARD,
            .score = 22.0f
        };
    }

    // Rule 7: Sacred Reflection Pool (Sunken basin inside courtyard)
    if (plan->has_courtyard && !plan->has_pool && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = -1,
            .partner_site_index = -1,
            .rule = RULE_EXCAVATE_POOL,
            .score = 25.0f
        };
    }

    // Rule 8: Sacred Olive Grove (Symmetrically paired trees)
    if (plan->has_courtyard && !plan->has_trees && count < max_opps) {
        out_opps[count++] = (GrowthOpportunity){
            .site_index = -1,
            .partner_site_index = -1,
            .rule = RULE_PLANT_GROVE,
            .score = 20.0f
        };
    }

    // Garden Precinct Rules (Sanctuary Stage)
    if (plan->stage == TEMPLE_STAGE_SANCTUARY && plan->has_courtyard) {
        if (!plan->has_garden_paths && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_LAYOUT_GARDEN_PATHS,
                .score = 30.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_tholos && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_ERECT_THOLOS,
                .score = 28.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_exedrae && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_INSTALL_EXEDRAE,
                .score = 26.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_naiskoi && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_BUILD_NAISKOI,
                .score = 25.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_fountains && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_CONSTRUCT_FOUNTAINS,
                .score = 27.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_pergolas && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_BUILD_PERGOLAS,
                .score = 24.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_braziers && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_ERECT_BRAZIERS,
                .score = 23.0f
            };
        }
        if (plan->has_garden_paths && !plan->has_garden_flora && count < max_opps) {
            out_opps[count++] = (GrowthOpportunity){
                .site_index = -1,
                .partner_site_index = -1,
                .rule = RULE_POPULATE_GARDEN_FLORA,
                .score = 22.0f
            };
        }
    }

    return count;
}

// Algebra: Applies winning rule to assemble components
bool greek_algebra_apply(TemplePlan *plan, const GrowthOpportunity *opp) {
    if (!opp) return false;

    int half_w = plan->cella_width_m / 2;
    int half_l = plan->cella_length_m / 2;
    int h = plan->column_height_m;
    int spacing = plan->column_spacing_m;

    switch (opp->rule) {
        case RULE_ATTACH_PRONAOS: {
            // Front porch extends +Z by spacing
            int z0 = half_l;
            int z1 = half_l + spacing;

            // Side wall antae projections
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-half_w, 0, z0, -half_w + 1, h, z1), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(half_w - 1, 0, z0, half_w, h, z1), 0, 0, true);

            // Columns in-antis across front line Z=z1
            // 2 symmetrical columns flanked by antae
            int col_offset = spacing;
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(-col_offset, 0, z1 - 1, -col_offset + 1, h, z1), 1, 0, true);
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(col_offset - 1, 0, z1 - 1, col_offset, h, z1), 1, 0, true);

            // Lintel across front columns and antae
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(-half_w, h, z1 - 1, half_w, h + 1, z1), 0, 0, true);

            // Update boundary site
            plan->sites[0].boundary.min_z = z1;
            plan->sites[0].boundary.max_z = z1;
            plan->sites[4].boundary.max_z = z1; // expand roof boundary

            plan->stage = TEMPLE_STAGE_PROSTYLE;
            plan->complexity++;
            return true;
        }

        case RULE_ATTACH_OPISTHODOMOS: {
            // Rear porch extends -Z by spacing (perfect reflection of Pronaos)
            int z0 = -half_l - spacing;
            int z1 = -half_l;

            // Side wall antae projections
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-half_w, 0, z0, -half_w + 1, h, z1), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(half_w - 1, 0, z0, half_w, h, z1), 0, 0, true);

            // 2 columns in-antis across rear line Z=z0
            int col_offset = spacing;
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(-col_offset, 0, z0, -col_offset + 1, h, z0 + 1), 1, 0, true);
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(col_offset - 1, 0, z0, col_offset, h, z0 + 1), 1, 0, true);

            // Lintel across rear columns
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(-half_w, h, z0, half_w, h + 1, z0 + 1), 0, 0, true);

            // Update boundary site
            plan->sites[1].boundary.min_z = z0;
            plan->sites[1].boundary.max_z = z0;
            plan->sites[4].boundary.min_z = z0; // expand roof boundary

            plan->stage = TEMPLE_STAGE_AMPHIPROSTYLE;
            plan->complexity++;
            return true;
        }

        case RULE_WRAP_PERISTYLE: {
            // Colonnade surrounds existing core
            ModBox3D core = compute_total_footprint(plan);
            int span_x = abs(core.min_x) > abs(core.max_x) ? abs(core.min_x) : abs(core.max_x);
            span_x += spacing;
            // Round up to multiple of spacing
            span_x = ((span_x + spacing - 1) / spacing) * spacing;

            int peri_min_z = core.min_z - spacing;
            int peri_max_z = core.max_z + spacing;

            // 1. Front (+Z) and Back (-Z) Colonnade: Symmetrical pairs around X=0
            // Columns are placed at +/- (offset) so an intercolumniation sits on X=0
            for (int col_i = 0; col_i * spacing < span_x; ++col_i) {
                int x_right = col_i * spacing + 1;
                int x_left = -x_right - 1; // [-x_right - 1, -x_right]

                // Front columns (+Z)
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(x_left, 0, peri_max_z - 1, x_left + 1, h, peri_max_z), 1, 0, true);
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(x_right, 0, peri_max_z - 1, x_right + 1, h, peri_max_z), 1, 0, true);

                // Back columns (-Z)
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(x_left, 0, peri_min_z, x_left + 1, h, peri_min_z + 1), 1, 0, true);
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(x_right, 0, peri_min_z, x_right + 1, h, peri_min_z + 1), 1, 0, true);
            }

            // 2. Flank Colonnades: Left (-X) and Right (+X)
            int flank_right = span_x;
            int flank_left = -flank_right - 1;
            for (int z = peri_min_z + spacing; z <= peri_max_z - spacing; z += spacing) {
                // Left flank column
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(flank_left, 0, z, flank_left + 1, h, z + 1), 1, 0, true);
                // Right flank column
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(flank_right, 0, z, flank_right + 1, h, z + 1), 1, 0, true);
            }

            // 3. Continuous Entablature / Architrave Lintel around perimeter
            int lintel_left = flank_left;
            int lintel_right = flank_right + 1; // [-span, span] is self-symmetric
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(lintel_left, h, peri_min_z, lintel_right, h + 1, peri_min_z + 1), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(lintel_left, h, peri_max_z - 1, lintel_right, h + 1, peri_max_z), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(lintel_left, h, peri_min_z, lintel_left + 1, h + 1, peri_max_z), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(flank_right, h, peri_min_z, flank_right + 1, h + 1, peri_max_z), 0, 0, true);

            // Update boundaries
            plan->sites[4].boundary = make_box(lintel_left, h, peri_min_z, lintel_right, h, peri_max_z);
            plan->stage = TEMPLE_STAGE_PERIPTERAL;
            plan->complexity += 2;
            return true;
        }

        case RULE_EXPAND_STYLOBATE: {
            ModBox3D foot = compute_total_footprint(plan);
            // Add 3 stepped tiers (Crepidoma)
            for (int tier = 0; tier < 3; ++tier) {
                int expand = (3 - tier);
                int y_level = -tier - 1;
                add_node(plan, ARCH_PRIMITIVE_STEP,
                         make_box(foot.min_x - expand, y_level, foot.min_z - expand,
                                  foot.max_x + expand, y_level + 1, foot.max_z + expand),
                         0, 0, true);
            }
            plan->has_stylobate = true;
            plan->complexity++;
            return true;
        }

        case RULE_RAISE_PEDIMENT: {
            ModBox3D roof_box = plan->sites[4].boundary;
            int y_base = plan->column_height_m + 1;
            int pediment_height = (roof_box.max_x - roof_box.min_x) / 4;
            if (pediment_height < 2) pediment_height = 2;

            // Add triangular pediment / gable roof node
            add_node(plan, ARCH_PRIMITIVE_PEDIMENT,
                     make_box(roof_box.min_x, y_base, roof_box.min_z,
                              roof_box.max_x, y_base + pediment_height, roof_box.max_z),
                     0, 0, true);

            plan->has_pediment = true;
            plan->sites[4].active = false;
            plan->complexity++;
            return true;
        }

        case RULE_ATTACH_COURTYARD: {
            ModBox3D core = compute_total_footprint(plan);
            int span_x = abs(core.min_x) > abs(core.max_x) ? abs(core.min_x) : abs(core.max_x);
            int court_min_z = core.max_z + 1;
            int court_len = 10;
            int court_max_z = court_min_z + court_len;

            // 1. Paved Courtyard Ground
            add_node(plan, ARCH_PRIMITIVE_COURTYARD, make_box(-span_x, 0, court_min_z, span_x, 1, court_max_z), 0, 0, true);

            // 2. Colonnaded Stoa Wings (Left & Right flanks)
            int flank_left = -span_x;
            int flank_right = span_x - 1;
            for (int z = court_min_z + 2; z <= court_max_z - 2; z += spacing) {
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(flank_left, 0, z, flank_left + 1, h, z + 1), 1, 0, true);
                add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(flank_right, 0, z, flank_right + 1, h, z + 1), 1, 0, true);
            }
            // Stoa Entablatures
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(flank_left, h, court_min_z, flank_left + 1, h + 1, court_max_z), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(flank_right, h, court_min_z, flank_right + 1, h + 1, court_max_z), 0, 0, true);

            // 3. Perimeter Enclosure Walls along front
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(-span_x, 0, court_max_z - 1, -2, 2, court_max_z), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_WALL, make_box(2, 0, court_max_z - 1, span_x, 2, court_max_z), 0, 0, true);

            // 4. Propylaea Gateway Columns flanking entrance gap
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(-2, 0, court_max_z - 1, -1, h, court_max_z), 1, 0, true);
            add_node(plan, ARCH_PRIMITIVE_COLUMN, make_box(1, 0, court_max_z - 1, 2, h, court_max_z), 1, 0, true);
            add_node(plan, ARCH_PRIMITIVE_LINTEL, make_box(-2, h, court_max_z - 1, 2, h + 1, court_max_z), 0, 0, true);

            plan->has_courtyard = true;
            plan->stage = TEMPLE_STAGE_SANCTUARY;
            plan->complexity += 3;
            return true;
        }

        case RULE_EXCAVATE_POOL: {
            int court_min_z = 0;
            for (int i = 0; i < plan->node_count; ++i) {
                if (plan->nodes[i].type == ARCH_PRIMITIVE_COURTYARD) {
                    court_min_z = plan->nodes[i].box.min_z;
                    break;
                }
            }
            int z_mid = court_min_z + 5;
            // Sunken basin node (centered across X=0)
            add_node(plan, ARCH_PRIMITIVE_POOL, make_box(-3, -1, z_mid - 2, 3, 0, z_mid + 2), 0, 0, true);
            plan->has_pool = true;
            plan->complexity++;
            return true;
        }

        case RULE_PLANT_GROVE: {
            int court_min_z = 0;
            for (int i = 0; i < plan->node_count; ++i) {
                if (plan->nodes[i].type == ARCH_PRIMITIVE_COURTYARD) {
                    court_min_z = plan->nodes[i].box.min_z;
                    break;
                }
            }
            int z_mid = court_min_z + 5;
            int x_tree = 4; // Center of tree at -4 and +4
            int tree_h = 4;
            // Pair 1 (Z = z_mid - 3)
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-x_tree - 1, 0, z_mid - 4, -x_tree + 1, tree_h, z_mid - 2), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(x_tree - 1, 0, z_mid - 4, x_tree + 1, tree_h, z_mid - 2), 0, 0, true);
            // Pair 2 (Z = z_mid + 3)
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-x_tree - 1, 0, z_mid + 2, -x_tree + 1, tree_h, z_mid + 4), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(x_tree - 1, 0, z_mid + 2, x_tree + 1, tree_h, z_mid + 4), 0, 0, true);

            plan->has_trees = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_LAYOUT_GARDEN_PATHS: {
            // 1. Outer lateral avenues along left and right flanks
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(-14, 0, -17, -12, 1, 18), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(12, 0, -17, 14, 1, 18), 0, 0, true);

            // 2. Transverse connecting paths
            // Rear walkway behind temple (leading to Tholos)
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(-12, 0, -14, 12, 1, -12), 0, 0, true);
            // Front walkway in front of Propylaea gate
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(-12, 0, 18, 12, 1, 20), 0, 0, true);
            // Mid-flank connectors to ambulatory
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(-12, 0, -1, -8, 1, 1), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_PATH, make_box(8, 0, -1, 12, 1, 1), 0, 0, true);

            plan->has_garden_paths = true;
            plan->complexity += 3;
            return true;
        }

        case RULE_ERECT_THOLOS: {
            // Circular monopteros rotunda in rear glade (centered on X=0)
            add_node(plan, ARCH_PRIMITIVE_THOLOS, make_box(-3, 0, -18, 3, 5, -12), 0, 0, true);
            plan->has_tholos = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_INSTALL_EXEDRAE: {
            // Semicircular marble philosopher benches along flank avenues
            add_node(plan, ARCH_PRIMITIVE_EXEDRA, make_box(-17, 0, -2, -14, 2, 2), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_EXEDRA, make_box(14, 0, -2, 17, 2, 2), 0, 0, true);
            plan->has_exedrae = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_BUILD_NAISKOI: {
            // Miniature votive temple shrines along flank avenues
            add_node(plan, ARCH_PRIMITIVE_NAISKOS, make_box(-17, 0, 4, -14, 4, 8), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_NAISKOS, make_box(14, 0, 4, 17, 4, 8), 0, 0, true);
            plan->has_naiskoi = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_CONSTRUCT_FOUNTAINS: {
            // Secondary marble fountain basins with PBF fluid
            add_node(plan, ARCH_PRIMITIVE_FOUNTAIN, make_box(-17, 0, -8, -14, 2, -4), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_FOUNTAIN, make_box(14, 0, -8, 17, 2, -4), 0, 0, true);
            plan->has_fountains = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_BUILD_PERGOLAS: {
            // Shaded post-and-beam vine trellises along forward avenues
            add_node(plan, ARCH_PRIMITIVE_PERGOLA, make_box(-14, 0, 10, -12, 3, 16), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_PERGOLA, make_box(12, 0, 10, 14, 3, 16), 0, 0, true);
            plan->has_pergolas = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_ERECT_BRAZIERS: {
            // Ceremonial fire altars
            // Pair 1: Flanking Propylaea front portal
            add_node(plan, ARCH_PRIMITIVE_BRAZIER, make_box(-4, 0, 19, -3, 2, 20), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_BRAZIER, make_box(3, 0, 19, 4, 2, 20), 0, 0, true);
            // Pair 2: Flanking Tholos rear entrance
            add_node(plan, ARCH_PRIMITIVE_BRAZIER, make_box(-2, 0, -12, -1, 2, -11), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_BRAZIER, make_box(1, 0, -12, 2, 2, -11), 0, 0, true);
            plan->has_braziers = true;
            plan->complexity += 2;
            return true;
        }

        case RULE_POPULATE_GARDEN_FLORA: {
            // 1. Tall Cypress Trees (slender column form, height 6)
            // Rear glade behind Tholos
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-6, 0, -17, -4, 6, -15), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(4, 0, -17, 6, 6, -15), 0, 0, true);
            // Outer corner sentinels
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-18, 0, -17, -16, 6, -15), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(16, 0, -17, 18, 6, -15), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-18, 0, 16, -16, 6, 18), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(16, 0, 16, 18, 6, 18), 0, 0, true);
            // Flank sentinels
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-18, 0, -1, -16, 6, 1), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(16, 0, -1, 18, 6, 1), 0, 0, true);

            // 2. Ancient Sprawling Olive Trees (wide canopy, height 4)
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-10, 0, -7, -8, 4, -5), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(8, 0, -7, 10, 4, -5), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-10, 0, 4, -8, 4, 6), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(8, 0, 4, 10, 4, 6), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-18, 0, -11, -16, 4, -9), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(16, 0, -11, 18, 4, -9), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(-18, 0, 9, -16, 4, 11), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_TREE, make_box(16, 0, 9, 18, 4, 11), 0, 0, true);

            // 3. Flowering Shrub Beds & Boxwood Borders
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(-11, 0, -6, -10, 1, -2), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(10, 0, -6, 11, 1, -2), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(-11, 0, 2, -10, 1, 6), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(10, 0, 2, 11, 1, 6), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(-11, 0, -16, -7, 1, -15), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(7, 0, -16, 11, 1, -15), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(-11, 0, 16, -7, 1, 17), 0, 0, true);
            add_node(plan, ARCH_PRIMITIVE_SHRUB, make_box(7, 0, 16, 11, 1, 17), 0, 0, true);

            plan->has_garden_flora = true;
            plan->complexity += 4;
            return true;
        }

        default:
            return false;
    }
}

// Bi-Algebra Developmental Loop
TemplePlan generate_greek_temple(uint32_t seed, TempleStage target_stage, int module_voxels) {
    TemplePlan plan;
    memset(&plan, 0, sizeof(plan));
    plan.seed = seed;
    plan.module_voxels = (module_voxels > 0) ? module_voxels : 2;

    // Proportions adhering to classical ratios:
    // Cella width: 6M, length: 12M (2:1 ratio)
    plan.column_spacing_m = 2;
    plan.column_height_m = 6;
    plan.cella_width_m = 6;
    plan.cella_length_m = 12;
    plan.cella_height_m = 6;

    uint32_t prng = seed;

    // 1. Initial seed S_0: sacred rectangular cella
    make_initial_cella(&plan);

    // 2. Developmental stepping loop
    int max_iterations = 45;
    while (max_iterations-- > 0) {
        // Stop condition: reached target stage, has stylobate and pediment
        if (plan.stage >= target_stage && plan.has_pediment && plan.has_stylobate) {
            if (target_stage < TEMPLE_STAGE_SANCTUARY || 
                (plan.has_courtyard && plan.has_pool && plan.has_trees &&
                 plan.has_garden_paths && plan.has_tholos && plan.has_exedrae &&
                 plan.has_naiskoi && plan.has_fountains && plan.has_pergolas &&
                 plan.has_braziers && plan.has_garden_flora)) {
                break;
            }
        }

        // Coalgebra: observe current state and generate opportunities
        GrowthOpportunity opps[MAX_OPPORTUNITIES];
        int count = greek_coalgebra_inspect(&plan, opps, MAX_OPPORTUNITIES);
        if (count == 0) break;

        // If a rule would push the building past the requested target_stage, filter it out
        int chosen_idx = -1;
        float total_weight = 0.0f;
        for (int i = 0; i < count; ++i) {
            GrammarRule r = opps[i].rule;
            if (target_stage < TEMPLE_STAGE_SANCTUARY && (
                r == RULE_ATTACH_COURTYARD || r == RULE_EXCAVATE_POOL || r == RULE_PLANT_GROVE ||
                r == RULE_LAYOUT_GARDEN_PATHS || r == RULE_ERECT_THOLOS || r == RULE_INSTALL_EXEDRAE ||
                r == RULE_BUILD_NAISKOI || r == RULE_CONSTRUCT_FOUNTAINS || r == RULE_BUILD_PERGOLAS ||
                r == RULE_ERECT_BRAZIERS || r == RULE_POPULATE_GARDEN_FLORA)) {
                opps[i].score = 0.0f;
            } else if (target_stage < TEMPLE_STAGE_PERIPTERAL && r == RULE_WRAP_PERISTYLE) {
                opps[i].score = 0.0f;
            } else if (target_stage < TEMPLE_STAGE_AMPHIPROSTYLE && r == RULE_ATTACH_OPISTHODOMOS) {
                opps[i].score = 0.0f;
            } else if (target_stage < TEMPLE_STAGE_PROSTYLE && r == RULE_ATTACH_PRONAOS) {
                opps[i].score = 0.0f;
            } else if (plan.stage >= target_stage && 
                       r != RULE_RAISE_PEDIMENT && r != RULE_EXPAND_STYLOBATE && 
                       r != RULE_EXCAVATE_POOL && r != RULE_PLANT_GROVE &&
                       r != RULE_LAYOUT_GARDEN_PATHS && r != RULE_ERECT_THOLOS &&
                       r != RULE_INSTALL_EXEDRAE && r != RULE_BUILD_NAISKOI &&
                       r != RULE_CONSTRUCT_FOUNTAINS && r != RULE_BUILD_PERGOLAS &&
                       r != RULE_ERECT_BRAZIERS && r != RULE_POPULATE_GARDEN_FLORA) {
                opps[i].score = 0.0f;
            }
            total_weight += opps[i].score;
        }

        if (total_weight <= 0.0f) break;

        // Weighted random selection
        float r = prng_uniform(&prng) * total_weight;
        float accum = 0.0f;
        for (int i = 0; i < count; ++i) {
            if (opps[i].score <= 0.0f) continue;
            accum += opps[i].score;
            if (r <= accum) {
                chosen_idx = i;
                break;
            }
        }

        if (chosen_idx >= 0) {
            // Algebra: apply chosen rule to synthesize new state S_{t+1}
            greek_algebra_apply(&plan, &opps[chosen_idx]);
        }
    }

    // Guarantee essential completing elements if loop terminated early
    if (!plan.has_stylobate) {
        GrowthOpportunity opp = { -1, -1, RULE_EXPAND_STYLOBATE, 1.0f };
        greek_algebra_apply(&plan, &opp);
    }
    if (!plan.has_pediment) {
        GrowthOpportunity opp = { 4, -1, RULE_RAISE_PEDIMENT, 1.0f };
        greek_algebra_apply(&plan, &opp);
    }
    if (target_stage == TEMPLE_STAGE_SANCTUARY) {
        if (!plan.has_courtyard) {
            GrowthOpportunity opp = { 0, -1, RULE_ATTACH_COURTYARD, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_pool) {
            GrowthOpportunity opp = { -1, -1, RULE_EXCAVATE_POOL, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_trees) {
            GrowthOpportunity opp = { -1, -1, RULE_PLANT_GROVE, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_garden_paths) {
            GrowthOpportunity opp = { -1, -1, RULE_LAYOUT_GARDEN_PATHS, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_tholos) {
            GrowthOpportunity opp = { -1, -1, RULE_ERECT_THOLOS, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_exedrae) {
            GrowthOpportunity opp = { -1, -1, RULE_INSTALL_EXEDRAE, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_naiskoi) {
            GrowthOpportunity opp = { -1, -1, RULE_BUILD_NAISKOI, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_fountains) {
            GrowthOpportunity opp = { -1, -1, RULE_CONSTRUCT_FOUNTAINS, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_pergolas) {
            GrowthOpportunity opp = { -1, -1, RULE_BUILD_PERGOLAS, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_braziers) {
            GrowthOpportunity opp = { -1, -1, RULE_ERECT_BRAZIERS, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
        if (!plan.has_garden_flora) {
            GrowthOpportunity opp = { -1, -1, RULE_POPULATE_GARDEN_FLORA, 1.0f };
            greek_algebra_apply(&plan, &opp);
        }
    }

    return plan;
}

// ASCII Blueprint Diagnostic Renderer
void print_temple_blueprint(const TemplePlan *plan) {
    ModBox3D bounds = compute_total_footprint(plan);
    int pad = 3;
    int min_x = bounds.min_x - pad;
    int max_x = bounds.max_x + pad;
    int min_z = bounds.min_z - pad;
    int max_z = bounds.max_z + pad;

    int w = max_x - min_x + 1;
    int h = max_z - min_z + 1;
    if (w > 120 || h > 120) {
        printf("Temple too large for ASCII preview (%dx%d)\n", w, h);
        return;
    }

    char *grid = (char *)malloc(w * h);
    memset(grid, ' ', w * h);

    // Rasterize 2D symbols
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        char sym = ' ';
        switch (n->type) {
            case ARCH_PRIMITIVE_STEP:      sym = '='; break;
            case ARCH_PRIMITIVE_CELL:      sym = '.'; break;
            case ARCH_PRIMITIVE_WALL:      sym = '#'; break;
            case ARCH_PRIMITIVE_COLUMN:    sym = 'O'; break;
            case ARCH_PRIMITIVE_LINTEL:    sym = '-'; break;
            case ARCH_PRIMITIVE_PEDIMENT:  sym = '^'; break;
            case ARCH_PRIMITIVE_COURTYARD: sym = '.'; break;
            case ARCH_PRIMITIVE_POOL:      sym = '~'; break;
            case ARCH_PRIMITIVE_TREE:      sym = '*'; break;
            case ARCH_PRIMITIVE_THOLOS:    sym = '@'; break;
            case ARCH_PRIMITIVE_EXEDRA:    sym = 'C'; break;
            case ARCH_PRIMITIVE_NAISKOS:   sym = 'n'; break;
            case ARCH_PRIMITIVE_FOUNTAIN:  sym = '$'; break;
            case ARCH_PRIMITIVE_BRAZIER:   sym = '!'; break;
            case ARCH_PRIMITIVE_PERGOLA:   sym = '%'; break;
            case ARCH_PRIMITIVE_PATH:      sym = '+'; break;
            case ARCH_PRIMITIVE_SHRUB:     sym = ','; break;
            default: continue;
        }

        for (int z = n->box.min_z; z < n->box.max_z; ++z) {
            for (int x = n->box.min_x; x < n->box.max_x; ++x) {
                int gx = x - min_x;
                int gz = z - min_z;
                if (gx >= 0 && gx < w && gz >= 0 && gz < h) {
                    char prev = grid[gz * w + gx];
                    // Precedence: Tree > Tholos > Column > Wall > Fountain > Brazier > Pool > Bench > Path > Floor
                    if (prev == '*' || prev == '@' || prev == 'O') continue;
                    if (prev == '#' && sym != 'O' && sym != '*' && sym != '@') continue;
                    grid[gz * w + gx] = sym;
                }
            }
        }
    }

    const char *stage_name = "Cella";
    if (plan->stage == TEMPLE_STAGE_PROSTYLE) stage_name = "Prostyle (Front Porch)";
    else if (plan->stage == TEMPLE_STAGE_AMPHIPROSTYLE) stage_name = "Amphiprostyle (Front & Rear Porches)";
    else if (plan->stage == TEMPLE_STAGE_PERIPTERAL) stage_name = "Peripteral (Full Colonnade)";
    else if (plan->stage == TEMPLE_STAGE_SANCTUARY) stage_name = "Sanctuary (Temenos, Pool & Grove)";

    printf("\n=======================================================\n");
    printf(" TEMPLE BLUEPRINT: Stage = %s, Nodes = %d\n", stage_name, plan->node_count);
    printf(" Key: [O] Column, [#] Wall, [.] Floor, [=] Steps, [~] Pool, [*] Tree,\n");
    printf("      [@] Tholos, [C] Exedra, [n] Naiskos, [$] Fountain, [!] Brazier, [%%] Pergola, [+] Path\n");
    printf("=======================================================\n");

    for (int z = h - 1; z >= 0; --z) {
        printf("%2d | ", z + min_z);
        for (int x = 0; x < w; ++x) {
            putchar(grid[z * w + x]);
        }
        putchar('\n');
    }
    printf("     ");
    for (int x = 0; x < w; ++x) putchar('-');
    printf("\n     ");
    for (int x = 0; x < w; x += 4) printf("%-4d", x + min_x);
    printf("\n\n");

    free(grid);
}
