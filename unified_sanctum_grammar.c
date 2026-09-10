#include "unified_sanctum_grammar.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// ---------------------------------------------------------------------------
// Internal Helper: Add Node to Citadel Plan
// ---------------------------------------------------------------------------
static int sanctum_add_node(SanctumCitadelPlan *plan,
                            SanctumPrimitiveType prim_type,
                            SanctumMotif motif,
                            int x, int y, int z,
                            int w, int h, int d,
                            float cant_x, float cant_z,
                            Color color,
                            bool is_emissive,
                            bool is_void,
                            bool is_fluid) {
    if (!plan || plan->node_count >= MAX_SANCTUM_NODES) return -1;
    int idx = plan->node_count++;
    SanctumNode *n = &plan->nodes[idx];
    n->prim_type = prim_type;
    n->motif = motif;
    n->x = x; n->y = y; n->z = z;
    n->w = w; n->h = h; n->d = d;
    n->cant_angle_x = cant_x;
    n->cant_angle_z = cant_z;
    n->color = color;
    n->is_emissive = is_emissive;
    n->is_void = is_void;
    n->is_fluid = is_fluid;

    // Update spatial bounds
    if (x - w/2 < plan->min_x) plan->min_x = x - w/2;
    if (x + w/2 > plan->max_x) plan->max_x = x + w/2;
    if (y < plan->min_y) plan->min_y = y;
    if (y + h > plan->max_y) plan->max_y = y + h;
    if (z - d/2 < plan->min_z) plan->min_z = z - d/2;
    if (z + d/2 > plan->max_z) plan->max_z = z + d/2;

    // Census tracking
    if (is_void) plan->count_void++;
    else if (is_fluid) plan->count_fluid++;
    else if (prim_type == SANCTUM_PRIM_HARDLIGHT_BRIDGE) plan->count_hardlight++;
    else if (prim_type == SANCTUM_PRIM_STONE_ORTHOSTAT || prim_type == SANCTUM_PRIM_STONE_LINTEL) plan->count_stone++;
    else if (prim_type == SANCTUM_PRIM_MARBLE_STYLOBATE || prim_type == SANCTUM_PRIM_MARBLE_COLUMN ||
             prim_type == SANCTUM_PRIM_MARBLE_ARCHITRAVE || prim_type == SANCTUM_PRIM_MARBLE_PEDIMENT) plan->count_marble++;
    else if (prim_type == SANCTUM_PRIM_TITANIUM_PYLON || prim_type == SANCTUM_PRIM_TITANIUM_BUTTRESS) plan->count_titanium++;

    return idx;
}

// ---------------------------------------------------------------------------
// Internal Helper: Register Open Connection Socket
// ---------------------------------------------------------------------------
static int sanctum_add_socket(SanctumCitadelPlan *plan,
                              Vector3 pos, Vector3 dir,
                              SanctumSocketType type,
                              int parent_node_idx,
                              int elevation) {
    if (!plan || plan->socket_count >= MAX_SANCTUM_SOCKETS) return -1;
    int idx = plan->socket_count++;
    SanctumSocket *s = &plan->sockets[idx];
    s->pos = pos;
    s->dir = dir;
    s->type = type;
    s->parent_node_idx = parent_node_idx;
    s->elevation = elevation;
    s->is_occupied = false;
    return idx;
}

// ---------------------------------------------------------------------------
// Plan Initialization: Root Propylaea & Core Nexus
// ---------------------------------------------------------------------------
void sanctum_plan_init(SanctumCitadelPlan *plan, uint32_t seed) {
    if (!plan) return;
    memset(plan, 0, sizeof(SanctumCitadelPlan));
    plan->seed = seed;

    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Weathered granite sarsen
    Color c_marble   = (Color){ 232, 230, 222, 255 }; // Pentelic white marble
    Color c_titanium = (Color){ 62, 68, 76, 255 };    // Brutalist dark pewter titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Emissive cyan conduit
    Color c_water    = (Color){ 50, 160, 220, 255 };  // Sacred PBF water

    // 1. Central Stepped Marble Stylobate Podium (14 x 2 x 14)
    int root = sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_CORE_NEXUS,
                                0, 0, 0, 14, 2, 14, 0.0f, 0.0f, c_marble, false, false, false);

    // 2. Central Sunken Reflecting Water Basin (6 x 1 x 6)
    sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_CORE_NEXUS,
                     0, 1, 0, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

    // 3. Four Concentric Fluted Doric Columns at diagonal corners
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_CORE_NEXUS,
                     -4, 2, -4, 2, 6, 2, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_CORE_NEXUS,
                      4, 2, -4, 2, 6, 2, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_CORE_NEXUS,
                     -4, 2,  4, 2, 6, 2, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_CORE_NEXUS,
                      4, 2,  4, 2, 6, 2, 0.0f, 0.0f, c_marble, false, false, false);

    // 4. Four Cyclopean Sarsen Uprights at cardinal edges with Titanium Clamps
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_CORE_NEXUS,
                     0, 2, -6, 4, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_CORE_NEXUS,
                     0, 2,  6, 4, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_CORE_NEXUS,
                     -6, 2, 0, 2, 7, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_CORE_NEXUS,
                      6, 2, 0, 2, 7, 4, 0.0f, 0.0f, c_stone, false, false, false);

    // Titanium Ring Architrave binding the tops of the megaliths
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CORE_NEXUS,
                     0, 9, 0, 12, 1, 12, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_CORE_NEXUS,
                     0, 9, 0, 8, 1, 8, 0.0f, 0.0f, c_cyan, true, false, false);

    // Initial Boundary Connection Sockets
    sanctum_add_socket(plan, (Vector3){  0.0f, 0.0f,  9.0f }, (Vector3){  0.0f, 0.0f,  1.0f }, SANCTUM_SOCKET_AXIAL_PATH, root, 0);
    sanctum_add_socket(plan, (Vector3){  0.0f, 0.0f, -9.0f }, (Vector3){  0.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_AXIAL_PATH, root, 0);
    sanctum_add_socket(plan, (Vector3){  9.0f, 0.0f,  0.0f }, (Vector3){  1.0f, 0.0f,  0.0f }, SANCTUM_SOCKET_AXIAL_PATH, root, 0);
    sanctum_add_socket(plan, (Vector3){ -9.0f, 0.0f,  0.0f }, (Vector3){ -1.0f, 0.0f,  0.0f }, SANCTUM_SOCKET_AXIAL_PATH, root, 0);
    sanctum_add_socket(plan, (Vector3){  0.0f, 10.0f, 0.0f }, (Vector3){  0.0f, 1.0f,  0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, root, 10);
    sanctum_add_socket(plan, (Vector3){  0.0f, -1.0f, 0.0f }, (Vector3){  0.0f, -1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, root, -1);
}

// ---------------------------------------------------------------------------
// Coalgebra γ(S, F): Inspect Open Sockets and Afford Continuations
// ---------------------------------------------------------------------------
int sanctum_coalgebra_frontier(const SanctumCitadelPlan *plan, SanctumOpportunity *out_opps, int max_opps) {
    if (!plan || !out_opps || max_opps <= 0) return 0;
    int count = 0;

    for (int i = 0; i < plan->socket_count && count < max_opps; ++i) {
        const SanctumSocket *s = &plan->sockets[i];
        if (s->is_occupied) continue;

        // High entropy pseudo-random scoring hash based on seed and socket index
        uint32_t h = (plan->seed + (uint32_t)i * 2654435761u) ^ (uint32_t)(s->pos.x * 37 + s->pos.z * 101);
        h ^= (h << 13); h ^= (h >> 17); h ^= (h << 5);

        switch (s->type) {
            case SANCTUM_SOCKET_AXIAL_PATH: {
                // Afford 1: Megalithic Titanium Trilithon Gateway
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_MEGALITHIC_TRILITHON;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.85f + (float)(h % 15) / 100.0f;
                }
                // Afford 2: Doric Colonnade Terrace Avenue
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_DORIC_COLONNADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.90f + (float)((h >> 4) % 15) / 100.0f;
                }
                // Afford 3: Canted Titanium Pylon Portal
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_CANTED_PYLON_PORTAL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.88f + (float)((h >> 8) % 15) / 100.0f;
                }
                // Afford 4: Abyssal Chasm with Hard-Light Bridge
                if (count < max_opps && plan->count_hardlight < (plan->node_count / 4 + 1)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_HARDLIGHT_CHASM;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.92f + (float)((h >> 12) % 15) / 100.0f;
                }
                // Afford 5: Sacred PBF Cascade
                if (count < max_opps && (h % 3 == 0)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_SACRED_PBF_CASCADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.86f;
                }
                break;
            }

            case SANCTUM_SOCKET_COLONNADE_FLANK: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_DORIC_COLONNADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.85f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_CANTED_PYLON_PORTAL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.80f;
                }
                break;
            }

            case SANCTUM_SOCKET_PORTAL_GATE: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_THOLOS_GRAVITY_PIT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.95f;
                }
                break;
            }

            case SANCTUM_SOCKET_VERTICAL_APEX: {
                if (s->dir.y > 0.0f && count < max_opps) {
                    // Zenith: Pinnacle Spire with Levitating Core
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_APEX_SPIRE_MATRIX;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.98f;
                } else if (s->dir.y < 0.0f && count < max_opps) {
                    // Nadir: Subterranean Corbelled Crypt Vault
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_CORBELLED_CRYPT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.94f;
                }
                break;
            }

            case SANCTUM_SOCKET_HYDRAULIC_GUTTER: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_SACRED_PBF_CASCADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.90f;
                }
                break;
            }

            default:
                break;
        }
    }

    return count;
}

// ---------------------------------------------------------------------------
// Algebra α(S, F, opp): Snap Motif to Socket and Expand Structure
// ---------------------------------------------------------------------------
bool sanctum_algebra_expand(SanctumCitadelPlan *plan, const SanctumOpportunity *opp) {
    if (!plan || !opp || opp->socket_idx < 0 || opp->socket_idx >= plan->socket_count) return false;

    SanctumSocket *target_socket = &plan->sockets[opp->socket_idx];
    if (target_socket->is_occupied) return false;
    target_socket->is_occupied = true;

    Vector3 pos = opp->spawn_pos;
    Vector3 dir = opp->spawn_dir;
    int px = (int)roundf(pos.x);
    int py = (int)roundf(pos.y);
    int pz = (int)roundf(pos.z);

    // Directional orientation: is forward along Z or X?
    bool along_z = fabsf(dir.z) >= fabsf(dir.x);
    int step_x = along_z ? 0 : (dir.x > 0.0f ? 1 : -1);
    int step_z = along_z ? (dir.z > 0.0f ? 1 : -1) : 0;

    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Sarsen stone
    Color c_marble   = (Color){ 235, 232, 222, 255 }; // Pentelic marble
    Color c_titanium = (Color){ 58, 64, 72, 255 };    // Pewter titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Emissive cyan
    Color c_gold     = (Color){ 255, 185, 45, 255 };  // Amber/Gold
    Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Fluid

    switch (opp->motif) {
        case SANCTUM_MOTIF_MEGALITHIC_TRILITHON: {
            int span = 8;
            int length = 6;
            int h_post = 8;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Left and Right Orthostats (Rough Sarsen Megaliths)
            int ox1 = along_z ? center_x - 3 : center_x;
            int oz1 = along_z ? center_z : center_z - 3;
            int ox2 = along_z ? center_x + 3 : center_x;
            int oz2 = along_z ? center_z : center_z + 3;

            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             ox1, py, oz1, 2, h_post, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             ox2, py, oz2, 2, h_post, 2, 0.0f, 0.0f, c_stone, false, false, false);

            // Titanium Clamping Collars
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             ox1, py + h_post - 2, oz1, 3, 2, 3, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             ox2, py + h_post - 2, oz2, 3, 2, 3, 0.0f, 0.0f, c_titanium, false, false, false);

            // Megalithic Capstone Lintel
            int lw = along_z ? span : 3;
            int ld = along_z ? 3 : span;
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_LINTEL, opp->motif,
                             center_x, py + h_post, center_z, lw, 2, ld, 0.0f, 0.0f, c_stone, false, false, false);

            // Cyan Light Inlay across lintel
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + h_post + 1, center_z, lw - 2, 1, ld - 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // Forward Socket
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 2)), (float)py, (float)(pz + step_z * (length + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_DORIC_COLONNADE: {
            int length = 12;
            int width = 8;
            int col_h = 7;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Stepped Pentelic Marble Stylobate
            int sw = along_z ? width : length;
            int sd = along_z ? length : width;
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, sw, 1, sd, 0.0f, 0.0f, c_marble, false, false, false);

            // Fluted Columns on Left and Right flanks
            for (int k = -1; k <= 1; ++k) {
                int cx1 = along_z ? (center_x - 3) : (center_x + k * 4);
                int cz1 = along_z ? (center_z + k * 4) : (center_z - 3);
                int cx2 = along_z ? (center_x + 3) : (center_x + k * 4);
                int cz2 = along_z ? (center_z + k * 4) : (center_z + 3);

                sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                 cx1, py + 1, cz1, 1, col_h, 1, 0.0f, 0.0f, c_marble, false, false, false);
                sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                 cx2, py + 1, cz2, 1, col_h, 1, 0.0f, 0.0f, c_marble, false, false, false);
            }

            // Continuous Architrave Frieze
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, opp->motif,
                             center_x, py + col_h + 1, center_z, sw, 1, sd, 0.0f, 0.0f, c_marble, false, false, false);

            // Glowing conduit running down the center aisle
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py, center_z, along_z ? 1 : length, 1, along_z ? length : 1, 0.0f, 0.0f, c_cyan, true, false, false);

            // Forward Socket & Lateral Flank Sockets
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 1)), (float)py, (float)(pz + step_z * (length + 1)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);

            Vector3 l_dir = along_z ? (Vector3){ -1.0f, 0.0f, 0.0f } : (Vector3){ 0.0f, 0.0f, -1.0f };
            Vector3 r_dir = along_z ? (Vector3){  1.0f, 0.0f, 0.0f } : (Vector3){ 0.0f, 0.0f,  1.0f };
            sanctum_add_socket(plan, (Vector3){ (float)(center_x + l_dir.x * 5), (float)py, (float)(center_z + l_dir.z * 5) },
                               l_dir, SANCTUM_SOCKET_COLONNADE_FLANK, plan->node_count - 1, py);
            sanctum_add_socket(plan, (Vector3){ (float)(center_x + r_dir.x * 5), (float)py, (float)(center_z + r_dir.z * 5) },
                               r_dir, SANCTUM_SOCKET_COLONNADE_FLANK, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_CANTED_PYLON_PORTAL: {
            int length = 8;
            int span = 10;
            int p_height = 10;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Twin Canted Pewter Pylons
            int px1 = along_z ? center_x - 4 : center_x;
            int pz1 = along_z ? center_z : center_z - 4;
            int px2 = along_z ? center_x + 4 : center_x;
            int pz2 = along_z ? center_z : center_z + 4;

            float cant_x1 = along_z ? -22.0f : 0.0f;
            float cant_z1 = along_z ? 0.0f : -22.0f;
            float cant_x2 = along_z ?  22.0f : 0.0f;
            float cant_z2 = along_z ? 0.0f :  22.0f;

            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             px1, py, pz1, 3, p_height, 3, cant_x1, cant_z1, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             px2, py, pz2, 3, p_height, 3, cant_x2, cant_z2, c_titanium, false, false, false);

            // Overhead Chevron Lintel Span
            int cw = along_z ? span : 3;
            int cd = along_z ? 3 : span;
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py + p_height, center_z, cw, 2, cd, 0.0f, 0.0f, c_titanium, false, false, false);

            // Emissive Glyph Seams
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + p_height + 1, center_z, cw - 2, 1, cd - 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // Forward Socket & Zenith Socket
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 2)), (float)py, (float)(pz + step_z * (length + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            sanctum_add_socket(plan,
                               (Vector3){ (float)center_x, (float)(py + p_height + 3), (float)center_z },
                               (Vector3){ 0.0f, 1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, plan->node_count - 1, py + p_height);
            break;
        }

        case SANCTUM_MOTIF_HARDLIGHT_CHASM: {
            int chasm_len = 16;
            int chasm_w = 12;
            int chasm_depth = 10;
            int center_x = px + step_x * (chasm_len / 2);
            int center_z = pz + step_z * (chasm_len / 2);

            // 1. Excavate Abyssal Chasm Void (Pass 1 in rasterizer)
            int vw = along_z ? chasm_w : chasm_len;
            int vd = along_z ? chasm_len : chasm_w;
            sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, opp->motif,
                             center_x, py - chasm_depth, center_z, vw, chasm_depth, vd, 0.0f, 0.0f, BLACK, false, true, false);

            // 2. Twin Solid Abutments at near and far edge
            int ax_near = along_z ? center_x : px + step_x;
            int az_near = along_z ? pz + step_z : center_z;
            int ax_far  = along_z ? center_x : px + step_x * (chasm_len - 1);
            int az_far  = along_z ? pz + step_z * (chasm_len - 1) : center_z;

            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             ax_near, py - 2, az_near, along_z ? 6 : 2, 3, along_z ? 2 : 6, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             ax_far, py - 2, az_far, along_z ? 6 : 2, 3, along_z ? 2 : 6, 0.0f, 0.0f, c_titanium, false, false, false);

            // 3. Suspended Cyan Luminescent Hard-Light Bridge Span
            int bw = along_z ? 4 : chasm_len;
            int bd = along_z ? chasm_len : 4;
            sanctum_add_node(plan, SANCTUM_PRIM_HARDLIGHT_BRIDGE, opp->motif,
                             center_x, py, center_z, bw, 1, bd, 0.0f, 0.0f, (Color){ 45, 225, 255, 230 }, true, false, false);

            // Forward Socket on other side of chasm
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (chasm_len + 2)), (float)py, (float)(pz + step_z * (chasm_len + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_THOLOS_GRAVITY_PIT: {
            int radius = 8;
            int pit_depth = 8;
            int center_x = px + step_x * (radius + 2);
            int center_z = pz + step_z * (radius + 2);

            // 1. Excavate Central Subterranean Gravity Pit
            sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, opp->motif,
                             center_x, py - pit_depth, center_z, 10, pit_depth, 10, 0.0f, 0.0f, BLACK, false, true, false);

            // 2. Annular PBF Reflection Fluid Moat around rim
            sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
                             center_x, py, center_z, 14, 1, 14, 0.0f, 0.0f, c_water, false, false, true);

            // 3. Concentric Fluted Doric Columns & Megalith Uprights in circle
            for (int a = 0; a < 8; ++a) {
                float rad = (float)a * (2.0f * PI / 8.0f);
                int cx = center_x + (int)roundf(cosf(rad) * 6.0f);
                int cz = center_z + (int)roundf(sinf(rad) * 6.0f);
                if (a % 2 == 0) {
                    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                     cx, py + 1, cz, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
                } else {
                    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                                     cx, py + 1, cz, 2, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
                }
            }

            // 4. Levitating Octahedral Gravity Core at pit center
            sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, opp->motif,
                             center_x, py + 2, center_z, 2, 2, 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // 5. Overhead Titanium Dome Architrave
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py + 8, center_z, 12, 1, 12, 0.0f, 0.0f, c_titanium, false, false, false);

            // Outward Sockets in 3 remaining cardinal directions
            Vector3 left_dir  = (Vector3){ -dir.z, 0.0f,  dir.x };
            Vector3 right_dir = (Vector3){  dir.z, 0.0f, -dir.x };
            sanctum_add_socket(plan,
                               (Vector3){ (float)(center_x + dir.x * (radius + 2)), (float)py, (float)(center_z + dir.z * (radius + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            sanctum_add_socket(plan,
                               (Vector3){ (float)(center_x + left_dir.x * (radius + 2)), (float)py, (float)(center_z + left_dir.z * (radius + 2)) },
                               left_dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            sanctum_add_socket(plan,
                               (Vector3){ (float)(center_x + right_dir.x * (radius + 2)), (float)py, (float)(center_z + right_dir.z * (radius + 2)) },
                               right_dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            sanctum_add_socket(plan,
                               (Vector3){ (float)center_x, (float)(py + 9), (float)center_z },
                               (Vector3){ 0.0f, 1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, plan->node_count - 1, py + 8);
            break;
        }

        case SANCTUM_MOTIF_APEX_SPIRE_MATRIX: {
            int spire_h = 14;
            int center_x = px;
            int center_z = pz;

            // 1. Towering Central Spire Pinnacle Needle
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py, center_z, 3, spire_h, 3, 0.0f, 0.0f, c_titanium, false, false, false);

            // 2. Quad Flying Buttresses Canted at 28 degrees
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, opp->motif,
                             center_x - 3, py, center_z, 2, 6, 2, -28.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, opp->motif,
                             center_x + 3, py, center_z, 2, 6, 2,  28.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, opp->motif,
                             center_x, py, center_z - 3, 2, 6, 2, 0.0f, -28.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, opp->motif,
                             center_x, py, center_z + 3, 2, 6, 2, 0.0f,  28.0f, c_titanium, false, false, false);

            // 3. Cantilevered Observation Balcony Terrace (Pentelic Marble)
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py + 4, center_z, 8, 1, 8, 0.0f, 0.0f, c_marble, false, false, false);

            // 4. Apex Amber Pulse Beacon & Levitating Core
            sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, opp->motif,
                             center_x, py + spire_h, center_z, 2, 3, 2, 0.0f, 0.0f, c_gold, true, false, false);
            break;
        }

        case SANCTUM_MOTIF_CORBELLED_CRYPT: {
            int crypt_h = 6;
            int crypt_w = 8;
            int crypt_d = 8;
            int center_x = px;
            int center_z = pz;
            int cy = py - crypt_h;

            // Excavate Subterranean Crypt Cavern
            sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, opp->motif,
                             center_x, cy, center_z, crypt_w, crypt_h, crypt_d, 0.0f, 0.0f, BLACK, false, true, false);

            // Megalithic Corbelled Walls
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - 3, cy, center_z, 2, crypt_h, crypt_d, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + 3, cy, center_z, 2, crypt_h, crypt_d, 0.0f, 0.0f, c_stone, false, false, false);

            // Central Reliquary Altar with Glowing Cyan Conduit
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_LINTEL, opp->motif,
                             center_x, cy, center_z, 2, 2, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, cy + 2, center_z, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);
            break;
        }

        case SANCTUM_MOTIF_SACRED_PBF_CASCADE: {
            int length = 10;
            int width = 6;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Stepped Marble Cascade Podium
            int mw = along_z ? width : length;
            int md = along_z ? length : width;
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, mw, 1, md, 0.0f, 0.0f, c_marble, false, false, false);

            // Dynamic PBF Water Flume
            int ww = along_z ? 2 : length;
            int wd = along_z ? length : 2;
            sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
                             center_x, py + 1, center_z, ww, 1, wd, 0.0f, 0.0f, c_water, false, false, true);

            // Flanking Sarsen Boundary Stele
            int sx1 = along_z ? center_x - 2 : center_x;
            int sz1 = along_z ? center_z : center_z - 2;
            int sx2 = along_z ? center_x + 2 : center_x;
            int sz2 = along_z ? center_z : center_z + 2;
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             sx1, py + 1, sz1, 1, 3, 1, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             sx2, py + 1, sz2, 1, 3, 1, 0.0f, 0.0f, c_stone, false, false, false);

            // Forward Socket
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 1)), (float)py, (float)(pz + step_z * (length + 1)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            break;
        }

        default:
            return false;
    }

    return true;
}

// ---------------------------------------------------------------------------
// Open-Ended Generation Loop: Step-by-Step Frontier Expansion
// ---------------------------------------------------------------------------
SanctumCitadelPlan generate_unified_sanctum(uint32_t seed, int growth_steps) {
    SanctumCitadelPlan plan;
    sanctum_plan_init(&plan, seed);
    plan.growth_steps = growth_steps;

    SanctumOpportunity opps[MAX_SANCTUM_OPPS];

    for (int step = 0; step < growth_steps; ++step) {
        int opp_count = sanctum_coalgebra_frontier(&plan, opps, MAX_SANCTUM_OPPS);
        if (opp_count <= 0) break;

        // Choose opportunity with highest score (or weighted pseudo-random selection)
        int best_idx = 0;
        float best_score = -1.0f;
        for (int i = 0; i < opp_count; ++i) {
            if (opps[i].score > best_score) {
                best_score = opps[i].score;
                best_idx = i;
            }
        }

        // Apply selected continuation
        sanctum_algebra_expand(&plan, &opps[best_idx]);
    }

    return plan;
}

// ---------------------------------------------------------------------------
// Invariant Verification
// ---------------------------------------------------------------------------
bool verify_sanctum_invariants(const SanctumCitadelPlan *plan) {
    if (!plan || plan->node_count <= 0) return false;

    // Invariant 1: Structural diversity — Must contain members from all 3 architectural families
    if (plan->count_stone <= 0) {
        printf("Verification failed: Missing megalithic sarsen stone primitives!\n");
        return false;
    }
    if (plan->count_marble <= 0) {
        printf("Verification failed: Missing Classical Greek Pentelic marble primitives!\n");
        return false;
    }
    if (plan->count_titanium <= 0) {
        printf("Verification failed: Missing Forerunner titanium alloy primitives!\n");
        return false;
    }

    // Invariant 2: Active connection sockets must have valid direction normals
    for (int i = 0; i < plan->socket_count; ++i) {
        const SanctumSocket *s = &plan->sockets[i];
        float len = sqrtf(s->dir.x * s->dir.x + s->dir.y * s->dir.y + s->dir.z * s->dir.z);
        if (len < 0.9f || len > 1.1f) {
            printf("Verification failed: Malformed socket direction vector at idx %d\n", i);
            return false;
        }
    }

    // Invariant 3: Bounding box integrity
    if (plan->min_x > plan->max_x || plan->min_y > plan->max_y || plan->min_z > plan->max_z) {
        printf("Verification failed: Corrupted bounding box!\n");
        return false;
    }

    return true;
}
