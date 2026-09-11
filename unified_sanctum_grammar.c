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
    else if (prim_type == SANCTUM_PRIM_STONE_ORTHOSTAT || prim_type == SANCTUM_PRIM_STONE_LINTEL ||
             prim_type == SANCTUM_PRIM_OBELISK || prim_type == SANCTUM_PRIM_ALTAR) plan->count_stone++;
    else if (prim_type == SANCTUM_PRIM_MARBLE_STYLOBATE || prim_type == SANCTUM_PRIM_MARBLE_COLUMN ||
             prim_type == SANCTUM_PRIM_MARBLE_ARCHITRAVE || prim_type == SANCTUM_PRIM_MARBLE_PEDIMENT ||
             prim_type == SANCTUM_PRIM_BALUSTRADE) plan->count_marble++;
    else if (prim_type == SANCTUM_PRIM_TITANIUM_PYLON || prim_type == SANCTUM_PRIM_TITANIUM_BUTTRESS ||
             prim_type == SANCTUM_PRIM_COFFERED_CEILING || prim_type == SANCTUM_PRIM_BRAZIER) plan->count_titanium++;

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
    float len = sqrtf(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    if (len > 0.0001f) {
        dir.x /= len;
        dir.y /= len;
        dir.z /= len;
    }
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
// Tactical Gameplay Registration Helpers
// ---------------------------------------------------------------------------
static void sanctum_add_supply(SanctumCitadelPlan *plan, Vector3 pos, SanctumSupplyType type) {
    if (!plan || plan->supply_count >= MAX_SANCTUM_SUPPLIES) return;
    plan->supply_points[plan->supply_count++] = (SanctumSupplyPoint){ pos, type };
}

static void sanctum_add_spawn(SanctumCitadelPlan *plan, Vector3 pos, float yaw, bool is_player) {
    if (!plan || plan->spawn_point_count >= MAX_SANCTUM_SPAWNS) return;
    plan->spawn_points[plan->spawn_point_count++] = (SanctumCombatSpawn){ pos, yaw, is_player };
}

static void sanctum_add_jump_pad(SanctumCitadelPlan *plan, Vector3 pos, float power) {
    if (!plan || plan->jump_pad_count >= MAX_SANCTUM_JUMP_PADS) return;
    plan->jump_pads[plan->jump_pad_count++] = (SanctumJumpPad){ pos, power };
}

// ---------------------------------------------------------------------------
// Macro-Topology 0: The Sovereign Stronghold (Central Hilltop Holdout)
// ---------------------------------------------------------------------------
static void sanctum_init_stronghold(SanctumCitadelPlan *plan) {
    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Sarsen
    Color c_marble   = (Color){ 232, 230, 222, 255 }; // Marble
    Color c_titanium = (Color){ 62, 68, 76, 255 };    // Titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Cyan
    // Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Water (disabled for now)
    Color c_gold     = (Color){ 255, 185, 45, 255 };  // Gold

    // Tier 0 Foundation Podium: 52 x 4 x 52
    int base_podium = sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                                       0, 0, 0, 52, 4, 52, 0.0f, 0.0f, c_stone, false, false, false);

    // Four Fortified Pewter Corner Bastions with beacon braziers
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 0, -22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 0, -22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 0,  22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 0,  22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);

    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 8, -22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 8, -22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 8,  22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 8,  22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);

    // Tier 1 Grand Stylobate Terrace: 46 x 2 x 46
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                     0, 4, 0, 46, 2, 46, 0.0f, 0.0f, c_marble, false, false, false);

    // Cardinal Stairs
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  0, 1,  25, 14, 3, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  0, 1, -25, 14, 3, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  25, 1,  0, 6, 3, 14, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS, -25, 1,  0, 6, 3, 14, 0.0f, 0.0f, c_marble, false, false, false);

    // North Covered Stoa
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIMETER_STOA, 0, 6, 21, 38, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    for (int k = -3; k <= 3; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_PERIMETER_STOA, k * 5, 6, 17, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_PERIMETER_STOA, 0, 12, 19, 38, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_PERIMETER_STOA, 0, 13, 19, 38, 1, 6, 0.0f, 0.0f, c_titanium, false, false, false);

    // Central Holdout Sanctuary (Y=6..38)
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 6, 0, 28, 2, 28, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 8, 11, 26, 12, 3, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -12, 8, 0, 3, 12, 22, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  12, 8, 0, 3, 12, 22, 0.0f, 0.0f, c_stone, false, false, false);

    // Low Chest-High Parapets (Half-Cover)
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -9, 8, -12, 6, 2, 2, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  9, 8, -12, 6, 2, 2, 0.0f, 0.0f, c_marble, false, false, false);

    // Propylaea Portal
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -7, 8, -12, 4, 14, 4, -18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  7, 8, -12, 4, 14, 4,  18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 20, -12, 18, 3, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 21, -12, 14, 1, 2, 0.0f, 0.0f, c_cyan, true, false, false);

    // Interior Nave
    for (int k = -1; k <= 2; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -5, 8, k * 4 - 2, 2, 10, 2, 0.0f, 0.0f, c_marble, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  5, 8, k * 4 - 2, 2, 10, 2, 0.0f, 0.0f, c_marble, false, false, false);
    }
    // Sacred Healing Well inside Holdout (disabled for now)
    // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_RECHARGER_WELL, 0, 8, 2, 6, 1, 8, 0.0f, 0.0f, c_water, false, false, true);

    // Coffered Ceiling & Roof Gallery
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 20, 0, 26, 2, 26, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 20, 2, 8, 2, 8, 0.0f, 0.0f, BLACK, false, true, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 22, 0, 20, 2, 20, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 24, -9, 18, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 24,  9, 18, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

    // Apex Spire Needle (Y=26..38)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_APEX_SPIRE_MATRIX, 0, 26, 0, 4, 12, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX, -4, 24, 0, 2, 6, 2, -28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,  4, 24, 0, 2, 6, 2,  28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX, 0, 24, -4, 2, 6, 2, 0.0f, -28.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX, 0, 24,  4, 2, 6, 2, 0.0f,  28.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, SANCTUM_MOTIF_APEX_SPIRE_MATRIX, 0, 35, 0, 3, 3, 3, 0.0f, 0.0f, c_gold, true, false, false);

    // Power Weapon Altar at South Approach
    sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, SANCTUM_MOTIF_WEAPON_ALTAR, 0, 6, -11, 4, 2, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_WEAPON_ALTAR, 0, 8, -11, 2, 1, 2, 0.0f, 0.0f, c_gold, true, false, false);

    // Four Peripheral Enemy Spawn Crypts in Quadrants
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT,  18, 6,  18, 4, 5, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT, -18, 6,  18, 4, 5, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT,  18, 6, -18, 4, 5, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT, -18, 6, -18, 4, 5, 4, 0.0f, 0.0f, c_stone, false, false, false);

    // East & West Ammo Obelisk Plazas
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_OBELISK_PLAZA,  17, 6, 0, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_OBELISK, SANCTUM_MOTIF_OBELISK_PLAZA,  17, 7, 0, 2, 12, 2, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_OBELISK_PLAZA, -17, 6, 0, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_OBELISK, SANCTUM_MOTIF_OBELISK_PLAZA, -17, 7, 0, 2, 12, 2, 0.0f, 0.0f, c_stone, false, false, false);

    // Tactical Spawns & Supplies Registration
    sanctum_add_spawn(plan, (Vector3){ 0.0f, 8.0f, 2.0f }, 180.0f, true); // Player Holdout
    sanctum_add_spawn(plan, (Vector3){  18.0f, 6.0f,  18.0f }, 225.0f, false); // Enemy NE
    sanctum_add_spawn(plan, (Vector3){ -18.0f, 6.0f,  18.0f }, 135.0f, false); // Enemy NW
    sanctum_add_spawn(plan, (Vector3){  18.0f, 6.0f, -18.0f }, 315.0f, false); // Enemy SE
    sanctum_add_spawn(plan, (Vector3){ -18.0f, 6.0f, -18.0f },  45.0f, false); // Enemy SW

    sanctum_add_supply(plan, (Vector3){ 0.0f, 8.0f, 2.0f }, SANCTUM_SUPPLY_HEALTH);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 8.0f, -11.0f }, SANCTUM_SUPPLY_DYNAMIC_SHOT);
    sanctum_add_supply(plan, (Vector3){  17.0f, 6.0f, 0.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){ -17.0f, 6.0f, 0.0f }, SANCTUM_SUPPLY_AMMO);

    // Gravity Jump Lifts
    sanctum_add_jump_pad(plan, (Vector3){ 0.0f, 1.0f, -25.0f }, 18.0f);
    sanctum_add_jump_pad(plan, (Vector3){ 0.0f, 1.0f,  25.0f }, 18.0f);

    // Outward Frontier Sockets
    sanctum_add_socket(plan, (Vector3){ 0.0f, 0.0f, 29.0f }, (Vector3){ 0.0f, 0.0f, 1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);
    sanctum_add_socket(plan, (Vector3){ 0.0f, 0.0f, -29.0f }, (Vector3){ 0.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);
    sanctum_add_socket(plan, (Vector3){ 20.0f, 6.0f, 19.0f }, (Vector3){ 1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_STOA_CONTINUATION, base_podium, 6);
    sanctum_add_socket(plan, (Vector3){ -20.0f, 6.0f, 19.0f }, (Vector3){ -1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_STOA_CONTINUATION, base_podium, 6);
}

// ---------------------------------------------------------------------------
// Macro-Topology 1: The Abyssal Rift (Dual-Spire Chasm Arena)
// ---------------------------------------------------------------------------
static void sanctum_init_abyssal_rift(SanctumCitadelPlan *plan) {
    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Sarsen
    Color c_marble   = (Color){ 232, 230, 222, 255 }; // Marble
    Color c_titanium = (Color){ 62, 68, 76, 255 };    // Titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Cyan
    // Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Water (disabled for now)
    Color c_gold     = (Color){ 255, 185, 45, 255 };  // Gold

    // Central Abyssal Void Canyon bisecting the arena along X
    sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     0, -2, 0, 52, 12, 16, 0.0f, 0.0f, BLACK, false, true, false);

    // ================= SOUTH FORTRESS PLATEAU (Player Holdout Bluff) =================
    int base_south = sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                                      0, 0, -18, 48, 6, 20, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     0, 6, -18, 44, 2, 16, 0.0f, 0.0f, c_marble, false, false, false);

    // Fortified Megaron Bunker
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     0, 8, -20, 24, 10, 12, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     0, 18, -20, 24, 2, 12, 0.0f, 0.0f, c_titanium, false, false, false);

    // Low Chest-High Breastworks along the cliff edge overlooking the chasm
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     0, 8, -9, 38, 2, 1, 0.0f, 0.0f, c_marble, false, false, false);

    // Transmission Pinnacle Needle & Flying Buttresses
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     0, 18, -20, 4, 16, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     -4, 18, -20, 2, 6, 2, -28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                      4, 18, -20, 2, 6, 2,  28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     0, 32, -20, 3, 3, 3, 0.0f, 0.0f, c_gold, true, false, false);

    // South Sacred Healing Well (disabled for now)
    // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_RECHARGER_WELL,
    //                  0, 8, -22, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

    // ================= THREE CHASM BRIDGES =================
    // Central Luminescent Hard-Light Bridge
    sanctum_add_node(plan, SANCTUM_PRIM_HARDLIGHT_BRIDGE, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     0, 7, 0, 4, 1, 16, 0.0f, 0.0f, c_cyan, true, false, false);

    // West Cyclopean Megalithic Arch Bridge
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_LINTEL, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     -16, 7, 0, 4, 2, 16, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     -16, 6, -8, 5, 3, 2, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     -16, 6,  8, 5, 3, 2, 0.0f, 0.0f, c_titanium, false, false, false);

    // East Hard-Light Bridge
    sanctum_add_node(plan, SANCTUM_PRIM_HARDLIGHT_BRIDGE, SANCTUM_MOTIF_HARDLIGHT_CHASM,
                     16, 7, 0, 4, 1, 16, 0.0f, 0.0f, c_cyan, true, false, false);

    // Central Suspended Power Weapon Shrine
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_WEAPON_ALTAR,
                     0, -2, 0, 4, 7, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, SANCTUM_MOTIF_WEAPON_ALTAR,
                     0, 5, 0, 4, 2, 4, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_WEAPON_ALTAR,
                     0, 7, 0, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

    // ================= NORTH FORTRESS PLATEAU (Enemy Ingress Grounds) =================
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIBOLOS_RAMPART,
                     0, 0, 18, 48, 6, 20, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_PERIBOLOS_RAMPART,
                     0, 6, 18, 44, 2, 16, 0.0f, 0.0f, c_marble, false, false, false);

    // North Twin Pylon Portals
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CANTED_PYLON_PORTAL,
                     -10, 8, 14, 4, 14, 4, -18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CANTED_PYLON_PORTAL,
                      10, 8, 14, 4, 14, 4,  18.0f, 0.0f, c_titanium, false, false, false);

    // North Spawn Crypts
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT, 0, 8, 23, 6, 6, 6, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT, -16, 8, 21, 5, 6, 5, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT,  16, 8, 21, 5, 6, 5, 0.0f, 0.0f, c_stone, false, false, false);

    // Tactical Registrations
    sanctum_add_spawn(plan, (Vector3){ 0.0f, 8.0f, -18.0f }, 0.0f, true); // Player Holdout
    sanctum_add_spawn(plan, (Vector3){   0.0f, 8.0f, 23.0f }, 180.0f, false); // Enemy North
    sanctum_add_spawn(plan, (Vector3){ -16.0f, 8.0f, 21.0f }, 135.0f, false); // Enemy NW
    sanctum_add_spawn(plan, (Vector3){  16.0f, 8.0f, 21.0f }, 225.0f, false); // Enemy NE

    sanctum_add_supply(plan, (Vector3){ 0.0f, 8.0f, -22.0f }, SANCTUM_SUPPLY_HEALTH);
    sanctum_add_supply(plan, (Vector3){ -12.0f, 8.0f, -14.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 8.0f, 9.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 6.0f, 0.0f }, SANCTUM_SUPPLY_VOID);

    // Jump Pads inside Chasm Void launching to bridge deck
    sanctum_add_jump_pad(plan, (Vector3){ -10.0f, -1.0f, 0.0f }, 22.0f);
    sanctum_add_jump_pad(plan, (Vector3){  10.0f, -1.0f, 0.0f }, 22.0f);

    // Outward Sockets
    sanctum_add_socket(plan, (Vector3){ -24.0f, 6.0f, -18.0f }, (Vector3){ -1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_south, 6);
    sanctum_add_socket(plan, (Vector3){  24.0f, 6.0f, -18.0f }, (Vector3){  1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_south, 6);
}

// ---------------------------------------------------------------------------
// Macro-Topology 2: The Sunken Crucible (Inverted Colosseum Arena)
// ---------------------------------------------------------------------------
static void sanctum_init_sunken_crucible(SanctumCitadelPlan *plan) {
    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Sarsen
    Color c_marble   = (Color){ 232, 230, 222, 255 }; // Marble
    Color c_titanium = (Color){ 62, 68, 76, 255 };    // Titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Cyan
    // Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Water (disabled for now)

    // Outer Colosseum Ring: 52 x 8 x 52 retaining wall at Y=4
    int base_colosseum = sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                                          0, 4, 0, 52, 8, 52, 0.0f, 0.0f, c_stone, false, false, false);

    // Excavated Central Battle Bowl: 38 x 8 x 38 void down to Y=2
    sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, SANCTUM_MOTIF_THOLOS_GRAVITY_PIT,
                     0, 2, 0, 38, 8, 38, 0.0f, 0.0f, BLACK, false, true, false);

    // Promenade Colonnades lining the upper rim at Y=10
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_DORIC_COLONNADE, 0, 14,  21, 38, 1, 4, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_DORIC_COLONNADE, 0, 14, -21, 38, 1, 4, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_DORIC_COLONNADE,  21, 14, 0, 4, 1, 38, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_DORIC_COLONNADE, -21, 14, 0, 4, 1, 38, 0.0f, 0.0f, c_marble, false, false, false);

    // Colonnade Columns along upper promenade
    for (int k = -3; k <= 3; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_DORIC_COLONNADE, k * 5, 10,  21, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_DORIC_COLONNADE, k * 5, 10, -21, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }

    // Low Balustrades overlooking pit
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_DORIC_COLONNADE, 0, 11,  19, 36, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_DORIC_COLONNADE, 0, 11, -19, 36, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

    // Four Ingress Gatehouses built into outer walls
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_SPAWN_CRYPT,   0, 10,  23, 8, 6, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_SPAWN_CRYPT,   0, 10, -23, 8, 6, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_SPAWN_CRYPT,  23, 10,   0, 4, 6, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_SPAWN_CRYPT, -23, 10,   0, 4, 6, 8, 0.0f, 0.0f, c_titanium, false, false, false);

    // Sunken Arena Floor at Y=2 (disabled for now)
    // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_THOLOS_GRAVITY_PIT, 0, 2, 0, 30, 1, 30, 0.0f, 0.0f, c_water, false, false, true);

    // Central Holdout Sanctuary Island: 14 x 2 x 14 at Y=2
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 2, 0, 14, 2, 14, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 4, 0, 12, 1, 12, 0.0f, 0.0f, c_marble, false, false, false);

    // 4 Doric Columns & Titanium Canopy over central island
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -4, 5, -4, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  4, 5, -4, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -4, 5,  4, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  4, 5,  4, 1, 4, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, 0, 9, 0, 10, 1, 10, 0.0f, 0.0f, c_titanium, false, false, false);

    // Stepped Power Weapon Altar on Island
    sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, SANCTUM_MOTIF_WEAPON_ALTAR, 0, 5, 0, 4, 2, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, SANCTUM_MOTIF_WEAPON_ALTAR, 0, 7, 0, 2, 2, 2, 0.0f, 0.0f, c_cyan, true, false, false);

    // 4 Sarsen Half-Cover Orthostats in the sunken pit
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -10, 2, -10, 2, 4, 3, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  10, 2, -10, 2, 4, 3, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT, -10, 2,  10, 2, 4, 3, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,  10, 2,  10, 2, 4, 3, 0.0f, 0.0f, c_stone, false, false, false);

    // Tactical Spawns & Supplies
    sanctum_add_spawn(plan, (Vector3){ 0.0f, 5.0f, 0.0f }, 0.0f, true); // Player Holdout
    sanctum_add_spawn(plan, (Vector3){   0.0f, 10.0f,  23.0f }, 180.0f, false); // Enemy North
    sanctum_add_spawn(plan, (Vector3){   0.0f, 10.0f, -23.0f },   0.0f, false); // Enemy South
    sanctum_add_spawn(plan, (Vector3){  23.0f, 10.0f,   0.0f }, 270.0f, false); // Enemy East
    sanctum_add_spawn(plan, (Vector3){ -23.0f, 10.0f,   0.0f },  90.0f, false); // Enemy West

    sanctum_add_supply(plan, (Vector3){ 0.0f, 5.0f, 0.0f }, SANCTUM_SUPPLY_HEALTH);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 6.0f, 0.0f }, SANCTUM_SUPPLY_DYNAMIC_SHOT);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 2.0f,  12.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){ 0.0f, 2.0f, -12.0f }, SANCTUM_SUPPLY_AMMO);

    // Four Gravity Jump Lifts launching players from sunken floor onto upper colonnade
    sanctum_add_jump_pad(plan, (Vector3){ -14.0f, 2.0f, -14.0f }, 20.0f);
    sanctum_add_jump_pad(plan, (Vector3){  14.0f, 2.0f, -14.0f }, 20.0f);
    sanctum_add_jump_pad(plan, (Vector3){ -14.0f, 2.0f,  14.0f }, 20.0f);
    sanctum_add_jump_pad(plan, (Vector3){  14.0f, 2.0f,  14.0f }, 20.0f);

    // Outward Sockets
    sanctum_add_socket(plan, (Vector3){ 0.0f, 10.0f, 27.0f }, (Vector3){ 0.0f, 0.0f, 1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_colosseum, 10);
    sanctum_add_socket(plan, (Vector3){ 0.0f, 10.0f, -27.0f }, (Vector3){ 0.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_colosseum, 10);
}

// ---------------------------------------------------------------------------
// Macro-Topology 3: The Asymmetric Outpost (Urban Bunker Citadel)
// ---------------------------------------------------------------------------
static void sanctum_init_asymmetric_outpost(SanctumCitadelPlan *plan) {
    Color c_stone    = (Color){ 105, 102, 95, 255 };  // Sarsen
    Color c_marble   = (Color){ 232, 230, 222, 255 }; // Marble
    Color c_titanium = (Color){ 62, 68, 76, 255 };    // Titanium
    Color c_cyan     = (Color){ 45, 225, 255, 255 };  // Cyan
    // Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Water (disabled for now)
    Color c_gold     = (Color){ 255, 185, 45, 255 };  // Gold

    // ================= SOUTH-WEST SECTOR (Player Holdout Bunker) =================
    int base_sw = sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                                   -12, 0, -12, 28, 4, 28, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     -12, 4, -12, 26, 1, 26, 0.0f, 0.0f, c_marble, false, false, false);

    // Enclosed Megaron Bunker
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     -12, 5, -12, 20, 9, 20, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     -12, 14, -12, 22, 2, 22, 0.0f, 0.0f, c_titanium, false, false, false);

    // Watchtower in SW Corner (rising to Y=28)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     -20, 5, -20, 5, 22, 5, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     -16, 16, -20, 2, 6, 2, 28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_APEX_SPIRE_MATRIX,
                     -20, 27, -20, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);

    // Low Chest-High Parapet Breastworks at bunker exit
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_FIREFIGHT_HOLDOUT,
                     -2, 5, -6, 2, 2, 8, 0.0f, 0.0f, c_marble, false, false, false);

    // Bunker Healing Well (disabled for now)
    // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_RECHARGER_WELL,
    //                  -14, 5, -14, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

    // ================= CENTRAL ZIGZAGGING PERIBOLOS CURTAIN WALL =================
    // Wall segment 1
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIBOLOS_RAMPART,
                     -2, 4, 14, 4, 8, 14, 0.0f, 0.0f, c_stone, false, false, false);
    // Wall segment 2
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIBOLOS_RAMPART,
                     4, 4, 7, 10, 8, 4, 0.0f, 0.0f, c_stone, false, false, false);

    // Fortified Gatehouse Chokepoint (Twin Canted Pylons)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CANTED_PYLON_PORTAL,
                     4, 4, -1, 4, 12, 4, -18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CANTED_PYLON_PORTAL,
                     8, 4, -1, 4, 12, 4,  18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_CANTED_PYLON_PORTAL,
                     6, 14, -1, 10, 2, 4, 0.0f, 0.0f, c_titanium, false, false, false);

    // Power Weapon Altar at Chokepoint
    sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, SANCTUM_MOTIF_WEAPON_ALTAR,
                     6, 4, 0, 4, 2, 4, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, SANCTUM_MOTIF_WEAPON_ALTAR,
                     6, 6, 0, 2, 2, 2, 0.0f, 0.0f, c_cyan, true, false, false);

    // ================= NORTH-EAST SECTOR (Enemy Staging Grounds) =================
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_COURTYARD_PLAZA,
                     14, 4, 14, 26, 2, 26, 0.0f, 0.0f, c_marble, false, false, false);

    // Ruined Peristyle Colonnade in NE Plaza
    for (int k = 0; k < 4; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_COURTYARD_PLAZA,
                         6 + k * 5, 6, 22, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }

    // Staggered Enemy Spawn Crypts
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT,  20, 5, 22, 6, 6, 6, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT,  22, 5,  6, 6, 6, 6, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_SPAWN_CRYPT, -16, 5, 20, 6, 6, 6, 0.0f, 0.0f, c_stone, false, false, false);

    // Tactical Registrations
    sanctum_add_spawn(plan, (Vector3){ -12.0f, 5.0f, -12.0f }, 45.0f, true); // Player Holdout
    sanctum_add_spawn(plan, (Vector3){  20.0f, 5.0f,  22.0f }, 225.0f, false); // Enemy NE
    sanctum_add_spawn(plan, (Vector3){  22.0f, 5.0f,   6.0f }, 270.0f, false); // Enemy East
    sanctum_add_spawn(plan, (Vector3){ -16.0f, 5.0f,  20.0f }, 180.0f, false); // Enemy North Flank

    sanctum_add_supply(plan, (Vector3){ -14.0f, 5.0f, -14.0f }, SANCTUM_SUPPLY_HEALTH);
    sanctum_add_supply(plan, (Vector3){  -8.0f, 5.0f, -16.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){  14.0f, 5.0f,  14.0f }, SANCTUM_SUPPLY_AMMO);
    sanctum_add_supply(plan, (Vector3){   6.0f, 5.0f,   0.0f }, SANCTUM_SUPPLY_VOID);

    // Jump Pad launching up to curtain wall walkway
    sanctum_add_jump_pad(plan, (Vector3){ -4.0f, 4.0f, -2.0f }, 16.0f);

    // Outward Sockets
    sanctum_add_socket(plan, (Vector3){ 24.0f, 4.0f, 14.0f }, (Vector3){ 1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_sw, 4);
}

// ---------------------------------------------------------------------------
// Plan Initialization: Seed-Driven Macro-Topology Selection
// ---------------------------------------------------------------------------
void sanctum_plan_init_ex(SanctumCitadelPlan *plan, uint32_t seed, int forced_topology) {
    if (!plan) return;
    memset(plan, 0, sizeof(SanctumCitadelPlan));
    plan->seed = seed;

    // Seed-Driven Macro-Topology: forced override or seed % SANCTUM_TOPO_COUNT
    if (forced_topology >= 0 && forced_topology < SANCTUM_TOPO_COUNT) {
        plan->topology = (SanctumTopology)forced_topology;
    } else {
        plan->topology = (SanctumTopology)(seed % SANCTUM_TOPO_COUNT);
    }

    switch (plan->topology) {
        case SANCTUM_TOPO_STRONGHOLD:
            sanctum_init_stronghold(plan);
            break;
        case SANCTUM_TOPO_ABYSSAL_RIFT:
            sanctum_init_abyssal_rift(plan);
            break;
        case SANCTUM_TOPO_SUNKEN_CRUCIBLE:
            sanctum_init_sunken_crucible(plan);
            break;
        case SANCTUM_TOPO_ASYMMETRIC_OUTPOST:
            sanctum_init_asymmetric_outpost(plan);
            break;
        default:
            sanctum_init_stronghold(plan);
            break;
    }
}

void sanctum_plan_init(SanctumCitadelPlan *plan, uint32_t seed) {
    sanctum_plan_init_ex(plan, seed, -1);
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
            case SANCTUM_SOCKET_QUADRANT_INFILL:
            case SANCTUM_SOCKET_PLAZA_INFILL: {
                // High-priority infill for diagonal quadrants and open plaza space!
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_OBELISK_PLAZA;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.97f + (float)(h % 10) / 100.0f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_STELAE_AVENUE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.95f + (float)((h >> 2) % 10) / 100.0f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_COURTYARD_PLAZA;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.94f + (float)((h >> 4) % 10) / 100.0f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_HYPOSTYLE_HALL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.92f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_PERIBOLOS_RAMPART;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.88f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_AMMO_CACHE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.91f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_RECHARGER_WELL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.89f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_GRAVITY_LIFT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.87f;
                }
                break;
            }

            case SANCTUM_SOCKET_STOA_CONTINUATION: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_PERIMETER_STOA;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.98f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_PERIBOLOS_RAMPART;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.89f;
                }
                break;
            }

            case SANCTUM_SOCKET_AXIAL_PATH: {
                // Monumental Macro-Structure: Colossal Hypostyle Basilica
                if (count < max_opps && (h % 3 == 0)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_HYPOSTYLE_HALL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.95f;
                }
                // Afford: Doric Colonnade Terrace Avenue
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_DORIC_COLONNADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.90f + (float)((h >> 4) % 15) / 100.0f;
                }
                // Afford: Megalithic Titanium Trilithon Gateway
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_MEGALITHIC_TRILITHON;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.85f + (float)(h % 15) / 100.0f;
                }
                // Afford: Canted Titanium Pylon Portal
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_CANTED_PYLON_PORTAL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.88f + (float)((h >> 8) % 15) / 100.0f;
                }
                // Afford: Abyssal Chasm with Hard-Light Bridge
                if (count < max_opps && plan->count_hardlight < (plan->node_count / 3 + 1)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_HARDLIGHT_CHASM;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.91f + (float)((h >> 12) % 15) / 100.0f;
                }
                // Afford: Sacred PBF Cascade
                if (count < max_opps && (h % 4 == 0)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_SACRED_PBF_CASCADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.86f;
                }
                // Afford: Firefight Holdout Redoubt
                if (count < max_opps && (h % 3 == 1)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_FIREFIGHT_HOLDOUT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.93f;
                }
                // Afford: Enemy Wave Spawn Crypt
                if (count < max_opps && (h % 3 == 2)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_SPAWN_CRYPT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.92f;
                }
                // Afford: Power Weapon Altar
                if (count < max_opps && (h % 5 == 0)) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_WEAPON_ALTAR;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.94f;
                }
                break;
            }

            case SANCTUM_SOCKET_COLONNADE_FLANK: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_COURTYARD_PLAZA;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.92f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_PERIBOLOS_RAMPART;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.86f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_DORIC_COLONNADE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.82f;
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

            case SANCTUM_SOCKET_TACTICAL_CACHE: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_AMMO_CACHE;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.98f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_RECHARGER_WELL;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.95f;
                }
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_WEAPON_ALTAR;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.92f;
                }
                break;
            }

            case SANCTUM_SOCKET_ENEMY_INLET: {
                if (count < max_opps) {
                    SanctumOpportunity *opp = &out_opps[count++];
                    opp->socket_idx = i;
                    opp->motif = SANCTUM_MOTIF_SPAWN_CRYPT;
                    opp->spawn_pos = s->pos;
                    opp->spawn_dir = s->dir;
                    opp->score = 0.99f;
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
    // Color c_water    = (Color){ 50, 160, 220, 255 };  // PBF Fluid (disabled for now)

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

            // Fluted Columns on Left and Right flanks + Vaulted Hard-Light Ribs Overhead
            for (int k = -1; k <= 1; ++k) {
                int cx1 = along_z ? (center_x - 3) : (center_x + k * 4);
                int cz1 = along_z ? (center_z + k * 4) : (center_z - 3);
                int cx2 = along_z ? (center_x + 3) : (center_x + k * 4);
                int cz2 = along_z ? (center_z + k * 4) : (center_z + 3);

                sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                 cx1, py + 1, cz1, 1, col_h, 1, 0.0f, 0.0f, c_marble, false, false, false);
                sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                 cx2, py + 1, cz2, 1, col_h, 1, 0.0f, 0.0f, c_marble, false, false, false);

                // Vaulted Hard-Light Energy Ribs arching between the paired columns
                sanctum_add_node(plan, SANCTUM_PRIM_HARDLIGHT_BRIDGE, opp->motif,
                                 along_z ? center_x : (center_x + k * 4),
                                 py + col_h + 1,
                                 along_z ? (center_z + k * 4) : center_z,
                                 along_z ? 6 : 1, 1, along_z ? 1 : 6,
                                 0.0f, 0.0f, (Color){ 45, 225, 255, 220 }, true, false, false);
            }

            // Continuous Architrave Frieze
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, opp->motif,
                             center_x, py + col_h + 1, center_z, sw, 1, sd, 0.0f, 0.0f, c_marble, false, false, false);

            // Glowing conduit running down the center aisle
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py, center_z, along_z ? 1 : length, 1, along_z ? length : 1, 0.0f, 0.0f, c_cyan, true, false, false);

            // Entrance Bronze Braziers on Left and Right Plinths
            int bx1 = along_z ? (center_x - 4) : px;
            int bz1 = along_z ? pz : (center_z - 4);
            int bx2 = along_z ? (center_x + 4) : px;
            int bz2 = along_z ? pz : (center_z + 4);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             bx1, py + 1, bz1, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             bx2, py + 1, bz2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

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

            // 2. Annular PBF Reflection Fluid Moat around rim (disabled for now)
            // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
            //                  center_x, py, center_z, 14, 1, 14, 0.0f, 0.0f, c_water, false, false, true);

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

            // Dynamic PBF Water Flume (disabled for now)
            // int ww = along_z ? 2 : length;
            // int wd = along_z ? length : 2;
            // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
            //                  center_x, py + 1, center_z, ww, 1, wd, 0.0f, 0.0f, c_water, false, false, true);

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

        case SANCTUM_MOTIF_HYPOSTYLE_HALL: {
            int span = 18;
            int length = 18;
            int hall_h = 9;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // 1. Monumental Stepped Marble Stylobate Floor
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, span, 2, length, 0.0f, 0.0f, c_marble, false, false, false);

            // 2. 4x4 Grid of 16 Colossal Fluted Marble Columns
            for (int r = -2; r <= 1; ++r) {
                for (int c = -2; c <= 1; ++c) {
                    int col_x = center_x + (r * 4 + 2);
                    int col_z = center_z + (c * 4 + 2);
                    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                     col_x, py + 2, col_z, 2, hall_h - 2, 2, 0.0f, 0.0f, c_marble, false, false, false);
                }
            }

            // 3. Heavy Coffered Ceiling Slab in Dark Titanium Alloy
            sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, opp->motif,
                             center_x, py + hall_h, center_z, span, 2, length, 0.0f, 0.0f, c_titanium, false, false, false);

            // 4. Central Clerestory Skylight Void in the ceiling
            sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, opp->motif,
                             center_x, py + hall_h, center_z, 6, 2, 6, 0.0f, 0.0f, BLACK, false, true, false);

            // 5. Four Corner Sarsen Buttress Piers
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - 8, py + 2, center_z - 8, 3, hall_h - 2, 3, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + 8, py + 2, center_z - 8, 3, hall_h - 2, 3, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - 8, py + 2, center_z + 8, 3, hall_h - 2, 3, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + 8, py + 2, center_z + 8, 3, hall_h - 2, 3, 0.0f, 0.0f, c_stone, false, false, false);

            // 6. Cyan Light Frieze around interior entablature
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + hall_h - 1, center_z, span - 2, 1, length - 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // 7. Four Bronze Braziers at interior corners
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             center_x - 6, py + 2, center_z - 6, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             center_x + 6, py + 2, center_z + 6, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

            // Sockets: Forward axial path, lateral flank exits, zenith roof apex
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 2)), (float)py, (float)(pz + step_z * (length + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);

            Vector3 l_dir = along_z ? (Vector3){ -1.0f, 0.0f, 0.0f } : (Vector3){ 0.0f, 0.0f, -1.0f };
            Vector3 r_dir = along_z ? (Vector3){  1.0f, 0.0f, 0.0f } : (Vector3){ 0.0f, 0.0f,  1.0f };
            sanctum_add_socket(plan, (Vector3){ (float)(center_x + l_dir.x * 10), (float)py, (float)(center_z + l_dir.z * 10) },
                               l_dir, SANCTUM_SOCKET_COLONNADE_FLANK, plan->node_count - 1, py);
            sanctum_add_socket(plan, (Vector3){ (float)(center_x + r_dir.x * 10), (float)py, (float)(center_z + r_dir.z * 10) },
                               r_dir, SANCTUM_SOCKET_COLONNADE_FLANK, plan->node_count - 1, py);
            sanctum_add_socket(plan, (Vector3){ (float)center_x, (float)(py + hall_h + 2), (float)center_z },
                               (Vector3){ 0.0f, 1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, plan->node_count - 1, py + hall_h);
            break;
        }

        case SANCTUM_MOTIF_COURTYARD_PLAZA: {
            int size = 16;
            int center_x = px + (int)roundf(dir.x * 8.0f);
            int center_z = pz + (int)roundf(dir.z * 8.0f);

            // 1. Paved Marble Plaza Floor
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, size, 1, size, 0.0f, 0.0f, c_marble, false, false, false);

            // 2. Central Sunken Reflecting Pool with PBF Water (disabled for now)
            // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
            //                  center_x, py + 1, center_z, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

            // 3. Four Corner Sarsen Stelae with Titanium Collars
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - 6, py + 1, center_z - 6, 2, 5, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + 6, py + 1, center_z - 6, 2, 5, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - 6, py + 1, center_z + 6, 2, 5, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + 6, py + 1, center_z + 6, 2, 5, 2, 0.0f, 0.0f, c_stone, false, false, false);

            // 4. Low Decorative Marble Balustrades along outer edges
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z - 7, 12, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z + 7, 12, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

            // 5. Bronze Tripod Braziers at fountain corners
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             center_x - 4, py + 1, center_z - 4, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             center_x + 4, py + 1, center_z + 4, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

            // Sockets: Continue quadrant expansion & flank connections
            sanctum_add_socket(plan,
                               (Vector3){ (float)(center_x + dir.x * 9.0f), (float)py, (float)(center_z + dir.z * 9.0f) },
                               dir, SANCTUM_SOCKET_QUADRANT_INFILL, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_PERIBOLOS_RAMPART: {
            int r_len = 16;
            int r_w = 4;
            int r_h = 6;
            int center_x = px + step_x * (r_len / 2);
            int center_z = pz + step_z * (r_len / 2);

            int rw = along_z ? r_w : r_len;
            int rd = along_z ? r_len : r_w;

            // 1. Cyclopean Stone Curtain Wall
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x, py, center_z, rw, r_h, rd, 0.0f, 0.0f, c_stone, false, false, false);

            // 2. Titanium Parapet Walkway and Crenellations
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py + r_h, center_z, rw, 1, rd, 0.0f, 0.0f, c_titanium, false, false, false);

            // 3. Corner Bastion Watchtower
            int tx = along_z ? center_x : px + step_x * r_len;
            int tz = along_z ? pz + step_z * r_len : center_z;
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             tx, py, tz, 5, r_h + 3, 5, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             tx, py + r_h + 4, tz, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);

            // Socket: continue wall or turn corner
            Vector3 next_dir = along_z ? (Vector3){ 1.0f, 0.0f, 0.0f } : (Vector3){ 0.0f, 0.0f, 1.0f };
            sanctum_add_socket(plan, (Vector3){ (float)tx, (float)py, (float)tz },
                               next_dir, SANCTUM_SOCKET_COLONNADE_FLANK, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_GRAND_STAIRS: {
            int length = 8;
            int width = 12;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            int sw = along_z ? width : length;
            int sd = along_z ? length : width;

            // Stepped Pentelic Marble tiers
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, sw, 2, sd, 0.0f, 0.0f, c_marble, false, false, false);

            // Flanking Guardian Megaliths with glowing light stelae
            int gx1 = along_z ? (center_x - 7) : center_x;
            int gz1 = along_z ? center_z : (center_z - 7);
            int gx2 = along_z ? (center_x + 7) : center_x;
            int gz2 = along_z ? center_z : (center_z + 7);

            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             gx1, py, gz1, 2, 6, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             gx2, py, gz2, 2, 6, 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             gx1, py + 6, gz1, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             gx2, py + 6, gz2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

            // Forward Socket
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (length + 2)), (float)py, (float)(pz + step_z * (length + 2)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_PERIMETER_STOA: {
            int s_len = 16;
            int s_w = 6;
            int s_h = 7;
            int center_x = px + step_x * (s_len / 2);
            int center_z = pz + step_z * (s_len / 2);

            int sw = along_z ? s_w : s_len;
            int sd = along_z ? s_len : s_w;

            // 1. Stepped Marble Stylobate Base
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, sw, 1, sd, 0.0f, 0.0f, c_marble, false, false, false);

            // 2. Back Cyclopean Retaining Wall
            int wx = along_z ? (center_x + (dir.x > 0 ? 2 : -2)) : center_x;
            int wz = along_z ? center_z : (center_z + (dir.z > 0 ? 2 : -2));
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             wx, py + 1, wz, along_z ? 2 : s_len, s_h - 1, along_z ? s_len : 2, 0.0f, 0.0f, c_stone, false, false, false);

            // 3. Colonnade Front (4 Doric Marble Columns)
            int fx = along_z ? (center_x - (dir.x > 0 ? 2 : -2)) : center_x;
            int fz = along_z ? center_z : (center_z - (dir.z > 0 ? 2 : -2));
            for (int k = -2; k <= 1; ++k) {
                int cx = along_z ? fx : (center_x + k * 4 + 2);
                int cz = along_z ? (center_z + k * 4 + 2) : fz;
                sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, opp->motif,
                                 cx, py + 1, cz, 1, s_h - 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
            }

            // 4. Architrave & Coffered Ceiling
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, opp->motif,
                             center_x, py + s_h, center_z, sw, 1, sd, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, opp->motif,
                             center_x, py + s_h + 1, center_z, sw, 1, sd, 0.0f, 0.0f, c_titanium, false, false, false);

            // 5. Emissive Light Channel running down aisle
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + 1, center_z, along_z ? 1 : s_len - 2, 1, along_z ? s_len - 2 : 1, 0.0f, 0.0f, c_cyan, true, false, false);

            // Forward Socket to continue stoa
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (s_len + 1)), (float)py, (float)(pz + step_z * (s_len + 1)) },
                               dir, SANCTUM_SOCKET_STOA_CONTINUATION, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_OBELISK_PLAZA: {
            int p_size = 14;
            int center_x = px + (int)roundf(dir.x * 7.0f);
            int center_z = pz + (int)roundf(dir.z * 7.0f);

            // 1. Paved Marble Stylobate Plaza
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, p_size, 1, p_size, 0.0f, 0.0f, c_marble, false, false, false);

            // 2. Central Sunken Pool with PBF water (disabled for now)
            // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
            //                  center_x, py + 1, center_z, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

            // 3. Central Monolithic Tapered Obelisk
            sanctum_add_node(plan, SANCTUM_PRIM_OBELISK, opp->motif,
                             center_x, py + 1, center_z, 2, 11, 2, 0.0f, 0.0f, c_stone, false, false, false);

            // 4. Four Corner Braziers
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif, center_x - 5, py + 1, center_z - 5, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif, center_x + 5, py + 1, center_z - 5, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif, center_x - 5, py + 1, center_z + 5, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif, center_x + 5, py + 1, center_z + 5, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

            // 5. Low Decorative Balustrades
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z - 6, 10, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z + 6, 10, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

            // Sockets outward
            sanctum_add_socket(plan,
                               (Vector3){ (float)(center_x + dir.x * 8.0f), (float)py, (float)(center_z + dir.z * 8.0f) },
                               dir, SANCTUM_SOCKET_QUADRANT_INFILL, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_STELAE_AVENUE: {
            int a_len = 14;
            int a_w = 8;
            int center_x = px + step_x * (a_len / 2);
            int center_z = pz + step_z * (a_len / 2);

            int aw = along_z ? a_w : a_len;
            int ad = along_z ? a_len : a_w;

            // 1. Paved Marble Avenue Floor
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, aw, 1, ad, 0.0f, 0.0f, c_marble, false, false, false);

            // 2. Central Cyan Light Conduit
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + 1, center_z, along_z ? 1 : a_len, 1, along_z ? a_len : 1, 0.0f, 0.0f, c_cyan, true, false, false);

            // 3. Symmetrical Rows of 3 Megalithic Stelae with cyan glyph inserts
            for (int k = -1; k <= 1; ++k) {
                int sx1 = along_z ? (center_x - 3) : (center_x + k * 4);
                int sz1 = along_z ? (center_z + k * 4) : (center_z - 3);
                int sx2 = along_z ? (center_x + 3) : (center_x + k * 4);
                int sz2 = along_z ? (center_z + k * 4) : (center_z + 3);

                sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                                 sx1, py + 1, sz1, 1, 5, 1, 0.0f, 0.0f, c_stone, false, false, false);
                sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                                 sx1, py + 3, sz1, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);

                sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                                 sx2, py + 1, sz2, 1, 5, 1, 0.0f, 0.0f, c_stone, false, false, false);
                sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                                 sx2, py + 3, sz2, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);
            }

            // 4. Stepped Altar at avenue midpoint
            sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, opp->motif,
                             center_x, py + 1, center_z, 3, 2, 3, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, opp->motif,
                             center_x, py + 3, center_z, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

            // Forward Socket
            sanctum_add_socket(plan,
                               (Vector3){ (float)(px + step_x * (a_len + 1)), (float)py, (float)(pz + step_z * (a_len + 1)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py);
            break;
        }

        case SANCTUM_MOTIF_FIREFIGHT_HOLDOUT: {
            int span = 12;
            int length = 12;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Stepped Platform foundation
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x, py, center_z, span, 2, length, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py + 2, center_z, span - 2, 1, length - 2, 0.0f, 0.0f, c_marble, false, false, false);

            // Perimeter half-cover balustrades (waist-high shooting positions)
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 3, center_z + (length / 2 - 1), span - 4, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 3, center_z - (length / 2 - 1), span - 4, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

            // Corner full-cover titanium bastions
            int cx_off = span / 2 - 2;
            int cz_off = length / 2 - 2;
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif, center_x + cx_off, py + 3, center_z + cz_off, 2, 4, 2, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif, center_x - cx_off, py + 3, center_z + cz_off, 2, 4, 2, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif, center_x + cx_off, py + 3, center_z - cz_off, 2, 4, 2, 0.0f, 0.0f, c_titanium, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif, center_x - cx_off, py + 3, center_z - cz_off, 2, 4, 2, 0.0f, 0.0f, c_titanium, false, false, false);

            // Tactical Resupply Ammo & Health in the holdout
            sanctum_add_supply(plan, (Vector3){ (float)center_x, (float)(py + 3), (float)center_z }, SANCTUM_SUPPLY_AMMO);
            sanctum_add_supply(plan, (Vector3){ (float)(center_x + 2), (float)(py + 3), (float)center_z }, SANCTUM_SUPPLY_HEALTH);

            // Tactical sockets: outward flanking fire sockets
            sanctum_add_socket(plan, (Vector3){ (float)(px + step_x * (length + 1)), (float)(py + 2), (float)(pz + step_z * (length + 1)) },
                               dir, SANCTUM_SOCKET_AXIAL_PATH, plan->node_count - 1, py + 2);
            break;
        }

        case SANCTUM_MOTIF_SPAWN_CRYPT: {
            int span = 8;
            int length = 8;
            int center_x = px + step_x * (length / 2);
            int center_z = pz + step_z * (length / 2);

            // Bedrock chamber floor
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x, py, center_z, span, 1, length, 0.0f, 0.0f, c_stone, false, false, false);

            // Heavy cyclopean side walls
            int side_x = along_z ? 3 : 0;
            int side_z = along_z ? 0 : 3;
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x + side_x, py + 1, center_z + side_z, along_z ? 2 : span - 2, 5, along_z ? length : 2, 0.0f, 0.0f, c_stone, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, opp->motif,
                             center_x - side_x, py + 1, center_z - side_z, along_z ? 2 : span - 2, 5, along_z ? length : 2, 0.0f, 0.0f, c_stone, false, false, false);

            // Massive Megalithic Capstone / Lintel
            sanctum_add_node(plan, SANCTUM_PRIM_STONE_LINTEL, opp->motif,
                             center_x, py + 6, center_z, span + 1, 2, length + 1, 0.0f, 0.0f, c_titanium, false, false, false);

            // Glowing cyan wave emitter / teleporter aperture
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + 1, center_z, 2, 1, 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // Register Enemy Spawn Point inside the crypt chamber facing outward!
            float spawn_yaw = along_z ? (dir.z > 0 ? 0.0f : 180.0f) : (dir.x > 0 ? 90.0f : 270.0f);
            sanctum_add_spawn(plan, (Vector3){ (float)center_x, (float)(py + 1), (float)center_z }, spawn_yaw, false);
            break;
        }

        case SANCTUM_MOTIF_AMMO_CACHE: {
            int center_x = px + step_x * 2;
            int center_z = pz + step_z * 2;

            // Marble plinth
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, 4, 1, 4, 0.0f, 0.0f, c_marble, false, false, false);
            // Titanium locker frame
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py + 1, center_z, 2, 2, 2, 0.0f, 0.0f, c_titanium, false, false, false);
            // Glowing ammo node
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + 2, center_z, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);

            // Register Ammo Supply Point
            sanctum_add_supply(plan, (Vector3){ (float)center_x, (float)(py + 2), (float)center_z }, SANCTUM_SUPPLY_AMMO);
            break;
        }

        case SANCTUM_MOTIF_WEAPON_ALTAR: {
            int center_x = px + step_x * 3;
            int center_z = pz + step_z * 3;

            // Stepped sacrificial altar plinth
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, opp->motif,
                             center_x, py + 1, center_z, 4, 2, 4, 0.0f, 0.0f, c_stone, false, false, false);
            // Levitating Gravity Core / Power Weapon Plinth
            sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, opp->motif,
                             center_x, py + 3, center_z, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);

            // Register High-Value Supply Point (Dynamic Shot or Void)
            SanctumSupplyType ptype = (plan->seed % 2 == 0) ? SANCTUM_SUPPLY_DYNAMIC_SHOT : SANCTUM_SUPPLY_VOID;
            sanctum_add_supply(plan, (Vector3){ (float)center_x, (float)(py + 3), (float)center_z }, ptype);
            break;
        }

        case SANCTUM_MOTIF_RECHARGER_WELL: {
            int center_x = px + step_x * 3;
            int center_z = pz + step_z * 3;

            // Marble basin frame
            sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, opp->motif,
                             center_x, py, center_z, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z + 2, 6, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
            sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, opp->motif,
                             center_x, py + 1, center_z - 2, 6, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);

            // Sacred PBF restorative pool (disabled for now)
            // sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
            //                  center_x, py + 1, center_z, 4, 1, 4, 0.0f, 0.0f, c_water, false, false, true);

            // Register Health Supply Point
            sanctum_add_supply(plan, (Vector3){ (float)center_x, (float)(py + 2), (float)center_z }, SANCTUM_SUPPLY_HEALTH);
            break;
        }

        case SANCTUM_MOTIF_GRAVITY_LIFT: {
            int center_x = px + step_x * 2;
            int center_z = pz + step_z * 2;

            // Octagonal / square titanium jump pad base
            sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, opp->motif,
                             center_x, py, center_z, 4, 1, 4, 0.0f, 0.0f, c_titanium, false, false, false);
            // Glowing cyan gravitational lift emitter
            sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, opp->motif,
                             center_x, py + 1, center_z, 2, 1, 2, 0.0f, 0.0f, c_cyan, true, false, false);

            // Register Jump Pad
            sanctum_add_jump_pad(plan, (Vector3){ (float)center_x, (float)(py + 1), (float)center_z }, 18.0f);
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
SanctumCitadelPlan generate_unified_sanctum_ex(uint32_t seed, int growth_steps, int forced_topology) {
    SanctumCitadelPlan plan;
    sanctum_plan_init_ex(&plan, seed, forced_topology);
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

SanctumCitadelPlan generate_unified_sanctum(uint32_t seed, int growth_steps) {
    return generate_unified_sanctum_ex(seed, growth_steps, -1);
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

    // Invariant 4: Firefight Combat Invariants
    if (plan->spawn_point_count < 3) {
        printf("Verification failed: Firefight map has insufficient spawns (%d < 3)!\n", plan->spawn_point_count);
        return false;
    }
    bool has_player = false;
    int enemy_count = 0;
    for (int i = 0; i < plan->spawn_point_count; ++i) {
        if (plan->spawn_points[i].is_player) has_player = true;
        else enemy_count++;
    }
    if (!has_player) {
        printf("Verification failed: Firefight map is missing player holdout spawn!\n");
        return false;
    }
    if (enemy_count < 2) {
        printf("Verification failed: Firefight map has insufficient enemy wave ingress points (%d < 2)!\n", enemy_count);
        return false;
    }
    if (plan->supply_count < 2) {
        printf("Verification failed: Firefight map has insufficient tactical supply points (%d < 2)!\n", plan->supply_count);
        return false;
    }

    return true;
}
