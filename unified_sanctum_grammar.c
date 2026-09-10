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
// Plan Initialization: Monumental Acropolis Foundation & Central Sanctuary
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
    Color c_gold     = (Color){ 255, 185, 45, 255 };  // Sacred brazier fire / electrum

    // =======================================================================
    // 1. THE MONUMENTAL ACROPOLIS PODIUM (Multi-Tiered Grounding Fortress)
    // =======================================================================
    // Tier 0 Lower Cyclopean Bedrock Base: 52 x 4 x 52 megalith mass at Y=0
    int base_podium = sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                                       0, 0, 0, 52, 4, 52, 0.0f, 0.0f, c_stone, false, false, false);

    // Four Massive Fortified Pewter Titanium Corner Bastions: 8 x 8 x 8 at Y=0..8
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                     -22, 0, -22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                      22, 0, -22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                     -22, 0,  22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                      22, 0,  22, 8, 8, 8, 0.0f, 0.0f, c_titanium, false, false, false);

    // Bastion Beacon Braziers atop corner bastions at Y=8
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 8, -22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 8, -22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM, -22, 8,  22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_ACROPOLIS_PODIUM,  22, 8,  22, 2, 2, 2, 0.0f, 0.0f, c_gold, true, false, false);

    // Tier 1 Grand Stylobate Terrace: 46 x 2 x 46 polished Pentelic marble at Y=4
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_ACROPOLIS_PODIUM,
                     0, 4, 0, 46, 2, 46, 0.0f, 0.0f, c_marble, false, false, false);

    // Four Cardinal Processional Staircases descending from terrace Y=4 to ground Y=0
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  0, 1,  25, 14, 3, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  0, 1, -25, 14, 3, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS,  25, 1,  0, 6, 3, 14, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_GRAND_STAIRS, -25, 1,  0, 6, 3, 14, 0.0f, 0.0f, c_marble, false, false, false);

    // =======================================================================
    // 2. CONTINUOUS PERIMETER ENCLOSURES (Grand Doric Stoas & Curtain Walls)
    // =======================================================================
    // North Perimeter Stoa: Colonnaded covered gallery along Z = 19, Y = 6..13
    // Back Cyclopean Wall
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIMETER_STOA,
                     0, 6, 21, 38, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    // Colonnade Front (8 Doric Columns)
    for (int k = -3; k <= 3; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_PERIMETER_STOA,
                         k * 5, 6, 17, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }
    // Stoa Architrave & Coffered Roof
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_ARCHITRAVE, SANCTUM_MOTIF_PERIMETER_STOA,
                     0, 12, 19, 38, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_PERIMETER_STOA,
                     0, 13, 19, 38, 1, 6, 0.0f, 0.0f, c_titanium, false, false, false);

    // South Perimeter Stoa Wings (Flanking the Grand Propylaea Avenue):
    // West Wing (X: -19 to -7)
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIMETER_STOA,
                     -13, 6, -21, 14, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    for (int k = 0; k < 3; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_PERIMETER_STOA,
                         -17 + k * 4, 6, -17, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_PERIMETER_STOA,
                     -13, 12, -19, 14, 2, 6, 0.0f, 0.0f, c_titanium, false, false, false);

    // East Wing (X: 7 to 19)
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_PERIMETER_STOA,
                      13, 6, -21, 14, 7, 2, 0.0f, 0.0f, c_stone, false, false, false);
    for (int k = 0; k < 3; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_PERIMETER_STOA,
                          9 + k * 4, 6, -17, 1, 6, 1, 0.0f, 0.0f, c_marble, false, false, false);
    }
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_PERIMETER_STOA,
                      13, 12, -19, 14, 2, 6, 0.0f, 0.0f, c_titanium, false, false, false);

    // =======================================================================
    // 3. THE COLOSSAL MEGARON MEGASTRUCTURE (Commanding Ziggurat-Citadel, Y=6..38)
    // =======================================================================
    // A. Raised Temple Stylobate (28 x 2 x 28 at Y = 6)
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 6, 0, 28, 2, 28, 0.0f, 0.0f, c_marble, false, false, false);

    // B. Enclosed Cella Outer Walls: Cyclopean Stone Masses (Y=8..20, Height 12)
    // North Back Wall (26 x 12 x 3 at Z = 11)
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 8, 11, 26, 12, 3, 0.0f, 0.0f, c_stone, false, false, false);
    // West Flank Wall (3 x 12 x 22 at X = -12)
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     -12, 8, 0, 3, 12, 22, 0.0f, 0.0f, c_stone, false, false, false);
    // East Flank Wall (3 x 12 x 22 at X = 12)
    sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                      12, 8, 0, 3, 12, 22, 0.0f, 0.0f, c_stone, false, false, false);

    // Engaged Fluted Doric Columns along East and West Exterior Walls
    for (int k = -2; k <= 2; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                         -14, 8, k * 5, 2, 12, 2, 0.0f, 0.0f, c_marble, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                          14, 8, k * 5, 2, 12, 2, 0.0f, 0.0f, c_marble, false, false, false);
    }

    // C. South Facade Monumental Propylaea Gatehouse (Grand Entrance)
    // Twin Canted Titanium Pylons (canted at 18 degrees)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     -7, 8, -12, 4, 14, 4, -18.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                      7, 8, -12, 4, 14, 4,  18.0f, 0.0f, c_titanium, false, false, false);
    // Monumental Inscribed Titanium Lintel Beam spanning the portal at Y=20
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 20, -12, 18, 3, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    // Cyan Emissive Glyph Seam across the entrance lintel
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 21, -12, 14, 1, 2, 0.0f, 0.0f, c_cyan, true, false, false);

    // D. Interior Great Hall (Cella Interior)
    // Two rows of 4 fluted marble columns lining the central nave
    for (int k = -1; k <= 2; ++k) {
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                         -5, 8, k * 4 - 2, 2, 10, 2, 0.0f, 0.0f, c_marble, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_COLUMN, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                          5, 8, k * 4 - 2, 2, 10, 2, 0.0f, 0.0f, c_marble, false, false, false);
    }
    // Sacred Central Reflecting Basin
    sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 8, 2, 6, 1, 10, 0.0f, 0.0f, c_water, false, false, true);
    // Glowing cyan floor energy conduit
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 8, 2, 1, 1, 14, 0.0f, 0.0f, c_cyan, true, false, false);
    // Bronze altar braziers in the sanctuary
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_MONUMENTAL_CITADEL, -3, 8, 8, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_MONUMENTAL_CITADEL,  3, 8, 8, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

    // E. Heavy Coffered Ceiling Slab & Clerestory Skylight at Y=20
    sanctum_add_node(plan, SANCTUM_PRIM_COFFERED_CEILING, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 20, 0, 26, 2, 26, 0.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_CHASM_VOID, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 20, 2, 8, 2, 8, 0.0f, 0.0f, BLACK, false, true, false);

    // F. Upper Sky Cella & Cantilevered Observation Gallery (Y=22..26)
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 22, 0, 20, 2, 20, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 24, -9, 18, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 24,  9, 18, 1, 1, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     -9, 24, 0, 1, 1, 18, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BALUSTRADE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                      9, 24, 0, 1, 1, 18, 0.0f, 0.0f, c_marble, false, false, false);

    // Upper Stepped Pyramidal Tier (14 x 2 x 14 at Y=24)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 24, 0, 14, 2, 14, 0.0f, 0.0f, c_titanium, false, false, false);

    // G. Apex Transmission Spire Needle & Canted Flying Buttresses (Y=26..38)
    // Soaring Central Titanium Needle (4 x 12 x 4)
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_PYLON, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 26, 0, 4, 12, 4, 0.0f, 0.0f, c_titanium, false, false, false);
    // Four Canted Flying Buttresses anchoring the needle to the upper terrace
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     -4, 24, 0, 2, 6, 2, -28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                      4, 24, 0, 2, 6, 2,  28.0f, 0.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 24, -4, 2, 6, 2, 0.0f, -28.0f, c_titanium, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_TITANIUM_BUTTRESS, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 24,  4, 2, 6, 2, 0.0f,  28.0f, c_titanium, false, false, false);
    // Levitating Promethean Oracle Gravity Core at Y=35
    sanctum_add_node(plan, SANCTUM_PRIM_GRAVITY_CORE, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 35, 0, 3, 3, 3, 0.0f, 0.0f, c_gold, true, false, false);
    // Luminescent apex conduit ring
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_MONUMENTAL_CITADEL,
                     0, 34, 0, 6, 1, 6, 0.0f, 0.0f, c_cyan, true, false, false);

    // =======================================================================
    // 4. DENSE DETAIL INFILL STRUCTURES (Filling Intermediate Spaces)
    // =======================================================================
    // A. South Processional Stelae Avenue: Rows of orthostats with cyan glyphs
    for (int s = 0; s < 4; ++s) {
        int sz = -15 - s * 3;
        sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_STELAE_AVENUE,
                         -5, 6, sz, 1, 4, 1, 0.0f, 0.0f, c_stone, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_STELAE_AVENUE,
                         -5, 8, sz, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_STONE_ORTHOSTAT, SANCTUM_MOTIF_STELAE_AVENUE,
                          5, 6, sz, 1, 4, 1, 0.0f, 0.0f, c_stone, false, false, false);
        sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_STELAE_AVENUE,
                          5, 8, sz, 1, 1, 1, 0.0f, 0.0f, c_cyan, true, false, false);
    }
    // Processional light conduit along south avenue
    sanctum_add_node(plan, SANCTUM_PRIM_LIGHT_CHANNEL, SANCTUM_MOTIF_STELAE_AVENUE,
                     0, 6, -20, 1, 1, 12, 0.0f, 0.0f, c_cyan, true, false, false);

    // B. Monumental Sacrificial Hearth / Altar Plinth at South Approach (0, 6, -11)
    sanctum_add_node(plan, SANCTUM_PRIM_ALTAR, SANCTUM_MOTIF_STELAE_AVENUE,
                     0, 6, -11, 4, 2, 4, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_STELAE_AVENUE,
                     0, 8, -11, 2, 1, 2, 0.0f, 0.0f, c_gold, true, false, false);

    // C. West Monumental Obelisk Plaza (X = -17, Z = 0)
    // Stepped Marble Base
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_OBELISK_PLAZA,
                     -17, 6, 0, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    // Tapered Sarsen Obelisk rising 12 voxels to Y=19
    sanctum_add_node(plan, SANCTUM_PRIM_OBELISK, SANCTUM_MOTIF_OBELISK_PLAZA,
                     -17, 7, 0, 2, 12, 2, 0.0f, 0.0f, c_stone, false, false, false);
    // Four Perimeter Bronze Braziers
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA, -19, 7, -2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA, -15, 7, -2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA, -19, 7,  2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA, -15, 7,  2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

    // D. East Monumental Obelisk Plaza (X = 17, Z = 0)
    sanctum_add_node(plan, SANCTUM_PRIM_MARBLE_STYLOBATE, SANCTUM_MOTIF_OBELISK_PLAZA,
                      17, 6, 0, 6, 1, 6, 0.0f, 0.0f, c_marble, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_OBELISK, SANCTUM_MOTIF_OBELISK_PLAZA,
                      17, 7, 0, 2, 12, 2, 0.0f, 0.0f, c_stone, false, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA,  15, 7, -2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA,  19, 7, -2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA,  15, 7,  2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);
    sanctum_add_node(plan, SANCTUM_PRIM_BRAZIER, SANCTUM_MOTIF_OBELISK_PLAZA,  19, 7,  2, 1, 1, 1, 0.0f, 0.0f, c_gold, true, false, false);

    // E. Hydraulic Water Rills connecting the central megaron to the obelisk basins
    sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_OBELISK_PLAZA,
                     -15, 6, 0, 4, 1, 2, 0.0f, 0.0f, c_water, false, false, true);
    sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, SANCTUM_MOTIF_OBELISK_PLAZA,
                      15, 6, 0, 4, 1, 2, 0.0f, 0.0f, c_water, false, false, true);

    // =======================================================================
    // 5. EXPANSION FRONTIER SOCKETS
    // =======================================================================
    // Cardinal Processional Avenues extending outward at ground level
    sanctum_add_socket(plan, (Vector3){   0.0f, 0.0f,  29.0f }, (Vector3){  0.0f, 0.0f,  1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);
    sanctum_add_socket(plan, (Vector3){   0.0f, 0.0f, -29.0f }, (Vector3){  0.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);
    sanctum_add_socket(plan, (Vector3){  29.0f, 0.0f,   0.0f }, (Vector3){  1.0f, 0.0f,  0.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);
    sanctum_add_socket(plan, (Vector3){ -29.0f, 0.0f,   0.0f }, (Vector3){ -1.0f, 0.0f,  0.0f }, SANCTUM_SOCKET_AXIAL_PATH, base_podium, 0);

    // Terrace Quadrant Infill Sockets (North-East, North-West, South-East, South-West)
    sanctum_add_socket(plan, (Vector3){  16.0f, 6.0f,  12.0f }, (Vector3){  1.0f, 0.0f,  1.0f }, SANCTUM_SOCKET_QUADRANT_INFILL, base_podium, 6);
    sanctum_add_socket(plan, (Vector3){ -16.0f, 6.0f,  12.0f }, (Vector3){ -1.0f, 0.0f,  1.0f }, SANCTUM_SOCKET_QUADRANT_INFILL, base_podium, 6);
    sanctum_add_socket(plan, (Vector3){  16.0f, 6.0f, -12.0f }, (Vector3){  1.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_QUADRANT_INFILL, base_podium, 6);
    sanctum_add_socket(plan, (Vector3){ -16.0f, 6.0f, -12.0f }, (Vector3){ -1.0f, 0.0f, -1.0f }, SANCTUM_SOCKET_QUADRANT_INFILL, base_podium, 6);

    // Stoa Continuation Sockets (East and West ends of stoas)
    sanctum_add_socket(plan, (Vector3){  20.0f, 6.0f,  19.0f }, (Vector3){  1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_STOA_CONTINUATION, base_podium, 6);
    sanctum_add_socket(plan, (Vector3){ -20.0f, 6.0f,  19.0f }, (Vector3){ -1.0f, 0.0f, 0.0f }, SANCTUM_SOCKET_STOA_CONTINUATION, base_podium, 6);

    // Vertical Apex Sockets (Apex Needle Zenith and Subterranean Crypt)
    sanctum_add_socket(plan, (Vector3){ 0.0f, 38.0f, 0.0f }, (Vector3){ 0.0f,  1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, base_podium, 38);
    sanctum_add_socket(plan, (Vector3){ 0.0f, -1.0f, 0.0f }, (Vector3){ 0.0f, -1.0f, 0.0f }, SANCTUM_SOCKET_VERTICAL_APEX, base_podium, -1);
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

            // 2. Central Sunken Reflecting Pool with PBF Water
            sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
                             center_x, py + 1, center_z, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

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

            // 2. Central Sunken Pool with PBF water
            sanctum_add_node(plan, SANCTUM_PRIM_PBF_WATER, opp->motif,
                             center_x, py + 1, center_z, 6, 1, 6, 0.0f, 0.0f, c_water, false, false, true);

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
