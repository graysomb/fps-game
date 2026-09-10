#ifndef MEGALITH_GRAMMAR_H
#define MEGALITH_GRAMMAR_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Coarse geological & structural primitive types
typedef enum {
    MEGALITH_PRIMITIVE_NONE = 0,
    MEGALITH_PRIMITIVE_ORTHOSTAT,   // Upright standing megalith / menhir (bearing pillar)
    MEGALITH_PRIMITIVE_CAPSTONE,    // Massive horizontal lintel / table slab
    MEGALITH_PRIMITIVE_WALL_SLAB,   // Drystone / upright revetment slab for chambers
    MEGALITH_PRIMITIVE_CORBEL,      // Inward-stepping cantilevered vault slab
    MEGALITH_PRIMITIVE_CAIRN,       // Packed stone rubble fill / burial cairn
    MEGALITH_PRIMITIVE_MOUND,       // Earth / turf mantle covering chamber and passage
    MEGALITH_PRIMITIVE_FLOOR_SLAB,  // Paved flagstone passage / threshold
    MEGALITH_PRIMITIVE_HEARTH       // Central ritual fire pit / offering stone
} MegalithPrimitiveType;

// High-level architectural archetypes
typedef enum {
    MEGALITH_ARCHETYPE_PASSAGE_GRAVE = 0, // Sepulchral: Dolmen -> Chamber -> Passage -> Tumulus
    MEGALITH_ARCHETYPE_STONE_CIRCLE        // Astronomical: Menhir -> Arc -> Cromlech -> Henge -> Avenue
} MegalithArchetype;

// Developmental maturation stages
typedef enum {
    MEGALITH_STAGE_MENHIR = 0,        // Solitary standing stone
    MEGALITH_STAGE_DOLMEN,            // 2-3 orthostats + massive capstone (portal tomb)
    MEGALITH_STAGE_CHAMBER,           // Polygonal enclosed vault with corbelled roof
    MEGALITH_STAGE_PASSAGE_GRAVE,     // Elongated solstice passage + chamber
    MEGALITH_STAGE_TUMULUS,           // Earthen barrow mantle with peristalith kerb ring
    MEGALITH_STAGE_HENGE              // Concentric megalithic ring + trilithons + avenue
} MegalithStage;

// Rules available in the prehistoric grammar
typedef enum {
    MEGALITH_RULE_NONE = 0,
    MEGALITH_RULE_ERECT_ORTHOSTAT,    // Place standing stone on terrain
    MEGALITH_RULE_PAIR_ORTHOSTATS,     // Erect opposing bearing stone to form portal
    MEGALITH_RULE_BRIDGE_CAPSTONE,     // Place massive slab atop bearing supports
    MEGALITH_RULE_ENCLOSE_CHAMBER,     // Curve orthostats into polygonal vault
    MEGALITH_RULE_CORBEL_ROOF,         // Step corbel slabs inward to dome chamber
    MEGALITH_RULE_EXTEND_PASSAGE,      // Add paired orthostats along astronomical axis
    MEGALITH_RULE_EXPAND_CIRCLE,       // Add standing stone along circular arc
    MEGALITH_RULE_RAISE_TRILITHON,     // Construct inner monumental trilithon horseshoe
    MEGALITH_RULE_ALIGN_AVENUE,        // Extend parallel megalithic processional rows
    MEGALITH_RULE_ACCUMULATE_MOUND     // Mound earth, turf, and kerb stones over vault
} MegalithRule;

// 3D Box in modular geological units (1 unit ~ 1 meter)
typedef struct {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
} MegaBox3D;

// Semantic stone node placed in the structure
typedef struct {
    int node_id;
    MegalithPrimitiveType type;
    MegaBox3D box;             // Extents in modular units
    int support_parents[4];   // Node IDs of stones supporting this node
    int support_count;
    float tilt_x;             // Weathered natural lean (-0.2 to +0.2)
    float tilt_z;
    float roughness;          // Surface noise & erosion intensity
    uint32_t flags;
} MegalithNode;

// Candidate developmental opportunity emitted by coalgebra gamma(S)
typedef struct {
    MegalithRule rule;
    MegaBox3D target_box;
    int support_parents[4];
    int support_count;
    float stability_score;    // Center-of-mass balance & contact area
    float alignment_score;    // Solstice / astronomical axis bias
    float enclosure_score;    // Negative-space boundary closure
    float total_score;
} MegalithOpportunity;

#define MAX_MEGALITH_NODES 512
#define MAX_MEGALITH_OPPORTUNITIES 64

// Full Megalithic Building State S
typedef struct {
    MegalithNode nodes[MAX_MEGALITH_NODES];
    int node_count;

    MegalithArchetype archetype;
    MegalithStage stage;
    int target_mass_nodes;    // Stopping condition: accumulated stone mass
    uint32_t seed;

    // Astronomical alignment vector (e.g. Winter Solstice sunrise azimuth)
    float axis_dx;
    float axis_dz;

    // Chamber center and dimensions
    int chamber_cx;
    int chamber_cz;
    int chamber_radius;
    int passage_length;
    int circle_radius;
    int circle_stones_placed;
    int circle_stones_total;

    bool has_capstone;
    bool has_corbel;
    bool has_passage;
    bool has_mound;
    bool has_hearth;
    bool has_avenue;
} MegalithPlan;

// --- Coalgebra (gamma): Inspects S for support, void, boundary, axis, and terrain ---
int megalith_coalgebra_inspect(const MegalithPlan *plan, MegalithOpportunity *out_opps, int max_opps);

// --- Algebra (alpha): Places stone nodes according to winning opportunity ---
bool megalith_algebra_apply(MegalithPlan *plan, const MegalithOpportunity *opp);

// --- Developmental Stepper: Runs developmental growth loop from seed to target mass ---
MegalithPlan generate_megalith_structure(uint32_t seed, MegalithArchetype archetype, MegalithStage target_stage);

// --- Diagnostics & Blueprint Rendering ---
void print_megalith_blueprint(const MegalithPlan *plan);

#ifdef __cplusplus
}
#endif

#endif // MEGALITH_GRAMMAR_H
