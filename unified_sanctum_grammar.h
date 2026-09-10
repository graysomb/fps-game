#ifndef UNIFIED_SANCTUM_GRAMMAR_H
#define UNIFIED_SANCTUM_GRAMMAR_H

#include "raylib.h"
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_SANCTUM_NODES 2048
#define MAX_SANCTUM_SOCKETS 512
#define MAX_SANCTUM_OPPS 96

// ---------------------------------------------------------------------------
// Architectural Socket & Frontier Typology
// ---------------------------------------------------------------------------
typedef enum {
    SANCTUM_SOCKET_AXIAL_PATH,       // North/South/East/West avenue or skybridge
    SANCTUM_SOCKET_COLONNADE_FLANK,  // Lateral extension for peristyles and pylons
    SANCTUM_SOCKET_PORTAL_GATE,      // Gateways, trilithons, and chevron portals
    SANCTUM_SOCKET_ROTUNDA_RADIAL,   // Circular perimeter connections (Tholos/Stone Circle)
    SANCTUM_SOCKET_VERTICAL_APEX,    // Upward (spires, pediments) or downward (crypts, pits)
    SANCTUM_SOCKET_HYDRAULIC_GUTTER, // Fluid aqueduct and coolant channels
    SANCTUM_SOCKET_QUADRANT_INFILL,  // Diagonal quadrant corner (plaza, hypostyle, rampart)
    SANCTUM_SOCKET_STOA_CONTINUATION,// Perimeter stoa colonnade run
    SANCTUM_SOCKET_PLAZA_INFILL      // Courtyard plaza detail cluster
} SanctumSocketType;

typedef struct {
    Vector3 pos;             // Grid center (x, y, z)
    Vector3 dir;             // Outward normal direction vector
    SanctumSocketType type;   // Connection type
    int parent_node_idx;     // Index of generating node
    int elevation;           // Grid height Y
    bool is_occupied;        // True if connected to a continuation motif
} SanctumSocket;

// ---------------------------------------------------------------------------
// Syncretic Primitive Categories
// ---------------------------------------------------------------------------
typedef enum {
    SANCTUM_PRIM_STONE_ORTHOSTAT,    // Weathered cyclopean sarsen monolith
    SANCTUM_PRIM_STONE_LINTEL,       // Massive megalithic capstone
    SANCTUM_PRIM_MARBLE_STYLOBATE,   // Stepped Pentelic white marble platform
    SANCTUM_PRIM_MARBLE_COLUMN,      // Fluted Doric/Ionic column
    SANCTUM_PRIM_MARBLE_ARCHITRAVE,  // Classical entablature & frieze
    SANCTUM_PRIM_MARBLE_PEDIMENT,    // Triangular tympanum with gold accents
    SANCTUM_PRIM_TITANIUM_PYLON,     // Canted metallic alloy bastion
    SANCTUM_PRIM_TITANIUM_BUTTRESS,  // Sloped exterior flying buttress
    SANCTUM_PRIM_HARDLIGHT_BRIDGE,   // Cyan luminescent forcefield bridge
    SANCTUM_PRIM_LIGHT_CHANNEL,      // Glowing emissive circuit groove
    SANCTUM_PRIM_GRAVITY_CORE,       // Levitating octahedral matrix
    SANCTUM_PRIM_PBF_WATER,          // Physical PBF fluid spring/channel
    SANCTUM_PRIM_CHASM_VOID,         // Bedrock excavation void
    SANCTUM_PRIM_COFFERED_CEILING,   // Vaulted coffered roof slab
    SANCTUM_PRIM_BALUSTRADE,         // Low decorative marble parapet
    SANCTUM_PRIM_BRAZIER,            // Sacred bronze flame / telemetry beacon
    SANCTUM_PRIM_OBELISK,            // Monolithic tapered needle stela
    SANCTUM_PRIM_ALTAR               // Stepped sacrificial hearth plinth
} SanctumPrimitiveType;

// ---------------------------------------------------------------------------
// Syncretic Architectural Motifs
// ---------------------------------------------------------------------------
typedef enum {
    SANCTUM_MOTIF_CORE_NEXUS,            // Central hybrid altar & tholos
    SANCTUM_MOTIF_ACROPOLIS_PODIUM,      // Colossal stepped bedrock foundation terrace
    SANCTUM_MOTIF_HYPOSTYLE_HALL,        // Enclosed 4x4 monumental column basilica
    SANCTUM_MOTIF_COURTYARD_PLAZA,       // Paved peristyle quadrangle with fountain
    SANCTUM_MOTIF_PERIBOLOS_RAMPART,     // Fortified cyclopean curtain wall & battlements
    SANCTUM_MOTIF_GRAND_STAIRS,          // Wide monumental stepped flight
    SANCTUM_MOTIF_MEGALITHIC_TRILITHON,  // Titanium-clamped megalithic gateway
    SANCTUM_MOTIF_DORIC_COLONNADE,       // Fluted marble colonnade with light frieze
    SANCTUM_MOTIF_CANTED_PYLON_PORTAL,   // Canted twin pylons with chevron span
    SANCTUM_MOTIF_THOLOS_GRAVITY_PIT,    // Sunken circular rotunda over gravity well
    SANCTUM_MOTIF_HARDLIGHT_CHASM,       // Abyssal canyon with hard-light crossing
    SANCTUM_MOTIF_CORBELLED_CRYPT,       // Subterranean cyclopean vault
    SANCTUM_MOTIF_APEX_SPIRE_MATRIX,     // Towering pinnacle with levitating core
    SANCTUM_MOTIF_SACRED_PBF_CASCADE,    // Stepped water cascade & coolant canal
    SANCTUM_MOTIF_MONUMENTAL_CITADEL,    // Colossal 3-tiered Ziggurat-Megaron Citadel (Y=0..36)
    SANCTUM_MOTIF_PERIMETER_STOA,        // Continuous colonnaded covered portico
    SANCTUM_MOTIF_STELAE_AVENUE,         // Processional double row of megalithic stelae
    SANCTUM_MOTIF_OBELISK_PLAZA,         // Tapered needle obelisk plaza with braziers
    SANCTUM_MOTIF_COUNT
} SanctumMotif;

// ---------------------------------------------------------------------------
// Structural Node Definition
// ---------------------------------------------------------------------------
typedef struct {
    SanctumPrimitiveType prim_type;
    SanctumMotif motif;
    int x, y, z;             // Voxel grid coordinates
    int w, h, d;             // Voxel volume extent (width, height, depth)
    float cant_angle_x;      // Slope angle in X (degrees)
    float cant_angle_z;      // Slope angle in Z (degrees)
    Color color;             // Material tint
    bool is_emissive;        // Glows in dark / hard-light
    bool is_void;            // True if excavates existing voxels
    bool is_fluid;           // Spawns active physical PBF fluid voxels
} SanctumNode;

// ---------------------------------------------------------------------------
// Generative Opportunity (Coalgebra Continuation)
// ---------------------------------------------------------------------------
typedef struct {
    int socket_idx;
    SanctumMotif motif;
    Vector3 spawn_pos;
    Vector3 spawn_dir;
    float score;
} SanctumOpportunity;

// ---------------------------------------------------------------------------
// Citadel Plan (Complete Structural State)
// ---------------------------------------------------------------------------
typedef struct {
    uint32_t seed;
    int growth_steps;
    int node_count;
    int socket_count;
    SanctumNode nodes[MAX_SANCTUM_NODES];
    SanctumSocket sockets[MAX_SANCTUM_SOCKETS];

    // Bounding metrics
    int min_x, max_x;
    int min_y, max_y;
    int min_z, max_z;

    // Architectural census
    int count_stone;
    int count_marble;
    int count_titanium;
    int count_hardlight;
    int count_fluid;
    int count_void;
} SanctumCitadelPlan;

// ---------------------------------------------------------------------------
// Bi-Algebra API
// ---------------------------------------------------------------------------

// Initialize citadel plan with central core nexus
void sanctum_plan_init(SanctumCitadelPlan *plan, uint32_t seed);

// Coalgebra γ(S, F): Inspect open sockets and afford valid continuation motifs
int sanctum_coalgebra_frontier(const SanctumCitadelPlan *plan, SanctumOpportunity *out_opps, int max_opps);

// Algebra α(S, F, opp): Snap chosen motif to socket, spawn nodes, register new sockets
bool sanctum_algebra_expand(SanctumCitadelPlan *plan, const SanctumOpportunity *opp);

// Complete open-ended generation loop for K growth steps
SanctumCitadelPlan generate_unified_sanctum(uint32_t seed, int growth_steps);

// Architectural invariant verifier
bool verify_sanctum_invariants(const SanctumCitadelPlan *plan);

#ifdef __cplusplus
}
#endif

#endif // UNIFIED_SANCTUM_GRAMMAR_H
