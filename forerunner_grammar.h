#ifndef FORERUNNER_GRAMMAR_H
#define FORERUNNER_GRAMMAR_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_FORERUNNER_NODES 256
#define MAX_FORERUNNER_PARENTS 8
#define MAX_FORERUNNER_OPPORTUNITIES 32

// Architectural archetypes representing distinct Halo Forerunner typologies
typedef enum {
    FORERUNNER_ARCHETYPE_CARTOGRAPHER = 0, // Monumental axial chasm, canted chevron gateway & floating core
    FORERUNNER_ARCHETYPE_CROSSROADS,       // Intersecting cruciform chasm with bi-level overlapping cross-bridges
    FORERUNNER_ARCHETYPE_CRUCIBLE,         // Radial ring of inward-canted pylons around a subterranean gravity pit
    FORERUNNER_ARCHETYPE_SPIRE,            // Towering central beam pinnacle with quad flying buttress cantilevers
    FORERUNNER_ARCHETYPE_COUNT
} ForerunnerArchetype;

// Architectural materials for Forerunner brutalist megastructures
typedef enum {
    FORERUNNER_MAT_PEWTER = 0,       // Brushed gunmetal / titanium alloy hull plating
    FORERUNNER_MAT_DARK_BASALT,      // Subterranean trench bedrock / trench floor
    FORERUNNER_MAT_HARDLIGHT_CYAN,   // Emissive solid-light energy bridge & conduit
    FORERUNNER_MAT_SOLAR_AMBER,      // Telemetry beam & power coupling
    FORERUNNER_MAT_BRONZE_ACCENT     // Warm metallic alloy on glyph panels
} ForerunnerMaterial;

// Architectural primitives in the Forerunner vocabulary
typedef enum {
    FORERUNNER_PRIM_NONE = 0,
    FORERUNNER_PRIM_CHASM,           // Intentional negative void excavation
    FORERUNNER_PRIM_PYLON,           // Canted inward-leaning structural buttress
    FORERUNNER_PRIM_LINTEL,          // Heavy angled chevron lintel beam / collar
    FORERUNNER_PRIM_BRIDGE_SPAN,     // Cantilevered solid metal deck
    FORERUNNER_PRIM_HARDLIGHT_BRIDGE,// Pure energy crossing over abyssal void
    FORERUNNER_PRIM_SPIRE,           // Towering vertical energy conduit spire
    FORERUNNER_PRIM_LIGHT_CHANNEL,   // Recessed glowing emissive conduit groove
    FORERUNNER_PRIM_GRAVITY_CORE,    // Hovering central installation matrix
    FORERUNNER_PRIM_TERRACE          // Stepped side gallery / wall tier
} ForerunnerPrimitiveType;

// Developmental maturation stages
typedef enum {
    FORERUNNER_STAGE_FOUNDATION = 0, // Abyssal excavation (chasm, cruciform rift, or gravity pit)
    FORERUNNER_STAGE_SUPERSTRUCTURE, // Monumental pylons, chevron portal, or radial ring
    FORERUNNER_STAGE_TRANSIT,        // Hard-light bridges, bi-level cross spans, or perimeter catwalks
    FORERUNNER_STAGE_TERRACES,       // Subterranean vault galleries, observation wings, or side buttresses
    FORERUNNER_STAGE_APEX            // Levitating gravity core, telemetry conduits, or apex solar emitter
} ForerunnerStage;

// Backward-compatible aliases
#define FORERUNNER_STAGE_CHASM        FORERUNNER_STAGE_FOUNDATION
#define FORERUNNER_STAGE_GATEWAY      FORERUNNER_STAGE_SUPERSTRUCTURE
#define FORERUNNER_STAGE_SKYBRIDGE    FORERUNNER_STAGE_TRANSIT
#define FORERUNNER_STAGE_VAULT        FORERUNNER_STAGE_TERRACES
#define FORERUNNER_STAGE_CARTOGRAPHER FORERUNNER_STAGE_APEX

// Rules in the mass-void-path generative grammar
typedef enum {
    FORERUNNER_RULE_NONE = 0,
    // Axial Cartographer Rules
    FORERUNNER_RULE_EXCAVATE_CHASM,
    FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY,
    FORERUNNER_RULE_SPAN_SKYBRIDGE,
    FORERUNNER_RULE_TERRACE_VAULT,
    FORERUNNER_RULE_RAISE_GRAVITY_CORE,
    FORERUNNER_RULE_CARVE_LIGHT_CHANNELS,

    // Crossroads Rules
    FORERUNNER_RULE_EXCAVATE_CRUCIFORM_CHASM,
    FORERUNNER_RULE_ERECT_CORNER_PYLONS,
    FORERUNNER_RULE_SPAN_CROSS_BRIDGES,
    FORERUNNER_RULE_ADD_INTERSECTION_NEXUS,

    // Crucible Rules
    FORERUNNER_RULE_EXCAVATE_GRAVITY_PIT,
    FORERUNNER_RULE_ERECT_RADIAL_PYLON_RING,
    FORERUNNER_RULE_SPAN_RADIAL_CATWALK,
    FORERUNNER_RULE_SUSPEND_CENTRIFUGE_CORE,

    // Spire Citadel Rules
    FORERUNNER_RULE_EXCAVATE_ANNULAR_MOAT,
    FORERUNNER_RULE_ERECT_CENTRAL_PINNACLE,
    FORERUNNER_RULE_ADD_FLYING_BUTTRESSES,
    FORERUNNER_RULE_ATTACH_OBSERVATION_DECKS
} ForerunnerRule;

// Integer 3D bounding box (in modular half-meter coordinates)
typedef struct {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
} ForerunnerBox3D;

// Structural node in the Forerunner graph
typedef struct {
    ForerunnerPrimitiveType type;
    ForerunnerMaterial material;
    ForerunnerBox3D box;
    int parent_ids[MAX_FORERUNNER_PARENTS];
    int parent_count;
    float cant_angle_deg; // Inward rake angle X
    float cant_angle_z;   // Inward rake angle Z
    bool is_void;         // If true, represents an excavated void box
    float mass_weight;    // Mass rating for physics calculations
} ForerunnerNode;

// Complete generative plan
typedef struct {
    uint32_t seed;
    ForerunnerArchetype archetype;
    ForerunnerStage stage;
    ForerunnerNode nodes[MAX_FORERUNNER_NODES];
    int node_count;

    // Geometric Metrics
    int chasm_width;
    int chasm_depth;
    int chasm_length;

    // State Tracking
    bool has_chasm;
    bool has_gateway;
    bool has_bridge;
    bool has_vault;
    bool has_core;
    int pylon_count;
    int light_channel_count;
} ForerunnerPlan;

// Coalgebraic affordance opportunity
typedef struct {
    ForerunnerRule rule;
    ForerunnerPrimitiveType primitive_type;
    ForerunnerMaterial material;
    ForerunnerBox3D target_box;
    int parent_ids[MAX_FORERUNNER_PARENTS];
    int parent_count;
    float score;
} ForerunnerOpportunity;

// Core Bi-Algebra API
int forerunner_coalgebra_inspect(const ForerunnerPlan *plan, ForerunnerOpportunity *out_opps, int max_opps);
bool forerunner_algebra_apply(ForerunnerPlan *plan, const ForerunnerOpportunity *opp);
ForerunnerPlan generate_forerunner_structure(uint32_t seed, ForerunnerArchetype archetype, ForerunnerStage target_stage);
void print_forerunner_blueprint(const ForerunnerPlan *plan);

#ifdef __cplusplus
}
#endif

#endif // FORERUNNER_GRAMMAR_H
