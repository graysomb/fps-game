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
    FORERUNNER_PRIM_PYLON,           // 25-degree canted inward-leaning structural buttress
    FORERUNNER_PRIM_LINTEL,          // Heavy angled chevron lintel beam
    FORERUNNER_PRIM_BRIDGE_SPAN,     // Cantilevered solid metal deck
    FORERUNNER_PRIM_HARDLIGHT_BRIDGE,// Pure energy crossing over abyssal void
    FORERUNNER_PRIM_SPIRE,           // Towering vertical energy conduit spire
    FORERUNNER_PRIM_LIGHT_CHANNEL,   // Recessed glowing emissive conduit groove
    FORERUNNER_PRIM_GRAVITY_CORE,    // Hovering central installation matrix
    FORERUNNER_PRIM_TERRACE          // Stepped side gallery / wall tier
} ForerunnerPrimitiveType;

// Developmental maturation stages
typedef enum {
    FORERUNNER_STAGE_CHASM = 0,      // Abyssal trench & solitary sentinel beacon
    FORERUNNER_STAGE_GATEWAY,        // Monumental chevron entrance portal
    FORERUNNER_STAGE_SKYBRIDGE,      // High-altitude cantilever & hard-light crossing
    FORERUNNER_STAGE_VAULT,          // Terraced chamber galleries & flanking buttresses
    FORERUNNER_STAGE_CARTOGRAPHER    // Apex monument: hovering core, telemetry spires & full chasm
} ForerunnerStage;

// Rules in the mass-void-path generative grammar
typedef enum {
    FORERUNNER_RULE_NONE = 0,
    FORERUNNER_RULE_EXCAVATE_CHASM,
    FORERUNNER_RULE_ERECT_CHEVRON_GATEWAY,
    FORERUNNER_RULE_SPAN_SKYBRIDGE,
    FORERUNNER_RULE_TERRACE_VAULT,
    FORERUNNER_RULE_RAISE_GRAVITY_CORE,
    FORERUNNER_RULE_CARVE_LIGHT_CHANNELS
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
    float cant_angle_deg; // Inward rake angle (typically 20-30 deg)
    bool is_void;         // If true, represents an excavated void box
    float mass_weight;    // Mass rating for physics calculations
} ForerunnerNode;

// Complete generative plan
typedef struct {
    uint32_t seed;
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
ForerunnerPlan generate_forerunner_structure(uint32_t seed, ForerunnerStage target_stage);
void print_forerunner_blueprint(const ForerunnerPlan *plan);

#ifdef __cplusplus
}
#endif

#endif // FORERUNNER_GRAMMAR_H
