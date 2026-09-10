#ifndef HYPERBOREAN_GRAMMAR_H
#define HYPERBOREAN_GRAMMAR_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Maximum capacity constants
#define MAX_HYPER_NODES 256
#define MAX_HYPER_PARENTS 8
#define MAX_HYPER_OPPORTUNITIES 32

// Material categories contrasting raw megalithic mass with refined Hellenic marble
typedef enum {
    HYPER_MAT_WEATHERED_SARSEN = 0, // Ancient dark/grey megalithic stone with moss/lichen
    HYPER_MAT_PENTELIC_MARBLE,      // Gleaming white fluted marble
    HYPER_MAT_GOLD_ACCENT,          // Chased gold leaf on friezes, tripods & finials
    HYPER_MAT_TERRACOTTA,           // Classical tiled roof elements
    HYPER_MAT_PBF_WATER,            // Sacred spring / reflection moat fluid
    HYPER_MAT_EMBER,                // Eternal ceremonial flame / fire brazier
    HYPER_MAT_CYPRESS,              // Sacred grove trees
    HYPER_MAT_BRONZE                // Sacrificial tripod & statues
} HyperMaterial;

// Primitive types in the Syncretic Hyperborean vocabulary
typedef enum {
    HYPER_PRIM_NONE = 0,
    // Megalithic Chthonic
    HYPER_PRIM_MENHIR,              // Rough standing stone / orthostat
    HYPER_PRIM_ROUGH_LINTEL,        // Cyclopean bridging beam
    HYPER_PRIM_HEEL_STONE,          // Colossal solar marker monolith
    // Classical Apollonian
    HYPER_PRIM_FLUTED_COLUMN,       // Carved Doric column with flutes & echinus capital
    HYPER_PRIM_ARCHITRAVE,          // Polished epistyle beam & triglyph-metope frieze
    HYPER_PRIM_PEDIMENT,            // Triangular pediment with carved relief tympanum
    HYPER_PRIM_THOLOS_PODIUM,       // Concentric tiered circular marble steps
    HYPER_PRIM_FLUID_BASIN,         // Sunken circular cistern for sacred PBF spring
    HYPER_PRIM_BRAZIER,             // Bronze ritual tripod with perpetual flame
    HYPER_PRIM_STATUE,              // Votive marble herm / Apollo statue
    HYPER_PRIM_CYPRESS,             // Slender Mediterranean cypress tree
    // Syncretic Hybrid
    HYPER_PRIM_HYBRID_TRILITHON,    // Colossal sarsen orthostats carrying Classical Doric entablature
    HYPER_PRIM_CYCLOPEAN_ARCHITRAVE // Massive stone lintel carved with classical guttae & triglyphs
} HyperPrimitiveType;

// Developmental maturation stages
typedef enum {
    HYPER_STAGE_AVENUE = 0,         // Solstice dromos corridor & Heel Stone marker
    HYPER_STAGE_OUTER_HENGE,        // Perimeter cromlech ring of 16 sarsen menhirs + lintels
    HYPER_STAGE_CLASSICAL_PERISTYLE,// Middle concentric ring of 12 fluted marble Doric columns
    HYPER_STAGE_GREAT_TRILITHONS,   // Inner horseshoe of colossal hybrid trilithons with Doric pediments
    HYPER_STAGE_FULL_SANCTUM        // Complete complex: Sunken Tholos, PBF spring moat, flame, and cypress alley
} HyperStage;

// Rules in the syncretic generative grammar
typedef enum {
    HYPER_RULE_NONE = 0,
    HYPER_RULE_ALIGN_SOLSTICE_AVENUE,
    HYPER_RULE_ERECT_OUTER_HENGE,
    HYPER_RULE_ERECT_MARBLE_PERISTYLE,
    HYPER_RULE_RAISE_HYBRID_TRILITHONS,
    HYPER_RULE_EXCAVATE_THOLOS_BASIN,
    HYPER_RULE_POPULATE_SACRED_GROVE
} HyperRule;

// Integer 3D bounding box (in modular half-meter coordinate space)
typedef struct {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
} HyperBox3D;

// Structural node in the synthesized complex
typedef struct {
    HyperPrimitiveType type;
    HyperMaterial material;
    HyperBox3D box;
    int parent_ids[MAX_HYPER_PARENTS];
    int parent_count;
    float stability_weight;
} HyperNode;

// Complete blueprint plan for the Hyperborean Sun-Henge
typedef struct {
    uint32_t seed;
    HyperStage stage;
    HyperNode nodes[MAX_HYPER_NODES];
    int node_count;

    // Astronomical Alignment (Midsummer Sunrise axis)
    float axis_dx, axis_dz;

    // Concentric Ring Radii
    int r_outer_henge;      // ~17 modules
    int r_mid_peristyle;    // ~11 modules
    int r_inner_trilithon;  // ~6 modules
    int r_tholos_core;      // ~3 modules

    // State Tracking
    bool has_avenue;
    int outer_stones_placed;
    int peristyle_cols_placed;
    int trilithons_placed;
    bool has_tholos;
    bool has_moat;
    bool has_brazier;
} HyperPlan;

// Coalgebraic affordance opportunity
typedef struct {
    HyperRule rule;
    HyperPrimitiveType primitive_type;
    HyperMaterial material;
    HyperBox3D target_box;
    int parent_ids[MAX_HYPER_PARENTS];
    int parent_count;
    float score;
} HyperOpportunity;

// Core Bi-Algebra API
int hyperborean_coalgebra_inspect(const HyperPlan *plan, HyperOpportunity *out_opps, int max_opps);
bool hyperborean_algebra_apply(HyperPlan *plan, const HyperOpportunity *opp);
HyperPlan generate_hyperborean_structure(uint32_t seed, HyperStage target_stage);
void print_hyperborean_blueprint(const HyperPlan *plan);

#ifdef __cplusplus
}
#endif

#endif // HYPERBOREAN_GRAMMAR_H
