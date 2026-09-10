#ifndef GREEK_GRAMMAR_H
#define GREEK_GRAMMAR_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Coarse architectural primitive types
typedef enum {
    ARCH_PRIMITIVE_NONE = 0,
    ARCH_PRIMITIVE_CELL,        // Enclosed chamber volume (Cella / Naos)
    ARCH_PRIMITIVE_WALL,        // Masonry boundary (with or without portal)
    ARCH_PRIMITIVE_COLUMN,      // Vertical support with base, shaft drums, capital
    ARCH_PRIMITIVE_STEP,        // Stepped stylobate / stereobate tier
    ARCH_PRIMITIVE_LINTEL,      // Architrave and entablature beam
    ARCH_PRIMITIVE_PEDIMENT,    // Triangular gable roof element
    ARCH_PRIMITIVE_COURTYARD,   // Open paved temenos enclosure
    ARCH_PRIMITIVE_POOL,        // Sunken basin containing liquid / PBF fluid
    ARCH_PRIMITIVE_TREE,        // Timber trunk with organic foliage canopy
    ARCH_PRIMITIVE_THOLOS,      // Circular monopteros rotunda / gazebo
    ARCH_PRIMITIVE_EXEDRA,      // Curved or rectangular marble philosopher bench
    ARCH_PRIMITIVE_NAISKOS,     // Miniature votive temple shrine with pediment
    ARCH_PRIMITIVE_FOUNTAIN,    // Garden fountain basin with water spout and PBF fluid
    ARCH_PRIMITIVE_BRAZIER,     // Stone pedestal fire altar / tripod with burning coals
    ARCH_PRIMITIVE_PERGOLA,     // Timber post-and-beam arbor trellis with vines
    ARCH_PRIMITIVE_PATH,        // Flagstone and gravel processional walkway
    ARCH_PRIMITIVE_SHRUB        // Low flowering shrubs, lavender, and planter beds
} ArchPrimitiveType;

// Cardinal faces of an architectural interface / growth site
typedef enum {
    SITE_FACE_FRONT = 0, // Facing sacred entrance (+Z)
    SITE_FACE_BACK  = 1, // Facing rear treasury (-Z)
    SITE_FACE_LEFT  = 2, // Flank (-X)
    SITE_FACE_RIGHT = 3, // Flank (+X)
    SITE_FACE_TOP   = 4  // Superstructure / Entablature (+Y)
} SiteFace;

// High-level historical maturation stage
typedef enum {
    TEMPLE_STAGE_CELLA = 0,        // Simple Naos cell
    TEMPLE_STAGE_PROSTYLE,         // Cella + front porch (Pronaos)
    TEMPLE_STAGE_AMPHIPROSTYLE,    // Cella + front & rear porches (Opisthodomos)
    TEMPLE_STAGE_PERIPTERAL,       // Surrounded by full colonnade (Peristyle)
    TEMPLE_STAGE_SANCTUARY         // Grand sacred precinct with Courtyard, Pool, and Grove
} TempleStage;

// Rules available in the grammar
typedef enum {
    RULE_NONE = 0,
    RULE_ATTACH_PRONAOS,          // Add front entrance porch with columns in-antis
    RULE_ATTACH_OPISTHODOMOS,     // Add mirrored rear porch across symmetry plane
    RULE_WRAP_PERISTYLE,          // Add full outer ambulatory colonnade
    RULE_EXPAND_STYLOBATE,        // Add tiered stepping podium below base
    RULE_RAISE_PEDIMENT,          // Add entablature and triangular pediment
    RULE_ATTACH_COURTYARD,        // Add peristyle forecourt with stoa wings
    RULE_EXCAVATE_POOL,           // Add sunken stone reflection basin with PBF fluid
    RULE_PLANT_GROVE,             // Add symmetrically paired sacred olive trees
    RULE_LAYOUT_GARDEN_PATHS,     // Add processional flagstone walkways encircling the precinct
    RULE_ERECT_THOLOS,            // Add sacred circular monopteros rotunda in rear glade
    RULE_INSTALL_EXEDRAE,         // Add paired curved marble philosopher benches
    RULE_BUILD_NAISKOI,           // Add paired miniature votive shrines
    RULE_CONSTRUCT_FOUNTAINS,     // Add paired garden fountains with PBF fluid
    RULE_BUILD_PERGOLAS,          // Add post-and-beam shaded vine trellises
    RULE_ERECT_BRAZIERS,          // Add ceremonial fire altars with glowing coals
    RULE_POPULATE_GARDEN_FLORA    // Add dense groves of cypresses, olives, and flowering shrubs
} GrammarRule;

// 3D Integer box in modular units
typedef struct {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
} ModBox3D;

// Semantic node placed in the building
typedef struct {
    ArchPrimitiveType type;
    ModBox3D box;             // Extents in modular units
    int column_count;         // If colonnade / row
    int column_spacing;       // Spacing between columns in modular units
    bool is_exterior;
    uint32_t flags;
} ArchNode;

// Boundary interface exposed by the current building state S
typedef struct {
    int site_id;
    SiteFace face;
    ModBox3D boundary;
    bool is_on_symmetry_axis;
    int symmetry_partner_id;  // Index of reflected site across X=0, or -1 if self
    bool active;
} GrowthSite;

// Candidate developmental opportunity emitted by coalgebra γ
typedef struct {
    int site_index;
    int partner_site_index;   // If paired by symmetry (-1 if none)
    GrammarRule rule;
    float score;
} GrowthOpportunity;

#define MAX_ARCH_NODES 1024
#define MAX_GROWTH_SITES 128
#define MAX_OPPORTUNITIES 64

// Full building state S
typedef struct {
    ArchNode nodes[MAX_ARCH_NODES];
    int node_count;

    GrowthSite sites[MAX_GROWTH_SITES];
    int site_count;

    // Proportional & modular geometry parameters
    int module_voxels;        // Number of voxels per module M (e.g. 2)
    int cella_width_m;        // Cella width in M (lateral, X)
    int cella_length_m;       // Cella length in M (longitudinal, Z)
    int cella_height_m;       // Cella wall height in M (vertical, Y)
    int column_spacing_m;     // Intercolumniation in M
    int column_height_m;      // Column height in M

    TempleStage stage;
    int complexity;
    uint32_t seed;
    bool has_pediment;
    bool has_stylobate;
    bool has_courtyard;
    bool has_pool;
    bool has_trees;
    bool has_garden_paths;
    bool has_tholos;
    bool has_exedrae;
    bool has_naiskoi;
    bool has_fountains;
    bool has_pergolas;
    bool has_braziers;
    bool has_garden_flora;
} TemplePlan;

// --- Coalgebra (γ): Observes S, finds exposed boundaries, pairs symmetries, scores options ---
int greek_coalgebra_inspect(const TemplePlan *plan, GrowthOpportunity *out_opps, int max_opps);

// --- Algebra (α): Applies winning rule r to site e, assembling parts into composite building ---
bool greek_algebra_apply(TemplePlan *plan, const GrowthOpportunity *opp);

// --- Bi-Algebra Stepper: Runs the developmental growth loop from seed to target complexity ---
TemplePlan generate_greek_temple(uint32_t seed, TempleStage target_stage, int module_voxels);

// --- Diagnostics & Blueprint Rendering ---
void print_temple_blueprint(const TemplePlan *plan);

#ifdef __cplusplus
}
#endif

#endif // GREEK_GRAMMAR_H
