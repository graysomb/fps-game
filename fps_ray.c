// Raylib split-screen FPS prototype (port of fps_game.c)
#include "raylib.h"
#include "rlgl.h" // for rlBegin/rlEnd
#include "raymath.h" // for MatrixIdentity()
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <stdarg.h>

// Game state enum
typedef enum {
    GAME_STATE_MENU,
    GAME_STATE_PLAYING,
    GAME_STATE_PAUSED,
    GAME_STATE_SETTINGS
} GameState;
static GameState gameState = GAME_STATE_MENU;

// Input type enum
#define MAX_PLAYERS 4
typedef enum {
    INPUT_TYPE_KEYBOARD,
    INPUT_TYPE_GAMEPAD
} InputType;
static InputType playerInput[MAX_PLAYERS] = {
    INPUT_TYPE_KEYBOARD,
    INPUT_TYPE_KEYBOARD,
    INPUT_TYPE_GAMEPAD,
    INPUT_TYPE_GAMEPAD
};

//gcc fps_ray.c -o fps_ray  $(pkg-config --cflags --libs raylib) -lm -Wl,-rpath,/usr/local/lib
// Screen and game constants
#define SCREEN_WIDTH    1600
#define SCREEN_HEIGHT    900
#define MOVE_SPEED      5.0f    // units per second (unused: using acceleration)
#define TURN_SPEED     90.0f    // degrees per second
#define JUMP_SPEED     7.5f    // initial jump velocity
#define GRAVITY         9.8f*1.0f    // gravity acceleration
#define BASE_EYE_HEIGHT 1.0f    // player eye height above floor
#define HUD_FONT_SIZE      28
#define HUD_BAR_HEIGHT      6
#define HUD_PADDING_X      10
#define HUD_PADDING_Y       8
#define HUD_BAR_WIDTH     320
#define HUD_BAR_THICKNESS  30
#define HUD_BAR_SPACING     8
#define HUD_BAR_TEXT_SIZE  22
#define HUD_BULLET_TEXT_SIZE 78
#define AMMO_MAX 6
#define AMMO_RECHARGE_SECONDS 2.0f
#define HUD_AMMO_SIZE 12
#define HUD_AMMO_GAP 4
#define BULLET_COOLDOWN_SECONDS 0.25f
#define ACCELERATION   400.0f    // horizontal acceleration
#define FREEZE_GROUND_WEIGHT       0.9f
#define FREEZE_NEIGHBOR_WEIGHT     0.4f*0.0f
#define FREEZE_GROUND_REF_HEIGHT   6.0f
#define FREEZE_PROPAGATION_ITERATIONS 100
#define FREEZE_PROPAGATION_ATTENUATION 1.0f
#define FREEZE_PROPAGATION_EPSILON 1e-6f
#define FREEZE_PATH_DECAY 0.85f
#define FREEZE_OVERHANG_DECAY 0.7f
#define ACTIVATION_VELOCITY_WEIGHT 0.6f
#define ACTIVATION_VELOCITY_REF_SPEED 6.0f
#define ACTIVATION_STRAIN_WEIGHT   0.1f
#define ACTIVATION_STRAIN_REF      0.1f
#define ACTIVATION_GLUE_WEIGHT     0.1f
#define ACTIVATION_GLUE_REF_SPEED  3.0f
#define ACTIVATION_DYNAMIC_WEIGHT  0.5f
#define FREEZE_BELIEF_IMPORTANCE   0.9f
#define ACTIVATION_HYSTERESIS      0.1f
#define FRICTION       400.0f    // ground friction deceleration
#define TURN_ACCELERATION 400.0f
#define TURN_FRICTION 400.0f
#define PLAYER_RADIUS   0.5f
#define FLOOR_SIZE     20.0f*1.0f    // half-size of floor in world units
#define MAX_STATIC_COLLISION_NEIGHBORS 64
#define PLAYER_SIZE 0.5f
#define PARTICLE_RADIUS (VOXEL_SIZE * 0.5f)
#define VGS_ALPHA 0.75f
#define VGS_BETA 0.35f
#define VGS_ITERS 3
#define VGS_EPS 1e-6f
#define VGS_EARLY_OUT_EPS 0.0002f
#define VOXEL_CORNER_COUNT 8
#define VOXEL_CENTER_INDEX 8
#define VOXEL_PARTICLE_COUNT 9
#define PBD_MAX_STEP_DT 0.005f
#define PBD_SUBSTEPS 2
#define PBD_CONSTRAINT_ITERS 6
#define COLLISION_RELAXATION 0.99f
#define COLLISION_CENTROID_ONLY_DT 0.02f
#define SPLIT_VELOCITY_DAMP 0.1f
#define CENTER_RELAXATION 0.99f
#define VELOCITY_DAMPING .999f
#define GLUE_RELAXATION 0.9f*0.4f
//#define GLUE_EPS 1e-6f
#define GLUE_EPS 0.0002f
#define GLUE_BREAK_STRAIN 0.4f
#define GLUE_BREAK_HINGE_ANGLE_DEG 20.0f
#define GLUE_BREAK_VELOCITY_SKIP_FRAMES 30
#define GLUE_VIRTUAL_EDGE_STRENGTH 0.4f*0.4f
#define GLUE_VIRTUAL_CENTER_STRENGTH 0.2f*0.4f
#define RECYCLE_DYNAMIC_MAX_FRAMES (60 * 10)
#define RECYCLE_STATIC_RESTORE_INTERVAL 1
#define RECYCLE_STATIC_RESTORE_DELAY (60 * 10)
#define RECYCLE_OWNED_STATIC_MAX_FRAMES (60 * 10)
#define STATIC_DEBRIS_OWNER (-2)
#define VOXEL_SPLIT_STRAIN_THRESHOLD 0.9f
#define VOXEL_SPLIT_SHEAR_THRESHOLD 0.9f
#define VOXEL_DUST_STRAIN_THRESHOLD 2.0f
#define VOXEL_DUST_SHEAR_THRESHOLD 0.98f
#define VOXEL_HASH_REBUILD_INTERVAL 2
#define TABLE_CACHE_SIZE 4
#define FACE_BLOCK_MIN_OVERLAP (VOXEL_SIZE * 0.25f)

#define MAX_NEIGHBOR_VOXELS 128
#define MAX_FACE_NEIGHBORS   64
#define GLUE_NEIGHBOR_HASH_SIZE 128
#define MAX_SPLIT_CHILDREN    8
static const float GRID_EPSILON = 1e-4f;
#define VOXEL_ACTIVATION_RADIUS 2*1
//#define VOXEL_ACTIVATION_UNIT_BUDGET 128
#define VOXEL_ACTIVATION_UNIT_BUDGET 128*5
#define VOXEL_DEACTIVATION_VELOCITY_THRESHOLD 1.0f
#define VOXEL_DEACTIVATION_STRAIN_THRESHOLD 0.2f
#define VOXEL_DEACTIVATION_SHEAR_THRESHOLD 0.2f
#define VOXEL_DEACTIVATION_FRAMES 10
#define VOXEL_MAX_DEACTIVATIONS_PER_FRAME 128*5
#define STATIC_RESTORE_SEARCH_RADIUS 2*1
#define DEBRIS_ACTIVATION_COOLDOWN_FRAMES (60 * 10)
#define STATIC_REBUILD_ACTIVATION_COOLDOWN_FRAMES (60 * 10)

// KD-stats constants
#define BASE_HEALTH 100
#define BASE_SHIELD 150
#define SHIELD_REGEN_DELAY 5.0f
#define VOXEL_DAMAGE 50

// Voxel physics constants
#define MAX_VOXELS    131072
#define HASH_SIZE     131072    // must be power of two
#define VOXEL_SIZE     0.5f    // size of each voxel cube

// Tunable voxel edit brush (per-axis span of the add/remove operation)
static int voxelBrushSpan = 3;

// Player structure
typedef struct {
    Vector3 pos;
    float yaw;
    float pitch;
    float yaw_vel;
    float pitch_vel;
    Vector3 vel;
    bool onGround;
    int vType;
    int kills;
    int debrisKills;
    int deaths;
    float kd_ratio;
    int health;
    int shield;
    float last_damage_time;
    int ammo;
    float ammo_recharge_timer;
    float last_shot_time;
} Player;
static Player players[MAX_PLAYERS];
static int activePlayers = 2;
static bool randomSpawnEnabled = true;
static const Vector3 playerSpawnPositions[MAX_PLAYERS] = {
    { 0.0f,  BASE_EYE_HEIGHT, -9.0f },
    { 0.0f,  BASE_EYE_HEIGHT,  9.0f },
    { -6.0f, BASE_EYE_HEIGHT, -6.0f },
    {  6.0f, BASE_EYE_HEIGHT,  6.0f }
};
static const float playerSpawnYaw[MAX_PLAYERS] = { 0.0f, 180.0f, 45.0f, -135.0f };

typedef struct {
    Vector3 pos;
    Vector3 prev_pos;
    Vector3 predicted_pos;
    Vector3 vel;
    float inv_mass;
} Particle;

// Voxel structure
typedef struct {
    Vector3 pos;
    Vector3 vel;
    bool simulate;
    bool fixed;
    bool glueEligible;
    bool pendingActivation;
    bool isBullet;
    int activationCooldownFrames;
    Color color;
    int type;
    int owner;
    int activator;
    /* Visible faces mask: surface[i]=true if face i is visible:
       0=+X,1=-X,2=+Y,3=-Y,4=+Z,5=-Z */
    bool surface[6];
    int  gx, gy, gz;
    int  span;
    int  min_gx, max_gx;
    int  min_gy, max_gy;
    int  min_gz, max_gz;
    int  rest_min_gx, rest_max_gx;
    int  rest_min_gy, rest_max_gy;
    int  rest_min_gz, rest_max_gz;
    int  orig_min_gx, orig_max_gx;
    int  orig_min_gy, orig_max_gy;
    int  orig_min_gz, orig_max_gz;
    Particle particles[8];
    Particle center_particle;
    float rest_volume;
    float rest_edge;
    float particle_radius;
    int sleepFrames;
    float freezeBelief;
    float activationBelief;
    float groundSupport;
    uint8_t neighborSupport;
    uint8_t supportMask;
    uint8_t skipCollisionVelocityFrames;
    int lifeFrames;
    int debugClusterTag;
    int prevGlueClusterId;
    int prevGlueClusterSize;
    int prevGlueClusterTag;
    uint8_t prevGlueClusterValid;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;
static int glueClusterIndices[MAX_VOXELS];
static unsigned char glueClusterVisited[MAX_VOXELS];
static bool dynamicGlueClustersInitialized = false;
#define DEBUG_CLUSTER_TAG_MAX 64
static int debugTagOffset[DEBUG_CLUSTER_TAG_MAX][3];
static float freezeBeliefScratch[MAX_VOXELS];
static int freezeDistance[MAX_VOXELS];
static int overhangDistance[MAX_VOXELS];
static uint8_t freezeBoundaryFlags[MAX_VOXELS];
static uint8_t staticBeliefDirty[MAX_VOXELS];
static int staticBeliefDirtyList[MAX_VOXELS];
static int staticBeliefDirtyCount = 0;
static uint8_t staticBeliefQueued[MAX_VOXELS];
static int staticBeliefQueue[MAX_VOXELS];
static bool staticBeliefsInitialized = false;
static bool staticBeliefsForceFullRefresh = false;
static const int DEBUG_SPAN2_FACE_LOG_INIT = 32;
static int debugSpan2FaceLogBudget = 0;
static int debugSmushLogBudget = 0;
static bool debugLogSpan2Faces = false;
static Voxel recycleQueue[MAX_VOXELS];
static int recycleQueueHead = 0;
static int recycleQueueTail = 0;
static int recycleQueueCount = 0;
static int recycleFrameCounter = 0;
static int addVoxel(float px, float py, float pz, bool fixed, bool simulate, Color color, int type);
static void rebuild_voxel_hash(void);
static void voxel_measure_strain(const Voxel *voxel,
                                 float *out_max_strain,
                                 float *out_max_shear);
static int build_glue_cluster_indices(int start_idx, int *out_indices);
static int gather_glued_neighbors(int voxel_idx, int *out, int max_out);
static bool restore_glue_cluster_to_static(const int *cluster, int cluster_count);
static void solve_span_voxel_collisions(void);
static void log_dynamic_glue_cluster_breaks(void);
static void compute_voxel_center_and_mass(const Voxel *voxel, Vector3 *center, float *inv_mass_sum);
typedef struct {
    int gx, gy, gz;
    Color color;
    int type;
    bool fixed;
    int voxelIndex;
    int debugTag;
    int activator;
} UnitVoxelSeed;

typedef struct {
    UnitVoxelSeed voxels[MAX_VOXELS];
    int count;
} UnitVoxelBuffer;

static void refresh_static_voxel_beliefs(void);
static void update_static_voxel_belief(int idx);
static void mark_static_beliefs_dirty_for_voxel(const Voxel *voxel);
static void mark_static_beliefs_dirty_column_above(int gx, int gz, int gy);
static void update_dynamic_activation_beliefs(void);
static float compute_cluster_freeze_belief(const UnitVoxelBuffer *buffer, int startIndex);
static void rollback_activation_buffer(UnitVoxelBuffer *buffer, int startIndex);
static bool dynamic_belief_overcomes_static(float dynamicBelief, float frozenBelief);
static uint8_t compute_static_support_mask(const Voxel *voxel);
static bool list_contains_index(const int *list, int count, int value);
static int gather_static_face_neighbors(const Voxel *voxel, int *out, int max_out);
static int gather_static_voxels_near_point(Vector3 point, float radius, int *out, int max_out);
static bool push_particle_out_of_static(const Voxel *static_voxel, Particle *particle, float radius);
static int gather_static_voxels_near_voxel(const Voxel *voxel, int *out, int max_out);
static bool resolve_span_static_overlap(int dynamic_idx, Voxel *dynamic,
                                        int static_idx, const Voxel *static_voxel);
static bool nudge_voxel_bottom_above_static(int voxel_idx, Voxel *voxel);
static void glue_dynamic_voxel_to_static_neighbors(void);
static void glue_dynamic_voxel_to_static_neighbors_for_voxel(int voxel_idx);
static bool recycle_queue_push(const Voxel *voxel);
static bool recycle_queue_pop(Voxel *out);
static bool voxel_is_at_rest_location(const Voxel *voxel);
static bool voxel_outside_world_bounds(const Voxel *voxel);
static bool spawn_static_at_rest(const Voxel *snapshot);
static void recycle_dead_voxels(void);
static void update_projectiles(float dt);
static void handle_pbd_projectile_hits(void);
static Vector3 v_add(Vector3 a, Vector3 b);
static Vector3 v_sub(Vector3 a, Vector3 b);

static void unit_voxel_buffer_clear(UnitVoxelBuffer *buffer) {
    if (!buffer) {
        return;
    }
    int existing = buffer->count;
    if (existing < 0 || existing > MAX_VOXELS) {
        existing = 0;
    }
    for (int i = 0; i < existing; ++i) {
        int idx = buffer->voxels[i].voxelIndex;
        if (idx >= 0 && idx < voxel_count) {
            voxels[idx].pendingActivation = false;
        }
    }
    buffer->count = 0;
}

static bool unit_voxel_buffer_push(UnitVoxelBuffer *buffer,
                                   int gx, int gy, int gz,
                                   Color color, int type,
                                   bool fixed, int voxelIndex,
                                   int debugTag, int activator)
{
    if (!buffer || buffer->count >= MAX_VOXELS) {
        return false;
    }
    buffer->voxels[buffer->count++] = (UnitVoxelSeed){
        .gx = gx,
        .gy = gy,
        .gz = gz,
        .color = color,
        .type = type,
        .fixed = fixed,
        .voxelIndex = voxelIndex,
        .debugTag = debugTag,
        .activator = activator
    };
    return true;
}

static void rollback_activation_buffer(UnitVoxelBuffer *buffer, int startIndex)
{
    if (!buffer) {
        return;
    }
    if (startIndex < 0) {
        startIndex = 0;
    }
    if (startIndex > buffer->count) {
        return;
    }
    for (int i = startIndex; i < buffer->count; ++i) {
        int idx = buffer->voxels[i].voxelIndex;
        if (idx >= 0 && idx < voxel_count) {
            voxels[idx].pendingActivation = false;
        }
    }
    buffer->count = startIndex;
}

static inline float voxel_particle_radius(const Voxel *v) {
    float r = v->particle_radius;
    if (r <= 0.0f) {
        r = 0.5f * v->rest_edge;
    }
    float clamp = PARTICLE_RADIUS * (float)(v->span > 0 ? v->span : 1);
    if (clamp > 0.0f) {
        r = fminf(r, clamp);
    }
    return r;
}

static bool debugDrawParticles = false;
static bool debugColorParticlesByVelocity = false;
static const float PARTICLE_DEBUG_MARKER_RADIUS = 0.6f;
static const float PARTICLE_DEBUG_MAX_SPEED = 20.0f;
static bool debugLogSpanCollisions = false;
static int debugSpanEdgeLogBudget = 0;
static int debugSpanCollisionLogBudget = 0;
static bool debugLogGlue = false;
static bool debugLogGlueClusters = false;
static bool skipGlueClusterCollisions = true;
static bool renderAllDynamicFaces = false;
static bool debugLogDynamicVoxels = false;
static bool debugLogVoxelRecycle = false;
static bool debugLogActivationFailures = true;
static bool debugLogRestoreFailures = true;
static bool debugLogVoxelBlowup = true;
static bool debugLogSmush = true;
static bool debugLogSmushSpawns = false;
static bool debugLogSmushHits = false;
static bool debugLogSmushDeaths = false;
static bool debugLogMultiscale = false;
static bool debugLogRestoreClusters = false;
static int debugBlowupLogBudget = 0;
static bool debugLogClusterBreaksOnly = true;
static int debugGlueBuildLogBudget = 0;
static int debugGlueSolveLogBudget = 0;
static int debugGlueBreakLogBudget = 0;
static const int DEBUG_GLUE_BUILD_LOG_INIT = 64;
static const int DEBUG_GLUE_SOLVE_LOG_INIT = 64;
static const int DEBUG_GLUE_BREAK_LOG_INIT = 16;
static bool debugLogActivation = false;
static const int DEBUG_ACTIVATION_LOG_INIT = 64;
static int activationGlueLogBudget = 0;
static int activationLogStart = -1;
static int activationLogEnd = -1;
static bool debugLogVoxelDeactivation = false;
static FILE *debugLogFile = NULL;
static int debugSupportLogBudget = 512;
static const float DEBUG_FALL_LOG_THRESHOLD = -5.0f;
static const int DEBUG_FALL_LOG_BUDGET = 32;
static bool debugLogFall = false;
static bool debugShowBeliefColors = false;
static unsigned char debugTagBreakLogged[DEBUG_CLUSTER_TAG_MAX];

static const char *trace_level_label(int level) {
    switch (level) {
        case LOG_TRACE: return "TRACE";
        case LOG_DEBUG: return "DEBUG";
        case LOG_INFO: return "INFO";
        case LOG_WARNING: return "WARN";
        case LOG_ERROR: return "ERROR";
        case LOG_FATAL: return "FATAL";
        default: return "LOG";
    }
}

static const char *debug_cluster_tag_label(int tag) {
    switch (tag) {
        case 1: return "corner-chunk";
        case 2: return "stacked-pillar";
        case 3: return "three-span-bar";
        case 4: return "plate-cap";
        case 5: return "step-chain";
        case 6: return "solid-cube";
        case 7: return "long-beam";
        case 8: return "brace-frame";
        case 9: return "staggered-stack";
        case 10: return "unit-cube-control";
        case 11: return "single-span-control";
        case 12: return "thin-slab";
        case 13: return "flat-cross";
        case 14: return "vertical-column";
        default: return "unlabeled";
    }
}

static bool debug_should_log_tag_break(int tag)
{
    if (tag <= 0 || tag >= DEBUG_CLUSTER_TAG_MAX) {
        return false;
    }
    if (debugTagBreakLogged[tag]) {
        return false;
    }
    return true;
}

static bool debug_should_log_message(const char *text) {
    if (!debugLogClusterBreaksOnly) {
        return true;
    }
    if (!text) {
        return false;
    }
    if (strstr(text, "[GlueClusterBreak]")) {
        return true;
    }
    if (strstr(text, "[GlueTagBreak]")) {
        return true;
    }
    if (strstr(text, "[Recycle]")) {
        return true;
    }
    if (strstr(text, "[ActivationFail]")) {
        return true;
    }
    if (strstr(text, "[RestoreFail]")) {
        return true;
    }
    if (strstr(text, "[VoxelBlowup]")) {
        return true;
    }
    if (strstr(text, "[DynamicVoxel]")) {
        return true;
    }
    if (strstr(text, "[Smush]")) {
        return true;
    }
    if (text[0] == '[') {
        return true;
    }
    if (strstr(text, "Logging TraceLog output to")) {
        return true;
    }
    if (strstr(text, "Failed to open log file")) {
        return true;
    }
    return false;
}

static void DebugTraceLogCallback(int logLevel, const char *text, va_list args) {
    if (!debug_should_log_message(text)) {
        return;
    }
    FILE *console = (logLevel >= LOG_WARNING) ? stderr : stdout;
    //FILE *console = (logLevel >= LOG_FATAL) ? stderr : stdout;
    va_list console_args;
    va_copy(console_args, args);
    vfprintf(console, text, console_args);
    fprintf(console, "\n");
    va_end(console_args);

    if (!debugLogFile) {
        return;
    }

    va_list file_args;
    va_copy(file_args, args);
    fprintf(debugLogFile, "[%s] ", trace_level_label(logLevel));
    vfprintf(debugLogFile, text, file_args);
    fprintf(debugLogFile, "\n");
    fflush(debugLogFile);
    va_end(file_args);
}

static void InitDebugLogging(void) {
    const char *log_path = getenv("FPS_RAY_LOG");
    if (!log_path || log_path[0] == '\0') {
        log_path = "fps_ray.log";
    }
    debugLogFile = fopen(log_path, "w");
    SetTraceLogCallback(DebugTraceLogCallback);
    if (debugLogFile) {
        TraceLog(LOG_INFO, "Logging TraceLog output to %s", log_path);
    } else {
        TraceLog(LOG_WARNING, "Failed to open log file at %s", log_path);
    }
}

static bool debug_should_log_span_pair(const Voxel *a, const Voxel *b, int *budget) {
    if (!debugLogSpanCollisions || budget == NULL || *budget <= 0) {
        return false;
    }
    if (a->span <= 1 && b->span <= 1) {
        return false;
    }
    --(*budget);
    return true;
}

static const int corner_signs[8][3] = {
    { -1, -1, -1 }, {  1, -1, -1 },
    { -1,  1, -1 }, {  1,  1, -1 },
    { -1, -1,  1 }, {  1, -1,  1 },
    { -1,  1,  1 }, {  1,  1,  1 }
};

static const int voxel_edge_pairs[12][2] = {
    {0,1},{1,3},{3,2},{2,0},
    {4,5},{5,7},{7,6},{6,4},
    {0,4},{1,5},{3,7},{2,6}
};

typedef struct {
    float minx, maxx;
    float miny, maxy;
    float minz, maxz;
} VoxelWorldBounds;

static void voxel_world_bounds(const Voxel *v, VoxelWorldBounds *out) {
    if (!out || !v) {
        return;
    }
    float half = 0.5f * v->rest_edge;
    out->minx = v->pos.x - half;
    out->maxx = v->pos.x + half;
    out->miny = v->pos.y - half;
    out->maxy = v->pos.y + half;
    out->minz = v->pos.z - half;
    out->maxz = v->pos.z + half;
}

static void voxel_particle_world_bounds(const Voxel *v, VoxelWorldBounds *out)
{
    if (!out || !v) {
        return;
    }
    Vector3 p = v->particles[0].pos;
    float minx = p.x, maxx = p.x;
    float miny = p.y, maxy = p.y;
    float minz = p.z, maxz = p.z;
    for (int i = 1; i < 8; ++i) {
        p = v->particles[i].pos;
        if (p.x < minx) minx = p.x;
        if (p.x > maxx) maxx = p.x;
        if (p.y < miny) miny = p.y;
        if (p.y > maxy) maxy = p.y;
        if (p.z < minz) minz = p.z;
        if (p.z > maxz) maxz = p.z;
    }
    out->minx = minx;
    out->maxx = maxx;
    out->miny = miny;
    out->maxy = maxy;
    out->minz = minz;
    out->maxz = maxz;
}

static void voxel_predicted_bounds(const Voxel *v, VoxelWorldBounds *out)
{
    if (!out || !v) {
        return;
    }
    Vector3 p = v->particles[0].predicted_pos;
    float minx = p.x, maxx = p.x;
    float miny = p.y, maxy = p.y;
    float minz = p.z, maxz = p.z;
    for (int i = 1; i < 8; ++i) {
        p = v->particles[i].predicted_pos;
        if (p.x < minx) minx = p.x;
        if (p.x > maxx) maxx = p.x;
        if (p.y < miny) miny = p.y;
        if (p.y > maxy) maxy = p.y;
        if (p.z < minz) minz = p.z;
        if (p.z > maxz) maxz = p.z;
    }
    out->minx = minx;
    out->maxx = maxx;
    out->miny = miny;
    out->maxy = maxy;
    out->minz = minz;
    out->maxz = maxz;
}

static void voxel_visibility_bounds(const Voxel *v, VoxelWorldBounds *out)
{
    if (!v || !out) {
        return;
    }
    if (v->simulate) {
        voxel_particle_world_bounds(v, out);
    } else {
        voxel_world_bounds(v, out);
    }
}

static bool axis_contact_state(float minA, float maxA,
                               float minB, float maxB,
                               float eps, bool *touch, bool *overlap)
{
    if (maxA < minB - eps || maxB < minA - eps) {
        if (touch) *touch = false;
        if (overlap) *overlap = false;
        return false;
    }

    float overlapMin = fmaxf(minA, minB);
    float overlapMax = fminf(maxA, maxB);
    bool touching = (overlapMax - overlapMin) <= eps;

    if (touch) *touch = touching;
    if (overlap) *overlap = !touching;
    return true;
}

static bool bounds_overlap(float minA, float maxA,
                           float minB, float maxB, float eps)
{
    return !(maxA < minB - eps || maxB < minA - eps);
}

static float bounds_overlap_length(float minA, float maxA,
                                   float minB, float maxB)
{
    float lo = fmaxf(minA, minB);
    float hi = fminf(maxA, maxB);
    float len = hi - lo;
    return (len > 0.0f) ? len : 0.0f;
}

static void translate_voxel_particles(Voxel *voxel, Vector3 delta)
{
    if (!voxel) {
        return;
    }
    for (int i = 0; i < 8; ++i) {
        voxel->particles[i].pos = v_add(voxel->particles[i].pos, delta);
        voxel->particles[i].prev_pos = v_add(voxel->particles[i].prev_pos, delta);
        voxel->particles[i].predicted_pos = v_add(voxel->particles[i].predicted_pos, delta);
    }
}

static bool segment_intersects_aabb(Vector3 start, Vector3 end,
                                    Vector3 box_min, Vector3 box_max)
{
    Vector3 dir = v_sub(end, start);
    float tmin = 0.0f;
    float tmax = 1.0f;
    const float eps = 1e-6f;

    float min_vals[3] = { box_min.x, box_min.y, box_min.z };
    float max_vals[3] = { box_max.x, box_max.y, box_max.z };
    float start_vals[3] = { start.x, start.y, start.z };
    float dir_vals[3] = { dir.x, dir.y, dir.z };

    for (int axis = 0; axis < 3; ++axis) {
        float d = dir_vals[axis];
        float s = start_vals[axis];
        if (fabsf(d) < eps) {
            if (s < min_vals[axis] || s > max_vals[axis]) {
                return false;
            }
            continue;
        }
        float inv_d = 1.0f / d;
        float t1 = (min_vals[axis] - s) * inv_d;
        float t2 = (max_vals[axis] - s) * inv_d;
        if (t1 > t2) {
            float tmp = t1;
            t1 = t2;
            t2 = tmp;
        }
        if (t1 > tmin) tmin = t1;
        if (t2 < tmax) tmax = t2;
        if (tmin > tmax) {
            return false;
        }
    }
    return true;
}

static bool face_blocked_by_voxel(const Voxel *self, const Voxel *neighbor, int face)
{
    if (!self || !neighbor) {
        return false;
    }
    VoxelWorldBounds a, b;
    voxel_visibility_bounds(self, &a);
    voxel_visibility_bounds(neighbor, &b);
    const float eps = VOXEL_SIZE * 0.05f;

    switch (face) {
        case 0: // +X
            if (b.minx > a.maxx + eps || b.maxx < a.maxx - eps) return false;
            return bounds_overlap(a.miny, a.maxy, b.miny, b.maxy, eps) &&
                   bounds_overlap(a.minz, a.maxz, b.minz, b.maxz, eps) &&
                   bounds_overlap_length(a.miny, a.maxy, b.miny, b.maxy) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.minz, a.maxz, b.minz, b.maxz) >= FACE_BLOCK_MIN_OVERLAP;
        case 1: // -X
            if (b.maxx < a.minx - eps || b.minx > a.minx + eps) return false;
            return bounds_overlap(a.miny, a.maxy, b.miny, b.maxy, eps) &&
                   bounds_overlap(a.minz, a.maxz, b.minz, b.maxz, eps) &&
                   bounds_overlap_length(a.miny, a.maxy, b.miny, b.maxy) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.minz, a.maxz, b.minz, b.maxz) >= FACE_BLOCK_MIN_OVERLAP;
        case 2: // +Y
            if (b.miny > a.maxy + eps || b.maxy < a.maxy - eps) return false;
            return bounds_overlap(a.minx, a.maxx, b.minx, b.maxx, eps) &&
                   bounds_overlap(a.minz, a.maxz, b.minz, b.maxz, eps) &&
                   bounds_overlap_length(a.minx, a.maxx, b.minx, b.maxx) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.minz, a.maxz, b.minz, b.maxz) >= FACE_BLOCK_MIN_OVERLAP;
        case 3: // -Y
            if (b.maxy < a.miny - eps || b.miny > a.miny + eps) return false;
            return bounds_overlap(a.minx, a.maxx, b.minx, b.maxx, eps) &&
                   bounds_overlap(a.minz, a.maxz, b.minz, b.maxz, eps) &&
                   bounds_overlap_length(a.minx, a.maxx, b.minx, b.maxx) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.minz, a.maxz, b.minz, b.maxz) >= FACE_BLOCK_MIN_OVERLAP;
        case 4: // +Z
            if (b.minz > a.maxz + eps || b.maxz < a.maxz - eps) return false;
            return bounds_overlap(a.minx, a.maxx, b.minx, b.maxx, eps) &&
                   bounds_overlap(a.miny, a.maxy, b.miny, b.maxy, eps) &&
                   bounds_overlap_length(a.minx, a.maxx, b.minx, b.maxx) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.miny, a.maxy, b.miny, b.maxy) >= FACE_BLOCK_MIN_OVERLAP;
        case 5: // -Z
            if (b.maxz < a.minz - eps || b.minz > a.minz + eps) return false;
            return bounds_overlap(a.minx, a.maxx, b.minx, b.maxx, eps) &&
                   bounds_overlap(a.miny, a.maxy, b.miny, b.maxy, eps) &&
                   bounds_overlap_length(a.minx, a.maxx, b.minx, b.maxx) >= FACE_BLOCK_MIN_OVERLAP &&
                   bounds_overlap_length(a.miny, a.maxy, b.miny, b.maxy) >= FACE_BLOCK_MIN_OVERLAP;
        default:
            return false;
    }
}

static void log_span2_face_cull(const Voxel *self, int self_idx,
                                const Voxel *neighbor, int neighbor_idx,
                                int face)
{
    if (!debugLogSpan2Faces || debugSpan2FaceLogBudget <= 0) {
        return;
    }
    if (!self || !neighbor) {
        return;
    }
    if (self->span != 2 && neighbor->span != 2) {
        return;
    }
    VoxelWorldBounds a, b;
    voxel_visibility_bounds(self, &a);
    voxel_visibility_bounds(neighbor, &b);
    TraceLog(LOG_INFO,
             "[FaceCull] face=%d self=%d span=%d neighbor=%d span=%d "
             "self=(%.2f,%.2f,%.2f)-(%.2f,%.2f,%.2f) "
             "neighbor=(%.2f,%.2f,%.2f)-(%.2f,%.2f,%.2f)",
             face, self_idx, self->span, neighbor_idx, neighbor->span,
             a.minx, a.miny, a.minz, a.maxx, a.maxy, a.maxz,
             b.minx, b.miny, b.minz, b.maxx, b.maxy, b.maxz);
    --debugSpan2FaceLogBudget;
}

static int voxel_touching_axes(const Voxel *a, const Voxel *b,
                               int *overlap_axes, float eps)
{
    if (!a || !b) {
        if (overlap_axes) *overlap_axes = 0;
        return 0;
    }

    VoxelWorldBounds boundsA, boundsB;
    voxel_world_bounds(a, &boundsA);
    voxel_world_bounds(b, &boundsB);

    bool touchX = false, overlapX = false;
    bool touchY = false, overlapY = false;
    bool touchZ = false, overlapZ = false;

    if (!axis_contact_state(boundsA.minx, boundsA.maxx,
                            boundsB.minx, boundsB.maxx, eps,
                            &touchX, &overlapX)) {
        if (overlap_axes) *overlap_axes = 0;
        return 0;
    }

    if (!axis_contact_state(boundsA.miny, boundsA.maxy,
                            boundsB.miny, boundsB.maxy, eps,
                            &touchY, &overlapY)) {
        if (overlap_axes) *overlap_axes = 0;
        return 0;
    }

    if (!axis_contact_state(boundsA.minz, boundsA.maxz,
                            boundsB.minz, boundsB.maxz, eps,
                            &touchZ, &overlapZ)) {
        if (overlap_axes) *overlap_axes = 0;
        return 0;
    }

    if (overlap_axes) {
        *overlap_axes = (overlapX ? 1 : 0) +
                        (overlapY ? 1 : 0) +
                        (overlapZ ? 1 : 0);
    }

    return (touchX ? 1 : 0) +
           (touchY ? 1 : 0) +
           (touchZ ? 1 : 0);
}

static bool ranges_touch_int(int minA, int maxA, int minB, int maxB) {
    return (maxA + 1 == minB) || (maxB + 1 == minA);
}

static bool ranges_overlap_int(int minA, int maxA, int minB, int maxB) {
    return !(maxA < minB || maxB < minA);
}

static bool voxels_share_edge_or_corner_rest(const Voxel *voxel_a, const Voxel *voxel_b) {
    int aMinX = voxel_a->rest_min_gx, aMaxX = voxel_a->rest_max_gx;
    int aMinY = voxel_a->rest_min_gy, aMaxY = voxel_a->rest_max_gy;
    int aMinZ = voxel_a->rest_min_gz, aMaxZ = voxel_a->rest_max_gz;
    int bMinX = voxel_b->rest_min_gx, bMaxX = voxel_b->rest_max_gx;
    int bMinY = voxel_b->rest_min_gy, bMaxY = voxel_b->rest_max_gy;
    int bMinZ = voxel_b->rest_min_gz, bMaxZ = voxel_b->rest_max_gz;

    bool touchX = ranges_touch_int(aMinX, aMaxX, bMinX, bMaxX);
    bool touchY = ranges_touch_int(aMinY, aMaxY, bMinY, bMaxY);
    bool touchZ = ranges_touch_int(aMinZ, aMaxZ, bMinZ, bMaxZ);

    bool overlapX = ranges_overlap_int(aMinX, aMaxX, bMinX, bMaxX);
    bool overlapY = ranges_overlap_int(aMinY, aMaxY, bMinY, bMaxY);
    bool overlapZ = ranges_overlap_int(aMinZ, aMaxZ, bMinZ, bMaxZ);

    if ((!touchX && !overlapX) ||
        (!touchY && !overlapY) ||
        (!touchZ && !overlapZ)) {
        return false;
    }

    int touching_axes = (touchX ? 1 : 0) + (touchY ? 1 : 0) + (touchZ ? 1 : 0);
    return touching_axes >= 2;
}

static bool voxels_share_face_rest(const Voxel *voxel_a, const Voxel *voxel_b) {
    int aMinX = voxel_a->rest_min_gx, aMaxX = voxel_a->rest_max_gx;
    int aMinY = voxel_a->rest_min_gy, aMaxY = voxel_a->rest_max_gy;
    int aMinZ = voxel_a->rest_min_gz, aMaxZ = voxel_a->rest_max_gz;
    int bMinX = voxel_b->rest_min_gx, bMaxX = voxel_b->rest_max_gx;
    int bMinY = voxel_b->rest_min_gy, bMaxY = voxel_b->rest_max_gy;
    int bMinZ = voxel_b->rest_min_gz, bMaxZ = voxel_b->rest_max_gz;

    bool overlapX = ranges_overlap_int(aMinX, aMaxX, bMinX, bMaxX);
    bool overlapY = ranges_overlap_int(aMinY, aMaxY, bMinY, bMaxY);
    bool overlapZ = ranges_overlap_int(aMinZ, aMaxZ, bMinZ, bMaxZ);

    bool faceX = overlapY && overlapZ &&
                 ((aMaxX + 1 == bMinX) || (bMaxX + 1 == aMinX));
    bool faceY = overlapX && overlapZ &&
                 ((aMaxY + 1 == bMinY) || (bMaxY + 1 == aMinY));
    bool faceZ = overlapX && overlapY &&
                 ((aMaxZ + 1 == bMinZ) || (bMaxZ + 1 == aMinZ));
    return faceX || faceY || faceZ;
}

static float voxel_rest_axis_min(const Voxel *v, int axis) {
    switch (axis) {
        case 0: return (float)v->rest_min_gx * VOXEL_SIZE;
        case 1: return (float)v->rest_min_gy * VOXEL_SIZE;
        default: return (float)v->rest_min_gz * VOXEL_SIZE;
    }
}

static float voxel_rest_axis_max(const Voxel *v, int axis) {
    switch (axis) {
        case 0: return (float)(v->rest_max_gx + 1) * VOXEL_SIZE;
        case 1: return (float)(v->rest_max_gy + 1) * VOXEL_SIZE;
        default: return (float)(v->rest_max_gz + 1) * VOXEL_SIZE;
    }
}

static float voxel_rest_corner_axis_coord(const Voxel *v, int axis, int corner_idx) {
    float min = voxel_rest_axis_min(v, axis);
    float max = voxel_rest_axis_max(v, axis);
    int sign = corner_signs[corner_idx][axis];
    return (sign >= 0) ? max : min;
}

static Vector3 voxel_rest_corner_world(const Voxel *v, int corner_idx) {
    return (Vector3){
        voxel_rest_corner_axis_coord(v, 0, corner_idx),
        voxel_rest_corner_axis_coord(v, 1, corner_idx),
        voxel_rest_corner_axis_coord(v, 2, corner_idx)
    };
}

static int table_get(int x, int y, int z);
static void rebuild_glue_constraints(void);
static void solve_voxel_glue(bool allow_break);
static bool voxels_share_edge_or_corner(const Voxel *voxel_a, const Voxel *voxel_b);
static bool voxels_share_face(const Voxel *voxel_a, const Voxel *voxel_b);

typedef struct {
    int dx, dy, dz;
    int faceA[4];
    int faceB[4];
} GlueDirection;

static const GlueDirection glueDirections[3] = {
    { 1, 0, 0, { 1, 3, 5, 7 }, { 0, 2, 4, 6 } },
    { 0, 1, 0, { 2, 3, 6, 7 }, { 0, 1, 4, 5 } },
    { 0, 0, 1, { 4, 5, 6, 7 }, { 0, 1, 2, 3 } },
};

static const char *glue_direction_label_from_delta(int dx, int dy, int dz)
{
    if (dx > 0) return "+X";
    if (dx < 0) return "-X";
    if (dy > 0) return "+Y";
    if (dy < 0) return "-Y";
    if (dz > 0) return "+Z";
    if (dz < 0) return "-Z";
    return "unknown";
}

static void voxel_compute_bounds(const Voxel *v,
                                 int *minx, int *maxx,
                                 int *miny, int *maxy,
                                 int *minz, int *maxz)
{
    int span = (v->span > 0) ? v->span : 1;
    if (span == 1) {
        if (minx) *minx = v->gx;
        if (maxx) *maxx = v->gx;
        if (miny) *miny = v->gy;
        if (maxy) *maxy = v->gy;
        if (minz) *minz = v->gz;
        if (maxz) *maxz = v->gz;
        return;
    }

    const float half = 0.5f * v->rest_edge;
    if (minx) *minx = (int)floorf((v->pos.x - half + GRID_EPSILON) / VOXEL_SIZE);
    if (maxx) *maxx = (int)floorf((v->pos.x + half - GRID_EPSILON) / VOXEL_SIZE);
    if (miny) *miny = (int)floorf((v->pos.y - half + GRID_EPSILON) / VOXEL_SIZE);
    if (maxy) *maxy = (int)floorf((v->pos.y + half - GRID_EPSILON) / VOXEL_SIZE);
    if (minz) *minz = (int)floorf((v->pos.z - half + GRID_EPSILON) / VOXEL_SIZE);
    if (maxz) *maxz = (int)floorf((v->pos.z + half - GRID_EPSILON) / VOXEL_SIZE);
}

static void voxel_grid_bounds(const Voxel *v,
                              int *minx, int *maxx,
                              int *miny, int *maxy,
                              int *minz, int *maxz)
{
    if (v->max_gx >= v->min_gx && v->max_gy >= v->min_gy && v->max_gz >= v->min_gz) {
        if (minx) *minx = v->min_gx;
        if (maxx) *maxx = v->max_gx;
        if (miny) *miny = v->min_gy;
        if (maxy) *maxy = v->max_gy;
        if (minz) *minz = v->min_gz;
        if (maxz) *maxz = v->max_gz;
        return;
    }

    voxel_compute_bounds(v, minx, maxx, miny, maxy, minz, maxz);
}


// Utility functions
static float randomInRange(float min, float max) {
    return min + ((float)rand()/ (float)RAND_MAX) * (max - min);
}
static float clampf(float v, float lo, float hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static float saturatef(float v) {
    return clampf(v, 0.0f, 1.0f);
}

static float mixf(float a, float b, float t) {
    return a * (1.0f - t) + b * t;
}
static Vector3 v_add(Vector3 a, Vector3 b) {
    return (Vector3){ a.x + b.x, a.y + b.y, a.z + b.z };
}
static Vector3 v_sub(Vector3 a, Vector3 b) {
    return (Vector3){ a.x - b.x, a.y - b.y, a.z - b.z };
}
static Vector3 v_mul(Vector3 v, float s) {
    return (Vector3){ v.x*s, v.y*s, v.z*s };
}
static float v_dot(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static Vector3 v_cross(Vector3 a, Vector3 b) {
    return (Vector3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static bool v_isfinite(Vector3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

static Particle *voxel_particle_at(Voxel *voxel, int idx) {
    return (idx == VOXEL_CENTER_INDEX) ? &voxel->center_particle : &voxel->particles[idx];
}

static float v_length(Vector3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

static Vector3 v_norm(Vector3 v) {
    float len = v_length(v);
    if (len == 0.0f) {
        // Return zero vector or handle however you want
        return (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    return v_mul(v, 1.0f / len);
}

static bool face_local_coords(const Vector3 origin,
                              const Vector3 U,
                              const Vector3 V,
                              const Vector3 normal,
                              float UU,
                              float VV,
                              float UV,
                              float invDet,
                              Vector3 point,
                              float *outU,
                              float *outV)
{
    Vector3 d = v_sub(point, origin);
    float signedDistance = v_dot(normal, d);
    Vector3 projected = v_sub(point, v_mul(normal, signedDistance));
    Vector3 p = v_sub(projected, origin);

    float du = v_dot(p, U);
    float dv = v_dot(p, V);

    float u = (du * VV - dv * UV) * invDet;
    float v = (dv * UU - du * UV) * invDet;
    if (outU) *outU = u;
    if (outV) *outV = v;
    return true;
}

// Patch for a merged quad in one of the principal planes
typedef struct {
    int   plane;     // 0=XY,1=XZ,2=YZ
    int   layer;     // z for XY, y for XZ, x for YZ
    int   i0, j0;    // lower-left grid cell in plane coords
    int   di, dj;    // width/height in voxels along in-plane axes
    bool  positive;  // true => positive-facing side (+Z,+Y,+X), else negative
    Color col;       // face color
    int   voxelIndex; // representative voxel driving this patch's debug color
} Patch;

static Patch patches[MAX_VOXELS];
static int   patchCount = 0;
static Mesh  greedyMesh = { 0 };
static Material greedyMaterial;
static bool greedyMaterialInit = false;
static bool  meshDirty  = true;

static Color voxel_belief_debug_color(const Voxel *voxel)
{
    if (!voxel) {
        return (Color){ 255, 255, 255, 255 };
    }
    float belief = voxel->simulate ? voxel->activationBelief : voxel->freezeBelief;
    belief = clampf(belief, 0.0f, 1.0f);
    float hue = voxel->simulate ? 5.0f : 210.0f;
    float saturation = 0.35f + 0.65f * belief;
    float value = 0.25f + 0.75f * belief;
    Color col = ColorFromHSV(hue, saturation, value);
    col.a = 255;
    return col;
}

static Color voxel_display_color(const Voxel *voxel)
{
    if (!voxel) {
        return (Color){ 160, 160, 160, 255 };
    }
    if (debugShowBeliefColors) {
        return voxel_belief_debug_color(voxel);
    }
    return voxel->color;
}

typedef struct {
    int      coarseVoxel;
    int      fineVoxel;
    int      coarseCorner[4];
    int      fineCornerFace[4];
    float    w[4];
    int      fineCorner;
    uint8_t  coarseMask;
    uint8_t  fineMask;
    float    restLocalU;
    float    restLocalV;
    float    restLocalN;
    float    restNormalAngle;
    float    strength;
    int      dirX;
    int      dirY;
    int      dirZ;
    bool     active;
} GlueConstraint;

static GlueConstraint glueConstraints[MAX_VOXELS * 48];
static int glueConstraintCount = 0;
static uint8_t gluedNeighborCounts[MAX_VOXELS];
static int gluedNeighborList[MAX_VOXELS][MAX_FACE_NEIGHBORS];
static uint16_t gluedNeighborRefCounts[MAX_VOXELS][MAX_FACE_NEIGHBORS];
static int gluedNeighborHashKeys[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE];
static uint8_t gluedNeighborHashIndex[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE];
static uint32_t gluedNeighborHashStamp[MAX_VOXELS][GLUE_NEIGHBOR_HASH_SIZE];
static uint32_t gluedNeighborEpoch[MAX_VOXELS];
static uint8_t glueAdjacencyDirtyFlags[MAX_VOXELS];
static int glueAdjacencyDirtyList[MAX_VOXELS];
static int glueAdjacencyDirtyCount = 0;
static bool glueAdjacencyDirtyAll = true;

static void mark_glue_adjacency_dirty_for_voxel(int voxel_idx) {
    if (glueAdjacencyDirtyAll) {
        return;
    }
    if (voxel_idx < 0 || voxel_idx >= MAX_VOXELS) {
        return;
    }
    if (glueAdjacencyDirtyFlags[voxel_idx]) {
        return;
    }
    glueAdjacencyDirtyFlags[voxel_idx] = 1;
    if (glueAdjacencyDirtyCount < MAX_VOXELS) {
        glueAdjacencyDirtyList[glueAdjacencyDirtyCount++] = voxel_idx;
    } else {
        glueAdjacencyDirtyAll = true;
        glueAdjacencyDirtyCount = 0;
    }
}

static void glue_neighbor_hash_reset(int voxel_idx) {
    uint32_t epoch = ++gluedNeighborEpoch[voxel_idx];
    if (epoch == 0) {
        memset(gluedNeighborHashStamp[voxel_idx], 0, sizeof(gluedNeighborHashStamp[voxel_idx]));
        gluedNeighborEpoch[voxel_idx] = 1;
    }
}

static unsigned glue_neighbor_hash_seed(int neighbor_idx) {
    return (unsigned)neighbor_idx * 2654435761u;
}

static bool glue_neighbor_hash_find(int voxel_idx, int neighbor_idx,
                                    int *slot_out, int *list_index_out) {
    uint32_t epoch = gluedNeighborEpoch[voxel_idx];
    if (epoch == 0) {
        glue_neighbor_hash_reset(voxel_idx);
        epoch = gluedNeighborEpoch[voxel_idx];
    }
    unsigned mask = GLUE_NEIGHBOR_HASH_SIZE - 1;
    unsigned h = glue_neighbor_hash_seed(neighbor_idx);
    for (unsigned probe = 0; probe < GLUE_NEIGHBOR_HASH_SIZE; ++probe) {
        unsigned slot = (h + probe) & mask;
        if (gluedNeighborHashStamp[voxel_idx][slot] != epoch) {
            if (slot_out) {
                *slot_out = (int)slot;
            }
            if (list_index_out) {
                *list_index_out = -1;
            }
            return false;
        }
        if (gluedNeighborHashKeys[voxel_idx][slot] == neighbor_idx) {
            if (slot_out) {
                *slot_out = (int)slot;
            }
            if (list_index_out) {
                *list_index_out = (int)gluedNeighborHashIndex[voxel_idx][slot];
            }
            return true;
        }
    }
    if (slot_out) {
        *slot_out = -1;
    }
    if (list_index_out) {
        *list_index_out = -1;
    }
    return false;
}

static bool glue_neighbor_hash_insert(int voxel_idx, int neighbor_idx, int list_index) {
    int slot = -1;
    int ignored = -1;
    if (glue_neighbor_hash_find(voxel_idx, neighbor_idx, &slot, &ignored)) {
        return true;
    }
    if (slot < 0) {
        return false;
    }
    gluedNeighborHashStamp[voxel_idx][slot] = gluedNeighborEpoch[voxel_idx];
    gluedNeighborHashKeys[voxel_idx][slot] = neighbor_idx;
    gluedNeighborHashIndex[voxel_idx][slot] = (uint8_t)list_index;
    return true;
}

static void glue_neighbor_hash_rebuild(int voxel_idx) {
    glue_neighbor_hash_reset(voxel_idx);
    int count = gluedNeighborCounts[voxel_idx];
    for (int i = 0; i < count; ++i) {
        int neighbor_idx = gluedNeighborList[voxel_idx][i];
        glue_neighbor_hash_insert(voxel_idx, neighbor_idx, i);
    }
}

static void glue_adjacency_reset_voxel(int voxel_idx) {
    gluedNeighborCounts[voxel_idx] = 0;
    for (int i = 0; i < MAX_FACE_NEIGHBORS; ++i) {
        gluedNeighborList[voxel_idx][i] = -1;
        gluedNeighborRefCounts[voxel_idx][i] = 0;
    }
    glue_neighbor_hash_reset(voxel_idx);
}

static void glue_adjacency_clear_all(void) {
    for (int i = 0; i < voxel_count; ++i) {
        glue_adjacency_reset_voxel(i);
    }
    memset(glueAdjacencyDirtyFlags, 0, sizeof(glueAdjacencyDirtyFlags));
    glueAdjacencyDirtyCount = 0;
    glueAdjacencyDirtyAll = false;
}

static void glue_adjacency_add_ref_oneway(int voxel_idx, int neighbor_idx, int ref_count) {
    if (voxel_idx < 0 || voxel_idx >= MAX_VOXELS ||
        neighbor_idx < 0 || neighbor_idx >= MAX_VOXELS) {
        return;
    }
    if (ref_count <= 0) {
        return;
    }
    int slot = -1;
    int list_index = -1;
    if (glue_neighbor_hash_find(voxel_idx, neighbor_idx, &slot, &list_index)) {
        if (list_index >= 0 && list_index < MAX_FACE_NEIGHBORS) {
            uint32_t sum = (uint32_t)gluedNeighborRefCounts[voxel_idx][list_index] +
                           (uint32_t)ref_count;
            if (sum > UINT16_MAX) {
                sum = UINT16_MAX;
            }
            gluedNeighborRefCounts[voxel_idx][list_index] = (uint16_t)sum;
        }
        return;
    }
    int count = gluedNeighborCounts[voxel_idx];
    if (count >= MAX_FACE_NEIGHBORS) {
        mark_glue_adjacency_dirty_for_voxel(voxel_idx);
        return;
    }
    gluedNeighborList[voxel_idx][count] = neighbor_idx;
    gluedNeighborRefCounts[voxel_idx][count] = (uint16_t)ref_count;
    gluedNeighborCounts[voxel_idx] = (uint8_t)(count + 1);
    if (!glue_neighbor_hash_insert(voxel_idx, neighbor_idx, count)) {
        mark_glue_adjacency_dirty_for_voxel(voxel_idx);
    }
}

static void glue_adjacency_add_ref_pair(int voxel_a, int voxel_b, int ref_count) {
    glue_adjacency_add_ref_oneway(voxel_a, voxel_b, ref_count);
    glue_adjacency_add_ref_oneway(voxel_b, voxel_a, ref_count);
}

static void glue_adjacency_remove_ref_oneway(int voxel_idx, int neighbor_idx, int ref_count) {
    if (voxel_idx < 0 || voxel_idx >= MAX_VOXELS ||
        neighbor_idx < 0 || neighbor_idx >= MAX_VOXELS) {
        return;
    }
    if (ref_count <= 0) {
        return;
    }
    int slot = -1;
    int list_index = -1;
    if (!glue_neighbor_hash_find(voxel_idx, neighbor_idx, &slot, &list_index)) {
        mark_glue_adjacency_dirty_for_voxel(voxel_idx);
        return;
    }
    if (list_index < 0 || list_index >= MAX_FACE_NEIGHBORS) {
        mark_glue_adjacency_dirty_for_voxel(voxel_idx);
        return;
    }
    uint16_t current = gluedNeighborRefCounts[voxel_idx][list_index];
    if (ref_count >= current) {
        int last = gluedNeighborCounts[voxel_idx] - 1;
        if (last < 0) {
            gluedNeighborCounts[voxel_idx] = 0;
        } else if (list_index != last) {
            gluedNeighborList[voxel_idx][list_index] = gluedNeighborList[voxel_idx][last];
            gluedNeighborRefCounts[voxel_idx][list_index] = gluedNeighborRefCounts[voxel_idx][last];
        }
        if (last >= 0) {
            gluedNeighborList[voxel_idx][last] = -1;
            gluedNeighborRefCounts[voxel_idx][last] = 0;
        }
        if (gluedNeighborCounts[voxel_idx] > 0) {
            gluedNeighborCounts[voxel_idx] = (uint8_t)last;
        }
        glue_neighbor_hash_rebuild(voxel_idx);
    } else {
        gluedNeighborRefCounts[voxel_idx][list_index] = (uint16_t)(current - ref_count);
    }
}

static void glue_adjacency_remove_ref_pair(int voxel_a, int voxel_b, int ref_count) {
    glue_adjacency_remove_ref_oneway(voxel_a, voxel_b, ref_count);
    glue_adjacency_remove_ref_oneway(voxel_b, voxel_a, ref_count);
}

static void rebuild_glue_adjacency_if_dirty(void) {
    if (!glueAdjacencyDirtyAll && glueAdjacencyDirtyCount == 0) {
        return;
    }
    bool rebuild_all = glueAdjacencyDirtyAll;
    if (rebuild_all) {
        glue_adjacency_clear_all();
    } else {
        for (int i = 0; i < glueAdjacencyDirtyCount; ++i) {
            int voxel_idx = glueAdjacencyDirtyList[i];
            glue_adjacency_reset_voxel(voxel_idx);
        }
        glueAdjacencyDirtyCount = 0;
    }
    for (int g = 0; g < glueConstraintCount; ++g) {
        const GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }
        if (rebuild_all) {
            glue_adjacency_add_ref_pair(gc->coarseVoxel, gc->fineVoxel, 1);
        } else {
            if (gc->coarseVoxel >= 0 && gc->coarseVoxel < MAX_VOXELS &&
                glueAdjacencyDirtyFlags[gc->coarseVoxel]) {
                glue_adjacency_add_ref_oneway(gc->coarseVoxel, gc->fineVoxel, 1);
            }
            if (gc->fineVoxel >= 0 && gc->fineVoxel < MAX_VOXELS &&
                glueAdjacencyDirtyFlags[gc->fineVoxel]) {
                glue_adjacency_add_ref_oneway(gc->fineVoxel, gc->coarseVoxel, 1);
            }
        }
    }
    if (!rebuild_all) {
        memset(glueAdjacencyDirtyFlags, 0, sizeof(glueAdjacencyDirtyFlags));
    }
    glueAdjacencyDirtyAll = false;
}

static bool face_normal_predicted(const Voxel *voxel, const int corners[4],
                                  Vector3 *out_normal) {
    Vector3 p0 = voxel->particles[corners[0]].predicted_pos;
    Vector3 p1 = voxel->particles[corners[1]].predicted_pos;
    Vector3 p2 = voxel->particles[corners[2]].predicted_pos;
    if (!v_isfinite(p0) || !v_isfinite(p1) || !v_isfinite(p2)) {
        return false;
    }
    Vector3 u = v_sub(p1, p0);
    Vector3 v = v_sub(p2, p0);
    Vector3 n = v_cross(u, v);
    float len = v_length(n);
    if (len < 1e-6f) {
        return false;
    }
    *out_normal = v_mul(n, 1.0f / len);
    return true;
}

static bool face_normal_rest(const Voxel *voxel, const int corners[4],
                             Vector3 *out_normal) {
    Vector3 p0 = voxel_rest_corner_world(voxel, corners[0]);
    Vector3 p1 = voxel_rest_corner_world(voxel, corners[1]);
    Vector3 p2 = voxel_rest_corner_world(voxel, corners[2]);
    Vector3 u = v_sub(p1, p0);
    Vector3 v = v_sub(p2, p0);
    Vector3 n = v_cross(u, v);
    float len = v_length(n);
    if (len < 1e-6f) {
        return false;
    }
    *out_normal = v_mul(n, 1.0f / len);
    return true;
}

static inline int voxel_span_for_glue(const Voxel *v) {
    return (v && v->span > 0) ? v->span : 1;
}

static void get_face_corners_for_direction(const GlueDirection *dir,
                                           bool positive_side,
                                           int outCorners[4])
{
    const int *src = positive_side ? dir->faceA : dir->faceB;
    for (int i = 0; i < 4; ++i) {
        outCorners[i] = src[i];
    }
}

static void order_coarse_fine_pair(int negativeIdx, int positiveIdx,
                                   int *coarseIdx, bool *coarseIsPositive,
                                   int *fineIdx, bool *fineIsPositive)
{
    const Voxel *neg = &voxels[negativeIdx];
    const Voxel *pos = &voxels[positiveIdx];
    int spanNeg = voxel_span_for_glue(neg);
    int spanPos = voxel_span_for_glue(pos);

    if (spanNeg >= spanPos) {
        *coarseIdx = negativeIdx;
        *fineIdx = positiveIdx;
    } else {
        *coarseIdx = positiveIdx;
        *fineIdx = negativeIdx;
    }

    *coarseIsPositive = (*coarseIdx == positiveIdx);
    *fineIsPositive = (*fineIdx == positiveIdx);
}



// Spatial hash table for voxels
//static struct { int key, idx; } table[HASH_SIZE];
// Spatial hash helpers
// static int hashVoxel(int x, int y, int z) {
//     unsigned int h = (unsigned int)(x*73856093 ^ y*19349663 ^ z*83492791);
//     return (int)(h & (HASH_SIZE - 1));
// }
// static void table_set(int x, int y, int z, int idx) {
//     int h = hashVoxel(x,y,z);
//     while (table[h].key) h = (h + 1) & (HASH_SIZE - 1);
//     table[h].key = 1;
//     table[h].idx = idx;
// }
// static int table_get(int x, int y, int z) {
//     int h = hashVoxel(x,y,z);
//     while (table[h].key) {
//         Voxel *v = &voxels[table[h].idx];
//         if (v->gx==x && v->gy==y && v->gz==z) return table[h].idx;
//         h = (h + 1) & (HASH_SIZE - 1);
//     }
//     return -1;
// }

/*----------------------------------------------------------*/
/*  Morton helpers                                           */
/*----------------------------------------------------------*/

static inline uint64_t part1by2(uint32_t n)
/* Interleave lower 21 bits of n with two 0-bits. */
{
    uint64_t x = n & 0x1FFFFFu;            // keep 21 bits, promote to 64-bit

    x = (x | (x << 32)) & 0x1F00000000FFFFULL;
    x = (x | (x << 16)) & 0x1F0000FF0000FFULL;
    x = (x | (x <<  8)) & 0x100F00F00F00F00FULL;
    x = (x | (x <<  4)) & 0x10C30C30C30C30C3ULL;
    x = (x | (x <<  2)) & 0x1249249249249249ULL;
    return x;
}

static inline uint64_t mortonKey(int x, int y, int z)
/* Pack signed coords into an **unsigned** Morton key.
 * We bias them by +2^20 so −1 048 575..+1 048 575 map into 0..2^21-1. */
{
    const uint32_t bias = 1u << 20;       // 1 048 576
    return  (part1by2((uint32_t)(x + bias))      ) |   // bits 0,3,6…
            (part1by2((uint32_t)(y + bias)) << 1 ) |   // bits 1,4,7…
            (part1by2((uint32_t)(z + bias)) << 2 );    // bits 2,5,8…
}

/*----------------------------------------------------------*/
/*  Hash table bucket & hash function                       */
/*----------------------------------------------------------*/

typedef struct {
    uint64_t key;   // 0 = empty
    int      idx;   // index in `voxels[]`
} Bucket;

static Bucket table[HASH_SIZE];
typedef struct {
    uint64_t key;
    int idx;
    uint32_t gen;
} TableCacheEntry;
static TableCacheEntry tableCache[TABLE_CACHE_SIZE];
static uint32_t tableCacheGeneration = 1;
static uint32_t tableCacheCursor = 0;
static int voxelHashFramesSinceRebuild = 0;

static inline void table_cache_invalidate(void)
{
    ++tableCacheGeneration;
    if (tableCacheGeneration == 0) {
        tableCacheGeneration = 1;
    }
}

static inline size_t hashVoxelKey(uint64_t k)
/*  SplitMix64 finalizer—fast avalanche, single multiply.
 *  HASH_SIZE **must** be a power of two so `&` is equivalent to % */
{
    k ^= k >> 30;  k *= 0xbf58476d1ce4e5b9ULL;
    k ^= k >> 27;  k *= 0x94d049bb133111ebULL;
    k ^= k >> 31;
    return (size_t)(k & (HASH_SIZE - 1));
}

/*----------------------------------------------------------*/
/*  Public helpers                                          */
/*----------------------------------------------------------*/

static int hashVoxel(int x, int y, int z)
/*  **Only** used by table_set and table_get; keep it for API parity. */
{
    return (int)hashVoxelKey(mortonKey(x, y, z));
}

static void table_set(int x, int y, int z, int idx)
{
    uint64_t k = mortonKey(x, y, z);
    size_t   h = hashVoxelKey(k);

    /* linear probe until we find an empty slot or an exact key match */
    while (table[h].key && table[h].key != k)
        h = (h + 1) & (HASH_SIZE - 1);

    table[h].key = k;      // stores the full Morton key
    table[h].idx = idx;
}

static int table_get(int x, int y, int z)
/*  Returns voxel index or –1 if empty */
{
    uint64_t k = mortonKey(x, y, z);
    for (int i = 0; i < TABLE_CACHE_SIZE; ++i) {
        if (tableCache[i].gen == tableCacheGeneration && tableCache[i].key == k) {
            return tableCache[i].idx;
        }
    }
    size_t   h = hashVoxelKey(k);

    while (1) {
        uint64_t bk = table[h].key;
        if (bk == 0) {
            tableCache[tableCacheCursor].key = k;
            tableCache[tableCacheCursor].idx = -1;
            tableCache[tableCacheCursor].gen = tableCacheGeneration;
            tableCacheCursor = (tableCacheCursor + 1) & (TABLE_CACHE_SIZE - 1);
            return -1;          // empty bucket ⇒ miss
        }
        if (bk == k) {
            int idx = table[h].idx;
            tableCache[tableCacheCursor].key = k;
            tableCache[tableCacheCursor].idx = idx;
            tableCache[tableCacheCursor].gen = tableCacheGeneration;
            tableCacheCursor = (tableCacheCursor + 1) & (TABLE_CACHE_SIZE - 1);
            return idx;  // exact key ⇒ hit
        }
        h = (h + 1) & (HASH_SIZE - 1);        // step to next bucket
    }
}
// Remove voxel entry from spatial hash and rehash subsequent cluster entries
static void table_remove(int x, int y, int z) {
    uint64_t k = mortonKey(x, y, z);
    size_t mask = HASH_SIZE - 1;
    size_t h = hashVoxelKey(k);
    // Locate the bucket for this key
    while (table[h].key) {
        if (table[h].key == k) break;
        h = (h + 1) & mask;
    }
    if (!table[h].key) return; // Key not found
    // Remove the entry
    table[h].key = 0;
    // Rehash any following entries in this probe cluster
    size_t j = (h + 1) & mask;
    while (table[j].key) {
        uint64_t k2 = table[j].key;
        int idx2 = table[j].idx;
        table[j].key = 0;
        // Find new home for this entry
        size_t h2 = hashVoxelKey(k2);
        while (table[h2].key) {
            h2 = (h2 + 1) & mask;
        }
        table[h2].key = k2;
        table[h2].idx = idx2;
        j = (j + 1) & mask;
    }
}

static void voxel_table_register(Voxel *v, int idx)
{
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_compute_bounds(v, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    v->min_gx = minx; v->max_gx = maxx;
    v->min_gy = miny; v->max_gy = maxy;
    v->min_gz = minz; v->max_gz = maxz;
    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            for (int z = minz; z <= maxz; ++z) {
                table_set(x, y, z, idx);
            }
        }
    }
    mark_glue_adjacency_dirty_for_voxel(idx);
}

static void voxel_table_unregister(const Voxel *v)
{
    if (!v) {
        return;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(v, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            for (int z = minz; z <= maxz; ++z) {
                table_remove(x, y, z);
            }
        }
    }
}

static void remove_voxel_index(int idx)
{
    if (idx < 0 || idx >= voxel_count) {
        return;
    }
    table_cache_invalidate();
    int voxel_count_before = voxel_count;
    Voxel *victim = &voxels[idx];
    if (!victim->simulate) {
        mark_static_beliefs_dirty_for_voxel(victim);
    }
    voxel_table_unregister(victim);

    int last = voxel_count - 1;
    if (idx != last) {
        Voxel *moved = &voxels[last];
        voxel_table_unregister(moved);
        voxels[idx] = *moved;
        voxel_table_register(&voxels[idx], idx);
        if (!voxels[idx].simulate) {
            mark_static_beliefs_dirty_for_voxel(&voxels[idx]);
        }
    }
    --voxel_count;
    if (debugLogVoxelDeactivation) {
        TraceLog(LOG_INFO,
                 "[RemoveVoxel] idx=%d voxel_count_before=%d moved_from=%d voxel_count_after=%d",
                 idx, voxel_count_before, (idx != last) ? last : -1, voxel_count);
    }
}

static bool occupied(int x, int y, int z) {
    return table_get(x, y, z) >= 0;
}

static int add_static_voxel_at_grid(int gx, int gy, int gz, Color color, int type)
{
    float px = ((float)gx + 0.5f) * VOXEL_SIZE;
    float py = ((float)gy + 0.5f) * VOXEL_SIZE;
    float pz = ((float)gz + 0.5f) * VOXEL_SIZE;
    return addVoxel(px, py, pz, true, false, color, type);
}

static void tag_owned_static_voxel(int idx, int owner)
{
    if (idx < 0 || idx >= voxel_count) {
        return;
    }
    Voxel *voxel = &voxels[idx];
    if (voxel->simulate) {
        return;
    }
    voxel->owner = owner;
    voxel->lifeFrames = 0;
}

static bool grid_region_is_free(int minx, int maxx,
                                int miny, int maxy,
                                int minz, int maxz)
{
    if (minx > maxx || miny > maxy || minz > maxz) {
        return false;
    }
    if (miny < 0) {
        return false;
    }
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx >= voxel_count) {
                    continue;
                }
                Voxel *blocker = &voxels[idx];
                if (blocker->simulate) {
                    return false;
                }
            }
        }
    }
    return true;
}

static void remove_static_voxels_in_region(int minx, int maxx,
                                           int miny, int maxy,
                                           int minz, int maxz)
{
    if (minx > maxx || miny > maxy || minz > maxz) {
        return;
    }
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                while (1) {
                    int idx = table_get(x, y, z);
                    if (idx < 0 || idx >= voxel_count) {
                        break;
                    }
                    Voxel *candidate = &voxels[idx];
                    if (candidate->simulate) {
                        break;
                    }
                    remove_voxel_index(idx);
                }
            }
        }
    }
}

static void remove_static_voxels_in_region_recycle(int minx, int maxx,
                                                   int miny, int maxy,
                                                   int minz, int maxz)
{
    if (minx > maxx || miny > maxy || minz > maxz) {
        return;
    }
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                while (1) {
                    int idx = table_get(x, y, z);
                    if (idx < 0 || idx >= voxel_count) {
                        break;
                    }
                    Voxel *candidate = &voxels[idx];
                    if (candidate->simulate) {
                        break;
                    }
                    if (candidate->owner == -1) {
                        if (recycle_queue_push(candidate) && debugLogVoxelRecycle) {
                            TraceLog(LOG_INFO,
                                     "[Recycle] enqueue-bullet voxel=%d span=%d rest=(%d..%d,%d..%d,%d..%d)",
                                     idx, candidate->span,
                                     candidate->rest_min_gx, candidate->rest_max_gx,
                                     candidate->rest_min_gy, candidate->rest_max_gy,
                                     candidate->rest_min_gz, candidate->rest_max_gz);
                        }
                    }
                    remove_voxel_index(idx);
                }
            }
        }
    }
}

static void remove_unowned_static_voxels_in_region(int minx, int maxx,
                                                   int miny, int maxy,
                                                   int minz, int maxz)
{
    if (minx > maxx || miny > maxy || minz > maxz) {
        return;
    }
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                while (1) {
                    int idx = table_get(x, y, z);
                    if (idx < 0 || idx >= voxel_count) {
                        break;
                    }
                    Voxel *candidate = &voxels[idx];
                    if (candidate->simulate) {
                        break;
                    }
                    if (candidate->owner >= 0) {
                        break;
                    }
                    remove_voxel_index(idx);
                }
            }
        }
    }
}

static void remove_debris_static_in_columns(int minx, int maxx,
                                            int minz, int maxz,
                                            int miny, int maxy)
{
    if (minx > maxx || minz > maxz || miny > maxy) {
        return;
    }
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                while (1) {
                    int idx = table_get(x, y, z);
                    if (idx < 0 || idx >= voxel_count) {
                        break;
                    }
                    Voxel *candidate = &voxels[idx];
                    if (candidate->simulate) {
                        break;
                    }
                    if (candidate->owner != STATIC_DEBRIS_OWNER) {
                        break;
                    }
                    remove_voxel_index(idx);
                }
            }
        }
    }
}

static bool remove_all_debris_static(void)
{
    bool removed = false;
    int i = 0;
    while (i < voxel_count) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate && voxel->owner == STATIC_DEBRIS_OWNER) {
            remove_voxel_index(i);
            removed = true;
            continue;
        }
        ++i;
    }
    return removed;
}

static bool find_nearest_free_static_region(int base_minx, int base_maxx,
                                            int base_miny, int base_maxy,
                                            int base_minz, int base_maxz,
                                            int min_y_limit,
                                            int *out_minx, int *out_maxx,
                                            int *out_miny, int *out_maxy,
                                            int *out_minz, int *out_maxz)
{
    if (!out_minx || !out_maxx || !out_miny || !out_maxy || !out_minz || !out_maxz) {
        return false;
    }
    if (base_miny < min_y_limit) {
        base_maxy += (min_y_limit - base_miny);
        base_miny = min_y_limit;
    }
    if (grid_region_is_free(base_minx, base_maxx,
                            base_miny, base_maxy,
                            base_minz, base_maxz))
    {
        *out_minx = base_minx; *out_maxx = base_maxx;
        *out_miny = base_miny; *out_maxy = base_maxy;
        *out_minz = base_minz; *out_maxz = base_maxz;
        return true;
    }

    const int search = STATIC_RESTORE_SEARCH_RADIUS;
    for (int radius = 1; radius <= search; ++radius) {
        bool found = false;
        int best_dx = 0, best_dy = 0, best_dz = 0;
        int best_dist_sq = INT_MAX;
        int limit = radius;
        for (int dz = -limit; dz <= limit; ++dz) {
            for (int dy = -limit; dy <= limit; ++dy) {
                for (int dx = -limit; dx <= limit; ++dx) {
                    if (dx == 0 && dy == 0 && dz == 0) {
                        continue;
                    }
                    int dist_sq = dx*dx + dy*dy + dz*dz;
                    if (dist_sq > radius * radius) {
                        continue;
                    }
                    int test_minx = base_minx + dx;
                    int test_maxx = base_maxx + dx;
                    int test_miny = base_miny + dy;
                    int test_maxy = base_maxy + dy;
                    int test_minz = base_minz + dz;
                    int test_maxz = base_maxz + dz;
                    if (test_miny < min_y_limit) {
                        continue;
                    }
                    if (!grid_region_is_free(test_minx, test_maxx,
                                             test_miny, test_maxy,
                                             test_minz, test_maxz)) {
                        continue;
                    }
                    if (dist_sq < best_dist_sq) {
                        best_dist_sq = dist_sq;
                        best_dx = dx;
                        best_dy = dy;
                        best_dz = dz;
                        found = true;
                    }
                }
            }
        }

        if (found) {
            *out_minx = base_minx + best_dx;
            *out_maxx = base_maxx + best_dx;
            *out_miny = base_miny + best_dy;
            *out_maxy = base_maxy + best_dy;
            *out_minz = base_minz + best_dz;
            *out_maxz = base_maxz + best_dz;
            return true;
        }
    }

    return false;
}

static bool spawn_static_covering_voxel(const Voxel *voxel)
{
    if (!voxel) {
        TraceLog(LOG_WARNING, "[Multiscale] Static restore failed: null voxel");
        return false;
    }
    int base_minx, base_maxx, base_miny, base_maxy, base_minz, base_maxz;
    bool has_rest_bounds =
        (voxel->rest_min_gx <= voxel->rest_max_gx) &&
        (voxel->rest_min_gy <= voxel->rest_max_gy) &&
        (voxel->rest_min_gz <= voxel->rest_max_gz);

    if (has_rest_bounds) {
        float rest_center_x = 0.5f * ((float)voxel->rest_min_gx + (float)voxel->rest_max_gx + 1.0f);
        float rest_center_y = 0.5f * ((float)voxel->rest_min_gy + (float)voxel->rest_max_gy + 1.0f);
        float rest_center_z = 0.5f * ((float)voxel->rest_min_gz + (float)voxel->rest_max_gz + 1.0f);
        float current_center_x = voxel->pos.x / VOXEL_SIZE;
        float current_center_y = voxel->pos.y / VOXEL_SIZE;
        float current_center_z = voxel->pos.z / VOXEL_SIZE;

        int shift_x = (int)roundf(current_center_x - rest_center_x);
        int shift_y = (int)roundf(current_center_y - rest_center_y);
        int shift_z = (int)roundf(current_center_z - rest_center_z);

        base_minx = voxel->rest_min_gx + shift_x;
        base_maxx = voxel->rest_max_gx + shift_x;
        base_miny = voxel->rest_min_gy + shift_y;
        base_maxy = voxel->rest_max_gy + shift_y;
        base_minz = voxel->rest_min_gz + shift_z;
        base_maxz = voxel->rest_max_gz + shift_z;
    } else {
        voxel_grid_bounds(voxel, &base_minx, &base_maxx, &base_miny, &base_maxy, &base_minz, &base_maxz);
    }

    if (base_minx > base_maxx || base_miny > base_maxy || base_minz > base_maxz) {
        TraceLog(LOG_WARNING,
                 "[Multiscale] Static restore failed: invalid bounds (%d..%d,%d..%d,%d..%d)",
                 base_minx, base_maxx, base_miny, base_maxy, base_minz, base_maxz);
        return false;
    }

    int minx = base_minx;
    int maxx = base_maxx;
    int miny = base_miny;
    int maxy = base_maxy;
    int minz = base_minz;
    int maxz = base_maxz;

    bool placed = find_nearest_free_static_region(base_minx, base_maxx,
                                                  base_miny, base_maxy,
                                                  base_minz, base_maxz,
                                                  base_miny,
                                                  &minx, &maxx,
                                                  &miny, &maxy,
                                                  &minz, &maxz);
    if (!placed && has_rest_bounds) {
        int curr_minx, curr_maxx, curr_miny, curr_maxy, curr_minz, curr_maxz;
        voxel_grid_bounds(voxel, &curr_minx, &curr_maxx, &curr_miny, &curr_maxy, &curr_minz, &curr_maxz);
        placed = find_nearest_free_static_region(curr_minx, curr_maxx,
                                                 curr_miny, curr_maxy,
                                                 curr_minz, curr_maxz,
                                                 curr_miny,
                                                 &minx, &maxx,
                                                 &miny, &maxy,
                                                 &minz, &maxz);
    }
    if (!placed) {
        TraceLog(LOG_WARNING,
                 "[Multiscale] Static restore failed: no free cells near (%.2f, %.2f, %.2f)",
                 voxel->pos.x, voxel->pos.y, voxel->pos.z);
        return false;
    }

    int desired_span = (voxel->span > 0) ? voxel->span : 1;
    if (desired_span > 1) {
        float half_span = 0.5f * (float)desired_span;
        float center_x_cells = voxel->pos.x / VOXEL_SIZE;
        float center_y_cells = voxel->pos.y / VOXEL_SIZE;
        float center_z_cells = voxel->pos.z / VOXEL_SIZE;

        int width_x = maxx - minx + 1;
        if (width_x != desired_span) {
            int forced_minx = (int)floorf(center_x_cells - half_span + GRID_EPSILON);
            minx = forced_minx;
            maxx = forced_minx + desired_span - 1;
        }

        int width_y = maxy - miny + 1;
        if (width_y != desired_span) {
            int forced_miny = (int)floorf(center_y_cells - half_span + GRID_EPSILON);
            miny = forced_miny;
            maxy = forced_miny + desired_span - 1;
            if (miny < 0) {
                maxy += -miny;
                miny = 0;
            }
        }

        int width_z = maxz - minz + 1;
        if (width_z != desired_span) {
            int forced_minz = (int)floorf(center_z_cells - half_span + GRID_EPSILON);
            minz = forced_minz;
            maxz = forced_minz + desired_span - 1;
        }
    }

    int span_count_x = maxx - minx + 1;
    int span_count_y = maxy - miny + 1;
    int span_count_z = maxz - minz + 1;
    if (span_count_x <= 0 || span_count_y <= 0 || span_count_z <= 0) {
        TraceLog(LOG_WARNING,
                 "[Multiscale] Static restore failed: empty span (%d,%d,%d)",
                 span_count_x, span_count_y, span_count_z);
        return false;
    }
    int cell_count = span_count_x * span_count_y * span_count_z;
    if (cell_count <= 0) {
        TraceLog(LOG_WARNING,
                 "[Multiscale] Static restore failed: invalid cell count %d",
                 cell_count);
        return false;
    }
    if (voxel_count + cell_count > MAX_VOXELS) {
        TraceLog(LOG_WARNING, "[Multiscale] Static restore skipped: insufficient capacity (%d needed)", cell_count);
        return false;
    }

    remove_static_voxels_in_region(minx, maxx, miny, maxy, minz, maxz);

    for (int gz = minz; gz <= maxz; ++gz) {
        for (int gy = miny; gy <= maxy; ++gy) {
            for (int gx = minx; gx <= maxx; ++gx) {
                int idx = add_static_voxel_at_grid(gx, gy, gz, voxel->color, voxel->type);
                if (idx < 0) {
                    TraceLog(LOG_WARNING,
                             "[Multiscale] Static restore aborted: failed to spawn voxel at (%d,%d,%d)",
                             gx, gy, gz);
                    remove_static_voxels_in_region(minx, maxx, miny, maxy, minz, maxz);
                    return false;
                }
                if (voxel->owner == -1) {
                    voxels[idx].owner = STATIC_DEBRIS_OWNER;
                    voxels[idx].lifeFrames = 0;
                } else {
                    voxels[idx].owner = voxel->owner;
                    voxels[idx].lifeFrames = 0;
                }
            }
        }
    }
    return true;
}


static void mark_surface(int idx) {
    Voxel *v = &voxels[idx];
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(v, &minx, &maxx, &miny, &maxy, &minz, &maxz);

    v->surface[0] = true; // +X
    for (int y = miny; y <= maxy && v->surface[0]; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            if (occupied(maxx + 1, y, z)) {
                v->surface[0] = false;
                break;
            }
        }
    }

    v->surface[1] = true; // -X
    for (int y = miny; y <= maxy && v->surface[1]; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            if (occupied(minx - 1, y, z)) {
                v->surface[1] = false;
                break;
            }
        }
    }

    v->surface[2] = true; // +Y
    for (int x = minx; x <= maxx && v->surface[2]; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            if (occupied(x, maxy + 1, z)) {
                v->surface[2] = false;
                break;
            }
        }
    }

    v->surface[3] = true; // -Y
    for (int x = minx; x <= maxx && v->surface[3]; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            if (occupied(x, miny - 1, z)) {
                v->surface[3] = false;
                break;
            }
        }
    }

    v->surface[4] = true; // +Z
    for (int x = minx; x <= maxx && v->surface[4]; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            if (occupied(x, y, maxz + 1)) {
                v->surface[4] = false;
                break;
            }
        }
    }

    v->surface[5] = true; // -Z
    for (int x = minx; x <= maxx && v->surface[5]; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            if (occupied(x, y, minz - 1)) {
                v->surface[5] = false;
                break;
            }
        }
    }
}
// Mark surfaces of voxels adjacent to a given world position
static void mark_surface_neighbors(Vector3 pos) {
    // Compute grid coordinates of the given position
    int gx = (int)floorf(pos.x / VOXEL_SIZE);
    int gy = (int)floorf(pos.y / VOXEL_SIZE);
    int gz = (int)floorf(pos.z / VOXEL_SIZE);

    // Offsets for the 6 face-adjacent neighbors
    const int offs[6][3] = {
        { 1,  0,  0 }, { -1,  0,  0 },
        { 0,  1,  0 }, {  0, -1,  0 },
        { 0,  0,  1 }, {  0,  0, -1 }
    };

    // Lookup each neighbor and mark its surface
    for (int i = 0; i < 6; i++) {
        int nx = gx + offs[i][0];
        int ny = gy + offs[i][1];
        int nz = gz + offs[i][2];
        int neighbor_idx = table_get(nx, ny, nz);
        if (neighbor_idx >= 0) {
            mark_surface(neighbor_idx);
        }
    }
}

static void rebuild_all_voxel_surfaces(void) {
    for (int i = 0; i < voxel_count; ++i) {
        mark_surface(i);
    }
}
//-----------------------------------------------------------------------------
// Check for voxels adjacent to a given world position
// neighbors[0] = +X, neighbors[1] = -X,
// neighbors[2] = +Y, neighbors[3] = -Y,
// neighbors[4] = +Z, neighbors[5] = -Z
static void get_adjacent_voxel_directions(Vector3 pos, bool neighbors[6]) {
    // Compute grid coordinates of the position
    int gx = (int)floorf(pos.x / VOXEL_SIZE);
    int gy = (int)floorf(pos.y / VOXEL_SIZE-VOXEL_SIZE);
    int gz = (int)floorf(pos.z / VOXEL_SIZE);
    // Offsets for the 6 face-adjacent neighbors
    const int offs[6][3] = {
        { 1,  0,  0 }, { -1,  0,  0 },
        { 0,  1,  0 }, {  0, -1,  0 },
        { 0,  0,  1 }, {  0,  0, -1 }
    };
    // Test occupancy at each neighboring cell
    for (int i = 0; i < 6; i++) {
        neighbors[i] = occupied(
            gx + offs[i][0],
            gy + offs[i][1],
            gz + offs[i][2]
        );
    }
}

static void init_voxel_struct(Voxel *v,
                              float px, float py, float pz,
                              bool fixed, bool simulate,
                              Color color, int type,
                              int span, int owner)
{
    if (!v) {
        return;
    }
    if (span < 1) span = 1;
    float edge = VOXEL_SIZE * (float)span;
    float half = 0.5f * edge;

    v->pos = (Vector3){ px, py, pz };
    v->vel = (Vector3){ 0,0,0 };
    v->fixed = fixed;
    v->simulate = simulate;
    v->glueEligible = true;
    v->pendingActivation = false;
    v->isBullet = false;
    v->activationCooldownFrames = 0;
    v->color = color;
    v->type = type;
    v->owner = owner;
    v->activator = -1;
    memset(v->surface, 0, sizeof v->surface);
    v->min_gx = v->min_gy = v->min_gz = INT_MAX;
    v->max_gx = v->max_gy = v->max_gz = INT_MIN;
    v->rest_min_gx = v->rest_min_gy = v->rest_min_gz = INT_MAX;
    v->rest_max_gx = v->rest_max_gy = v->rest_max_gz = INT_MIN;
    v->orig_min_gx = v->orig_min_gy = v->orig_min_gz = INT_MAX;
    v->orig_max_gx = v->orig_max_gy = v->orig_max_gz = INT_MIN;
    v->gx = (int)floorf(px / VOXEL_SIZE);
    v->gy = (int)floorf(py / VOXEL_SIZE);
    v->gz = (int)floorf(pz / VOXEL_SIZE);
    v->span = span;
    v->rest_edge = edge;
    v->rest_volume = edge * edge * edge;
    v->particle_radius = 0.5f * edge;
    for (int i = 0; i < VOXEL_CORNER_COUNT; ++i) {
        Particle *p = &v->particles[i];
        p->pos = (Vector3){
            px + corner_signs[i][0] * half,
            py + corner_signs[i][1] * half,
            pz + corner_signs[i][2] * half
        };
        p->prev_pos = p->pos;
        p->predicted_pos = p->pos;
        p->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
        if (fixed || !simulate) {
            p->inv_mass = 0.0f;
        } else {
            float mass_scale = (float)(span * span * span);
            if (mass_scale <= 0.0f) mass_scale = 1.0f;
            p->inv_mass = 1.0f / mass_scale;
        }
    }
    v->center_particle.pos = v->pos;
    v->center_particle.prev_pos = v->pos;
    v->center_particle.predicted_pos = v->pos;
    v->center_particle.vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    if (fixed || !simulate) {
        v->center_particle.inv_mass = 0.0f;
    } else {
        float mass_scale = (float)(span * span * span);
        if (mass_scale <= 0.0f) mass_scale = 1.0f;
        v->center_particle.inv_mass = 1.0f / mass_scale;
    }
    int rest_minx, rest_maxx, rest_miny, rest_maxy, rest_minz, rest_maxz;
    voxel_compute_bounds(v, &rest_minx, &rest_maxx, &rest_miny, &rest_maxy, &rest_minz, &rest_maxz);
    v->rest_min_gx = rest_minx;
    v->rest_max_gx = rest_maxx;
    v->rest_min_gy = rest_miny;
    v->rest_max_gy = rest_maxy;
    v->rest_min_gz = rest_minz;
    v->rest_max_gz = rest_maxz;
    v->orig_min_gx = rest_minx;
    v->orig_max_gx = rest_maxx;
    v->orig_min_gy = rest_miny;
    v->orig_max_gy = rest_maxy;
    v->orig_min_gz = rest_minz;
    v->orig_max_gz = rest_maxz;
    v->sleepFrames = 0;
    v->freezeBelief = 0.0f;
    v->activationBelief = 0.0f;
    v->groundSupport = 0.0f;
    v->neighborSupport = 0;
    v->supportMask = 0;
    v->skipCollisionVelocityFrames = 0;
    v->lifeFrames = 0;
    v->debugClusterTag = 0;
    v->prevGlueClusterId = -1;
    v->prevGlueClusterSize = 0;
    v->prevGlueClusterTag = 0;
    v->prevGlueClusterValid = 0;
}

static int addVoxelSized(float px, float py, float pz, bool fixed, bool simulate,
                         Color color, int type, int span) {
    if (voxel_count >= MAX_VOXELS) {
        if (debugLogActivationFailures) {
            TraceLog(LOG_WARNING,
                     "[ActivationFail] addVoxelSized capacity full simulate=%d span=%d pos=(%.2f,%.2f,%.2f)",
                     simulate ? 1 : 0, span, px, py, pz);
        }
        return -1;
    }
    int idx = voxel_count++;
    Voxel *v = &voxels[idx];
    init_voxel_struct(v, px, py, pz, fixed, simulate, color, type, span, -1);
    table_cache_invalidate();
    voxel_table_register(v, idx);
    if (!simulate) {
        mark_static_beliefs_dirty_for_voxel(v);
    }
    return idx;
}

// Add a voxel (static or dynamic)
static int addVoxel(float px, float py, float pz, bool fixed, bool simulate, Color color, int type) {
    int idx = addVoxelSized(px, py, pz, fixed, simulate, color, type, 1);
    if (idx >= 0) {
        voxels[idx].glueEligible = false;
    }
    return idx;
}

static inline size_t unit_voxel_grid_index(int x, int y, int z,
                                           int dimx, int dimy, int dimz)
{
    (void)dimz;
    return ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
}

static int maximal_cube_span_at(int sx, int sy, int sz,
                                const unsigned char *occupied,
                                const unsigned char *consumed,
                                int dimx, int dimy, int dimz)
{
    int max_span = dimx - sx;
    int limit_y = dimy - sy;
    int limit_z = dimz - sz;
    if (limit_y < max_span) max_span = limit_y;
    if (limit_z < max_span) max_span = limit_z;
    if (max_span < 1) {
        return 1;
    }

    int best_span = 1;
    for (int span = 1; span <= max_span; ++span) {
        bool ok = true;
        for (int dz = 0; dz < span && ok; ++dz) {
            for (int dy = 0; dy < span && ok; ++dy) {
                for (int dx = 0; dx < span; ++dx) {
                    size_t idx = unit_voxel_grid_index(sx + dx, sy + dy, sz + dz,
                                                       dimx, dimy, dimz);
                    if (!occupied[idx] || consumed[idx]) {
                        ok = false;
                        break;
                    }
                }
            }
        }
        if (ok) {
            best_span = span;
        } else {
            break;
        }
    }
    return best_span;
}

static int emit_multiscale_voxels_from_units(const UnitVoxelBuffer *buffer,
                                             bool fixed, bool simulate,
                                             int type_override)
{
    if (!buffer || buffer->count <= 0) {
        return 0;
    }

    int minx = INT_MAX, miny = INT_MAX, minz = INT_MAX;
    int maxx = INT_MIN, maxy = INT_MIN, maxz = INT_MIN;
    for (int i = 0; i < buffer->count; ++i) {
        const UnitVoxelSeed *seed = &buffer->voxels[i];
        if (seed->gx < minx) minx = seed->gx;
        if (seed->gy < miny) miny = seed->gy;
        if (seed->gz < minz) minz = seed->gz;
        if (seed->gx > maxx) maxx = seed->gx;
        if (seed->gy > maxy) maxy = seed->gy;
        if (seed->gz > maxz) maxz = seed->gz;
    }

    if (minx > maxx || miny > maxy || minz > maxz) {
        return 0;
    }

    int dimx = maxx - minx + 1;
    int dimy = maxy - miny + 1;
    int dimz = maxz - minz + 1;

    size_t cell_count = (size_t)dimx * (size_t)dimy * (size_t)dimz;
    if (cell_count == 0) {
        return 0;
    }

    unsigned char *occupied = (unsigned char *)calloc(cell_count, sizeof(unsigned char));
    unsigned char *consumed = (unsigned char *)calloc(cell_count, sizeof(unsigned char));
    Color *color_grid = (Color *)malloc(cell_count * sizeof(Color));
    int *type_grid = (int *)malloc(cell_count * sizeof(int));
    int *tag_grid = (int *)malloc(cell_count * sizeof(int));
    int *activator_grid = (int *)malloc(cell_count * sizeof(int));
    if (!occupied || !consumed || !color_grid || !type_grid || !tag_grid || !activator_grid) {
        free(occupied);
        free(consumed);
        free(color_grid);
        free(type_grid);
        free(tag_grid);
        free(activator_grid);
        TraceLog(LOG_WARNING, "[Multiscale] Failed to allocate voxel mask (%dx%dx%d)", dimx, dimy, dimz);
        return 0;
    }

    for (int i = 0; i < buffer->count; ++i) {
        const UnitVoxelSeed *seed = &buffer->voxels[i];
        int local_x = seed->gx - minx;
        int local_y = seed->gy - miny;
        int local_z = seed->gz - minz;
        size_t idx = unit_voxel_grid_index(local_x, local_y, local_z,
                                           dimx, dimy, dimz);
        occupied[idx] = 1;
        color_grid[idx] = seed->color;
        type_grid[idx] = (type_override >= 0) ? type_override : seed->type;
        tag_grid[idx] = seed->debugTag;
        activator_grid[idx] = seed->activator;
    }

    int spawned = 0;
    for (int z = 0; z < dimz; ++z) {
        for (int y = 0; y < dimy; ++y) {
            for (int x = 0; x < dimx; ++x) {
                size_t cell_idx = unit_voxel_grid_index(x, y, z, dimx, dimy, dimz);
                if (!occupied[cell_idx] || consumed[cell_idx]) {
                    continue;
                }

                int span = maximal_cube_span_at(x, y, z, occupied, consumed, dimx, dimy, dimz);
                Color color = color_grid[cell_idx];
                int spawn_type = type_grid[cell_idx];
                int spawn_tag = tag_grid[cell_idx];
                int spawn_activator = activator_grid[cell_idx];

                for (int dz = 0; dz < span; ++dz) {
                    for (int dy = 0; dy < span; ++dy) {
                        for (int dx = 0; dx < span; ++dx) {
                            size_t idx = unit_voxel_grid_index(x + dx, y + dy, z + dz,
                                                               dimx, dimy, dimz);
                            consumed[idx] = 1;
                        }
                    }
                }

                int gx = minx + x;
                int gy = miny + y;
                int gz = minz + z;
                float px = ((float)gx + 0.5f * (float)span) * VOXEL_SIZE;
                float py = ((float)gy + 0.5f * (float)span) * VOXEL_SIZE;
                float pz = ((float)gz + 0.5f * (float)span) * VOXEL_SIZE;
                int new_idx = addVoxelSized(px, py, pz, fixed, simulate, color, spawn_type, span);
                if (new_idx >= 0) {
                    voxels[new_idx].debugClusterTag = spawn_tag;
                    voxels[new_idx].activator = spawn_activator;
                    if (debugLogSmush && debugLogSmushSpawns &&
                        spawn_activator >= 0 && debugSmushLogBudget > 0) {
                        TraceLog(LOG_INFO,
                                 "[Smush] spawn dynamic voxel=%d span=%d activator=%d",
                                 new_idx, span, spawn_activator);
                        --debugSmushLogBudget;
                    }
                } else if (debugLogSmush && debugLogSmushSpawns && debugSmushLogBudget > 0) {
                    TraceLog(LOG_INFO,
                             "[Smush] spawn failed span=%d activator=%d",
                             span, spawn_activator);
                    --debugSmushLogBudget;
                }
                if (debugLogMultiscale) {
                    TraceLog(LOG_INFO,
                             "[Multiscale] Emit span=%d grid=(%d..%d,%d..%d,%d..%d) pos=(%.2f,%.2f,%.2f)",
                             span,
                             gx, gx + span - 1,
                             gy, gy + span - 1,
                             gz, gz + span - 1,
                             px, py, pz);
                }
                ++spawned;
            }
        }
    }

    free(occupied);
    free(consumed);
    free(color_grid);
    free(type_grid);
    free(tag_grid);
    free(activator_grid);

    if (debugLogMultiscale) {
        TraceLog(LOG_INFO,
                 "[Multiscale] Converted %d unit voxels into %d span clusters",
                 buffer->count, spawned);
    }
    return spawned;
}

static void emit_static_voxels_from_units(const UnitVoxelBuffer *buffer)
{
    if (!buffer) {
        return;
    }
    for (int i = 0; i < buffer->count; ++i) {
        const UnitVoxelSeed *seed = &buffer->voxels[i];
        int idx = add_static_voxel_at_grid(seed->gx, seed->gy, seed->gz, seed->color, seed->type);
        if (idx >= 0) {
            voxels[idx].debugClusterTag = seed->debugTag;
        }
    }
}

static void remove_buffered_static_voxels(const UnitVoxelBuffer *buffer)
{
    if (!buffer) {
        return;
    }
    int indices[VOXEL_ACTIVATION_UNIT_BUDGET];
    int idx_count = 0;
    for (int i = 0; i < buffer->count; ++i) {
        const UnitVoxelSeed *seed = &buffer->voxels[i];
        if (seed->voxelIndex >= 0 && idx_count < VOXEL_ACTIVATION_UNIT_BUDGET) {
            indices[idx_count++] = seed->voxelIndex;
        }
    }
    for (int i = 0; i < idx_count - 1; ++i) {
        for (int j = i + 1; j < idx_count; ++j) {
            if (indices[i] < indices[j]) {
                int tmp = indices[i];
                indices[i] = indices[j];
                indices[j] = tmp;
            }
        }
    }
    for (int i = 0; i < idx_count; ++i) {
        int idx = indices[i];
        if (idx >= 0 && idx < voxel_count) {
            Voxel *voxel = &voxels[idx];
            if (!voxel->simulate && voxel->owner == -1) {
                if (recycle_queue_push(voxel) && debugLogVoxelRecycle) {
                    TraceLog(LOG_INFO,
                             "[Recycle] enqueue-activation voxel=%d span=%d rest=(%d..%d,%d..%d,%d..%d)",
                             idx, voxel->span,
                             voxel->orig_min_gx, voxel->orig_max_gx,
                             voxel->orig_min_gy, voxel->orig_max_gy,
                             voxel->orig_min_gz, voxel->orig_max_gz);
                }
            }
            remove_voxel_index(idx);
        }
    }
}

static bool activation_try_enqueue(int voxel_idx,
                                   int activator,
                                   int center_gx, int center_gy, int center_gz,
                                   float radius_sq,
                                   UnitVoxelBuffer *buffer,
                                   int queue[],
                                   int *tail)
{
    if (!buffer || buffer->count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
        return false;
    }
    if (voxel_idx < 0 || voxel_idx >= voxel_count) {
        return false;
    }

    Voxel *candidate = &voxels[voxel_idx];
    if (candidate->simulate || candidate->span != 1 || candidate->pendingActivation) {
        return false;
    }

    if (radius_sq >= 0.0f) {
        float dx = (float)(candidate->gx - center_gx);
        float dy = (float)(candidate->gy - center_gy);
        float dz = (float)(candidate->gz - center_gz);
        if ((dx*dx + dy*dy + dz*dz) > radius_sq) {
            return false;
        }
    }

    candidate->pendingActivation = true;
    if (!unit_voxel_buffer_push(buffer,
                                candidate->gx, candidate->gy, candidate->gz,
                                candidate->color, candidate->type,
                                candidate->fixed, voxel_idx,
                                candidate->debugClusterTag, activator))
    {
        candidate->pendingActivation = false;
        return false;
    }

    if (queue && tail && *tail < VOXEL_ACTIVATION_UNIT_BUDGET) {
        queue[(*tail)++] = voxel_idx;
    }
    return true;
}

static int collect_static_activation_cluster(int seed_idx,
                                             int activator,
                                             int center_gx, int center_gy, int center_gz,
                                             float seed_radius_sq,
                                             UnitVoxelBuffer *buffer)
{
    if (!buffer || buffer->count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
        return 0;
    }

    int queue[VOXEL_ACTIVATION_UNIT_BUDGET];
    int head = 0;
    int tail = 0;
    int added = 0;

    if (!activation_try_enqueue(seed_idx, activator, center_gx, center_gy, center_gz,
                                seed_radius_sq, buffer, queue, &tail))
    {
        return 0;
    }
    added++;

    static const int neighbor_dirs[6][3] = {
        { 1,  0,  0 }, { -1,  0,  0 },
        { 0,  1,  0 }, {  0, -1,  0 },
        { 0,  0,  1 }, {  0,  0, -1 }
    };

    while (head < tail && buffer->count < VOXEL_ACTIVATION_UNIT_BUDGET) {
        int current_idx = queue[head++];
        if (current_idx < 0 || current_idx >= voxel_count) {
            continue;
        }
        const Voxel *current = &voxels[current_idx];
        int base_gx = current->gx;
        int base_gy = current->gy;
        int base_gz = current->gz;

        for (int n = 0; n < 6 && buffer->count < VOXEL_ACTIVATION_UNIT_BUDGET; ++n) {
            int nx = base_gx + neighbor_dirs[n][0];
            int ny = base_gy + neighbor_dirs[n][1];
            int nz = base_gz + neighbor_dirs[n][2];
            int neighbor_idx = table_get(nx, ny, nz);
            if (neighbor_idx < 0) {
                continue;
            }
            if (activation_try_enqueue(neighbor_idx, activator, center_gx, center_gy, center_gz,
                                       seed_radius_sq, buffer, queue, &tail))
            {
                added++;
                if (buffer->count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
                    break;
                }
            }
        }
    }

    return added;
}

static int expand_activation_cluster_unbounded(UnitVoxelBuffer *buffer, int startIndex,
                                               float dynamicBelief, int activator)
{
    if (!buffer) {
        return 0;
    }
    if (startIndex < 0) {
        startIndex = 0;
    }
    if (startIndex >= buffer->count) {
        return 0;
    }

    int queue[VOXEL_ACTIVATION_UNIT_BUDGET];
    int head = 0;
    int tail = 0;
    for (int i = startIndex; i < buffer->count && tail < VOXEL_ACTIVATION_UNIT_BUDGET; ++i) {
        int idx = buffer->voxels[i].voxelIndex;
        if (idx >= 0 && idx < voxel_count) {
            queue[tail++] = idx;
        }
    }

    static const int neighbor_dirs[6][3] = {
        { 1,  0,  0 }, { -1,  0,  0 },
        { 0,  1,  0 }, {  0, -1,  0 },
        { 0,  0,  1 }, {  0,  0, -1 }
    };

    int added = 0;
    while (head < tail && buffer->count < VOXEL_ACTIVATION_UNIT_BUDGET) {
        int current_idx = queue[head++];
        if (current_idx < 0 || current_idx >= voxel_count) {
            continue;
        }
        const Voxel *current = &voxels[current_idx];
        int base_gx = current->gx;
        int base_gy = current->gy;
        int base_gz = current->gz;

        for (int n = 0; n < 6 && buffer->count < VOXEL_ACTIVATION_UNIT_BUDGET; ++n) {
            int nx = base_gx + neighbor_dirs[n][0];
            int ny = base_gy + neighbor_dirs[n][1];
            int nz = base_gz + neighbor_dirs[n][2];
            int neighbor_idx = table_get(nx, ny, nz);
            if (neighbor_idx < 0) {
                continue;
            }
            if (neighbor_idx >= 0 && neighbor_idx < voxel_count) {
                Voxel *neighbor = &voxels[neighbor_idx];
                if (!dynamic_belief_overcomes_static(dynamicBelief, neighbor->freezeBelief)) {
                    continue;
                }
            }
            if (activation_try_enqueue(neighbor_idx, activator, 0, 0, 0,
                                       -1.0f, buffer, queue, &tail))
            {
                ++added;
                if (buffer->count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
                    break;
                }
            }
        }
    }

    return added;
}

static float compute_cluster_freeze_belief(const UnitVoxelBuffer *buffer, int startIndex)
{
    if (!buffer) {
        return 0.0f;
    }
    if (startIndex < 0) {
        startIndex = 0;
    }
    if (startIndex >= buffer->count) {
        return 0.0f;
    }
    float accum = 0.0f;
    float minBelief = 1.0f;
    int count = 0;
    for (int i = startIndex; i < buffer->count; ++i) {
        int idx = buffer->voxels[i].voxelIndex;
        if (idx < 0 || idx >= voxel_count) {
            continue;
        }
        float belief = voxels[idx].freezeBelief;
        accum += belief;
        if (belief < minBelief) {
            minBelief = belief;
        }
        ++count;
    }
    if (count <= 0) {
        return 0.0f;
    }
    float average = accum / (float)count;
    return 0.5f * (average + minBelief);
}

static bool dynamic_belief_overcomes_static(float dynamicBelief, float frozenBelief)
{
    float weightedDynamic = dynamicBelief * ACTIVATION_DYNAMIC_WEIGHT;
    float weightedFrozen = frozenBelief * FREEZE_BELIEF_IMPORTANCE + ACTIVATION_HYSTERESIS;
    return weightedDynamic >= weightedFrozen;
}

static bool activation_range_contains(int idx)
{
    if (activationLogStart < 0 || activationLogEnd <= activationLogStart) {
        return false;
    }
    return idx >= activationLogStart && idx < activationLogEnd;
}

static void log_activation_new_spans(int start_idx, int end_idx)
{
    int budget = DEBUG_ACTIVATION_LOG_INIT;
    for (int i = start_idx; i < end_idx && budget > 0; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate) {
            continue;
        }
        if (v->span <= 1) {
            continue;
        }
        int minx, maxx, miny, maxy, minz, maxz;
        voxel_grid_bounds(v, &minx, &maxx, &miny, &maxy, &minz, &maxz);
        float target_x = 0.5f * ((float)minx + (float)maxx + 1.0f) * VOXEL_SIZE;
        float target_y = 0.5f * ((float)miny + (float)maxy + 1.0f) * VOXEL_SIZE;
        float target_z = 0.5f * ((float)minz + (float)maxz + 1.0f) * VOXEL_SIZE;
        float dx = v->pos.x - target_x;
        float dy = v->pos.y - target_y;
        float dz = v->pos.z - target_z;
        TraceLog(LOG_INFO,
                 "[Activation] new span=%d idx=%d grid=(%d..%d,%d..%d,%d..%d) pos=(%.2f,%.2f,%.2f) "
                 "offset=(%.3f,%.3f,%.3f)",
                 v->span, i,
                 minx, maxx, miny, maxy, minz, maxz,
                 v->pos.x, v->pos.y, v->pos.z,
                 dx, dy, dz);
        --budget;
    }
}

static void log_activation_glue_mismatches(int start_idx, int end_idx)
{
    int budget = DEBUG_ACTIVATION_LOG_INIT;
    for (int i = start_idx; i < end_idx && budget > 0; ++i) {
        if (i < 0 || i >= voxel_count) {
            continue;
        }
        Voxel *v = &voxels[i];
        if (!v->simulate) {
            continue;
        }
        int neighbors[MAX_FACE_NEIGHBORS];
        int neighbor_count = gather_glued_neighbors(i, neighbors, MAX_FACE_NEIGHBORS);
        for (int n = 0; n < neighbor_count && budget > 0; ++n) {
            int nidx = neighbors[n];
            if (nidx < 0 || nidx >= voxel_count) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (!neighbor->simulate) {
                continue;
            }
            if (neighbor->span == v->span) {
                continue;
            }
            TraceLog(LOG_INFO,
                     "[Activation] glue span mismatch a=%d span=%d b=%d span=%d "
                     "posA=(%.2f,%.2f,%.2f) posB=(%.2f,%.2f,%.2f)",
                     i, v->span, nidx, neighbor->span,
                     v->pos.x, v->pos.y, v->pos.z,
                     neighbor->pos.x, neighbor->pos.y, neighbor->pos.z);
            --budget;
        }
    }
}

static bool activate_static_voxels_near_dynamic(void)
{
    UnitVoxelBuffer buffer = { 0 };
    unit_voxel_buffer_clear(&buffer);
    refresh_static_voxel_beliefs();
    update_dynamic_activation_beliefs();
    float radius = (float)VOXEL_ACTIVATION_RADIUS;
    float radius_sq = radius * radius;

    for (int i = 0; i < voxel_count; ++i) {
        if (buffer.count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
            break;
        }
        Voxel *dynamic = &voxels[i];
        if (!dynamic->simulate || dynamic->type == 1 || dynamic->type == 2 || dynamic->isBullet) {
            continue;
        }
        if (dynamic->activationCooldownFrames > 0) {
            continue;
        }
        int center_gx = (int)floorf(dynamic->pos.x / VOXEL_SIZE);
        int center_gy = (int)floorf(dynamic->pos.y / VOXEL_SIZE);
        int center_gz = (int)floorf(dynamic->pos.z / VOXEL_SIZE);

        int activator = dynamic->activator;
        if (activator < 0) {
            activator = dynamic->owner;
        }
        for (int gz = center_gz - VOXEL_ACTIVATION_RADIUS;
             gz <= center_gz + VOXEL_ACTIVATION_RADIUS && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET;
             ++gz)
        {
            for (int gy = center_gy - VOXEL_ACTIVATION_RADIUS;
                 gy <= center_gy + VOXEL_ACTIVATION_RADIUS && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET;
                 ++gy)
            {
                for (int gx = center_gx - VOXEL_ACTIVATION_RADIUS;
                     gx <= center_gx + VOXEL_ACTIVATION_RADIUS && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET;
                     ++gx)
                {
                    float dx = (float)(gx - center_gx);
                    float dy = (float)(gy - center_gy);
                    float dz = (float)(gz - center_gz);
                    if ((dx*dx + dy*dy + dz*dz) > radius_sq) {
                        continue;
                    }

                    int idx = table_get(gx, gy, gz);
                    if (idx < 0 || idx >= voxel_count) {
                        continue;
                    }
                    Voxel *candidate = &voxels[idx];
                    if (candidate->simulate || candidate->span != 1 ||
                        candidate->pendingActivation || candidate->activationCooldownFrames > 0) {
                        continue;
                    }

                    int previousCount = buffer.count;
                    int added = collect_static_activation_cluster(idx, activator,
                                                                  center_gx, center_gy, center_gz,
                                                                  radius_sq,
                                                                  &buffer);
                    if (added <= 0) {
                        rollback_activation_buffer(&buffer, previousCount);
                        continue;
                    }
                    float clusterBelief = compute_cluster_freeze_belief(&buffer, previousCount);
                    if (!dynamic_belief_overcomes_static(dynamic->activationBelief, clusterBelief)) {
                        rollback_activation_buffer(&buffer, previousCount);
                        continue;
                    }
                    expand_activation_cluster_unbounded(&buffer, previousCount,
                                                        dynamic->activationBelief, activator);
                    if (buffer.count >= VOXEL_ACTIVATION_UNIT_BUDGET) {
                        break;
                    }
                }
            }
        }
    }

    if (buffer.count <= 0) {
        return false;
    }

    remove_buffered_static_voxels(&buffer);
    int activation_base = voxel_count;
    int spawned = emit_multiscale_voxels_from_units(&buffer, false, true, -1);
    if (debugLogActivation) {
        TraceLog(LOG_INFO,
                 "[Activation] static->dynamic units=%d spawned=%d base=%d after=%d",
                 buffer.count, spawned, activation_base, voxel_count);
    }

    activationLogStart = activation_base;
    activationLogEnd = voxel_count;
    activationGlueLogBudget = DEBUG_ACTIVATION_LOG_INIT;
    rebuild_voxel_hash();
    rebuild_all_voxel_surfaces();
    rebuild_glue_constraints();
    if (debugLogActivation) {
        log_activation_new_spans(activation_base, voxel_count);
        log_activation_glue_mismatches(activation_base, voxel_count);
    }
    activationLogStart = -1;
    activationLogEnd = -1;
    refresh_static_voxel_beliefs();
    meshDirty = true;
    return true;
}

static bool restore_dynamic_voxel_to_static(int idx)
{
    if (idx < 0 || idx >= voxel_count) {
        return false;
    }
    Voxel voxel = voxels[idx];
    remove_voxel_index(idx);
    return spawn_static_covering_voxel(&voxel);
}

static bool restore_dynamic_snapshot(const Voxel *snapshot)
{
    if (!snapshot) {
        return false;
    }
    if (voxel_count >= MAX_VOXELS) {
        if (debugLogRestoreFailures) {
            TraceLog(LOG_WARNING,
                     "[RestoreFail] dynamic restore capacity exceeded span=%d pos=(%.2f,%.2f,%.2f)",
                     snapshot->span,
                     snapshot->pos.x, snapshot->pos.y, snapshot->pos.z);
        }
        return false;
    }
    voxels[voxel_count] = *snapshot;
    if (voxels[voxel_count].simulate) {
        voxels[voxel_count].lifeFrames = 0;
    }
    table_cache_invalidate();
    voxel_table_register(&voxels[voxel_count], voxel_count);
    ++voxel_count;
    return true;
}

static bool cell_contains_static_voxel(int x, int y, int z)
{
    int idx = table_get(x, y, z);
    if (idx < 0 || idx >= voxel_count) {
        return false;
    }
    return !voxels[idx].simulate;
}

static bool recycle_queue_push(const Voxel *voxel) {
    if (!voxel || recycleQueueCount >= MAX_VOXELS) {
        return false;
    }
    Voxel snapshot = *voxel;
    snapshot.lifeFrames = recycleFrameCounter;
    recycleQueue[recycleQueueTail] = snapshot;
    recycleQueueTail = (recycleQueueTail + 1) % MAX_VOXELS;
    ++recycleQueueCount;
    return true;
}

static bool recycle_queue_pop(Voxel *out) {
    if (!out || recycleQueueCount <= 0) {
        return false;
    }
    *out = recycleQueue[recycleQueueHead];
    recycleQueueHead = (recycleQueueHead + 1) % MAX_VOXELS;
    --recycleQueueCount;
    return true;
}

static bool voxel_is_at_rest_location(const Voxel *voxel) {
    if (!voxel) {
        return true;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    if (voxel->orig_min_gx <= voxel->orig_max_gx &&
        voxel->orig_min_gy <= voxel->orig_max_gy &&
        voxel->orig_min_gz <= voxel->orig_max_gz) {
        return (minx == voxel->orig_min_gx && maxx == voxel->orig_max_gx &&
                miny == voxel->orig_min_gy && maxy == voxel->orig_max_gy &&
                minz == voxel->orig_min_gz && maxz == voxel->orig_max_gz);
    }
    return (minx == voxel->rest_min_gx && maxx == voxel->rest_max_gx &&
            miny == voxel->rest_min_gy && maxy == voxel->rest_max_gy &&
            minz == voxel->rest_min_gz && maxz == voxel->rest_max_gz);
}

static bool voxel_outside_world_bounds(const Voxel *voxel) {
    if (!voxel) {
        return false;
    }
    VoxelWorldBounds bounds;
    voxel_world_bounds(voxel, &bounds);
    if (bounds.maxx < -FLOOR_SIZE || bounds.minx > FLOOR_SIZE) {
        return true;
    }
    if (bounds.maxz < -FLOOR_SIZE || bounds.minz > FLOOR_SIZE) {
        return true;
    }
    if (bounds.maxy < -VOXEL_SIZE || bounds.miny > (FLOOR_SIZE * 2.0f)) {
        return true;
    }
    return false;
}

static bool spawn_static_at_rest(const Voxel *snapshot) {
    if (!snapshot) {
        return false;
    }
    int minx = snapshot->rest_min_gx;
    int maxx = snapshot->rest_max_gx;
    int miny = snapshot->rest_min_gy;
    int maxy = snapshot->rest_max_gy;
    int minz = snapshot->rest_min_gz;
    int maxz = snapshot->rest_max_gz;
    if (snapshot->orig_min_gx <= snapshot->orig_max_gx &&
        snapshot->orig_min_gy <= snapshot->orig_max_gy &&
        snapshot->orig_min_gz <= snapshot->orig_max_gz) {
        minx = snapshot->orig_min_gx;
        maxx = snapshot->orig_max_gx;
        miny = snapshot->orig_min_gy;
        maxy = snapshot->orig_max_gy;
        minz = snapshot->orig_min_gz;
        maxz = snapshot->orig_max_gz;
    }
    if (minx > maxx || miny > maxy || minz > maxz) {
        return false;
    }
    int max_world_y = (int)ceilf((FLOOR_SIZE * 2.0f) / VOXEL_SIZE);
    remove_debris_static_in_columns(minx, maxx, minz, maxz, 0, max_world_y);
    remove_unowned_static_voxels_in_region(minx, maxx, miny, maxy, minz, maxz);
    if (!grid_region_is_free(minx, maxx, miny, maxy, minz, maxz)) {
        return false;
    }

    float center_x = 0.5f * ((float)minx + (float)maxx + 1.0f);
    float center_y = 0.5f * ((float)miny + (float)maxy + 1.0f);
    float center_z = 0.5f * ((float)minz + (float)maxz + 1.0f);
    float px = center_x * VOXEL_SIZE;
    float py = center_y * VOXEL_SIZE;
    float pz = center_z * VOXEL_SIZE;
    int span = (snapshot->span > 0) ? snapshot->span : 1;

    int idx = addVoxelSized(px, py, pz, true, false, snapshot->color, snapshot->type, span);
    if (idx < 0) {
        return false;
    }
    if (STATIC_REBUILD_ACTIVATION_COOLDOWN_FRAMES > 0) {
        voxels[idx].activationCooldownFrames = STATIC_REBUILD_ACTIVATION_COOLDOWN_FRAMES;
    }
    voxels[idx].owner = snapshot->owner;
    if (snapshot->orig_min_gx <= snapshot->orig_max_gx &&
        snapshot->orig_min_gy <= snapshot->orig_max_gy &&
        snapshot->orig_min_gz <= snapshot->orig_max_gz) {
        voxels[idx].orig_min_gx = snapshot->orig_min_gx;
        voxels[idx].orig_max_gx = snapshot->orig_max_gx;
        voxels[idx].orig_min_gy = snapshot->orig_min_gy;
        voxels[idx].orig_max_gy = snapshot->orig_max_gy;
        voxels[idx].orig_min_gz = snapshot->orig_min_gz;
        voxels[idx].orig_max_gz = snapshot->orig_max_gz;
    }
    return true;
}

static void recycle_dead_voxels(void) {
    recycleFrameCounter++;
    bool changed = false;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate) {
            continue;
        }
        voxel->lifeFrames++;
        if (voxel->lifeFrames > RECYCLE_DYNAMIC_MAX_FRAMES ||
            voxel_outside_world_bounds(voxel))
        {
            if (debugLogVoxelRecycle) {
                TraceLog(LOG_INFO,
                         "[Recycle] dynamic->static voxel=%d owner=%d life=%d",
                         i, voxel->owner, voxel->lifeFrames);
            }
            if (voxel->owner != -1) {
                remove_voxel_index(i);
                changed = true;
                --i;
                continue;
            }
            restore_dynamic_voxel_to_static(i);
            changed = true;
            --i;
        }
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->simulate) {
            continue;
        }
        if (voxel->owner >= 0) {
            voxel->lifeFrames++;
            if (voxel->lifeFrames > RECYCLE_OWNED_STATIC_MAX_FRAMES) {
                if (debugLogVoxelRecycle) {
                    TraceLog(LOG_INFO,
                             "[Recycle] remove-owned-static voxel=%d owner=%d life=%d",
                             i, voxel->owner, voxel->lifeFrames);
                }
                remove_voxel_index(i);
                changed = true;
                --i;
            }
            continue;
        }
        if (!voxel_is_at_rest_location(voxel)) {
            if (recycle_queue_push(voxel)) {
                if (debugLogVoxelRecycle) {
                    TraceLog(LOG_INFO,
                             "[Recycle] enqueue-static voxel=%d span=%d rest=(%d..%d,%d..%d,%d..%d)",
                             i, voxel->span,
                             voxel->rest_min_gx, voxel->rest_max_gx,
                             voxel->rest_min_gy, voxel->rest_max_gy,
                             voxel->rest_min_gz, voxel->rest_max_gz);
                }
                remove_voxel_index(i);
                changed = true;
                --i;
            }
        }
    }

    if ((recycleFrameCounter % RECYCLE_STATIC_RESTORE_INTERVAL) == 0 &&
        recycleQueueCount > 0 &&
        (recycleFrameCounter - recycleQueue[recycleQueueHead].lifeFrames) >= RECYCLE_STATIC_RESTORE_DELAY)
    {
        Voxel snapshot;
        if (recycle_queue_pop(&snapshot)) {
            if (!spawn_static_at_rest(&snapshot)) {
                if (debugLogVoxelRecycle) {
                    TraceLog(LOG_INFO,
                             "[Recycle] restore-failed span=%d rest=(%d..%d,%d..%d,%d..%d) queue=%d",
                             snapshot.span,
                             snapshot.rest_min_gx, snapshot.rest_max_gx,
                             snapshot.rest_min_gy, snapshot.rest_max_gy,
                             snapshot.rest_min_gz, snapshot.rest_max_gz,
                             recycleQueueCount);
                }
                recycle_queue_push(&snapshot);
            } else {
                if (debugLogVoxelRecycle) {
                    TraceLog(LOG_INFO,
                             "[Recycle] restore-ok span=%d rest=(%d..%d,%d..%d,%d..%d) queue=%d",
                             snapshot.span,
                             snapshot.rest_min_gx, snapshot.rest_max_gx,
                             snapshot.rest_min_gy, snapshot.rest_max_gy,
                             snapshot.rest_min_gz, snapshot.rest_max_gz,
                             recycleQueueCount);
                }
                if (remove_all_debris_static()) {
                    changed = true;
                }
                changed = true;
            }
        }
    }

    if (changed) {
        rebuild_voxel_hash();
        rebuild_all_voxel_surfaces();
        meshDirty = true;
    }
}

static bool line_contains_static_along_x(int x, int miny, int maxy, int minz, int maxz)
{
    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            if (cell_contains_static_voxel(x, y, z)) {
                return true;
            }
        }
    }
    return false;
}

static bool line_contains_static_along_y(int y, int minx, int maxx, int minz, int maxz)
{
    for (int x = minx; x <= maxx; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            if (cell_contains_static_voxel(x, y, z)) {
                return true;
            }
        }
    }
    return false;
}

static bool line_contains_static_along_z(int z, int minx, int maxx, int miny, int maxy)
{
    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            if (cell_contains_static_voxel(x, y, z)) {
                return true;
            }
        }
    }
    return false;
}

static uint8_t compute_static_support_mask(const Voxel *voxel)
{
    if (!voxel) {
        return 0;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    if (minx > maxx || miny > maxy || minz > maxz) {
        return 0;
    }

    uint8_t mask = 0;
    if (line_contains_static_along_x(maxx + 1, miny, maxy, minz, maxz)) mask |= (1u << 0);
    if (line_contains_static_along_x(minx - 1, miny, maxy, minz, maxz)) mask |= (1u << 1);
    if (line_contains_static_along_y(maxy + 1, minx, maxx, minz, maxz)) mask |= (1u << 2);
    if (line_contains_static_along_y(miny - 1, minx, maxx, minz, maxz)) mask |= (1u << 3);
    if (line_contains_static_along_z(maxz + 1, minx, maxx, miny, maxy)) mask |= (1u << 4);
    if (line_contains_static_along_z(minz - 1, minx, maxx, miny, maxy)) mask |= (1u << 5);

    return mask;
}

static bool append_static_neighbor_index(int neighbor_idx,
                                         const Voxel *self,
                                         int *out, int *count, int max_out)
{
    if (!out || !count || !self) {
        return false;
    }
    if (neighbor_idx < 0 || neighbor_idx >= voxel_count) {
        return false;
    }
    const Voxel *neighbor = &voxels[neighbor_idx];
    if (neighbor == self || neighbor->simulate) {
        return false;
    }
    if (list_contains_index(out, *count, neighbor_idx)) {
        return false;
    }
    if (*count >= max_out) {
        return false;
    }
    out[(*count)++] = neighbor_idx;
    return true;
}

static int gather_static_face_neighbors(const Voxel *voxel, int *out, int max_out)
{
    if (!voxel || !out || max_out <= 0 || voxel->simulate) {
        return 0;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    if (minx > maxx || miny > maxy || minz > maxz) {
        return 0;
    }
    int count = 0;
    for (int y = miny; y <= maxy && count < max_out; ++y) {
        for (int z = minz; z <= maxz && count < max_out; ++z) {
            int idx = table_get(maxx + 1, y, z);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    for (int y = miny; y <= maxy && count < max_out; ++y) {
        for (int z = minz; z <= maxz && count < max_out; ++z) {
            int idx = table_get(minx - 1, y, z);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    for (int x = minx; x <= maxx && count < max_out; ++x) {
        for (int z = minz; z <= maxz && count < max_out; ++z) {
            int idx = table_get(x, maxy + 1, z);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    for (int x = minx; x <= maxx && count < max_out; ++x) {
        for (int z = minz; z <= maxz && count < max_out; ++z) {
            int idx = table_get(x, miny - 1, z);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    for (int x = minx; x <= maxx && count < max_out; ++x) {
        for (int y = miny; y <= maxy && count < max_out; ++y) {
            int idx = table_get(x, y, maxz + 1);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    for (int x = minx; x <= maxx && count < max_out; ++x) {
        for (int y = miny; y <= maxy && count < max_out; ++y) {
            int idx = table_get(x, y, minz - 1);
            append_static_neighbor_index(idx, voxel, out, &count, max_out);
        }
    }
    return count;
}

static int bitcount_u8(uint8_t value)
{
    int c = 0;
    while (value) {
        c += (value & 1u);
        value >>= 1;
    }
    return c;
}

static void clear_static_belief_dirty(void)
{
    for (int i = 0; i < staticBeliefDirtyCount; ++i) {
        int idx = staticBeliefDirtyList[i];
        if (idx >= 0 && idx < MAX_VOXELS) {
            staticBeliefDirty[idx] = 0;
        }
    }
    staticBeliefDirtyCount = 0;
}

static void mark_static_belief_dirty_index(int idx)
{
    if (idx < 0 || idx >= voxel_count) {
        return;
    }
    Voxel *voxel = &voxels[idx];
    if (voxel->simulate) {
        return;
    }
    if (staticBeliefDirty[idx]) {
        return;
    }
    if (staticBeliefDirtyCount >= MAX_VOXELS) {
        staticBeliefsForceFullRefresh = true;
        return;
    }
    staticBeliefDirty[idx] = 1;
    staticBeliefDirtyList[staticBeliefDirtyCount++] = idx;
}

static void mark_static_beliefs_dirty_column_above(int gx, int gz, int gy)
{
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->simulate) {
            continue;
        }
        if (voxel->gx != gx || voxel->gz != gz) {
            continue;
        }
        if (voxel->gy >= gy) {
            mark_static_belief_dirty_index(i);
        }
    }
}

static void mark_static_beliefs_dirty_for_voxel(const Voxel *voxel)
{
    if (!voxel || voxel->simulate) {
        return;
    }
    int idx = (int)(voxel - voxels);
    mark_static_belief_dirty_index(idx);
    int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
    int neighborCount = gather_static_face_neighbors(voxel, neighbors, MAX_STATIC_COLLISION_NEIGHBORS);
    for (int n = 0; n < neighborCount; ++n) {
        mark_static_belief_dirty_index(neighbors[n]);
    }
    mark_static_beliefs_dirty_column_above(voxel->gx, voxel->gz, voxel->gy);
}

static void update_static_voxel_belief(int idx)
{
    if (idx < 0 || idx >= voxel_count) {
        return;
    }
    Voxel *voxel = &voxels[idx];
    if (voxel->simulate) {
        voxel->supportMask = 0;
        voxel->neighborSupport = 0;
        voxel->groundSupport = 0.0f;
        freezeBoundaryFlags[idx] = 0;
        voxel->freezeBelief = 0.0f;
        freezeBeliefScratch[idx] = 0.0f;
        return;
    }

    uint8_t supportMask = compute_static_support_mask(voxel);
    voxel->supportMask = supportMask;
    voxel->neighborSupport = (uint8_t)bitcount_u8(supportMask);

    VoxelWorldBounds bounds;
    voxel_world_bounds(voxel, &bounds);
    bool touchesGround = (bounds.miny <= GRID_EPSILON);
    //bool zeroBoundary = (!touchesGround) && (voxel->gy >=3.0f) && (voxel->surface[0] || voxel->surface[1] || voxel->surface[2] || voxel->surface[3] || voxel->surface[4] || voxel->surface[5]);
    bool zeroBoundary = (!touchesGround) && ( voxel->surface[2] || voxel->surface[3]);

    voxel->groundSupport = touchesGround ? 1.0f : 0.0f;

    uint8_t flags = 0;
    if (touchesGround) flags |= 1u;
    if (zeroBoundary)  flags |= 2u;
    freezeBoundaryFlags[idx] = flags;

    if (touchesGround) {
        voxel->freezeBelief = 1.0f;
    } else if (zeroBoundary) {
        voxel->freezeBelief = 0.0f;
    } else {
        voxel->freezeBelief = 0.5f;
    }
    freezeBeliefScratch[idx] = voxel->freezeBelief;
}

static void recompute_static_freeze_beliefs_path_length(void)
{
    int queue[MAX_VOXELS];
    int head = 0;
    int tail = 0;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->simulate) {
            voxel->supportMask = 0;
            voxel->neighborSupport = 0;
            voxel->groundSupport = 0.0f;
            voxel->freezeBelief = 0.0f;
            freezeBeliefScratch[i] = 0.0f;
            freezeBoundaryFlags[i] = 0;
            freezeDistance[i] = INT_MAX;
            overhangDistance[i] = INT_MAX;
            continue;
        }

        uint8_t supportMask = compute_static_support_mask(voxel);
        voxel->supportMask = supportMask;
        voxel->neighborSupport = (uint8_t)bitcount_u8(supportMask);

        VoxelWorldBounds bounds;
        voxel_world_bounds(voxel, &bounds);
        bool touchesGround = (bounds.miny <= GRID_EPSILON);
        voxel->groundSupport = touchesGround ? 1.0f : 0.0f;

        freezeBoundaryFlags[i] = touchesGround ? 1u : 0u;
        if (touchesGround) {
            freezeDistance[i] = 0;
            if (tail < MAX_VOXELS) {
                queue[tail++] = i;
            }
        } else {
            freezeDistance[i] = INT_MAX;
        }
        overhangDistance[i] = INT_MAX;
    }

    while (head < tail) {
        int idx = queue[head++];
        if (idx < 0 || idx >= voxel_count) {
            continue;
        }
        const Voxel *voxel = &voxels[idx];
        if (voxel->simulate) {
            continue;
        }
        int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
        int neighborCount = gather_static_face_neighbors(voxel, neighbors, MAX_STATIC_COLLISION_NEIGHBORS);
        for (int n = 0; n < neighborCount; ++n) {
            int nidx = neighbors[n];
            if (nidx < 0 || nidx >= voxel_count) {
                continue;
            }
            if (voxels[nidx].simulate) {
                continue;
            }
            int nextDist = freezeDistance[idx] + 1;
            if (nextDist < freezeDistance[nidx]) {
                freezeDistance[nidx] = nextDist;
                if (tail < MAX_VOXELS) {
                    queue[tail++] = nidx;
                }
            }
        }
    }

    head = 0;
    tail = 0;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].simulate) {
            continue;
        }
        bool touchesGround = (freezeBoundaryFlags[i] & 1u) != 0;
        bool hasBelow = (voxels[i].supportMask & (1u << 3)) != 0;
        if (touchesGround || hasBelow) {
            continue;
        }
        overhangDistance[i] = 0;
        if (tail < MAX_VOXELS) {
            queue[tail++] = i;
        }
    }

    while (head < tail) {
        int idx = queue[head++];
        if (idx < 0 || idx >= voxel_count) {
            continue;
        }
        const Voxel *voxel = &voxels[idx];
        if (voxel->simulate) {
            continue;
        }
        int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
        int neighborCount = gather_static_face_neighbors(voxel, neighbors, MAX_STATIC_COLLISION_NEIGHBORS);
        for (int n = 0; n < neighborCount; ++n) {
            int nidx = neighbors[n];
            if (nidx < 0 || nidx >= voxel_count) {
                continue;
            }
            if (voxels[nidx].simulate) {
                continue;
            }
            bool touchesGround = (freezeBoundaryFlags[nidx] & 1u) != 0;
            bool hasBelow = (voxels[nidx].supportMask & (1u << 3)) != 0;
            if (touchesGround || hasBelow) {
                continue;
            }
            int nextDist = overhangDistance[idx] + 1;
            if (nextDist < overhangDistance[nidx]) {
                overhangDistance[nidx] = nextDist;
                if (tail < MAX_VOXELS) {
                    queue[tail++] = nidx;
                }
            }
        }
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->simulate) {
            continue;
        }
        if (freezeDistance[i] == INT_MAX) {
            voxel->freezeBelief = 0.0f;
        } else {
            voxel->freezeBelief = powf(FREEZE_PATH_DECAY, (float)freezeDistance[i]);
        }
        if (overhangDistance[i] != INT_MAX) {
            voxel->freezeBelief *= powf(FREEZE_OVERHANG_DECAY,
                                        (float)(overhangDistance[i] + 1));
        }
        freezeBeliefScratch[i] = voxel->freezeBelief;
    }
}

static void full_refresh_static_voxel_beliefs(void)
{
    recompute_static_freeze_beliefs_path_length();
}

static void refresh_static_voxel_beliefs(void)
{
    if (!staticBeliefsInitialized || staticBeliefsForceFullRefresh) {
        full_refresh_static_voxel_beliefs();
        staticBeliefsInitialized = true;
        staticBeliefsForceFullRefresh = false;
        clear_static_belief_dirty();
        return;
    }

    if (staticBeliefDirtyCount <= 0) {
        return;
    }
    full_refresh_static_voxel_beliefs();
    staticBeliefsInitialized = true;
    staticBeliefsForceFullRefresh = false;
    clear_static_belief_dirty();
}

static bool voxel_connected_to_static_world(const Voxel *voxel)
{
    if (!voxel) {
        return false;
    }

    VoxelWorldBounds bounds;
    voxel_world_bounds(voxel, &bounds);
    if (bounds.miny <= GRID_EPSILON) {
        return true;
    }

    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    if (minx > maxx || miny > maxy || minz > maxz) {
        if (debugSupportLogBudget > 0) {
        TraceLog(LOG_TRACE,
                 "[Support] voxel bounds invalid pos=(%.2f,%.2f,%.2f) bounds=(%d..%d,%d..%d,%d..%d)",
                 voxel->pos.x, voxel->pos.y, voxel->pos.z,
                 minx, maxx, miny, maxy, minz, maxz);
            --debugSupportLogBudget;
        }
        return false;
    }

    bool adjacent_static = false;
    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            if (cell_contains_static_voxel(maxx + 1, y, z) ||
                cell_contains_static_voxel(minx - 1, y, z)) {
                adjacent_static = true;
                break;
            }
        }
        if (adjacent_static) {
            break;
        }
    }
    if (!adjacent_static) {
        for (int x = minx; x <= maxx; ++x) {
            for (int z = minz; z <= maxz; ++z) {
                if (cell_contains_static_voxel(x, maxy + 1, z) ||
                    cell_contains_static_voxel(x, miny - 1, z)) {
                    adjacent_static = true;
                    break;
                }
            }
            if (adjacent_static) {
                break;
            }
        }
    }
    if (!adjacent_static) {
        for (int x = minx; x <= maxx; ++x) {
            for (int y = miny; y <= maxy; ++y) {
                if (cell_contains_static_voxel(x, y, maxz + 1) ||
                    cell_contains_static_voxel(x, y, minz - 1)) {
                    adjacent_static = true;
                    break;
                }
            }
            if (adjacent_static) {
                break;
            }
        }
    }
    if (!adjacent_static && debugSupportLogBudget > 0) {
        TraceLog(LOG_TRACE,
                 "[Support] no-static-connection pos=(%.2f,%.2f,%.2f) bounds=(%d..%d,%d..%d,%d..%d)",
                 voxel->pos.x, voxel->pos.y, voxel->pos.z,
                 minx, maxx, miny, maxy, minz, maxz);
        --debugSupportLogBudget;
    }
    return adjacent_static;
}

static bool glue_cluster_has_static_support(const int *cluster, int cluster_count)
{
    if (!cluster || cluster_count <= 0) {
        return false;
    }
    int invalid_indices = 0;
    int checked = 0;
    for (int i = 0; i < cluster_count; ++i) {
        int idx = cluster[i];
        if (idx < 0 || idx >= voxel_count) {
            ++invalid_indices;
            continue;
        }
        ++checked;
        if (voxel_connected_to_static_world(&voxels[idx])) {
            return true;
        }
    }
    TraceLog(LOG_TRACE,
             "[Deactivate] cluster_support=false count=%d checked=%d invalid=%d reason=no-static-connection",
             cluster_count, checked, invalid_indices);
    return false;
}


static void keep_cluster_awake(const int *cluster, int cluster_count)
{
    if (!cluster || cluster_count <= 0) {
        return;
    }
    for (int i = 0; i < cluster_count; ++i) {
        int idx = cluster[i];
        if (idx < 0 || idx >= voxel_count) {
            continue;
        }
        Voxel *member = &voxels[idx];
        if (!member->simulate) {
            continue;
        }
        member->sleepFrames = (int)fminf(
            (float)member->sleepFrames,
            (float)(VOXEL_DEACTIVATION_FRAMES - 1));
    }
}

static bool deactivate_sleeping_voxels(void)
{
    bool changed = false;
    int remaining_budget = VOXEL_MAX_DEACTIVATIONS_PER_FRAME;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate) {
            voxel->sleepFrames = 0;
            continue;
        }
        if (voxel->type != 0 || voxel->isBullet) {
            voxel->sleepFrames = 0;
            continue;
        }

        float speed = v_length(voxel->vel);
        float strain = 0.0f;
        float shear = 0.0f;
        voxel_measure_strain(voxel, &strain, &shear);
        bool calm = (speed < VOXEL_DEACTIVATION_VELOCITY_THRESHOLD) &&
                    (strain < VOXEL_DEACTIVATION_STRAIN_THRESHOLD) &&
                    (shear < VOXEL_DEACTIVATION_SHEAR_THRESHOLD);
        if (calm) {
            if (voxel->sleepFrames < INT_MAX) {
                voxel->sleepFrames++;
            }
        } else {
            voxel->sleepFrames = 0;
        }
    }
    int idx = 0;
    while (idx < voxel_count && remaining_budget > 0) {
        Voxel *voxel = &voxels[idx];
        if (!voxel->simulate) {
            ++idx;
            continue;
        }
        if (voxel->sleepFrames < VOXEL_DEACTIVATION_FRAMES) {
            if (debugLogVoxelDeactivation && voxel->sleepFrames > 0) {
                TraceLog(LOG_INFO,
                         "[Deactivate] voxel=%d span=%d sleep=%d/%d waiting",
                         idx, voxel->span, voxel->sleepFrames, VOXEL_DEACTIVATION_FRAMES);
            }
            ++idx;
            continue;
        }

        int cluster_count = build_glue_cluster_indices(idx, glueClusterIndices);
        if (cluster_count <= 0) {
            if (debugLogVoxelDeactivation) {
                TraceLog(LOG_INFO,
                         "[Deactivate] voxel=%d span=%d sleep=%d cluster-empty",
                         idx, voxel->span, voxel->sleepFrames);
            }
            ++idx;
            continue;
        }

        bool cluster_ready = true;
        int cluster_fail_idx = -1;
        int cluster_fail_sleep = 0;
        const char *cluster_fail_reason = NULL;
        for (int c = 0; c < cluster_count; ++c) {
            int cidx = glueClusterIndices[c];
            if (cidx < 0 || cidx >= voxel_count) {
                cluster_ready = false;
                cluster_fail_idx = cidx;
                cluster_fail_reason = "invalid-index";
                break;
            }
            Voxel *member = &voxels[cidx];
            if (!member->simulate || member->sleepFrames < VOXEL_DEACTIVATION_FRAMES) {
                cluster_ready = false;
                cluster_fail_idx = cidx;
                cluster_fail_sleep = member->sleepFrames;
                cluster_fail_reason = member->simulate ? "insufficient-sleep" : "static-member";
                break;
            }
        }

        if (!cluster_ready) {
            if (debugLogVoxelDeactivation) {
                TraceLog(LOG_INFO,
                         "[Deactivate] voxel=%d cluster_ready=false fail_idx=%d sleep=%d reason=%s",
                         idx, cluster_fail_idx, cluster_fail_sleep,
                         cluster_fail_reason ? cluster_fail_reason : "unknown");
            }
            keep_cluster_awake(glueClusterIndices, cluster_count);
            ++idx;
            continue;
        }

        if (!glue_cluster_has_static_support(glueClusterIndices, cluster_count)) {
            if (debugLogVoxelDeactivation) {
                TraceLog(LOG_INFO,
                         "[Deactivate] voxel=%d cluster=%d lacks static support",
                         idx, cluster_count);
            }
            keep_cluster_awake(glueClusterIndices, cluster_count);
            ++idx;
            continue;
        }

        if (!restore_glue_cluster_to_static(glueClusterIndices, cluster_count)) {
            if (debugLogVoxelDeactivation) {
                TraceLog(LOG_WARNING,
                         "[Deactivate] voxel=%d cluster=%d restore failed",
                         idx, cluster_count);
            }
            ++idx;
            continue;
        }

        if (debugLogVoxelDeactivation) {
            TraceLog(LOG_INFO,
                     "[Deactivate] restored cluster starting=%d count=%d remaining_budget_before=%d",
                     idx, cluster_count, remaining_budget);
        }
        remaining_budget -= cluster_count;
        if (remaining_budget < 0) {
            remaining_budget = 0;
        }
        changed = true;
        continue;
    }

    if (changed) {
        rebuild_voxel_hash();
        rebuild_all_voxel_surfaces();
        rebuild_glue_constraints();
        refresh_static_voxel_beliefs();
        meshDirty = true;
    }
    return changed;
}

static unsigned char mix_color_channel(float a, float b, float t) {
    return (unsigned char)clampf(mixf(a, b, t), 0.0f, 255.0f);
}

static void build_oblique_voxel_pyramid(UnitVoxelBuffer *buffer) {
    if (!buffer) {
        return;
    }

    const int pyramid_height = 6;
    const int base_length = 6;
    const int base_width = 6;
    const int origin_x = -8;
    const int origin_y = 1;
    const int origin_z = -4;

    for (int level = 0; level < pyramid_height; ++level) {
        int layer_y = origin_y + level;
        int shrink_forward = level;
        int shrink_z = (level + 1) / 2;
        int min_x = origin_x;
        int max_x = origin_x + base_length - 1 - shrink_forward;
        int min_z = origin_z + shrink_z;
        int max_z = origin_z + base_width - 1 - shrink_z;
        if (max_x < min_x || max_z < min_z) {
            continue;
        }

        float t = (pyramid_height > 1)
            ? ((float)level / (float)(pyramid_height - 1))
            : 0.0f;
        Color col = {
            mix_color_channel(170.0f, 230.0f, t),
            mix_color_channel(100.0f, 160.0f, t),
            mix_color_channel(140.0f, 210.0f, t),
            255
        };

        for (int gx = min_x; gx <= max_x; ++gx) {
            for (int gz = min_z; gz <= max_z; ++gz) {
                if (!unit_voxel_buffer_push(buffer, gx, layer_y, gz, col, 0, true, -1, 0, -1)) {
                    TraceLog(LOG_WARNING, "[Pyramid] Unit voxel buffer full");
                    return;
                }
            }
        }
    }
}

static float grid_span_center(int min_g, int span) {
    return (min_g + 0.5f * (float)span) * VOXEL_SIZE;
}

static void apply_debug_tag_offset(int tag, int *x, int *y, int *z)
{
    if (tag <= 0 || tag >= DEBUG_CLUSTER_TAG_MAX) {
        return;
    }
    if (x) *x += debugTagOffset[tag][0];
    if (y) *y += debugTagOffset[tag][1];
    if (z) *z += debugTagOffset[tag][2];
}

static void init_debug_tag_offsets(void)
{
    memset(debugTagOffset, 0, sizeof(debugTagOffset));
    const int spacing = 16;
    const int cols = 5;
    const int max_tag = 19;
    for (int tag = 1; tag <= max_tag; ++tag) {
        int idx = tag - 1;
        int col = idx % cols;
        int row = idx / cols;
        debugTagOffset[tag][0] = col * spacing;
        debugTagOffset[tag][1] = 0;
        debugTagOffset[tag][2] = row * spacing;
    }
}

static int add_dynamic_span_voxel_at_grid_tag(int minx, int miny, int minz,
                                              int span, Color color, int debugTag)
{
    apply_debug_tag_offset(debugTag, &minx, &miny, &minz);
    float px = grid_span_center(minx, span);
    float py = grid_span_center(miny, span);
    float pz = grid_span_center(minz, span);
    int idx = addVoxelSized(px, py, pz, false, true, color, 0, span);
    if (idx >= 0) {
        voxels[idx].debugClusterTag = debugTag;
    }
    return idx;
}

static void add_dynamic_span_voxel_at_grid(int minx, int miny, int minz,
                                           int span, Color color)
{
    add_dynamic_span_voxel_at_grid_tag(minx, miny, minz, span, color, 0);
}

static void add_dynamic_unit_block_tag(int minx, int miny, int minz,
                                       int sx, int sy, int sz, Color color,
                                       int debugTag)
{
    apply_debug_tag_offset(debugTag, &minx, &miny, &minz);
    for (int x = 0; x < sx; ++x) {
        for (int y = 0; y < sy; ++y) {
            for (int z = 0; z < sz; ++z) {
                add_dynamic_span_voxel_at_grid_tag(minx + x, miny + y, minz + z,
                                                   1, color, debugTag);
            }
        }
    }
}

static void add_dynamic_unit_block(int minx, int miny, int minz,
                                   int sx, int sy, int sz, Color color)
{
    add_dynamic_unit_block_tag(minx, miny, minz, sx, sy, sz, color, 0);
}

static inline void addVoxelAt(int gx, int gy, int gz, Color c) {
    float px = (gx + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    float py = (gy + 0.5f) * VOXEL_SIZE;
    float pz = (gz + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    addVoxel(px, py, pz, true, false, c, 0);
}

static inline int grid_to_world_g(int g) {
    return g - (int)floorf(FLOOR_SIZE / VOXEL_SIZE);
}

static void add_static_box_at_grid(int minx, int maxx,
                                   int miny, int maxy,
                                   int minz, int maxz,
                                   Color color)
{
    for (int y = miny; y <= maxy; ++y) {
        for (int x = minx; x <= maxx; ++x) {
            for (int z = minz; z <= maxz; ++z) {
                add_static_voxel_at_grid(x, y, z, color, 0);
            }
        }
    }
}

static void add_static_disc_at_grid(int cx, int cz, int radius, int y, Color color)
{
    int r2 = radius * radius;
    for (int dx = -radius; dx <= radius; ++dx) {
        for (int dz = -radius; dz <= radius; ++dz) {
            if (dx * dx + dz * dz > r2) {
                continue;
            }
            add_static_voxel_at_grid(cx + dx, y, cz + dz, color, 0);
        }
    }
}

static void add_static_cylinder_at_grid(int cx, int cz, int radius,
                                        int miny, int maxy,
                                        Color color)
{
    for (int y = miny; y <= maxy; ++y) {
        add_static_disc_at_grid(cx, cz, radius, y, color);
    }
}

static void add_static_ring_at_grid(int cx, int cz,
                                    int radius_outer, int radius_inner,
                                    int miny, int maxy,
                                    Color color)
{
    int outer2 = radius_outer * radius_outer;
    int inner2 = radius_inner * radius_inner;
    for (int y = miny; y <= maxy; ++y) {
        for (int dx = -radius_outer; dx <= radius_outer; ++dx) {
            for (int dz = -radius_outer; dz <= radius_outer; ++dz) {
                int dist2 = dx * dx + dz * dz;
                if (dist2 > outer2 || dist2 < inner2) {
                    continue;
                }
                add_static_voxel_at_grid(cx + dx, y, cz + dz, color, 0);
            }
        }
    }
}

static void add_static_boulder_at_grid(int cx, int cz, int base_y,
                                       int radius, Color color)
{
    static const int offsets[5][2] = {
        { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }
    };
    int layers = radius;
    for (int y = 0; y <= layers; ++y) {
        int layer_radius = radius - (y / 2);
        if (layer_radius < 1) {
            layer_radius = 1;
        }
        int oi = y % 5;
        add_static_disc_at_grid(cx + offsets[oi][0],
                                cz + offsets[oi][1],
                                layer_radius, base_y + y, color);
    }
}

static void add_static_ring_with_door_at_grid(int cx, int cz,
                                              int radius_outer, int radius_inner,
                                              int miny, int maxy,
                                              Color color,
                                              int door_dir_x,
                                              int door_width, int door_height)
{
    int outer2 = radius_outer * radius_outer;
    int inner2 = radius_inner * radius_inner;
    for (int y = miny; y <= maxy; ++y) {
        bool door_row = (y <= door_height);
        for (int dx = -radius_outer; dx <= radius_outer; ++dx) {
            for (int dz = -radius_outer; dz <= radius_outer; ++dz) {
                int dist2 = dx * dx + dz * dz;
                if (dist2 > outer2 || dist2 < inner2) {
                    continue;
                }
                if (door_row && abs(dz) <= door_width) {
                    if ((door_dir_x > 0 && dx >= radius_inner) ||
                        (door_dir_x < 0 && dx <= -radius_inner)) {
                        continue;
                    }
                }
                add_static_voxel_at_grid(cx + dx, y, cz + dz, color, 0);
            }
        }
    }
}

static void add_static_hollow_cylinder_at_grid(int cx, int cz,
                                               int radius_outer, int radius_inner,
                                               int miny, int maxy,
                                               Color color)
{
    if (radius_inner < 0) {
        radius_inner = 0;
    }
    if (radius_inner >= radius_outer) {
        radius_inner = radius_outer - 1;
    }
    add_static_ring_at_grid(cx, cz, radius_outer, radius_inner, miny, maxy, color);
}

static void add_static_arch_at_grid(int cx, int cz, int base_y,
                                    int span, int height, int thickness,
                                    Color color)
{
    int radius = span / 2;
    if (radius < 2) {
        return;
    }
    float radius_f = (float)radius;
    int half = span / 2;
    for (int dx = -half; dx <= half; ++dx) {
        float xf = (float)dx;
        float y_arc = 0.0f;
        float inside = radius_f * radius_f - xf * xf;
        if (inside > 0.0f) {
            y_arc = radius_f - sqrtf(inside);
        }
        int y = base_y + height - (int)roundf(y_arc);
        for (int t = 0; t < thickness; ++t) {
            add_static_box_at_grid(cx + dx, cx + dx,
                                   y + t, y + t,
                                   cz - thickness, cz + thickness,
                                   color);
        }
    }
}

static void buildStackedPillar(int cx, int cz, int r,
                               int segH, int segCount, int gapH,
                               int faultY, Color c) {
    int yBase = 0;
    for (int s = 0; s < segCount; ++s) {
        for (int y = 0; y < segH; ++y) {
            for (int dx = -r; dx <= r; ++dx) {
                for (int dz = -r; dz <= r; ++dz) {
                    if (abs(dx) > r || abs(dz) > r) continue;
                    if (abs(dx) == r && abs(dz) == r) continue;
                    int worldY = yBase + y;
                    if (faultY >= 0 && worldY == faultY) {
                        if (abs(dx) <= r - 1 && abs(dz) <= r - 1) {
                            continue;
                        }
                    }
                    addVoxelAt(cx + dx, worldY, cz + dz, c);
                }
            }
        }
        yBase += segH + gapH;
    }
}

static void buildSplitPlatform(int centerX, int centerZ,
                               int baseY, int size, int seamHalf,
                               Color deckC, Color ribC) {
    int half = size / 2;
    int x0 = centerX - half;
    int z0 = centerZ - half;
    int x1 = x0 + size;
    int z1 = z0 + size;

    int midX = centerX;
    int midZ = centerZ;

    for (int x = x0; x < x1; ++x) {
        for (int z = z0; z < z1; ++z) {
            if (seamHalf > 0 &&
                (abs(x - midX) <= seamHalf || abs(z - midZ) <= seamHalf)) {
                continue;
            }
            addVoxelAt(x, baseY, z, deckC);
        }
    }

    int ribY = baseY - 1;
    if (ribY >= 0) {
        for (int x = x0; x < x1; ++x) {
            if (x == midX) continue;
            addVoxelAt(x, ribY, centerZ - 2, ribC);
            addVoxelAt(x, ribY, centerZ + 0, ribC);
            addVoxelAt(x, ribY, centerZ + 2, ribC);
        }
    }

    if (seamHalf > 0) {
        int capY = baseY + 2;
        for (int dz = -2; dz <= 2; ++dz) {
            if (dz == 0) {
                addVoxelAt(centerX, capY, centerZ + dz, deckC);
            }
        }
        for (int dx = -2; dx <= 2; ++dx) {
            if (dx == 0) {
                addVoxelAt(centerX + dx, capY, centerZ, deckC);
            }
        }
    }
}

static void buildLeg2x2Notched(int gx, int gz, int y0, int y1, int topY, Color c) {
    for (int y = y0; y <= y1; ++y) {
        for (int dx = 0; dx < 2; ++dx) {
            for (int dz = 0; dz < 2; ++dz) {
                if (y == topY && dx == 1 && dz == 1) continue;
                addVoxelAt(gx + dx, y, gz + dz, c);
            }
        }
    }
}

static void buildBox(int cx, int cy, int cz, int sx, int sy, int sz, Color c) {
    int x0 = cx - sx / 2, x1 = x0 + sx;
    int y0 = cy,       y1 = y0 + sy;
    int z0 = cz - sz / 2, z1 = z0 + sz;
    for (int x = x0; x < x1; ++x) {
        for (int y = y0; y < y1; ++y) {
            for (int z = z0; z < z1; ++z) {
                addVoxelAt(x, y, z, c);
            }
        }
    }
}

static void buildFortWalls(int cx, int cz, int w, int d, int h, Color c) {
    int x0 = cx - w/2, x1 = x0 + w;
    int z0 = cz - d/2, z1 = z0 + d;
    for (int y = 0; y < h; ++y) {
        for (int x = x0; x < x1; ++x) {
            addVoxelAt(x, y, z0, c);
            addVoxelAt(x, y, z1 - 1, c);
        }
        for (int z = z0; z < z1; ++z) {
            addVoxelAt(x0, y, z, c);
            addVoxelAt(x1 - 1, y, z, c);
        }
    }
}

static void carveFortGateOnX(int cx, int cz, int w, int d,
                             int gate_w, int gate_h, int side) {
    int x0 = cx - w/2, x1 = x0 + w;
    int z0 = cz - d/2, z1 = z0 + d;
    int midZ = cz;
    int xWall = (side < 0) ? x0 : (x1 - 1);
    int halfGate = gate_w / 2;
    for (int y = 0; y <= gate_h; ++y) {
        for (int z = midZ - halfGate; z <= midZ + halfGate; ++z) {
            int gx = grid_to_world_g(xWall);
            int gz = grid_to_world_g(z);
            int idx = table_get(gx, y, gz);
            if (idx >= 0 && idx < voxel_count) {
                remove_voxel_index(idx);
            }
        }
    }
}

static void carveFortWindows(int cx, int cz, int w, int d,
                             int y0, int y1, int window_w) {
    int x0 = cx - w/2, x1 = x0 + w;
    int z0 = cz - d/2, z1 = z0 + d;
    int midX = cx;
    int midZ = cz;
    int halfWin = window_w / 2;
    for (int y = y0; y <= y1; ++y) {
        for (int x = midX - halfWin; x <= midX + halfWin; ++x) {
            int idx = table_get(grid_to_world_g(x), y, grid_to_world_g(z0));
            if (idx >= 0 && idx < voxel_count) remove_voxel_index(idx);
            idx = table_get(grid_to_world_g(x), y, grid_to_world_g(z1 - 1));
            if (idx >= 0 && idx < voxel_count) remove_voxel_index(idx);
        }
        for (int z = midZ - halfWin; z <= midZ + halfWin; ++z) {
            int idx = table_get(grid_to_world_g(x0), y, grid_to_world_g(z));
            if (idx >= 0 && idx < voxel_count) remove_voxel_index(idx);
            idx = table_get(grid_to_world_g(x1 - 1), y, grid_to_world_g(z));
            if (idx >= 0 && idx < voxel_count) remove_voxel_index(idx);
        }
    }
}

static void buildFortRoof(int cx, int cz, int w, int d, int y, Color c) {
    int x0 = cx - w/2, x1 = x0 + w;
    int z0 = cz - d/2, z1 = z0 + d;
    for (int x = x0; x < x1; ++x) {
        for (int z = z0; z < z1; ++z) {
            addVoxelAt(x, y, z, c);
        }
    }
}

static void buildGateFrame(int cx, int cz, int w, int h, int depth, Color c) {
    int x0 = cx - w / 2;
    int x1 = x0 + w;
    int z0 = cz - depth / 2;
    int z1 = z0 + depth;
    for (int y = 0; y < h; ++y) {
        addVoxelAt(x0, y, z0, c);
        addVoxelAt(x1 - 1, y, z0, c);
        addVoxelAt(x0, y, z1 - 1, c);
        addVoxelAt(x1 - 1, y, z1 - 1, c);
    }
    for (int x = x0; x < x1; ++x) {
        for (int z = z0; z < z1; ++z) {
            addVoxelAt(x, h, z, c);
        }
    }
}

static void buildTestWorld(void) {
    // Floor
    int M = (int)(2.0f * FLOOR_SIZE / VOXEL_SIZE);

    // Pillars
    int pillar_height = 15; // 45 - 10
    int pillar_radius = 3;
    int pillar_positions[4][2] = {
        { M / 4, M / 4 },
        { M / 4, 3 * M / 4 },
        { 3 * M / 4, M / 4 },
        { 3 * M / 4, 3 * M / 4 }
    };

    for (int p = 0; p < 4; p++) {
        int cx = pillar_positions[p][0];
        int cz = pillar_positions[p][1];
        for (int y = 0; y <= pillar_height; y++) {
            for (int dx = -pillar_radius; dx <= pillar_radius; dx++) {
                for (int dz = -pillar_radius; dz <= pillar_radius; dz++) {
                    if (dx*dx + dz*dz > pillar_radius*pillar_radius) continue; // circular pillar
                    float px = (cx + dx + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                    float py = (y + 0.5f) * VOXEL_SIZE;
                    float pz = (cz + dz + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                    addVoxel(px, py, pz, true, false, (Color){ 200, 100, 50, 255 }, 0);
                }
            }
        }
    }

    // Central platform (n=1)
    int platform_size = 10;
    int platform_height = 1; // 15 / 3
    int platform_base_height = 5; // to keep top at same level (21)
    for (int y = platform_base_height; y <= platform_base_height + platform_height; y++) {
        for (int x = M/2 - platform_size/2; x <= M/2 + platform_size/2; x++) {
            for (int z = M/2 - platform_size/2; z <= M/2 + platform_size/2; z++) {
                float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                float py = (y + 0.5f) * VOXEL_SIZE;
                float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                addVoxel(px, py, pz, true, false, (Color){ 100, 200, 100, 255 }, 0);
                //addVoxelSized(px, py, pz, false, true, (Color){ 100, 200, 100, 255 }, 0, 1);
            }
        }
    }

    // Platform legs: 2x2 columns at each corner down to the floor.
    int platform_min = M/2 - platform_size/2;
    int platform_max = M/2 + platform_size/2;
    int leg_min_y = 0;
    int leg_max_y = platform_base_height;
    if (leg_max_y >= leg_min_y) {
        int leg_x[4] = { platform_min, platform_min, platform_max - 1, platform_max - 1 };
        int leg_z[4] = { platform_min, platform_max - 1, platform_min, platform_max - 1 };
        for (int corner = 0; corner < 4; ++corner) {
            for (int y = leg_min_y; y <= leg_max_y; ++y) {
                for (int dx = 0; dx < 2; ++dx) {
                    for (int dz = 0; dz < 2; ++dz) {
                        int x = leg_x[corner] + dx;
                        int z = leg_z[corner] + dz;
                        float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                        float py = (y + 0.5f) * VOXEL_SIZE;
                        float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                        addVoxel(px, py, pz, true, false, (Color){ 120, 160, 90, 255 }, 0);
                    }
                }
            }
        }
    }
}

static void buildBloodWorld(void) {
    int M = (int)(2.0f * FLOOR_SIZE / VOXEL_SIZE);
    int center = M / 2;
    int map_radius_cells = center - 2;

    int pillar_radius = 3;
    int pillar_seg_height = 6;
    int pillar_seg_count = 3;
    int pillar_gap = 0;

    int platform_size = (int)roundf(0.25f * (float)map_radius_cells);
    int platform_base_height = (int)roundf(0.10f * (float)map_radius_cells);
    if (platform_size < 12) platform_size = 12;
    if (platform_base_height < 5) platform_base_height = 5;

    int pillar_offset = (int)roundf(0.50f * (float)map_radius_cells);
    if (pillar_offset < 26) pillar_offset = 26;

    int base_offset = (int)roundf(0.62f * (float)map_radius_cells);
    if (base_offset < 30) base_offset = 30;

    Color pillar_color = (Color){ 230, 160, 70, 255 };
    Color deck_color = (Color){ 90, 220, 150, 255 };
    Color rib_color = (Color){ 70, 190, 230, 255 };
    Color leg_color = (Color){ 200, 210, 70, 255 };
    Color cover_color = (Color){ 170, 110, 220, 255 };
    Color wall_color = (Color){ 220, 120, 120, 255 };

    int pillar_fault = pillar_seg_height;
    buildStackedPillar(center - pillar_offset, center - pillar_offset,
                       pillar_radius, pillar_seg_height, pillar_seg_count, pillar_gap,
                       pillar_fault, pillar_color);
    buildStackedPillar(center - pillar_offset, center + pillar_offset,
                       pillar_radius, pillar_seg_height, pillar_seg_count, pillar_gap,
                       pillar_fault, pillar_color);
    buildStackedPillar(center + pillar_offset, center - pillar_offset,
                       pillar_radius, pillar_seg_height, pillar_seg_count, pillar_gap,
                       pillar_fault, pillar_color);
    buildStackedPillar(center + pillar_offset, center + pillar_offset,
                       pillar_radius, pillar_seg_height, pillar_seg_count, pillar_gap,
                       pillar_fault, pillar_color);

    buildSplitPlatform(center, center, platform_base_height,
                       platform_size, 0, deck_color, rib_color);

    int half = platform_size / 2;
    int leg_min_y = 0;
    int leg_max_y = platform_base_height - 1;
    int topY = leg_max_y;
    buildLeg2x2Notched(center - half, center - half, leg_min_y, leg_max_y, topY, leg_color);
    buildLeg2x2Notched(center - half, center + half - 1, leg_min_y, leg_max_y, topY, leg_color);
    buildLeg2x2Notched(center + half - 1, center - half, leg_min_y, leg_max_y, topY, leg_color);
    buildLeg2x2Notched(center + half - 1, center + half - 1, leg_min_y, leg_max_y, topY, leg_color);

    int base_height = 9;
    int base_w = 16;
    int base_d = 12;
    int gate_w = 5;
    int gate_h = 6;
    int window_w = 3;
    int window_y0 = 3;
    int window_y1 = 5;
    buildFortWalls(center - base_offset, center, base_w, base_d, base_height, wall_color);
    buildFortWalls(center + base_offset, center, base_w, base_d, base_height, wall_color);
    buildFortRoof(center - base_offset, center, base_w, base_d, base_height, wall_color);
    buildFortRoof(center + base_offset, center, base_w, base_d, base_height, wall_color);
    carveFortGateOnX(center - base_offset, center, base_w, base_d, gate_w, gate_h, 1);
    carveFortGateOnX(center + base_offset, center, base_w, base_d, gate_w, gate_h, -1);
    carveFortWindows(center - base_offset, center, base_w, base_d, window_y0, window_y1, window_w);
    carveFortWindows(center + base_offset, center, base_w, base_d, window_y0, window_y1, window_w);

    buildBox(center - pillar_offset + 4, 0, center + 14, 4, 4, 4, cover_color);
    buildBox(center + pillar_offset - 4, 0, center - 14, 4, 4, 4, cover_color);
    buildBox(center - pillar_offset + 4, 0, center - 14, 4, 4, 4, cover_color);
    buildBox(center + pillar_offset - 4, 0, center + 14, 4, 4, 4, cover_color);

    int lane_len = 34;
    int lane_height = 5;
    for (int y = 0; y < lane_height; ++y) {
        for (int x = center - lane_len; x <= center + lane_len; ++x) {
            addVoxelAt(x, y, center + 18, wall_color);
            addVoxelAt(x, y, center - 18, wall_color);
        }
    }

    buildGateFrame(center - pillar_offset - 8, center + 18, 5, 6, 3, wall_color);
    buildGateFrame(center + pillar_offset + 8, center + 18, 5, 6, 3, wall_color);
    buildGateFrame(center - pillar_offset - 8, center - 18, 5, 6, 3, wall_color);
    buildGateFrame(center + pillar_offset + 8, center - 18, 5, 6, 3, wall_color);

    buildBox(center - pillar_offset - 10, 0, center + 18, 4, 4, 4, cover_color);
    buildBox(center + pillar_offset + 10, 0, center + 18, 4, 4, 4, cover_color);
    buildBox(center - pillar_offset - 10, 0, center - 18, 4, 4, 4, cover_color);
    buildBox(center + pillar_offset + 10, 0, center - 18, 4, 4, 4, cover_color);
}

static void buildDebugWorld(void) {
    //init_debug_tag_offsets();
    // Floating dynamic test clusters (mixed span sizes, zero initial glue stress).
    // Debug tags: 1=corner chunk, 2=stacked pillar, 3=three-span bar, 4=plate+cap,
    // 5=step chain, 6=solid cube, 7=long beam, 8=brace frame, 9=staggered stack,
    // 10=unit cube control, 11=single span control, 12=thin slab, 13=flat cross, 14=vertical column,
    // 15=slab as units, 16=slab as span-4, 17=beam as units, 18=beam as span-4, 19=checker slab,
    // 20=vertical beam (Y), 21=depth beam (Z), 22=YZ wall, 23=XZ slab, 24=stacked span-2 (Y).
    // {
    //     int bx = -16, by = 10, bz = -12;
    //     int tag = 1;
    //     add_dynamic_unit_block_tag(bx, by, bz, 2, 2, 2, (Color){ 210, 120, 90, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by, bz, 2, (Color){ 240, 160, 100, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx, by + 2, bz, 2, (Color){ 200, 90, 140, 255 }, tag);
    // }
    // {
    //     int bx = 6, by = 12, bz = -10;
    //     int tag = 2;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 4, (Color){ 180, 150, 80, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 4, by, bz, 2, (Color){ 220, 180, 90, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 4, bz + 1, 2, (Color){ 250, 210, 130, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 6, bz + 1, 1, (Color){ 240, 200, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 6, bz + 2, 1, (Color){ 240, 200, 120, 255 }, tag);
    // }
    // {
    //     int bx = -8, by = 16, bz = 6;
    //     int tag = 3;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 2, (Color){ 90, 170, 230, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by, bz, 2, (Color){ 70, 140, 210, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 4, by, bz, 2, (Color){ 60, 120, 190, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 2, bz, 1, (Color){ 120, 200, 250, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 2, bz + 1, 1, (Color){ 120, 200, 250, 255 }, tag);
    // }
    // {
    //     int bx = 10, by = 8, bz = 8;
    //     int tag = 4;
    //     add_dynamic_unit_block_tag(bx, by, bz, 3, 1, 3, (Color){ 140, 200, 140, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 1, bz + 1, 2, (Color){ 90, 170, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 3, bz + 1, 1, (Color){ 60, 140, 90, 255 }, tag);
    // }
    // {
    //     int bx = -2, by = 20, bz = -2;
    //     int tag = 5;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 2, (Color){ 200, 120, 210, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by, bz, 2, (Color){ 170, 90, 180, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx, by + 2, bz, 1, (Color){ 210, 140, 230, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 2, bz, 1, (Color){ 210, 140, 230, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 2, bz, 1, (Color){ 210, 140, 230, 255 }, tag);
    // }
    // {
    //     int bx = -14, by = 24, bz = 6;
    //     int tag = 6;
    //     add_dynamic_unit_block_tag(bx, by, bz, 3, 3, 3, (Color){ 200, 170, 110, 255 }, tag);
    // }
    // {
    //     int bx = 2, by = 22, bz = 12;
    //     int tag = 7;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 1, (Color){ 160, 210, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by, bz, 1, (Color){ 160, 210, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by, bz, 1, (Color){ 160, 210, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 3, by, bz, 1, (Color){ 160, 210, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 4, by, bz, 1, (Color){ 160, 210, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 1, bz, 1, (Color){ 120, 170, 90, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by - 1, bz, 1, (Color){ 120, 170, 90, 255 }, tag);
    // }
    // {
    //     int bx = 12, by = 18, bz = -6;
    //     int tag = 8;
    //     add_dynamic_unit_block_tag(bx, by, bz, 3, 1, 3, (Color){ 120, 150, 220, 255 }, tag);
    //     add_dynamic_unit_block_tag(bx, by + 2, bz, 3, 1, 3, (Color){ 120, 150, 220, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx, by + 1, bz, 1, (Color){ 80, 110, 200, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 1, bz, 1, (Color){ 80, 110, 200, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx, by + 1, bz + 2, 1, (Color){ 80, 110, 200, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 1, bz + 2, 1, (Color){ 80, 110, 200, 255 }, tag);
    // }
    // {
    //     int bx = -6, by = 26, bz = -14;
    //     int tag = 9;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 2, (Color){ 210, 140, 120, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 2, by + 1, bz, 2, (Color){ 190, 120, 100, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 3, bz + 1, 1, (Color){ 230, 160, 130, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 3, by + 3, bz + 1, 1, (Color){ 230, 160, 130, 255 }, tag);
    // }
    // {
    //     int bx = -20, by = 30, bz = 12;
    //     int tag = 10;
    //     add_dynamic_unit_block_tag(bx, by, bz, 2, 2, 2, (Color){ 200, 200, 200, 255 }, tag);
    // }
    // {
    //     int bx = -10, by = 30, bz = 12;
    //     int tag = 11;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 2, (Color){ 180, 180, 200, 255 }, tag);
    // }
    // {
    //     int bx = 0, by = 30, bz = 12;
    //     int tag = 12;
    //     add_dynamic_unit_block_tag(bx, by, bz, 4, 1, 4, (Color){ 160, 180, 220, 255 }, tag);
    // }
    // {
    //     int bx = 8, by = 30, bz = 14;
    //     int tag = 13;
    //     add_dynamic_unit_block_tag(bx, by, bz, 3, 1, 3, (Color){ 150, 160, 230, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by + 1, bz + 1, 1, (Color){ 110, 130, 210, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 1, by - 1, bz + 1, 1, (Color){ 110, 130, 210, 255 }, tag);
    // }
    // {
    //     int bx = 18, by = 30, bz = 12;
    //     int tag = 14;
    //     add_dynamic_unit_block_tag(bx, by, bz, 1, 5, 1, (Color){ 200, 150, 130, 255 }, tag);
    // }
    // {
    //     int bx = -22, by = 34, bz = -2;
    //     int tag = 15;
    //     add_dynamic_unit_block_tag(bx, by, bz, 4, 1, 4, (Color){ 130, 170, 210, 255 }, tag);
    // }
    // {
    //     int bx = -14, by = 34, bz = -2;
    //     int tag = 16;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 4, (Color){ 110, 150, 190, 255 }, tag);
    // }
    // {
    //     int bx = -6, by = 34, bz = -2;
    //     int tag = 17;
    //     add_dynamic_unit_block_tag(bx, by, bz, 6, 1, 1, (Color){ 140, 200, 160, 255 }, tag);
    // }
    // {
    //     int bx = 2, by = 34, bz = -2;
    //     int tag = 18;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 4, (Color){ 120, 180, 140, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx + 4, by, bz, 2, (Color){ 120, 180, 140, 255 }, tag);
    // }
    // {
    //     int bx = 10, by = 34, bz = -4;
    //     int tag = 19;
    //     for (int x = 0; x < 4; ++x) {
    //         for (int z = 0; z < 4; ++z) {
    //             if (((x + z) & 1) == 0) {
    //                 add_dynamic_span_voxel_at_grid_tag(bx + x, by, bz + z, 1,
    //                                                    (Color){ 180, 140, 160, 255 }, tag);
    //             }
    //         }
    //     }
    // }
    // {
    //     int bx = -22, by = 34, bz = -12;
    //     int tag = 20;
    //     add_dynamic_unit_block_tag(bx, by, bz, 1, 6, 1, (Color){ 160, 190, 230, 255 }, tag);
    // }
    // {
    //     int bx = -14, by = 34, bz = -12;
    //     int tag = 21;
    //     add_dynamic_unit_block_tag(bx, by, bz, 1, 1, 6, (Color){ 180, 210, 170, 255 }, tag);
    // }
    // {
    //     int bx = -6, by = 34, bz = -12;
    //     int tag = 22;
    //     add_dynamic_unit_block_tag(bx, by, bz, 1, 4, 4, (Color){ 140, 170, 200, 255 }, tag);
    // }
    // {
    //     int bx = 2, by = 34, bz = -12;
    //     int tag = 23;
    //     add_dynamic_unit_block_tag(bx, by, bz, 4, 1, 4, (Color){ 200, 170, 140, 255 }, tag);
    // }
    // {
    //     int bx = 12, by = 34, bz = -12;
    //     int tag = 24;
    //     add_dynamic_span_voxel_at_grid_tag(bx, by, bz, 2, (Color){ 170, 140, 200, 255 }, tag);
    //     add_dynamic_span_voxel_at_grid_tag(bx, by + 2, bz, 2, (Color){ 170, 140, 200, 255 }, tag);
    // }
    // for (int x = 0; x <= M; x++) {
    //     for (int z = 0; z <= M; z++) {
    //         float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //         float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //         addVoxel(px, 0, pz, true, false, (Color){ 150, 150, 150, 255 }, 0);
    //     }
    // }
    // UnitVoxelBuffer pyramid_units;
    // unit_voxel_buffer_clear(&pyramid_units);
    // build_oblique_voxel_pyramid(&pyramid_units);
    // emit_static_voxels_from_units(&pyramid_units);
    //Span-2 dynamic voxel near origin for floor collision testing
    // {
    //     int span = 4;
    //     float px = -2.0f * VOXEL_SIZE;
    //     float pz = 2.0f * VOXEL_SIZE;
    //     float py = 2.0f+0.5f * (float)span * VOXEL_SIZE;
    //     addVoxelSized(px, py, pz, false, true, (Color){ 240, 160, 60, 255 }, 0, span);
    // }
    // {
    //     int span = 7;
    //     float px = -2.0f * VOXEL_SIZE;
    //     float pz = 2.0f * VOXEL_SIZE;
    //     float py = 2.0f+0.5f * (float)span * VOXEL_SIZE;
    //     addVoxelSized(px, py*2.0f, pz, false, true, (Color){ 240, 160, 60, 255 }, 0, span);
    // }
    // // Static 1x4x4 pad with a span-2 block hovering above for collision testing
    // {
    //     const int layer_size = 4;
    //     const int layer_origin_x = 4;
    //     const int layer_origin_z = -6;
    //     const int layer_y = 1;
    //     Color pad_color = (Color){ 80, 140, 210, 255 };
    //     for (int lx = 0; lx < layer_size; ++lx) {
    //         for (int lz = 0; lz < layer_size; ++lz) {
    //             int gx = layer_origin_x + lx;
    //             int gz = layer_origin_z + lz;
    //             add_static_voxel_at_grid(gx, layer_y, gz, pad_color, 0);
    //         }
    //     }
    //     float center_x = (layer_origin_x + 0.5f * (float)layer_size) * VOXEL_SIZE;
    //     float center_z = (layer_origin_z + 0.5f * (float)layer_size) * VOXEL_SIZE;
    //     float center_y = ((float)layer_y + 0.5f) * VOXEL_SIZE;
    //     int span = 2;
    //     float static_half = 0.5f * VOXEL_SIZE;
    //     float vertical_gap = 2.0f * VOXEL_SIZE;
    //     float dynamic_half = 0.5f * VOXEL_SIZE * (float)span;
    //     float py = center_y + static_half + vertical_gap + dynamic_half;
    //     addVoxelSized(center_x, py, center_z, false, true, (Color){ 230, 80, 120, 255 }, 0, span);
    // }
}

// Build static demo cube of voxels
static void buildDemo(void) {
    buildTestWorld();
    //buildBloodWorld();
    //buildDebugWorld();
    rebuild_glue_constraints();
}



static int first_voxel_hit(Ray ray, float t_max, int ignore_id);
static void UpdateKdRatio(int player_index);
static Vector3 pick_player_spawn(int player_index);
static int brush_extent_for_voxel(const Voxel *v);
static int player_max_health(const Player *p);
static int player_max_shield(const Player *p);
static void update_player_ammo(float dt);

// Reset game: players and voxels
static void ResetGame(void) {
    // init players
    for (int i = 0; i < MAX_PLAYERS; i++) {
        players[i].pos = pick_player_spawn(i);
        players[i].yaw = playerSpawnYaw[i];
        players[i].pitch = 0;
        players[i].yaw_vel = 0;
        players[i].pitch_vel = 0;
        players[i].vel = (Vector3){0,0,0};
        players[i].onGround = true;
        players[i].vType = 0;
        players[i].kills = 0;
        players[i].debrisKills = 0;
        players[i].deaths = 0;
        UpdateKdRatio(i);
        players[i].health = BASE_HEALTH;
        players[i].shield = BASE_SHIELD;
        players[i].last_damage_time = 0.0f;
        players[i].ammo = AMMO_MAX;
        players[i].ammo_recharge_timer = 0.0f;
        players[i].last_shot_time = -1000.0f;
    }
    // clear voxels
    voxel_count = 0;
    staticBeliefsInitialized = false;
    staticBeliefsForceFullRefresh = false;
    staticBeliefDirtyCount = 0;
    dynamicGlueClustersInitialized = false;
    memset(debugTagBreakLogged, 0, sizeof(debugTagBreakLogged));
    memset(staticBeliefDirty, 0, sizeof(staticBeliefDirty));
    memset(staticBeliefQueued, 0, sizeof(staticBeliefQueued));
    // clear hash
    memset(table, 0, sizeof(table));
    // build static blocks
    buildDemo();
    rebuild_all_voxel_surfaces();
    //rebuild_glue_constraints();
    meshDirty = true;

}

static void update_player_ammo(float dt) {
    for (int i = 0; i < activePlayers; ++i) {
        Player *p = &players[i];
        if (p->ammo >= AMMO_MAX) {
            p->ammo = AMMO_MAX;
            p->ammo_recharge_timer = 0.0f;
            continue;
        }
        p->ammo_recharge_timer += dt;
        while (p->ammo < AMMO_MAX && p->ammo_recharge_timer >= AMMO_RECHARGE_SECONDS) {
            p->ammo++;
            p->ammo_recharge_timer -= AMMO_RECHARGE_SECONDS;
        }
    }
}

static void UpdateKdRatio(int player_index) {
    Player *p = &players[player_index];
    p->kd_ratio = (float)(p->kills + 1) / (p->deaths + 1);
}

static void apply_damage_to_player(int player_index, int attacker_index, int damage,
                                   bool award_kill, bool award_debris)
{
    if (player_index < 0 || player_index >= activePlayers) {
        return;
    }
    Player *player = &players[player_index];
    player->last_damage_time = (float)GetTime();
    if (player->shield > 0) {
        player->shield -= damage;
        if (player->shield < 0) {
            player->health += player->shield;
            player->shield = 0;
        }
    } else {
        player->health -= damage;
    }

    if (player->health <= 0) {
        if (debugLogSmush && debugLogSmushDeaths) {
            TraceLog(LOG_WARNING,
                     "[Smush] death victim=%d attacker=%d awardKill=%d awardDebris=%d hp=%d shield=%d",
                     player_index, attacker_index,
                     award_kill ? 1 : 0, award_debris ? 1 : 0,
                     player->health, player->shield);
        }
        player->deaths++;
        UpdateKdRatio(player_index);
        if (attacker_index >= 0 && attacker_index < activePlayers && attacker_index != player_index) {
            if (award_kill) {
                players[attacker_index].kills++;
                UpdateKdRatio(attacker_index);
            }
            if (award_debris) {
                players[attacker_index].debrisKills++;
                if (debugLogSmush && debugLogSmushDeaths) {
                    TraceLog(LOG_INFO,
                             "[Smush] credit attacker=%d victim=%d debrisKills=%d",
                             attacker_index, player_index, players[attacker_index].debrisKills);
                }
            }
        }
        player->pos = pick_player_spawn(player_index);
        player->vel = (Vector3){ 0, 0, 0 };
        player->onGround = true;
        if (player_index >= 0 && player_index < MAX_PLAYERS) {
            player->yaw = playerSpawnYaw[player_index];
        }
        player->pitch = 0;
        player->health = player_max_health(player);
        player->shield = player_max_shield(player);
    }
}

static bool activate_static_neighbors_of_region(int minx, int maxx,
                                                int miny, int maxy,
                                                int minz, int maxz,
                                                int activator)
{
    UnitVoxelBuffer buffer = { 0 };
    unit_voxel_buffer_clear(&buffer);

    int ex_minx = minx - 1;
    int ex_maxx = maxx + 1;
    int ex_miny = miny - 1;
    int ex_maxy = maxy + 1;
    int ex_minz = minz - 1;
    int ex_maxz = maxz + 1;

    if (ex_miny < 0) {
        ex_miny = 0;
    }

    for (int z = ex_minz; z <= ex_maxz && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET; ++z) {
        for (int y = ex_miny; y <= ex_maxy && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET; ++y) {
            for (int x = ex_minx; x <= ex_maxx && buffer.count < VOXEL_ACTIVATION_UNIT_BUDGET; ++x) {
                if (x >= minx && x <= maxx &&
                    y >= miny && y <= maxy &&
                    z >= minz && z <= maxz) {
                    continue;
                }
                int idx = table_get(x, y, z);
                if (idx < 0 || idx >= voxel_count) {
                    continue;
                }
                Voxel *candidate = &voxels[idx];
                if (candidate->simulate || candidate->span != 1 || candidate->pendingActivation) {
                    continue;
                }
                candidate->pendingActivation = true;
                if (!unit_voxel_buffer_push(&buffer,
                                            candidate->gx, candidate->gy, candidate->gz,
                                            candidate->color, candidate->type,
                                            candidate->fixed, idx,
                                            candidate->debugClusterTag, activator))
                {
                    candidate->pendingActivation = false;
                    break;
                }
            }
        }
    }

    if (buffer.count <= 0) {
        return false;
    }

    remove_buffered_static_voxels(&buffer);
    if (debugLogSmush) {
        debugSmushLogBudget = 32;
    }
    emit_multiscale_voxels_from_units(&buffer, false, true, -1);
    if (debugLogSmush && debugSmushLogBudget > 0) {
        TraceLog(LOG_INFO,
                 "[Smush] activated units=%d activator=%d region=(%d..%d,%d..%d,%d..%d)",
                 buffer.count, activator, minx, maxx, miny, maxy, minz, maxz);
        --debugSmushLogBudget;
    }
    return true;
}

static void rebuild_voxel_hash(void) {
    table_cache_invalidate();
    memset(table, 0, sizeof(table));
    table_cache_invalidate();
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        v->gx = (int)floorf(v->pos.x / VOXEL_SIZE);
        v->gy = (int)floorf(v->pos.y / VOXEL_SIZE);
        v->gz = (int)floorf(v->pos.z / VOXEL_SIZE);
        voxel_table_register(v, i);
    }
    voxelHashFramesSinceRebuild = 0;
}

// Physics step for voxels
static void physics_step(float dt) {    // Rebuild spatial hash
    (void)dt;
    rebuild_voxel_hash();
    // // Simulate dynamic voxels
    // for (int i = 0; i < voxel_count; i++) {
    //     Voxel *v = &voxels[i];
    //     if (!v->simulate) continue;

    //     // Apply gravity
    //     v->vel.y -= GRAVITY * dt;

    //     // Continuous collision detection
    //     Vector3 displacement = v_mul(v->vel, dt);
    //     float distance = v_length(displacement);
    //     bool hit_voxel = false;

    //     if (distance > 0.0001f) { // only cast if moving
    //         Ray ray = { v->pos, v_norm(v->vel) };
    //         int hit_id = first_voxel_hit(ray, distance, i);

    //         if (hit_id >= 0) {
    //             hit_voxel = true;

    //             // Stop the bullet
    //             v->simulate = false;
    //             v->fixed = true;
    //             v->pos = (Vector3){-999.0f, -999.0f, -999.0f};
    //             deactivate_constraints_for_voxel(i);

    //             int brushExtent = (voxelBrushSpan < 1) ? 1 : voxelBrushSpan;

    //             if (v->type == 1) { // DESTRUCTION
    //                 Voxel *u = &voxels[hit_id];
    //                 int anchorX = u->gx;
    //                 int anchorY = u->gy;
    //                 int anchorZ = u->gz;

    //                 for (int dx = 0; dx < brushExtent; dx++) {
    //                     for (int dy = 0; dy < brushExtent; dy++) {
    //                         for (int dz = 0; dz < brushExtent; dz++) {
    //                             int victim_idx = table_get(anchorX + dx, anchorY + dy, anchorZ + dz);
    //                             if (victim_idx >= 0) {
    //                                 Voxel *victim = &voxels[victim_idx];
    //                                 table_remove(victim->gx, victim->gy, victim->gz);
    //                                 mark_surface_neighbors(victim->pos);
    //                                 victim->simulate = false;
    //                                 victim->fixed = true;
    //                                 victim->pos = (Vector3){-999.0f, -999.0f, -999.0f};
    //                                 deactivate_constraints_for_voxel(victim_idx);
    //                             }
    //                         }
    //                     }
    //                 }
    //             } else { // CONSTRUCTION
    //                 int anchorX = v->gx;
    //                 int anchorY = v->gy;
    //                 int anchorZ = v->gz;

    //                 for (int dx = 0; dx < brushExtent; dx++) {
    //                     for (int dy = 0; dy < brushExtent; dy++) {
    //                         for (int dz = 0; dz < brushExtent; dz++) {
    //                             int targetX = anchorX + dx;
    //                             int targetY = anchorY + dy;
    //                             int targetZ = anchorZ + dz;
    //                             if (!occupied(targetX, targetY, targetZ)) {
    //                                 float px = (targetX + 0.5f) * VOXEL_SIZE;
    //                                 float py = (targetY + 0.5f) * VOXEL_SIZE;
    //                                 float pz = (targetZ + 0.5f) * VOXEL_SIZE;
    //                                 int new_idx = addVoxel(px, py, pz, true, false, v->color, 0);
    //                                 if (new_idx >= 0) {
    //                                     mark_surface(new_idx);
    //                                     mark_surface_neighbors(voxels[new_idx].pos);
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //             //reset mesh
    //             meshDirty = true;
    //         }
    //     }

    //     if (!hit_voxel) {
    //         // Move
    //         v->pos = v_add(v->pos, displacement);
    //         for (int j = 0; j < 2; j++) {
    //             float dx = v->pos.x - players[j].pos.x;
    //             float dy = v->pos.y - players[j].pos.y;
    //             float dz = v->pos.z - players[j].pos.z;
    //             if (fabsf(dx) < PLAYER_SIZE && fabsf(dy) < PLAYER_SIZE && fabsf(dz) < PLAYER_SIZE) {
    //                 players[j].last_damage_time = (float)GetTime();
    //                 if (players[j].shield > 0) {
    //                     players[j].shield -= VOXEL_DAMAGE;
    //                     if (players[j].shield < 0) {
    //                         players[j].health += players[j].shield;
    //                         players[j].shield = 0;
    //                     }
    //                 } else {
    //                     players[j].health -= VOXEL_DAMAGE;
    //                 }

    //                 v->simulate = false;
    //                 v->fixed    = true;
    //                 v->pos      = (Vector3){ -999.0f, -999.0f, -999.0f };
    //                 deactivate_constraints_for_voxel(i);

    //                 if (players[j].health <= 0) {
    //                     if (v->owner >= 0 && v->owner != j) {
    //                         players[v->owner].kills++;
    //                         players[j].deaths++;
    //                         UpdateKdRatio(v->owner);
    //                         UpdateKdRatio(j);
    //                     }
    //                     players[j].pos     = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
    //                     players[j].vel     = (Vector3){0,0,0};
    //                     players[j].onGround= true;
    //                     players[j].yaw     = (j == 0 ? 0 : 180);
    //                     players[j].pitch   = 0;
    //                     players[j].health = BASE_HEALTH / players[j].kd_ratio;
    //                     players[j].shield = BASE_SHIELD;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }
}

static void update_projectiles(float dt)
{
    if (debugLogSmush) {
        debugSmushLogBudget = 32;
    }
    bool static_changed = false;
    int i = 0;
    while (i < voxel_count) {
        Voxel *v = &voxels[i];
        if (!v->simulate || (v->type == 0 && !v->isBullet)) {
            ++i;
            continue;
        }

        Vector3 gravity = { 0.0f, -GRAVITY, 0.0f };
        v->vel = v_mul(v->vel, VELOCITY_DAMPING);
        Vector3 displacement = v_add(v_mul(v->vel, dt), v_mul(gravity, dt * dt));
        v->vel = v_add(v->vel, v_mul(gravity, dt));
        float distance = v_length(displacement);
        if (distance <= 1e-6f) {
            ++i;
            continue;
        }

        Vector3 start = v->pos;
        Vector3 end = v_add(start, displacement);

        bool handled = false;
        for (int j = 0; j < activePlayers; ++j) {
            Player *pl = &players[j];
            Vector3 box_min = {
                pl->pos.x - PLAYER_SIZE * 0.5f,
                pl->pos.y - PLAYER_SIZE * 0.5f,
                pl->pos.z - PLAYER_SIZE * 0.5f
            };
            Vector3 box_max = {
                pl->pos.x + PLAYER_SIZE * 0.5f,
                pl->pos.y + PLAYER_SIZE * 0.5f,
                pl->pos.z + PLAYER_SIZE * 0.5f
            };
            if (segment_intersects_aabb(start, end, box_min, box_max)) {
                apply_damage_to_player(j, v->owner, VOXEL_DAMAGE, true, false);
                remove_voxel_index(i);
                handled = true;
                break;
            }
        }
        if (handled) {
            continue;
        }

        Ray ray = { start, v_norm(v->vel) };
        int hit_id = first_voxel_hit(ray, distance, i);
        if (hit_id >= 0 && hit_id < voxel_count) {
            if (!voxels[hit_id].simulate) {
                int brushExtent = brush_extent_for_voxel(v);
                Voxel *hit_voxel = &voxels[hit_id];
                int halfBrush = brushExtent / 2;
                int anchorX = hit_voxel->gx - halfBrush;
                int anchorY = hit_voxel->gy - halfBrush;
                int anchorZ = hit_voxel->gz - halfBrush;
                if (v->type == 2 && hit_voxel->gy <= 0) {
                    anchorY = hit_voxel->gy + 1;
                }
                int minx = anchorX;
                int maxx = anchorX + brushExtent - 1;
                int miny = anchorY;
                int maxy = anchorY + brushExtent - 1;
                int minz = anchorZ;
                int maxz = anchorZ + brushExtent - 1;

                if (v->type == 1) {
                    remove_static_voxels_in_region_recycle(minx, maxx, miny, maxy, minz, maxz);
                    static_changed = true;
                } else if (v->type == 2) {
                    int min_g = (int)ceilf((-FLOOR_SIZE / VOXEL_SIZE) - 0.5f);
                    int max_g = (int)floorf((FLOOR_SIZE / VOXEL_SIZE) - 0.5f);
                    for (int dx = 0; dx < brushExtent; ++dx) {
                        for (int dy = 0; dy < brushExtent; ++dy) {
                            int targetY = anchorY + dy;
                            if (targetY < 0) {
                                continue;
                            }
                            for (int dz = 0; dz < brushExtent; ++dz) {
                                int targetX = anchorX + dx;
                                int targetZ = anchorZ + dz;
                                if (targetX < min_g || targetX > max_g ||
                                    targetZ < min_g || targetZ > max_g) {
                                    continue;
                                }
                                if (!occupied(targetX, targetY, targetZ)) {
                                    int idx = add_static_voxel_at_grid(targetX, targetY, targetZ, v->color, 0);
                                    if (idx >= 0) {
                                        tag_owned_static_voxel(idx, v->owner);
                                    }
                                    static_changed = true;
                                }
                            }
                        }
                    }
                }

                if (v->type != 1 && v->type != 2) {
                    int activator = v->activator;
                    if (activator < 0) {
                        activator = v->owner;
                    }
                    if (activate_static_neighbors_of_region(minx, maxx, miny, maxy, minz, maxz, activator)) {
                        static_changed = true;
                    }
                }
            }

            remove_voxel_index(i);
            continue;
        }

        if (v->type == 2 && start.y > 0.0f && end.y <= 0.0f) {
            float t = start.y / (start.y - end.y);
            Vector3 hit_pos = v_add(start, v_mul(displacement, t));
            int brushExtent = brush_extent_for_voxel(v);
            int halfBrush = brushExtent / 2;
            int hit_gx = (int)floorf(hit_pos.x / VOXEL_SIZE);
            int hit_gz = (int)floorf(hit_pos.z / VOXEL_SIZE);
            int anchorX = hit_gx - halfBrush;
            int anchorY = 0;
            int anchorZ = hit_gz - halfBrush;
            int minx = anchorX;
            int maxx = anchorX + brushExtent - 1;
            int miny = anchorY;
            int maxy = anchorY + brushExtent - 1;
            int minz = anchorZ;
            int maxz = anchorZ + brushExtent - 1;
            int min_g = (int)ceilf((-FLOOR_SIZE / VOXEL_SIZE) - 0.5f);
            int max_g = (int)floorf((FLOOR_SIZE / VOXEL_SIZE) - 0.5f);

            for (int dx = 0; dx < brushExtent; ++dx) {
                for (int dy = 0; dy < brushExtent; ++dy) {
                    int targetY = anchorY + dy;
                    if (targetY < 0) {
                        continue;
                    }
                    for (int dz = 0; dz < brushExtent; ++dz) {
                        int targetX = anchorX + dx;
                        int targetZ = anchorZ + dz;
                        if (targetX < min_g || targetX > max_g ||
                            targetZ < min_g || targetZ > max_g) {
                            continue;
                        }
                        if (!occupied(targetX, targetY, targetZ)) {
                            int idx = add_static_voxel_at_grid(targetX, targetY, targetZ, v->color, 0);
                            if (idx >= 0) {
                                tag_owned_static_voxel(idx, v->owner);
                            }
                            static_changed = true;
                        }
                    }
                }
            }

            remove_voxel_index(i);
            continue;
        }

        v->pos = end;
        translate_voxel_particles(v, displacement);
        ++i;
    }

    if (static_changed) {
        rebuild_voxel_hash();
        rebuild_all_voxel_surfaces();
        //rebuild_glue_constraints();
        meshDirty = true;
    }
}

static int resolve_smush_activator(int voxel_idx)
{
    if (voxel_idx < 0 || voxel_idx >= voxel_count) {
        return -1;
    }
    Voxel *voxel = &voxels[voxel_idx];
    if (voxel->activator >= 0) {
        return voxel->activator;
    }
    if (voxel->owner >= 0) {
        return voxel->owner;
    }
    int cluster_count = build_glue_cluster_indices(voxel_idx, glueClusterIndices);
    for (int c = 0; c < cluster_count; ++c) {
        int idx = glueClusterIndices[c];
        if (idx < 0 || idx >= voxel_count) {
            continue;
        }
        Voxel *member = &voxels[idx];
        if (member->activator >= 0) {
            return member->activator;
        }
        if (member->owner >= 0) {
            return member->owner;
        }
    }
    return -1;
}

static bool dynamic_voxel_glued_to_static(int voxel_idx)
{
    if (voxel_idx < 0 || voxel_idx >= voxel_count) {
        return false;
    }
    Voxel *voxel = &voxels[voxel_idx];
    if (!voxel->simulate || !voxel->glueEligible) {
        return false;
    }
    int cluster_count = build_glue_cluster_indices(voxel_idx, glueClusterIndices);
    if (cluster_count <= 0) {
        return false;
    }
    return glue_cluster_has_static_support(glueClusterIndices, cluster_count);
}

static void handle_pbd_projectile_hits(void)
{
    if (debugLogSmush) {
        debugSmushLogBudget = 32;
    }
    int i = 0;
    while (i < voxel_count) {
        Voxel *v = &voxels[i];
        if (!v->simulate || v->type != 0 || v->isBullet) {
            ++i;
            continue;
        }
        bool glued_to_static = dynamic_voxel_glued_to_static(i);
        bool removed = false;
        for (int j = 0; j < activePlayers; ++j) {
            if (v->owner == j) {
                continue;
            }
            float dx = v->pos.x - players[j].pos.x;
            float dy = v->pos.y - players[j].pos.y;
            float dz = v->pos.z - players[j].pos.z;
            if (fabsf(dx) < PLAYER_SIZE && fabsf(dy) < PLAYER_SIZE && fabsf(dz) < PLAYER_SIZE) {
                if (glued_to_static) {
                    ++i;
                    removed = true;
                    break;
                }
                bool is_projectile = !v->glueEligible;
                int activator = is_projectile ? -1 : resolve_smush_activator(i);
                bool debris = (!is_projectile && activator >= 0);
                int attacker = debris ? activator : v->owner;
                bool award_kill = !debris;
                bool award_debris = debris;
                if (debugLogSmush && debugLogSmushHits && debugSmushLogBudget > 0) {
                    TraceLog(LOG_INFO,
                             "[Smush] hit-check victim=%d voxel=%d owner=%d activator=%d simulate=%d span=%d projectile=%d debris=%d",
                             j, i, v->owner, activator, v->simulate ? 1 : 0, v->span,
                             is_projectile ? 1 : 0, debris ? 1 : 0);
                    --debugSmushLogBudget;
                }
                if (debris && v->activator < 0) {
                    v->activator = activator;
                }
                if (debugLogSmush && debugLogSmushHits && debris && debugSmushLogBudget > 0) {
                    TraceLog(LOG_INFO,
                             "[Smush] hit victim=%d attacker=%d voxel=%d span=%d",
                             j, attacker, i, v->span);
                    --debugSmushLogBudget;
                }
                apply_damage_to_player(j, attacker, VOXEL_DAMAGE, award_kill, award_debris);
                remove_voxel_index(i);
                removed = true;
                break;
            }
        }
        if (!removed) {
            ++i;
        }
    }
}

// Predict positions for the next step (equivalent to the GPU PredictPositions kernel).
static void integrate_particles(float dt) {
    const Vector3 gravity = { 0.0f, -GRAVITY*1.0f, 0.0f };
    const float dt_sq = dt * dt;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->type != 0 || voxel->isBullet) {
            continue;
        }

        for (int j = 0; j < VOXEL_PARTICLE_COUNT; ++j) {
            Particle *p = voxel_particle_at(voxel, j);

            p->prev_pos = p->pos;
            p->predicted_pos = p->pos;

            if (p->inv_mass == 0.0f) {
                continue;
            }

            // Damped semi-implicit Euler: carry over velocity and integrate constant gravity.
            p->vel = v_mul(p->vel, VELOCITY_DAMPING);
            Vector3 step = v_mul(p->vel, dt);
            Vector3 accel = v_mul(gravity, dt_sq);
            p->predicted_pos = v_add(p->predicted_pos, v_add(step, accel));
        }
    }
}

static Vector3 vgs_project(Vector3 onto, Vector3 vec) {
    float denom = v_dot(onto, onto);
    if (denom < VGS_EPS) {
        return (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    float scale = v_dot(onto, vec) / denom;
    return v_mul(onto, scale);
}

static void reset_voxel_shape_to_rest(const Voxel *voxel, Vector3 *p, const float *w,
                                      Vector3 centroid) {
    float edge = voxel->rest_edge;
    if (edge <= 0.0f) {
        edge = VOXEL_SIZE;
    }
    if (!v_isfinite(centroid)) {
        centroid = voxel->pos;
    }
    float half = 0.5f * edge;
    for (int i = 0; i < 8; ++i) {
        if (w[i] == 0.0f) {
            continue;
        }
        p[i] = (Vector3){
            centroid.x + corner_signs[i][0] * half,
            centroid.y + corner_signs[i][1] * half,
            centroid.z + corner_signs[i][2] * half
        };
    }
}

// Voxel Gram-Schmidt shape matching (Algorithm 1 in the paper) keeps each cell near-rest.
static void solve_voxel_shape(Voxel *voxel) {
    bool has_dynamic = false;
    bool needs_reset = false;
    Vector3 p[8];
    float w[8];

    for (int i = 0; i < 8; ++i) {
        Particle *part = &voxel->particles[i];
        if (!v_isfinite(part->predicted_pos)) {
            part->predicted_pos = part->pos;
            needs_reset = true;
        }
        p[i] = part->predicted_pos;
        w[i] = part->inv_mass;
        if (w[i] > 0.0f) {
            has_dynamic = true;
        }
    }

    if (!has_dynamic) {
        return;
    }

    const float rest_volume = voxel->rest_volume;
    const float rest_edge = voxel->rest_edge;
    Vector3 centroid = { 0.0f, 0.0f, 0.0f };

    if (needs_reset) {
        for (int i = 0; i < 8; ++i) {
            centroid = v_add(centroid, p[i]);
        }
        centroid = v_mul(centroid, 1.0f / 8.0f);
        reset_voxel_shape_to_rest(voxel, p, w, centroid);
        voxel->pos = centroid;
        for (int i = 0; i < 8; ++i) {
            voxel->particles[i].predicted_pos = p[i];
        }
        return;
    }

    for (int iter = 0; iter < VGS_ITERS; ++iter) {
        centroid = (Vector3){ 0.0f, 0.0f, 0.0f };
        for (int i = 0; i < 8; ++i) {
            centroid = v_add(centroid, p[i]);
        }
        centroid = v_mul(centroid, 1.0f / 8.0f);

        // Compute principal axes (v0..v2) and damp them toward orthogonality via Gram-Schmidt.
        Vector3 v0 = v_add(v_add(v_sub(p[1], p[0]), v_sub(p[3], p[2])),
                           v_add(v_sub(p[5], p[4]), v_sub(p[7], p[6])));
        v0 = v_mul(v0, 0.25f);

        Vector3 v1 = v_add(v_add(v_sub(p[2], p[0]), v_sub(p[3], p[1])),
                           v_add(v_sub(p[6], p[4]), v_sub(p[7], p[5])));
        v1 = v_mul(v1, 0.25f);

        Vector3 v2 = v_add(v_add(v_sub(p[4], p[0]), v_sub(p[5], p[1])),
                           v_add(v_sub(p[6], p[2]), v_sub(p[7], p[3])));
        v2 = v_mul(v2, 0.25f);

        Vector3 u0 = v_sub(v0, v_mul(v_add(vgs_project(v1, v0), vgs_project(v2, v0)), VGS_ALPHA));
        Vector3 u1 = v_sub(v1, v_mul(v_add(vgs_project(v2, v1), vgs_project(v0, v1)), VGS_ALPHA));
        Vector3 u2 = v_sub(v2, v_mul(v_add(vgs_project(v0, v2), vgs_project(v1, v2)), VGS_ALPHA));

        float len0 = v_length(u0);
        float len1 = v_length(u1);
        float len2 = v_length(u2);
        if (!isfinite(len0) || !isfinite(len1) || !isfinite(len2) ||
            len0 < VGS_EPS || len1 < VGS_EPS || len2 < VGS_EPS) {
            needs_reset = true;
            break;
        }

        float lenp0 = 4.0f * v_length(v0);
        float lenp1 = 4.0f * v_length(v1);
        float lenp2 = 4.0f * v_length(v2);
        float r_v = 1.0f;
        float denom = lenp0 * lenp1 * lenp2;
        float rest_demom = rest_edge * rest_edge * rest_edge;
        if (fabsf(denom) < VGS_EPS || !isfinite(denom)) {
            needs_reset = true;
            break;
        }
        if (fabs(denom-rest_demom) > VGS_EPS) {
            float ratio = (rest_edge * rest_edge * rest_edge) / denom;
            float root = cbrtf(fabsf(ratio));
            r_v = (ratio < 0.0f) ? -root : root;
        }

        float target0 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp0 * r_v));
        float target1 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp1 * r_v));
        float target2 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp2 * r_v));

        float edge_eps = VGS_EARLY_OUT_EPS * rest_edge;
        float volume_eps = VGS_EARLY_OUT_EPS * rest_volume;
        float d0 = fabsf(len0 - target0);
        float d1 = fabsf(len1 - target1);
        float d2 = fabsf(len2 - target2);
        float raw_volume = v_dot(v_cross(u0, u1), u2);
        if (!isfinite(raw_volume) || fabsf(raw_volume) < VGS_EPS) {
            needs_reset = true;
            break;
        }
        if (d0 <= edge_eps && d1 <= edge_eps && d2 <= edge_eps &&
            fabsf(raw_volume - rest_volume) <= volume_eps) {
            break;
        }

        if (fabs(len0-target0) > VGS_EPS) u0 = v_mul(u0, target0 / len0);
        if (fabs(len1-target1) > VGS_EPS) u1 = v_mul(u1, target1 / len1);
        if (fabs(len2-target2) > VGS_EPS) u2 = v_mul(u2, target2 / len2);

        // Volume correction mirrors the GPU "ResizeVoxelBasis" stage.
        float volume = v_dot(v_cross(u0, u1), u2);
        if (!isfinite(volume) || fabsf(volume) < VGS_EPS) {
            needs_reset = true;
            break;
        }
        if (fabsf(volume) > VGS_EPS && fabsf(volume-rest_volume) > VGS_EPS) {
            float scale = rest_volume / volume;
            float root = cbrtf(fabsf(scale));
            if (scale < 0.0f) {
                root = -root;
            }
            u0 = v_mul(u0, root);
            u1 = v_mul(u1, root);
            u2 = v_mul(u2, root);
        }

        // Rebuild the voxel corners from the orthogonal frame and push dynamic particles only.
        u0 = v_mul(u0, 0.5f);
        u1 = v_mul(u1, 0.5f);
        u2 = v_mul(u2, 0.5f);
        Vector3 new_p[8];
        new_p[0] = v_sub(v_sub(v_sub(centroid, u0), u1), u2);
        new_p[1] = v_sub(v_sub(v_add(centroid, u0), u1), u2);
        new_p[2] = v_sub(v_add(v_sub(centroid, u0), u1), u2);
        new_p[3] = v_sub(v_add(v_add(centroid, u0), u1), u2);
        new_p[4] = v_add(v_sub(v_sub(centroid, u0), u1), u2);
        new_p[5] = v_add(v_sub(v_add(centroid, u0), u1), u2);
        new_p[6] = v_add(v_add(v_sub(centroid, u0), u1), u2);
        new_p[7] = v_add(v_add(v_add(centroid, u0), u1), u2);

        for (int i = 0; i < 8; ++i) {
            if (w[i] == 0.0f)
                continue;
            p[i] = new_p[i];
        }

    }

    if (needs_reset) {
        centroid = (Vector3){ 0.0f, 0.0f, 0.0f };
        for (int i = 0; i < 8; ++i) {
            centroid = v_add(centroid, p[i]);
        }
        centroid = v_mul(centroid, 1.0f / 8.0f);
        reset_voxel_shape_to_rest(voxel, p, w, centroid);
    }

    voxel->pos = centroid;
    for (int i = 0; i < 8; ++i) {
        voxel->particles[i].predicted_pos = p[i];
    }
}

static void solve_voxel_glue(bool allow_break) {
    bool glue_break = false;
    const float hinge_break_angle = GLUE_BREAK_HINGE_ANGLE_DEG * DEG2RAD;
    if (debugLogGlue) {
        debugGlueSolveLogBudget = DEBUG_GLUE_SOLVE_LOG_INIT;
        debugGlueBreakLogBudget = DEBUG_GLUE_BREAK_LOG_INIT;
    }
    for (int i = 0; i < glueConstraintCount; ++i) {
        GlueConstraint *gc = &glueConstraints[i];
        if (!gc->active) {
            continue;
        }

        Voxel *coarse = &voxels[gc->coarseVoxel];
        Voxel *fine   = &voxels[gc->fineVoxel];
        if (!coarse->simulate && !fine->simulate) {
            continue;
        }

        Particle *coarseParticles[4];
        float weights[4];
        for (int k = 0; k < 4; ++k) {
            coarseParticles[k] = &coarse->particles[gc->coarseCorner[k]];
            weights[k] = gc->w[k];
        }
        Particle *fineParticle = &fine->particles[gc->fineCorner];

        Vector3 coarsePred[4];
        for (int k = 0; k < 4; ++k) {
            coarsePred[k] = coarseParticles[k]->predicted_pos;
        }
        if (!v_isfinite(coarsePred[0]) || !v_isfinite(coarsePred[1]) ||
            !v_isfinite(coarsePred[2]) || !v_isfinite(coarsePred[3]) ||
            !v_isfinite(fineParticle->predicted_pos)) {
            continue;
        }
        Vector3 coarseU = v_sub(coarsePred[1], coarsePred[0]);
        Vector3 coarseV = v_sub(coarsePred[2], coarsePred[0]);
        Vector3 coarseNormal = v_cross(coarseU, coarseV);
        float coarseNormalLenSq = v_dot(coarseNormal, coarseNormal);
        float baryUU = v_dot(coarseU, coarseU);
        float baryVV = v_dot(coarseV, coarseV);
        float baryUV = v_dot(coarseU, coarseV);
        float baryDet = baryUU * baryVV - baryUV * baryUV;
        bool baryValid = (coarseNormalLenSq > 1e-12f) && (fabsf(baryDet) > 1e-12f);
        float solveRawU = 0.0f;
        float solveRawV = 0.0f;
        float solveU = 0.0f;
        float solveV = 0.0f;
        if (baryValid) {
            Vector3 normalDir = v_mul(coarseNormal, 1.0f / sqrtf(coarseNormalLenSq));
            baryValid = face_local_coords(coarsePred[0], coarseU, coarseV, normalDir,
                                          baryUU, baryVV, baryUV, 1.0f / baryDet,
                                          fineParticle->predicted_pos,
                                          &solveRawU, &solveRawV);
            if (baryValid) {
                solveU = clampf(solveRawU, 0.0f, 1.0f);
                solveV = clampf(solveRawV, 0.0f, 1.0f);
            }
        }

        Vector3 basisU = v_norm(coarseU);
        Vector3 basisV = v_sub(coarseV, v_mul(basisU, v_dot(coarseV, basisU)));
        float basisVLen = v_length(basisV);
        if (basisVLen < 1e-6f) {
            continue;
        }
        basisV = v_mul(basisV, 1.0f / basisVLen);
        Vector3 basisN = v_norm(v_cross(basisU, basisV));

        Vector3 anchor = { 0.0f, 0.0f, 0.0f };
        for (int k = 0; k < 4; ++k) {
            anchor = v_add(anchor, v_mul(coarseParticles[k]->predicted_pos, weights[k]));
        }
        Vector3 target = v_add(anchor,
                               v_add(v_mul(basisU, gc->restLocalU),
                                     v_add(v_mul(basisV, gc->restLocalV),
                                           v_mul(basisN, gc->restLocalN))));
        Vector3 C = v_sub(target, fineParticle->predicted_pos);
        float violation = v_length(C);
        float rest_edge_min = fminf(coarse->rest_edge, fine->rest_edge);
        float break_distance = GLUE_BREAK_STRAIN * rest_edge_min;
        bool hinge_break = false;
        if (allow_break && hinge_break_angle > 0.0f) {
            Vector3 n_coarse = { 0.0f, 0.0f, 0.0f };
            Vector3 n_fine = { 0.0f, 0.0f, 0.0f };
            if (face_normal_predicted(coarse, gc->coarseCorner, &n_coarse) &&
                face_normal_predicted(fine, gc->fineCornerFace, &n_fine)) {
                float cur_dot = clampf(v_dot(n_coarse, n_fine), -1.0f, 1.0f);
                float cur_angle = acosf(cur_dot);
                if (fabsf(cur_angle - gc->restNormalAngle) > hinge_break_angle) {
                    hinge_break = true;
                }
            }
        }
        if (allow_break && (violation > break_distance || hinge_break)) {
            if (debugLogGlue && debugGlueBreakLogBudget > 0) {
                TraceLog(LOG_DEBUG,
                         "[GlueBreak] pair=(%d,%d) spans=(%d,%d) violation=%.5f break=%.5f hinge=%d uvPred=(%.3f,%.3f) rawUV=(%.3f,%.3f) baryValid=%s weights=(%.3f,%.3f,%.3f,%.3f) coarsePos=(%.2f,%.2f,%.2f) finePos=(%.2f,%.2f,%.2f)",
                         gc->coarseVoxel, gc->fineVoxel,
                         voxel_span_for_glue(coarse), voxel_span_for_glue(fine),
                         violation, break_distance,
                         hinge_break ? 1 : 0,
                         solveU, solveV,
                         solveRawU, solveRawV,
                         baryValid ? "true" : "false",
                         weights[0], weights[1], weights[2], weights[3],
                         coarse->pos.x, coarse->pos.y, coarse->pos.z,
                         fine->pos.x, fine->pos.y, fine->pos.z);
                --debugGlueBreakLogBudget;
            }
            if (debugLogGlueClusters) {
                int tagA = coarse->debugClusterTag;
                int tagB = fine->debugClusterTag;
                if (debug_should_log_tag_break(tagA)) {
                    int neighborsA[MAX_FACE_NEIGHBORS];
                    int neighborCountA = gather_glued_neighbors(gc->coarseVoxel,
                                                               neighborsA,
                                                               MAX_FACE_NEIGHBORS);
                    Vector3 rel = v_sub(coarse->vel, fine->vel);
                    const int uvOutU = (solveRawU < 0.0f || solveRawU > 1.0f) ? 1 : 0;
                    const int uvOutV = (solveRawV < 0.0f || solveRawV > 1.0f) ? 1 : 0;
                    TraceLog(LOG_INFO,
                             "[GlueTagBreak] tag=%d label=%s voxel=%d span=%d neighborCount=%d "
                             "violation=%.5f break=%.5f relVel=%.4f glueConstraints=%d "
                             "coarseEdge=%.4f fineEdge=%.4f rawUV=(%.4f,%.4f) uv=(%.4f,%.4f) "
                             "uvOut=(%d,%d) baryValid=%d dir=%s normalLen=%.6f baryDet=%.6f",
                             tagA, debug_cluster_tag_label(tagA),
                             gc->coarseVoxel, voxel_span_for_glue(coarse),
                             neighborCountA, violation, break_distance,
                             v_length(rel), glueConstraintCount,
                             coarse->rest_edge, fine->rest_edge,
                             solveRawU, solveRawV, solveU, solveV,
                             uvOutU, uvOutV, baryValid ? 1 : 0,
                             glue_direction_label_from_delta(gc->dirX, gc->dirY, gc->dirZ),
                             sqrtf(coarseNormalLenSq), baryDet);
                    debugTagBreakLogged[tagA] = 1;
                }
                if (debug_should_log_tag_break(tagB)) {
                    int neighborsB[MAX_FACE_NEIGHBORS];
                    int neighborCountB = gather_glued_neighbors(gc->fineVoxel,
                                                               neighborsB,
                                                               MAX_FACE_NEIGHBORS);
                    Vector3 rel = v_sub(fine->vel, coarse->vel);
                    const int uvOutU = (solveRawU < 0.0f || solveRawU > 1.0f) ? 1 : 0;
                    const int uvOutV = (solveRawV < 0.0f || solveRawV > 1.0f) ? 1 : 0;
                    TraceLog(LOG_INFO,
                             "[GlueTagBreak] tag=%d label=%s voxel=%d span=%d neighborCount=%d "
                             "violation=%.5f break=%.5f relVel=%.4f glueConstraints=%d "
                             "coarseEdge=%.4f fineEdge=%.4f rawUV=(%.4f,%.4f) uv=(%.4f,%.4f) "
                             "uvOut=(%d,%d) baryValid=%d dir=%s normalLen=%.6f baryDet=%.6f",
                             tagB, debug_cluster_tag_label(tagB),
                             gc->fineVoxel, voxel_span_for_glue(fine),
                             neighborCountB, violation, break_distance,
                             v_length(rel), glueConstraintCount,
                             coarse->rest_edge, fine->rest_edge,
                             solveRawU, solveRawV, solveU, solveV,
                             uvOutU, uvOutV, baryValid ? 1 : 0,
                             glue_direction_label_from_delta(gc->dirX, gc->dirY, gc->dirZ),
                             sqrtf(coarseNormalLenSq), baryDet);
                    debugTagBreakLogged[tagB] = 1;
                }
            }
            if (coarse->simulate) {
                coarse->skipCollisionVelocityFrames = GLUE_BREAK_VELOCITY_SKIP_FRAMES;
            }
            if (fine->simulate) {
                fine->skipCollisionVelocityFrames = GLUE_BREAK_VELOCITY_SKIP_FRAMES;
            }
            gc->active = false;
            glue_adjacency_remove_ref_pair(gc->coarseVoxel, gc->fineVoxel, 1);
            glue_break = true;
            if (DEBRIS_ACTIVATION_COOLDOWN_FRAMES > 0) {
                int neighbors[MAX_FACE_NEIGHBORS];
                if (coarse->simulate) {
                    int count = gather_glued_neighbors(gc->coarseVoxel,
                                                       neighbors,
                                                       MAX_FACE_NEIGHBORS);
                    if (count == 0) {
                        coarse->activationCooldownFrames = DEBRIS_ACTIVATION_COOLDOWN_FRAMES;
                    }
                }
                if (fine->simulate) {
                    int count = gather_glued_neighbors(gc->fineVoxel,
                                                       neighbors,
                                                       MAX_FACE_NEIGHBORS);
                    if (count == 0) {
                        fine->activationCooldownFrames = DEBRIS_ACTIVATION_COOLDOWN_FRAMES;
                    }
                }
            }
            continue;
        }
        if (violation < GLUE_EPS) {
            continue;
        }

        float invMassSum = fineParticle->inv_mass;
        for (int k = 0; k < 4; ++k) {
            float inv_m = coarseParticles[k]->inv_mass;
            float wk = weights[k];
            invMassSum += wk * wk * inv_m;
        }
        if (invMassSum <= 0.0f) {
            continue;
        }

        Vector3 lambda = v_mul(C, -GLUE_RELAXATION * gc->strength / invMassSum);
        float lambda_mag = v_length(lambda);
        Vector3 coarseDeltaAccum = { 0.0f, 0.0f, 0.0f };
        float coarseDeltaMax = 0.0f;
        for (int k = 0; k < 4; ++k) {
            float inv_m = coarseParticles[k]->inv_mass;
            if (inv_m <= 0.0f) {
                continue;
            }
            float scale = weights[k] * inv_m;
            Vector3 deltaMove = v_mul(lambda, scale);
            coarseParticles[k]->predicted_pos = v_add(coarseParticles[k]->predicted_pos,
                                                      deltaMove);
            coarseDeltaAccum = v_add(coarseDeltaAccum, deltaMove);
            float deltaLen = v_length(deltaMove);
            if (deltaLen > coarseDeltaMax) {
                coarseDeltaMax = deltaLen;
            }
        }
        Vector3 fineDelta = { 0.0f, 0.0f, 0.0f };
        if (fineParticle->inv_mass > 0.0f) {
            fineDelta = v_mul(lambda, fineParticle->inv_mass);
            fineParticle->predicted_pos = v_sub(fineParticle->predicted_pos, fineDelta);
        }

        if (debugLogGlue && debugGlueSolveLogBudget > 0) {
            TraceLog(LOG_DEBUG,
                     "[GlueSolve] pair=(%d,%d) spans=(%d,%d) violation=%.5f break=%.5f invMass=%.5f lambda=%.5f coarseDeltaMax=%.5f fineDelta=%.5f baryValid=%s uvPred=(%.3f,%.3f) rawUV=(%.3f,%.3f) weights=(%.3f,%.3f,%.3f,%.3f) coarsePosY=%.2f finePosY=%.2f coarseVelY=%.2f fineVelY=%.2f",
                     gc->coarseVoxel, gc->fineVoxel,
                     voxel_span_for_glue(coarse), voxel_span_for_glue(fine),
                     violation, break_distance,
                     invMassSum, lambda_mag,
                     coarseDeltaMax, v_length(fineDelta),
                     baryValid ? "true" : "false",
                     solveU, solveV,
                     solveRawU, solveRawV,
                     weights[0], weights[1], weights[2], weights[3],
                     coarse->pos.y, fine->pos.y,
                     coarse->vel.y, fine->vel.y);
            --debugGlueSolveLogBudget;
        }
    }
    if (glue_break) {
        /* adjacency updates applied incrementally per constraint */
    }
}

static void log_dynamic_voxel_positions(void) {
    if (!debugLogDynamicVoxels) {
        return;
    }
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (!voxel->simulate) {
            continue;
        }
        TraceLog(LOG_INFO,
                 "[DynamicVoxel] voxel=%d pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f) span=%d",
                 i,
                 voxel->pos.x, voxel->pos.y, voxel->pos.z,
                 voxel->vel.x, voxel->vel.y, voxel->vel.z,
                 voxel->span);
    }
}

static int gather_face_neighbors(int voxel_idx, const GlueDirection *dir,
                                 int out[], int max_out)
{
    const Voxel *voxel = &voxels[voxel_idx];
    if (!voxel->glueEligible) {
        return 0;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);

    int count = 0;
    if (dir->dx != 0) {
        int x = (dir->dx > 0) ? (maxx + 1) : (minx - 1);
        for (int y = miny; y <= maxy; ++y) {
            for (int z = minz; z <= maxz; ++z) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx == voxel_idx) continue;
                if (!voxels[idx].glueEligible) continue;
                bool seen = false;
                for (int n = 0; n < count; ++n) {
                    if (out[n] == idx) { seen = true; break; }
                }
                if (!seen && count < max_out) {
                    out[count++] = idx;
                }
            }
        }
    } else if (dir->dy != 0) {
        int y = (dir->dy > 0) ? (maxy + 1) : (miny - 1);
        for (int x = minx; x <= maxx; ++x) {
            for (int z = minz; z <= maxz; ++z) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx == voxel_idx) continue;
                if (!voxels[idx].glueEligible) continue;
                bool seen = false;
                for (int n = 0; n < count; ++n) {
                    if (out[n] == idx) { seen = true; break; }
                }
                if (!seen && count < max_out) {
                    out[count++] = idx;
                }
            }
        }
    } else {
        int z = (dir->dz > 0) ? (maxz + 1) : (minz - 1);
        for (int x = minx; x <= maxx; ++x) {
            for (int y = miny; y <= maxy; ++y) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx == voxel_idx) continue;
                if (!voxels[idx].glueEligible) continue;
                bool seen = false;
                for (int n = 0; n < count; ++n) {
                    if (out[n] == idx) { seen = true; break; }
                }
                if (!seen && count < max_out) {
                    out[count++] = idx;
                }
            }
        }
    }
    return count;
}

static void add_bilinear_glue_constraints_for_pair(int negativeIdx, int positiveIdx,
                                                   const GlueDirection *dir)
{
    int coarseIdx, fineIdx;
    bool coarsePositive, finePositive;
    order_coarse_fine_pair(negativeIdx, positiveIdx,
                           &coarseIdx, &coarsePositive,
                           &fineIdx, &finePositive);

    Voxel *coarse = &voxels[coarseIdx];
    Voxel *fine   = &voxels[fineIdx];

    int coarseFace[4];
    int fineFace[4];
    get_face_corners_for_direction(dir, !coarsePositive, coarseFace);
    get_face_corners_for_direction(dir, !finePositive, fineFace);

    uint8_t coarseMask = 0;
    for (int k = 0; k < 4; ++k) {
        coarseMask |= (uint8_t)(1u << coarseFace[k]);
    }

    int normalAxis = (dir->dx != 0) ? 0 : ((dir->dy != 0) ? 1 : 2);
    int axisU, axisV;
    if (normalAxis == 0) {
        axisU = 1; axisV = 2;
    } else if (normalAxis == 1) {
        axisU = 0; axisV = 2;
    } else {
        axisU = 0; axisV = 1;
    }

    float coarseMinU = voxel_rest_axis_min(coarse, axisU);
    float coarseMaxU = voxel_rest_axis_max(coarse, axisU);
    float coarseMinV = voxel_rest_axis_min(coarse, axisV);
    float coarseMaxV = voxel_rest_axis_max(coarse, axisV);
    float coarseExtentU = coarseMaxU - coarseMinU;
    float coarseExtentV = coarseMaxV - coarseMinV;
    if (coarseExtentU <= 0.0f || coarseExtentV <= 0.0f) {
        return;
    }

    Vector3 coarseRest[4];
    for (int k = 0; k < 4; ++k) {
        coarseRest[k] = voxel_rest_corner_world(coarse, coarseFace[k]);
    }
    Vector3 restUVec = v_sub(coarseRest[1], coarseRest[0]);
    Vector3 restVVec = v_sub(coarseRest[2], coarseRest[0]);
    Vector3 restBasisU = v_norm(restUVec);
    Vector3 restBasisV = v_sub(restVVec, v_mul(restBasisU, v_dot(restVVec, restBasisU)));
    float restBasisVLen = v_length(restBasisV);
    if (restBasisVLen < 1e-6f) {
        return;
    }
    restBasisV = v_mul(restBasisV, 1.0f / restBasisVLen);
    Vector3 restBasisN = v_norm(v_cross(restBasisU, restBasisV));
    float restNormalAngle = PI;
    Vector3 restNormalCoarse = { 0.0f, 0.0f, 0.0f };
    Vector3 restNormalFine = { 0.0f, 0.0f, 0.0f };
    if (face_normal_rest(coarse, coarseFace, &restNormalCoarse) &&
        face_normal_rest(fine, fineFace, &restNormalFine)) {
        float rest_dot = clampf(v_dot(restNormalCoarse, restNormalFine), -1.0f, 1.0f);
        restNormalAngle = acosf(rest_dot);
    }

    for (int c = 0; c < 4; ++c) {
        if (glueConstraintCount >= (int)(sizeof(glueConstraints) / sizeof(glueConstraints[0]))) {
            return;
        }

        int fineCorner = fineFace[c];
        float fineCoordU = voxel_rest_corner_axis_coord(fine, axisU, fineCorner);
        float fineCoordV = voxel_rest_corner_axis_coord(fine, axisV, fineCorner);
        float rawU = (fineCoordU - coarseMinU) / coarseExtentU;
        float rawV = (fineCoordV - coarseMinV) / coarseExtentV;
        float u = clampf(rawU, 0.0f, 1.0f);
        float v = clampf(rawV, 0.0f, 1.0f);

        float w0 = (1.0f - u) * (1.0f - v);
        float w1 =        u  * (1.0f - v);
        float w2 = (1.0f - u) *        v;
        float w3 =        u  *        v;
        float w_sum = w0 + w1 + w2 + w3;
        if (w_sum <= 0.0f) {
            continue;
        }
        float inv_w_sum = 1.0f / w_sum;
        w0 *= inv_w_sum;
        w1 *= inv_w_sum;
        w2 *= inv_w_sum;
        w3 *= inv_w_sum;

        Vector3 fineRest = voxel_rest_corner_world(fine, fineCorner);
        Vector3 restAnchor = v_add(v_add(v_mul(coarseRest[0], w0), v_mul(coarseRest[1], w1)),
                                   v_add(v_mul(coarseRest[2], w2), v_mul(coarseRest[3], w3)));
        Vector3 restDelta = v_sub(fineRest, restAnchor);

        GlueConstraint *gc = &glueConstraints[glueConstraintCount++];
        gc->coarseVoxel = coarseIdx;
        gc->fineVoxel   = fineIdx;
        for (int k = 0; k < 4; ++k) {
            gc->coarseCorner[k] = coarseFace[k];
        }
        gc->w[0] = w0;
        gc->w[1] = w1;
        gc->w[2] = w2;
        gc->w[3] = w3;
        gc->fineCorner = fineCorner;
        gc->coarseMask = coarseMask;
        gc->fineMask = (uint8_t)(1u << fineCorner);
        gc->restLocalU = v_dot(restDelta, restBasisU);
        gc->restLocalV = v_dot(restDelta, restBasisV);
        gc->restLocalN = v_dot(restDelta, restBasisN);
        gc->restNormalAngle = restNormalAngle;
        gc->strength = 1.0f;
        gc->dirX = dir->dx;
        gc->dirY = dir->dy;
        gc->dirZ = dir->dz;
        for (int k = 0; k < 4; ++k) {
            gc->fineCornerFace[k] = fineFace[k];
        }
        gc->active = true;
        glue_adjacency_add_ref_pair(coarseIdx, fineIdx, 1);

        float edgeStrengthU = (u < 0.5f) ? (1.0f - u) : u;
        float edgeStrengthV = (v < 0.5f) ? (1.0f - v) : v;
        float centerStrength = 0.25f + 0.75f *
                               (1.0f - fabsf(u - 0.5f) * 2.0f) *
                               (1.0f - fabsf(v - 0.5f) * 2.0f);
        edgeStrengthU = clampf(edgeStrengthU, 0.0f, 1.0f) * GLUE_VIRTUAL_EDGE_STRENGTH;
        edgeStrengthV = clampf(edgeStrengthV, 0.0f, 1.0f) * GLUE_VIRTUAL_EDGE_STRENGTH;
        centerStrength = clampf(centerStrength, 0.0f, 1.0f) * GLUE_VIRTUAL_CENTER_STRENGTH;

        if (edgeStrengthU > 0.0f &&
            glueConstraintCount < (int)(sizeof(glueConstraints) / sizeof(glueConstraints[0])))
        {
            float e0 = 0.0f, e1 = 0.0f, e2 = 0.0f, e3 = 0.0f;
            if (u < 0.5f) {
                e0 = 0.5f;
                e2 = 0.5f;
            } else {
                e1 = 0.5f;
                e3 = 0.5f;
            }
            Vector3 edgeAnchor = v_add(v_add(v_mul(coarseRest[0], e0), v_mul(coarseRest[1], e1)),
                                       v_add(v_mul(coarseRest[2], e2), v_mul(coarseRest[3], e3)));
            Vector3 edgeDelta = v_sub(fineRest, edgeAnchor);
            GlueConstraint *edgeGc = &glueConstraints[glueConstraintCount++];
            edgeGc->coarseVoxel = coarseIdx;
            edgeGc->fineVoxel   = fineIdx;
            for (int k = 0; k < 4; ++k) {
                edgeGc->coarseCorner[k] = coarseFace[k];
            }
            edgeGc->w[0] = e0;
            edgeGc->w[1] = e1;
            edgeGc->w[2] = e2;
            edgeGc->w[3] = e3;
            edgeGc->fineCorner = fineCorner;
            edgeGc->coarseMask = coarseMask;
            edgeGc->fineMask = (uint8_t)(1u << fineCorner);
            edgeGc->restLocalU = v_dot(edgeDelta, restBasisU);
            edgeGc->restLocalV = v_dot(edgeDelta, restBasisV);
            edgeGc->restLocalN = v_dot(edgeDelta, restBasisN);
            edgeGc->restNormalAngle = restNormalAngle;
            edgeGc->strength = edgeStrengthU;
            edgeGc->dirX = dir->dx;
            edgeGc->dirY = dir->dy;
            edgeGc->dirZ = dir->dz;
            for (int k = 0; k < 4; ++k) {
                edgeGc->fineCornerFace[k] = fineFace[k];
            }
            edgeGc->active = true;
            glue_adjacency_add_ref_pair(coarseIdx, fineIdx, 1);
        }

        if (edgeStrengthV > 0.0f &&
            glueConstraintCount < (int)(sizeof(glueConstraints) / sizeof(glueConstraints[0])))
        {
            float e0 = 0.0f, e1 = 0.0f, e2 = 0.0f, e3 = 0.0f;
            if (v < 0.5f) {
                e0 = 0.5f;
                e1 = 0.5f;
            } else {
                e2 = 0.5f;
                e3 = 0.5f;
            }
            Vector3 edgeAnchor = v_add(v_add(v_mul(coarseRest[0], e0), v_mul(coarseRest[1], e1)),
                                       v_add(v_mul(coarseRest[2], e2), v_mul(coarseRest[3], e3)));
            Vector3 edgeDelta = v_sub(fineRest, edgeAnchor);
            GlueConstraint *edgeGc = &glueConstraints[glueConstraintCount++];
            edgeGc->coarseVoxel = coarseIdx;
            edgeGc->fineVoxel   = fineIdx;
            for (int k = 0; k < 4; ++k) {
                edgeGc->coarseCorner[k] = coarseFace[k];
            }
            edgeGc->w[0] = e0;
            edgeGc->w[1] = e1;
            edgeGc->w[2] = e2;
            edgeGc->w[3] = e3;
            edgeGc->fineCorner = fineCorner;
            edgeGc->coarseMask = coarseMask;
            edgeGc->fineMask = (uint8_t)(1u << fineCorner);
            edgeGc->restLocalU = v_dot(edgeDelta, restBasisU);
            edgeGc->restLocalV = v_dot(edgeDelta, restBasisV);
            edgeGc->restLocalN = v_dot(edgeDelta, restBasisN);
            edgeGc->restNormalAngle = restNormalAngle;
            edgeGc->strength = edgeStrengthV;
            edgeGc->dirX = dir->dx;
            edgeGc->dirY = dir->dy;
            edgeGc->dirZ = dir->dz;
            for (int k = 0; k < 4; ++k) {
                edgeGc->fineCornerFace[k] = fineFace[k];
            }
            edgeGc->active = true;
            glue_adjacency_add_ref_pair(coarseIdx, fineIdx, 1);
        }

        if (centerStrength > 0.0f &&
            glueConstraintCount < (int)(sizeof(glueConstraints) / sizeof(glueConstraints[0])))
        {
            float c0 = 0.25f, c1 = 0.25f, c2 = 0.25f, c3 = 0.25f;
            Vector3 centerAnchor = v_add(v_add(v_mul(coarseRest[0], c0), v_mul(coarseRest[1], c1)),
                                         v_add(v_mul(coarseRest[2], c2), v_mul(coarseRest[3], c3)));
            Vector3 centerDelta = v_sub(fineRest, centerAnchor);
            GlueConstraint *centerGc = &glueConstraints[glueConstraintCount++];
            centerGc->coarseVoxel = coarseIdx;
            centerGc->fineVoxel   = fineIdx;
            for (int k = 0; k < 4; ++k) {
                centerGc->coarseCorner[k] = coarseFace[k];
            }
            centerGc->w[0] = c0;
            centerGc->w[1] = c1;
            centerGc->w[2] = c2;
            centerGc->w[3] = c3;
            centerGc->fineCorner = fineCorner;
            centerGc->coarseMask = coarseMask;
            centerGc->fineMask = (uint8_t)(1u << fineCorner);
            centerGc->restLocalU = v_dot(centerDelta, restBasisU);
            centerGc->restLocalV = v_dot(centerDelta, restBasisV);
            centerGc->restLocalN = v_dot(centerDelta, restBasisN);
            centerGc->restNormalAngle = restNormalAngle;
            centerGc->strength = centerStrength;
            centerGc->dirX = dir->dx;
            centerGc->dirY = dir->dy;
            centerGc->dirZ = dir->dz;
            for (int k = 0; k < 4; ++k) {
                centerGc->fineCornerFace[k] = fineFace[k];
            }
            centerGc->active = true;
            glue_adjacency_add_ref_pair(coarseIdx, fineIdx, 1);
        }

        if (debugLogGlue && debugGlueBuildLogBudget > 0) {
            TraceLog(LOG_DEBUG,
                     "[GlueBuild] pair=(%d,%d) spans=(%d,%d) dir=(%d,%d,%d) coarsePosSide=%s finePosSide=%s fineCorner=%d rawUV=(%.3f,%.3f) uv=(%.3f,%.3f) weights=(%.3f,%.3f,%.3f,%.3f) coarseCenter=(%.2f,%.2f,%.2f) fineCenter=(%.2f,%.2f,%.2f)",
                     coarseIdx, fineIdx,
                     voxel_span_for_glue(coarse), voxel_span_for_glue(fine),
                     dir->dx, dir->dy, dir->dz,
                     coarsePositive ? "+face" : "-face",
                     finePositive ? "+face" : "-face",
                     fineCorner,
                     rawU, rawV,
                     u, v,
                     w0, w1, w2, w3,
                     coarse->pos.x, coarse->pos.y, coarse->pos.z,
                     fine->pos.x, fine->pos.y, fine->pos.z);
            --debugGlueBuildLogBudget;
        }
        if (debugLogActivation && activationGlueLogBudget > 0) {
            int spanA = voxel_span_for_glue(coarse);
            int spanB = voxel_span_for_glue(fine);
            if (spanA != spanB &&
                (activation_range_contains(coarseIdx) ||
                 activation_range_contains(fineIdx)))
            {
                TraceLog(LOG_INFO,
                         "[Activation] glue build mismatch pair=(%d,%d) spans=(%d,%d) dir=(%d,%d,%d) "
                         "rawUV=(%.3f,%.3f) uv=(%.3f,%.3f) weights=(%.3f,%.3f,%.3f,%.3f) "
                         "coarseCenter=(%.2f,%.2f,%.2f) fineCenter=(%.2f,%.2f,%.2f)",
                         coarseIdx, fineIdx,
                         spanA, spanB,
                         dir->dx, dir->dy, dir->dz,
                         rawU, rawV,
                         u, v,
                         w0, w1, w2, w3,
                         coarse->pos.x, coarse->pos.y, coarse->pos.z,
                         fine->pos.x, fine->pos.y, fine->pos.z);
                --activationGlueLogBudget;
            }
        }
    }
}

static void rebuild_glue_constraints(void) {
    glueConstraintCount = 0;
    glue_adjacency_clear_all();
    if (debugLogGlue) {
        debugGlueBuildLogBudget = DEBUG_GLUE_BUILD_LOG_INIT;
    }

    for (int i = 0; i < voxel_count; ++i) {
        if (!voxels[i].glueEligible) {
            continue;
        }
        for (int d = 0; d < 3; ++d) {
            const GlueDirection *dir = &glueDirections[d];
            int neighbors[MAX_FACE_NEIGHBORS] = {0};
            int neighborCount = gather_face_neighbors(i, dir, neighbors, MAX_FACE_NEIGHBORS);
            for (int n = 0; n < neighborCount; ++n) {
                int neighbor_idx = neighbors[n];
                add_bilinear_glue_constraints_for_pair(i, neighbor_idx, dir);
            }
        }
    }
    glue_dynamic_voxel_to_static_neighbors();
}

static void deactivate_glue_constraints_between(int a, int b) {
    int removed = 0;
    for (int g = 0; g < glueConstraintCount; ++g) {
        GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }
        if ((gc->coarseVoxel == a && gc->fineVoxel == b) ||
            (gc->coarseVoxel == b && gc->fineVoxel == a))
        {
            gc->active = false;
            ++removed;
        }
    }
    if (removed > 0) {
        glue_adjacency_remove_ref_pair(a, b, removed);
    }
}

static int gather_glued_neighbors(int voxel_idx, int *out, int max_out) {
    if (!out || max_out <= 0 || voxel_idx < 0 || voxel_idx >= voxel_count) {
        return 0;
    }
    rebuild_glue_adjacency_if_dirty();
    int count = gluedNeighborCounts[voxel_idx];
    if (count > max_out) {
        count = max_out;
    }
    for (int i = 0; i < count; ++i) {
        out[i] = gluedNeighborList[voxel_idx][i];
    }
    return count;
}

static int build_glue_cluster_indices(int start_idx, int *out_indices)
{
    if (start_idx < 0 || start_idx >= voxel_count || !out_indices) {
        return 0;
    }
    if (!voxels[start_idx].simulate) {
        return 0;
    }

    memset(glueClusterVisited, 0, sizeof(unsigned char) * (size_t)voxel_count);

    int head = 0;
    int tail = 0;
    out_indices[tail++] = start_idx;
    glueClusterVisited[start_idx] = 1;

    while (head < tail) {
        int idx = out_indices[head++];
        int neighbors[MAX_FACE_NEIGHBORS];
        int neighbor_count = gather_glued_neighbors(idx, neighbors, MAX_FACE_NEIGHBORS);
        for (int n = 0; n < neighbor_count; ++n) {
            int neighbor = neighbors[n];
            if (neighbor < 0 || neighbor >= voxel_count) {
                continue;
            }
             if (!voxels[neighbor].simulate) {
                 continue;
             }
            if (glueClusterVisited[neighbor]) {
                continue;
            }
            if (tail >= MAX_VOXELS) {
                break;
            }
            glueClusterVisited[neighbor] = 1;
            out_indices[tail++] = neighbor;
        }
        if (tail >= MAX_VOXELS) {
            break;
        }
    }

    return tail;
}

static void build_glue_cluster_ids(int *out_cluster_id)
{
    if (!out_cluster_id || voxel_count <= 0) {
        return;
    }

    static unsigned char cluster_visited[MAX_VOXELS];
    static int cluster_members[MAX_VOXELS];

    for (int i = 0; i < voxel_count; ++i) {
        out_cluster_id[i] = -1;
    }
    memset(cluster_visited, 0, sizeof(unsigned char) * (size_t)voxel_count);

    for (int i = 0; i < voxel_count; ++i) {
        if (!voxels[i].simulate || !voxels[i].glueEligible) {
            continue;
        }
        if (cluster_visited[i]) {
            continue;
        }

        int head = 0;
        int tail = 0;
        cluster_members[tail++] = i;
        cluster_visited[i] = 1;

        while (head < tail) {
            int idx = cluster_members[head++];
            int neighbors[MAX_FACE_NEIGHBORS];
            int neighbor_count = gather_glued_neighbors(idx, neighbors, MAX_FACE_NEIGHBORS);
            for (int n = 0; n < neighbor_count; ++n) {
                int neighbor = neighbors[n];
                if (neighbor < 0 || neighbor >= voxel_count) {
                    continue;
                }
                if (!voxels[neighbor].simulate || !voxels[neighbor].glueEligible) {
                    continue;
                }
                if (cluster_visited[neighbor]) {
                    continue;
                }
                if (tail >= MAX_VOXELS) {
                    break;
                }
                cluster_visited[neighbor] = 1;
                cluster_members[tail++] = neighbor;
            }
            if (tail >= MAX_VOXELS) {
                break;
            }
        }

        int rep = cluster_members[0];
        for (int m = 1; m < tail; ++m) {
            if (cluster_members[m] < rep) {
                rep = cluster_members[m];
            }
        }
        for (int m = 0; m < tail; ++m) {
            out_cluster_id[cluster_members[m]] = rep;
        }
    }
}

static void log_dynamic_glue_cluster_breaks(void)
{
    if (!debugLogGlueClusters || voxel_count <= 0) {
        return;
    }

    static int currentClusterId[MAX_VOXELS];
    static int currentClusterSize[MAX_VOXELS];
    static int clusterMembers[MAX_VOXELS];
    static int oldClusterPrimary[MAX_VOXELS];
    static int oldClusterSecondary[MAX_VOXELS];
    static int oldClusterNewCount[MAX_VOXELS];
    static int oldClusterTag[MAX_VOXELS];
    static unsigned char oldClusterLogged[MAX_VOXELS];

    for (int i = 0; i < voxel_count; ++i) {
        currentClusterId[i] = -1;
        currentClusterSize[i] = 0;
    }

    for (int i = 0; i < voxel_count; ++i) {
        if (!voxels[i].simulate || !voxels[i].glueEligible) {
            continue;
        }
        if (currentClusterId[i] != -1) {
            continue;
        }

        int head = 0;
        int tail = 0;
        clusterMembers[tail++] = i;
        currentClusterId[i] = -2;

        while (head < tail) {
            int idx = clusterMembers[head++];
            int neighbors[MAX_FACE_NEIGHBORS];
            int neighbor_count = gather_glued_neighbors(idx, neighbors, MAX_FACE_NEIGHBORS);
            for (int n = 0; n < neighbor_count; ++n) {
                int neighbor = neighbors[n];
                if (neighbor < 0 || neighbor >= voxel_count) {
                    continue;
                }
                if (!voxels[neighbor].simulate || !voxels[neighbor].glueEligible) {
                    continue;
                }
                if (currentClusterId[neighbor] != -1) {
                    continue;
                }
                currentClusterId[neighbor] = -2;
                clusterMembers[tail++] = neighbor;
                if (tail >= MAX_VOXELS) {
                    break;
                }
            }
            if (tail >= MAX_VOXELS) {
                break;
            }
        }

        int rep = clusterMembers[0];
        for (int m = 1; m < tail; ++m) {
            if (clusterMembers[m] < rep) {
                rep = clusterMembers[m];
            }
        }
        for (int m = 0; m < tail; ++m) {
            currentClusterId[clusterMembers[m]] = rep;
        }
        currentClusterSize[rep] = tail;
    }

    if (dynamicGlueClustersInitialized) {
        for (int i = 0; i < voxel_count; ++i) {
            oldClusterPrimary[i] = -1;
            oldClusterSecondary[i] = -1;
            oldClusterNewCount[i] = 0;
            oldClusterTag[i] = 0;
            oldClusterLogged[i] = 0;
        }

        for (int i = 0; i < voxel_count; ++i) {
            if (!voxels[i].simulate || !voxels[i].glueEligible) {
                continue;
            }
            Voxel *voxel = &voxels[i];
            if (!voxel->prevGlueClusterValid) {
                continue;
            }
            int oldId = voxel->prevGlueClusterId;
            int newId = currentClusterId[i];
            if (oldId < 0 || newId < 0) {
                continue;
            }

            if (oldClusterTag[oldId] == 0) {
                oldClusterTag[oldId] = voxel->prevGlueClusterTag;
            }
            if (oldClusterPrimary[oldId] == -1) {
                oldClusterPrimary[oldId] = newId;
                oldClusterNewCount[oldId] = 1;
            } else if (oldClusterPrimary[oldId] != newId) {
                if (oldClusterSecondary[oldId] == -1) {
                    oldClusterSecondary[oldId] = newId;
                    oldClusterNewCount[oldId] = 2;
                } else if (oldClusterSecondary[oldId] != newId) {
                    oldClusterNewCount[oldId] = 3;
                }
            }
        }

        for (int i = 0; i < voxel_count; ++i) {
            const Voxel *voxel = &voxels[i];
            if (!voxel->prevGlueClusterValid) {
                continue;
            }
            int oldId = voxel->prevGlueClusterId;
            if (oldId < 0 || oldId >= voxel_count) {
                continue;
            }
            if (oldClusterLogged[oldId]) {
                continue;
            }
            oldClusterLogged[oldId] = 1;
            int oldSize = voxel->prevGlueClusterSize;
            if (oldSize <= 1) {
                continue;
            }
            int newCount = oldClusterNewCount[oldId];
            if (newCount > 1) {
                int primary = oldClusterPrimary[oldId];
                int secondary = oldClusterSecondary[oldId];
                int primarySize = (primary >= 0) ? currentClusterSize[primary] : 0;
                int secondarySize = (secondary >= 0) ? currentClusterSize[secondary] : 0;
                int primaryTag = (primary >= 0 && primary < voxel_count) ? voxels[primary].debugClusterTag : 0;
                int secondaryTag = (secondary >= 0 && secondary < voxel_count) ? voxels[secondary].debugClusterTag : 0;
                TraceLog(LOG_INFO,
                         "[GlueClusterBreak] cluster=%d tag=%d size=%d new_clusters=%d primary=%d tag=%d size=%d secondary=%d tag=%d size=%d",
                         oldId, oldClusterTag[oldId], oldSize, newCount,
                         primary, primaryTag, primarySize,
                         secondary, secondaryTag, secondarySize);
            }
        }
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        voxel->prevGlueClusterId = currentClusterId[i];
        if (currentClusterId[i] >= 0) {
            voxel->prevGlueClusterSize = currentClusterSize[currentClusterId[i]];
            voxel->prevGlueClusterTag = voxel->debugClusterTag;
            voxel->prevGlueClusterValid = 1;
        } else {
            voxel->prevGlueClusterSize = 0;
            voxel->prevGlueClusterTag = 0;
            voxel->prevGlueClusterValid = 0;
        }
    }
    dynamicGlueClustersInitialized = true;
}

static void set_voxel_velocity(Voxel *voxel, Vector3 vel)
{
    if (!voxel) {
        return;
    }
    voxel->vel = vel;
    for (int i = 0; i < 8; ++i) {
        Particle *p = &voxel->particles[i];
        if (p->inv_mass > 0.0f) {
            p->vel = vel;
        } else {
            p->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
        }
    }
}

static bool voxel_center_near_grid(const Voxel *voxel, float epsilon)
{
    if (!voxel) {
        return false;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    float target_x = 0.5f * ((float)minx + (float)maxx + 1.0f) * VOXEL_SIZE;
    float target_y = 0.5f * ((float)miny + (float)maxy + 1.0f) * VOXEL_SIZE;
    float target_z = 0.5f * ((float)minz + (float)maxz + 1.0f) * VOXEL_SIZE;
    return (fabsf(voxel->pos.x - target_x) <= epsilon) &&
           (fabsf(voxel->pos.y - target_y) <= epsilon) &&
           (fabsf(voxel->pos.z - target_z) <= epsilon);
}

static bool batch_glued_dynamic_voxels(void)
{
    if (voxel_count <= 0) {
        return false;
    }

    unsigned char clusterVisited[MAX_VOXELS] = { 0 };
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *seed = &voxels[i];
        if (!seed->simulate || clusterVisited[i]) {
            continue;
        }

        int cluster_count = build_glue_cluster_indices(i, glueClusterIndices);
        if (cluster_count <= 1) {
            clusterVisited[i] = 1;
            continue;
        }

        for (int c = 0; c < cluster_count; ++c) {
            int idx = glueClusterIndices[c];
            if (idx >= 0 && idx < voxel_count) {
                clusterVisited[idx] = 1;
            }
        }

        UnitVoxelBuffer buffer = { 0 };
        buffer.count = 0;
        Vector3 sum_vel = { 0.0f, 0.0f, 0.0f };
        float vel_weight = 0.0f;
        int min_sleep = INT_MAX;
        bool eligible = true;
        float align_epsilon = VOXEL_SIZE * 0.1f;
        int cluster_activator = -1;

        for (int c = 0; c < cluster_count; ++c) {
            int idx = glueClusterIndices[c];
            if (idx < 0 || idx >= voxel_count) {
                eligible = false;
                break;
            }
            Voxel *voxel = &voxels[idx];
            if (!voxel->simulate) {
                eligible = false;
                break;
            }
            if (cluster_activator < 0 && voxel->activator >= 0) {
                cluster_activator = voxel->activator;
            }
            if (!voxel_center_near_grid(voxel, align_epsilon)) {
                eligible = false;
                break;
            }

            int minx, maxx, miny, maxy, minz, maxz;
            voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
            int span_count_x = maxx - minx + 1;
            int span_count_y = maxy - miny + 1;
            int span_count_z = maxz - minz + 1;
            int cell_count = span_count_x * span_count_y * span_count_z;
            if (cell_count <= 0) {
                eligible = false;
                break;
            }

            for (int gx = minx; gx <= maxx; ++gx) {
                for (int gy = miny; gy <= maxy; ++gy) {
                    for (int gz = minz; gz <= maxz; ++gz) {
                        if (!unit_voxel_buffer_push(&buffer, gx, gy, gz,
                                                    voxel->color, voxel->type, voxel->fixed, -1,
                                                    voxel->debugClusterTag, -1)) {
                            eligible = false;
                            break;
                        }
                    }
                    if (!eligible) {
                        break;
                    }
                }
                if (!eligible) {
                    break;
                }
            }
            if (!eligible) {
                break;
            }

            sum_vel = v_add(sum_vel, v_mul(voxel->vel, (float)cell_count));
            vel_weight += (float)cell_count;
            if (voxel->sleepFrames < min_sleep) {
                min_sleep = voxel->sleepFrames;
            }
        }

        if (!eligible || buffer.count <= 1) {
            continue;
        }
        if (cluster_activator >= 0) {
            for (int b = 0; b < buffer.count; ++b) {
                buffer.voxels[b].activator = cluster_activator;
            }
        }

        if (vel_weight <= 0.0f) {
            vel_weight = 1.0f;
        }
        Vector3 avg_vel = v_mul(sum_vel, 1.0f / vel_weight);
        if (min_sleep == INT_MAX) {
            min_sleep = 0;
        }

        Voxel *snapshots = (Voxel *)malloc(sizeof(Voxel) * (size_t)cluster_count);
        int *sorted = (int *)malloc(sizeof(int) * (size_t)cluster_count);
        if (!snapshots || !sorted) {
            free(snapshots);
            free(sorted);
            ++i;
            continue;
        }
        for (int c = 0; c < cluster_count; ++c) {
            int idx = glueClusterIndices[c];
            snapshots[c] = voxels[idx];
            sorted[c] = idx;
        }
        for (int a = 0; a < cluster_count - 1; ++a) {
            for (int b = a + 1; b < cluster_count; ++b) {
                if (sorted[a] < sorted[b]) {
                    int tmp = sorted[a];
                    sorted[a] = sorted[b];
                    sorted[b] = tmp;
                }
            }
        }
        for (int c = 0; c < cluster_count; ++c) {
            remove_voxel_index(sorted[c]);
        }

        int before = voxel_count;
        int spawned = emit_multiscale_voxels_from_units(&buffer, false, true, -1);
        if (spawned <= 0) {
            for (int c = 0; c < cluster_count; ++c) {
                restore_dynamic_snapshot(&snapshots[c]);
            }
            free(sorted);
            free(snapshots);
            continue;
        }

        for (int v = before; v < voxel_count; ++v) {
            set_voxel_velocity(&voxels[v], avg_vel);
            voxels[v].sleepFrames = min_sleep;
        }

        free(sorted);
        free(snapshots);
        rebuild_voxel_hash();
        rebuild_all_voxel_surfaces();
        rebuild_glue_constraints();
        meshDirty = true;
        return true;
    }

    return false;
}

static bool restore_glue_cluster_to_static(const int *cluster, int cluster_count)
{
    if (!cluster || cluster_count <= 0) {
        return false;
    }

    int voxel_count_before = voxel_count;
    Voxel *snapshots = (Voxel *)malloc(sizeof(Voxel) * (size_t)cluster_count);
    int *sorted = (int *)malloc(sizeof(int) * (size_t)cluster_count);
    if (!snapshots || !sorted) {
        free(snapshots);
        free(sorted);
        return false;
    }

    for (int i = 0; i < cluster_count; ++i) {
        int idx = cluster[i];
        if (idx < 0 || idx >= voxel_count) {
            free(snapshots);
            free(sorted);
            return false;
        }
        snapshots[i] = voxels[idx];
        sorted[i] = idx;
    }

    for (int i = 0; i < cluster_count - 1; ++i) {
        for (int j = i + 1; j < cluster_count; ++j) {
            if (sorted[i] < sorted[j]) {
                int tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    for (int i = 0; i < cluster_count; ++i) {
        remove_voxel_index(sorted[i]);
    }
    int voxel_count_after_removal = voxel_count;
    int static_success = 0;
    bool converted = false;
    for (int i = 0; i < cluster_count; ++i) {
        if (spawn_static_covering_voxel(&snapshots[i])) {
            converted = true;
            ++static_success;
        } else {
            if (debugLogRestoreFailures) {
                TraceLog(LOG_WARNING,
                         "[RestoreFail] static restore failed voxel=%d span=%d rest=(%d..%d,%d..%d,%d..%d)",
                         i, snapshots[i].span,
                         snapshots[i].rest_min_gx, snapshots[i].rest_max_gx,
                         snapshots[i].rest_min_gy, snapshots[i].rest_max_gy,
                         snapshots[i].rest_min_gz, snapshots[i].rest_max_gz);
            }
            restore_dynamic_snapshot(&snapshots[i]);
        }
    }

    if (debugLogRestoreClusters) {
        TraceLog(LOG_INFO,
                 "[RestoreCluster] count=%d voxel_count before=%d after_removal=%d after_restore=%d static_success=%d",
                 cluster_count, voxel_count_before, voxel_count_after_removal, voxel_count,
                 static_success);
    }
    free(sorted);
    free(snapshots);
    return converted;
}

static float compute_glue_relative_velocity_score(int voxel_idx, const Voxel *voxel)
{
    if (voxel_idx < 0 || voxel_idx >= voxel_count || !voxel) {
        return 0.0f;
    }
    int neighbors[MAX_FACE_NEIGHBORS];
    int neighbor_count = gather_glued_neighbors(voxel_idx, neighbors, MAX_FACE_NEIGHBORS);
    float maxRelativeSpeed = 0.0f;
    for (int i = 0; i < neighbor_count; ++i) {
        int nidx = neighbors[i];
        if (nidx < 0 || nidx >= voxel_count) {
            continue;
        }
        const Voxel *neighbor = &voxels[nidx];
        if (!neighbor->simulate) {
            continue;
        }
        Vector3 rel = v_sub(voxel->vel, neighbor->vel);
        float relSpeed = v_length(rel);
        if (relSpeed > maxRelativeSpeed) {
            maxRelativeSpeed = relSpeed;
        }
    }
    return saturatef(maxRelativeSpeed / ACTIVATION_GLUE_REF_SPEED);
}

static void update_dynamic_activation_beliefs(void)
{
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->activationCooldownFrames > 0) {
            voxel->activationCooldownFrames--;
        }
        if (!voxel->simulate) {
            voxel->activationBelief = 0.0f;
            continue;
        }
        if (voxel->isBullet) {
            voxel->activationBelief = 0.0f;
            continue;
        }
        float speed = v_length(voxel->vel);
        float velocityScore = saturatef(speed / ACTIVATION_VELOCITY_REF_SPEED);

        float strain = 0.0f;
        float shear = 0.0f;
        voxel_measure_strain(voxel, &strain, &shear);
        float strainScore = saturatef(fmaxf(strain, shear) / ACTIVATION_STRAIN_REF);

        float glueScore = compute_glue_relative_velocity_score(i, voxel);

        float weightSum = ACTIVATION_VELOCITY_WEIGHT +
                          ACTIVATION_STRAIN_WEIGHT +
                          ACTIVATION_GLUE_WEIGHT;
        if (weightSum <= 0.0f) {
            weightSum = 1.0f;
        }
        float weighted = velocityScore * ACTIVATION_VELOCITY_WEIGHT +
                         strainScore * ACTIVATION_STRAIN_WEIGHT +
                         glueScore * ACTIVATION_GLUE_WEIGHT;
        voxel->activationBelief = saturatef(weighted / weightSum);
    }
}

static bool ranges_overlap(int minA, int maxA, int minB, int maxB) {
    if (maxA < minB) return false;
    if (maxB < minA) return false;
    return true;
}

static bool voxels_face_direction(int idxA, int idxB,
                                  const GlueDirection **out_dir,
                                  bool *a_is_negative)
{
    const Voxel *a = &voxels[idxA];
    const Voxel *b = &voxels[idxB];
    int a_minx, a_maxx, a_miny, a_maxy, a_minz, a_maxz;
    int b_minx, b_maxx, b_miny, b_maxy, b_minz, b_maxz;
    voxel_grid_bounds(a, &a_minx, &a_maxx, &a_miny, &a_maxy, &a_minz, &a_maxz);
    voxel_grid_bounds(b, &b_minx, &b_maxx, &b_miny, &b_maxy, &b_minz, &b_maxz);

    if (ranges_overlap(a_miny, a_maxy, b_miny, b_maxy) &&
        ranges_overlap(a_minz, a_maxz, b_minz, b_maxz))
    {
        if (a_maxx + 1 == b_minx) {
            if (out_dir) *out_dir = &glueDirections[0];
            if (a_is_negative) *a_is_negative = true;
            return true;
        }
        if (b_maxx + 1 == a_minx) {
            if (out_dir) *out_dir = &glueDirections[0];
            if (a_is_negative) *a_is_negative = false;
            return true;
        }
    }

    if (ranges_overlap(a_minx, a_maxx, b_minx, b_maxx) &&
        ranges_overlap(a_minz, a_maxz, b_minz, b_maxz))
    {
        if (a_maxy + 1 == b_miny) {
            if (out_dir) *out_dir = &glueDirections[1];
            if (a_is_negative) *a_is_negative = true;
            return true;
        }
        if (b_maxy + 1 == a_miny) {
            if (out_dir) *out_dir = &glueDirections[1];
            if (a_is_negative) *a_is_negative = false;
            return true;
        }
    }

    if (ranges_overlap(a_minx, a_maxx, b_minx, b_maxx) &&
        ranges_overlap(a_miny, a_maxy, b_miny, b_maxy))
    {
        if (a_maxz + 1 == b_minz) {
            if (out_dir) *out_dir = &glueDirections[2];
            if (a_is_negative) *a_is_negative = true;
            return true;
        }
        if (b_maxz + 1 == a_minz) {
            if (out_dir) *out_dir = &glueDirections[2];
            if (a_is_negative) *a_is_negative = false;
            return true;
        }
    }

    return false;
}

static void rebuild_glue_between_pair(int idxA, int idxB) {
    if (!voxels[idxA].glueEligible || !voxels[idxB].glueEligible) {
        return;
    }
    const GlueDirection *dir = NULL;
    bool a_is_negative = false;
    if (!voxels_face_direction(idxA, idxB, &dir, &a_is_negative) || !dir) {
        return;
    }
    if (a_is_negative) {
        add_bilinear_glue_constraints_for_pair(idxA, idxB, dir);
    } else {
        add_bilinear_glue_constraints_for_pair(idxB, idxA, dir);
    }
}

static void rebuild_glue_children_with_neighbors(const int *children, int child_count,
                                                 const int *neighbors, int neighbor_count)
{
    if (!children || !neighbors) {
        return;
    }
    for (int n = 0; n < neighbor_count; ++n) {
        int neighbor_idx = neighbors[n];
        if (neighbor_idx < 0 || neighbor_idx >= voxel_count) {
            continue;
        }
        for (int c = 0; c < child_count; ++c) {
            int child_idx = children[c];
            if (child_idx < 0 || child_idx >= voxel_count || child_idx == neighbor_idx) {
                continue;
            }
            rebuild_glue_between_pair(child_idx, neighbor_idx);
        }
    }
}

static void compact_glue_constraints(void) {
    int write = 0;
    for (int g = 0; g < glueConstraintCount; ++g) {
        if (!glueConstraints[g].active) {
            continue;
        }
        if (write != g) {
            glueConstraints[write] = glueConstraints[g];
        }
        ++write;
    }
    glueConstraintCount = write;
}

// Returns true when the provided voxel/corner pair is part of an active glue constraint.
static bool particles_are_glued_pair(int voxel_idx_a, int corner_idx_a,
                                     int voxel_idx_b, int corner_idx_b) {
    uint8_t bitA = (uint8_t)(1u << corner_idx_a);
    uint8_t bitB = (uint8_t)(1u << corner_idx_b);

    for (int g = 0; g < glueConstraintCount; ++g) {
        const GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }

        if (gc->coarseVoxel == voxel_idx_a && gc->fineVoxel == voxel_idx_b) {
            if ((gc->coarseMask & bitA) && (gc->fineMask & bitB)) {
                return true;
            }
        } else if (gc->coarseVoxel == voxel_idx_b && gc->fineVoxel == voxel_idx_a) {
            if ((gc->coarseMask & bitB) && (gc->fineMask & bitA)) {
                return true;
            }
        }
    }

    return false;
}

static int gather_neighbor_voxels(const Voxel *voxel, int voxel_idx, int *out, int max_out)
{
    if (max_out <= 0) {
        return 0;
    }

    int count = 0;
    out[count++] = voxel_idx;

    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);

    for (int x = minx - 1; x <= maxx + 1; ++x) {
        for (int y = miny - 1; y <= maxy + 1; ++y) {
            for (int z = minz - 1; z <= maxz + 1; ++z) {
                int idx = table_get(x, y, z);
                if (idx < 0) {
                    continue;
                }
                bool seen = false;
                for (int n = 0; n < count; ++n) {
                    if (out[n] == idx) {
                        seen = true;
                        break;
                    }
                }
                if (!seen && count < max_out) {
                    out[count++] = idx;
                }
            }
        }
    }

    return count;
}

static bool list_contains_index(const int *list, int count, int value)
{
    if (!list || count <= 0) {
        return false;
    }
    for (int i = 0; i < count; ++i) {
        if (list[i] == value) {
            return true;
        }
    }
    return false;
}

static int gather_static_voxels_near_point(Vector3 point, float radius, int *out, int max_out)
{
    if (!out || max_out <= 0) {
        return 0;
    }
    float minx = point.x - radius;
    float maxx = point.x + radius;
    float miny = point.y - radius;
    float maxy = point.y + radius;
    float minz = point.z - radius;
    float maxz = point.z + radius;

    int gx0 = (int)floorf(minx / VOXEL_SIZE);
    int gx1 = (int)floorf(maxx / VOXEL_SIZE);
    int gy0 = (int)floorf(miny / VOXEL_SIZE);
    int gy1 = (int)floorf(maxy / VOXEL_SIZE);
    int gz0 = (int)floorf(minz / VOXEL_SIZE);
    int gz1 = (int)floorf(maxz / VOXEL_SIZE);

    int count = 0;
    for (int z = gz0; z <= gz1; ++z) {
        for (int y = gy0; y <= gy1; ++y) {
            for (int x = gx0; x <= gx1; ++x) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx >= voxel_count) {
                    continue;
                }
                const Voxel *candidate = &voxels[idx];
                if (candidate->simulate) {
                    continue;
                }
                if (list_contains_index(out, count, idx)) {
                    continue;
                }
                if (count < max_out) {
                    out[count++] = idx;
                } else {
                    return count;
                }
            }
        }
    }
    return count;
}

static int gather_static_voxels_near_voxel(const Voxel *voxel, int *out, int max_out)
{
    if (!voxel || !out || max_out <= 0) {
        return 0;
    }
    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    if (minx > maxx || miny > maxy || minz > maxz) {
        return 0;
    }
    minx -= 1; maxx += 1;
    miny -= 1; if (miny < 0) miny = 0;
    maxy += 1;
    minz -= 1; maxz += 1;

    int count = 0;
    for (int z = minz; z <= maxz; ++z) {
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                int idx = table_get(x, y, z);
                if (idx < 0 || idx >= voxel_count) {
                    continue;
                }
                const Voxel *candidate = &voxels[idx];
                if (candidate->simulate) {
                    continue;
                }
                if (list_contains_index(out, count, idx)) {
                    continue;
                }
                if (count < max_out) {
                    out[count++] = idx;
                } else {
                    return count;
                }
            }
        }
    }
    return count;
}

static bool push_particle_out_of_static(const Voxel *static_voxel, Particle *particle, float radius)
{
    if (!static_voxel || !particle || radius <= 0.0f) {
        return false;
    }
    VoxelWorldBounds bounds;
    voxel_world_bounds(static_voxel, &bounds);
    Vector3 pos = particle->predicted_pos;

    bool inside =
        (pos.x >= bounds.minx && pos.x <= bounds.maxx) &&
        (pos.y >= bounds.miny && pos.y <= bounds.maxy) &&
        (pos.z >= bounds.minz && pos.z <= bounds.maxz);

    const float eps = 1e-6f;
    if (!inside) {
        Vector3 closest = {
            clampf(pos.x, bounds.minx, bounds.maxx),
            clampf(pos.y, bounds.miny, bounds.maxy),
            clampf(pos.z, bounds.minz, bounds.maxz)
        };
        Vector3 delta = v_sub(pos, closest);
        float dist_sq = v_dot(delta, delta);
        float radius_sq = radius * radius;
        if (dist_sq >= radius_sq) {
            return false;
        }
        float dist = sqrtf(fmaxf(dist_sq, eps));
        float penetration = radius - dist;
        if (penetration <= 0.0f) {
            return false;
        }
        Vector3 normal = (dist > eps)
            ? v_mul(delta, 1.0f / dist)
            : (Vector3){ 0.0f, 1.0f, 0.0f };
        particle->predicted_pos = v_add(particle->predicted_pos,
                                        v_mul(normal, penetration));
        return true;
    } else {
        float topGap = bounds.maxy - pos.y;
        float bottomGap = pos.y - bounds.miny;
        const float verticalBias = VOXEL_SIZE * 0.25f;
        bool preferTop = (pos.y >= static_voxel->pos.y) &&
                         (topGap <= radius + verticalBias);
        bool preferBottom = (pos.y < static_voxel->pos.y) &&
                            (bottomGap <= radius + verticalBias);
        if (preferTop || preferBottom) {
            bool top = preferTop && (!preferBottom || topGap <= bottomGap);
            Vector3 normal = top ? (Vector3){ 0.0f,  1.0f, 0.0f }
                                 : (Vector3){ 0.0f, -1.0f, 0.0f };
            float gap = top ? topGap : bottomGap;
            float penetration = radius - gap;
            if (penetration < 0.0f) {
                penetration = radius;
            }
            particle->predicted_pos = v_add(particle->predicted_pos,
                                            v_mul(normal, penetration));
            return true;
        }

        float distances[6] = {
            pos.x - bounds.minx,
            bounds.maxx - pos.x,
            pos.y - bounds.miny,
            bounds.maxy - pos.y,
            pos.z - bounds.minz,
            bounds.maxz - pos.z
        };
        int min_idx = 0;
        float min_dist = distances[0];
        for (int k = 1; k < 6; ++k) {
            if (distances[k] < min_dist) {
                min_dist = distances[k];
                min_idx = k;
            }
        }
        Vector3 normal = { 0.0f, 0.0f, 0.0f };
        switch (min_idx) {
            case 0: normal.x = -1.0f; break;
            case 1: normal.x =  1.0f; break;
            case 2: normal.y = -1.0f; break;
            case 3: normal.y =  1.0f; break;
            case 4: normal.z = -1.0f; break;
            case 5: normal.z =  1.0f; break;
            default: normal.y = 1.0f; break;
        }
        float penetration = radius + min_dist;
        if (penetration <= 0.0f) {
            return false;
        }
        particle->predicted_pos = v_add(particle->predicted_pos,
                                        v_mul(normal, penetration));
        return true;
    }
}

static bool resolve_span_static_overlap(int dynamic_idx, Voxel *dynamic,
                                        int static_idx, const Voxel *static_voxel)
{
    if (!dynamic || !static_voxel || !dynamic->simulate || dynamic->span <= 1) {
        return false;
    }
    VoxelWorldBounds boundsA, boundsB;
    voxel_world_bounds(dynamic, &boundsA);
    voxel_world_bounds(static_voxel, &boundsB);

    float overlapX = fminf(boundsA.maxx, boundsB.maxx) - fmaxf(boundsA.minx, boundsB.minx);
    float overlapY = fminf(boundsA.maxy, boundsB.maxy) - fmaxf(boundsA.miny, boundsB.miny);
    float overlapZ = fminf(boundsA.maxz, boundsB.maxz) - fmaxf(boundsA.minz, boundsB.minz);
    if (overlapX <= 0.0f || overlapY <= 0.0f || overlapZ <= 0.0f) {
        return false;
    }

    float overlaps[3] = { overlapX, overlapY, overlapZ };
    int axis = 0;
    if (overlapY < overlaps[axis]) axis = 1;
    if (overlapZ < overlaps[axis]) axis = 2;
    float penetration = overlaps[axis];
    if (penetration <= 0.0f) {
        return false;
    }

    Vector3 centerA = dynamic->pos;
    Vector3 centerB = static_voxel->pos;
    float inv_mass_sum = 0.0f;
    compute_voxel_center_and_mass(dynamic, &centerA, &inv_mass_sum);
    if (inv_mass_sum <= 0.0f) {
        return false;
    }

    Vector3 normal = { 0.0f, 0.0f, 0.0f };
    if (axis == 0) {
        normal.x = (centerA.x >= centerB.x) ? 1.0f : -1.0f;
    } else if (axis == 1) {
        normal.y = (centerA.y >= centerB.y) ? 1.0f : -1.0f;
    } else {
        normal.z = (centerA.z >= centerB.z) ? 1.0f : -1.0f;
    }

    if (debugLogSpanCollisions) {
        TraceLog(LOG_INFO,
                 "[SpanStatic] overlap dyn=%d span=%d static=%d axis=%d pen=%.4f overlaps=(%.4f,%.4f,%.4f) "
                 "centersA=(%.2f,%.2f,%.2f) centersB=(%.2f,%.2f,%.2f) normal=(%.1f,%.1f,%.1f)",
                 dynamic_idx, dynamic->span, static_idx, axis, penetration,
                 overlapX, overlapY, overlapZ,
                 centerA.x, centerA.y, centerA.z,
                 centerB.x, centerB.y, centerB.z,
                 normal.x, normal.y, normal.z);
    }

    const float omega = COLLISION_RELAXATION;
    Vector3 delta = v_mul(normal, penetration * omega);
    for (int j = 0; j < 8; ++j) {
        Particle *p = &dynamic->particles[j];
        if (p->inv_mass <= 0.0f) {
            continue;
        }
        float weight = p->inv_mass / inv_mass_sum;
        p->predicted_pos = v_add(p->predicted_pos, v_mul(delta, weight));
    }

    return true;
}

static bool nudge_voxel_bottom_above_static(int voxel_idx, Voxel *voxel)
{
    if (!voxel || !voxel->simulate || voxel->span <= 1) {
        return false;
    }

    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    int support_y = miny - 1;
    if (support_y < 0) {
        return false;
    }

    float half_edge = 0.5f * voxel->rest_edge;
    Vector3 center = voxel->pos;
    compute_voxel_center_and_mass(voxel, &center, NULL);
    float current_bottom = center.y - half_edge;
    float required_bottom = current_bottom;
    bool found_support = false;

    for (int z = minz; z <= maxz; ++z) {
        for (int x = minx; x <= maxx; ++x) {
            int idx = table_get(x, support_y, z);
            if (idx < 0 || idx >= voxel_count) {
                continue;
            }
            const Voxel *support = &voxels[idx];
            if (support->simulate) {
                continue;
            }
            float top = ((float)(support_y + 1)) * VOXEL_SIZE;
            if (top > required_bottom) {
                required_bottom = top;
            }
            found_support = true;
        }
    }

    float lift = required_bottom - current_bottom;
    if (!found_support || lift <= 0.0f) {
        return false;
    }

    Vector3 delta = { 0.0f, lift, 0.0f };
    voxel->pos = v_add(voxel->pos, delta);
    for (int j = 0; j < 8; ++j) {
        voxel->particles[j].predicted_pos = v_add(voxel->particles[j].predicted_pos, delta);
    }

    if (debugLogSpanCollisions) {
        TraceLog(LOG_INFO,
                 "[SpanStaticLift] voxel=%d span=%d lift=%.4f bottom=%.4f target=%.4f",
                 voxel_idx, voxel->span, lift, current_bottom, required_bottom);
    }
    return true;
}

static void glue_dynamic_voxel_to_static_neighbors_for_voxel(int voxel_idx)
{
    if (voxel_idx < 0 || voxel_idx >= voxel_count) {
        return;
    }
    Voxel *dynamic = &voxels[voxel_idx];
    if (!dynamic->simulate || !dynamic->glueEligible) {
        return;
    }

    int processed[MAX_FACE_NEIGHBORS];
    int processed_count = 0;

    int minx, maxx, miny, maxy, minz, maxz;
    voxel_grid_bounds(dynamic, &minx, &maxx, &miny, &maxy, &minz, &maxz);

    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            int idx = table_get(maxx + 1, y, z);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
    for (int y = miny; y <= maxy; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            int idx = table_get(minx - 1, y, z);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
    for (int x = minx; x <= maxx; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            int idx = table_get(x, maxy + 1, z);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
    for (int x = minx; x <= maxx; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            int idx = table_get(x, miny - 1, z);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            int idx = table_get(x, y, maxz + 1);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            int idx = table_get(x, y, minz - 1);
            if (idx < 0 || idx >= voxel_count) continue;
            Voxel *neighbor = &voxels[idx];
            if (neighbor->simulate || list_contains_index(processed, processed_count, idx)) continue;
            bool prevEligible = neighbor->glueEligible;
            neighbor->glueEligible = true;
            rebuild_glue_between_pair(voxel_idx, idx);
            neighbor->glueEligible = prevEligible;
            if (processed_count < MAX_FACE_NEIGHBORS) processed[processed_count++] = idx;
        }
    }
}

static void glue_dynamic_voxel_to_static_neighbors(void)
{
    for (int i = 0; i < voxel_count; ++i) {
        if (!voxels[i].simulate) {
            continue;
        }
        glue_dynamic_voxel_to_static_neighbors_for_voxel(i);
    }
}

// True when two voxels share any active glue constraint (face coupling).
static bool voxels_are_glued(int voxel_idx_a, int voxel_idx_b) {
    if (voxel_idx_a == voxel_idx_b) {
        return false;
    }

    rebuild_glue_adjacency_if_dirty();
    if (voxel_idx_a < 0 || voxel_idx_a >= voxel_count) {
        return false;
    }
    int count = gluedNeighborCounts[voxel_idx_a];
    for (int i = 0; i < count; ++i) {
        if (gluedNeighborList[voxel_idx_a][i] == voxel_idx_b) {
            return true;
        }
    }
    return false;
}

// Edge/corner adjacency causes constraints to overlap; skip collisions for those pairs as well.
// Returns false when intervals are separated by more than eps.
// Sets *touch when they meet (within eps) but do not overlap.
static bool voxels_share_edge_or_corner(const Voxel *voxel_a, const Voxel *voxel_b) {
    const float eps = VOXEL_SIZE * 0.05f;
    int overlap_axes = 0;
    int touching_axes = voxel_touching_axes(voxel_a, voxel_b, &overlap_axes, eps);
    bool share_edge = (touching_axes >= 2);
    bool rest_share = voxels_share_edge_or_corner_rest(voxel_a, voxel_b);
    bool rest_face = voxels_share_face_rest(voxel_a, voxel_b);
    

    if (debug_should_log_span_pair(voxel_a, voxel_b, &debugSpanEdgeLogBudget)) {
        TraceLog(LOG_INFO,
                 "[SpanDebug] adjacency shareEdge=%s spans=(%d,%d) touchingAxes=%d overlapAxes=%d "
                 "restShare=%s restFace=%s posA=(%.2f,%.2f,%.2f) posB=(%.2f,%.2f,%.2f)",
                 share_edge ? "true" : "false",
                 voxel_a->span, voxel_b->span,
                 touching_axes, overlap_axes,
                 rest_share ? "true" : "false",
                 rest_face ? "true" : "false",
                 voxel_a->pos.x, voxel_a->pos.y, voxel_a->pos.z,
                 voxel_b->pos.x, voxel_b->pos.y, voxel_b->pos.z);
    }

    return share_edge || rest_share || rest_face;
}


static void compute_voxel_center_and_mass(const Voxel *voxel, Vector3 *center, float *inv_mass_sum) {
    Vector3 c = { 0.0f, 0.0f, 0.0f };
    float sum = 0.0f;

    for (int i = 0; i < 8; ++i) {
        const Particle *p = &voxel->particles[i];
        c = v_add(c, p->predicted_pos);
        sum += p->inv_mass;
    }

    if (center) {
        *center = v_mul(c, 1.0f / 8.0f);
    }
    if (inv_mass_sum) {
        *inv_mass_sum = sum;
    }
}

// Resolve collisions against the scene and neighbouring voxels (mirrors ResolveCollisions compute pass). need to filter based on glue
static void solve_particle_collisions(float dt) {
    const float half_player = PLAYER_SIZE * 0.5f;
    const float omega = COLLISION_RELAXATION;
    const float eps = 1e-6f;
    const bool centroid_only = (dt > COLLISION_CENTROID_ONLY_DT);
    const int particle_start = centroid_only ? VOXEL_CENTER_INDEX : 0;
    const int particle_end = centroid_only ? (VOXEL_CENTER_INDEX + 1) : VOXEL_PARTICLE_COUNT;
    static int glue_cluster_id[MAX_VOXELS];
    if (skipGlueClusterCollisions) {
        build_glue_cluster_ids(glue_cluster_id);
    }
    if (debugLogSpanCollisions) {
        debugSpanCollisionLogBudget = 32;
        debugSpanEdgeLogBudget = 32;
    } else {
        debugSpanCollisionLogBudget = 0;
        debugSpanEdgeLogBudget = 0;
    }

    // First, clamp predictions against static scene bounds and player capsules.
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!(voxel->simulate) || voxel->type != 0 || voxel->isBullet){
            continue;
        }
        float voxel_radius = voxel_particle_radius(voxel);
        float terrain_limit = FLOOR_SIZE - voxel_radius;
        int span = (voxel->span > 0) ? voxel->span : 1;
        float span_extent = 0.5f * VOXEL_SIZE * (float)(span - 1);
        float static_collision_radius = voxel_radius - span_extent;
        float min_static_radius = 0.5f * VOXEL_SIZE;
        if (static_collision_radius < min_static_radius) {
            static_collision_radius = min_static_radius;
        }

        for (int j = particle_start; j < particle_end; ++j) {
            Particle *p = voxel_particle_at(voxel, j);

            Vector3 pos = p->predicted_pos;

            float floor_offset = 0.5f * VOXEL_SIZE * (float)span;
            float floor_limit = fmaxf(0.0f, floor_offset - voxel_radius);
            if (pos.y < floor_limit) {
                pos.y = floor_limit;
            }

            pos.x = clampf(pos.x, -terrain_limit, terrain_limit);
            pos.z = clampf(pos.z, -terrain_limit, terrain_limit);

            // Interactions with player bounding boxes (kept for gameplay parity).
            for (int player_idx = 0; player_idx < activePlayers; ++player_idx) {
                Player *pl = &players[player_idx];
                Vector3 box_min = {
                    pl->pos.x - half_player,
                    pl->pos.y - half_player,
                    pl->pos.z - half_player
                };
                Vector3 box_max = {
                    pl->pos.x + half_player,
                    pl->pos.y + half_player,
                    pl->pos.z + half_player
                };

                Vector3 nearest = {
                    clampf(pos.x, box_min.x, box_max.x),
                    clampf(pos.y, box_min.y, box_max.y),
                    clampf(pos.z, box_min.z, box_max.z)
                };

                Vector3 delta = v_sub(pos, nearest);
                float dist_sq = v_dot(delta, delta);
                float radius_sq = voxel_radius * voxel_radius;
                if (dist_sq < radius_sq) {
                    float dist = sqrtf(fmaxf(dist_sq, eps));
                    float penetration = voxel_radius - dist;
                    Vector3 normal = (dist > eps)
                        ? v_mul(delta, -1.0f / dist)
                        : (Vector3){ 0.0f, 1.0f, 0.0f };
                    pos = v_add(pos, v_mul(normal, penetration));
                }
            }

            p->predicted_pos = pos;

            int static_neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
            int static_count = gather_static_voxels_near_point(
                p->predicted_pos, voxel_radius,
                static_neighbors, MAX_STATIC_COLLISION_NEIGHBORS);
            for (int s = 0; s < static_count; ++s) {
                int static_idx = static_neighbors[s];
                if (static_idx < 0 || static_idx >= voxel_count) {
                    continue;
                }
                const Voxel *static_voxel = &voxels[static_idx];
                if (static_voxel->simulate) {
                    continue;
                }
                push_particle_out_of_static(static_voxel, p, static_collision_radius);
            }
        }

        if (span > 1) {
            int static_indices[MAX_STATIC_COLLISION_NEIGHBORS];
            int static_count = gather_static_voxels_near_voxel(voxel, static_indices,
                                                               MAX_STATIC_COLLISION_NEIGHBORS);
            for (int s = 0; s < static_count; ++s) {
                int static_idx = static_indices[s];
                if (static_idx < 0 || static_idx >= voxel_count) {
                    continue;
                }
                const Voxel *static_voxel = &voxels[static_idx];
                if (static_voxel->simulate) {
                    continue;
                }
                resolve_span_static_overlap(i, voxel, static_idx, static_voxel);
            }
            nudge_voxel_bottom_above_static(i, voxel);
        }
    }

    // Particle-particle collisions using a symmetric correction identical to the compute shader.
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxelA = &voxels[i];
        if (!(voxelA->simulate) || voxelA->type != 0 || voxelA->isBullet){
            continue;
        }
        float radiusA = voxel_particle_radius(voxelA);
        int spanA = (voxelA->span > 0) ? voxelA->span : 1;
        float spanA_extent = 0.5f * VOXEL_SIZE * (float)(spanA - 1);
        VoxelWorldBounds boundsA;
        voxel_predicted_bounds(voxelA, &boundsA);
        int neighbor_ids[MAX_NEIGHBOR_VOXELS];
        int neighbor_count = gather_neighbor_voxels(voxelA, i, neighbor_ids, MAX_NEIGHBOR_VOXELS);

        for (int j = particle_start; j < particle_end; ++j) {
            Particle *pa = voxel_particle_at(voxelA, j);
            float wa = pa->inv_mass;

            for (int n = 0; n < neighbor_count; ++n) {
                int neighbor_idx = neighbor_ids[n];
                if (neighbor_idx < i) {
                    continue;
                }
                if (neighbor_idx == i) {
                    continue;
                }

                Voxel *voxelB = &voxels[neighbor_idx];
                float radiusB = voxel_particle_radius(voxelB);
                int spanB = (voxelB->span > 0) ? voxelB->span : 1;
                float spanB_extent = 0.5f * VOXEL_SIZE * (float)(spanB - 1);
                VoxelWorldBounds boundsB;
                voxel_predicted_bounds(voxelB, &boundsB);

                if (boundsA.maxx + radiusA < boundsB.minx - radiusB ||
                    boundsB.maxx + radiusB < boundsA.minx - radiusA ||
                    boundsA.maxy + radiusA < boundsB.miny - radiusB ||
                    boundsB.maxy + radiusB < boundsA.miny - radiusA ||
                    boundsA.maxz + radiusA < boundsB.minz - radiusB ||
                    boundsB.maxz + radiusB < boundsA.minz - radiusA) {
                    continue;
                }

                if (!centroid_only && (spanA > 1 || spanB > 1)) {
                    continue;
                }

                bool glued = voxels_are_glued(i, neighbor_idx);
                bool share_edge_corner = voxels_share_edge_or_corner(voxelA, voxelB);
                if (glued || share_edge_corner) {
                    continue;
                }
                if (skipGlueClusterCollisions &&
                    glue_cluster_id[i] >= 0 && glue_cluster_id[i] == glue_cluster_id[neighbor_idx]) {
                    continue;
                }

                for (int q = particle_start; q < particle_end; ++q) {
                    if (neighbor_idx == i && q <= j) {
                        continue;
                    }

                    if (j < VOXEL_CORNER_COUNT && q < VOXEL_CORNER_COUNT) {
                        if (particles_are_glued_pair(i, j, neighbor_idx, q)) {
                            continue;
                        }
                    }

                    Particle *pb = voxel_particle_at(voxelB, q);
                    float wb = pb->inv_mass;

                    float w_sum = wa + wb;
                    if (w_sum <= 0.0f) {
                        continue;
                    }

                    Vector3 delta = v_sub(pa->predicted_pos, pb->predicted_pos);
                    float dist_sq = v_dot(delta, delta);
                    float span_offset = spanA_extent + spanB_extent;
                    float target_dist = radiusA + radiusB - span_offset;
                    if (target_dist < 0.0f) {
                        target_dist = 0.0f;
                    }

                    if (dist_sq >= (target_dist * target_dist)) {
                        continue;
                    }

                    float dist = sqrtf(fmaxf(dist_sq, eps));
                    float penetration = target_dist - dist;
                    if (penetration <= 0.0f) {
                        continue;
                    }

                    if (debug_should_log_span_pair(voxelA, voxelB, &debugSpanCollisionLogBudget)) {
                        int overlap_axes = 0;
                        int touching_axes = voxel_touching_axes(
                            voxelA, voxelB, &overlap_axes, VOXEL_SIZE * 0.05f);
                        TraceLog(LOG_INFO,
                                 "[SpanDebug] collision voxels=(%d,%d) particles=(%d,%d) spans=(%d,%d) "
                                 "touchingAxes=%d overlapAxes=%d glued=%d shareEdge=%d dist=%.4f target=%.4f pen=%.4f",
                                 i, neighbor_idx, j, q,
                                 voxelA->span, voxelB->span,
                                 touching_axes, overlap_axes,
                                 glued ? 1 : 0, share_edge_corner ? 1 : 0,
                                 dist, target_dist, penetration);
                    }

                    Vector3 normal = (dist > eps)
                        ? v_mul(delta, 1.0f / dist)
                        : (Vector3){ 1.0f, 0.0f, 0.0f };

                    float h = 0.5f * penetration;
                    float scale = omega * h / w_sum;

                    if (wa > 0.0f) {
                        pa->predicted_pos = v_add(pa->predicted_pos, v_mul(normal, scale * wa));
                    }
                    if (wb > 0.0f) {
                        pb->predicted_pos = v_sub(pb->predicted_pos, v_mul(normal, scale * wb));
                    }
                }
            }
        }
    }
}

static void solve_span_voxel_collisions(void) {
    const float contact_eps = 0.0f;
    const float omega = COLLISION_RELAXATION;
    static int glue_cluster_id[MAX_VOXELS];
    if (skipGlueClusterCollisions) {
        build_glue_cluster_ids(glue_cluster_id);
    }

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxelA = &voxels[i];
        if (!voxelA->simulate || voxelA->type != 0 || voxelA->isBullet) {
            continue;
        }
        VoxelWorldBounds boundsA_pred;
        voxel_predicted_bounds(voxelA, &boundsA_pred);
        int neighbor_ids[MAX_NEIGHBOR_VOXELS];
        int neighbor_count = gather_neighbor_voxels(voxelA, i, neighbor_ids, MAX_NEIGHBOR_VOXELS);

        for (int n = 0; n < neighbor_count; ++n) {
            int neighbor_idx = neighbor_ids[n];
            if (neighbor_idx <= i) {
                continue;
            }

            Voxel *voxelB = &voxels[neighbor_idx];
            int spanA = (voxelA->span > 0) ? voxelA->span : 1;
            int spanB = (voxelB->span > 0) ? voxelB->span : 1;
            if (spanA <= 1 && spanB <= 1) {
                continue;
            }
            if (voxels_are_glued(i, neighbor_idx)) {
                continue;
            }
            if (voxels_share_edge_or_corner(voxelA, voxelB)) {
                continue;
            }
            if (skipGlueClusterCollisions &&
                glue_cluster_id[i] >= 0 && glue_cluster_id[i] == glue_cluster_id[neighbor_idx]) {
                continue;
            }

            VoxelWorldBounds boundsB_pred;
            voxel_predicted_bounds(voxelB, &boundsB_pred);
            if (boundsA_pred.maxx < boundsB_pred.minx - contact_eps ||
                boundsB_pred.maxx < boundsA_pred.minx - contact_eps ||
                boundsA_pred.maxy < boundsB_pred.miny - contact_eps ||
                boundsB_pred.maxy < boundsA_pred.miny - contact_eps ||
                boundsA_pred.maxz < boundsB_pred.minz - contact_eps ||
                boundsB_pred.maxz < boundsA_pred.minz - contact_eps) {
                continue;
            }

            VoxelWorldBounds boundsA, boundsB;
            voxel_world_bounds(voxelA, &boundsA);
            voxel_world_bounds(voxelB, &boundsB);

            float overlapX = fminf(boundsA.maxx, boundsB.maxx) - fmaxf(boundsA.minx, boundsB.minx);
            float overlapY = fminf(boundsA.maxy, boundsB.maxy) - fmaxf(boundsA.miny, boundsB.miny);
            float overlapZ = fminf(boundsA.maxz, boundsB.maxz) - fmaxf(boundsA.minz, boundsB.minz);
            if (overlapX <= contact_eps || overlapY <= contact_eps || overlapZ <= contact_eps) {
                continue;
            }

            float overlaps[3] = { overlapX, overlapY, overlapZ };
            int axis = 0;
            if (overlapY < overlaps[axis]) {
                axis = 1;
            }
            if (overlapZ < overlaps[axis]) {
                axis = 2;
            }

            float penetration = overlaps[axis];
            if (penetration <= 0.0f) {
                continue;
            }

            Vector3 centerA = voxelA->pos;
            Vector3 centerB = voxelB->pos;
            float inv_massA = 0.0f;
            float inv_massB = 0.0f;
            compute_voxel_center_and_mass(voxelA, &centerA, &inv_massA);
            compute_voxel_center_and_mass(voxelB, &centerB, &inv_massB);

            Vector3 normal = { 0.0f, 0.0f, 0.0f };
            if (axis == 0) {
                normal.x = (centerA.x >= centerB.x) ? 1.0f : -1.0f;
            } else if (axis == 1) {
                normal.y = (centerA.y >= centerB.y) ? 1.0f : -1.0f;
            } else {
                normal.z = (centerA.z >= centerB.z) ? 1.0f : -1.0f;
            }

            float w_sum = inv_massA + inv_massB;
            if (w_sum <= 0.0f) {
                continue;
            }

            float correction = penetration * omega;
            float moveA = (inv_massA > 0.0f) ? correction * (inv_massA / w_sum) : 0.0f;
            float moveB = (inv_massB > 0.0f) ? correction * (inv_massB / w_sum) : 0.0f;

            if (moveA > 0.0f) {
                Vector3 delta = v_mul(normal, moveA);
                for (int j = 0; j < 8; ++j) {
                    Particle *pa = &voxelA->particles[j];
                    if (pa->inv_mass <= 0.0f) {
                        continue;
                    }
                    float weight = (inv_massA > 0.0f) ? (pa->inv_mass / inv_massA) : 0.0f;
                    pa->predicted_pos = v_add(pa->predicted_pos, v_mul(delta, weight));
                }
            }

            if (moveB > 0.0f) {
                Vector3 delta = v_mul(normal, moveB);
                for (int j = 0; j < 8; ++j) {
                    Particle *pb = &voxelB->particles[j];
                    if (pb->inv_mass <= 0.0f) {
                        continue;
                    }
                    float weight = (inv_massB > 0.0f) ? (pb->inv_mass / inv_massB) : 0.0f;
                    pb->predicted_pos = v_sub(pb->predicted_pos, v_mul(delta, weight));
                }
            }
        }
    }
}

static void update_particle_velocities(float dt) {
    float inv_dt = (dt > 0.0f) ? 1.0f / dt : 0.0f;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->type != 0 || voxel->isBullet) {
            continue;
        }
        bool skipVelocity = (voxel->skipCollisionVelocityFrames > 0);
        if (skipVelocity && voxel->skipCollisionVelocityFrames > 0) {
            voxel->skipCollisionVelocityFrames--;
        }
        Vector3 centroid = { 0.0f, 0.0f, 0.0f };
        Vector3 prev_centroid = { 0.0f, 0.0f, 0.0f };

        for (int j = 0; j < 8; ++j) {
            Particle *p = &voxel->particles[j];
            centroid = v_add(centroid, p->predicted_pos);
            prev_centroid = v_add(prev_centroid, p->prev_pos);

            Vector3 new_pos = p->predicted_pos;
            Vector3 delta = v_sub(new_pos, p->prev_pos);

            if (!skipVelocity) {
                if (p->inv_mass > 0.0f) {
                    p->vel = v_mul(delta, inv_dt);
                } else {
                    p->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
                }
            }

            p->pos = new_pos;
        }

        centroid = v_mul(centroid, 1.0f / 8.0f);
        prev_centroid = v_mul(prev_centroid, 1.0f / 8.0f);
        if (!skipVelocity) {
            voxel->vel = v_mul(v_sub(centroid, prev_centroid), inv_dt);
        }
        voxel->pos = centroid;

        Particle *center = &voxel->center_particle;
        center->prev_pos = prev_centroid;
        center->predicted_pos = centroid;
        center->pos = centroid;
        center->vel = voxel->vel;
    }
}

static void voxel_measure_strain(const Voxel *voxel,
                                 float *out_max_strain,
                                 float *out_max_shear)
{
    float max_strain = 0.0f;
    float max_shear = 0.0f;
    float rest_edge = voxel->rest_edge;
    if (rest_edge <= 0.0f) {
        rest_edge = VOXEL_SIZE;
    }
    float inv_rest_edge = (rest_edge > 0.0f) ? (1.0f / rest_edge) : 0.0f;

    for (int e = 0; e < 12; ++e) {
        int a_idx = voxel_edge_pairs[e][0];
        int b_idx = voxel_edge_pairs[e][1];
        Vector3 a = voxel->particles[a_idx].predicted_pos;
        Vector3 b = voxel->particles[b_idx].predicted_pos;
        float len = v_length(v_sub(a, b));
        float strain = fabsf(len - rest_edge) * inv_rest_edge;
        if (strain > max_strain) {
            max_strain = strain;
        }
    }

    Vector3 axis_vecs[3];
    axis_vecs[0] = v_sub(voxel->particles[1].predicted_pos,
                         voxel->particles[0].predicted_pos);
    axis_vecs[1] = v_sub(voxel->particles[2].predicted_pos,
                         voxel->particles[0].predicted_pos);
    axis_vecs[2] = v_sub(voxel->particles[4].predicted_pos,
                         voxel->particles[0].predicted_pos);
    for (int a = 0; a < 3; ++a) {
        float len = v_length(axis_vecs[a]);
        if (len > 1e-6f) {
            axis_vecs[a] = v_mul(axis_vecs[a], 1.0f / len);
        } else {
            axis_vecs[a] = (Vector3){ 0.0f, 0.0f, 0.0f };
        }
    }
    const int shear_pairs[3][2] = { {0,1}, {0,2}, {1,2} };
    for (int s = 0; s < 3; ++s) {
        int i = shear_pairs[s][0];
        int j = shear_pairs[s][1];
        float dot_abs = fabsf(v_dot(axis_vecs[i], axis_vecs[j]));
        if (dot_abs > max_shear) {
            max_shear = dot_abs;
        }
    }

    if (out_max_strain) {
        *out_max_strain = max_strain;
    }
    if (out_max_shear) {
        *out_max_shear = max_shear;
    }
}

static void apply_uniform_velocity(Voxel *v, Vector3 vel, float dt) {
    v->vel = vel;
    for (int i = 0; i < 8; ++i) {
        v->particles[i].vel = vel;
        if (dt > 0.0f) {
            Vector3 offset = v_mul(vel, dt);
            v->particles[i].prev_pos = v_sub(v->particles[i].pos, offset);
        }
        v->particles[i].predicted_pos = v->particles[i].pos;
    }
    v->center_particle.vel = vel;
    if (dt > 0.0f) {
        Vector3 offset = v_mul(vel, dt);
        v->center_particle.prev_pos = v_sub(v->pos, offset);
    }
    v->center_particle.predicted_pos = v->pos;
    v->center_particle.pos = v->pos;
}

static int split_voxel_at(int idx, float dt, int *out_children, int max_children) {
    if (idx < 0 || idx >= voxel_count || !out_children || max_children < MAX_SPLIT_CHILDREN) {
        return 0;
    }

    Voxel parent = voxels[idx];
    if (parent.span < 2) {
        return 0;
    }
    const int additional_children = MAX_SPLIT_CHILDREN - 1;
    if (voxel_count + additional_children > MAX_VOXELS) {
        return 0;
    }

    int child_span = parent.span / 2;
    if (child_span < 1) {
        return 0;
    }
    float child_edge = VOXEL_SIZE * (float)child_span;
    float parent_half = 0.5f * parent.rest_edge;
    float child_offset = parent_half - 0.5f * child_edge;
    if (child_offset < 0.0f) {
        child_offset = 0.0f;
    }
    Vector3 split_vel = v_mul(parent.vel, SPLIT_VELOCITY_DAMP);

    int child_counter = 0;
    int next_index = voxel_count;
    for (int iz = 0; iz < 2; ++iz) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int ix = 0; ix < 2; ++ix) {
                float cx = parent.pos.x + (ix ? child_offset : -child_offset);
                float cy = parent.pos.y + (iy ? child_offset : -child_offset);
                float cz = parent.pos.z + (iz ? child_offset : -child_offset);

                int child_idx;
                if (child_counter == 0) {
                    child_idx = idx;
                } else {
                    child_idx = next_index++;
                }

                Voxel *child = &voxels[child_idx];
                init_voxel_struct(child,
                                  cx, cy, cz,
                                  parent.fixed, parent.simulate,
                                  parent.color, parent.type,
                                  child_span, parent.owner);
                child->activator = parent.activator;
                child->debugClusterTag = parent.debugClusterTag;
                apply_uniform_velocity(child, split_vel, dt);
                out_children[child_counter] = child_idx;
                child_counter++;
            }
        }
    }

    voxel_count = next_index;
    return child_counter;
}

static bool split_strained_voxels(float dt) {
    bool split_any = false;
    int i = 0;
    while (i < voxel_count) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->type != 0 || voxel->isBullet) {
            ++i;
            continue;
        }

        float strain = 0.0f;
        float shear = 0.0f;
        voxel_measure_strain(voxel, &strain, &shear);
        if (strain > VOXEL_DUST_STRAIN_THRESHOLD ||
            shear > VOXEL_DUST_SHEAR_THRESHOLD)
        {
            int glued_neighbors[MAX_FACE_NEIGHBORS];
            int glued_neighbor_count = gather_glued_neighbors(i, glued_neighbors, MAX_FACE_NEIGHBORS);
            for (int n = 0; n < glued_neighbor_count; ++n) {
                deactivate_glue_constraints_between(i, glued_neighbors[n]);
            }
            if (recycle_queue_push(voxel) && debugLogVoxelRecycle) {
                TraceLog(LOG_INFO,
                         "[Dust] enqueue-voxel idx=%d span=%d strain=%.3f shear=%.3f queue=%d",
                         i, voxel->span, strain, shear, recycleQueueCount);
            }
            remove_voxel_index(i);
            split_any = true;
            continue;
        }
        if (voxel->span <= 1) {
            ++i;
            continue;
        }
        if (strain > VOXEL_SPLIT_STRAIN_THRESHOLD ||
            shear > VOXEL_SPLIT_SHEAR_THRESHOLD)
        {
            int glued_neighbors[MAX_FACE_NEIGHBORS];
            int glued_neighbor_count = gather_glued_neighbors(i, glued_neighbors, MAX_FACE_NEIGHBORS);
            for (int n = 0; n < glued_neighbor_count; ++n) {
                deactivate_glue_constraints_between(i, glued_neighbors[n]);
            }

            int children[MAX_SPLIT_CHILDREN];
            int child_count = split_voxel_at(i, dt, children, MAX_SPLIT_CHILDREN);
            if (child_count > 0) {
                split_any = true;
                rebuild_glue_children_with_neighbors(children, child_count,
                                                     glued_neighbors, glued_neighbor_count);
                ++i;
                continue;
            } else {
                for (int n = 0; n < glued_neighbor_count; ++n) {
                    rebuild_glue_between_pair(i, glued_neighbors[n]);
                }
            }
        }
        ++i;
    }
    if (split_any) {
        compact_glue_constraints();
    }
    return split_any;
}

void simulate_voxel_pbd(float dt) {
    const int substeps = PBD_SUBSTEPS;
    const int constraint_iterations = PBD_CONSTRAINT_ITERS;
    const float sub_dt = (substeps > 0) ? dt / (float)substeps : dt;

    for (int step = 0; step < substeps; ++step) {
        if (debugLogVoxelBlowup) {
            debugBlowupLogBudget = 32;
        }
        integrate_particles(sub_dt);
        solve_particle_collisions(sub_dt);
        solve_span_voxel_collisions();
        bool split_this_step = split_strained_voxels(sub_dt);
        if (split_this_step) {
            rebuild_voxel_hash();
            rebuild_all_voxel_surfaces();
            rebuild_glue_constraints();
            meshDirty = true;
        }

        for (int it = 0; it < constraint_iterations; ++it) {
            for (int i = 0; i < voxel_count; ++i) {
                Voxel *voxel = &voxels[i];
                if (!voxel->simulate || voxel->type != 0 || voxel->isBullet)
                    continue;
                solve_voxel_shape(voxel);
            }
        solve_voxel_glue(false);
        }
        solve_voxel_glue(true);
        update_particle_velocities(sub_dt);
        if (debugLogVoxelBlowup && debugBlowupLogBudget > 0) {
            for (int i = 0; i < voxel_count && debugBlowupLogBudget > 0; ++i) {
                Voxel *voxel = &voxels[i];
                if (!voxel->simulate) {
                    continue;
                }
                Vector3 pos = voxel->pos;
                if (!isfinite(pos.x) || !isfinite(pos.y) || !isfinite(pos.z) ||
                    fabsf(pos.x) > 1e4f || fabsf(pos.y) > 1e4f || fabsf(pos.z) > 1e4f)
                {
                    TraceLog(LOG_WARNING,
                             "[VoxelBlowup] voxel=%d span=%d pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f)",
                             i, voxel->span,
                             pos.x, pos.y, pos.z,
                             voxel->vel.x, voxel->vel.y, voxel->vel.z);
                    --debugBlowupLogBudget;
                }
            }
        }
    }

    log_dynamic_voxel_positions();
}

/*
 * Sketch of the CPU voxel PBD loop, mirroring the paper's GPU pipeline.
 *
 * void simulate_voxel_pbd(float dt) {
 *     const int substeps = get_substep_count();
 *     const float sub_dt = dt / (float)substeps;
 *
 *     for (int step = 0; step < substeps; ++step) {
 *         // 1. Kinematics + collisions for all particles (predict positions).
 *         integrate_particles(sub_dt);
 *         solve_particle_collisions(sub_dt);
 *
 *         // 2. Constraint iterations (Gauss-Seidel order).
 *         for (int it = 0; it < get_constraint_iterations(); ++it) {
 *             // 2a. Voxel Gram-Schmidt shape matching (Algorithm 1).
 *             for_each_voxel(voxel_count, [&](Voxel *v) {
 *                 solve_voxel_shape(v, sub_dt);
 *             });
 *             // Face-constraint pass removed in this CPU prototype.
 *         }
 *
 *         // 3. Velocity update from corrected positions.
 *         update_particle_velocities(sub_dt);
 *     }
 * }
 */

// Fire a voxel bullet
static void FireVoxel(int idx) {
    Player *p = &players[idx];
    float now = (float)GetTime();
    if (now - p->last_shot_time < BULLET_COOLDOWN_SECONDS) {
        return;
    }
    if (p->vType == 1 || p->vType == 2) {
        if (p->ammo <= 0) {
            return;
        }
        p->ammo--;
    }
    float yawRad = DEG2RAD * p->yaw;
    float pitchRad = DEG2RAD * p->pitch;
    Vector3 dir = { sinf(-yawRad)*cosf(pitchRad), sinf(pitchRad), -cosf(yawRad)*cosf(pitchRad) };
    Vector3 start = v_add(p->pos, v_mul(dir, 0.8f));
    Color col = (p->vType==0? RED : BLUE);
    int vix = addVoxel(start.x, start.y, start.z, false, true, col, p->vType);
    if (vix >= 0) {
        Voxel *shot = &voxels[vix];
        Vector3 vel = v_mul(dir, 60.0f);
        shot->vel = vel;
        shot->owner = idx;
        shot->activator = idx;
        if (p->vType == 0) {
            shot->isBullet = true;
            shot->glueEligible = false;
        }
        for (int i = 0; i < 8; ++i) {
            shot->particles[i].vel = vel;
        }
    }
    p->last_shot_time = now;
}

static const char *bullet_type_symbol(int type) {
    switch (type) {
        case 0: return "O";
        case 1: return "-";
        case 2: return "+";
        default: return "?";
    }
}

static int clamp_active_players(int count) {
    if (count < 1) {
        return 1;
    }
    if (count > MAX_PLAYERS) {
        return MAX_PLAYERS;
    }
    return count;
}

static bool keyboard_player_index(int index) {
    return index == 0 || index == 1;
}

static void get_view_layout(int player_count, int *w, int *h) {
    if (player_count <= 1) {
        *w = SCREEN_WIDTH;
        *h = SCREEN_HEIGHT;
    } else if (player_count == 2) {
        *w = SCREEN_WIDTH / 2;
        *h = SCREEN_HEIGHT;
    } else {
        *w = SCREEN_WIDTH / 2;
        *h = SCREEN_HEIGHT / 2;
    }
}

static void get_viewport(int index, int player_count, int *x, int *y, int *w, int *h) {
    get_view_layout(player_count, w, h);
    if (player_count <= 1) {
        *x = 0;
        *y = 0;
    } else if (player_count == 2) {
        *x = index * (*w);
        *y = 0;
    } else {
        *x = (index % 2) * (*w);
        *y = (index / 2) * (*h);
    }
}

static void ensure_render_targets(RenderTexture2D *screens,
                                  int *current_players,
                                  int *current_w,
                                  int *current_h,
                                  int player_count)
{
    int view_w = 0;
    int view_h = 0;
    get_view_layout(player_count, &view_w, &view_h);
    if (player_count == *current_players &&
        view_w == *current_w &&
        view_h == *current_h) {
        return;
    }

    for (int i = 0; i < *current_players; ++i) {
        if (screens[i].id != 0) {
            UnloadRenderTexture(screens[i]);
        }
        screens[i] = (RenderTexture2D){ 0 };
    }

    for (int i = 0; i < player_count; ++i) {
        screens[i] = LoadRenderTexture(view_w, view_h);
    }
    *current_players = player_count;
    *current_w = view_w;
    *current_h = view_h;
}

static bool spawn_position_clear(Vector3 pos, float min_dist) {
    float half = PLAYER_SIZE * 0.5f;
    int minx = (int)floorf((pos.x - half) / VOXEL_SIZE);
    int maxx = (int)floorf((pos.x + half) / VOXEL_SIZE);
    int miny = (int)floorf((pos.y - half) / VOXEL_SIZE);
    int maxy = (int)floorf((pos.y + half) / VOXEL_SIZE);
    int minz = (int)floorf((pos.z - half) / VOXEL_SIZE);
    int maxz = (int)floorf((pos.z + half) / VOXEL_SIZE);

    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            for (int z = minz; z <= maxz; ++z) {
                if (occupied(x, y, z)) {
                    return false;
                }
            }
        }
    }

    if (min_dist > 0.0f) {
        float min_dist_sq = min_dist * min_dist;
        for (int i = 0; i < activePlayers; ++i) {
            Vector3 delta = v_sub(pos, players[i].pos);
            if (v_dot(delta, delta) < min_dist_sq) {
                return false;
            }
        }
    }
    return true;
}

static Vector3 pick_player_spawn(int player_index) {
    if (!randomSpawnEnabled) {
        return playerSpawnPositions[0];
    }

    float min_distance = PLAYER_SIZE * 2.5f;
    float min_pos = -FLOOR_SIZE + PLAYER_RADIUS;
    float max_pos = FLOOR_SIZE - PLAYER_RADIUS;
    for (int attempt = 0; attempt < 64; ++attempt) {
        Vector3 pos = {
            randomInRange(min_pos, max_pos),
            BASE_EYE_HEIGHT,
            randomInRange(min_pos, max_pos)
        };
        if (spawn_position_clear(pos, min_distance)) {
            return pos;
        }
    }
    if (player_index >= 0 && player_index < MAX_PLAYERS) {
        return playerSpawnPositions[player_index];
    }
    return (Vector3){ 0.0f, BASE_EYE_HEIGHT, 0.0f };
}

static int brush_extent_for_voxel(const Voxel *v) {
    int base = (voxelBrushSpan < 1) ? 1 : voxelBrushSpan;
    if (!v || (v->type != 1 && v->type != 2)) {
        return base;
    }
    int owner = v->owner;
    if (owner < 0 || owner >= activePlayers) {
        return base;
    }
    float ratio = fmaxf(players[owner].kd_ratio, 0.1f);
    float scaled = (float)base / ratio;
    int extent = (int)roundf(scaled);
    if (extent < 1) {
        extent = 1;
    } else if (extent > 5) {
        extent = 5;
    }
    return extent;
}

static int player_max_health(const Player *p) {
    float ratio = fmaxf(p ? p->kd_ratio : 1.0f, 0.1f);
    int value = (int)roundf((float)BASE_HEALTH / ratio);
    return (value < 1) ? 1 : value;
}

static int player_max_shield(const Player *p) {
    float ratio = fmaxf(p ? p->kd_ratio : 1.0f, 0.1f);
    int value = (int)roundf((float)BASE_SHIELD / ratio);
    return (value < 1) ? 1 : value;
}

static void draw_hud_bars(const Player *p, int viewport_w, int viewport_h) {
    if (!p || viewport_w <= 0 || viewport_h <= 0) {
        return;
    }

    int max_w = viewport_w - HUD_PADDING_X * 2;
    int bar_w = HUD_BAR_WIDTH;
    float health_max_ratio = (float)player_max_health(p) / (float)BASE_HEALTH;
    float shield_max_ratio = (float)player_max_shield(p) / (float)BASE_SHIELD;
    float max_ratio = fmaxf(health_max_ratio, shield_max_ratio);
    if (max_ratio < 0.1f) {
        max_ratio = 0.1f;
    }
    if ((float)bar_w * max_ratio > (float)max_w) {
        bar_w = (int)floorf((float)max_w / max_ratio);
    }
    if (bar_w < 1) {
        bar_w = 1;
    }

    int bullet_gap = HUD_BAR_SPACING;
    int total_h = HUD_BAR_THICKNESS * 2 + HUD_BAR_SPACING + HUD_BULLET_TEXT_SIZE + bullet_gap;
    int y = viewport_h - HUD_PADDING_Y - total_h;

    int health_max = player_max_health(p);
    int shield_max = player_max_shield(p);
    float health_ratio = clampf((float)p->health / (float)health_max, 0.0f, 1.0f);
    float shield_ratio = clampf((float)p->shield / (float)shield_max, 0.0f, 1.0f);
    int health_bar_w = (int)roundf((float)bar_w * health_max_ratio);
    int shield_bar_w = (int)roundf((float)bar_w * shield_max_ratio);
    int max_bar_w = (health_bar_w > shield_bar_w) ? health_bar_w : shield_bar_w;
    int x_max = viewport_w - HUD_PADDING_X - max_bar_w;

    Color bar_bg = (Color){ 25, 25, 25, 200 };
    Color health_color = (Color){ 210, 70, 70, 230 };
    Color shield_color = (Color){ 70, 140, 220, 230 };

    Color bullet_color = BLACK;
    const char *bullet_symbol = bullet_type_symbol(p->vType);
    if (p->vType == 1) {
        bullet_color = RED;
    } else if (p->vType == 2) {
        bullet_color = BLUE;
    }
    int bullet_w = MeasureText(bullet_symbol, HUD_BULLET_TEXT_SIZE);
    int bullet_x = x_max + (max_bar_w - bullet_w) / 2;
    DrawText(bullet_symbol, bullet_x, y, HUD_BULLET_TEXT_SIZE, bullet_color);

    int bar_y = y + HUD_BULLET_TEXT_SIZE + bullet_gap;
    int health_x = viewport_w - HUD_PADDING_X - health_bar_w;
    DrawRectangle(health_x, bar_y, health_bar_w, HUD_BAR_THICKNESS, bar_bg);
    DrawRectangle(health_x, bar_y, (int)roundf(health_bar_w * health_ratio), HUD_BAR_THICKNESS, health_color);
    DrawRectangleLines(health_x, bar_y, health_bar_w, HUD_BAR_THICKNESS, BLACK);

    int shield_y = bar_y + HUD_BAR_THICKNESS + HUD_BAR_SPACING;
    int shield_x = viewport_w - HUD_PADDING_X - shield_bar_w;
    DrawRectangle(shield_x, shield_y, shield_bar_w, HUD_BAR_THICKNESS, bar_bg);
    DrawRectangle(shield_x, shield_y, (int)roundf(shield_bar_w * shield_ratio), HUD_BAR_THICKNESS, shield_color);
    DrawRectangleLines(shield_x, shield_y, shield_bar_w, HUD_BAR_THICKNESS, BLACK);

    int health_display = p->health;
    if (health_display > health_max) {
        health_display = health_max;
    }
    int shield_display = p->shield;
    if (shield_display > shield_max) {
        shield_display = shield_max;
    }
    const char *health_text = TextFormat("HP %d/%d", health_display, health_max);
    const char *shield_text = TextFormat("SH %d/%d", shield_display, shield_max);
    int health_text_w = MeasureText(health_text, HUD_BAR_TEXT_SIZE);
    int shield_text_w = MeasureText(shield_text, HUD_BAR_TEXT_SIZE);
    int health_text_x = health_x + (health_bar_w - health_text_w) / 2;
    int shield_text_x = shield_x + (shield_bar_w - shield_text_w) / 2;
    int text_y = bar_y + (HUD_BAR_THICKNESS - HUD_BAR_TEXT_SIZE) / 2;
    int shield_text_y = shield_y + (HUD_BAR_THICKNESS - HUD_BAR_TEXT_SIZE) / 2;
    DrawText(health_text, health_text_x, text_y, HUD_BAR_TEXT_SIZE, WHITE);
    DrawText(shield_text, shield_text_x, shield_text_y, HUD_BAR_TEXT_SIZE, WHITE);

    int ammo_size = HUD_AMMO_SIZE * 3;
    int ammo_gap = HUD_AMMO_GAP * 3;
    int ammo_w = AMMO_MAX * ammo_size + (AMMO_MAX - 1) * ammo_gap;
    int ammo_x = HUD_PADDING_X;
    int ammo_y = viewport_h - HUD_PADDING_Y - ammo_size;
    Color ammo_full = (Color){ 230, 230, 230, 220 };
    Color ammo_empty = (Color){ 60, 60, 60, 200 };
    for (int i = 0; i < AMMO_MAX; ++i) {
        int x = ammo_x + i * (ammo_size + ammo_gap);
        Color fill = (i < p->ammo) ? ammo_full : ammo_empty;
        DrawRectangle(x, ammo_y, ammo_size, ammo_size, fill);
        DrawRectangleLines(x, ammo_y, ammo_size, ammo_size, BLACK);
    }
}
// Append the 12 edges (24 vertices) of a cube to the current RL_LINES batch
static void drawCubeEdges(const Voxel *voxel)
{
    Vector3 v[8];
    for (int i = 0; i < 8; ++i) v[i] = voxel->particles[i].pos;

    static const int edge_indices[12][2] = {
        {0,1},{1,3},{3,2},{2,0},
        {4,5},{5,7},{7,6},{6,4},
        {0,4},{1,5},{3,7},{2,6}
    };

    for (int e = 0; e < 12; ++e) {
        const Vector3 a = v[edge_indices[e][0]];
        const Vector3 b = v[edge_indices[e][1]];
        rlVertex3f(a.x, a.y, a.z);
        rlVertex3f(b.x, b.y, b.z);
    }
}

static void compute_dynamic_face_visibility(int idx, bool faces[6])
{
    if (!faces || idx < 0 || idx >= voxel_count) {
        return;
    }
    Voxel *v = &voxels[idx];
    VoxelWorldBounds bounds;
    voxel_particle_world_bounds(v, &bounds);
    int minx = (int)floorf((bounds.minx + GRID_EPSILON) / VOXEL_SIZE);
    int maxx = (int)floorf((bounds.maxx - GRID_EPSILON) / VOXEL_SIZE);
    int miny = (int)floorf((bounds.miny + GRID_EPSILON) / VOXEL_SIZE);
    int maxy = (int)floorf((bounds.maxy - GRID_EPSILON) / VOXEL_SIZE);
    int minz = (int)floorf((bounds.minz + GRID_EPSILON) / VOXEL_SIZE);
    int maxz = (int)floorf((bounds.maxz - GRID_EPSILON) / VOXEL_SIZE);

    for (int f = 0; f < 6; ++f) {
        faces[f] = true;
    }

    for (int y = miny; y <= maxy && faces[0]; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            int nidx = table_get(maxx + 1, y, z);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 0)) {
                faces[0] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 0);
                break;
            }
        }
    }

    for (int y = miny; y <= maxy && faces[1]; ++y) {
        for (int z = minz; z <= maxz; ++z) {
            int nidx = table_get(minx - 1, y, z);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 1)) {
                faces[1] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 1);
                break;
            }
        }
    }

    for (int x = minx; x <= maxx && faces[2]; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            int nidx = table_get(x, maxy + 1, z);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 2)) {
                faces[2] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 2);
                break;
            }
        }
    }

    for (int x = minx; x <= maxx && faces[3]; ++x) {
        for (int z = minz; z <= maxz; ++z) {
            int nidx = table_get(x, miny - 1, z);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 3)) {
                faces[3] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 3);
                break;
            }
        }
    }

    for (int x = minx; x <= maxx && faces[4]; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            int nidx = table_get(x, y, maxz + 1);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 4)) {
                faces[4] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 4);
                break;
            }
        }
    }

    for (int x = minx; x <= maxx && faces[5]; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            int nidx = table_get(x, y, minz - 1);
            if (nidx < 0 || nidx >= voxel_count || nidx == idx) {
                continue;
            }
            Voxel *neighbor = &voxels[nidx];
            if (face_blocked_by_voxel(v, neighbor, 5)) {
                faces[5] = false;
                log_span2_face_cull(v, idx, neighbor, nidx, 5);
                break;
            }
        }
    }
}

static void compute_voxel_face_visibility(int idx, bool faces[6])
{
    if (!faces || idx < 0 || idx >= voxel_count) {
        return;
    }
    if (voxels[idx].simulate) {
        if (renderAllDynamicFaces) {
            for (int f = 0; f < 6; ++f) {
                faces[f] = true;
            }
            return;
        }
        compute_dynamic_face_visibility(idx, faces);
        return;
    }
    mark_surface(idx);
    memcpy(faces, voxels[idx].surface, sizeof(voxels[idx].surface));
}

static void drawCubeMan(const Voxel *voxel, const bool faces[6])
{
    Vector3 v[8];
    for (int i = 0; i < 8; ++i) v[i] = voxel->particles[i].pos;

    Color displayColor = voxel_display_color(voxel);
    rlColor4ub(displayColor.r, displayColor.g, displayColor.b, displayColor.a);

    if (!faces || faces[4]) {
        rlNormal3f(0.0f, 0.0f, 1.0f);
        rlVertex3f(v[4].x, v[4].y, v[4].z);
        rlVertex3f(v[5].x, v[5].y, v[5].z);
        rlVertex3f(v[7].x, v[7].y, v[7].z);
        rlVertex3f(v[4].x, v[4].y, v[4].z);
        rlVertex3f(v[7].x, v[7].y, v[7].z);
        rlVertex3f(v[6].x, v[6].y, v[6].z);
    }

    if (!faces || faces[5]) {
        rlNormal3f(0.0f, 0.0f, -1.0f);
        rlVertex3f(v[1].x, v[1].y, v[1].z);
        rlVertex3f(v[0].x, v[0].y, v[0].z);
        rlVertex3f(v[2].x, v[2].y, v[2].z);
        rlVertex3f(v[1].x, v[1].y, v[1].z);
        rlVertex3f(v[2].x, v[2].y, v[2].z);
        rlVertex3f(v[3].x, v[3].y, v[3].z);
    }

    if (!faces || faces[2]) {
        rlNormal3f(0.0f, 1.0f, 0.0f);
        rlVertex3f(v[6].x, v[6].y, v[6].z);
        rlVertex3f(v[7].x, v[7].y, v[7].z);
        rlVertex3f(v[3].x, v[3].y, v[3].z);
        rlVertex3f(v[6].x, v[6].y, v[6].z);
        rlVertex3f(v[3].x, v[3].y, v[3].z);
        rlVertex3f(v[2].x, v[2].y, v[2].z);
    }

    if (!faces || faces[3]) {
        rlNormal3f(0.0f, -1.0f, 0.0f);
        rlVertex3f(v[0].x, v[0].y, v[0].z);
        rlVertex3f(v[1].x, v[1].y, v[1].z);
        rlVertex3f(v[5].x, v[5].y, v[5].z);
        rlVertex3f(v[0].x, v[0].y, v[0].z);
        rlVertex3f(v[5].x, v[5].y, v[5].z);
        rlVertex3f(v[4].x, v[4].y, v[4].z);
    }

    if (!faces || faces[0]) {
        rlNormal3f(1.0f, 0.0f, 0.0f);
        rlVertex3f(v[5].x, v[5].y, v[5].z);
        rlVertex3f(v[1].x, v[1].y, v[1].z);
        rlVertex3f(v[3].x, v[3].y, v[3].z);
        rlVertex3f(v[5].x, v[5].y, v[5].z);
        rlVertex3f(v[3].x, v[3].y, v[3].z);
        rlVertex3f(v[7].x, v[7].y, v[7].z);
    }

    if (!faces || faces[1]) {
        rlNormal3f(-1.0f, 0.0f, 0.0f);
        rlVertex3f(v[0].x, v[0].y, v[0].z);
        rlVertex3f(v[4].x, v[4].y, v[4].z);
        rlVertex3f(v[6].x, v[6].y, v[6].z);
        rlVertex3f(v[0].x, v[0].y, v[0].z);
        rlVertex3f(v[6].x, v[6].y, v[6].z);
        rlVertex3f(v[2].x, v[2].y, v[2].z);
    }
}

static Color particle_velocity_color(float speed)
{
    float t = clampf(speed / PARTICLE_DEBUG_MAX_SPEED, 0.0f, 1.0f);
    float hue = 240.0f - 240.0f * t; // 240=blue (slow) to 0=red (fast)
    Color c = ColorFromHSV(hue, 0.85f, 1.0f);
    c.a = 255;
    return c;
}

static void draw_particle_debug(void)
{
    const Color baseColor = { 255, 200, 80, 255 };

    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (!voxel->simulate) continue;
        float baseRadius = voxel_particle_radius(voxel);
        float markerRadius = fmaxf(baseRadius * PARTICLE_DEBUG_MARKER_RADIUS, 0.02f);
        for (int j = 0; j < 8; ++j) {
            const Particle *p = &voxel->particles[j];
            Color markerColor = baseColor;
            if (debugColorParticlesByVelocity && p->inv_mass > 0.0f) {
                float speed = v_length(p->vel);
                markerColor = particle_velocity_color(speed);
            }
            DrawSphere(p->pos, markerRadius, markerColor);
        }
    }
}

// Return the parametric distance t (along the ray) to the first voxel boundary
// we’ll cross on this axis.  If the ray is exactly parallel to the axis (d == 0),
// we return +INF so the DDA never steps on this axis.
static inline float ray_t_to_next_plane(float p,      // ray.position.x  (or y / z)
                                        float d,      // ray.direction.x (or y / z) – assumed normalised-ish
                                        int   cell,   // floor(p / VOXEL_SIZE)
                                        float voxel)  // VOXEL_SIZE
{
    if (d == 0.0f) return INFINITY;          // will never cross a plane on this axis

    // Which side of the voxel are we travelling toward?
    // If d > 0 we leave via the "high" face at (cell+1)*voxel,
    // else we leave via the "low" face at  cell   *voxel.
    float nextPlane = (d > 0.0f) ? (cell + 1) * voxel
                                 :  cell      * voxel;

    // Distance along the ray parameter t = (planePos - pos) / dir
    return (nextPlane - p) / d;              // d has sign, so t is positive
}

typedef struct {
    int id;       // voxel index in your array
    int x, y, z;  // voxel grid coords (optional, handy for debugging)
    float t;      // parametric distance along the ray
} VoxelHit;

static int first_voxel_hit(Ray ray, float t_max, int ignore_id) {
    // 1. Normalise once
    ray.direction = v_norm(ray.direction);
    Vector3 dir   = ray.direction;

    // 2. Starting voxel indices
    int x = (int)floorf(ray.position.x / VOXEL_SIZE);
    int y = (int)floorf(ray.position.y / VOXEL_SIZE);
    int z = (int)floorf(ray.position.z / VOXEL_SIZE);

    // 3. Step direction per axis
    int stepX = (dir.x > 0.0f) ?  1 : (dir.x < 0.0f ? -1 : 0);
    int stepY = (dir.y > 0.0f) ?  1 : (dir.y < 0.0f ? -1 : 0);
    int stepZ = (dir.z > 0.0f) ?  1 : (dir.z < 0.0f ? -1 : 0);

    // 4. t delta per axis (distance to cross one voxel)
    float txDelta = (dir.x == 0.0f) ? INFINITY : fabsf(VOXEL_SIZE / dir.x);
    float tyDelta = (dir.y == 0.0f) ? INFINITY : fabsf(VOXEL_SIZE / dir.y);
    float tzDelta = (dir.z == 0.0f) ? INFINITY : fabsf(VOXEL_SIZE / dir.z);

    // 5. first plane crossings
    float txNext = ray_t_to_next_plane(ray.position.x, dir.x, x, VOXEL_SIZE);
    float tyNext = ray_t_to_next_plane(ray.position.y, dir.y, y, VOXEL_SIZE);
    float tzNext = ray_t_to_next_plane(ray.position.z, dir.z, z, VOXEL_SIZE);

    // 6. DDA walk
    while (fminf(txNext, fminf(tyNext, tzNext)) <= t_max) {
        int id = table_get(x, y, z);
        if (id >= 0 && id != ignore_id) {                // hit!
            return id;
        }

        if (txNext < tyNext && txNext < tzNext) {
            x += stepX;
            txNext += txDelta;
        } else if (tyNext < tzNext) {
            y += stepY;
            tyNext += tyDelta;
        } else {
            z += stepZ;
            tzNext += tzDelta;
        }
    }
    return -1;                        // no voxel found within t_max
}



// Generate a greedy mesh of all visible voxels ( i think the bug where single voxels are not drawn right is somewhere in here)
// Only one layer is drawn per voxel, this makes layers disappear on the individual voxel level

static void merge_rects_on_plane(int count, int *list, int plane, bool positive) {
    if (!count) return;

    // Group voxels by layer
    int layers[count], layerCount = 0;
    for (int i = 0; i < count; i++) {
        int idx = list[i];
        int layer = 0;
        switch (plane) {
            case 0: layer = voxels[idx].gz; break; // XY-plane, group by z
            case 1: layer = voxels[idx].gy; break; // XZ-plane, group by y
            case 2: layer = voxels[idx].gx; break; // YZ-plane, group by x
        }
        bool seen = false;
        for (int j = 0; j < layerCount; j++) {
            if (layers[j] == layer) {
                seen = true;
                break;
            }
        }
        if (!seen) layers[layerCount++] = layer;
    }

    // For each layer, find and merge rectangles
    for (int l = 0; l < layerCount; l++) {
        int layer = layers[l];

        // Find the bounds of the current layer
        int minI = INT_MAX, maxI = INT_MIN, minJ = INT_MAX, maxJ = INT_MIN;
        for (int i = 0; i < count; i++) {
            Voxel *v = &voxels[list[i]];
            int vLayer = 0, vI = 0, vJ = 0;
            switch (plane) {
                case 0: vLayer = v->gz; vI = v->gx; vJ = v->gy; break;
                case 1: vLayer = v->gy; vI = v->gx; vJ = v->gz; break;
                case 2: vLayer = v->gx; vI = v->gy; vJ = v->gz; break;
            }

            if (vLayer != layer) continue;
            if (vI < minI) minI = vI;
            if (vI > maxI) maxI = vI;
            if (vJ < minJ) minJ = vJ;
            if (vJ > maxJ) maxJ = vJ;
        }

        int w = maxI - minI + 1;
        int h = maxJ - minJ + 1;
        bool mask[w][h];
        memset(mask, 0, sizeof(mask));

        for (int i = 0; i < count; i++) {
            Voxel *v = &voxels[list[i]];
            int vLayer = 0, vI = 0, vJ = 0;
            switch (plane) {
                case 0: vLayer = v->gz; vI = v->gx; vJ = v->gy; break;
                case 1: vLayer = v->gy; vI = v->gx; vJ = v->gz; break;
                case 2: vLayer = v->gx; vI = v->gy; vJ = v->gz; break;
            }
            if (vLayer != layer) continue;
            mask[vI - minI][vJ - minJ] = true;
        }

        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                if (!mask[i][j]) continue;

                int si = i, sj = j;
                int ww = 1;
                while (si + ww < w && mask[si + ww][sj]) ww++;

                int hh = 1;
                for (;;) {
                    bool block = false;
                    for (int k = 0; k < ww; k++) {
                        if (sj + hh >= h || !mask[si + k][sj + hh]) {
                            block = true;
                            break;
                        }
                    }
                    if (block) break;
                    hh++;
                }
                
                for (int dj = 0; dj < hh; dj++) {
                    for (int di = 0; di < ww; di++) {
                        mask[si + di][sj + dj] = false;
                    }
                }
                
                i = si + ww - 1;
                
                Patch *pt = &patches[patchCount++];
                pt->plane = plane;
                pt->layer = layer;
                pt->i0 = minI + si;
                pt->j0 = minJ + sj;
                pt->di = ww;
                pt->dj = hh;
                pt->positive = positive;
                
                int baseIdx = -1;
                switch (plane) {
                    case 0: baseIdx = table_get(pt->i0, pt->j0, layer); break;
                    case 1: baseIdx = table_get(pt->i0, layer, pt->j0); break;
                    case 2: baseIdx = table_get(layer, pt->i0, pt->j0); break;
                }
                pt->voxelIndex = baseIdx;
                Color baseColor = { 160, 160, 160, 255 };
                if (baseIdx >= 0 && baseIdx < voxel_count) {
                    baseColor = voxel_display_color(&voxels[baseIdx]);
                }
                pt->col = baseColor;
            }
        }
    }
}


static int yzPosList[MAX_VOXELS], yzNegList[MAX_VOXELS];
static int xzPosList[MAX_VOXELS], xzNegList[MAX_VOXELS];
static int xyPosList[MAX_VOXELS], xyNegList[MAX_VOXELS];

static Mesh gen_greedy_mesh(void) {
    Mesh mesh = { 0 };
    patchCount = 0;

    // Collect all static (non-simulated) voxels by visible faces (principal planes)
    int yzPosCount = 0, yzNegCount = 0; // +X, -X
    int xzPosCount = 0, xzNegCount = 0; // +Y, -Y
    int xyPosCount = 0, xyNegCount = 0; // +Z, -Z
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (v->simulate) continue;
        if (v->surface[0]) yzPosList[yzPosCount++] = i; // +X face
        if (v->surface[1]) yzNegList[yzNegCount++] = i; // -X face
        if (v->surface[2]) xzPosList[xzPosCount++] = i; // +Y face
        if (v->surface[3]) xzNegList[xzNegCount++] = i; // -Y face
        if (v->surface[4]) xyPosList[xyPosCount++] = i; // +Z face
        if (v->surface[5]) xyNegList[xyNegCount++] = i; // -Z face
    }

    // Merge rectangles in each principal plane
    merge_rects_on_plane(xyPosCount, xyPosList, 0, true); // XY-plane, +Z
    merge_rects_on_plane(xyNegCount, xyNegList, 0, false); // XY-plane, -Z
    merge_rects_on_plane(xzPosCount, xzPosList, 1, true); // XZ-plane, +Y
    merge_rects_on_plane(xzNegCount, xzNegList, 1, false); // XZ-plane, -Y
    merge_rects_on_plane(yzPosCount, yzPosList, 2, true); // YZ-plane, +X
    merge_rects_on_plane(yzNegCount, yzNegList, 2, false); // YZ-plane, -X

    // Build Raylib mesh from collected patches
    int vCount = patchCount * 4;
    int iCount = patchCount * 6;
    float *verts  = calloc(vCount*3, sizeof(float));
    float *norms  = calloc(vCount*3, sizeof(float));
    float *uvs    = calloc(vCount*2, sizeof(float));
    unsigned char *cols = calloc(vCount*4, sizeof(unsigned char));
    unsigned short *inds = calloc(iCount, sizeof(unsigned short));

    for (int p = 0; p < patchCount; p++) {
        Patch *pt = &patches[p];
        // Determine plane basis vectors and position
        Vector3 origin = {0};
        Vector3 iu, ju, norm;
        switch (pt->plane) {
        case 0: // XY plane at Z=layer
            origin = (Vector3){ (pt->i0 + 0.0f)*VOXEL_SIZE,
                                (pt->j0 + 0.0f)*VOXEL_SIZE,
                                (pt->layer + (pt->positive?1:0))*VOXEL_SIZE };
            iu = (Vector3){ VOXEL_SIZE*pt->di, 0, 0 };
            ju = (Vector3){ 0, VOXEL_SIZE*pt->dj, 0 };
            norm = pt->positive? (Vector3){0,0,1} : (Vector3){0,0,-1};
            break;
        case 1: // XZ plane at Y=layer
            origin = (Vector3){ (pt->i0 + 0.0f)*VOXEL_SIZE,
                                (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                (pt->j0 + 0.0f)*VOXEL_SIZE };
            iu = (Vector3){ VOXEL_SIZE*pt->di, 0, 0 };
            ju = (Vector3){ 0, 0, VOXEL_SIZE*pt->dj };
            norm = pt->positive? (Vector3){0,1,0} : (Vector3){0,-1,0};
            break;
        default: // YZ plane at X=layer
            origin = (Vector3){ (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                (pt->i0 + 0.0f)*VOXEL_SIZE,
                                (pt->j0 + 0.0f)*VOXEL_SIZE };
            iu = (Vector3){ 0, VOXEL_SIZE*pt->di, 0 };
            ju = (Vector3){ 0, 0, VOXEL_SIZE*pt->dj };
            norm = pt->positive? (Vector3){1,0,0} : (Vector3){-1,0,0};
            break;
        }
        // build quad vertices
        int vofs = p*4;
        Vector3 corners[4] = {
            origin,
            v_add(origin, iu),
            v_add(v_add(origin, iu), ju),
            v_add(origin, ju)
        };
        for (int k = 0; k < 4; k++) {
            verts[(vofs+k)*3 + 0] = corners[k].x;
            verts[(vofs+k)*3 + 1] = corners[k].y;
            verts[(vofs+k)*3 + 2] = corners[k].z;
            norms[(vofs+k)*3 + 0] = norm.x;
            norms[(vofs+k)*3 + 1] = norm.y;
            norms[(vofs+k)*3 + 2] = norm.z;
            uvs[(vofs+k)*2 + 0] = (k==1||k==2);
            uvs[(vofs+k)*2 + 1] = (k>=2);
            cols[(vofs+k)*4 + 0] = pt->col.r;
            cols[(vofs+k)*4 + 1] = pt->col.g;
            cols[(vofs+k)*4 + 2] = pt->col.b;
            cols[(vofs+k)*4 + 3] = pt->col.a;
        }
        int iofs = p*6;
        inds[iofs + 0] = vofs + 0;
        inds[iofs + 1] = vofs + 1;
        inds[iofs + 2] = vofs + 2;
        inds[iofs + 3] = vofs + 0;
        inds[iofs + 4] = vofs + 2;
        inds[iofs + 5] = vofs + 3;
    }
    mesh.vertexCount   = vCount;
    mesh.triangleCount = patchCount * 2;
    mesh.vertices      = verts;
    mesh.normals       = norms;
    mesh.texcoords     = uvs;
    mesh.colors        = cols;
    mesh.indices       = inds;
    UploadMesh(&mesh, false);

    return mesh;
}

static void update_static_patch_colors(void)
{
    if (!greedyMesh.vertices || !greedyMesh.colors || patchCount <= 0) {
        return;
    }

    bool anyChange = false;
    unsigned char *colorBuffer = (unsigned char *)greedyMesh.colors;
    for (int p = 0; p < patchCount; ++p) {
        Patch *pt = &patches[p];
        Color newColor = pt->col;
        if (pt->voxelIndex >= 0 && pt->voxelIndex < voxel_count) {
            newColor = voxel_display_color(&voxels[pt->voxelIndex]);
        }
        if (newColor.r != pt->col.r || newColor.g != pt->col.g ||
            newColor.b != pt->col.b || newColor.a != pt->col.a) {
            pt->col = newColor;
            anyChange = true;
        }

        int vofs = p * 4;
        for (int k = 0; k < 4; ++k) {
            unsigned char *dst = &colorBuffer[(vofs + k) * 4];
            dst[0] = pt->col.r;
            dst[1] = pt->col.g;
            dst[2] = pt->col.b;
            dst[3] = pt->col.a;
        }
    }

    if (anyChange) {
        UpdateMeshBuffer(greedyMesh, 3, greedyMesh.colors,
                         greedyMesh.vertexCount * 4 * sizeof(unsigned char), 0);
    }
}

// Draw all voxels via greedy mesh instead of per-voxel raycasting
static void DrawVoxels(Camera3D cam) {
    (void)cam;

    if (debugLogSpan2Faces) {
        debugSpan2FaceLogBudget = DEBUG_SPAN2_FACE_LOG_INIT;
    }

    if (meshDirty) {
        if (greedyMesh.vertices) {
            UnloadMesh(greedyMesh);
            greedyMesh = (Mesh){ 0 };
        }
        greedyMesh = gen_greedy_mesh();
        meshDirty = false;
    }
    if (!greedyMaterialInit) {
        greedyMaterial = LoadMaterialDefault();
        greedyMaterialInit = true;
    }

    if (greedyMesh.vertices) {
        update_static_patch_colors();
    }

    rlDisableBackfaceCulling();
    if (greedyMesh.vertices) {
        DrawMesh(greedyMesh, greedyMaterial, MatrixIdentity());
    }
    rlBegin(RL_LINES);
    rlColor4ub(0, 0, 0, 60);
    for (int p = 0; p < patchCount; p++) {
        Patch *pt = &patches[p];
        Vector3 origin, iu, ju;
        switch (pt->plane) {
            case 0:
                origin = (Vector3){ pt->i0*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE, 0, 0 };
                ju = (Vector3){ 0, VOXEL_SIZE, 0 };
                break;
            case 1:
                origin = (Vector3){ pt->i0*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE, 0, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE };
                break;
            default:
                origin = (Vector3){ (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    pt->i0*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE };
                iu = (Vector3){ 0, VOXEL_SIZE, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE };
                break;
        }
        for (int ix = 0; ix <= pt->di; ++ix) {
            Vector3 a = v_add(origin, v_mul(iu, ix));
            Vector3 b = v_add(a, v_mul(ju, pt->dj));
            rlVertex3f(a.x, a.y, a.z); rlVertex3f(b.x, b.y, b.z);
        }
        for (int iy = 0; iy <= pt->dj; ++iy) {
            Vector3 a = v_add(origin, v_mul(ju, iy));
            Vector3 b = v_add(a, v_mul(iu, pt->di));
            rlVertex3f(a.x, a.y, a.z); rlVertex3f(b.x, b.y, b.z);
        }
    }
    rlEnd();

    rlBegin(RL_TRIANGLES);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!v->simulate) {
            continue;
        }
        bool faces[6];
        compute_voxel_face_visibility(i, faces);
        drawCubeMan(v, faces);
    }
    rlEnd();

    rlBegin(RL_LINES);
    rlColor4ub(0, 0, 0, 255);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!v->simulate) {
            continue;
        }
        drawCubeEdges(v);
    }
    rlEnd();
    rlEnableBackfaceCulling();

    if (debugDrawParticles) {
        draw_particle_debug();
    }
}

static Color scale_color(Color c, float scale, unsigned char alpha) {
    Color out = {
        (unsigned char)clampf(c.r * scale, 0.0f, 255.0f),
        (unsigned char)clampf(c.g * scale, 0.0f, 255.0f),
        (unsigned char)clampf(c.b * scale, 0.0f, 255.0f),
        alpha
    };
    return out;
}

static Color player_palette_color(int index) {
    static const Color palette[MAX_PLAYERS] = {
        { 230, 80, 80, 255 },
        { 80, 200, 120, 255 },
        { 80, 140, 230, 255 },
        { 230, 180, 70, 255 }
    };
    if (index < 0 || index >= MAX_PLAYERS) {
        return (Color){ 200, 200, 200, 255 };
    }
    return palette[index];
}

static void draw_players(void) {
    for (int i = 0; i < activePlayers; i++) {
        Player *p = &players[i];
        Color base = player_palette_color(i);
        Color base_dark = scale_color(base, 0.35f, 255);
        DrawCube(p->pos, PLAYER_SIZE,PLAYER_SIZE, PLAYER_SIZE, base_dark);
        DrawCubeWires(p->pos, PLAYER_SIZE,PLAYER_SIZE,PLAYER_SIZE, base_dark);
        if (p->shield > 0) {
            int max_shield = player_max_shield(p);
            float shield_percentage = (float)p->shield / (float)max_shield;
            unsigned char alpha = (unsigned char)clampf(90.0f + 90.0f * shield_percentage, 80.0f, 180.0f);
            Color shield_color = scale_color(base, 1.0f, alpha);
            DrawCube(p->pos, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, shield_color);
            DrawCubeWires(p->pos, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, Fade(shield_color, 0.8f));
        }
    }
}

static void HandleKeyboardInput(int i, float dt) {
    if (!keyboard_player_index(i)) {
        return;
    }
    Player *p = &players[i];
    // turn
    float yaw_accel = 0.0f;
    if ((i==0 && IsKeyDown(KEY_F)) || (i==1 && IsKeyDown(KEY_LEFT)))  yaw_accel += TURN_ACCELERATION * fmaxf(p->kd_ratio,1);
    if ((i==0 && IsKeyDown(KEY_H)) || (i==1 && IsKeyDown(KEY_RIGHT))) yaw_accel -= TURN_ACCELERATION * fmaxf(p->kd_ratio,1);

    if (yaw_accel != 0.0f) {
        p->yaw_vel += yaw_accel * dt;
    } else {
        // friction
        if (p->yaw_vel > 0) {
            p->yaw_vel -= TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->yaw_vel < 0) p->yaw_vel = 0;
        } else if (p->yaw_vel < 0) {
            p->yaw_vel += TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->yaw_vel > 0) p->yaw_vel = 0;
        }
    }
    p->yaw_vel = clampf(p->yaw_vel, -TURN_SPEED * fmaxf(p->kd_ratio,1), TURN_SPEED * fmaxf(p->kd_ratio,1));
    p->yaw += p->yaw_vel * dt;

    // look up/down
    float pitch_accel = 0.0f;
    if ((i==0 && IsKeyDown(KEY_T)) || (i==1 && IsKeyDown(KEY_UP)))   pitch_accel += TURN_ACCELERATION * fmaxf(p->kd_ratio,1);
    if ((i==0 && IsKeyDown(KEY_G)) || (i==1 && IsKeyDown(KEY_DOWN))) pitch_accel -= TURN_ACCELERATION * fmaxf(p->kd_ratio,1);

    if (pitch_accel != 0.0f) {
        p->pitch_vel += pitch_accel * dt;
    } else {
        // friction
        if (p->pitch_vel > 0) {
            p->pitch_vel -= TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->pitch_vel < 0) p->pitch_vel = 0;
        } else if (p->pitch_vel < 0) {
            p->pitch_vel += TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->pitch_vel > 0) p->pitch_vel = 0;
        }
    }
    p->pitch_vel = clampf(p->pitch_vel, -TURN_SPEED * fmaxf(p->kd_ratio,1), TURN_SPEED * fmaxf(p->kd_ratio,1));
    p->pitch += p->pitch_vel * dt;
    p->pitch = clampf(p->pitch, -89, 89);
    // compute forward/right
    float yr = DEG2RAD * p->yaw;
    Vector3 forward = { sinf(-yr), 0, -cosf(yr) };
    Vector3 right   = { -forward.z, 0, forward.x };
    // acceleration
    Vector3 accel = {0,0,0};
    if ((i==0 && IsKeyDown(KEY_W)) || (i==1 && IsKeyDown(KEY_I))) accel = v_add(accel, forward);
    if ((i==0 && IsKeyDown(KEY_S)) || (i==1 && IsKeyDown(KEY_K))) accel = v_add(accel, v_mul(forward, -1));
    if ((i==0 && IsKeyDown(KEY_A)) || (i==1 && IsKeyDown(KEY_J))) accel = v_add(accel, v_mul(right, -1));
    if ((i==0 && IsKeyDown(KEY_D)) || (i==1 && IsKeyDown(KEY_L))) accel = v_add(accel, right);
    if (accel.x!=0 || accel.z!=0) {
        float len = sqrtf(accel.x*accel.x + accel.z*accel.z);
        accel = v_mul(accel, 1/len);
        p->vel = v_add(p->vel, v_mul(accel, ACCELERATION * fmaxf(p->kd_ratio,1) * dt));
    } else {
        // friction
        float sp = sqrtf(p->vel.x*p->vel.x + p->vel.z*p->vel.z);
        if (sp > 0) {
            float dec = FRICTION * fmaxf(p->kd_ratio,1) * dt;
            float ns = sp - dec; if (ns < 0) ns = 0;
            p->vel.x *= ns/sp;
            p->vel.z *= ns/sp;
        }
    }
}

static void HandleGamepadInput(int i, float dt) {
    if (!IsGamepadAvailable(i)) return;

    Player *p = &players[i];

    // turn (right stick horizontal)
    float yaw_accel = 0.0f;
    float yaw_axis = GetGamepadAxisMovement(i, GAMEPAD_AXIS_RIGHT_X);
    if (fabsf(yaw_axis) > 0.1f) {
        yaw_accel = -yaw_axis * TURN_ACCELERATION * fmaxf(p->kd_ratio,1);
    }

    if (yaw_accel != 0.0f) {
        p->yaw_vel += yaw_accel * dt;
    } else {
        // friction
        if (p->yaw_vel > 0) {
            p->yaw_vel -= TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->yaw_vel < 0) p->yaw_vel = 0;
        } else if (p->yaw_vel < 0) {
            p->yaw_vel += TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->yaw_vel > 0) p->yaw_vel = 0;
        }
    }
    p->yaw_vel = clampf(p->yaw_vel, -TURN_SPEED * fmaxf(p->kd_ratio,1), TURN_SPEED * fmaxf(p->kd_ratio,1));
    p->yaw += p->yaw_vel * dt;

    // look up/down (right stick vertical)
    float pitch_accel = 0.0f;
    float pitch_axis = GetGamepadAxisMovement(i, GAMEPAD_AXIS_RIGHT_Y);
    if (fabsf(pitch_axis) > 0.1f) {
        pitch_accel = -pitch_axis * TURN_ACCELERATION * fmaxf(p->kd_ratio,1);
    }

    if (pitch_accel != 0.0f) {
        p->pitch_vel += pitch_accel * dt;
    } else {
        // friction
        if (p->pitch_vel > 0) {
            p->pitch_vel -= TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->pitch_vel < 0) p->pitch_vel = 0;
        } else if (p->pitch_vel < 0) {
            p->pitch_vel += TURN_FRICTION * fmaxf(p->kd_ratio,1) * dt;
            if (p->pitch_vel > 0) p->pitch_vel = 0;
        }
    }
    p->pitch_vel = clampf(p->pitch_vel, -TURN_SPEED * fmaxf(p->kd_ratio,1), TURN_SPEED * fmaxf(p->kd_ratio,1));
    p->pitch += p->pitch_vel * dt;
    p->pitch = clampf(p->pitch, -89, 89);

    // compute forward/right
    float yr = DEG2RAD * p->yaw;
    Vector3 forward = { sinf(-yr), 0, -cosf(yr) };
    Vector3 right   = { -forward.z, 0, forward.x };

    // acceleration (left stick)
    Vector3 accel = {0,0,0};
    float accel_x = GetGamepadAxisMovement(i, GAMEPAD_AXIS_LEFT_X);
    float accel_y = GetGamepadAxisMovement(i, GAMEPAD_AXIS_LEFT_Y);

    if (fabsf(accel_y) > 0.1f) {
        accel = v_add(accel, v_mul(forward, -accel_y));
    }
    if (fabsf(accel_x) > 0.1f) {
        accel = v_add(accel, v_mul(right, accel_x));
    }

    if (accel.x!=0 || accel.z!=0) {
        float len = sqrtf(accel.x*accel.x + accel.z*accel.z);
        accel = v_mul(accel, 1/len);
        p->vel = v_add(p->vel, v_mul(accel, ACCELERATION * fmaxf(p->kd_ratio,1) * dt));
    } else {
        // friction
        float sp = sqrtf(p->vel.x*p->vel.x + p->vel.z*p->vel.z);
        if (sp > 0) {
            float dec = FRICTION * fmaxf(p->kd_ratio,1) * dt;
            float ns = sp - dec; if (ns < 0) ns = 0;
            p->vel.x *= ns/sp;
            p->vel.z *= ns/sp;
        }
    }
}

static void HandleKeyboardInput(int i, float dt);
static void HandleGamepadInput(int i, float dt);

int main(void) {
    int countFrame = 0;
    InitDebugLogging();
    SetTraceLogLevel(LOG_DEBUG);
    // init window and render textures
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Split-Screen FPS (raylib)");
    SetTargetFPS(60);
    // seed RNG
    srand((unsigned)time(NULL));
    // reset game state
    ResetGame();
    // prepare split-screen render textures
    RenderTexture2D screens[MAX_PLAYERS] = { 0 };
    int renderPlayers = 0;
    int renderW = 0;
    int renderH = 0;
    // main loop
    while (!WindowShouldClose()) {
        switch (gameState) {
            case GAME_STATE_MENU:
                // Draw menu
                BeginDrawing();
                    ClearBackground(RAYWHITE);
                    DrawText("Main Menu", SCREEN_WIDTH / 2 - MeasureText("Main Menu", 40) / 2, 100, 40, BLACK);
                    DrawText("Press ENTER to Start", SCREEN_WIDTH / 2 - MeasureText("Press ENTER to Start", 20) / 2, 200, 20, DARKGRAY);
                    DrawText("Press S for Settings", SCREEN_WIDTH / 2 - MeasureText("Press S for Settings", 20) / 2, 250, 20, DARKGRAY);
                EndDrawing();

                if (IsKeyPressed(KEY_ENTER)) {
                    gameState = GAME_STATE_PLAYING;
                }
                if (IsKeyPressed(KEY_S)) {
                    gameState = GAME_STATE_SETTINGS;
                }
                break;
            case GAME_STATE_PLAYING:
        float dt = GetFrameTime();
        update_player_ammo(dt);
        // input: shooting, bullet type, jump
        for (int i = 0; i < activePlayers; ++i) {
            if (playerInput[i] == INPUT_TYPE_GAMEPAD) {
                if (IsGamepadButtonPressed(i, GAMEPAD_BUTTON_RIGHT_TRIGGER_2)) {
                    FireVoxel(i);
                }
                if (IsGamepadButtonPressed(i, GAMEPAD_BUTTON_LEFT_TRIGGER_1)) {
                    players[i].vType = (players[i].vType + 1) % 3;
                }
                if (IsGamepadButtonPressed(i, GAMEPAD_BUTTON_RIGHT_FACE_DOWN) && players[i].onGround) {
                    players[i].vel.y = JUMP_SPEED;
                    players[i].onGround = false;
                }
            }
        }

        if (playerInput[0] == INPUT_TYPE_KEYBOARD) {
            if (IsKeyPressed(KEY_LEFT_CONTROL))  FireVoxel(0);
            if (IsKeyPressed(KEY_Q)) players[0].vType = (players[0].vType + 1) % 3;
            if (IsKeyPressed(KEY_SPACE) && players[0].onGround) {
                players[0].vel.y = JUMP_SPEED;
                players[0].onGround = false;
            }
        }
        if (playerInput[1] == INPUT_TYPE_KEYBOARD) {
            if (IsKeyPressed(KEY_RIGHT_CONTROL)) FireVoxel(1);
            if (IsKeyPressed(KEY_U)) players[1].vType = (players[1].vType + 1) % 3;
            if (IsKeyPressed(KEY_RIGHT_SHIFT) && players[1].onGround) {
                players[1].vel.y = JUMP_SPEED;
                players[1].onGround = false;
            }
        }

        if (IsKeyPressed(KEY_P)) {
            gameState = GAME_STATE_PAUSED;
        }

        if (IsKeyPressed(KEY_F3)) {
            debugDrawParticles = !debugDrawParticles;
            if (!debugDrawParticles) {
                debugColorParticlesByVelocity = false;
            }
        }
        if (IsKeyPressed(KEY_F4)) {
            debugColorParticlesByVelocity = !debugColorParticlesByVelocity;
            if (debugColorParticlesByVelocity) {
                debugDrawParticles = true;
            }
        }

        // Shield regeneration
        for (int i = 0; i < activePlayers; i++) {
            if ((float)GetTime() - players[i].last_damage_time > SHIELD_REGEN_DELAY) {
                int max_shield = player_max_shield(&players[i]);
                if (players[i].shield < max_shield) {
                    players[i].shield += 1; // Regenerate 1 shield point per frame
                }
                if (players[i].shield > max_shield) {
                    players[i].shield = max_shield;
                }
            }
        }

        // update players
        for (int i = 0; i < activePlayers; i++) {
            if (playerInput[i] == INPUT_TYPE_KEYBOARD) {
                HandleKeyboardInput(i, dt);
            } else if (playerInput[i] == INPUT_TYPE_GAMEPAD) {
                HandleGamepadInput(i, dt);
            }
            Player *p = &players[i];
            

            // clamp horizontal speed
            {
                float speed = sqrtf(p->vel.x*p->vel.x + p->vel.z*p->vel.z);
                if (speed > MOVE_SPEED * fmaxf(p->kd_ratio,1)) {
                    p->vel.x *= MOVE_SPEED * fmaxf(p->kd_ratio,1) / speed;
                    p->vel.z *= MOVE_SPEED * fmaxf(p->kd_ratio,1) / speed;
                }
            }

            // Block velocity where a neighbor voxel exists in movement direction
            {
                bool neigh[6];
                get_adjacent_voxel_directions(v_add(p->pos,v_mul(p->vel,dt)), neigh);
                // X-axis (+X/neigh[0], -X/neigh[1])
                if ((p->vel.x > 0 && neigh[0]) || (p->vel.x < 0 && neigh[1])) p->vel.x = 0;
                // Y-axis (+Y/neigh[2], -Y/neigh[3])
                if ((p->vel.y > 0 && neigh[2]) || (p->vel.y < 0 && neigh[3])) p->vel.y = 0;
                // Z-axis (+Z/neigh[4], -Z/neigh[5])
                if ((p->vel.z > 0 && neigh[4]) || (p->vel.z < 0 && neigh[5])) p->vel.z = 0;
                if (!neigh[3]){
                    // apply gravity
                    p->vel.y -= GRAVITY*dt;
                    p->onGround = false;
                }else{
                    p->onGround = true;
                }
            }
            

            // apply movement
            p->pos.x += p->vel.x*dt;
            p->pos.y += p->vel.y*dt;
            p->pos.z += p->vel.z*dt;

            // ground clamp
            if (p->pos.y <= BASE_EYE_HEIGHT) {
                p->pos.y = BASE_EYE_HEIGHT;
                p->vel.y = 0;
                p->onGround = true;
            }
            // world bounds clamp for X,Z
            p->pos.x = clampf(p->pos.x, -FLOOR_SIZE+PLAYER_RADIUS, FLOOR_SIZE-PLAYER_RADIUS);
            p->pos.z = clampf(p->pos.z, -FLOOR_SIZE+PLAYER_RADIUS, FLOOR_SIZE-PLAYER_RADIUS);
        }

        activate_static_voxels_near_dynamic();
        batch_glued_dynamic_voxels();
        //update voxel physics
        int subStep = 1;
        for( int i = 0; i < subStep; i++){
            physics_step(dt/subStep);
        }
        update_projectiles(dt);
        simulate_voxel_pbd(dt);
        recycle_dead_voxels();
        log_dynamic_glue_cluster_breaks();
        handle_pbd_projectile_hits();
        voxelHashFramesSinceRebuild++;
        if (voxelHashFramesSinceRebuild >= VOXEL_HASH_REBUILD_INTERVAL) {
            rebuild_voxel_hash();
        }
        deactivate_sleeping_voxels();
        if (debugLogFall) {
            int fall_log_budget = DEBUG_FALL_LOG_BUDGET;
            for (int i = 0; i < voxel_count && fall_log_budget > 0; ++i) {
                Voxel *v = &voxels[i];
                if (!v->simulate) {
                    continue;
                }
                if (v->pos.y < DEBUG_FALL_LOG_THRESHOLD) {
                    TraceLog(LOG_WARNING,
                             "[Fall] voxel=%d span=%d pos=(%.2f,%.2f,%.2f) vel=(%.2f,%.2f,%.2f)",
                             i, v->span,
                             v->pos.x, v->pos.y, v->pos.z,
                             v->vel.x, v->vel.y, v->vel.z);
                    --fall_log_budget;
                }
            }
        }

        // frame by frame debug
        // int frac = 3;
        // if(countFrame==0){
        // integrate_particles(dt);
        // solve_particle_collisions(dt);
        // }

        // if(countFrame==1){
        // for (int i = 0; i < voxel_count; ++i) {
        //     Voxel *voxel = &voxels[i];
        //     if (!voxel->simulate)
        //         continue;
        //     solve_voxel_shape(voxel);
        // }
        // }

        // if(countFrame==2){
        // solve_voxel_glue();
        // }

        // if(countFrame==3){
        // update_particle_velocities(dt);
        // }

        // countFrame += 1;

        // if(countFrame>3){
        // countFrame = 0;
        // }


        activePlayers = clamp_active_players(activePlayers);
        ensure_render_targets(screens, &renderPlayers, &renderW, &renderH, activePlayers);
        Rectangle screenRec = { 0, 0, (float)renderW, (float)-renderH };

        Camera3D cams[MAX_PLAYERS] = { 0 };
        for (int i = 0; i < activePlayers; ++i) {
            cams[i].up = (Vector3){0,1,0};
            cams[i].fovy = 60;
            cams[i].projection = CAMERA_PERSPECTIVE;
            cams[i].position = players[i].pos;
            float yr = DEG2RAD*players[i].yaw, pr = DEG2RAD*players[i].pitch;
            cams[i].target = v_add(players[i].pos, (Vector3){ sinf(-yr)*cosf(pr), sinf(pr), -cosf(yr)*cosf(pr) });
        }

        for (int i = 0; i < activePlayers; ++i) {
            BeginTextureMode(screens[i]);
                ClearBackground(SKYBLUE);
                BeginMode3D(cams[i]);
                    DrawPlane((Vector3){0,0,0}, (Vector2){FLOOR_SIZE*2, FLOOR_SIZE*2}, DARKGRAY);
                    DrawVoxels(cams[i]);
                    draw_players();
                EndMode3D();
                int view_x = 0, view_y = 0, view_w = 0, view_h = 0;
                get_viewport(i, activePlayers, &view_x, &view_y, &view_w, &view_h);
                DrawRectangle(0, 0, view_w, HUD_BAR_HEIGHT, Fade(BLACK, 0.5f));
                DrawText(TextFormat("P%d | Shmush: %d | Kills: %d Deaths: %d",
                                    i + 1, players[i].debrisKills, players[i].kills, players[i].deaths),
                         HUD_PADDING_X, HUD_PADDING_Y, HUD_FONT_SIZE, WHITE);
                draw_hud_bars(&players[i], view_w, view_h);
                {
                    int cx = view_w / 2;
                    int cy = view_h / 2;
                    int cross = 9;
                    Color cross_color = Fade(RED, 0.45f);
                    DrawLine(cx - cross, cy, cx + cross, cy, cross_color);
                    DrawLine(cx, cy - cross, cx, cy + cross, cross_color);
                }
            EndTextureMode();
        }

        BeginDrawing();
            ClearBackground(BLACK);
            for (int i = 0; i < activePlayers; ++i) {
                int view_x = 0, view_y = 0, view_w = 0, view_h = 0;
                get_viewport(i, activePlayers, &view_x, &view_y, &view_w, &view_h);
                DrawTextureRec(screens[i].texture, screenRec, (Vector2){(float)view_x, (float)view_y}, WHITE);
            }
            if (activePlayers == 2) {
                DrawRectangle(SCREEN_WIDTH/2-2, 0, 4, SCREEN_HEIGHT, LIGHTGRAY);
            } else if (activePlayers > 2) {
                DrawRectangle(SCREEN_WIDTH/2-2, 0, 4, SCREEN_HEIGHT, LIGHTGRAY);
                DrawRectangle(0, SCREEN_HEIGHT/2-2, SCREEN_WIDTH, 4, LIGHTGRAY);
            }
            // particle debug text removed
        EndDrawing();
        break;
            case GAME_STATE_PAUSED:
                // Draw pause menu
                BeginDrawing();
                    ClearBackground(RAYWHITE);
                    DrawText("Paused", SCREEN_WIDTH / 2 - MeasureText("Paused", 40) / 2, 100, 40, BLACK);
                    DrawText("Press P to Resume", SCREEN_WIDTH / 2 - MeasureText("Press P to Resume", 20) / 2, 200, 20, DARKGRAY);
                    DrawText("Press Q to Quit", SCREEN_WIDTH / 2 - MeasureText("Press Q to Quit", 20) / 2, 250, 20, DARKGRAY);
                EndDrawing();

                if (IsKeyPressed(KEY_P)) {
                    gameState = GAME_STATE_PLAYING;
                }
                if (IsKeyPressed(KEY_Q)) {
                    gameState = GAME_STATE_MENU;
                }
                break;
            case GAME_STATE_SETTINGS:
                // Draw settings menu
                BeginDrawing();
                    ClearBackground(RAYWHITE);
                    DrawText("Settings", SCREEN_WIDTH / 2 - MeasureText("Settings", 40) / 2, 100, 40, BLACK);
                    DrawText(TextFormat("Active Players: %d", activePlayers), 100, 170, 20, DARKGRAY);
                    DrawText(TextFormat("Player 1 Input: %s", playerInput[0] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 220, 20, DARKGRAY);
                    DrawText(TextFormat("Player 2 Input: %s", playerInput[1] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 250, 20, DARKGRAY);
                    DrawText(TextFormat("Player 3 Input: %s", playerInput[2] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 280, 20, DARKGRAY);
                    DrawText(TextFormat("Player 4 Input: %s", playerInput[3] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 310, 20, DARKGRAY);
                    DrawText(TextFormat("Random Spawn: %s", randomSpawnEnabled ? "ON" : "OFF"), 100, 340, 20, DARKGRAY);
                    DrawText("Press R to toggle random spawns", 100, 370, 20, DARKGRAY);
                    DrawText("Press +/- to change active player count", 100, 400, 20, DARKGRAY);
                    DrawText("Press 1-4 to toggle player input", 100, 430, 20, DARKGRAY);
                    DrawText("Press M to return to Main Menu", 100, 480, 20, DARKGRAY);
                EndDrawing();

                if (IsKeyPressed(KEY_ONE)) {
                    playerInput[0] = 1 - playerInput[0];
                }
                if (IsKeyPressed(KEY_TWO)) {
                    playerInput[1] = 1 - playerInput[1];
                }
                if (IsKeyPressed(KEY_THREE)) {
                    playerInput[2] = 1 - playerInput[2];
                }
                if (IsKeyPressed(KEY_FOUR)) {
                    playerInput[3] = 1 - playerInput[3];
                }
                if (IsKeyPressed(KEY_R)) {
                    randomSpawnEnabled = !randomSpawnEnabled;
                }
                if (IsKeyPressed(KEY_EQUAL) || IsKeyPressed(KEY_KP_ADD)) {
                    activePlayers = clamp_active_players(activePlayers + 1);
                }
                if (IsKeyPressed(KEY_MINUS) || IsKeyPressed(KEY_KP_SUBTRACT)) {
                    activePlayers = clamp_active_players(activePlayers - 1);
                }
                if (IsKeyPressed(KEY_M)) {
                    gameState = GAME_STATE_MENU;
                }
                break;

        }
    }
    // cleanup
    for (int i = 0; i < renderPlayers; ++i) {
        if (screens[i].id != 0) {
            UnloadRenderTexture(screens[i]);
        }
    }
    if (greedyMesh.vertices) {
        UnloadMesh(greedyMesh);
        greedyMesh = (Mesh){ 0 };
    }
    if (greedyMaterialInit) {
        UnloadMaterial(greedyMaterial);
        greedyMaterialInit = false;
    }
    CloseWindow();
    if (debugLogFile) {
        fclose(debugLogFile);
        debugLogFile = NULL;
    }
    return 0;
}
