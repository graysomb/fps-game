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

// Game state enum
typedef enum {
    GAME_STATE_MENU,
    GAME_STATE_PLAYING,
    GAME_STATE_PAUSED,
    GAME_STATE_SETTINGS
} GameState;
static GameState gameState = GAME_STATE_MENU;

// Input type enum
typedef enum {
    INPUT_TYPE_KEYBOARD,
    INPUT_TYPE_GAMEPAD
} InputType;
static InputType playerInput[2] = { INPUT_TYPE_KEYBOARD, INPUT_TYPE_KEYBOARD };

//gcc fps_ray.c -o fps_ray  $(pkg-config --cflags --libs raylib) -lm -Wl,-rpath,/usr/local/lib
// Screen and game constants
#define SCREEN_WIDTH    1600
#define SCREEN_HEIGHT    900
#define MOVE_SPEED      5.0f    // units per second (unused: using acceleration)
#define TURN_SPEED     90.0f    // degrees per second
#define JUMP_SPEED     10.0f    // initial jump velocity
#define GRAVITY         9.8f    // gravity acceleration
#define BASE_EYE_HEIGHT 1.0f    // player eye height above floor
#define ACCELERATION   400.0f    // horizontal acceleration
#define FRICTION       400.0f    // ground friction deceleration
#define TURN_ACCELERATION 400.0f
#define TURN_FRICTION 400.0f
#define PLAYER_RADIUS   0.5f
#define FLOOR_SIZE     20.0f    // half-size of floor in world units
#define PLAYER_SIZE 0.5f
#define PARTICLE_RADIUS (VOXEL_SIZE * 0.5f)
#define VGS_ALPHA 0.5f
#define VGS_BETA 0.5f
#define VGS_ITERS 6
#define VGS_EPS 1e-6f
#define PBD_MAX_STEP_DT 0.005f
#define PBD_SUBSTEPS 3
#define PBD_CONSTRAINT_ITERS 12
#define COLLISION_RELAXATION 0.99f
#define CENTER_RELAXATION 0.9f
#define VELOCITY_DAMPING 0.99f
#define GLUE_RELAXATION 1.0f
#define GLUE_EPS 1e-6f
#define GLUE_BREAK_STRAIN 10.2f
#define VOXEL_SPLIT_STRAIN_THRESHOLD 10.05f
#define VOXEL_SPLIT_SHEAR_THRESHOLD 10.05f

#define MAX_NEIGHBOR_VOXELS 128
#define MAX_FACE_NEIGHBORS   64
#define MAX_SPLIT_CHILDREN    8
static const float GRID_EPSILON = 1e-4f;

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
static int voxelBrushSpan = 4;

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
    int deaths;
    float kd_ratio;
    int health;
    int shield;
    float last_damage_time;
} Player;
static Player players[2];

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
    Color color;
    int type;
    int owner;
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
    Particle particles[8];
    float rest_volume;
    float rest_edge;
    float particle_radius;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;

typedef struct {
    int gx, gy, gz;
    Color color;
} UnitVoxelSeed;

typedef struct {
    UnitVoxelSeed voxels[MAX_VOXELS];
    int count;
} UnitVoxelBuffer;

static void unit_voxel_buffer_clear(UnitVoxelBuffer *buffer) {
    if (!buffer) {
        return;
    }
    buffer->count = 0;
}

static bool unit_voxel_buffer_push(UnitVoxelBuffer *buffer,
                                   int gx, int gy, int gz,
                                   Color color)
{
    if (!buffer || buffer->count >= MAX_VOXELS) {
        return false;
    }
    buffer->voxels[buffer->count++] = (UnitVoxelSeed){
        .gx = gx,
        .gy = gy,
        .gz = gz,
        .color = color
    };
    return true;
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
static int debugGlueBuildLogBudget = 0;
static int debugGlueSolveLogBudget = 0;
static int debugGlueBreakLogBudget = 0;
static const int DEBUG_GLUE_BUILD_LOG_INIT = 64;
static const int DEBUG_GLUE_SOLVE_LOG_INIT = 64;
static const int DEBUG_GLUE_BREAK_LOG_INIT = 16;

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

static int table_get(int x, int y, int z);
static void rebuild_glue_constraints(void);
static void solve_voxel_glue(void);
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
} Patch;

static Patch patches[MAX_VOXELS];
static int   patchCount = 0;
static Mesh  greedyMesh = { 0 };
static bool  meshDirty  = true;

typedef struct {
    int      coarseVoxel;
    int      fineVoxel;
    int      coarseCorner[4];
    float    w[4];
    int      fineCorner;
    uint8_t  coarseMask;
    uint8_t  fineMask;
    bool     active;
} GlueConstraint;

static GlueConstraint glueConstraints[MAX_VOXELS * 12];
static int glueConstraintCount = 0;

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
    size_t   h = hashVoxelKey(k);

    while (1) {
        uint64_t bk = table[h].key;
        if (bk == 0)      return -1;          // empty bucket ⇒ miss
        if (bk == k)      return table[h].idx;  // exact key ⇒ hit
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
}


// Check occupancy
static bool occupied(int x,int y,int z){ return  table_get(x,y,z)>=0; }

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
    v->color = color;
    v->type = type;
    v->owner = owner;
    memset(v->surface, 0, sizeof v->surface);
    v->min_gx = v->min_gy = v->min_gz = INT_MAX;
    v->max_gx = v->max_gy = v->max_gz = INT_MIN;
    v->rest_min_gx = v->rest_min_gy = v->rest_min_gz = INT_MAX;
    v->rest_max_gx = v->rest_max_gy = v->rest_max_gz = INT_MIN;
    v->gx = (int)floorf(px / VOXEL_SIZE);
    v->gy = (int)floorf(py / VOXEL_SIZE);
    v->gz = (int)floorf(pz / VOXEL_SIZE);
    v->span = span;
    v->rest_edge = edge;
    v->rest_volume = edge * edge * edge;
    v->particle_radius = 0.5f * edge;
    for (int i = 0; i < 8; ++i) {
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
    int rest_minx, rest_maxx, rest_miny, rest_maxy, rest_minz, rest_maxz;
    voxel_compute_bounds(v, &rest_minx, &rest_maxx, &rest_miny, &rest_maxy, &rest_minz, &rest_maxz);
    v->rest_min_gx = rest_minx;
    v->rest_max_gx = rest_maxx;
    v->rest_min_gy = rest_miny;
    v->rest_max_gy = rest_maxy;
    v->rest_min_gz = rest_minz;
    v->rest_max_gz = rest_maxz;
}

static int addVoxelSized(float px, float py, float pz, bool fixed, bool simulate,
                         Color color, int type, int span) {
    if (voxel_count >= MAX_VOXELS) return -1;
    int idx = voxel_count++;
    Voxel *v = &voxels[idx];
    init_voxel_struct(v, px, py, pz, fixed, simulate, color, type, span, -1);
    voxel_table_register(v, idx);
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
                                             int type)
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
    if (!occupied || !consumed || !color_grid) {
        free(occupied);
        free(consumed);
        free(color_grid);
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
                addVoxelSized(px, py, pz, fixed, simulate, color, type, span);
                ++spawned;
            }
        }
    }

    free(occupied);
    free(consumed);
    free(color_grid);

    TraceLog(LOG_INFO,
             "[Multiscale] Converted %d unit voxels into %d span clusters",
             buffer->count, spawned);
    return spawned;
}

static unsigned char mix_color_channel(float a, float b, float t) {
    return (unsigned char)clampf(mixf(a, b, t), 0.0f, 255.0f);
}

static void build_oblique_voxel_pyramid(UnitVoxelBuffer *buffer) {
    if (!buffer) {
        return;
    }

    const int pyramid_height = 3;
    const int base_length = 4;
    const int base_width = 4;
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
                if (!unit_voxel_buffer_push(buffer, gx, layer_y, gz, col)) {
                    TraceLog(LOG_WARNING, "[Pyramid] Unit voxel buffer full");
                    return;
                }
            }
        }
    }
}

// Build static demo cube of voxels
static void buildDemo(void) {
    // Floor
    int M = (int)(2.0f * FLOOR_SIZE / VOXEL_SIZE);
    /* for (int x = 0; x <= M; x++) {
        for (int z = 0; z <= M; z++) {
            float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            addVoxel(px, 0, pz, true, false, (Color){ 150, 150, 150, 255 }, 0);
        }
    } */

    // Pillars
    // int pillar_height = 15; // 45 - 10
    // int pillar_radius = 3;
    // int pillar_positions[4][2] = {
    //     { M / 4, M / 4 },
    //     { M / 4, 3 * M / 4 },
    //     { 3 * M / 4, M / 4 },
    //     { 3 * M / 4, 3 * M / 4 }
    // };

    // for (int p = 0; p < 4; p++) {
    //     int cx = pillar_positions[p][0];
    //     int cz = pillar_positions[p][1];
    //     for (int y = 1; y <= pillar_height; y++) {
    //         for (int dx = -pillar_radius; dx <= pillar_radius; dx++) {
    //             for (int dz = -pillar_radius; dz <= pillar_radius; dz++) {
    //                 if (dx*dx + dz*dz > pillar_radius*pillar_radius) continue; // circular pillar
    //                 float px = (cx + dx + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //                 float py = (y + 0.5f) * VOXEL_SIZE;
    //                 float pz = (cz + dz + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //                 addVoxel(px, py, pz, false, true, (Color){ 200, 100, 50, 255 }, 0);
    //             }
    //         }
    //     }
    //}

    // Central platform (n=1)
    // int platform_size = 2;
    // int platform_height = 1; // 15 / 3
    // int platform_base_height = 10; // to keep top at same level (21)
    // for (int y = platform_base_height; y <= platform_base_height + platform_height; y++) {
    //     for (int x = M/2 - platform_size/2; x <= M/2 + platform_size/2; x++) {
    //         for (int z = M/2 - platform_size/2; z <= M/2 + platform_size/2; z++) {
    //             float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //             float py = (y + 0.5f) * VOXEL_SIZE;
    //             float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    //             addVoxel(px, py, pz, false, true, (Color){ 100, 200, 100, 255 }, 0);
    //             //addVoxelSized(px, py, pz, false, true, (Color){ 100, 200, 100, 255 }, 0, 1);
    //         }
    //     }
    // }
    UnitVoxelBuffer pyramid_units;
    unit_voxel_buffer_clear(&pyramid_units);
    build_oblique_voxel_pyramid(&pyramid_units);
    emit_multiscale_voxels_from_units(&pyramid_units, false, true, 0);
}



static int first_voxel_hit(Ray ray, float t_max, int ignore_id);
static void UpdateKdRatio(int player_index);

// Reset game: players and voxels
static void ResetGame(void) {
    // init players
    for (int i = 0; i < 2; i++) {
        //players[i].pos = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
        players[i].pos = (Vector3){ 0, BASE_EYE_HEIGHT, -9 };
        players[i].yaw = (i == 0) ? 0 : 180;
        players[i].pitch = 0;
        players[i].yaw_vel = 0;
        players[i].pitch_vel = 0;
        players[i].vel = (Vector3){0,0,0};
        players[i].onGround = true;
        players[i].vType = 0;
        players[i].kills = 0;
        players[i].deaths = 0;
        UpdateKdRatio(i);
        players[i].health = BASE_HEALTH;
        players[i].shield = BASE_SHIELD;
        players[i].last_damage_time = 0.0f;
    }
    // clear voxels
    voxel_count = 0;
    // clear hash
    memset(table, 0, sizeof(table));
    // build static blocks
    buildDemo();
    rebuild_all_voxel_surfaces();
    rebuild_glue_constraints();
    meshDirty = true;

}

static void UpdateKdRatio(int player_index) {
    Player *p = &players[player_index];
    p->kd_ratio = (float)(p->kills + 1) / (p->deaths + 1);
}

static void rebuild_voxel_hash(void) {
    memset(table, 0, sizeof(table));
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        v->gx = (int)floorf(v->pos.x / VOXEL_SIZE);
        v->gy = (int)floorf(v->pos.y / VOXEL_SIZE);
        v->gz = (int)floorf(v->pos.z / VOXEL_SIZE);
        voxel_table_register(v, i);
    }
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

// Predict positions for the next step (equivalent to the GPU PredictPositions kernel).
static void integrate_particles(float dt) {
    const Vector3 gravity = { 0.0f, -GRAVITY*1.0f, 0.0f };
    const float dt_sq = dt * dt;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];

        for (int j = 0; j < 8; ++j) {
            Particle *p = &voxel->particles[j];

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

// Voxel Gram-Schmidt shape matching (Algorithm 1 in the paper) keeps each cell near-rest.
static void solve_voxel_shape(Voxel *voxel) {
    bool has_dynamic = false;
    Vector3 p[8];
    float w[8];

    for (int i = 0; i < 8; ++i) {
        Particle *part = &voxel->particles[i];
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

        float lenp0 = 4.0f * v_length(v0);
        float lenp1 = 4.0f * v_length(v1);
        float lenp2 = 4.0f * v_length(v2);
        float r_v = 1.0f;
        float denom = lenp0 * lenp1 * lenp2;
        float rest_demom = rest_edge * rest_edge * rest_edge;
        if (fabs(denom-rest_demom) > VGS_EPS) {
            float ratio = (rest_edge * rest_edge * rest_edge) / denom;
            float root = cbrtf(fabsf(ratio));
            r_v = (ratio < 0.0f) ? -root : root;
        }

        float target0 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp0 * r_v));
        float target1 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp1 * r_v));
        float target2 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * (lenp2 * r_v));

        if (fabs(len0-target0) > VGS_EPS) u0 = v_mul(u0, target0 / len0);
        if (fabs(len1-target1) > VGS_EPS) u1 = v_mul(u1, target1 / len1);
        if (fabs(len2-target2) > VGS_EPS) u2 = v_mul(u2, target2 / len2);

        // Volume correction mirrors the GPU "ResizeVoxelBasis" stage.
        float volume = v_dot(v_cross(u0, u1), u2);
        if (fabsf(volume-rest_volume) > VGS_EPS) {
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

    voxel->pos = centroid;
    for (int i = 0; i < 8; ++i) {
        voxel->particles[i].predicted_pos = p[i];
    }
}

static void solve_voxel_glue(void) {
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

        Vector3 blended = { 0.0f, 0.0f, 0.0f };
        for (int k = 0; k < 4; ++k) {
            blended = v_add(blended, v_mul(coarseParticles[k]->predicted_pos, weights[k]));
        }
        Vector3 C = v_sub(blended, fineParticle->predicted_pos);
        float violation = v_length(C);
        float rest_edge_min = fminf(coarse->rest_edge, fine->rest_edge);
        float break_distance = GLUE_BREAK_STRAIN * rest_edge_min;
        if (violation > break_distance) {
            if (debugLogGlue && debugGlueBreakLogBudget > 0) {
                TraceLog(LOG_INFO,
                         "[GlueBreak] pair=(%d,%d) spans=(%d,%d) violation=%.5f break=%.5f uvPred=(%.3f,%.3f) rawUV=(%.3f,%.3f) baryValid=%s weights=(%.3f,%.3f,%.3f,%.3f) coarsePos=(%.2f,%.2f,%.2f) finePos=(%.2f,%.2f,%.2f)",
                         gc->coarseVoxel, gc->fineVoxel,
                         voxel_span_for_glue(coarse), voxel_span_for_glue(fine),
                         violation, break_distance,
                         solveU, solveV,
                         solveRawU, solveRawV,
                         baryValid ? "true" : "false",
                         weights[0], weights[1], weights[2], weights[3],
                         coarse->pos.x, coarse->pos.y, coarse->pos.z,
                         fine->pos.x, fine->pos.y, fine->pos.z);
                --debugGlueBreakLogBudget;
            }
            gc->active = false;
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

        Vector3 lambda = v_mul(C, -GLUE_RELAXATION / invMassSum);
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
            TraceLog(LOG_INFO,
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
        gc->active = true;

        if (debugLogGlue && debugGlueBuildLogBudget > 0) {
            TraceLog(LOG_INFO,
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
    }
}

static void rebuild_glue_constraints(void) {
    glueConstraintCount = 0;
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
}

static void deactivate_glue_constraints_between(int a, int b) {
    for (int g = 0; g < glueConstraintCount; ++g) {
        GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }
        if ((gc->coarseVoxel == a && gc->fineVoxel == b) ||
            (gc->coarseVoxel == b && gc->fineVoxel == a))
        {
            gc->active = false;
        }
    }
}

static int gather_glued_neighbors(int voxel_idx, int *out, int max_out) {
    if (!out || max_out <= 0) {
        return 0;
    }
    int count = 0;
    for (int g = 0; g < glueConstraintCount; ++g) {
        const GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }
        int neighbor = -1;
        if (gc->coarseVoxel == voxel_idx) {
            neighbor = gc->fineVoxel;
        } else if (gc->fineVoxel == voxel_idx) {
            neighbor = gc->coarseVoxel;
        }
        if (neighbor < 0) {
            continue;
        }
        bool seen = false;
        for (int n = 0; n < count; ++n) {
            if (out[n] == neighbor) {
                seen = true;
                break;
            }
        }
        if (!seen && count < max_out) {
            out[count++] = neighbor;
        }
    }
    return count;
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

// True when two voxels share any active glue constraint (face coupling).
static bool voxels_are_glued(int voxel_idx_a, int voxel_idx_b) {
    if (voxel_idx_a == voxel_idx_b) {
        return false;
    }

    for (int g = 0; g < glueConstraintCount; ++g) {
        const GlueConstraint *gc = &glueConstraints[g];
        if (!gc->active) {
            continue;
        }

        bool forward = (gc->coarseVoxel == voxel_idx_a && gc->fineVoxel == voxel_idx_b);
        bool reverse = (gc->coarseVoxel == voxel_idx_b && gc->fineVoxel == voxel_idx_a);
        if (forward || reverse) {
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
    (void)dt;

    const float half_player = PLAYER_SIZE * 0.5f;
    const float omega = COLLISION_RELAXATION;
    const float eps = 1e-6f;
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
        float voxel_radius = voxel_particle_radius(voxel);
        float terrain_limit = FLOOR_SIZE - voxel_radius;

        for (int j = 0; j < 8; ++j) {
            Particle *p = &voxel->particles[j];

            Vector3 pos = p->predicted_pos;

            if (pos.y < voxel_radius) {
                pos.y = voxel_radius;
            }

            pos.x = clampf(pos.x, -terrain_limit, terrain_limit);
            pos.z = clampf(pos.z, -terrain_limit, terrain_limit);

            // Interactions with player bounding boxes (kept for gameplay parity).
            for (int player_idx = 0; player_idx < 2; ++player_idx) {
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
        }
    }

    // Particle-particle collisions using a symmetric correction identical to the compute shader.
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxelA = &voxels[i];
        float radiusA = voxel_particle_radius(voxelA);
        int neighbor_ids[MAX_NEIGHBOR_VOXELS];
        int neighbor_count = gather_neighbor_voxels(voxelA, i, neighbor_ids, MAX_NEIGHBOR_VOXELS);

        for (int j = 0; j < 8; ++j) {
            Particle *pa = &voxelA->particles[j];
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

                bool glued = voxels_are_glued(i, neighbor_idx);
                bool share_edge_corner = voxels_share_edge_or_corner(voxelA, voxelB);
                if (glued || share_edge_corner) {
                    continue;
                }

                for (int q = 0; q < 8; ++q) {
                    if (neighbor_idx == i && q <= j) {
                        continue;
                    }

                    if (particles_are_glued_pair(i, j, neighbor_idx, q)) {
                        continue;
                    }

                    Particle *pb = &voxelB->particles[q];
                    float wb = pb->inv_mass;

                    float w_sum = wa + wb;
                    if (w_sum <= 0.0f) {
                        continue;
                    }

                    Vector3 delta = v_sub(pa->predicted_pos, pb->predicted_pos);
                    float dist_sq = v_dot(delta, delta);
                    float target_dist = radiusA + radiusB;

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

static void update_particle_velocities(float dt) {
    float inv_dt = (dt > 0.0f) ? 1.0f / dt : 0.0f;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        Vector3 centroid = { 0.0f, 0.0f, 0.0f };
        Vector3 prev_centroid = { 0.0f, 0.0f, 0.0f };

        for (int j = 0; j < 8; ++j) {
            Particle *p = &voxel->particles[j];
            centroid = v_add(centroid, p->predicted_pos);
            prev_centroid = v_add(prev_centroid, p->prev_pos);

            Vector3 new_pos = p->predicted_pos;
            Vector3 delta = v_sub(new_pos, p->prev_pos);

            if (p->inv_mass > 0.0f) {
                p->vel = v_mul(delta, inv_dt);
            } else {
                p->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
            }

            p->pos = new_pos;
        }

        centroid = v_mul(centroid, 1.0f / 8.0f);
        prev_centroid = v_mul(prev_centroid, 1.0f / 8.0f);
        voxel->vel = v_mul(v_sub(centroid, prev_centroid), inv_dt);
        voxel->pos = centroid;
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
}

static int split_voxel_at(int idx, float dt, int *out_children, int max_children) {
    if (idx < 0 || idx >= voxel_count || !out_children || max_children < MAX_SPLIT_CHILDREN) {
        return 0;
    }

    Voxel parent = voxels[idx];
    if (parent.span < 2 || (parent.span % 2) != 0) {
        return 0;
    }
    const int additional_children = MAX_SPLIT_CHILDREN - 1;
    if (voxel_count + additional_children > MAX_VOXELS) {
        return 0;
    }

    int child_span = parent.span / 2;
    float child_edge = VOXEL_SIZE * (float)child_span;
    float parent_half = 0.5f * parent.rest_edge;
    float start_x = parent.pos.x - parent_half + 0.5f * child_edge;
    float start_y = parent.pos.y - parent_half + 0.5f * child_edge;
    float start_z = parent.pos.z - parent_half + 0.5f * child_edge;

    int child_counter = 0;
    int next_index = voxel_count;
    for (int iz = 0; iz < 2; ++iz) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int ix = 0; ix < 2; ++ix) {
                float cx = start_x + ix * child_edge;
                float cy = start_y + iy * child_edge;
                float cz = start_z + iz * child_edge;

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
                apply_uniform_velocity(child, parent.vel, dt);
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
        if (!voxel->simulate || voxel->span <= 1 || (voxel->span % 2) != 0) {
            ++i;
            continue;
        }

        float strain = 0.0f;
        float shear = 0.0f;
        voxel_measure_strain(voxel, &strain, &shear);
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
        integrate_particles(sub_dt);
        solve_particle_collisions(sub_dt);
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
                if (!voxel->simulate)
                    continue;
                solve_voxel_shape(voxel);
            }
        solve_voxel_glue();
        }
        update_particle_velocities(sub_dt);
    }
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
    float yawRad = DEG2RAD * p->yaw;
    float pitchRad = DEG2RAD * p->pitch;
    Vector3 dir = { sinf(-yawRad)*cosf(pitchRad), sinf(pitchRad), -cosf(yawRad)*cosf(pitchRad) };
    Vector3 start = v_add(p->pos, v_mul(dir, 0.8f));
    Color col = (p->vType==0? RED : BLUE);
    int vix = addVoxel(start.x, start.y, start.z, false, true, col, p->vType);
    if (vix >= 0) {
        Voxel *shot = &voxels[vix];
        Vector3 vel = v_mul(dir, 50.0f);
        shot->vel = vel;
        shot->owner = idx;
        for (int i = 0; i < 8; ++i) {
            shot->particles[i].vel = vel;
        }
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

static void drawCubeMan(const Voxel *voxel)
{
    Vector3 v[8];
    for (int i = 0; i < 8; ++i) v[i] = voxel->particles[i].pos;

    rlColor4ub(voxel->color.r, voxel->color.g, voxel->color.b, voxel->color.a);

    rlNormal3f(0.0f, 0.0f, 1.0f);
    rlVertex3f(v[4].x, v[4].y, v[4].z);
    rlVertex3f(v[5].x, v[5].y, v[5].z);
    rlVertex3f(v[7].x, v[7].y, v[7].z);
    rlVertex3f(v[4].x, v[4].y, v[4].z);
    rlVertex3f(v[7].x, v[7].y, v[7].z);
    rlVertex3f(v[6].x, v[6].y, v[6].z);

    rlNormal3f(0.0f, 0.0f, -1.0f);
    rlVertex3f(v[1].x, v[1].y, v[1].z);
    rlVertex3f(v[0].x, v[0].y, v[0].z);
    rlVertex3f(v[2].x, v[2].y, v[2].z);
    rlVertex3f(v[1].x, v[1].y, v[1].z);
    rlVertex3f(v[2].x, v[2].y, v[2].z);
    rlVertex3f(v[3].x, v[3].y, v[3].z);

    rlNormal3f(0.0f, 1.0f, 0.0f);
    rlVertex3f(v[6].x, v[6].y, v[6].z);
    rlVertex3f(v[7].x, v[7].y, v[7].z);
    rlVertex3f(v[3].x, v[3].y, v[3].z);
    rlVertex3f(v[6].x, v[6].y, v[6].z);
    rlVertex3f(v[3].x, v[3].y, v[3].z);
    rlVertex3f(v[2].x, v[2].y, v[2].z);

    rlNormal3f(0.0f, -1.0f, 0.0f);
    rlVertex3f(v[0].x, v[0].y, v[0].z);
    rlVertex3f(v[1].x, v[1].y, v[1].z);
    rlVertex3f(v[5].x, v[5].y, v[5].z);
    rlVertex3f(v[0].x, v[0].y, v[0].z);
    rlVertex3f(v[5].x, v[5].y, v[5].z);
    rlVertex3f(v[4].x, v[4].y, v[4].z);

    rlNormal3f(1.0f, 0.0f, 0.0f);
    rlVertex3f(v[5].x, v[5].y, v[5].z);
    rlVertex3f(v[1].x, v[1].y, v[1].z);
    rlVertex3f(v[3].x, v[3].y, v[3].z);
    rlVertex3f(v[5].x, v[5].y, v[5].z);
    rlVertex3f(v[3].x, v[3].y, v[3].z);
    rlVertex3f(v[7].x, v[7].y, v[7].z);

    rlNormal3f(-1.0f, 0.0f, 0.0f);
    rlVertex3f(v[0].x, v[0].y, v[0].z);
    rlVertex3f(v[4].x, v[4].y, v[4].z);
    rlVertex3f(v[6].x, v[6].y, v[6].z);
    rlVertex3f(v[0].x, v[0].y, v[0].z);
    rlVertex3f(v[6].x, v[6].y, v[6].z);
    rlVertex3f(v[2].x, v[2].y, v[2].z);
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
                
                int baseIdx = 0;
                switch (plane) {
                    case 0: baseIdx = table_get(pt->i0, pt->j0, layer); break;
                    case 1: baseIdx = table_get(pt->i0, layer, pt->j0); break;
                    case 2: baseIdx = table_get(layer, pt->i0, pt->j0); break;
                }
                pt->col = voxels[baseIdx].color;
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

// Draw all voxels via greedy mesh instead of per-voxel raycasting
static void DrawVoxels(Camera3D cam) {
    /* (void)cam;
    if (meshDirty) {
        if (greedyMesh.vertices) UnloadMesh(greedyMesh);
        greedyMesh = gen_greedy_mesh();
        meshDirty = false;
    }
    rlDisableBackfaceCulling();
    DrawMesh(greedyMesh, LoadMaterialDefault(), MatrixIdentity());
    // Draw voxel grid lines atop merged quads
    rlBegin(RL_LINES);
    rlColor4ub(0, 0, 0, 255);
    for (int p = 0; p < patchCount; p++) {
        Patch *pt = &patches[p];
        Vector3 origin, iu, ju;
        switch (pt->plane) {
            case 0: // XY
                origin = (Vector3){ pt->i0*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE, 0, 0 };
                ju = (Vector3){ 0, VOXEL_SIZE, 0 };
                break;
            case 1: // XZ
                origin = (Vector3){ pt->i0*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE, 0, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE };
                break;
            default: // YZ
                origin = (Vector3){ (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    pt->i0*VOXEL_SIZE,
                                    pt->j0*VOXEL_SIZE };
                iu = (Vector3){ 0, VOXEL_SIZE, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE };
                break;
        }
        for (int ix = 0; ix <= pt->di; ix++) {
            Vector3 a = v_add(origin, v_mul(iu, ix));
            Vector3 b = v_add(a, v_mul(ju, pt->dj));
            rlVertex3f(a.x, a.y, a.z); rlVertex3f(b.x, b.y, b.z);
        }
        for (int jy = 0; jy <= pt->dj; jy++) {
            Vector3 a = v_add(origin, v_mul(ju, jy));
            Vector3 b = v_add(a, v_mul(iu, pt->di));
            rlVertex3f(a.x, a.y, a.z); rlVertex3f(b.x, b.y, b.z);
        }
    }
    rlEnd();
    rlEnableBackfaceCulling();*/

    //draw moving voxels
    rlBegin(RL_TRIANGLES);
    /*for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (v->simulate) {
            drawCubeMan(v);
        }
    }*/
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        drawCubeMan(v);
    }
    rlEnd();
    rlBegin(RL_LINES);
    rlColor4ub(0, 0, 0, 255);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        drawCubeEdges(v);
    }
    rlEnd();

    if (debugDrawParticles) {
        draw_particle_debug();
    }
}

static void draw_players(void) {
    for (int i = 0; i < 2; i++) {
        Player *p = &players[i];
        DrawCube(p->pos, PLAYER_SIZE,PLAYER_SIZE, PLAYER_SIZE, BLACK);
        DrawCubeWires(p->pos, PLAYER_SIZE,PLAYER_SIZE,PLAYER_SIZE, BLACK);
        if (p->shield > 0) {
            float shield_percentage = (float)p->shield / BASE_SHIELD;
            float fluctuation = (1.0f - shield_percentage) * 100.0f;
            Color shield_color = {
                (unsigned char)randomInRange(0, fluctuation),
                (unsigned char)randomInRange(0, fluctuation),
                (unsigned char)clampf(255 - randomInRange(0, fluctuation), 0, 255),
                128
            };
            DrawCube(p->pos, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, shield_color);
            DrawCubeWires(p->pos, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, PLAYER_SIZE + 0.2f, Fade(shield_color, 0.8f));
        }
    }
}

static void HandleKeyboardInput(int i, float dt) {
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
    // init window and render textures
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Split-Screen FPS (raylib)");
    SetTargetFPS(60);
    // seed RNG
    srand((unsigned)time(NULL));
    // reset game state
    ResetGame();
    // prepare split-screen render textures
    RenderTexture2D screen0 = LoadRenderTexture(SCREEN_WIDTH/2, SCREEN_HEIGHT);
    RenderTexture2D screen1 = LoadRenderTexture(SCREEN_WIDTH/2, SCREEN_HEIGHT);
    Rectangle screenRec = { 0, 0, (float)screen0.texture.width, (float)-screen0.texture.height };
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
        // input: shooting and jump
        if (playerInput[0] == INPUT_TYPE_KEYBOARD && IsKeyPressed(KEY_LEFT_CONTROL))  FireVoxel(0);
        if (playerInput[0] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(0, GAMEPAD_BUTTON_RIGHT_TRIGGER_2)) FireVoxel(0);
        if (playerInput[1] == INPUT_TYPE_KEYBOARD && IsKeyPressed(KEY_RIGHT_CONTROL)) FireVoxel(1);
        if (playerInput[1] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(1, GAMEPAD_BUTTON_RIGHT_TRIGGER_2)) FireVoxel(1);

        if (IsKeyPressed(KEY_Q)) players[0].vType = 1-players[0].vType;
        if (IsKeyPressed(KEY_U)) players[0].vType = 1-players[0].vType;
        if (playerInput[0] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(0, GAMEPAD_BUTTON_LEFT_TRIGGER_1)) players[0].vType = 1-players[0].vType;
        if (playerInput[1] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(1, GAMEPAD_BUTTON_LEFT_TRIGGER_1)) players[1].vType = 1-players[1].vType;

        if (playerInput[0] == INPUT_TYPE_KEYBOARD && IsKeyPressed(KEY_SPACE) && players[0].onGround) { players[0].vel.y = JUMP_SPEED; players[0].onGround = false; }
        if (playerInput[0] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(0, GAMEPAD_BUTTON_RIGHT_FACE_DOWN) && players[0].onGround) { players[0].vel.y = JUMP_SPEED; players[0].onGround = false; }
        if (playerInput[1] == INPUT_TYPE_KEYBOARD && IsKeyPressed(KEY_RIGHT_SHIFT) && players[1].onGround) { players[1].vel.y = JUMP_SPEED; players[1].onGround = false; }
        if (playerInput[1] == INPUT_TYPE_GAMEPAD && IsGamepadButtonPressed(1, GAMEPAD_BUTTON_RIGHT_FACE_DOWN) && players[1].onGround) { players[1].vel.y = JUMP_SPEED; players[1].onGround = false; }

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
        for (int i = 0; i < 2; i++) {
            if ((float)GetTime() - players[i].last_damage_time > SHIELD_REGEN_DELAY) {
                if (players[i].shield < BASE_SHIELD) {
                    players[i].shield += 1; // Regenerate 1 shield point per frame
                }
                if (players[i].shield > BASE_SHIELD) {
                    players[i].shield = BASE_SHIELD;
                }
            }
        }

        // update players
        for (int i = 0; i < 2; i++) {
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
        //update voxel physics
        int subStep = 3;
        for( int i = 0; i < subStep; i++){
            physics_step(dt/subStep);
        }
        simulate_voxel_pbd(dt);

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


        // setup cameras
        Camera3D cam0 = {0}, cam1 = {0};
        cam0.up = cam1.up = (Vector3){0,1,0};
        cam0.fovy = cam1.fovy = 60;
        cam0.projection = cam1.projection = CAMERA_PERSPECTIVE;
        // camera 0
        cam0.position = players[0].pos;
        {
            float yr = DEG2RAD*players[0].yaw, pr = DEG2RAD*players[0].pitch;
            cam0.target = v_add(players[0].pos, (Vector3){ sinf(-yr)*cosf(pr), sinf(pr), -cosf(yr)*cosf(pr) });
        }
        // camera 1
        cam1.position = players[1].pos;
        {
            float yr = DEG2RAD*players[1].yaw, pr = DEG2RAD*players[1].pitch;
            cam1.target = v_add(players[1].pos, (Vector3){ sinf(-yr)*cosf(pr), sinf(pr), -cosf(yr)*cosf(pr) });
        }
        // render to textures
        BeginTextureMode(screen0);
            ClearBackground(SKYBLUE);
            BeginMode3D(cam0);
                DrawPlane((Vector3){0,0,0}, (Vector2){FLOOR_SIZE*2, FLOOR_SIZE*2}, DARKGRAY);
                DrawVoxels(cam0);
                draw_players();
            EndMode3D();
            // UI p1
            DrawRectangle(0,0, SCREEN_WIDTH/2, 4.0f, Fade(BLACK, 0.5f));
            DrawText(TextFormat("P1 | Kills: %d Deaths: %d | Health: %d | Shield: %d", players[0].kills, players[0].deaths, players[0].health, players[0].shield), 10, 10, 20, WHITE);
            DrawLine(SCREEN_WIDTH/4-10, SCREEN_HEIGHT/2, SCREEN_WIDTH/4+10, SCREEN_HEIGHT/2, WHITE);
            DrawLine(SCREEN_WIDTH/4, SCREEN_HEIGHT/2-10, SCREEN_WIDTH/4, SCREEN_HEIGHT/2+10, WHITE);
        EndTextureMode();
        BeginTextureMode(screen1);
            ClearBackground(SKYBLUE);
            BeginMode3D(cam1);
                DrawPlane((Vector3){0,0,0}, (Vector2){FLOOR_SIZE*2, FLOOR_SIZE*2}, DARKGRAY);
                DrawVoxels(cam1);
                draw_players();
            EndMode3D();
            // UI p2
            DrawRectangle(0,0, SCREEN_WIDTH/2, 4.0f, Fade(BLACK, 0.5f));
            DrawText(TextFormat("P2 | Kills: %d Deaths: %d | Health: %d | Shield: %d", players[1].kills, players[1].deaths, players[1].health, players[1].shield), 10, 10, 20, WHITE);
            DrawLine(SCREEN_WIDTH/4-10, SCREEN_HEIGHT/2, SCREEN_WIDTH/4+10, SCREEN_HEIGHT/2, WHITE);
            DrawLine(SCREEN_WIDTH/4, SCREEN_HEIGHT/2-10, SCREEN_WIDTH/4, SCREEN_HEIGHT/2+10, WHITE);
        EndTextureMode();
        // draw both splits
        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureRec(screen0.texture, screenRec, (Vector2){0,0}, WHITE);
            DrawTextureRec(screen1.texture, screenRec, (Vector2){SCREEN_WIDTH/2,0}, WHITE);
            DrawRectangle(SCREEN_WIDTH/2-2, 0, 4, SCREEN_HEIGHT, LIGHTGRAY);
            DrawText(TextFormat("Particles (F3): %s", debugDrawParticles ? "ON" : "OFF"), 20, 20, 20, LIGHTGRAY);
            if (debugDrawParticles) {
                DrawText(TextFormat("Velocity Heatmap (F4): %s", debugColorParticlesByVelocity ? "ON" : "OFF"), 20, 44, 20, LIGHTGRAY);
            }
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
                    DrawText(TextFormat("Player 1 Input: %s", playerInput[0] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 200, 20, DARKGRAY);
                    DrawText(TextFormat("Player 2 Input: %s", playerInput[1] == INPUT_TYPE_KEYBOARD ? "Keyboard" : "Gamepad"), 100, 250, 20, DARKGRAY);
                    DrawText("Press 1 to toggle Player 1 input", 100, 350, 20, DARKGRAY);
                    DrawText("Press 2 to toggle Player 2 input", 100, 400, 20, DARKGRAY);
                    DrawText("Press M to return to Main Menu", 100, 500, 20, DARKGRAY);
                EndDrawing();

                if (IsKeyPressed(KEY_ONE)) {
                    playerInput[0] = 1 - playerInput[0];
                }
                if (IsKeyPressed(KEY_TWO)) {
                    playerInput[1] = 1 - playerInput[1];
                }
                if (IsKeyPressed(KEY_M)) {
                    gameState = GAME_STATE_MENU;
                }
                break;

        }
    }
    // cleanup
    UnloadRenderTexture(screen0);
    UnloadRenderTexture(screen1);
    CloseWindow();
    return 0;
}
