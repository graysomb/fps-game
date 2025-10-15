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

// KD-stats constants
#define BASE_HEALTH 100
#define BASE_SHIELD 150
#define SHIELD_REGEN_DELAY 5.0f
#define VOXEL_DAMAGE 50

// Voxel physics constants
#define MAX_VOXELS    131072
#define HASH_SIZE     131072    // must be power of two
#define VOXEL_SIZE     0.2f    // size of each voxel cube

#define VGS_ALPHA_DEFAULT       0.5f
#define VGS_BETA_DEFAULT        1.0f
#define VGS_LOCAL_ITERATIONS    3
#define CONSTRAINT_ITERATIONS   3

#define VOXEL_PARTICLE_RADIUS   (VOXEL_SIZE * 0.5f)
#define HALF_VOXEL_SIZE         (VOXEL_SIZE * 0.5f)

#define STRAIN_LIMIT_BASE_TENSION  0.35f
#define STRAIN_LIMIT_BASE_COMPRESSION (-0.25f)
#define STRAIN_LIMIT_SHELL_BONUS   1000.0f

#define MAX_VOXEL_LINEAR_SPEED   50.0f

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

// Voxel structure
typedef struct {
    Vector3 pos;
    Vector3 vel;
    bool simulate;
    bool fixed;
    Color color;
    int type;
    int owner;
    /* Visible faces mask: surface[i]=true if face i is visible:
       0=+X,1=-X,2=+Y,3=-Y,4=+Z,5=-Z */
    bool surface[6];
    int  gx, gy, gz;
    Vector3 corners[8];
    Vector3 corners_prev[8];
    Vector3 corner_vel[8];
    float  inv_mass[8];
    float  rest_volume;
    float  particle_radius;
    float  distance_to_surface;
    bool   skeleton;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;

static const Vector3 voxel_corner_offsets[8] = {
    { -HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE },
    {  HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE },
    { -HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE },
    {  HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE },
    { -HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE },
    {  HALF_VOXEL_SIZE, -HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE },
    { -HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE },
    {  HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE,  HALF_VOXEL_SIZE }
};

static const int face_axis_offsets[3][3] = {
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 }
};

static const int face_vertex_pairs[3][4][2] = {
    { {1,0}, {3,2}, {5,4}, {7,6} }, // +X face against neighbor's -X face
    { {2,0}, {3,1}, {6,4}, {7,5} }, // +Y
    { {4,0}, {5,1}, {6,2}, {7,3} }  // +Z
};

static const int axis_edge_pairs[3][4][2] = {
    { {0,1}, {2,3}, {4,5}, {6,7} }, // X edges
    { {0,2}, {1,3}, {4,6}, {5,7} }, // Y edges
    { {0,4}, {1,5}, {2,6}, {3,7} }  // Z edges
};

typedef struct {
    int   a;
    int   b;
    float strain_compressive;
    float strain_tensile;
    bool  active;
} FaceConstraint;

static FaceConstraint face_constraints[3][MAX_VOXELS];
static int face_constraint_counts[3] = { 0, 0, 0 };

static void clear_face_constraints(void) {
    for (int axis = 0; axis < 3; axis++) {
        face_constraint_counts[axis] = 0;
    }
}

static void add_face_constraint(int axis, int a, int b, float strain_c, float strain_t) {
    if (axis < 0 || axis >= 3) return;
    if (face_constraint_counts[axis] >= MAX_VOXELS) return;
    FaceConstraint *dst = &face_constraints[axis][face_constraint_counts[axis]++];
    dst->a = a;
    dst->b = b;
    dst->strain_compressive = strain_c;
    dst->strain_tensile = strain_t;
    dst->active = true;
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
static Vector3 v_add(Vector3 a, Vector3 b) {
    return (Vector3){ a.x + b.x, a.y + b.y, a.z + b.z };
}
static Vector3 v_mul(Vector3 v, float s) {
    return (Vector3){ v.x*s, v.y*s, v.z*s };
}
static Vector3 v_sub(Vector3 a, Vector3 b) {
    return (Vector3){ a.x - b.x, a.y - b.y, a.z - b.z };
}
static float v_dot(Vector3 a, Vector3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static Vector3 v_cross(Vector3 a, Vector3 b) {
    return (Vector3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

static Vector3 project_onto(Vector3 base, Vector3 vec) {
    float denom = v_dot(base, base);
    if (denom <= 1e-6f) return (Vector3){0.0f, 0.0f, 0.0f};
    float scale = v_dot(vec, base) / denom;
    return v_mul(base, scale);
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

static Vector3 v_normalize(Vector3 v) {
    float len = v_length(v);
    if (len <= 1e-6f) return (Vector3){0.0f, 0.0f, 0.0f};
    return v_mul(v, 1.0f / len);
}

static bool voxel_is_removed(const Voxel *v) {
    return (v->pos.x <= -900.0f && v->pos.y <= -900.0f && v->pos.z <= -900.0f);
}

static bool voxel_is_active(const Voxel *v) {
    return !voxel_is_removed(v);
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
static bool  voxelModelDirty = false;
static bool  voxelDebugLog = false;
static int   voxelDebugFrame = 0;
static float timeScale = 0.01f;



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


// Check occupancy
static bool occupied(int x,int y,int z){ return  table_get(x,y,z)>=0; }

static void rebuild_spatial_hash(void) {
    memset(table, 0, sizeof(table));
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_active(v)) continue;
        v->gx = (int)floorf(v->pos.x / VOXEL_SIZE);
        v->gy = (int)floorf(v->pos.y / VOXEL_SIZE);
        v->gz = (int)floorf(v->pos.z / VOXEL_SIZE);
        table_set(v->gx, v->gy, v->gz, i);
    }
}

static void mark_surface(int idx) {
    Voxel *v = &voxels[idx];
    int x = v->gx, y = v->gy, z = v->gz;
    v->surface[0] = !occupied(x+1, y,   z  );
    v->surface[1] = !occupied(x-1, y,   z  );
    v->surface[2] = !occupied(x,   y+1, z  );
    v->surface[3] = !occupied(x,   y-1, z  );
    v->surface[4] = !occupied(x,   y,   z+1);
    v->surface[5] = !occupied(x,   y,   z-1);
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

static void voxel_update_inverse_masses(Voxel *v) {
    float inv = (v->simulate && !v->fixed) ? 1.0f : 0.0f;
    for (int i = 0; i < 8; i++) {
        v->inv_mass[i] = inv;
    }
}

static void voxel_sync_corners_to_center(Voxel *v) {
    for (int i = 0; i < 8; i++) {
        v->corners[i] = v_add(v->pos, voxel_corner_offsets[i]);
        v->corner_vel[i] = v->vel;
    }
}

static void voxel_apply_translation(Voxel *v, Vector3 delta) {
    for (int i = 0; i < 8; i++) {
        v->corners[i] = v_add(v->corners[i], delta);
    }
}

static void voxel_store_prev_corners(Voxel *v) {
    for (int i = 0; i < 8; i++) {
        v->corners_prev[i] = v->corners[i];
    }
}

static Vector3 voxel_centroid_from_corners(const Vector3 corners[8]) {
    Vector3 sum = {0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 8; i++) sum = v_add(sum, corners[i]);
    return v_mul(sum, 1.0f / 8.0f);
}

static void voxel_update_center_from_corners(Voxel *v) {
    v->pos = voxel_centroid_from_corners(v->corners);
}

static void voxel_update_corner_velocities(Voxel *v, float inv_dt) {
    if (!isfinite(inv_dt)) return;
    for (int i = 0; i < 8; i++) {
        Vector3 delta = v_sub(v->corners[i], v->corners_prev[i]);
        v->corner_vel[i] = v_mul(delta, inv_dt);
    }
}

static void voxel_update_linear_velocity(Voxel *v, float inv_dt) {
    if (!isfinite(inv_dt)) return;
    Vector3 prev_c = voxel_centroid_from_corners(v->corners_prev);
    Vector3 curr_c = v->pos;
    Vector3 delta = v_sub(curr_c, prev_c);
    v->vel = v_mul(delta, inv_dt);
}

static void voxel_deactivate(Voxel *v) {
    v->simulate = false;
    v->fixed = true;
    voxel_update_inverse_masses(v);
    v->vel = (Vector3){0.0f, 0.0f, 0.0f};
    v->pos = (Vector3){ -999.0f, -999.0f, -999.0f };
    for (int i = 0; i < 8; i++) {
        v->corners[i] = v->pos;
        v->corner_vel[i] = (Vector3){0.0f, 0.0f, 0.0f};
        v->corners_prev[i] = v->pos;
    }
    v->distance_to_surface = 0.0f;
    v->skeleton = false;
}

static void voxel_set_velocity(Voxel *v, Vector3 vel) {
    v->vel = vel;
    for (int i = 0; i < 8; i++) v->corner_vel[i] = vel;
}

static void compute_voxel_distance_transform(void) {
    if (voxel_count <= 0) return;

    float *dist = (float *)malloc(sizeof(float) * (size_t)voxel_count);
    int   *queue = (int *)malloc(sizeof(int) * (size_t)voxel_count);
    if (!dist || !queue) {
        free(dist);
        free(queue);
        return;
    }

    int head = 0, tail = 0;
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_active(v)) {
            dist[i] = 0.0f;
            continue;
        }
        bool isSurface = false;
        for (int s = 0; s < 6; s++) {
            if (v->surface[s]) { isSurface = true; break; }
        }
        if (isSurface) {
            dist[i] = 0.0f;
            queue[tail++] = i;
        } else {
            dist[i] = FLT_MAX;
        }
    }

    if (tail == 0) {
        for (int i = 0; i < voxel_count; i++) {
            if (!voxel_is_active(&voxels[i])) continue;
            dist[i] = 0.0f;
            queue[tail++] = i;
        }
    }

    const int offs[6][3] = {
        { 1, 0, 0 }, { -1, 0, 0 },
        { 0, 1, 0 }, { 0,-1, 0 },
        { 0, 0, 1 }, { 0, 0,-1 }
    };

    while (head < tail) {
        int idx = queue[head++];
        float base = dist[idx];
        Voxel *v = &voxels[idx];
        if (!voxel_is_active(v)) continue;
        for (int n = 0; n < 6; n++) {
            int nx = v->gx + offs[n][0];
            int ny = v->gy + offs[n][1];
            int nz = v->gz + offs[n][2];
            int neighbor = table_get(nx, ny, nz);
            if (neighbor < 0) continue;
            if (dist[neighbor] > base + 1.0f) {
                dist[neighbor] = base + 1.0f;
                queue[tail++] = neighbor;
            }
        }
    }

    for (int i = 0; i < voxel_count; i++) {
        float voxDist = (dist[i] == FLT_MAX) ? 0.0f : dist[i] * VOXEL_SIZE;
        voxels[i].distance_to_surface = voxDist;
    }

    free(dist);
    free(queue);
}

static void compute_voxel_skeleton_from_distance(void) {
    const int offs[6][3] = {
        { 1, 0, 0 }, { -1, 0, 0 },
        { 0, 1, 0 }, { 0,-1, 0 },
        { 0, 0, 1 }, { 0, 0,-1 }
    };

    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_active(v)) {
            v->skeleton = false;
            continue;
        }
        float selfDist = v->distance_to_surface;
        bool localMax = true;
        for (int n = 0; n < 6; n++) {
            int nx = v->gx + offs[n][0];
            int ny = v->gy + offs[n][1];
            int nz = v->gz + offs[n][2];
            int neighbor = table_get(nx, ny, nz);
            if (neighbor < 0) continue;
            if (voxels[neighbor].distance_to_surface > selfDist + 1e-4f) {
                localMax = false;
                break;
            }
        }
        v->skeleton = (selfDist >= 2.0f * VOXEL_SIZE) && localMax;
    }
}

static void rebuild_face_constraints(void) {
    clear_face_constraints();

    const int offs[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 }
    };

    for (int i = 0; i < voxel_count; i++) {
        Voxel *a = &voxels[i];
        if (!voxel_is_active(a)) continue;
        for (int axis = 0; axis < 3; axis++) {
            int nx = a->gx + offs[axis][0];
            int ny = a->gy + offs[axis][1];
            int nz = a->gz + offs[axis][2];
            int neighbor = table_get(nx, ny, nz);
            if (neighbor < 0) continue;
            Voxel *b = &voxels[neighbor];
            if (!voxel_is_active(b)) continue;

            float da = a->distance_to_surface;
            float db = b->distance_to_surface;
            float diff = fabsf(da - db);

            float minDepth = fminf(da, db);
            float depthFactor = clampf(minDepth / (VOXEL_SIZE * 4.0f), 0.0f, 1.0f);
            float surfaceSoftness = 1.0f - clampf(minDepth / (VOXEL_SIZE * 1.25f), 0.0f, 1.0f);

            float strain_c = STRAIN_LIMIT_BASE_COMPRESSION - depthFactor * 0.3f;
            float strain_t = STRAIN_LIMIT_BASE_TENSION + depthFactor * 0.6f;

            strain_c *= (0.7f + 0.3f * (1.0f - surfaceSoftness));
            strain_t *= (0.4f + 0.6f * (1.0f - surfaceSoftness));

            if (axis == 1) { // vertical layers tend to delaminate
                strain_c *= 0.8f;
                strain_t *= 0.75f;
            } else if (axis == 2) { // Z-axis slightly stronger
                strain_c *= 1.1f;
                strain_t *= 1.1f;
            }

            if (a->skeleton && b->skeleton && diff < VOXEL_SIZE * 0.5f) {
                strain_c = -STRAIN_LIMIT_SHELL_BONUS;
                strain_t = STRAIN_LIMIT_SHELL_BONUS;
            } else if (diff > VOXEL_SIZE * 1.5f) {
                strain_c *= 0.5f;
                strain_t *= 0.5f;
            }

            add_face_constraint(axis, i, neighbor, strain_c, strain_t);
        }
    }
}

static void build_voxel_model_metadata(void) {
    rebuild_spatial_hash();
    for (int i = 0; i < voxel_count; i++) {
        if (voxel_is_active(&voxels[i])) mark_surface(i);
    }
    compute_voxel_distance_transform();
    compute_voxel_skeleton_from_distance();
    rebuild_face_constraints();
    voxelModelDirty = false;
}

static void voxel_apply_vgs(Voxel *v, float alpha, float beta, int iterations) {
    if (!v->simulate) return;

    Vector3 p[8];
    for (int i = 0; i < 8; i++) p[i] = v->corners[i];

    float r = v->particle_radius;
    float V0 = v->rest_volume;

    for (int it = 0; it < iterations; it++) {
        Vector3 c = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < 8; i++) c = v_add(c, p[i]);
        c = v_mul(c, 1.0f / 8.0f);

        Vector3 v0 = v_mul(v_add(v_add(v_add(v_sub(p[1], p[0]), v_sub(p[3], p[2])), v_sub(p[5], p[4])), v_sub(p[7], p[6])), 0.25f);
        Vector3 v1 = v_mul(v_add(v_add(v_add(v_sub(p[2], p[0]), v_sub(p[3], p[1])), v_sub(p[6], p[4])), v_sub(p[7], p[5])), 0.25f);
        Vector3 v2 = v_mul(v_add(v_add(v_add(v_sub(p[4], p[0]), v_sub(p[5], p[1])), v_sub(p[6], p[2])), v_sub(p[7], p[3])), 0.25f);

        Vector3 u0 = v_sub(v0, v_mul(v_add(project_onto(v1, v0), project_onto(v2, v0)), alpha));
        Vector3 u1 = v_sub(v1, v_mul(v_add(project_onto(v2, v1), project_onto(v0, v1)), alpha));
        Vector3 u2 = v_sub(v2, v_mul(v_add(project_onto(v0, v2), project_onto(v1, v2)), alpha));

        float len0 = v_length(u0);
        float len1 = v_length(u1);
        float len2 = v_length(u2);

        float target0 = (1.0f - beta) * r + beta * v_length(v0);
        float target1 = (1.0f - beta) * r + beta * v_length(v1);
        float target2 = (1.0f - beta) * r + beta * v_length(v2);

        if (len0 > 1e-6f) u0 = v_mul(u0, target0 / len0);
        if (len1 > 1e-6f) u1 = v_mul(u1, target1 / len1);
        if (len2 > 1e-6f) u2 = v_mul(u2, target2 / len2);

        float volume = v_dot(v_cross(u0, u1), u2);
        if (fabsf(volume) > 1e-6f) {
            float ratio = V0 / volume;
            ratio = clampf(ratio, -100.0f, 100.0f);
            float scaleMag = powf(fabsf(ratio), 1.0f / 3.0f);
            scaleMag = clampf(scaleMag, 0.25f, 4.0f);
            float scale = 0.5f * scaleMag;
            if (volume < 0.0f) scale = -scale;
            u0 = v_mul(u0, scale);
            u1 = v_mul(u1, scale);
            u2 = v_mul(u2, scale);
        } else {
            Vector3 center = c;
            for (int vid = 0; vid < 8; vid++) {
                p[vid] = v_add(center, voxel_corner_offsets[vid]);
            }
            break;
        }

        for (int vid = 0; vid < 8; vid++) {
            if (v->inv_mass[vid] <= 0.0f) continue;
            Vector3 offset = {0.0f, 0.0f, 0.0f};
            offset = v_add(offset, v_mul(u0, (vid & 1) ? 1.0f : -1.0f));
            offset = v_add(offset, v_mul(u1, (vid & 2) ? 1.0f : -1.0f));
            offset = v_add(offset, v_mul(u2, (vid & 4) ? 1.0f : -1.0f));
            p[vid] = v_add(c, offset);
        }
    }

    for (int i = 0; i < 8; i++) v->corners[i] = p[i];
}

static void solve_distance_constraint(Vector3 *pa, Vector3 *pb, float wa, float wb, float rest_length) {
    float wsum = wa + wb;
    if (wsum <= 0.0f) return;
    Vector3 diff = v_sub(*pb, *pa);
    float len = v_length(diff);
    if (len <= 1e-6f) return;
    Vector3 n = v_mul(diff, 1.0f / len);
    float constraint = len - rest_length;
    float lambda = -constraint / wsum;
    if (wa > 0.0f) *pa = v_add(*pa, v_mul(n, wa * lambda));
    if (wb > 0.0f) *pb = v_add(*pb, v_mul(n, -wb * lambda));
}

static void solve_face_constraints_axis(int axis) {
    const int (*pairs_surface)[2] = face_vertex_pairs[axis];
    const int (*axis_edges)[2] = axis_edge_pairs[axis];

    int count = face_constraint_counts[axis];
    for (int idx = 0; idx < count; idx++) {
        FaceConstraint *c = &face_constraints[axis][idx];
        if (!c->active) continue;
        if (c->a < 0 || c->a >= voxel_count || c->b < 0 || c->b >= voxel_count) {
            c->active = false;
            continue;
        }
        Voxel *a = &voxels[c->a];
        Voxel *b = &voxels[c->b];

        float maxStrain = -FLT_MAX;
        float minStrain = FLT_MAX;

        for (int e = 0; e < 4; e++) {
            int va0 = axis_edges[e][0];
            int va1 = axis_edges[e][1];
            Vector3 edgeA = v_sub(a->corners[va0], a->corners[va1]);
            float lenA = v_length(edgeA);
            float strainA = (lenA - VOXEL_SIZE) / VOXEL_SIZE;
            if (strainA > maxStrain) maxStrain = strainA;
            if (strainA < minStrain) minStrain = strainA;

            Vector3 edgeB = v_sub(b->corners[va0], b->corners[va1]);
            float lenB = v_length(edgeB);
            float strainB = (lenB - VOXEL_SIZE) / VOXEL_SIZE;
            if (strainB > maxStrain) maxStrain = strainB;
            if (strainB < minStrain) minStrain = strainB;
        }

        if (maxStrain > c->strain_tensile || minStrain < c->strain_compressive) {
            c->active = false;
            continue;
        }

        for (int p = 0; p < 4; p++) {
            int ia = pairs_surface[p][0];
            int ib = pairs_surface[p][1];
            solve_distance_constraint(&a->corners[ia], &b->corners[ib], a->inv_mass[ia], b->inv_mass[ib], 0.0f);
        }
    }
}

static void apply_voxel_constraints(float dt) {
    if (dt <= 0.0f) return;
    for (int iter = 0; iter < CONSTRAINT_ITERATIONS; iter++) {
        for (int i = 0; i < voxel_count; i++) {
            float alpha = VGS_ALPHA_DEFAULT;
            float beta  = VGS_BETA_DEFAULT;
            if (voxels[i].simulate) {
                alpha = 0.25f;
                beta  = 0.8f;
            }
            voxel_apply_vgs(&voxels[i], alpha, beta, VGS_LOCAL_ITERATIONS);
        }
        for (int axis = 0; axis < 3; axis++) {
            solve_face_constraints_axis(axis);
        }
    }

    float inv_dt = 1.0f / dt;
    for (int i = 0; i < voxel_count; i++) {
        if (!voxels[i].simulate) continue;
        voxel_update_center_from_corners(&voxels[i]);
        voxel_update_corner_velocities(&voxels[i], inv_dt);
        voxel_update_linear_velocity(&voxels[i], inv_dt);
        float speed = v_length(voxels[i].vel);
        if (speed > MAX_VOXEL_LINEAR_SPEED) {
            voxels[i].vel = v_mul(voxels[i].vel, MAX_VOXEL_LINEAR_SPEED / speed);
        }
    }
}


// Add a voxel (static or dynamic)
static int addVoxel(float px, float py, float pz, bool fixed, bool simulate, Color color, int type) {
    if (voxel_count >= MAX_VOXELS) return -1;
    int idx = voxel_count++;
    Voxel *v = &voxels[idx];
    v->pos = (Vector3){ px, py, pz };
    voxel_set_velocity(v, (Vector3){0.0f, 0.0f, 0.0f});
    v->fixed = fixed;
    v->simulate = simulate;
    v->color = color;
    v->type = type;
    v->owner   = -1;
    memset(v->surface, 0, sizeof v->surface);
    // compute grid coords
    v->gx = (int)floorf(px / VOXEL_SIZE);
    v->gy = (int)floorf(py / VOXEL_SIZE);
    v->gz = (int)floorf(pz / VOXEL_SIZE);
    v->rest_volume = VOXEL_SIZE * VOXEL_SIZE * VOXEL_SIZE;
    v->particle_radius = VOXEL_PARTICLE_RADIUS;
    voxel_update_inverse_masses(v);
    voxel_sync_corners_to_center(v);
    voxel_store_prev_corners(v);
    v->distance_to_surface = 0.0f;
    v->skeleton = false;
    table_set(v->gx, v->gy, v->gz, idx);
    return idx;
}

// Build static demo cube of voxels
static void buildDemo(void) {
    const int floorExtent = 55; // covers player spawn range (≈22 units)
    for (int gx = -floorExtent; gx <= floorExtent; gx++) {
        for (int gz = -floorExtent; gz <= floorExtent; gz++) {
            if (abs(gx) > 45 && abs(gz) > 45) continue; // trim corners for density reduction
            float px = (gx + 0.5f) * VOXEL_SIZE;
            float pz = (gz + 0.5f) * VOXEL_SIZE;
            addVoxel(px, 0.0f, pz, true, false, (Color){ 140, 140, 150, 255 }, 0);
        }
    }

    const int pillarHeight = 12;
    const int pillarRadius = 1;
    const int pillarPositions[4][2] = {
        { -20, -20 },
        { -20,  20 },
        {  20, -20 },
        {  20,  20 }
    };
    for (int p = 0; p < 4; p++) {
        int cx = pillarPositions[p][0];
        int cz = pillarPositions[p][1];
        for (int gy = 1; gy <= pillarHeight; gy++) {
            for (int dx = -pillarRadius; dx <= pillarRadius; dx++) {
                for (int dz = -pillarRadius; dz <= pillarRadius; dz++) {
                    if (abs(dx) + abs(dz) > pillarRadius + 1) continue;
                    float px = (cx + dx + 0.5f) * VOXEL_SIZE;
                    float py = (gy + 0.5f) * VOXEL_SIZE;
                    float pz = (cz + dz + 0.5f) * VOXEL_SIZE;
                    addVoxel(px, py, pz, true, false, (Color){ 165, 115, 90, 255 }, 0);
                }
            }
        }
    }
    
    const int blobRadius = 5;
    const int blobHeight = 10;
    Vector3 blobCenter = { 0.0f, (blobHeight + 4) * VOXEL_SIZE, 0.0f };
    for (int gx = -blobRadius; gx <= blobRadius; gx++) {
        for (int gy = 0; gy < blobHeight; gy++) {
            for (int gz = -blobRadius; gz <= blobRadius; gz++) {
                float dx = (float)gx;
                float dz = (float)gz;
                if (dx*dx + dz*dz > (blobRadius + 0.5f)*(blobRadius + 0.5f)) continue;
                float px = blobCenter.x + (gx + 0.5f) * VOXEL_SIZE;
                float py = blobCenter.y + (gy + 0.5f) * VOXEL_SIZE;
                float pz = blobCenter.z + (gz + 0.5f) * VOXEL_SIZE;
                Color c = (Color){ 80, (unsigned char)(130 + gy*4), 200, 255 };
                addVoxel(px, py, pz, false, true, c, 0);
            }
        }
    }
}



static int first_voxel_hit(Ray ray, float t_max, int ignore_id);
static void UpdateKdRatio(int player_index);

// Reset game: players and voxels
static void ResetGame(void) {
    // init players
    for (int i = 0; i < 2; i++) {
        players[i].pos = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
        players[i].yaw = (i == 0) ? 0 : 180;
        players[i].pitch = 0;
        players[i].yaw_vel = 0;
        players[i].pitch_vel = 0;
        players[i].vel = (Vector3){0,0,0};
        players[i].onGround = true;
        players[i].vType = 1;
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
    build_voxel_model_metadata();
    meshDirty = true;

} 

static void UpdateKdRatio(int player_index) {
    Player *p = &players[player_index];
    p->kd_ratio = (float)(p->kills + 1) / (p->deaths + 1);
}

// Physics step for voxels
static void physics_step(float dt) {    // Rebuild spatial hash
    rebuild_spatial_hash();
    // Simulate dynamic voxels
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        voxel_store_prev_corners(v);

        // Apply gravity
        v->vel.y -= GRAVITY * dt;

        // Continuous collision detection
        Vector3 displacement = v_mul(v->vel, dt);
        float distance = v_length(displacement);
        bool hit_voxel = false;

        if (distance > 0.0001f) { // only cast if moving
            Ray ray = { v->pos, v_norm(v->vel) };
            int hit_id = first_voxel_hit(ray, distance, i);

            if (hit_id >= 0) {
                hit_voxel = true;

                // Stop the bullet
                voxel_deactivate(v);

                int brushExtent = (voxelBrushSpan < 1) ? 1 : voxelBrushSpan;

                if (v->type == 1) { // DESTRUCTION
                    Voxel *u = &voxels[hit_id];
                    int anchorX = u->gx;
                    int anchorY = u->gy;
                    int anchorZ = u->gz;

                    for (int dx = 0; dx < brushExtent; dx++) {
                        for (int dy = 0; dy < brushExtent; dy++) {
                            for (int dz = 0; dz < brushExtent; dz++) {
                                int victim_idx = table_get(anchorX + dx, anchorY + dy, anchorZ + dz);
                                if (victim_idx >= 0) {
                                    Voxel *victim = &voxels[victim_idx];
                                    table_remove(victim->gx, victim->gy, victim->gz);
                                    mark_surface_neighbors(victim->pos);
                                    voxel_deactivate(victim);
                                    voxelModelDirty = true;
                                }
                            }
                        }
                    }
                } else { // CONSTRUCTION
                    int anchorX = v->gx;
                    int anchorY = v->gy;
                    int anchorZ = v->gz;

                    for (int dx = 0; dx < brushExtent; dx++) {
                        for (int dy = 0; dy < brushExtent; dy++) {
                            for (int dz = 0; dz < brushExtent; dz++) {
                                int targetX = anchorX + dx;
                                int targetY = anchorY + dy;
                                int targetZ = anchorZ + dz;
                                if (!occupied(targetX, targetY, targetZ)) {
                                    float px = (targetX + 0.5f) * VOXEL_SIZE;
                                    float py = (targetY + 0.5f) * VOXEL_SIZE;
                                    float pz = (targetZ + 0.5f) * VOXEL_SIZE;
                                    int new_idx = addVoxel(px, py, pz, true, false, v->color, 0);
                                    if (new_idx >= 0) {
                                        mark_surface(new_idx);
                                        mark_surface_neighbors(voxels[new_idx].pos);
                                        voxelModelDirty = true;
                                    }
                                }
                            }
                        }
                    }
                }
                //reset mesh
                meshDirty = true;
                voxelModelDirty = true;
            }
        }

        if (!hit_voxel) {
            // Move
            v->pos = v_add(v->pos, displacement);
            if (fabsf(v->pos.x) > FLOOR_SIZE * 8.0f ||
                fabsf(v->pos.y) > FLOOR_SIZE * 8.0f ||
                fabsf(v->pos.z) > FLOOR_SIZE * 8.0f) {
                voxel_deactivate(v);
            } else {
                voxel_apply_translation(v, displacement);
            }
            for (int j = 0; j < 2; j++) {
                float dx = v->pos.x - players[j].pos.x;
                float dy = v->pos.y - players[j].pos.y;
                float dz = v->pos.z - players[j].pos.z;
                if (fabsf(dx) < PLAYER_SIZE && fabsf(dy) < PLAYER_SIZE && fabsf(dz) < PLAYER_SIZE) {
                    players[j].last_damage_time = (float)GetTime();
                    if (players[j].shield > 0) {
                        players[j].shield -= VOXEL_DAMAGE;
                        if (players[j].shield < 0) {
                        players[j].health += players[j].shield;
                        players[j].shield = 0;
                    }
                } else {
                    players[j].health -= VOXEL_DAMAGE;
                }

                    voxel_deactivate(v);

                    if (players[j].health <= 0) {
                        if (v->owner >= 0 && v->owner != j) {
                            players[v->owner].kills++;
                            players[j].deaths++;
                            UpdateKdRatio(v->owner);
                            UpdateKdRatio(j);
                        }
                        players[j].pos     = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
                        players[j].vel     = (Vector3){0,0,0};
                        players[j].onGround= true;
                        players[j].yaw     = (j == 0 ? 0 : 180);
                        players[j].pitch   = 0;
                        players[j].health = BASE_HEALTH / players[j].kd_ratio;
                        players[j].shield = BASE_SHIELD;
                    }
                    break;
                }
            }
        }
    }

    apply_voxel_constraints(dt);

    if (voxelModelDirty) {
        build_voxel_model_metadata();
    }

    if (voxelDebugLog) {
        voxelDebugFrame++;
        if (voxelDebugFrame % 5 == 0) {
            int logged = 0;
            for (int i = 0; i < voxel_count && logged < 8; i++) {
                Voxel *v = &voxels[i];
                if (!v->simulate) continue;
                printf("[voxel-debug] id=%d pos=(%.2f,%.2f,%.2f) vel=(%.2f,%.2f,%.2f) dist=%.2f skeleton=%d\n",
                       i,
                       v->pos.x, v->pos.y, v->pos.z,
                       v->vel.x, v->vel.y, v->vel.z,
                       v->distance_to_surface,
                       v->skeleton);
                logged++;
            }
            if (logged == 0) {
                printf("[voxel-debug] no active voxels\n");
            }
        }
    }
}

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
        voxel_set_velocity(&voxels[vix], v_mul(dir, 50.0f));
        voxels[vix].owner = idx;
    }
}
// Append the 12 edges (24 vertices) of a cube to the current RL_LINES batch
static void drawCubeEdges(Vector3 pos, float w, float h, float d)
{
    float hw = w * 0.5f, hh = h * 0.5f, hd = d * 0.5f;
    float x = pos.x, y = pos.y, z = pos.z;

    // 4 bottom vertices
    Vector3 v000 = { x - hw, y - hh, z - hd };
    Vector3 v001 = { x - hw, y - hh, z + hd };
    Vector3 v101 = { x + hw, y - hh, z + hd };
    Vector3 v100 = { x + hw, y - hh, z - hd };

    // 4 top vertices
    Vector3 v010 = { x - hw, y + hh, z - hd };
    Vector3 v011 = { x - hw, y + hh, z + hd };
    Vector3 v111 = { x + hw, y + hh, z + hd };
    Vector3 v110 = { x + hw, y + hh, z - hd };

    // Bottom rectangle
    rlVertex3f(v000.x, v000.y, v000.z); rlVertex3f(v001.x, v001.y, v001.z);
    rlVertex3f(v001.x, v001.y, v001.z); rlVertex3f(v101.x, v101.y, v101.z);
    rlVertex3f(v101.x, v101.y, v101.z); rlVertex3f(v100.x, v100.y, v100.z);
    rlVertex3f(v100.x, v100.y, v100.z); rlVertex3f(v000.x, v000.y, v000.z);

    // Top rectangle
    rlVertex3f(v010.x, v010.y, v010.z); rlVertex3f(v011.x, v011.y, v011.z);
    rlVertex3f(v011.x, v011.y, v011.z); rlVertex3f(v111.x, v111.y, v111.z);
    rlVertex3f(v111.x, v111.y, v111.z); rlVertex3f(v110.x, v110.y, v110.z);
    rlVertex3f(v110.x, v110.y, v110.z); rlVertex3f(v010.x, v010.y, v010.z);

    // Vertical pillars
    rlVertex3f(v000.x, v000.y, v000.z); rlVertex3f(v010.x, v010.y, v010.z);
    rlVertex3f(v001.x, v001.y, v001.z); rlVertex3f(v011.x, v011.y, v011.z);
    rlVertex3f(v101.x, v101.y, v101.z); rlVertex3f(v111.x, v111.y, v111.z);
    rlVertex3f(v100.x, v100.y, v100.z); rlVertex3f(v110.x, v110.y, v110.z);
}

static void drawCubeMan(Vector3 pos,
                        float width,
                        float height,
                        float depth,
                        Color color)
{
    float hw = width  * 0.5f;
    float hh = height * 0.5f;
    float hd = depth  * 0.5f;

    float x = pos.x, y = pos.y, z = pos.z;

    rlColor4ub(color.r, color.g, color.b, color.a);

    /* ----------  FRONT  (+Z)  ---------- */
    rlNormal3f(0.0f, 0.0f, 1.0f);
    rlVertex3f(x - hw, y - hh, z + hd);   // CCW
    rlVertex3f(x + hw, y - hh, z + hd);
    rlVertex3f(x + hw, y + hh, z + hd);

    rlVertex3f(x - hw, y - hh, z + hd);
    rlVertex3f(x + hw, y + hh, z + hd);
    rlVertex3f(x - hw, y + hh, z + hd);

    /* ----------  BACK  (-Z)  ----------- */
    rlNormal3f(0.0f, 0.0f, -1.0f);
    rlVertex3f(x + hw, y - hh, z - hd);
    rlVertex3f(x - hw, y - hh, z - hd);
    rlVertex3f(x - hw, y + hh, z - hd);

    rlVertex3f(x + hw, y - hh, z - hd);
    rlVertex3f(x - hw, y + hh, z - hd);
    rlVertex3f(x + hw, y + hh, z - hd);

    /* ----------  TOP   (+Y)  ----------- */
    rlNormal3f(0.0f, 1.0f, 0.0f);
    rlVertex3f(x - hw, y + hh, z + hd);
    rlVertex3f(x + hw, y + hh, z + hd);
    rlVertex3f(x + hw, y + hh, z - hd);

    rlVertex3f(x - hw, y + hh, z + hd);
    rlVertex3f(x + hw, y + hh, z - hd);
    rlVertex3f(x - hw, y + hh, z - hd);

    /* ----------  BOTTOM (-Y) ----------- */
    rlNormal3f(0.0f, -1.0f, 0.0f);
    rlVertex3f(x - hw, y - hh, z + hd);
    rlVertex3f(x + hw, y - hh, z - hd);
    rlVertex3f(x + hw, y - hh, z + hd);

    rlVertex3f(x - hw, y - hh, z + hd);
    rlVertex3f(x - hw, y - hh, z - hd);
    rlVertex3f(x + hw, y - hh, z - hd);

    /* ----------  RIGHT (+X) ------------ */
    rlNormal3f(1.0f, 0.0f, 0.0f);
    rlVertex3f(x + hw, y - hh, z + hd);
    rlVertex3f(x + hw, y - hh, z - hd);
    rlVertex3f(x + hw, y + hh, z - hd);

    rlVertex3f(x + hw, y - hh, z + hd);
    rlVertex3f(x + hw, y + hh, z - hd);
    rlVertex3f(x + hw, y + hh, z + hd);

    /* ----------  LEFT  (-X) ------------ */
    rlNormal3f(-1.0f, 0.0f, 0.0f);
    rlVertex3f(x - hw, y - hh, z - hd);
    rlVertex3f(x - hw, y - hh, z + hd);
    rlVertex3f(x - hw, y + hh, z + hd);

    rlVertex3f(x - hw, y - hh, z - hd);
    rlVertex3f(x - hw, y + hh, z + hd);
    rlVertex3f(x - hw, y + hh, z - hd);
}

static void draw_deformed_voxel(const Voxel *v) {
    static const int faces[6][4] = {
        {4,5,7,6}, // +Z
        {0,2,3,1}, // -Z
        {2,6,7,3}, // +Y
        {0,1,5,4}, // -Y
        {1,3,7,5}, // +X
        {0,4,6,2}  // -X
    };

    rlColor4ub(v->color.r, v->color.g, v->color.b, v->color.a);
    if (fabsf(v->pos.x) > FLOOR_SIZE * 8.0f ||
        fabsf(v->pos.y) > FLOOR_SIZE * 8.0f ||
        fabsf(v->pos.z) > FLOOR_SIZE * 8.0f) {
        return;
    }
    for (int f = 0; f < 6; f++) {
        Vector3 a = v->corners[faces[f][0]];
        Vector3 b = v->corners[faces[f][1]];
        Vector3 c = v->corners[faces[f][2]];
        Vector3 d = v->corners[faces[f][3]];
        Vector3 normal = v_normalize(v_cross(v_sub(b, a), v_sub(c, a)));
        rlNormal3f(normal.x, normal.y, normal.z);
        rlVertex3f(a.x, a.y, a.z);
        rlVertex3f(b.x, b.y, b.z);
        rlVertex3f(c.x, c.y, c.z);

        rlVertex3f(a.x, a.y, a.z);
        rlVertex3f(c.x, c.y, c.z);
        rlVertex3f(d.x, d.y, d.z);
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
    (void)cam;
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
    rlEnableBackfaceCulling();
    //draw moving voxels
    rlBegin(RL_TRIANGLES);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (v->simulate){
            draw_deformed_voxel(v);
        }

    }
    rlEnd();
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
    // init window and render textures
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Split-Screen FPS (raylib)");
    SetTargetFPS(120);
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
            case GAME_STATE_PLAYING: {
        if (IsKeyPressed(KEY_ZERO)) {
            timeScale = 1.0f;
            printf("[time] scale reset to 1.0\n");
        }
        if (IsKeyPressed(KEY_MINUS)) {
            timeScale = fmaxf(0.1f, timeScale - 0.1f);
            printf("[time] scale = %.2f\n", timeScale);
        }
        if (IsKeyPressed(KEY_EQUAL)) {
            timeScale = fminf(4.0f, timeScale + 0.1f);
            printf("[time] scale = %.2f\n", timeScale);
        }

        float dt_raw = GetFrameTime();
        float dt = dt_raw * timeScale;
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
        if (IsKeyPressed(KEY_F4)) {
            voxelDebugLog = !voxelDebugLog;
            printf("[voxel-debug] logging %s\n", voxelDebugLog ? "enabled" : "disabled");
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
                HandleKeyboardInput(i, dt_raw);
            } else if (playerInput[i] == INPUT_TYPE_GAMEPAD) {
                HandleGamepadInput(i, dt_raw);
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
                get_adjacent_voxel_directions(v_add(p->pos,v_mul(p->vel,dt_raw)), neigh);
                // X-axis (+X/neigh[0], -X/neigh[1])
                if ((p->vel.x > 0 && neigh[0]) || (p->vel.x < 0 && neigh[1])) p->vel.x = 0;
                // Y-axis (+Y/neigh[2], -Y/neigh[3])
                if ((p->vel.y > 0 && neigh[2]) || (p->vel.y < 0 && neigh[3])) p->vel.y = 0;
                // Z-axis (+Z/neigh[4], -Z/neigh[5])
                if ((p->vel.z > 0 && neigh[4]) || (p->vel.z < 0 && neigh[5])) p->vel.z = 0;
                if (!neigh[3]){
                    // apply gravity
                    p->vel.y -= GRAVITY*dt_raw;
                    p->onGround = false;
                }else{
                    p->onGround = true;
                }
            }
            

            // apply movement
            p->pos.x += p->vel.x*dt_raw;
            p->pos.y += p->vel.y*dt_raw;
            p->pos.z += p->vel.z*dt_raw;

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
        // update voxel physics
        int subStep = 3;
        for( int i = 0; i < subStep; i++){
            physics_step(dt/subStep);
        }
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
        EndDrawing();
        break;
            }
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
