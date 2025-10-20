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
#define VGS_ALPHA 0.001f // 0 means no shear correction 1 means complete shear correction
#define VGS_BETA 0.001f // 0 means no tension correction 1 means total tension correction
#define VGS_ITERS 3
#define VGS_EPS 1e-6f
#define PBD_MAX_STEP_DT 0.005f
#define PBD_SUBSTEPS 3
#define PBD_CONSTRAINT_ITERS 5
#define COLLISION_RELAXATION 0.5f
#define CENTER_RELAXATION 0.9f
#define VELOCITY_DAMPING 0.99f

// KD-stats constants
#define BASE_HEALTH 100
#define BASE_SHIELD 150
#define SHIELD_REGEN_DELAY 5.0f
#define VOXEL_DAMAGE 50

// Voxel physics constants
#define MAX_VOXELS    131072
#define HASH_SIZE     131072    // must be power of two
#define VOXEL_SIZE     0.5f    // size of each voxel cube
#define MAX_PARTICLES (MAX_VOXELS * 8)
#define STRAIN_BREAK_THRESHOLD 1.0f
#define SHEAR_BREAK_THRESHOLD 0.5f

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
    int   refcount;
    bool  active;
    int   active_index;
} Particle;

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
    Particle *particles[8];
    float rest_volume;
    float rest_edge;
    float particle_radius;
    bool glued_faces[6];
    bool pending_full_break;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;

static Particle particles_pool[MAX_PARTICLES];
static Particle *active_particles[MAX_PARTICLES];
static int particle_pool_count = 0;
static int active_particle_count = 0;

static void reset_particle_pool(void) {
    particle_pool_count = 0;
    active_particle_count = 0;
    memset(particles_pool, 0, sizeof(particles_pool));
    memset(active_particles, 0, sizeof(active_particles));
}

static Particle *particle_create(Vector3 pos, float inv_mass) {
    if (particle_pool_count >= MAX_PARTICLES) {
        return NULL;
    }

    Particle *p = &particles_pool[particle_pool_count++];
    p->pos = pos;
    p->prev_pos = pos;
    p->predicted_pos = pos;
    p->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    p->inv_mass = inv_mass;
    p->refcount = 1;
    p->active = true;
    p->active_index = active_particle_count;
    active_particles[active_particle_count++] = p;
    return p;
}

static Particle *particle_clone(const Particle *src) {
    if (!src) {
        return NULL;
    }
    if (particle_pool_count >= MAX_PARTICLES) {
        return NULL;
    }

    Particle *p = &particles_pool[particle_pool_count++];
    *p = *src;
    p->refcount = 1;
    p->active = true;
    p->active_index = active_particle_count;
    active_particles[active_particle_count++] = p;
    return p;
}

static void particle_retain(Particle *p) {
    if (!p) {
        return;
    }
    p->refcount++;
}

static void particle_release(Particle *p) {
    if (!p || p->refcount <= 0) {
        return;
    }

    p->refcount--;
    if (p->refcount == 0 && p->active) {
        int idx = p->active_index;
        int last = --active_particle_count;
        if (last >= 0) {
            if (idx != last) {
                active_particles[idx] = active_particles[last];
                active_particles[idx]->active_index = idx;
            }
            active_particles[last] = NULL;
        }
        p->active = false;
        p->active_index = -1;
    }
}

static inline float voxel_particle_radius(const Voxel *v) {
    return (v->particle_radius > 0.0f) ? v->particle_radius : PARTICLE_RADIUS;
}

static bool debugDrawParticles = false;
static bool debugColorParticlesByVelocity = false;
static const float PARTICLE_DEBUG_MARKER_RADIUS = 0.6f;
static const float PARTICLE_DEBUG_MAX_SPEED = 20.0f;

static const int corner_signs[8][3] = {
    { -1, -1, -1 }, {  1, -1, -1 },
    { -1,  1, -1 }, {  1,  1, -1 },
    { -1, -1,  1 }, {  1, -1,  1 },
    { -1,  1,  1 }, {  1,  1,  1 }
};

static const int face_corner_indices[6][4] = {
    { 1, 3, 5, 7 }, // +X
    { 0, 2, 4, 6 }, // -X
    { 2, 3, 6, 7 }, // +Y
    { 0, 1, 4, 5 }, // -Y
    { 4, 5, 6, 7 }, // +Z
    { 0, 1, 2, 3 }  // -Z
};

static const int face_offsets[6][3] = {
    { 1, 0, 0 },
    { -1, 0, 0 },
    { 0, 1, 0 },
    { 0, -1, 0 },
    { 0, 0, 1 },
    { 0, 0, -1 }
};

static const int opposite_face[6] = { 1, 0, 3, 2, 5, 4 };

static int table_get(int x, int y, int z);


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

// Add a voxel (static or dynamic)
static int addVoxel(float px, float py, float pz, bool fixed, bool simulate, Color color, int type, float invm) {
    if (voxel_count >= MAX_VOXELS) return -1;
    int idx = voxel_count++;
    Voxel *v = &voxels[idx];
    v->pos = (Vector3){ px, py, pz };
    v->vel = (Vector3){ 0,0,0 };
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
    v->rest_edge = VOXEL_SIZE;
    v->particle_radius = PARTICLE_RADIUS;
    memset(v->glued_faces, 0, sizeof(v->glued_faces));
    v->pending_full_break = false;
    const float half = VOXEL_SIZE * 0.5f;
    for (int i = 0; i < 8; ++i) {
        Vector3 pos = {
            px + corner_signs[i][0] * half,
            py + corner_signs[i][1] * half,
            pz + corner_signs[i][2] * half
        };
        float inv_mass = (fixed || !simulate) ? 0.0f : invm;
        Particle *p = particle_create(pos, inv_mass);
        if (!p) {
            for (int j = 0; j < i; ++j) {
                particle_release(v->particles[j]);
                v->particles[j] = NULL;
            }
            voxel_count--;
            return -1;
        }
        v->particles[i] = p;
    }
    table_set(v->gx, v->gy, v->gz, idx);
    return idx;
}

static void glue_neighbor_faces(void) {
    static const int axis_offsets[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 }
    };
    static const int positive_face[3] = { 0, 2, 4 };
    static const int negative_face[3] = { 1, 3, 5 };

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *a = &voxels[i];
        for (int axis = 0; axis < 3; ++axis) {
            int nx = a->gx + axis_offsets[axis][0];
            int ny = a->gy + axis_offsets[axis][1];
            int nz = a->gz + axis_offsets[axis][2];

            int neighbor_idx = table_get(nx, ny, nz);
            if (neighbor_idx < 0) {
                continue;
            }

            Voxel *b = &voxels[neighbor_idx];
            int faceA = positive_face[axis];
            int faceB = negative_face[axis];

            for (int c = 0; c < 4; ++c) {
                int ia = face_corner_indices[faceA][c];
                int ib = face_corner_indices[faceB][c];

                Particle *shared = a->particles[ia];
                Particle *old = b->particles[ib];

                if (shared == old) {
                    continue;
                }

                particle_release(old);
                particle_retain(shared);
                b->particles[ib] = shared;
            }

            a->glued_faces[faceA] = true;
            b->glued_faces[faceB] = true;
        }
    }
}

static void detach_face_particles(Voxel *voxel, int face_index) {
    if (!voxel) {
        return;
    }

    for (int c = 0; c < 4; ++c) {
        int corner = face_corner_indices[face_index][c];
        Particle *p_old = voxel->particles[corner];
        if (!p_old || p_old->refcount <= 1) {
            continue;
        }

        Particle *p_new = particle_clone(p_old);
        if (!p_new) {
            continue;
        }

        voxel->particles[corner] = p_new;
        particle_release(p_old);
    }
}

static void break_face_link(Voxel *voxel, int face_index) {
    if (!voxel || !voxel->glued_faces[face_index]) {
        return;
    }

    int nx = voxel->gx + face_offsets[face_index][0];
    int ny = voxel->gy + face_offsets[face_index][1];
    int nz = voxel->gz + face_offsets[face_index][2];
    int neighbor_idx = table_get(nx, ny, nz);

    detach_face_particles(voxel, face_index);
    voxel->glued_faces[face_index] = false;

    if (neighbor_idx < 0) {
        return;
    }

    Voxel *neighbor = &voxels[neighbor_idx];
    int opposite = opposite_face[face_index];
    detach_face_particles(neighbor, opposite);
    neighbor->glued_faces[opposite] = false;
}

static void process_pending_breaks(void) {
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->pending_full_break) {
            continue;
        }

        for (int face = 0; face < 6; ++face) {
            break_face_link(voxel, face);
        }

        voxel->pending_full_break = false;
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

    // Central platform
    int platform_size = 2;
    int platform_height =2; // 15 / 3
    int platform_base_height = 1; // to keep top at same level (21)
    for (int y = platform_base_height; y <= platform_base_height + platform_height; y++) {
        for (int x = M/2 - platform_size/2; x <= M/2 + platform_size/2; x++) {
            for (int z = M/2 - platform_size/2; z <= M/2 + platform_size/2; z++) {
                float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                float py = (y + 0.5f) * VOXEL_SIZE;
                float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                if (y==platform_base_height){
                    addVoxel(px, py, pz, false, true, (Color){ 100, 200, 100, 255 }, 0, 0.0f);
                }else{
                addVoxel(px, py, pz, false, true, (Color){ 100, 200, 100, 255 }, 0, 1.0f);
                }
            }
        }
    }
    // float px1 = (M/2 + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    // float py1 = (0.5f + 0.5f) * VOXEL_SIZE;
    // float pz1 = (M/2 + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    // addVoxel(px1, py1, pz1, false, true, (Color){ 100, 200, 100, 255 }, 0, 0.0f);
    // float px2 = (M/2 + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    // float py2 = (1.5f + 0.5f) * VOXEL_SIZE;
    // float pz2 = (M/2 + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
    // addVoxel(px2, py2, pz2, false, true, (Color){ 100, 200, 100, 255 }, 0, 1.0f);
    glue_neighbor_faces();
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
    reset_particle_pool();
    // clear hash
    memset(table, 0, sizeof(table));
    // build static blocks
    buildDemo();
    for (int i = 0; i < voxel_count; i++) {
        mark_surface(i);
    }
    meshDirty = true;

} 

static void UpdateKdRatio(int player_index) {
    Player *p = &players[player_index];
    p->kd_ratio = (float)(p->kills + 1) / (p->deaths + 1);
}

// Physics step for voxels
static void physics_step(float dt) {    // Rebuild spatial hash
    memset(table, 0, sizeof(table));
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        int x = (int)floorf(v->pos.x / VOXEL_SIZE);
        int y = (int)floorf(v->pos.y / VOXEL_SIZE);
        int z = (int)floorf(v->pos.z / VOXEL_SIZE);
        v->gx = x; v->gy = y; v->gz = z;
        table_set(x, y, z, i);
    }
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

    for (int i = 0; i < active_particle_count; ++i) {
        Particle *p = active_particles[i];

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
        Particle *part = voxel->particles[i];
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

        float len_v0 = v_length(v0);
        float len_v1 = v_length(v1);
        float len_v2 = v_length(v2);

        if (!voxel->pending_full_break && voxel->rest_edge > 0.0f) {
            float strain_x = fabsf(len_v0 - voxel->rest_edge) / voxel->rest_edge;
            float strain_y = fabsf(len_v1 - voxel->rest_edge) / voxel->rest_edge;
            float strain_z = fabsf(len_v2 - voxel->rest_edge) / voxel->rest_edge;

            if ((voxel->glued_faces[0] || voxel->glued_faces[1]) && strain_x > STRAIN_BREAK_THRESHOLD) {
                voxel->pending_full_break = true;
            }
            if ((voxel->glued_faces[2] || voxel->glued_faces[3]) && strain_y > STRAIN_BREAK_THRESHOLD) {
                voxel->pending_full_break = true;
            }
            if ((voxel->glued_faces[4] || voxel->glued_faces[5]) && strain_z > STRAIN_BREAK_THRESHOLD) {
                voxel->pending_full_break = true;
            }

            if (!voxel->pending_full_break) {
                float inv_len0 = (len_v0 > VGS_EPS) ? 1.0f / len_v0 : 0.0f;
                float inv_len1 = (len_v1 > VGS_EPS) ? 1.0f / len_v1 : 0.0f;
                float inv_len2 = (len_v2 > VGS_EPS) ? 1.0f / len_v2 : 0.0f;

                float shear_xy = (inv_len0 > 0.0f && inv_len1 > 0.0f)
                    ? fabsf(v_dot(v0, v1)) * inv_len0 * inv_len1
                    : 0.0f;
                float shear_xz = (inv_len0 > 0.0f && inv_len2 > 0.0f)
                    ? fabsf(v_dot(v0, v2)) * inv_len0 * inv_len2
                    : 0.0f;
                float shear_yz = (inv_len1 > 0.0f && inv_len2 > 0.0f)
                    ? fabsf(v_dot(v1, v2)) * inv_len1 * inv_len2
                    : 0.0f;

                if ((voxel->glued_faces[4] || voxel->glued_faces[5]) && shear_xy > SHEAR_BREAK_THRESHOLD) {
                    voxel->pending_full_break = true;
                }
                if ((voxel->glued_faces[2] || voxel->glued_faces[3]) && shear_xz > SHEAR_BREAK_THRESHOLD) {
                    voxel->pending_full_break = true;
                }
                if ((voxel->glued_faces[0] || voxel->glued_faces[1]) && shear_yz > SHEAR_BREAK_THRESHOLD) {
                    voxel->pending_full_break = true;
                }
            }
        }

        Vector3 u0 = v_sub(v0, v_mul(v_add(vgs_project(v1, v0), vgs_project(v2, v0)), VGS_ALPHA));
        Vector3 u1 = v_sub(v1, v_mul(v_add(vgs_project(v2, v1), vgs_project(v0, v1)), VGS_ALPHA));
        Vector3 u2 = v_sub(v2, v_mul(v_add(vgs_project(v0, v2), vgs_project(v1, v2)), VGS_ALPHA));

        float len0 = v_length(u0);
        float len1 = v_length(u1);
        float len2 = v_length(u2);

        float target0 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * len_v0);
        float target1 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * len_v1);
        float target2 = ((1.0f - VGS_BETA) * rest_edge) + (VGS_BETA * len_v2);

        if (len0 > VGS_EPS) u0 = v_mul(u0, target0 / len0);
        if (len1 > VGS_EPS) u1 = v_mul(u1, target1 / len1);
        if (len2 > VGS_EPS) u2 = v_mul(u2, target2 / len2);

        // Volume correction mirrors the GPU "ResizeVoxelBasis" stage.
        float volume = v_dot(v_cross(u0, u1), u2);
        if (fabsf(volume) > VGS_EPS) {
            float scale = rest_volume / volume;
            float root = cbrtf(fabsf(scale));
            if (scale < 0.0f) {
                root = -root;
            }
            float factor = 0.5f * root;
            u0 = v_mul(u0, factor);
            u1 = v_mul(u1, factor);
            u2 = v_mul(u2, factor);
        }

        // Rebuild the voxel corners from the orthogonal frame and push dynamic particles only.
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
        voxel->particles[i]->predicted_pos = p[i];
    }
}

static void compute_voxel_center_and_mass(const Voxel *voxel, Vector3 *center, float *inv_mass_sum) {
    Vector3 c = { 0.0f, 0.0f, 0.0f };
    float sum = 0.0f;

    for (int i = 0; i < 8; ++i) {
        const Particle *p = voxel->particles[i];
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

// Resolve collisions against the scene and neighbouring voxels (mirrors ResolveCollisions compute pass).
static void solve_particle_collisions(float dt) {
    (void)dt;

    const float half_player = PLAYER_SIZE * 0.5f;
    const float omega = COLLISION_RELAXATION;
    const float eps = 1e-6f;

    // First, clamp predictions against static scene bounds and player capsules.
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        float voxel_radius = voxel_particle_radius(voxel);
        float terrain_limit = FLOOR_SIZE - voxel_radius;

        for (int j = 0; j < 8; ++j) {
            Particle *p = voxel->particles[j];

            Vector3 pos = p->predicted_pos;

            if (pos.y < voxel_radius) {
                pos.y = voxel_radius;
            }

            // pos.x = clampf(pos.x, -terrain_limit, terrain_limit);
            // pos.z = clampf(pos.z, -terrain_limit, terrain_limit);

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

        for (int j = 0; j < 8; ++j) {
            Particle *pa = voxelA->particles[j];
            float wa = pa->inv_mass;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        int nx = voxelA->gx + dx;
                        int ny = voxelA->gy + dy;
                        int nz = voxelA->gz + dz;
                        int neighbor_idx = table_get(nx, ny, nz);
                        if (neighbor_idx < 0) {
                            continue;
                        }

                        if (neighbor_idx < i) {
                            continue;
                        }

                        Voxel *voxelB = &voxels[neighbor_idx];
                        float radiusB = voxel_particle_radius(voxelB);

                        for (int q = 0; q < 8; ++q) {
                            if (neighbor_idx == i && q <= j) {
                                continue;
                            }

                            Particle *pb = voxelB->particles[q];
                            if (pa == pb) {
                                continue;
                            }

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

                            Vector3 normal = (dist > eps)
                                ? v_mul(delta, 1.0f / dist)
                                : (Vector3){ 1.0f, 0.0f, 0.0f };

                            float h = 0.5f * penetration;
                            float scale = omega * h / w_sum;

                            // Apply equal-and-opposite displacements weighted by inverse mass.
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
    }
}

static void update_particle_velocities(float dt) {
    float inv_dt = (dt > 0.0f) ? 1.0f / dt : 0.0f;

    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        Vector3 centroid = { 0.0f, 0.0f, 0.0f };
        Vector3 prev_centroid = { 0.0f, 0.0f, 0.0f };

        for (int j = 0; j < 8; ++j) {
            Particle *p = voxel->particles[j];
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

void simulate_voxel_pbd(float dt) {
    const int substeps = 3;
    const int constraint_iterations = 1;
    const float sub_dt = (substeps > 0) ? dt / (float)substeps : dt;

    for (int step = 0; step < substeps; ++step) {
        integrate_particles(sub_dt);
        solve_particle_collisions(sub_dt);

        for (int it = 0; it < constraint_iterations; ++it) {
            for (int i = 0; i < voxel_count; ++i) {
                Voxel *voxel = &voxels[i];
                if (!voxel->simulate)
                    continue;
                solve_voxel_shape(voxel);
            }
        }

        process_pending_breaks();

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
    int vix = addVoxel(start.x, start.y, start.z, false, true, col, p->vType,1.0f);
    if (vix >= 0) {
        Voxel *shot = &voxels[vix];
        Vector3 vel = v_mul(dir, 50.0f);
        shot->vel = vel;
        shot->owner = idx;
        for (int i = 0; i < 8; ++i) {
            shot->particles[i]->vel = vel;
        }
    }
}
// Append the 12 edges (24 vertices) of a cube to the current RL_LINES batch
static void drawCubeEdges(const Voxel *voxel)
{
    Vector3 v[8];
    for (int i = 0; i < 8; ++i) v[i] = voxel->particles[i]->pos;

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
    for (int i = 0; i < 8; ++i) v[i] = voxel->particles[i]->pos;

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
            const Particle *p = voxel->particles[j];
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
        // update voxel physics
        int subStep = 3;
        for( int i = 0; i < subStep; i++){
            physics_step(dt/subStep);
        }
         simulate_voxel_pbd(dt);
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
