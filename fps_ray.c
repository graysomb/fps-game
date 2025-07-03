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
#define PLAYER_RADIUS   0.5f
#define FLOOR_SIZE     20.0f    // half-size of floor in world units
#define PLAYER_SIZE 0.5f

// Voxel physics constants
#define MAX_VOXELS    131072
#define HASH_SIZE     131072    // must be power of two
#define VOXEL_SIZE     0.2f    // size of each voxel cube

// Player structure
typedef struct {
    Vector3 pos;
    float yaw;
    float pitch;
    Vector3 vel;
    bool onGround;
    int vType;
} Player;
static Player players[2];
static int scores[2] = { 0, 0 };

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
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;


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

// Add a voxel (static or dynamic)
static int addVoxel(float px, float py, float pz, bool fixed, bool simulate, Color color, int type) {
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
    table_set(v->gx, v->gy, v->gz, idx);
    return idx;
}

// Build static demo cube of voxels
static void buildDemo(void) {
    const int N = 10;
    for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) for (int z = 0; z < N; z++) {
        float px = (x + 0.5f) * VOXEL_SIZE;
        float py = (y + 0.5f) * VOXEL_SIZE;
        float pz = (z + 0.5f) * VOXEL_SIZE;
        addVoxel(px, py, pz, true, false, (Color){ 150,150,150,255 }, 0);
    }
    int M = (int)(2.0f*FLOOR_SIZE / VOXEL_SIZE);
    for (int x = 0; x <= M; x++) {
        for (int z = 0; z <= M; z++) {
            float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            addVoxel(px, 0, pz, true, false, (Color){ 150,150,150,255 }, 0);
        }
    }
}



// Reset game: players and voxels
static void ResetGame(void) {
    // init players
    players[0].pos = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
    players[1].pos = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
    players[0].yaw = 0; players[1].yaw = 180;
    players[0].pitch = players[1].pitch = 0;
    players[0].vel = players[1].vel = (Vector3){0,0,0};
    players[0].onGround = players[1].onGround = true;
    players[0].vType = players[1].vType = 0;
    // clear voxels
    voxel_count = 0;
    // clear hash
    memset(table, 0, sizeof(table));
    // build static blocks
    buildDemo();
    for (int i = 0; i < voxel_count; i++) {
        mark_surface(i);
    }
    meshDirty = true;

}

// Physics step for voxels
static void physics_step(float dt) {
    // Rebuild spatial hash
    memset(table, 0, sizeof(table));
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        int x = (int)floorf(v->pos.x / VOXEL_SIZE);
        int y = (int)floorf(v->pos.y / VOXEL_SIZE);
        int z = (int)floorf(v->pos.z / VOXEL_SIZE);
        v->gx = x; v->gy = y; v->gz = z;
        table_set(x, y, z, i);
    }
    // Simulate dynamic voxels
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;
        // Apply gravity
        v->vel.y -= GRAVITY * dt;
        // Move
        v->pos = v_add(v->pos, v_mul(v->vel, dt));
        for (int j = 0; j < 2; j++) {
            float dx = v->pos.x - players[j].pos.x;
            float dy = v->pos.y - players[j].pos.y;
            float dz = v->pos.z - players[j].pos.z;
            if (fabsf(dx) < PLAYER_SIZE && fabsf(dy) < PLAYER_SIZE && fabsf(dz) < PLAYER_SIZE) {
                if (v->owner >= 0 && v->owner != j) {
                    scores[v->owner]++;
                    printf("Player %d scored! total=%d\n", v->owner + 1, scores[v->owner]);
                    fflush(stdout);
                }
                players[j].pos     = (Vector3){ randomInRange(-9,9), BASE_EYE_HEIGHT, randomInRange(-9,9) };
                players[j].vel     = (Vector3){0,0,0};
                players[j].onGround= true;
                players[j].yaw     = (j == 0 ? 0 : 180);
                players[j].pitch   = 0;
                v->simulate = false;
                v->fixed    = true;
                v->pos      = (Vector3){ -999.0f, -999.0f, -999.0f };
                break;
            }
        }
        if (!v->simulate) continue;
        // Check grid collision
        int nx = (int)floorf(v->pos.x / VOXEL_SIZE);
        int ny = (int)floorf(v->pos.y / VOXEL_SIZE);
        int nz = (int)floorf(v->pos.z / VOXEL_SIZE);
        if (nx != v->gx || ny != v->gy || nz != v->gz) {
            int hit = table_get(nx, ny, nz);
            if (hit >= 0 && hit != i) {

                //reset mesh
                meshDirty = true;

                if (v->type == 1) {
                    // Destructive collision: remove bullet and hit block
                    v->simulate = false;
                    v->fixed = true;
                    v->pos = (Vector3){-999.0f, -999.0f, -999.0f};
                    // Remove existing block
                    Voxel *u = &voxels[hit];
                    u->simulate = false;
                    u->fixed = true;
                    u->pos = (Vector3){-999.0f, -999.0f, -999.0f};
                } else {
                    v->simulate = false;
                    v->fixed = true;
                    v->pos = (Vector3){(v->gx + 0.5f) * VOXEL_SIZE,
                                (v->gy + 0.5f) * VOXEL_SIZE,
                                (v->gz + 0.5f) * VOXEL_SIZE};
                }
            }
            mark_surface(i);
            if (hit >= 0) mark_surface(hit);
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
        voxels[vix].vel = v_mul(dir, 50.0f);
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

static int first_voxel_hit(Ray ray, float t_max) {
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
        if (id >= 0) {                // hit!
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



// Generate a greedy mesh of all visible voxels
static Mesh gen_greedy_mesh(void) {
    Mesh mesh = { 0 };
    patchCount = 0;

    // Collect all static (non-simulated) voxels by visible faces (principal planes)
    int xyCount = 0, xzCount = 0, yzCount = 0;
    int xyList[MAX_VOXELS], xzList[MAX_VOXELS], yzList[MAX_VOXELS];
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (v->simulate) continue;
        if (v->surface[4] || v->surface[5]) xyList[xyCount++] = i;  // XY-plane (+Z/-Z)
        if (v->surface[2] || v->surface[3]) xzList[xzCount++] = i;  // XZ-plane (+Y/-Y)
        if (v->surface[0] || v->surface[1]) yzList[yzCount++] = i;  // YZ-plane (+X/-X)
    }

    // Merge rectangles in each principal plane via a greedy 2D algorithm
    // XY-plane merging (group by z)
    if (xyCount) {
        int zLayers[xyCount], zLayerCount = 0;
        for (int ii = 0; ii < xyCount; ii++) {
            int idx = xyList[ii]; int z = voxels[idx].gz;
            bool seen = false;
            for (int ll = 0; ll < zLayerCount; ll++) if (zLayers[ll] == z) { seen = true; break; }
            if (!seen) zLayers[zLayerCount++] = z;
        }
        for (int zl = 0; zl < zLayerCount; zl++) {
            int z = zLayers[zl];
            // find bounds
            int minX = INT_MAX, maxX = INT_MIN, minY = INT_MAX, maxY = INT_MIN;
            for (int ii = 0; ii < xyCount; ii++) {
                Voxel *v = &voxels[xyList[ii]];
                if (v->gz != z) continue;
                if (v->gx < minX) minX = v->gx;
                if (v->gx > maxX) maxX = v->gx;
                if (v->gy < minY) minY = v->gy;
                if (v->gy > maxY) maxY = v->gy;
            }
            int w = maxX - minX + 1;
            int h = maxY - minY + 1;
            bool mask[w][h]; memset(mask, 0, sizeof mask);
            for (int ii = 0; ii < xyCount; ii++) {
                Voxel *v = &voxels[xyList[ii]];
                if (v->gz != z) continue;
                mask[v->gx - minX][v->gy - minY] = true;
            }
            for (int yy = 0; yy < h; yy++) {
                for (int xx = 0; xx < w; xx++) {
                    if (!mask[xx][yy]) continue;
                    // start cell
                    int sx = xx, sy = yy;
                    // determine width
                    int ww = 1;
                    while (sx + ww < w && mask[sx + ww][sy]) ww++;
                    // determine height
                    int hh = 1;
                    for (;;) {
                        bool block = false;
                        for (int k = 0; k < ww; k++) {
                            if (!mask[sx + k][sy + hh]) { block = true; break; }
                        }
                        if (block || sy + hh >= h) break;
                        hh++;
                    }
                    // mark used
                    for (int dy = 0; dy < hh; dy++)
                        for (int dx = 0; dx < ww; dx++)
                            mask[sx + dx][sy + dy] = false;
                    xx = sx + ww - 1;
                    // record merged patch in XY-plane at z
                    {
                        Patch *pt = &patches[patchCount++];
                        pt->plane    = 0;
                        pt->layer    = z;
                        pt->i0       = minX + sx;
                        pt->j0       = minY + sy;
                        pt->di       = ww;
                        pt->dj       = hh;
                        int baseIdx = table_get(minX + sx, minY + sy, z);
                        pt->positive = voxels[baseIdx].surface[4];
                        pt->col      = voxels[baseIdx].color;
                    }
                }
            }
        }
    }

    // XZ-plane merging (group by y)
    if (xzCount) {
        int yLayers[xzCount], yLayerCount = 0;
        for (int ii = 0; ii < xzCount; ii++) {
            int idx = xzList[ii]; int y = voxels[idx].gy;
            bool seen = false;
            for (int ll = 0; ll < yLayerCount; ll++) if (yLayers[ll] == y) { seen = true; break; }
            if (!seen) yLayers[yLayerCount++] = y;
        }
        for (int yl = 0; yl < yLayerCount; yl++) {
            int y = yLayers[yl];
            int minX = INT_MAX, maxX = INT_MIN, minZ = INT_MAX, maxZ = INT_MIN;
            for (int ii = 0; ii < xzCount; ii++) {
                Voxel *v = &voxels[xzList[ii]];
                if (v->gy != y) continue;
                if (v->gx < minX) minX = v->gx;
                if (v->gx > maxX) maxX = v->gx;
                if (v->gz < minZ) minZ = v->gz;
                if (v->gz > maxZ) maxZ = v->gz;
            }
            int w = maxX - minX + 1;
            int h = maxZ - minZ + 1;
            bool mask[w][h]; memset(mask, 0, sizeof mask);
            for (int ii = 0; ii < xzCount; ii++) {
                Voxel *v = &voxels[xzList[ii]];
                if (v->gy != y) continue;
                mask[v->gx - minX][v->gz - minZ] = true;
            }
            for (int zz = 0; zz < h; zz++) {
                for (int xx = 0; xx < w; xx++) {
                    if (!mask[xx][zz]) continue;
                    // start cell
                    int sx = xx, sz = zz;
                    // width
                    int ww = 1;
                    while (sx + ww < w && mask[sx + ww][sz]) ww++;
                    // height
                    int hh = 1;
                    for (;;) {
                        bool block = false;
                        for (int k = 0; k < ww; k++) {
                            if (!mask[sx + k][sz + hh]) { block = true; break; }
                        }
                        if (block || sz + hh >= h) break;
                        hh++;
                    }
                    // mark used
                    for (int dz = 0; dz < hh; dz++)
                        for (int dx = 0; dx < ww; dx++)
                            mask[sx + dx][sz + dz] = false;
                    xx = sx + ww - 1;
                    // record merged patch in XZ-plane at y
                    {
                        Patch *pt = &patches[patchCount++];
                        pt->plane    = 1;
                        pt->layer    = y;
                        pt->i0       = minX + sx;
                        pt->j0       = minZ + sz;
                        pt->di       = ww;
                        pt->dj       = hh;
                        int baseIdx = table_get(minX + sx, y, minZ + sz);
                        pt->positive = voxels[baseIdx].surface[2];
                        pt->col      = voxels[baseIdx].color;
                    }
                }
            }
        }
    }

    // YZ-plane merging (group by x)
    if (yzCount) {
        int xLayers[yzCount], xLayerCount = 0;
        for (int ii = 0; ii < yzCount; ii++) {
            int idx = yzList[ii]; int x = voxels[idx].gx;
            bool seen = false;
            for (int ll = 0; ll < xLayerCount; ll++) if (xLayers[ll] == x) { seen = true; break; }
            if (!seen) xLayers[xLayerCount++] = x;
        }
        for (int xl = 0; xl < xLayerCount; xl++) {
            int x = xLayers[xl];
            int minY = INT_MAX, maxY = INT_MIN, minZ = INT_MAX, maxZ = INT_MIN;
            for (int ii = 0; ii < yzCount; ii++) {
                Voxel *v = &voxels[yzList[ii]];
                if (v->gx != x) continue;
                if (v->gy < minY) minY = v->gy;
                if (v->gy > maxY) maxY = v->gy;
                if (v->gz < minZ) minZ = v->gz;
                if (v->gz > maxZ) maxZ = v->gz;
            }
            int w = maxY - minY + 1;
            int h = maxZ - minZ + 1;
            bool mask[w][h]; memset(mask, 0, sizeof mask);
            for (int ii = 0; ii < yzCount; ii++) {
                Voxel *v = &voxels[yzList[ii]];
                if (v->gx != x) continue;
                mask[v->gy - minY][v->gz - minZ] = true;
            }
            for (int zz = 0; zz < h; zz++) {
                for (int yy = 0; yy < w; yy++) {
                    if (!mask[yy][zz]) continue;
                    // start cell
                    int sy = yy, sz = zz;
                    // width
                    int ww = 1;
                    while (sy + ww < w && mask[sy + ww][sz]) ww++;
                    // height
                    int hh = 1;
                    for (;;) {
                        bool block = false;
                        for (int k = 0; k < ww; k++) {
                            if (!mask[sy + k][sz + hh]) { block = true; break; }
                        }
                        if (block || sz + hh >= h) break;
                        hh++;
                    }
                    // mark used
                    for (int dz = 0; dz < hh; dz++)
                        for (int dy = 0; dy < ww; dy++)
                            mask[sy + dy][sz + dz] = false;
                    yy = sy + ww - 1;
                    // record merged patch in YZ-plane at x
                    {
                        Patch *pt = &patches[patchCount++];
                        pt->plane    = 2;
                        pt->layer    = x;
                        pt->i0       = minY + sy;
                        pt->j0       = minZ + sz;
                        pt->di       = ww;
                        pt->dj       = hh;
                        int baseIdx = table_get(x, minY + sy, minZ + sz);
                        pt->positive = voxels[baseIdx].surface[0];
                        pt->col      = voxels[baseIdx].color;
                    }
                }
            }
        }
    }

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
        drawCubeMan(v->pos, VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE, v->color);
        }
        
    }
    rlEnd();
}

static void draw_players(void) {
    for (int i = 0; i < 2; i++) {
        Player *p = &players[i];
        DrawCube(p->pos, PLAYER_SIZE,PLAYER_SIZE, PLAYER_SIZE, BLACK);
        DrawCubeWires(p->pos, PLAYER_SIZE,PLAYER_SIZE,PLAYER_SIZE, BLACK);
    }
}

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
        float dt = GetFrameTime();
        // input: shooting and jump
        if (IsKeyPressed(KEY_LEFT_CONTROL))  FireVoxel(0);
        if (IsKeyPressed(KEY_RIGHT_CONTROL)) FireVoxel(1);
        if (IsKeyPressed(KEY_Q)) players[0].vType = 1-players[0].vType;
        if (IsKeyPressed(KEY_U)) players[0].vType = 1-players[0].vType;
        if (IsKeyPressed(KEY_SPACE) && players[0].onGround) { players[0].vel.y = JUMP_SPEED; players[0].onGround = false; }
        if (IsKeyPressed(KEY_RIGHT_SHIFT) && players[1].onGround) { players[1].vel.y = JUMP_SPEED; players[1].onGround = false; }
        // update players
        for (int i = 0; i < 2; i++) {
            Player *p = &players[i];
            // turn
            if ((i==0 && IsKeyDown(KEY_H)) || (i==1 && IsKeyDown(KEY_RIGHT))) p->yaw -= TURN_SPEED*dt;
            if ((i==0 && IsKeyDown(KEY_F)) || (i==1 && IsKeyDown(KEY_LEFT)))  p->yaw += TURN_SPEED*dt;
            // look up/down
            if ((i==0 && IsKeyDown(KEY_T)) || (i==1 && IsKeyDown(KEY_UP)))   p->pitch += TURN_SPEED*dt;
            if ((i==0 && IsKeyDown(KEY_G)) || (i==1 && IsKeyDown(KEY_DOWN))) p->pitch -= TURN_SPEED*dt;
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
                p->vel = v_add(p->vel, v_mul(accel, ACCELERATION*dt));
            } else {
                // friction
                float sp = sqrtf(p->vel.x*p->vel.x + p->vel.z*p->vel.z);
                if (sp > 0) {
                    float dec = FRICTION*dt;
                    float ns = sp - dec; if (ns < 0) ns = 0;
                    p->vel.x *= ns/sp;
                    p->vel.z *= ns/sp;
                }
            }
            // vertical
            p->vel.y -= GRAVITY*dt;
            p->pos.y += p->vel.y*dt;
            if (p->pos.y <= BASE_EYE_HEIGHT) {
                p->pos.y = BASE_EYE_HEIGHT;
                p->vel.y = 0;
                p->onGround = true;
            }
            //clamp speed
            {
                float speed = sqrtf(p->vel.x*p->vel.x + p->vel.z*p->vel.z);
                if (speed > MOVE_SPEED) {
                    p->vel.x *= MOVE_SPEED/speed;
                    p->vel.z *= MOVE_SPEED/speed;
                }
            }
            // horizontal move
            p->pos.x += p->vel.x*dt;
            p->pos.z += p->vel.z*dt;
            // bounds
            p->pos.x = clampf(p->pos.x, -FLOOR_SIZE+PLAYER_RADIUS, FLOOR_SIZE-PLAYER_RADIUS);
            p->pos.z = clampf(p->pos.z, -FLOOR_SIZE+PLAYER_RADIUS, FLOOR_SIZE-PLAYER_RADIUS);
        }
        // update voxel physics
        // int subStep = 8;
        // for( int i = 0; i < subStep; i++){
        //     physics_step(dt/subStep);
        // }
        physics_step(dt);
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
            DrawRectangle(0,0, SCREEN_WIDTH/2, 40, Fade(BLACK, 0.5f));
            DrawText("P1: WASD move, H/F turn, T/G look, LCTRL shoot", 10,10,20,WHITE);
            DrawText(TextFormat("Score: %d", scores[0]), SCREEN_WIDTH/2 -150,10,20,WHITE);
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
            DrawRectangle(0,0, SCREEN_WIDTH/2, 40, Fade(BLACK, 0.5f));
            DrawText("P2: I/K move, RIGHT/LEFT turn, UP/DOWN look, RCTRL shoot", 10,10,20,WHITE);
            DrawText(TextFormat("Score: %d", scores[1]), SCREEN_WIDTH/2 -150,10,20,WHITE);
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
    }
    // cleanup
    UnloadRenderTexture(screen0);
    UnloadRenderTexture(screen1);
    CloseWindow();
    return 0;
}
