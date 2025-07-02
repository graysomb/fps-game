// Raylib split-screen FPS prototype (port of fps_game.c)
#include "raylib.h"
#include "rlgl.h" // for rlBegin/rlEnd
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
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
#define MAX_VOXELS    20000
#define HASH_SIZE     32768    // must be power of two
#define VOXEL_SIZE     0.2f    // size of each voxel cube

// Octree constants
#define CHUNK_SIZE   16          // 16×16×16 voxels per leaf
#define MAX_LEVELS    6          // octree depth (world size = CHUNK_SIZE*2^MAX_LEVELS)


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
    bool  surface;
    int gx, gy, gz;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;
// Spatial hash table for voxels
static struct { int key, idx; } table[HASH_SIZE];


// ────────────────────────────────────────────────────────────────────────────
// Octree for greedy-meshed voxel chunks
// ────────────────────────────────────────────────────────────────────────────
typedef struct GreedyMeshChunk {
    // Spatial position (integer chunk grid, NOT voxel grid)
    int gx, gy, gz;          // chunk origin = (gx*CHUNK_SIZE)*VOXEL_SIZE …

    // GPU geometry
    Mesh  mesh;              // raw vertex data
    Model model;             // convenience wrapper returned by LoadModelFromMesh()

    // Dirty flag: 1 = voxel edits since last remesh
    bool  dirty;

    // Runtime bookkeeping  (optional: only if you want to iterate voxels quickly)
    Voxel *block;            // pointer into a contiguous CHUNK_SIZE³ slice
    int    blockCount;       // how many live voxels currently stored there
} GreedyMeshChunk;


// Each node occupies an axis-aligned cube in world space
typedef struct OctreeNode {
    BoundingBox bounds;              // min/max Vector3 in WORLD units

    // Children: NULL if this is a leaf
    struct OctreeNode *child[8];

    // Leaf payload (NULL for interior nodes)
    GreedyMeshChunk  *chunk;

    // Book-keeping
    unsigned char     level;         // 0 = root, grows toward leaves
    bool              isLeaf;        // handy boolean alias
} OctreeNode;


// Octree root and helper prototypes (Step 2 build phase)
static OctreeNode *root;
OctreeNode* BuildOctree(int chunkSize, int maxLevels);
void InsertVoxel(OctreeNode* root, int cx, int cy, int cz, Voxel* v);
// Locate leaf node for given chunk coords (Step 4 dirty propagation / Step 7 utilities)
OctreeNode* LocateLeaf(OctreeNode* n, int cx, int cy, int cz);

// Greedy-meshing hooks (TODO)
// TODO: implement BuildGreedyMesh to batch quads instead of cubes
Mesh BuildGreedyMesh(Voxel* block, int gx, int gy, int gz);
// TODO: remeshing not used in cube-batch mode
// Greedy-remeshing stubs (not used in cube-batch mode)
void GreedyRemeshLeaf(OctreeNode* leaf);
void GreedyRemeshLeaves(OctreeNode* root);

// Greedy-meshing implementation (TODO)
Mesh BuildGreedyMesh(Voxel* block, int gx, int gy, int gz) {
    // TODO: implement greedy-meshed quad batching here
    Mesh mesh = {0};
    return mesh;
}

void GreedyRemeshLeaf(OctreeNode* leaf) {
    GreedyMeshChunk* c = leaf->chunk;
    if (!c->dirty) return;
    Mesh m = BuildGreedyMesh(c->block, c->gx, c->gy, c->gz);
    if (c->model.meshCount == 0)
        c->model = LoadModelFromMesh(m);
    else
        UploadMesh(&c->model.meshes[0], false);
    c->dirty = false;
}

void GreedyRemeshLeaves(OctreeNode* n) {
    if (!n) return;
    if (n->chunk) {
        GreedyRemeshLeaf(n);
    } else {
        for (int i = 0; i < 8; i++) {
            GreedyRemeshLeaves(n->child[i]);
        }
    }
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

// Spatial hash helpers
static int hashVoxel(int x, int y, int z) {
    unsigned int h = (unsigned int)(x*73856093 ^ y*19349663 ^ z*83492791);
    return (int)(h & (HASH_SIZE - 1));
}
static void table_set(int x, int y, int z, int idx) {
    int h = hashVoxel(x,y,z);
    while (table[h].key) h = (h + 1) & (HASH_SIZE - 1);
    table[h].key = 1;
    table[h].idx = idx;
}
static int table_get(int x, int y, int z) {
    int h = hashVoxel(x,y,z);
    while (table[h].key) {
        Voxel *v = &voxels[table[h].idx];
        if (v->gx==x && v->gy==y && v->gz==z) return table[h].idx;
        h = (h + 1) & (HASH_SIZE - 1);
    }
    return -1;
}

// Check occupancy
static bool occupied(int x,int y,int z){ return  table_get(x,y,z)>=0; }

static void mark_surface(int idx){
    Voxel *v=&voxels[idx];
    int x=v->gx, y=v->gy, z=v->gz;
    v->surface = !occupied(x+1,y,z)||!occupied(x-1,y,z)||!occupied(x,y+1,z)||!occupied(x,y-1,z)||!occupied(x,y,z+1)||!occupied(x,y,z-1);
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
    v->surface = true;
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
    // int M = (int)(1.0f*FLOOR_SIZE / VOXEL_SIZE);
    // for (int x = 0; x <= M ; x++) {
    //    for (int z = 0; z <= M ; z++) {
    //        float px = (x + 0.5f) * VOXEL_SIZE-FLOOR_SIZE;
    //        float pz = (z + 0.5f) * VOXEL_SIZE-FLOOR_SIZE;
    //        addVoxel(px, 0, pz, true, false, (Color){ 150,150,150,255 }, 0);
    //    }
    // }
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

    // Build octree and meshing chunks for existing voxels
    root = BuildOctree(CHUNK_SIZE, MAX_LEVELS);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        int cx = v->gx / CHUNK_SIZE;
        int cy = v->gy / CHUNK_SIZE;
        int cz = v->gz / CHUNK_SIZE;
        InsertVoxel(root, cx, cy, cz, v);
    }
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
        // Check grid collision
        int nx = (int)floorf(v->pos.x / VOXEL_SIZE);
        int ny = (int)floorf(v->pos.y / VOXEL_SIZE);
        int nz = (int)floorf(v->pos.z / VOXEL_SIZE);
        if (nx != v->gx || ny != v->gy || nz != v->gz) {
            int hit = table_get(nx, ny, nz);
            if (hit >= 0 && hit != i) {
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
            mark_surface(hit);
            // Mark affected chunks dirty
            {
                int cxi = v->gx/CHUNK_SIZE;
                int cyi = v->gy/CHUNK_SIZE;
                int czi = v->gz/CHUNK_SIZE;
                OctreeNode *leaf = LocateLeaf(root, cxi, cyi, czi);
                if (leaf && leaf->chunk) leaf->chunk->dirty = true;
            }
            if (hit >= 0) {
                Voxel *u = &voxels[hit];
                int cxu = u->gx/CHUNK_SIZE;
                int cyu = u->gy/CHUNK_SIZE;
                int czu = u->gz/CHUNK_SIZE;
                OctreeNode *leaf2 = LocateLeaf(root, cxu, cyu, czu);
                if (leaf2 && leaf2->chunk) leaf2->chunk->dirty = true;
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
    if (vix >= 0) voxels[vix].vel = v_mul(dir, 50.0f);
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


// Draw all voxels as cubes
// Replace raw voxel draw with octree‑based frustum culling (Step 5)
// Frustum culling types & helpers (Step 5)
typedef struct Plane { Vector3 normal; float d; } Plane;
typedef struct Frustum { Plane planes[6]; } Frustum;
Frustum ExtractFrustumPlanes(Camera3D cam);
bool BoxInFrustum(const BoundingBox *box, const Frustum *fr);
Camera3D* GetCurrentCamera(void);
static void DrawOctreeNode(OctreeNode* n, Frustum* fr);

// Replace raw voxel draw with octree‑based frustum culling (Step 5)
static void DrawVoxels(void) {
    // Extract view frustum from current camera
    Camera3D* cam = GetCurrentCamera();
    Frustum fr = ExtractFrustumPlanes(*cam);
    // Batch triangles for all visible voxels
    rlBegin(RL_TRIANGLES);
        DrawOctreeNode(root, &fr);
    rlEnd();
    // Batch lines for all visible voxel edges
    rlBegin(RL_LINES);
        rlColor4ub(0, 0, 0, 255);
        DrawOctreeNode(root, &fr);
    rlEnd();
}

// Recursive frustum‑culling draw of octree
static void DrawOctreeNode(OctreeNode* n, Frustum* fr) {
    if (!n || !BoxInFrustum(&n->bounds, fr)) return;
    if (n->chunk) {
        for (int i = 0; i < n->chunk->blockCount; i++) {
            Voxel *v = &n->chunk->block[i];
            if (!v->surface) continue;
            drawCubeMan(v->pos, VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE, v->color);
        }
    } else {
        for (int i = 0; i < 8; i++) DrawOctreeNode(n->child[i], fr);
    }
}

// Stub implementations for frustum culling & camera helper
static Camera3D *currentCamera;
Frustum ExtractFrustumPlanes(Camera3D cam) { (void)cam; return (Frustum){ 0 }; }
bool BoxInFrustum(const BoundingBox *box, const Frustum *fr) { (void)box; (void)fr; return true; }
Camera3D* GetCurrentCamera(void) { return currentCamera; }

// Utility functions (Step 7)
OctreeNode* BuildOctree(int chunkSize, int maxLevels) {
    OctreeNode *n = calloc(1, sizeof(*n));
    n->level = 0;
    n->isLeaf = false;
    // World-space bounds: CHUNK_SIZE * 2^MAX_LEVELS voxels per axis
    float half = chunkSize * ((float)(1 << maxLevels)) * VOXEL_SIZE * 0.5f;
    n->bounds.min = (Vector3){ -half, -half, -half };
    n->bounds.max = (Vector3){  half,  half,  half };
    return n;
}

OctreeNode* LocateLeaf(OctreeNode* n, int cx, int cy, int cz) {
    while (n && !n->isLeaf) {
        int level = n->level;
        int shift = MAX_LEVELS - level - 1;
        int bx = (cx >> shift) & 1;
        int by = (cy >> shift) & 1;
        int bz = (cz >> shift) & 1;
        int idx = (bx << 2) | (by << 1) | bz;
        if (!n->child[idx]) {
            // Allocate child node
            OctreeNode *c = calloc(1, sizeof(*c));
            c->level = level + 1;
            c->isLeaf = (c->level == MAX_LEVELS);
            // Split bounds
            Vector3 lo = n->bounds.min;
            Vector3 hi = n->bounds.max;
            Vector3 mid = {(lo.x + hi.x)*0.5f, (lo.y + hi.y)*0.5f, (lo.z + hi.z)*0.5f};
            c->bounds.min.x = bx ? mid.x : lo.x;
            c->bounds.max.x = bx ? hi.x  : mid.x;
            c->bounds.min.y = by ? mid.y : lo.y;
            c->bounds.max.y = by ? hi.y  : mid.y;
            c->bounds.min.z = bz ? mid.z : lo.z;
            c->bounds.max.z = bz ? hi.z  : mid.z;
            if (c->isLeaf) {
                c->chunk = calloc(1, sizeof(*c->chunk));
                c->chunk->gx = cx;
                c->chunk->gy = cy;
                c->chunk->gz = cz;
                c->chunk->block = calloc(CHUNK_SIZE*CHUNK_SIZE*CHUNK_SIZE,
                                          sizeof(*c->chunk->block));
                c->chunk->blockCount = 0;
                c->chunk->dirty = true;
            }
            n->child[idx] = c;
        }
        n = n->child[idx];
    }
    return n;
}

void InsertVoxel(OctreeNode* root, int cx, int cy, int cz, Voxel* v) {
    OctreeNode *leaf = LocateLeaf(root, cx, cy, cz);
    if (!leaf || !leaf->chunk) return;
    GreedyMeshChunk *c = leaf->chunk;
    int lx = v->gx - c->gx*CHUNK_SIZE;
    int ly = v->gy - c->gy*CHUNK_SIZE;
    int lz = v->gz - c->gz*CHUNK_SIZE;
    if (lx < 0 || lx >= CHUNK_SIZE || ly < 0 || ly >= CHUNK_SIZE || lz < 0 || lz >= CHUNK_SIZE) return;
    int idx = (lz*CHUNK_SIZE + ly)*CHUNK_SIZE + lx;
    c->block[idx] = *v;
    if (c->block[idx].surface) c->blockCount++;
    c->dirty = true;
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
        currentCamera = &cam0;
        BeginMode3D(cam0);
                DrawPlane((Vector3){0,0,0}, (Vector2){FLOOR_SIZE*2, FLOOR_SIZE*2}, DARKGRAY);
                DrawVoxels();
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
        currentCamera = &cam1;
        BeginMode3D(cam1);
                DrawPlane((Vector3){0,0,0}, (Vector2){FLOOR_SIZE*2, FLOOR_SIZE*2}, DARKGRAY);
                DrawVoxels();
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
