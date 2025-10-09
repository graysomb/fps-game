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
#define BRIDGE_CLUSTER_THRESHOLD 3

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
    int groupId;
    /* Visible faces mask: surface[i]=true if face i is visible:
       0=+X,1=-X,2=+Y,3=-Y,4=+Z,5=-Z */
    bool surface[6];
    int  gx, gy, gz;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;
static int voxelGroupCount = 0;
static int voxelGroupQueue[MAX_VOXELS];
static int voxelGroupRenderCount = 0;

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

typedef struct {
    Mesh mesh;
    Patch *patches;
    int patchCount;
    bool dirty;
    Vector3 drawOffset;
} VoxelGroupRender;

#define MAX_VOXEL_GROUPS MAX_VOXELS
static VoxelGroupRender voxelGroupRender[MAX_VOXEL_GROUPS];
static bool groupsDirty = true;

typedef struct {
    bool active;
    bool falling;
    Vector3 velocity;
} VoxelGroupState;

static VoxelGroupState voxelGroupState[MAX_VOXEL_GROUPS];
static int groupMemberOffsets[MAX_VOXEL_GROUPS + 1];
static int groupMemberList[MAX_VOXELS];
static int groupMemberCounts[MAX_VOXEL_GROUPS];
static int groupMemberCursor[MAX_VOXEL_GROUPS];
static const float GROUP_GRAVITY = GRAVITY;
static int voxelLocalIndex[MAX_VOXELS];

static void unload_group_mesh(VoxelGroupRender *grp) {
    if (grp->mesh.vertices) {
        UnloadMesh(grp->mesh);
        grp->mesh = (Mesh){ 0 };
    }
    if (grp->patches) {
        free(grp->patches);
        grp->patches = NULL;
    }
    grp->patchCount = 0;
    grp->drawOffset = (Vector3){0,0,0};
}

static void mark_group_dirty(int groupId) {
    if (groupId >= 0 && groupId < MAX_VOXEL_GROUPS) {
        voxelGroupRender[groupId].dirty = true;
    }
}

static void mark_all_groups_dirty(void) {
    int limit = voxelGroupRenderCount;
    if (limit > MAX_VOXEL_GROUPS) limit = MAX_VOXEL_GROUPS;
    for (int i = 0; i < limit; i++) {
        mark_group_dirty(i);
    }
    groupsDirty = true;
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

static bool voxel_is_removed(const Voxel *v) {
    return (v->pos.x <= -900.0f && v->pos.y <= -900.0f && v->pos.z <= -900.0f);
}

static bool voxel_is_static_active(const Voxel *v) {
    return (!v->simulate && !voxel_is_removed(v));
}
static void refine_groups_break_bridges(void);

static Color color_for_group(int groupId) {
    if (groupId < 0) {
        return (Color){ 200, 200, 200, 255 };
    }
    float hue = fmodf((float)(groupId * 61), 360.0f);
    float saturation = 0.55f + 0.35f * fmodf((float)(groupId * 17) / 23.0f, 1.0f);
    float value = 0.85f;
    return ColorFromHSV(hue, saturation, value);
}

static Color voxel_group_color(const Voxel *v) {
    if (voxel_is_static_active(v) && v->groupId >= 0) {
        return color_for_group(v->groupId);
    }
    return v->color;
}

static void compute_static_voxel_groups(void) {
    // Flood-fill static voxels into 6-connected groups for visualization.
    static const int neighborOffsets[6][3] = {
        { 1,  0,  0 }, { -1,  0,  0 },
        { 0,  1,  0 }, {  0, -1,  0 },
        { 0,  0,  1 }, {  0,  0, -1 }
    };

    for (int i = 0; i < voxel_count; i++) {
        voxels[i].groupId = -1;
    }

    voxelGroupCount = 0;
    for (int i = 0; i < voxel_count; i++) {
        Voxel *root = &voxels[i];
        if (!voxel_is_static_active(root)) continue;
        if (root->groupId >= 0) continue;

        int groupId = voxelGroupCount++;
        int head = 0;
        int tail = 0;
        root->groupId = groupId;
        voxelGroupQueue[tail++] = i;

        while (head < tail) {
            int idx = voxelGroupQueue[head++];
            Voxel *current = &voxels[idx];
            for (int n = 0; n < 6; n++) {
                int nx = current->gx + neighborOffsets[n][0];
                int ny = current->gy + neighborOffsets[n][1];
                int nz = current->gz + neighborOffsets[n][2];
                int neighborIdx = table_get(nx, ny, nz);
                if (neighborIdx < 0) continue;
                Voxel *neighbor = &voxels[neighborIdx];
                if (!voxel_is_static_active(neighbor)) continue;
                if (neighbor->groupId >= 0) continue;
                neighbor->groupId = groupId;
                if (tail < MAX_VOXELS) {
                    voxelGroupQueue[tail++] = neighborIdx;
                }
            }
        }
    }
    refine_groups_break_bridges();
}

static void gather_group_members(void) {
    if (voxelGroupCount <= 0) {
        memset(groupMemberOffsets, 0, sizeof(groupMemberOffsets));
        return;
    }
    int limit = (voxelGroupCount > MAX_VOXEL_GROUPS) ? MAX_VOXEL_GROUPS : voxelGroupCount;
    memset(groupMemberCounts, 0, sizeof(int) * limit);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_static_active(v)) continue;
        if (v->groupId < 0 || v->groupId >= limit) continue;
        groupMemberCounts[v->groupId]++;
    }
    int total = 0;
    for (int g = 0; g < limit; g++) {
        groupMemberOffsets[g] = total;
        total += groupMemberCounts[g];
    }
    groupMemberOffsets[limit] = total;
    memcpy(groupMemberCursor, groupMemberOffsets, sizeof(int) * limit);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_static_active(v)) continue;
        if (v->groupId < 0 || v->groupId >= limit) continue;
        int cursor = groupMemberCursor[v->groupId];
        groupMemberList[cursor] = i;
        groupMemberCursor[v->groupId] = cursor + 1;
    }
    if (limit < MAX_VOXEL_GROUPS) {
        groupMemberOffsets[limit + 1] = total;
    }
}

static void refine_groups_break_bridges(void) {
    if (voxelGroupCount <= 0) return;
    int groupLimit = voxelGroupCount;
    int *counts = calloc(groupLimit, sizeof(int));
    if (!counts) return;
    int staticCount = 0;
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_static_active(v)) continue;
        if (v->groupId < 0 || v->groupId >= groupLimit) continue;
        counts[v->groupId]++;
        staticCount++;
    }
    int *offsets = calloc(groupLimit + 1, sizeof(int));
    if (!offsets) {
        free(counts);
        return;
    }
    int running = 0;
    for (int g = 0; g < groupLimit; g++) {
        offsets[g] = running;
        running += counts[g];
    }
    offsets[groupLimit] = running;
    int *cursor = malloc(sizeof(int) * groupLimit);
    int *members = malloc(sizeof(int) * running);
    if (!members) {
        free(cursor);
        free(offsets);
        free(counts);
        return;
    }
    if (!cursor) {
        free(members);
        free(offsets);
        free(counts);
        return;
    }
    memcpy(cursor, offsets, sizeof(int) * groupLimit);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_static_active(v)) continue;
        int gid = v->groupId;
        if (gid < 0 || gid >= groupLimit) continue;
        int pos = cursor[gid]++;
        members[pos] = i;
    }
    int newGroupCounter = 0;
    for (int g = 0; g < groupLimit; g++) {
        int start = offsets[g];
        int end = offsets[g + 1];
        int localCount = end - start;
        if (localCount <= 0) continue;
        for (int idx = start; idx < end; idx++) {
            voxelLocalIndex[members[idx]] = idx - start;
        }
        int localSize = localCount;
        int *degrees = malloc(sizeof(int) * localSize);
        unsigned char *candidate = calloc(localSize, sizeof(unsigned char));
        unsigned char *clusterVisited = calloc(localSize, sizeof(unsigned char));
        unsigned char *bridgeFlag = calloc(localSize, sizeof(unsigned char));
        int *queue = malloc(sizeof(int) * localSize);
        if (!degrees || !candidate || !clusterVisited || !bridgeFlag || !queue) {
            free(queue);
            free(bridgeFlag);
            free(clusterVisited);
            free(candidate);
            free(degrees);
            for (int idx = start; idx < end; idx++) {
                voxelLocalIndex[members[idx]] = -1;
            }
            continue;
        }
        const int dirs[6][3] = {
            {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}
        };
        for (int li = 0; li < localSize; li++) {
            degrees[li] = 0;
            int globalIdx = members[start + li];
            Voxel *v = &voxels[globalIdx];
            for (int n = 0; n < 6; n++) {
                int nx = v->gx + dirs[n][0];
                int ny = v->gy + dirs[n][1];
                int nz = v->gz + dirs[n][2];
                int neighborIdx = table_get(nx, ny, nz);
                if (neighborIdx < 0) continue;
                int localNeighbor = voxelLocalIndex[neighborIdx];
                if (localNeighbor >= 0 && localNeighbor < localSize) {
                    degrees[li]++;
                }
            }
            if (degrees[li] <= 2) {
                candidate[li] = 1;
            }
        }
        for (int li = 0; li < localSize; li++) {
            if (!candidate[li] || clusterVisited[li]) continue;
            int qh = 0, qt = 0;
            queue[qt++] = li;
            clusterVisited[li] = 1;
            int clusterSize = 0;
            int *clusterList = malloc(sizeof(int) * localSize);
            int adjCount = 0;
            int *adjList = malloc(sizeof(int) * localSize);
            if (!clusterList || !adjList) {
                free(adjList);
                free(clusterList);
                break;
            }
            unsigned char *adjMark = calloc(localSize, sizeof(unsigned char));
            if (!adjMark) {
                free(adjList);
                free(clusterList);
                break;
            }
            while (qh < qt) {
                int cur = queue[qh++];
                clusterList[clusterSize++] = cur;
                int globalIdx = members[start + cur];
                Voxel *v = &voxels[globalIdx];
                for (int n = 0; n < 6; n++) {
                    int nx = v->gx + dirs[n][0];
                    int ny = v->gy + dirs[n][1];
                    int nz = v->gz + dirs[n][2];
                    int neighborIdx = table_get(nx, ny, nz);
                    if (neighborIdx < 0) continue;
                    int localNeighbor = voxelLocalIndex[neighborIdx];
                    if (localNeighbor < 0 || localNeighbor >= localSize) continue;
                    if (candidate[localNeighbor]) {
                        if (!clusterVisited[localNeighbor]) {
                            clusterVisited[localNeighbor] = 1;
                            queue[qt++] = localNeighbor;
                        }
                    } else {
                        if (!adjMark[localNeighbor]) {
                            adjMark[localNeighbor] = 1;
                            adjList[adjCount++] = localNeighbor;
                        }
                    }
                }
            }
            if (clusterSize < BRIDGE_CLUSTER_THRESHOLD && adjCount >= 2) {
                for (int ci = 0; ci < clusterSize; ci++) {
                    bridgeFlag[clusterList[ci]] = 1;
                }
            }
            free(adjMark);
            free(adjList);
            free(clusterList);
        }
        int *componentId = malloc(sizeof(int) * localSize);
        if (!componentId) {
            free(queue);
            free(bridgeFlag);
            free(clusterVisited);
            free(candidate);
            free(degrees);
            for (int idx = start; idx < end; idx++) {
                voxelLocalIndex[members[idx]] = -1;
            }
            continue;
        }
        for (int li = 0; li < localSize; li++) componentId[li] = -1;
        int componentCount = 0;
        for (int li = 0; li < localSize; li++) {
            if (componentId[li] >= 0) continue;
            int qh = 0, qt = 0;
            queue[qt++] = li;
            componentId[li] = componentCount;
            unsigned char bridgeOnly = bridgeFlag[li];
            while (qh < qt) {
                int cur = queue[qh++];
                int globalIdx = members[start + cur];
                Voxel *v = &voxels[globalIdx];
                for (int n = 0; n < 6; n++) {
                    int nx = v->gx + dirs[n][0];
                    int ny = v->gy + dirs[n][1];
                    int nz = v->gz + dirs[n][2];
                    int neighborIdx = table_get(nx, ny, nz);
                    if (neighborIdx < 0) continue;
                    int localNeighbor = voxelLocalIndex[neighborIdx];
                    if (localNeighbor < 0 || localNeighbor >= localSize) continue;
                    if (componentId[localNeighbor] >= 0) continue;
                    if (bridgeFlag[cur] || bridgeFlag[localNeighbor]) {
                        if (!(bridgeFlag[cur] && bridgeFlag[localNeighbor])) continue;
                    }
                    componentId[localNeighbor] = componentCount;
                    queue[qt++] = localNeighbor;
                }
            }
            componentCount++;
        }
        int *componentMap = malloc(sizeof(int) * componentCount);
        if (!componentMap) {
            free(componentId);
            free(queue);
            free(bridgeFlag);
            free(clusterVisited);
            free(candidate);
            free(degrees);
            for (int idx = start; idx < end; idx++) {
                voxelLocalIndex[members[idx]] = -1;
            }
            continue;
        }
        for (int ci = 0; ci < componentCount; ci++) {
            componentMap[ci] = newGroupCounter++;
        }
        for (int li = 0; li < localSize; li++) {
            int globalIdx = members[start + li];
            int comp = componentId[li];
            if (comp >= 0) {
                voxels[globalIdx].groupId = componentMap[comp];
            }
        }
        for (int idx = start; idx < end; idx++) {
            voxelLocalIndex[members[idx]] = -1;
        }
        free(componentMap);
        free(componentId);
        free(queue);
        free(bridgeFlag);
        free(clusterVisited);
        free(candidate);
        free(degrees);
    }
    voxelGroupCount = newGroupCounter;
    free(cursor);
    free(members);
    free(offsets);
    free(counts);
}

static void translate_group(int groupId, float dy) {
    if (groupId < 0 || groupId >= voxelGroupCount) return;
    if (fabsf(dy) < 1e-6f) return;
    int start = groupMemberOffsets[groupId];
    int end = groupMemberOffsets[groupId + 1];
    if (start >= end) return;
    for (int idx = start; idx < end; idx++) {
        Voxel *v = &voxels[groupMemberList[idx]];
        table_remove(v->gx, v->gy, v->gz);
    }
    for (int idx = start; idx < end; idx++) {
        Voxel *v = &voxels[groupMemberList[idx]];
        v->pos.y += dy;
        v->gy = (int)floorf(v->pos.y / VOXEL_SIZE);
    }
    for (int idx = start; idx < end; idx++) {
        int vi = groupMemberList[idx];
        Voxel *v = &voxels[vi];
        table_set(v->gx, v->gy, v->gz, vi);
    }
    for (int idx = start; idx < end; idx++) {
        int vi = groupMemberList[idx];
        mark_surface(vi);
        mark_surface_neighbors(voxels[vi].pos);
    }
    if (groupId >= 0 && groupId < MAX_VOXEL_GROUPS) {
        voxelGroupRender[groupId].drawOffset.y += dy;
    }
}

static void ensure_voxel_groups_up_to_date(void) {
    if (!groupsDirty) return;
    int prevRenderCount = voxelGroupRenderCount;
    compute_static_voxel_groups();
    if (voxelGroupCount > MAX_VOXEL_GROUPS) {
        voxelGroupCount = MAX_VOXEL_GROUPS;
    }
    for (int i = 0; i < voxelGroupCount; i++) {
        voxelGroupRender[i].dirty = true;
        voxelGroupState[i].active = true;
        if (!voxelGroupState[i].falling) {
            voxelGroupState[i].velocity = (Vector3){0,0,0};
        }
    }
    for (int i = voxelGroupCount; i < prevRenderCount; i++) {
        unload_group_mesh(&voxelGroupRender[i]);
        voxelGroupRender[i].dirty = true;
        voxelGroupState[i].active = false;
        voxelGroupState[i].falling = false;
        voxelGroupState[i].velocity = (Vector3){0,0,0};
    }
    voxelGroupRenderCount = voxelGroupCount;
    groupsDirty = false;
}

static bool group_supported(int groupId) {
    if (groupId < 0 || groupId >= voxelGroupCount) return false;
    int start = groupMemberOffsets[groupId];
    int end = groupMemberOffsets[groupId + 1];
    if (start >= end) return true;
    const float half = VOXEL_SIZE * 0.5f;
    for (int idx = start; idx < end; idx++) {
        Voxel *v = &voxels[groupMemberList[idx]];
        float bottom = v->pos.y - half;
        if (bottom <= 0.001f) {
            return true;
        }
        int belowIdx = table_get(v->gx, v->gy - 1, v->gz);
        if (belowIdx >= 0) {
            Voxel *b = &voxels[belowIdx];
            if (voxel_is_static_active(b) && b->groupId != groupId) {
                return true;
            }
        }
    }
    return false;
}

static void update_group_physics(float dt) {
    ensure_voxel_groups_up_to_date();
    if (voxelGroupCount <= 0) return;
    gather_group_members();
    int limit = (voxelGroupCount > MAX_VOXEL_GROUPS) ? MAX_VOXEL_GROUPS : voxelGroupCount;
    const float half = VOXEL_SIZE * 0.5f;
    bool regroupNeeded = false;
    for (int g = 0; g < limit; g++) {
        int start = groupMemberOffsets[g];
        int end   = groupMemberOffsets[g + 1];
        VoxelGroupState *state = &voxelGroupState[g];
        if (start >= end) {
            state->active = false;
            state->falling = false;
            state->velocity = (Vector3){0,0,0};
            continue;
        }
        state->active = true;
        bool supported = group_supported(g);
        if (supported) {
            if (state->falling) {
                state->falling = false;
                state->velocity = (Vector3){0,0,0};
            }
            continue;
        }
        if (!state->falling) {
            state->falling = true;
            state->velocity = (Vector3){0,0,0};
        }
        state->velocity.y -= GROUP_GRAVITY * dt;
        float desired = state->velocity.y * dt;
        float allowed = desired;
        int collidedGroup = -1;
        for (int idx = start; idx < end; idx++) {
            Voxel *v = &voxels[groupMemberList[idx]];
            float bottom = v->pos.y - half;
            float groundLimit = -bottom;
            if (allowed < groundLimit) {
                allowed = groundLimit;
                collidedGroup = -1;
            }
            int belowIdx = table_get(v->gx, v->gy - 1, v->gz);
            if (belowIdx >= 0) {
                Voxel *b = &voxels[belowIdx];
                if (voxel_is_static_active(b) && b->groupId != g) {
                    float occupantTop = (b->gy + 1) * VOXEL_SIZE;
                    float dropLimit = occupantTop - bottom;
                    if (allowed < dropLimit) {
                        allowed = dropLimit;
                        collidedGroup = b->groupId;
                    }
                }
            }
        }
        if (fabsf(allowed) > 1e-6f) {
            translate_group(g, allowed);
        }
        if (allowed > desired + 1e-5f) {
            state->falling = false;
            state->velocity = (Vector3){0,0,0};
            if (collidedGroup >= 0) {
                mark_group_dirty(collidedGroup);
            }
            regroupNeeded = true;
        }
    }
    if (regroupNeeded) {
        groupsDirty = true;
    }
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
    v->groupId = -1;
    memset(v->surface, 0, sizeof v->surface);
    // compute grid coords
    v->gx = (int)floorf(px / VOXEL_SIZE);
    v->gy = (int)floorf(py / VOXEL_SIZE);
    v->gz = (int)floorf(pz / VOXEL_SIZE);
    table_set(v->gx, v->gy, v->gz, idx);
    if (!simulate) {
        groupsDirty = true;
    }
    return idx;
}

// Build static demo cube of voxels
static void buildDemo(void) {
    // Floor
    int M = (int)(2.0f * FLOOR_SIZE / VOXEL_SIZE);
    for (int x = 0; x <= M; x++) {
        for (int z = 0; z <= M; z++) {
            float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            addVoxel(px, 0, pz, true, false, (Color){ 150, 150, 150, 255 }, 0);
        }
    }

    // Pillars
    int pillar_height = 35; // 45 - 10
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
        for (int y = 1; y <= pillar_height; y++) {
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

    // Central platform
    int platform_size = M / 5;
    int platform_height = 5; // 15 / 3
    int platform_base_height = 16; // to keep top at same level (21)
    for (int y = platform_base_height; y <= platform_base_height + platform_height; y++) {
        for (int x = M/2 - platform_size/2; x <= M/2 + platform_size/2; x++) {
            for (int z = M/2 - platform_size/2; z <= M/2 + platform_size/2; z++) {
                float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                float py = (y + 0.5f) * VOXEL_SIZE;
                float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
                addVoxel(px, py, pz, true, false, (Color){ 100, 200, 100, 255 }, 0);
            }
        }
    }
}



static int first_voxel_hit(Ray ray, float t_max, int ignore_id);
static void UpdateKdRatio(int player_index);
static void gather_group_members(void);
static void translate_group(int groupId, float dy);
static void update_group_physics(float dt);

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
    for (int i = 0; i < voxel_count; i++) {
        mark_surface(i);
    }
    voxelGroupCount = 0;
    voxelGroupRenderCount = 0;
    mark_all_groups_dirty();

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
    // Simulate dynamic voxels
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!v->simulate) continue;

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
                v->simulate = false;
                v->fixed = true;
                v->pos = (Vector3){-999.0f, -999.0f, -999.0f};
                v->groupId = -1;

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
                                    victim->simulate = false;
                                    victim->fixed = true;
                                    victim->pos = (Vector3){-999.0f, -999.0f, -999.0f};
                                    victim->groupId = -1;
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
                                    }
                                }
                            }
                        }
                    }
                }
                mark_all_groups_dirty();
            }
        }

        if (!hit_voxel) {
            // Move
            v->pos = v_add(v->pos, displacement);
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

                    v->simulate = false;
                    v->fixed    = true;
                    v->pos      = (Vector3){ -999.0f, -999.0f, -999.0f };
                    v->groupId  = -1;

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

static void merge_rects_on_plane(int count, int *list, int plane, bool positive,
                                 Patch *patchBuffer, int patchCapacity, int *patchCount) {
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
                
                if (*patchCount >= patchCapacity) {
                    continue;
                }
                Patch *pt = &patchBuffer[(*patchCount)++];
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
                if (baseIdx >= 0) {
                    pt->col = voxel_group_color(&voxels[baseIdx]);
                } else {
                    pt->col = (Color){ 255, 255, 255, 255 };
                }
            }
        }
    }
}


static int yzPosList[MAX_VOXELS], yzNegList[MAX_VOXELS];
static int xzPosList[MAX_VOXELS], xzNegList[MAX_VOXELS];
static int xyPosList[MAX_VOXELS], xyNegList[MAX_VOXELS];

static Mesh gen_group_mesh(int groupId, Patch **outPatches, int *outPatchCount) {
    Mesh mesh = (Mesh){ 0 };
    Patch *patchTemp = malloc(sizeof(Patch) * MAX_VOXELS);
    if (!patchTemp) {
        *outPatchCount = 0;
        *outPatches = NULL;
        return mesh;
    }

    int yzPosCount = 0, yzNegCount = 0;
    int xzPosCount = 0, xzNegCount = 0;
    int xyPosCount = 0, xyNegCount = 0;
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (!voxel_is_static_active(v)) continue;
        if (v->groupId != groupId) continue;
        if (v->surface[0]) yzPosList[yzPosCount++] = i;
        if (v->surface[1]) yzNegList[yzNegCount++] = i;
        if (v->surface[2]) xzPosList[xzPosCount++] = i;
        if (v->surface[3]) xzNegList[xzNegCount++] = i;
        if (v->surface[4]) xyPosList[xyPosCount++] = i;
        if (v->surface[5]) xyNegList[xyNegCount++] = i;
    }

    int patchCount = 0;
    merge_rects_on_plane(xyPosCount, xyPosList, 0, true,  patchTemp, MAX_VOXELS, &patchCount);
    merge_rects_on_plane(xyNegCount, xyNegList, 0, false, patchTemp, MAX_VOXELS, &patchCount);
    merge_rects_on_plane(xzPosCount, xzPosList, 1, true,  patchTemp, MAX_VOXELS, &patchCount);
    merge_rects_on_plane(xzNegCount, xzNegList, 1, false, patchTemp, MAX_VOXELS, &patchCount);
    merge_rects_on_plane(yzPosCount, yzPosList, 2, true,  patchTemp, MAX_VOXELS, &patchCount);
    merge_rects_on_plane(yzNegCount, yzNegList, 2, false, patchTemp, MAX_VOXELS, &patchCount);

    if (patchCount == 0) {
        free(patchTemp);
        *outPatchCount = 0;
        *outPatches = NULL;
        return mesh;
    }

    int vCount = patchCount * 4;
    int iCount = patchCount * 6;
    float *verts  = calloc(vCount*3, sizeof(float));
    float *norms  = calloc(vCount*3, sizeof(float));
    float *uvs    = calloc(vCount*2, sizeof(float));
    unsigned char *cols = calloc(vCount*4, sizeof(unsigned char));
    unsigned short *inds = calloc(iCount, sizeof(unsigned short));

    for (int p = 0; p < patchCount; p++) {
        Patch *pt = &patchTemp[p];
        Vector3 origin = {0};
        Vector3 iu, ju, norm;
        switch (pt->plane) {
            case 0:
                origin = (Vector3){ (pt->i0 + 0.0f)*VOXEL_SIZE,
                                    (pt->j0 + 0.0f)*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE*pt->di, 0, 0 };
                ju = (Vector3){ 0, VOXEL_SIZE*pt->dj, 0 };
                norm = pt->positive? (Vector3){0,0,1} : (Vector3){0,0,-1};
                break;
            case 1:
                origin = (Vector3){ (pt->i0 + 0.0f)*VOXEL_SIZE,
                                    (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    (pt->j0 + 0.0f)*VOXEL_SIZE };
                iu = (Vector3){ VOXEL_SIZE*pt->di, 0, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE*pt->dj };
                norm = pt->positive? (Vector3){0,1,0} : (Vector3){0,-1,0};
                break;
            default:
                origin = (Vector3){ (pt->layer + (pt->positive?1:0))*VOXEL_SIZE,
                                    (pt->i0 + 0.0f)*VOXEL_SIZE,
                                    (pt->j0 + 0.0f)*VOXEL_SIZE };
                iu = (Vector3){ 0, VOXEL_SIZE*pt->di, 0 };
                ju = (Vector3){ 0, 0, VOXEL_SIZE*pt->dj };
                norm = pt->positive? (Vector3){1,0,0} : (Vector3){-1,0,0};
                break;
        }

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

    Patch *patchCopy = malloc(sizeof(Patch) * patchCount);
    if (patchCopy) {
        memcpy(patchCopy, patchTemp, sizeof(Patch) * patchCount);
    } else {
        patchCount = 0;
    }

    free(patchTemp);
    *outPatchCount = patchCount;
    *outPatches = patchCopy;
    return mesh;
}

static void ensure_group_mesh(int groupId) {
    if (groupId < 0 || groupId >= voxelGroupCount) return;
    VoxelGroupRender *grp = &voxelGroupRender[groupId];
    if (!grp->dirty) return;
    unload_group_mesh(grp);
    Patch *patchList = NULL;
    int patchCount = 0;
    grp->mesh = gen_group_mesh(groupId, &patchList, &patchCount);
    grp->patches = patchList;
    grp->patchCount = patchCount;
    grp->dirty = false;
    grp->drawOffset = (Vector3){0,0,0};
}

static void draw_group_grid_lines(const VoxelGroupRender *grp) {
    if (!grp->patches) return;
    Vector3 offset = grp->drawOffset;
    for (int p = 0; p < grp->patchCount; p++) {
        const Patch *pt = &grp->patches[p];
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
        origin = v_add(origin, offset);
        for (int ix = 0; ix <= pt->di; ix++) {
            Vector3 a = v_add(origin, v_mul(iu, ix));
            Vector3 b = v_add(a, v_mul(ju, pt->dj));
            rlVertex3f(a.x, a.y, a.z);
            rlVertex3f(b.x, b.y, b.z);
        }
        for (int jy = 0; jy <= pt->dj; jy++) {
            Vector3 a = v_add(origin, v_mul(ju, jy));
            Vector3 b = v_add(a, v_mul(iu, pt->di));
            rlVertex3f(a.x, a.y, a.z);
            rlVertex3f(b.x, b.y, b.z);
        }
    }
}

// Draw all voxels via per-group greedy meshes
static void DrawVoxels(Camera3D cam) {
    (void)cam;
    ensure_voxel_groups_up_to_date();
    static bool materialInit = false;
    static Material groupMaterial;
    if (!materialInit) {
        groupMaterial = LoadMaterialDefault();
        materialInit = true;
    }

    rlDisableBackfaceCulling();
    for (int g = 0; g < voxelGroupCount; g++) {
        ensure_group_mesh(g);
        VoxelGroupRender *grp = &voxelGroupRender[g];
        if (grp->mesh.vertexCount > 0) {
            Matrix transform = MatrixTranslate(grp->drawOffset.x, grp->drawOffset.y, grp->drawOffset.z);
            DrawMesh(grp->mesh, groupMaterial, transform);
        }
    }
    rlBegin(RL_LINES);
    rlColor4ub(0, 0, 0, 255);
    for (int g = 0; g < voxelGroupCount; g++) {
        draw_group_grid_lines(&voxelGroupRender[g]);
    }
    rlEnd();
    rlEnableBackfaceCulling();

    rlBegin(RL_TRIANGLES);
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if (v->simulate) {
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
        update_group_physics(dt);
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
