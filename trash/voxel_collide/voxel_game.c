/*
 * voxel_game.c – Minimal voxel world demo on Linux using SDL2 + OpenGL
 *
 * Features
 * --------
 *  • Tunable voxel size (VOXEL_SIZE)
 *  • Voxel struct with mass, gravity, velocity, position, color
 *  • Flags for fixed (static), simulate (active physics), surface (rendered)
 *  • Sparse physics – only voxels in `active` array are stepped each tick
 *  • Simple neighbor queries (hash map) to mark surface voxels
 *  • Player free‑fly camera with WASD (XZ) + T/G (Y) + F/H (yaw) controls
 *  • Demo scene containing a static cube 
 *  • Constant friction damping on all moving voxels
 *
 * Build (Debian/Ubuntu):
 *   sudo apt install build-essential libsdl2-dev libgl1-mesa-dev libglu1-mesa-dev
 *   gcc voxel_game.c -o voxel_game $(sdl2-config --cflags --libs) -lGL -lGLU -lm
 *
 * Run:
 *   ./voxel_game   (ESC quits)
 */

#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <float.h>  // for FLT_MAX
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

/* =================== CONFIG =================== */
#define MAX_VOXELS  20000
#define MAX_ACTIVE  20000
#define HASH_SIZE   32768  /* must be power‑of‑two */
#define GRAVITY     9.81f
#define FIXED_DT    0.0166667f  /* 60 Hz */
#define VOXEL_SIZE  0.1f
#define INTERACTION_K 1.0f   /* strength of inter-voxel forces */
#define EPSILON 1e-3f         /* minimum distance squared to avoid singularity */
// Limits for physics to prevent runaway acceleration/velocity
#define MAX_ACCELERATION 1000.0f  /* maximum acceleration magnitude */
#define MAX_VELOCITY     500.0f   /* maximum velocity magnitude */
#define FRICTION_COEFF   0.1f    /* velocity damping coefficient per second */

/* =================== MATH =================== */
typedef struct { float x, y, z; } Vec3;
static inline Vec3 v3(float x,float y,float z){ return (Vec3){x,y,z}; }
static inline Vec3 v_add(Vec3 a,Vec3 b){ return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 v_mul(Vec3 a,float s){ return v3(a.x*s, a.y*s, a.z*s); }
static inline Vec3 v_sub(Vec3 a,Vec3 b){ return v3(a.x-b.x, a.y-b.y, a.z-b.z); }

/* =================== VOXEL =================== */
typedef struct {
    int   gx, gy, gz;   /* grid coords */
    float mass;
    Vec3  vel;
    Vec3  acc;
    Vec3  pos;         /* continuous center position */
    bool  fixed;
    bool  simulate;
    bool  surface;
    float r, g, b;
    int   charge;      /* -1 repel, 1 attract, 0 none */
} Voxel;

static Voxel voxels[MAX_VOXELS];
static int    voxel_count=0;
static int    active[MAX_ACTIVE];
static int    active_count=0;

/* =================== HASH MAP (open‑address) =================== */
typedef struct { int key; int idx; } Slot;
static Slot table[HASH_SIZE];
static inline int hash(int x,int y,int z){ uint32_t h=(uint32_t)(x*73856093 ^ y*19349663 ^ z*83492791); return h&(HASH_SIZE-1);} 

static void hset(int x,int y,int z,int idx){
    int h=hash(x,y,z);
    while(table[h].key){ h=(h+1)&(HASH_SIZE-1);} /* empty */
    table[h].key=1; table[h].idx=idx;
}
// Retrieve the voxel index at grid position (x,y,z), or -1 if none
static int hget(int x,int y,int z){
    int h = hash(x,y,z);
    while(table[h].key){
        int i = table[h].idx;
        if(i >= 0){
            Voxel *v = &voxels[i];
            if(v->gx == x && v->gy == y && v->gz == z) return i;
        }
        h = (h + 1) & (HASH_SIZE - 1);
    }
    return -1;
}

/* =================== WORLD HELPERS =================== */
static bool occupied(int x,int y,int z){ return hget(x,y,z)>=0; }

static void mark_surface(int idx){
    Voxel *v=&voxels[idx];
    int x=v->gx, y=v->gy, z=v->gz;
    //v->surface = !occupied(x+1,y,z)||!occupied(x-1,y,z)||!occupied(x,y+1,z)||!occupied(x,y-1,z)||!occupied(x,y,z+1)||!occupied(x,y,z-1);
    v->surface = true;
}
static void update_around(int x,int y,int z){
    static const int nb[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    int idx=hget(x,y,z); if(idx>=0) mark_surface(idx);
    for(int i=0;i<6;i++){ idx=hget(x+nb[i][0],y+nb[i][1],z+nb[i][2]); if(idx>=0) mark_surface(idx);} }

static int add_voxel(int x,int y,int z,bool fixed,bool sim,float r,float g,float b){
    if(voxel_count>=MAX_VOXELS) return -1;
    int idx=voxel_count++;
    Vec3 init_acc=sim?v3(0,-GRAVITY,0):v3(0,0,0);
    // initialize voxel properties
    Voxel *v = &voxels[idx];
    v->gx = x; v->gy = y; v->gz = z;
    v->mass = 1.0f;
    v->vel = v3(0,0,0);
    v->acc = init_acc;
    v->pos = v3((x + 0.5f) * VOXEL_SIZE, (y + 0.5f) * VOXEL_SIZE, (z + 0.5f) * VOXEL_SIZE);
    v->fixed = fixed;
    v->simulate = sim;
    v->surface = false;
    v->r = r; v->g = g; v->b = b;
    v->charge = 0;
    // insert into spatial hash
    hset(x, y, z, idx);
    if(sim && active_count<MAX_ACTIVE) active[active_count++]=idx;
    return idx;
}

/* =================== Barnes-Hut Octree =================== */
#define BH_THETA 0.5f
#define BH_MAX_NODES (MAX_VOXELS * 2)

typedef struct {
    Vec3 center;        /* center of this cube */
    float halfSize;     /* half the width of the cube */
    float chargeSum;    /* total charge in this node */
    Vec3 com;           /* center-of-charge (weighted position) */
    int children[8];    /* child node indices, -1 if empty */
    int particleIndex;  /* voxel index if leaf, -1 otherwise */
} BHNode;

static BHNode bhNodes[BH_MAX_NODES];
static int bhNodeCount;

// forward prototypes
static void bh_insertNode(int nodeIndex, int vidx);
static void bh_initTree(void);
static Vec3 bh_computeForce(int nodeIndex, Vec3 pos, int selfIdx);

// Insert a voxel into the octree at nodeIndex
static void bh_insertNode(int nodeIndex, int vidx) {
    if(bhNodeCount >= BH_MAX_NODES) return;
    BHNode *node = &bhNodes[nodeIndex];
    Voxel *v = &voxels[vidx];
    int oct = 0;
    if(v->pos.x >= node->center.x) oct |= 4;
    if(v->pos.y >= node->center.y) oct |= 2;
    if(v->pos.z >= node->center.z) oct |= 1;
    // empty leaf?
    if(node->particleIndex < 0 && node->children[0] < 0) {
        node->particleIndex = vidx;
        node->chargeSum = (float)v->charge;
        node->com = v->pos;
        return;
    }
    // leaf with existing particle: subdivide
    if(node->children[0] < 0) {
        int existing = node->particleIndex;
        node->particleIndex = -1;
        for(int c = 0; c < 8; c++) node->children[c] = -1;
        bh_insertNode(nodeIndex, existing);
    }
    // internal node: ensure child exists
    if(node->children[oct] < 0) {
        int ci = bhNodeCount++;
        BHNode *child = &bhNodes[ci];
        float h = node->halfSize * 0.5f;
        child->halfSize = h;
        child->center = node->center;
        child->center.x += (oct & 4) ? h : -h;
        child->center.y += (oct & 2) ? h : -h;
        child->center.z += (oct & 1) ? h : -h;
        child->chargeSum = 0.0f;
        child->com = v3(0,0,0);
        child->particleIndex = -1;
        for(int c = 0; c < 8; c++) child->children[c] = -1;
        node->children[oct] = ci;
    }
    // insert recursively
    bh_insertNode(node->children[oct], vidx);
    // update aggregate
    node->chargeSum = 0.0f;
    node->com = v3(0,0,0);
    for(int c = 0; c < 8; c++){
        int ci = node->children[c];
        if(ci < 0) continue;
        BHNode *ch = &bhNodes[ci];
        if(ch->chargeSum != 0.0f){
            node->chargeSum += ch->chargeSum;
            node->com = v_add(node->com, v_mul(ch->com, ch->chargeSum));
        }
    }
    if(node->chargeSum != 0.0f) {
        node->com = v_mul(node->com, 1.0f / node->chargeSum);
    }
}

// Build the spatial octree from all voxels (for collision queries)
static void bh_initTree(void) {
    float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
    float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
    for(int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        Vec3 p = v->pos;
        if(p.x < minx) minx = p.x; if(p.y < miny) miny = p.y; if(p.z < minz) minz = p.z;
        if(p.x > maxx) maxx = p.x; if(p.y > maxy) maxy = p.y; if(p.z > maxz) maxz = p.z;
    }
    Vec3 center = v3((minx + maxx) * 0.5f, (miny + maxy) * 0.5f, (minz + maxz) * 0.5f);
    float half = fmaxf(fmaxf(maxx - minx, maxy - miny), maxz - minz) * 0.5f;
    if(half < 1e-3f) half = 1.0f;
    bhNodeCount = 1;
    BHNode *root = &bhNodes[0];
    root->center = center;
    root->halfSize = half;
    root->chargeSum = 0.0f;
    root->com = v3(0,0,0);
    root->particleIndex = -1;
    for(int c = 0; c < 8; c++) root->children[c] = -1;
    for(int i = 0; i < voxel_count; i++){
        // insert every voxel into the spatial tree
        bh_insertNode(0, i);
    }
}

// Traverse the tree to compute approximate Coulomb acceleration
static Vec3 bh_computeForce(int nodeIndex, Vec3 pos, int selfIdx) {
    BHNode *node = &bhNodes[nodeIndex];
    Vec3 force = v3(0,0,0);
    if(node->particleIndex == selfIdx && node->children[0] < 0) {
        return force;
    }
    Vec3 dpos = v_sub(node->com, pos);
    float dist2 = dpos.x*dpos.x + dpos.y*dpos.y + dpos.z*dpos.z;
    if(dist2 < EPSILON) return force;
    float dist = sqrtf(dist2);
    // leaf?
    if(node->children[0] < 0) {
        Voxel *v = &voxels[selfIdx];
        float Q = node->chargeSum;
        float fmag = INTERACTION_K * v->charge * Q / dist2;
        Vec3 dir = v_mul(dpos, 1.0f / dist);
        return v_mul(dir, fmag / v->mass);
    }
    // opening angle criterion
    float s = node->halfSize * 2.0f;
    if(s / dist < BH_THETA) {
        Voxel *v = &voxels[selfIdx];
        float Q = node->chargeSum;
        float fmag = INTERACTION_K * v->charge * Q / dist2;
        Vec3 dir = v_mul(dpos, 1.0f / dist);
        return v_mul(dir, fmag / v->mass);
    }
    // recurse
    for(int c = 0; c < 8; c++){
        int ci = node->children[c];
        if(ci >= 0) force = v_add(force, bh_computeForce(ci, pos, selfIdx));
    }
    return force;
}

// AABB intersection test between cubic node and query box
static bool aabb_intersect(const Vec3 center, float halfSize, Vec3 min, Vec3 max) {
    Vec3 nodeMin = v3(center.x - halfSize, center.y - halfSize, center.z - halfSize);
    Vec3 nodeMax = v3(center.x + halfSize, center.y + halfSize, center.z + halfSize);
    if (nodeMin.x > max.x || nodeMax.x < min.x) return false;
    if (nodeMin.y > max.y || nodeMax.y < min.y) return false;
    if (nodeMin.z > max.z || nodeMax.z < min.z) return false;
    return true;
}

// Query the octree for voxels whose centers lie within [min, max]
static void bh_queryRange(int nodeIndex, Vec3 min, Vec3 max, int *out, int *count) {
    BHNode *node = &bhNodes[nodeIndex];
    if (!aabb_intersect(node->center, node->halfSize, min, max)) return;
    // Leaf node with a particle
    if (node->children[0] < 0 && node->particleIndex >= 0) {
        Voxel *v = &voxels[node->particleIndex];
        if (v->pos.x >= min.x && v->pos.x <= max.x &&
            v->pos.y >= min.y && v->pos.y <= max.y &&
            v->pos.z >= min.z && v->pos.z <= max.z) {
            out[(*count)++] = node->particleIndex;
        }
        return;
    }
    // Recurse into children
    for (int i = 0; i < 8; i++) {
        int ci = node->children[i];
        if (ci >= 0) bh_queryRange(ci, min, max, out, count);
    }
}

/* =================== PHYSICS =================== */
static void physics_step(float dt){
    // Build spatial octree for collision broad-phase
    bh_initTree();
    for(int i=0;i<active_count;){
        int idx = active[i];
        Voxel *v = &voxels[idx];
        // skip non-simulating or fixed voxels
        if(!v->simulate || v->fixed){
            i++;
            continue;
        }
        // start with gravity only
        Vec3 acc = v3(0.0f, -GRAVITY, 0.0f);
        // clamp acceleration to prevent excessive forces
        {
            float acc_mag2 = acc.x*acc.x + acc.y*acc.y + acc.z*acc.z;
            if(acc_mag2 > MAX_ACCELERATION * MAX_ACCELERATION) {
                float acc_mag = sqrtf(acc_mag2);
                acc = v_mul(acc, MAX_ACCELERATION / acc_mag);
            }
        }
        v->vel = v_add(v->vel, v_mul(acc, dt));
        // clamp velocity to fixed maximum
        {
            float vel_mag2 = v->vel.x*v->vel.x + v->vel.y*v->vel.y + v->vel.z*v->vel.z;
            if(vel_mag2 > MAX_VELOCITY * MAX_VELOCITY) {
                float vel_mag = sqrtf(vel_mag2);
                v->vel = v_mul(v->vel, MAX_VELOCITY / vel_mag);
            }
        }
        // apply friction damping
        {
            float friction_factor = 1.0f - FRICTION_COEFF * dt;
            if(friction_factor < 0.0f) friction_factor = 0.0f;
            v->vel = v_mul(v->vel, friction_factor);
        }
        // update continuous position
        v->pos = v_add(v->pos, v_mul(v->vel, dt));
        // determine new grid cell based on position
        int nx = (int)floorf(v->pos.x / VOXEL_SIZE);
        int ny = (int)floorf(v->pos.y / VOXEL_SIZE);
        int nz = (int)floorf(v->pos.z / VOXEL_SIZE);
        /* floor collision: stop on the floor at y=0 (no bounce) */
        {
            const float half = 0.5f * VOXEL_SIZE;
            if(v->pos.y < half) {
                v->pos.y = half;
                // transfer vertical momentum to the floor (third-law reaction)
                v->vel.y = 0.0f;
                acc.y    = 0.0f;
            }
        }
        /* block-block collision: equal and opposite impulse transfer */
        if(nx != v->gx || ny != v->gy || nz != v->gz) {
            // Broad-phase via spatial octree query instead of direct hash lookup
            Vec3 qmin = v3(nx * VOXEL_SIZE, ny * VOXEL_SIZE, nz * VOXEL_SIZE);
            Vec3 qmax = v3((nx + 1) * VOXEL_SIZE, (ny + 1) * VOXEL_SIZE, (nz + 1) * VOXEL_SIZE);
            int cands[8], ccount = 0;
            bh_queryRange(0, qmin, qmax, cands, &ccount);
            int oidx = (ccount > 0) ? cands[0] : -1;
            if(oidx >= 0) {
                Voxel *u = &voxels[oidx];
                // break fixed voxel on impact: enable physics on it
                if(u->fixed) {
                    u->fixed = false;
                    u->simulate = true;
                    if(active_count < MAX_ACTIVE) active[active_count++] = oidx;
                    mark_surface(oidx);
                    update_around(u->gx, u->gy, u->gz);
                }
                /* compute impulse J = m1 * v1 */
                Vec3 v1 = v->vel;
                float m1 = v->mass;
                float m2 = u->mass;
                Vec3 J = v_mul(v1, m1);
                /* apply to occupant: Δv2 = J/m2 */
                if(!u->fixed) u->vel = v_add(u->vel, v_mul(J, 1.0f / m2));
                /* apply to self:   Δv1 = -J/m1 */
                v->vel = v_add(v->vel, v_mul(J, -1.0f / m1));
                /* snap back to current cell center */
                v->pos = v3((v->gx + 0.5f) * VOXEL_SIZE,
                           (v->gy + 0.5f) * VOXEL_SIZE,
                           (v->gz + 0.5f) * VOXEL_SIZE);
                continue;
            }
            /* free to move: update hash map */
            int ox=v->gx, oy=v->gy, oz=v->gz;
            int h=hash(ox,oy,oz);
            while(table[h].key) {
                if(table[h].idx==idx) { table[h].key=0; break; }
                h=(h+1)&(HASH_SIZE-1);
            }
            v->gx=nx; v->gy=ny; v->gz=nz;
            hset(nx,ny,nz,idx);
            update_around(nx,ny,nz);
        }
        i++;
    }
}

/* =================== RENDERING =================== */
static void draw_voxel(const Voxel *v){
    if(!v->surface) return;
    // compute cube corner from continuous center position
    float s = VOXEL_SIZE;
    float x = v->pos.x - 0.5f * s;
    float y = v->pos.y - 0.5f * s;
    float z = v->pos.z - 0.5f * s;
    glColor3f(v->r,v->g,v->b);
    if(!occupied(v->gx+1,v->gy,v->gz)){ glBegin(GL_QUADS); glVertex3f(x+s,y  ,z  ); glVertex3f(x+s,y  ,z+s); glVertex3f(x+s,y+s,z+s); glVertex3f(x+s,y+s,z  ); glEnd(); }
    if(!occupied(v->gx-1,v->gy,v->gz)){ glBegin(GL_QUADS); glVertex3f(x  ,y  ,z  ); glVertex3f(x  ,y+s,z  ); glVertex3f(x  ,y+s,z+s); glVertex3f(x  ,y  ,z+s); glEnd(); }
    if(!occupied(v->gx,v->gy+1,v->gz)){ glBegin(GL_QUADS); glVertex3f(x  ,y+s,z  ); glVertex3f(x+s,y+s,z  ); glVertex3f(x+s,y+s,z+s); glVertex3f(x  ,y+s,z+s); glEnd(); }
    if(!occupied(v->gx,v->gy-1,v->gz)){ glBegin(GL_QUADS); glVertex3f(x  ,y  ,z  ); glVertex3f(x  ,y  ,z+s); glVertex3f(x+s,y  ,z+s); glVertex3f(x+s,y  ,z  ); glEnd(); }
    if(!occupied(v->gx,v->gy,v->gz+1)){ glBegin(GL_QUADS); glVertex3f(x  ,y  ,z+s); glVertex3f(x  ,y+s,z+s); glVertex3f(x+s,y+s,z+s); glVertex3f(x+s,y  ,z+s); glEnd(); }
    if(!occupied(v->gx,v->gy,v->gz-1)){ glBegin(GL_QUADS); glVertex3f(x  ,y  ,z  ); glVertex3f(x+s,y  ,z  ); glVertex3f(x+s,y+s,z  ); glVertex3f(x  ,y+s,z  ); glEnd(); }
}
static void draw_bh_boxes(void) {
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(1.0f);
    glDisable(GL_DEPTH_TEST);
    glBegin(GL_LINES);
    for(int i = 0; i < bhNodeCount; i++){
        BHNode *n = &bhNodes[i];
        float minx = n->center.x - n->halfSize;
        float maxx = n->center.x + n->halfSize;
        float miny = n->center.y - n->halfSize;
        float maxy = n->center.y + n->halfSize;
        float minz = n->center.z - n->halfSize;
        float maxz = n->center.z + n->halfSize;
        glVertex3f(minx,miny,minz); glVertex3f(maxx,miny,minz);
        glVertex3f(maxx,miny,maxz); glVertex3f(minx,miny,maxz);
        glVertex3f(minx,maxy,minz); glVertex3f(maxx,maxy,minz);
        glVertex3f(maxx,maxy,maxz); glVertex3f(minx,maxy,maxz);
        glVertex3f(minx,miny,minz); glVertex3f(minx,maxy,minz);
        glVertex3f(maxx,miny,minz); glVertex3f(maxx,maxy,minz);
        glVertex3f(maxx,miny,maxz); glVertex3f(maxx,maxy,maxz);
        glVertex3f(minx,miny,maxz); glVertex3f(minx,maxy,maxz);

        // bottom face, front→back
        glVertex3f(minx, miny, minz);  
        glVertex3f(minx, miny, maxz);

        glVertex3f(maxx, miny, minz);  
        glVertex3f(maxx, miny, maxz);

        // top face, front→back
        glVertex3f(minx, maxy, minz);  
        glVertex3f(minx, maxy, maxz);

        glVertex3f(maxx, maxy, minz);  
        glVertex3f(maxx, maxy, maxz);
    }
    glEnd();
    glEnable(GL_DEPTH_TEST);
}
static void render_scene(void){ glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); for(int i=0;i<voxel_count;i++) draw_voxel(&voxels[i]); draw_bh_boxes(); }

/* =================== CAMERA & INPUT =================== */
/* position camera to face demo cube at startup */
static Vec3 camPos={3.0f,3.0f,15.0f};
static float camYaw=-1.5708f;  /* -90° yaw to look down -Z towards the cube center */
static void handle_input(float dt){ const Uint8*ks=SDL_GetKeyboardState(NULL); float mv=10.0f, rt=1.5f; Vec3 fwd=v3(cosf(camYaw),0,sinf(camYaw)); Vec3 right=v3(-fwd.z,0,fwd.x);
    if(ks[SDL_SCANCODE_W]) camPos=v_add(camPos,v_mul(fwd,mv*dt));
    if(ks[SDL_SCANCODE_S]) camPos=v_add(camPos,v_mul(fwd,-mv*dt));
    if(ks[SDL_SCANCODE_A]) camPos=v_add(camPos,v_mul(right,-mv*dt));
    if(ks[SDL_SCANCODE_D]) camPos=v_add(camPos,v_mul(right, mv*dt));
    if(ks[SDL_SCANCODE_T]) camPos.y+=mv*dt; if(ks[SDL_SCANCODE_G]) camPos.y-=mv*dt;
    if(ks[SDL_SCANCODE_F]) camYaw+=rt*dt;   if(ks[SDL_SCANCODE_H]) camYaw-=rt*dt; }
static void set_camera(void){ glMatrixMode(GL_PROJECTION); glLoadIdentity(); gluPerspective(60,4.0/3.0,0.1,1000); glMatrixMode(GL_MODELVIEW); glLoadIdentity(); Vec3 look=v_add(camPos,v3(cosf(camYaw),0,sinf(camYaw))); gluLookAt(camPos.x,camPos.y,camPos.z,look.x,look.y,look.z,0,1,0);} 

/* =================== DEMO SCENE =================== */
static void build_demo(void) {
    // Large static block that will break apart on impact
    const int N2 = 10;
    const int offset2x = 0, offset2y = 0, offset2z = 0;
    for(int x = 0; x < N2; x++) {
        for(int y = 0; y < N2; y++) {
            for(int z = 0; z < N2; z++) {
                add_voxel(x + offset2x, y + offset2y, z + offset2z,
                          /*fixed=*/true, /*sim=*/false,
                          0.6f, 0.6f, 0.6f);
            }
        }
    }
    // Moving block that will be propelled into the static block
    const int N1 = 8;
    const int start_z = offset2z + N2 + 500;
    Vec3 init_vel = v3(0.0f, 5.5f, -50.0f);
    for(int x = 0; x < N1; x++) {
        for(int y = 0; y < N1; y++) {
            for(int z = 0; z < N1; z++) {
                int idx = add_voxel(offset2x + (N2 - N1) / 2 + x,
                                    offset2y + (N2 - N1) / 2 + y,
                                    start_z + z,
                                    /*fixed=*/false, /*sim=*/true,
                                    1.0f, 0.2f, 0.2f);
                voxels[idx].vel = init_vel;
            }
        }
    }
    // mark all initial surfaces for rendering
    for(int i = 0; i < voxel_count; i++) mark_surface(i);
}

/* =================== MAIN =================== */
int main(){ SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
    SDL_Window*win=SDL_CreateWindow("Voxel Demo",SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,800,600,SDL_WINDOW_OPENGL);
    SDL_GLContext ctx=SDL_GL_CreateContext(win); SDL_GL_SetSwapInterval(1);
    glEnable(GL_DEPTH_TEST); glDisable(GL_CULL_FACE);
    build_demo();
    uint64_t prev=SDL_GetPerformanceCounter(); float acc=0;
    bool running=true; while(running){ 
        SDL_Event e; 
        while(SDL_PollEvent(&e)){ 
            if(e.type==SDL_QUIT) running=false; 
            if(e.type==SDL_KEYDOWN&&e.key.keysym.sym==SDLK_ESCAPE) running=false; }
        uint64_t now=SDL_GetPerformanceCounter();
        float dt=(now-prev)/(float) SDL_GetPerformanceFrequency();
        prev=now;
        handle_input(dt);
        acc+=dt;
        //float out = acc/dt;
        //printf( "Now = %f ticks\n", (float) out );
        while(acc>=FIXED_DT){
            physics_step(FIXED_DT);
            acc -= FIXED_DT;
        }
        set_camera();
        render_scene();
        SDL_GL_SwapWindow(win); }
    SDL_GL_DeleteContext(ctx); SDL_DestroyWindow(win); SDL_Quit(); return 0; }
