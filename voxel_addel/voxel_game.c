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
#define VOXEL_SIZE  0.2f
#define INTERACTION_K 1.0f   /* strength of inter-voxel forces */
#define EPSILON 1e-3f         /* minimum distance squared to avoid singularity */
// Limits for physics to prevent runaway acceleration/velocity
#define MAX_ACCELERATION 1000.0f  /* maximum acceleration magnitude */
#define MAX_VELOCITY     500.0f   /* maximum velocity magnitude */
#define FRICTION_COEFF   0.1f    /* velocity damping coefficient per second */
#define STICK   0.1f

/* =================== MATH =================== */
typedef struct { float x, y, z; } Vec3;
static inline Vec3 v3(float x,float y,float z){ return (Vec3){x,y,z}; }
static inline Vec3 v_add(Vec3 a,Vec3 b){ return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 v_mul(Vec3 a,float s){ return v3(a.x*s, a.y*s, a.z*s); }
static inline Vec3 v_sub(Vec3 a,Vec3 b){ return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline float v_prod(Vec3 a,Vec3 b){ return (a.x*b.x+ a.y*b.y + a.z*b.z); }


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
    int type;
    float r, g, b;
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

static int add_voxel(int x,int y,int z,bool fixed,bool sim,float r,float g,float b ,int type){
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
    v->type = type;
    v->r = r; v->g = g; v->b = b;
    // insert into spatial hash
    hset(x, y, z, idx);
    if(sim && active_count<MAX_ACTIVE) active[active_count++]=idx;
    return idx;
}

/* =================== PHYSICS =================== */
static void physics_step(float dt){
    // clear hash table
    memset(table, 0, sizeof(table));

    //rebuild
    for (int i = 0; i < voxel_count; i++) {
    Voxel* v = &voxels[i];
    int gx = (int)(v->pos.x / VOXEL_SIZE);
    int gy = (int)(v->pos.y / VOXEL_SIZE);
    int gz = (int)(v->pos.z / VOXEL_SIZE);
    v->gx = gx; v->gy = gy; v->gz = gz;
    hset(gx, gy, gz, i);
    }

    //iterate over voxels
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
        if (v_prod(v->vel,v->vel)<(10.0f*STICK)){
            v->vel = v3(0,0,0);
            v->simulate =false;
            v->fixed = true;
        }
        //update velocity
        v->vel = v_add(v->vel, v_mul(acc, dt));

        v->vel = v_mul(v->vel, 0.999f);
        // update continuous position
        v->pos = v_add(v->pos, v_mul(v->vel, dt));
        // determine new grid cell based on position
        int nx = (int)floorf(v->pos.x / VOXEL_SIZE);
        int ny = (int)floorf(v->pos.y / VOXEL_SIZE);
        int nz = (int)floorf(v->pos.z / VOXEL_SIZE);
        /* floor collision: stop on the floor at y=0 (no bounce) */
        {
            const float half = 0.5f * VOXEL_SIZE;
            if(v->pos.y < half ) {
                v->pos.y = half;
                // rebound with half velocity
                v->vel.y = -0.5f * v->vel.y;
            }
        }
        /* block-block collision: equal and opposite impulse transfer */
        if((nx != v->gx || ny != v->gy || nz != v->gz) && occupied(nx,ny,nz)) {
            // break fixed voxel on impact: enable physics on it
            int hit_index = hget(nx, ny, nz);
            if (hit_index >= 0) {
                Voxel* u = &voxels[hit_index];
                //delete both blocks apparently resizing is expensive
                if (v->type == 1){
                    // Remove both voxels by marking their hash entries as empty and disabling them
                    table[hash(v->gx, v->gy, v->gz)].key = 0;
                    table[hash(u->gx, u->gy, u->gz)].key = 0;

                    // Move them out of bounds and deactivate
                    v->simulate = false; v->fixed = true; v->pos = v3(-999, -999, -999);
                    u->simulate = false; u->fixed = true; u->pos = v3(-999, -999, -999);

                    // Optional: mark them as not surface
                    v->surface = false;
                    u->surface = false;

                    mark_surface(hit_index);
                    update_around(u->gx, u->gy, u->gz);

                    // Remove from active list (v first, then u)
                    for (int j = 0; j < active_count; j++) {
                        if (active[j] == idx) {
                            active[j] = active[--active_count];
                            break;
                        }
                    }
                    for (int j = 0; j < active_count; j++) {
                        if (active[j] == hit_index) {
                            active[j] = active[--active_count];
                            break;
                        }
                    }
                }else{
                    //stop particle and snap to grid
                    v->simulate = false; v->fixed = true;
                    v->pos = v3((v->gx + 0.5f) * VOXEL_SIZE,
                           (v->gy + 0.5f) * VOXEL_SIZE,
                           (v->gz + 0.5f) * VOXEL_SIZE);
                    mark_surface(hit_index);
                    update_around(u->gx, u->gy, u->gz);
                }

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

    float s = VOXEL_SIZE;
    float x = v->pos.x - 0.5f * s;
    float y = v->pos.y - 0.5f * s;
    float z = v->pos.z - 0.5f * s;

    typedef struct { float x, y, z; } Pt;
    Pt face[6][4] = {
        {{x+s,y  ,z  },{x+s,y  ,z+s},{x+s,y+s,z+s},{x+s,y+s,z  }}, // +X
        {{x  ,y  ,z  },{x  ,y+s,z  },{x  ,y+s,z+s},{x  ,y  ,z+s}}, // -X
        {{x  ,y+s,z  },{x+s,y+s,z  },{x+s,y+s,z+s},{x  ,y+s,z+s}}, // +Y
        {{x  ,y  ,z  },{x  ,y  ,z+s},{x+s,y  ,z+s},{x+s,y  ,z  }}, // -Y
        {{x  ,y  ,z+s},{x  ,y+s,z+s},{x+s,y+s,z+s},{x+s,y  ,z+s}}, // +Z
        {{x  ,y  ,z  },{x+s,y  ,z  },{x+s,y+s,z  },{x  ,y+s,z  }}  // -Z
    };
    int gx = v->gx, gy = v->gy, gz = v->gz;
    int neighbors[6][3] = {
        {gx+1, gy, gz}, {gx-1, gy, gz}, {gx, gy+1, gz},
        {gx, gy-1, gz}, {gx, gy, gz+1}, {gx, gy, gz-1}
    };

    for(int i = 0; i < 6; i++){
        if(occupied(neighbors[i][0], neighbors[i][1], neighbors[i][2]))
            continue;

        // Draw face color
        glColor3f(v->r, v->g, v->b);
        glBegin(GL_QUADS);
        for(int j = 0; j < 4; j++)
            glVertex3f(face[i][j].x, face[i][j].y, face[i][j].z);
        glEnd();

        // Draw black wireframe edge
        glColor3f(0, 0, 0);
        glBegin(GL_LINE_LOOP);
        for(int j = 0; j < 4; j++)
            glVertex3f(face[i][j].x, face[i][j].y, face[i][j].z);
        glEnd();
    }
}

static void render_scene(void){ glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); for(int i=0;i<voxel_count;i++) draw_voxel(&voxels[i]); }

/* =================== CAMERA & INPUT =================== */
/* position camera to face demo cube at startup */
static Vec3 camPos={3.0f,3.0f,10.0f};
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
                          0.6f, 0.6f, 0.6f, 0);
            }
        }
    }
    // Moving block that will be propelled into the static block
    const int N1 = 1;
    const int start_z = offset2z + N2 + 50;
    Vec3 init_vel = v3(0.0f, 7.0f, -5.0f);
    for(int x = 0; x < N1; x++) {
        for(int y = 0; y < N1; y++) {
            for(int z = 0; z < N1; z++) {
                int idx = add_voxel(offset2x + (N2 - N1) / 2 + x,
                                    offset2y + (N2 - N1) / 2 + y,
                                    start_z + z,
                                    /*fixed=*/false, /*sim=*/true,
                                    1.0f, 0.2f, 0.2f, 1);
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
