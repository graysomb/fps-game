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
 *  • Demo scene containing a static cube and a falling sand pile
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
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

/* =================== CONFIG =================== */
#define MAX_VOXELS  200
#define MAX_ACTIVE  200
#define HASH_SIZE   32768  /* must be power‑of‑two */
#define GRAVITY     9.81f
#define FIXED_DT    0.0166667f  /* 60 Hz */
#define VOXEL_SIZE  1.0f
#define INTERACTION_K 10.0f   /* strength of inter-voxel forces */
#define EPSILON 1e-4f         /* minimum distance squared to avoid singularity */
// Limits for physics to prevent runaway acceleration/velocity
#define MAX_ACCELERATION 100.0f  /* maximum acceleration magnitude */
#define MAX_VELOCITY     50.0f   /* maximum velocity magnitude */
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
static int hget(int x,int y,int z){
    int h=hash(x,y,z);
    while(table[h].key){
        int i=table[h].idx;
        if(i>=0){ Voxel *v=&voxels[i]; if(v->gx==x&&v->gy==y&&v->gz==z) return i; }
        h=(h+1)&(HASH_SIZE-1);
    }
    return -1;
}

/* =================== WORLD HELPERS =================== */
static bool occupied(int x,int y,int z){ return hget(x,y,z)>=0; }

static void mark_surface(int idx){
    Voxel *v=&voxels[idx];
    int x=v->gx, y=v->gy, z=v->gz;
    v->surface = !occupied(x+1,y,z)||!occupied(x-1,y,z)||!occupied(x,y+1,z)||!occupied(x,y-1,z)||!occupied(x,y,z+1)||!occupied(x,y,z-1);
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

/* =================== PHYSICS =================== */
static void physics_step(float dt){
    for(int i=0;i<active_count;){
        int idx=active[i]; Voxel *v=&voxels[idx]; if(!v->simulate){ i++; continue; }
        Vec3 acc=v->acc;
        /* inter-voxel attractive/repulsive forces */
        // continuous position for force calculations
        Vec3 pos_v = v->pos;
        for(int j=0;j<voxel_count;j++){
            if(j==idx) continue;
            Voxel *w=&voxels[j];
            if(v->charge==0 || w->charge==0) continue;
            Vec3 pos_w = w->pos;
            Vec3 dpos = v_sub(pos_w, pos_v);
            float dist2 = dpos.x*dpos.x + dpos.y*dpos.y + dpos.z*dpos.z;
            if(dist2 < EPSILON) continue;
            float dist = sqrtf(dist2);
            float force_mag = INTERACTION_K * v->charge * w->charge / dist2;
            Vec3 dir = v_mul(dpos, 1.0f / dist);
            acc = v_add(acc, v_mul(dir, force_mag / v->mass));
        }
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
        bool collide=false;
        /* collision with floor or occupied cell */
        if(ny<0){ ny=0; v->vel.y=0; collide=true; }
        if(occupied(nx,ny,nz) && !(nx==v->gx && ny==v->gy && nz==v->gz)){
            v->vel = v3(0,0,0);
            collide = true;
        }
        if(!collide){
            // only update hash map on actual cell move
            if(nx!=v->gx || ny!=v->gy || nz!=v->gz){
                // remove old mapping
                int ox=v->gx, oy=v->gy, oz=v->gz;
                int h=hash(ox,oy,oz);
                while(table[h].key){
                    if(table[h].idx==idx){ table[h].key=0; break; }
                    h=(h+1)&(HASH_SIZE-1);
                }
                // insert new mapping
                v->gx=nx; v->gy=ny; v->gz=nz;
                hset(nx,ny,nz,idx);
                update_around(nx,ny,nz);
            }
        } else if(fabsf(v->vel.x)<0.01f && fabsf(v->vel.y)<0.01f && fabsf(v->vel.z)<0.01f) {
            /* deactivate slow/stopped voxels */
            v->simulate = false;
            active[i] = active[--active_count];
            continue;
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
static void render_scene(void){ glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); for(int i=0;i<voxel_count;i++) draw_voxel(&voxels[i]); }

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
static void build_demo(void){ /* static cube 6×6×6 */
    /* cube 6×6×6 with random charge and color, enabled for simulation */
    for(int x=0;x<6;x++){
        for(int y=0;y<6;y++){
            for(int z=0;z<6;z++){
                /* make cube voxels simulating so they can move */
                int idx = add_voxel(x, y, z, /*fixed=*/false, /*sim=*/true, 0.2f, 0.7f, 0.9f);
                if(idx<0) continue;
                /* random charge: -1 repel, 0 neutral, 1 attract */
                int ch=(rand()%3)-1;
                voxels[idx].charge=ch;
                /* color based on charge */
                if(ch>0){
                    voxels[idx].r=1.0f; voxels[idx].g=0.0f; voxels[idx].b=0.0f;
                } else if(ch<0){
                    voxels[idx].r=0.0f; voxels[idx].g=0.0f; voxels[idx].b=1.0f;
                } else {
                    voxels[idx].r=0.5f; voxels[idx].g=0.5f; voxels[idx].b=0.5f;
                }
            }
        }
    }
    /* mark all initial surfaces */
    for(int i=0;i<voxel_count;i++){
        mark_surface(i);
    }
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
