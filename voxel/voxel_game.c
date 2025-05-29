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

/* =================== CONFIG =================== */
#define MAX_VOXELS  200
#define MAX_ACTIVE  200
#define HASH_SIZE   32768  /* must be power‑of‑two */
#define GRAVITY     9.81f
#define FIXED_DT    0.0166667f  /* 60 Hz */
#define VOXEL_SIZE  1.0f

/* =================== MATH =================== */
typedef struct { float x, y, z; } Vec3;
static inline Vec3 v3(float x,float y,float z){ return (Vec3){x,y,z}; }
static inline Vec3 v_add(Vec3 a,Vec3 b){ return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 v_mul(Vec3 a,float s){ return v3(a.x*s, a.y*s, a.z*s); }

/* =================== VOXEL =================== */
typedef struct {
    int   gx, gy, gz;   /* grid coords */
    float mass;
    Vec3  vel;
    bool  fixed;
    bool  simulate;
    bool  surface;
    float r,g,b;
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
    voxels[idx]=(Voxel){x,y,z,1.0f,v3(0,0,0),fixed,sim,false,r,g,b};
    hset(x,y,z,idx);
    if(sim && active_count<MAX_ACTIVE) active[active_count++]=idx;
    return idx;
}

/* =================== PHYSICS =================== */
static void physics_step(float dt){
    for(int i=0;i<active_count;){
        int idx=active[i]; Voxel *v=&voxels[idx]; if(!v->simulate){ i++; continue; }
        Vec3 acc=v3(0,-GRAVITY,0);
        v->vel=v_add(v->vel,v_mul(acc,dt));
        Vec3 pos=v_add(v3((v->gx+0.5f)*VOXEL_SIZE,(v->gy+0.5f)*VOXEL_SIZE,(v->gz+0.5f)*VOXEL_SIZE),v_mul(v->vel,dt));
        int nx=floorf(pos.x/VOXEL_SIZE), ny=floorf(pos.y/VOXEL_SIZE), nz=floorf(pos.z/VOXEL_SIZE);
        bool collide=false;
        if(ny<0){ ny=0; v->vel.y=0; collide=true; }
        if(occupied(nx,ny,nz) && !(nx==v->gx&&ny==v->gy&&nz==v->gz)){ v->vel=v3(0,0,0); collide=true; }
        if(!collide){
            hset(v->gx,v->gy,v->gz,-1); /* lazy delete */
            v->gx=nx; v->gy=ny; v->gz=nz; hset(nx,ny,nz,idx); update_around(nx,ny,nz);
        }else if(fabsf(v->vel.x)<0.01f&&fabsf(v->vel.y)<0.01f&&fabsf(v->vel.z)<0.01f){
            v->simulate=false; active[i]=active[--active_count]; continue; }
        i++; }
}

/* =================== RENDERING =================== */
static void draw_voxel(const Voxel *v){
    if(!v->surface) return;
    float x=v->gx*VOXEL_SIZE, y=v->gy*VOXEL_SIZE, z=v->gz*VOXEL_SIZE, s=VOXEL_SIZE;
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
static Vec3 camPos={10,10,20}; static float camYaw=-0.8f;
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
    for(int x=0;x<6;x++) for(int y=0;y<6;y++) for(int z=0;z<6;z++) add_voxel(x,y,z,true,false,0.2f,0.7f,0.9f);
    /* sand */
    for(int i=0;i<200;i++){ int sx=rand()%4+1, sz=rand()%4+1, sy=10+i/16; add_voxel(sx,sy,sz,false,true,0.9f,0.8f,0.2f);} 
    for(int i=0;i<voxel_count;i++) mark_surface(i);
}

/* =================== MAIN =================== */
int main(){ SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
    SDL_Window*win=SDL_CreateWindow("Voxel Demo",SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,800,600,SDL_WINDOW_OPENGL);
    SDL_GLContext ctx=SDL_GL_CreateContext(win); SDL_GL_SetSwapInterval(1);
    glEnable(GL_DEPTH_TEST); glDisable(GL_CULL_FACE);
    build_demo();
    uint64_t prev=SDL_GetPerformanceCounter(); float acc=0;
    bool running=true; while(running){ SDL_Event e; while(SDL_PollEvent(&e)){ if(e.type==SDL_QUIT) running=false; if(e.type==SDL_KEYDOWN&&e.key.keysym.sym==SDLK_ESCAPE) running=false; }
        uint64_t now=SDL_GetPerformanceCounter(); float dt=(now-prev)/(float)SDL_GetPerformanceFrequency(); prev=now;
        handle_input(dt); acc+=dt; while(acc>=FIXED_DT){ physics_step(FIXED_DT); acc-=FIXED_DT; }
        set_camera(); render_scene(); SDL_GL_SwapWindow(win); }
    SDL_GL_DeleteContext(ctx); SDL_DestroyWindow(win); SDL_Quit(); return 0; }
