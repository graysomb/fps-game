// Raylib split-screen FPS prototype (port of fps_game.c)
#include "raylib.h"
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
#define MAX_VOXELS    30000
#define HASH_SIZE     65536    // must be power of two
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
    bool  surface;
    int gx, gy, gz;
} Voxel;
static Voxel voxels[MAX_VOXELS];
static int voxel_count = 0;
// Spatial hash table for voxels
static struct { int key, idx; } table[HASH_SIZE];

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
    int M = (int)(2.0f*FLOOR_SIZE / VOXEL_SIZE);
    for (int x = 0; x <= M ; x++) {
        for (int z = 0; z <= M ; z++) {
            float px = (x + 0.5f) * VOXEL_SIZE-FLOOR_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE-FLOOR_SIZE;
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

// Draw all voxels as cubes
static void DrawVoxels(void) {
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        if(!v->surface) continue;
        DrawCube(v->pos, VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE, v->color);
        DrawCubeWires(v->pos, VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE, BLACK);
    }
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
