#include <SDL2/SDL.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_image.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
//gcc fps_game.c -o fps_game     $(sdl2-config --cflags --libs)     $(pkg-config --cflags --libs SDL2_ttf SDL2_image)     -lGL -lGLU -lm
// Forward declarations for procedural floor


// Constants for game mechanics
#define SCREEN_WIDTH  1600
#define SCREEN_HEIGHT 900
#define MOVE_SPEED    5.0f    // units per second
#define TURN_SPEED    90.0f   // degrees per second
#define JUMP_SPEED    10.0f    // initial upward velocity for jumps
#define GRAVITY       9.8f    // gravity acceleration (units/s^2)
#define PLAYER_RADIUS 0.5f    // radius for player collision (approximate)
#define BASE_EYE_HEIGHT 1.0f  // eye height when standing on ground
#define ACCELERATION   40.0f  // horizontal acceleration (units/s^2)
#define FRICTION       20.0f  // ground friction deceleration (units/s^2)

//constants for map
// Procedural floor parameters
#define FLOOR_SIZE           20.0f     // extent of floor in X/Z (+/- FLOOR_SIZE)

// Player 1 key mappings (WASD movement, F/H yaw, T/G pitch)
const SDL_Scancode P1_FORWARD    = SDL_SCANCODE_W;
const SDL_Scancode P1_BACKWARD   = SDL_SCANCODE_S;
const SDL_Scancode P1_STRAFE_L   = SDL_SCANCODE_A;
const SDL_Scancode P1_STRAFE_R   = SDL_SCANCODE_D;
// F/H for yaw (left/right), T/G for pitch (up/down)
// F/H swapped: H to turn left, F to turn right
const SDL_Scancode P1_TURN_L     = SDL_SCANCODE_H;    // yaw left
const SDL_Scancode P1_TURN_R     = SDL_SCANCODE_F;    // yaw right
const SDL_Scancode P1_PITCH_UP   = SDL_SCANCODE_T;    // look up
const SDL_Scancode P1_PITCH_DOWN = SDL_SCANCODE_G;    // look down
const SDL_Scancode P1_JUMP       = SDL_SCANCODE_SPACE;
const SDL_Scancode P1_SHOOT      = SDL_SCANCODE_LCTRL;

// Player 2 key mappings (IJKL movement)
const SDL_Scancode P2_FORWARD   = SDL_SCANCODE_I;
const SDL_Scancode P2_BACKWARD  = SDL_SCANCODE_K;
// Strafing: J to strafe left, L to strafe right (swapped)
const SDL_Scancode P2_STRAFE_L  = SDL_SCANCODE_J;    // strafe left
const SDL_Scancode P2_STRAFE_R  = SDL_SCANCODE_L;    // strafe right
// Yaw: RIGHT arrow to turn left, LEFT arrow to turn right
const SDL_Scancode P2_TURN_L    = SDL_SCANCODE_RIGHT; // yaw left
const SDL_Scancode P2_TURN_R    = SDL_SCANCODE_LEFT;  // yaw right
// Pitch: UP arrow to look up, DOWN arrow to look down
const SDL_Scancode P2_PITCH_UP  = SDL_SCANCODE_UP;     // look up
const SDL_Scancode P2_PITCH_DOWN= SDL_SCANCODE_DOWN;   // look down
const SDL_Scancode P2_JUMP      = SDL_SCANCODE_RSHIFT;
const SDL_Scancode P2_SHOOT     = SDL_SCANCODE_RCTRL;


int scores[2] = {0, 0}; // scores[0] for player 1, scores[1] for player 2

// Structure definitions
typedef struct {
    float x, y, z;
    float yaw;        // rotation around Y-axis in degrees
    float pitch;      // rotation around X-axis (look up/down) in degrees
    float vy;         // vertical velocity
    float vx, vz;     // horizontal velocity
    bool  onGround;
} Player;

typedef struct {
    float x, y, z;
    float vx, vy, vz;
    bool active;
    int owner;  // which player fired (0 or 1)
    float life; // remaining lifetime in seconds
} Bullet;
// === Voxel physics integration from voxel_game ===
#include <stdint.h>
#include <string.h>
#define MAX_VOXELS 20000
#define HASH_SIZE  32768
#define VOXEL_SIZE 0.2f

typedef struct { float x, y, z; } Vec3;
static inline Vec3 v3(float x, float y, float z) { return (Vec3){x, y, z}; }
static inline Vec3 v_add(Vec3 a, Vec3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline Vec3 v_mul(Vec3 a, float s) { return v3(a.x * s, a.y * s, a.z * s); }

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
static int voxel_count = 0;

typedef struct { int key; int idx; } Slot;
static Slot table[HASH_SIZE];

// Spatial hash for voxels
static inline int hash_voxel(int x, int y, int z) {
    uint32_t h = (uint32_t)(x*73856093 ^ y*19349663 ^ z*83492791);
    return h & (HASH_SIZE - 1);
}
// Insert voxel index into hash
static void table_set(int x, int y, int z, int idx) {
    int h = hash_voxel(x, y, z);
    while (table[h].key) h = (h + 1) & (HASH_SIZE - 1);
    table[h].key = 1;
    table[h].idx = idx;
}
// Retrieve voxel index at grid pos, or -1
static int table_get(int x, int y, int z) {
    int h = hash_voxel(x, y, z);
    while (table[h].key) {
        int i = table[h].idx;
        Voxel *v = &voxels[i];
        if (v->gx == x && v->gy == y && v->gz == z) return i;
        h = (h + 1) & (HASH_SIZE - 1);
    }
    return -1;
}
// Check occupancy
static bool occupied_voxel(int x, int y, int z) {
    return table_get(x, y, z) >= 0;
}

static bool occupied(int x,int y,int z){ return  table_get(x,y,z)>=0; }
// Add a voxel at continuous pos (px,py,pz)
static int add_voxel(float px, float py, float pz, bool fixed, bool sim,
                     float r, float g, float b, int type) {
    if (voxel_count >= MAX_VOXELS) return -1;
    int idx = voxel_count++;
    Voxel *v = &voxels[idx];
    v->pos = v3(px, py, pz);
    v->vel = v3(0.0f, 0.0f, 0.0f);
    v->fixed = fixed;
    v->simulate = sim;
    v->type = type;
    v->surface = true;
    v->r = r; v->g = g; v->b = b;
    int x = (int)floorf(px / VOXEL_SIZE);
    int y = (int)floorf(py / VOXEL_SIZE);
    int z = (int)floorf(pz / VOXEL_SIZE);
    v->gx = x; v->gy = y; v->gz = z;
    table_set(x, y, z, idx);
    return idx;
}

// Simple physics step: gravity & block-block collision
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
                    v->pos = v3(-999.0f, -999.0f, -999.0f);
                    // Remove existing block
                    Voxel *u = &voxels[hit];
                    u->simulate = false;
                    u->fixed = true;
                    u->pos = v3(-999.0f, -999.0f, -999.0f);
                } else {
                    v->simulate = false;
                    v->fixed = true;
                    v->pos = v3((v->gx + 0.5f) * VOXEL_SIZE,
                                (v->gy + 0.5f) * VOXEL_SIZE,
                                (v->gz + 0.5f) * VOXEL_SIZE);
                }
            }
        }
    }
}


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
// Draw all voxels as colored cubes
static void draw_voxels(void) {
    for (int i = 0; i < voxel_count; i++) {
        Voxel *v = &voxels[i];
        draw_voxel(v);
    }
    glColor3f(1.0f, 1.0f, 1.0f);
}

// Build a static cube of voxels like demo
static void build_demo(void) {
    const int N = 10;
    const int offx = 0, offy = 0, offz = 0;
    for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) for (int z = 0; z < N; z++) {
        add_voxel((x + 0.5f) * VOXEL_SIZE + offx,
                  (y + 0.5f) * VOXEL_SIZE + offy,
                  (z + 0.5f) * VOXEL_SIZE + offz,
                  true, false,
                  0.6f, 0.6f, 0.6f,
                  0);
    }
}

// Toggleable voxel type for player 1
static int p1_voxel_type = 0;

// Set a random position within a given range
float randomInRange(float min, float max) {
    return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}


// Call this once at startup to seed the random generator.
void SeedRandom() {
    srand((unsigned)time(NULL));
}

// Global game state
Player players[2];
#define MAX_BULLETS 50
#define BULLET_LIFETIME 3.0f        // bullet time-to-live in seconds
#define BULLET_WORLD_BOUND 20.0f    // world bounds for bullet Y coordinate
Bullet bullets[MAX_BULLETS];
int numBullets = 0;  // or we can reuse bullets in a pool


// Fire a voxel bullet from a player
static void FireVoxel(int playerIndex) {
    float yawRad = players[playerIndex].yaw * (M_PI / 180.0f);
    float pitchRad = players[playerIndex].pitch * (M_PI / 180.0f);
    float dirx = sinf(-yawRad) * cosf(pitchRad);
    float diry = sinf(pitchRad);
    float dirz = -cosf(yawRad) * cosf(pitchRad);
    float offset = 0.8f;
    float sx = players[playerIndex].x + dirx * offset;
    float sy = players[playerIndex].y + diry * offset;
    float sz = players[playerIndex].z + dirz * offset;
    int type = (playerIndex == 0 ? p1_voxel_type : 0);
    // choose color based on type
    float r = (type == 1 ? 1.0f : 1.0f);
    float g = (type == 1 ? 0.0f : 1.0f);
    float b = (type == 1 ? 0.0f : 1.0f);
    int idx = add_voxel(sx, sy, sz, false, true, r, g, b, type);
    if (idx >= 0) {
        voxels[idx].vel = v_mul(v3(dirx, diry, dirz), 20.0f);
    }
}
// Utility: Clamp a value between min and max
static float clamp(float value, float min, float max) {
    if(value < min) return min;
    if(value > max) return max;
    return value;
}


// Reset game (players positions, etc.)
void ResetGame() {
    // Initialize players at different positions
    players[0].x = randomInRange(-9.0f, 9.0f); players[0].y = BASE_EYE_HEIGHT; players[0].z = randomInRange(-9.0f, 9.0f);
    players[0].yaw   = 0.0f; players[0].pitch = 0.0f; players[0].vy = 0.0f; players[0].vx = 0.0f; players[0].vz = 0.0f; players[0].onGround = true;
    players[1].x = randomInRange(-9.0f, 9.0f);  players[1].y = BASE_EYE_HEIGHT; players[1].z = randomInRange(-9.0f, 9.0f);
    players[1].yaw   = 180.0f; players[1].pitch = 0.0f; players[1].vy = 0.0f; players[1].vx = 0.0f; players[1].vz = 0.0f; players[1].onGround = true;
    // Clear bullets
    numBullets = 0;
    // Reset bullets
    for(int i = 0; i < MAX_BULLETS; ++i) {
        bullets[i].active = false;
        bullets[i].life   = 0.0f;
    }

    // Define your overall map area based on procedural floor size
    float mapX = -FLOOR_SIZE;
    float mapZ = -FLOOR_SIZE;
    float mapWidth = FLOOR_SIZE * 2.0f;
    float mapDepth = FLOOR_SIZE * 2.0f;

}

// Spawn a bullet from a player
void FireBullet(int playerIndex) {
    // Find a free bullet slot
    for(int i=0; i<MAX_BULLETS; ++i) {
        if(!bullets[i].active) {
            // Set bullet initial position at player's eye
            // Compute direction with yaw and pitch
            float yawRad   = players[playerIndex].yaw   * (M_PI / 180.0f);
            float pitchRad = players[playerIndex].pitch * (M_PI / 180.0f);
            float dirx = sinf(-yawRad) * cosf(pitchRad);
            float diry = sinf(pitchRad);
            float dirz = -cosf(yawRad) * cosf(pitchRad);
            // Set bullet initial position slightly in front of the eye
            float offset = 0.8f;
            bullets[i].x = players[playerIndex].x + dirx * offset;
            bullets[i].y = players[playerIndex].y + diry * offset;
            bullets[i].z = players[playerIndex].z + dirz * offset;
            // Set bullet velocity
            float speed = 20.0f;
            bullets[i].vx = dirx * speed;
            bullets[i].vy = diry * speed;
            bullets[i].vz = dirz * speed;
            bullets[i].active = true;
            bullets[i].owner  = playerIndex;
            bullets[i].life   = BULLET_LIFETIME;
            break;
        }
    }
}


// Render helper: draw a solid colored cube centered at (cx, cy, cz) with given half-size
void DrawColoredCube(float cx, float cy, float cz, float halfSize) {
    float minx = cx - halfSize, maxx = cx + halfSize;
    float miny = cy - halfSize, maxy = cy + halfSize;
    float minz = cz - halfSize, maxz = cz + halfSize;
    // Each face with a fixed color
    // Front face (z = maxz) - red
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_QUADS);
      glVertex3f(minx, miny, maxz);
      glVertex3f(maxx, miny, maxz);
      glVertex3f(maxx, maxy, maxz);
      glVertex3f(minx, maxy, maxz);
    glEnd();
    // Back face (z = minz) - green
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_QUADS);
      glVertex3f(maxx, miny, minz);
      glVertex3f(minx, miny, minz);
      glVertex3f(minx, maxy, minz);
      glVertex3f(maxx, maxy, minz);
    glEnd();
    // Left face (x = minx) - blue
    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_QUADS);
      glVertex3f(minx, miny, minz);
      glVertex3f(minx, miny, maxz);
      glVertex3f(minx, maxy, maxz);
      glVertex3f(minx, maxy, minz);
    glEnd();
    // Right face (x = maxx) - yellow
    glColor3f(1.0f, 1.0f, 0.0f);
    glBegin(GL_QUADS);
      glVertex3f(maxx, miny, maxz);
      glVertex3f(maxx, miny, minz);
      glVertex3f(maxx, maxy, minz);
      glVertex3f(maxx, maxy, maxz);
    glEnd();
    // Bottom face (y = miny) - cyan
    glColor3f(0.0f, 1.0f, 1.0f);
    glBegin(GL_QUADS);
      glVertex3f(minx, miny, minz);
      glVertex3f(maxx, miny, minz);
      glVertex3f(maxx, miny, maxz);
      glVertex3f(minx, miny, maxz);
    glEnd();
    // Top face (y = maxy) - magenta
    glColor3f(1.0f, 0.0f, 1.0f);
    glBegin(GL_QUADS);
      glVertex3f(minx, maxy, maxz);
      glVertex3f(maxx, maxy, maxz);
      glVertex3f(maxx, maxy, minz);
      glVertex3f(minx, maxy, minz);
    glEnd();
    // Reset color to white (for other objects that might use texture)
    glColor3f(1.0f, 1.0f, 1.0f);
}

void ResolvePlayerCollisions(Player *p) {
    const int maxIterations = 5;
    for (int iter = 0; iter < maxIterations; iter++) {
        bool collisionFound = false;
        
        // Check against world bounds first (optional), using FLOOR_SIZE
        if (p->x < -FLOOR_SIZE + PLAYER_RADIUS) {
            p->x = -FLOOR_SIZE + PLAYER_RADIUS;
            collisionFound = true;
        } else if (p->x > FLOOR_SIZE - PLAYER_RADIUS) {
            p->x = FLOOR_SIZE - PLAYER_RADIUS;
            collisionFound = true;
        }
        if (p->z < -FLOOR_SIZE + PLAYER_RADIUS) {
            p->z = -FLOOR_SIZE + PLAYER_RADIUS;
            collisionFound = true;
        } else if (p->z > FLOOR_SIZE - PLAYER_RADIUS) {
            p->z = FLOOR_SIZE - PLAYER_RADIUS;
            collisionFound = true;
        }
        
        // If no collisions were found, exit early
        if (!collisionFound)
            break;
    }
}


int main(int argc, char *argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return 1;
    }
    // Request OpenGL context (at least 2.1 for fixed-function compatibility)
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    SDL_Window *window = SDL_CreateWindow("Split-Screen FPS Prototype",
                        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                        SCREEN_WIDTH, SCREEN_HEIGHT,
                        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if (!window) {
        fprintf(stderr, "Failed to create window: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_GLContext glContext = SDL_GL_CreateContext(window);
    if (!glContext) {
        fprintf(stderr, "Failed to create OpenGL context: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    SDL_GL_SetSwapInterval(1);  // Enable VSync (optional)

    // OpenGL initial state
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_TEXTURE_2D);
    glShadeModel(GL_SMOOTH);
    // Setup simple lighting if desired (or leave unlit)
    // glEnable(GL_LIGHTING); ... (skipped for simplicity)

    // Setup initial game state
    ResetGame();
    // Build static voxel cube in map
    build_demo();

    // Timing
    Uint32 lastTicks = SDL_GetTicks();

    bool running = true;
    SDL_Event event;
    printf("Starting main loop. Close window or press ESC to quit.\n");
    while (running) {
        // Event handling
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            } else if (event.type == SDL_KEYDOWN) {
                // Handle one-time key presses (shoot, jump)
                if (event.key.keysym.scancode == P1_SHOOT && !event.key.repeat) {
                    FireVoxel(0);
                }
                if (event.key.keysym.scancode == P2_SHOOT && !event.key.repeat) {
                    FireVoxel(1);
                }
                // Toggle voxel type for player 1
                if (event.key.keysym.scancode == SDL_SCANCODE_Q && !event.key.repeat) {
                    p1_voxel_type = 1 - p1_voxel_type;
                }
                if (event.key.keysym.scancode == P1_JUMP && !event.key.repeat) {
                    if (players[0].onGround) {
                        players[0].vy = JUMP_SPEED;
                        players[0].onGround = false;
                    }
                }
                if (event.key.keysym.scancode == P2_JUMP && !event.key.repeat) {
                    if (players[1].onGround) {
                        players[1].vy = JUMP_SPEED;
                        players[1].onGround = false;
                    }
                }
                // Quit on ESC key
                if (event.key.keysym.scancode == SDL_SCANCODE_ESCAPE) {
                    running = false;
                }
            }
        }

        // Calculate delta time
        Uint32 currentTicks = SDL_GetTicks();
        float dt = (currentTicks - lastTicks) / 1000.0f;
        if (dt > 0.1f) dt = 0.1f; // clamp delta (if debugger or freeze, avoid big jumps)
        lastTicks = currentTicks;

        // Get key state for continuous movement
        const Uint8 *keystate = SDL_GetKeyboardState(NULL);

        // Update each player
        for (int i = 0; i < 2; ++i) {
            // Turning (yaw)
            if ((i == 0 && keystate[P1_TURN_L]) || (i == 1 && keystate[P2_TURN_L])) {
                players[i].yaw -= TURN_SPEED * dt;
            }
            if ((i == 0 && keystate[P1_TURN_R]) || (i == 1 && keystate[P2_TURN_R])) {
                players[i].yaw += TURN_SPEED * dt;
            }
            // Keep yaw within [0,360)
            if (players[i].yaw < 0) players[i].yaw += 360.0f;
            if (players[i].yaw >= 360.0f) players[i].yaw -= 360.0f;
            // Pitch control for both players
            if ((i == 0 && keystate[P1_PITCH_UP]) || (i == 1 && keystate[P2_PITCH_UP])) {
                players[i].pitch += TURN_SPEED * dt;
            }
            if ((i == 0 && keystate[P1_PITCH_DOWN]) || (i == 1 && keystate[P2_PITCH_DOWN])) {
                players[i].pitch -= TURN_SPEED * dt;
            }
            // Clamp pitch to avoid flipping
            if (players[i].pitch > 89.0f)  players[i].pitch = 89.0f;
            if (players[i].pitch < -89.0f) players[i].pitch = -89.0f;
            // Movement with momentum and friction
            {
                float yawRad = players[i].yaw * (M_PI / 180.0f);
                float forwardX = sinf(-yawRad);
                float forwardZ = -cosf(yawRad);
                // Strafe direction (perpendicular)
                float rightX = -forwardZ;
                float rightZ =  forwardX;
                // Build input vector
                float inputX = 0.0f, inputZ = 0.0f;
                if ((i == 0 && keystate[P1_FORWARD]) || (i == 1 && keystate[P2_FORWARD])) {
                    inputX += forwardX; inputZ += forwardZ;
                }
                if ((i == 0 && keystate[P1_BACKWARD]) || (i == 1 && keystate[P2_BACKWARD])) {
                    inputX -= forwardX; inputZ -= forwardZ;
                }
                if ((i == 0 && keystate[P1_STRAFE_R]) || (i == 1 && keystate[P2_STRAFE_R])) {
                    inputX += rightX; inputZ += rightZ;
                }
                if ((i == 0 && keystate[P1_STRAFE_L]) || (i == 1 && keystate[P2_STRAFE_L])) {
                    inputX -= rightX; inputZ -= rightZ;
                }
                // Normalize input
                float inMag = sqrtf(inputX*inputX + inputZ*inputZ);
                if (inMag > 0.0f) {
                    inputX /= inMag; inputZ /= inMag;
                }
                // Apply acceleration
                players[i].vx += inputX * ACCELERATION * dt;
                players[i].vz += inputZ * ACCELERATION * dt;
                // Apply friction when no input
                if (inMag < 1e-6f) {
                    float speed = sqrtf(players[i].vx*players[i].vx + players[i].vz*players[i].vz);
                    if (speed > 0.0f) {
                        float decel = FRICTION * dt;
                        float newSpeed = speed - decel;
                        if (newSpeed < 0.0f) newSpeed = 0.0f;
                        players[i].vx *= newSpeed/speed;
                        players[i].vz *= newSpeed/speed;
                    }
                }
                // Clamp max speed
                {
                    float speed = sqrtf(players[i].vx*players[i].vx + players[i].vz*players[i].vz);
                    if (speed > MOVE_SPEED) {
                        players[i].vx *= MOVE_SPEED/speed;
                        players[i].vz *= MOVE_SPEED/speed;
                    }
                }
                // Move and resolve collisions
                players[i].x += players[i].vx * dt;
                players[i].z += players[i].vz * dt;
                ResolvePlayerCollisions(&players[i]);
            }
            
            // Apply gravity and jumping physics for player i
            {
                float oldY = players[i].y;
                // Apply gravity if airborne
                if (!players[i].onGround) {
                    players[i].vy -= GRAVITY * dt;
                }
                // Predict new vertical position
                float newY = oldY + players[i].vy * dt;
                // Determine support height (floor or obstacle tops) beneath player's feet
                float floorY = 0.0f;
                float supportHeight = floorY + BASE_EYE_HEIGHT;
                // Check for landing
                if (newY <= supportHeight) {
                    players[i].y = supportHeight;
                    players[i].vy = 0.0f;
                    players[i].onGround = true;
                } else {
                    players[i].y = newY;
                    players[i].onGround = false;
                }
            }
        } // end player loop

        // Update voxel physics
        physics_step(dt);

        // Rendering
        glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
        glClearColor(0.6f, 0.7f, 0.9f, 1.0f); // a light sky blue color
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Render two views: top (player 1) and bottom (player 2)
        for (int i = 0; i < 2; ++i) {
            if (i == 0) {
                // Top half for player 1
                glViewport(0, SCREEN_HEIGHT/2, SCREEN_WIDTH, SCREEN_HEIGHT/2);
            } else {
                // Bottom half for player 2
                glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT/2);
            }

            // Set perspective projection for this viewport
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            float aspect = (float)SCREEN_WIDTH / ((float)SCREEN_HEIGHT / 2.0f);
            gluPerspective(75.0, aspect, 0.1, 100.0);
            // Set camera view for player i
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            // Camera view with yaw (horizontal) and pitch (vertical)
            {
                float yawRad   = players[i].yaw   * (M_PI / 180.0f);
                float pitchRad = players[i].pitch * (M_PI / 180.0f);
                float dirX = sinf(-yawRad) * cosf(pitchRad);
                float dirY = sinf(pitchRad);
                float dirZ = -cosf(yawRad) * cosf(pitchRad);
                gluLookAt(players[i].x, players[i].y, players[i].z,
                          players[i].x + dirX, players[i].y + dirY, players[i].z + dirZ,
                          0.0f, 1.0f, 0.0f);
            }
            
            // Draw voxels
            glDisable(GL_TEXTURE_2D);

            draw_voxels();

            // Draw other player (cube to represent opponent)
            int other = (i == 0 ? 1 : 0);
            // We draw the opponent only if within view (optional frustum cull – omitted for simplicity)
            // Disable texture, draw a taller cube for player
            glDisable(GL_TEXTURE_2D);
            glColor3f(0.8f, 0.2f, 0.2f); // opponent color (reddish)
            // Draw player model (e.g., cube of height ~1.5 centered at their eye level)
            DrawColoredCube(players[other].x, players[other].y, players[other].z, 0.5f);

             // For example, if you store player scores in an array 'scores':
            char scoreStr[64];
            sprintf(scoreStr, "Score: %d", scores[i]);

        } // end for each viewport

        SDL_GL_SwapWindow(window);
    } // end main loop

    // Cleanup
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
};
