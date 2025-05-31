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


// Utility: Clamp a value between min and max
static float clamp(float value, float min, float max) {
    if(value < min) return min;
    if(value > max) return max;
    return value;
}

// Check collision between a point (px, pz) (player position on XZ plane) 
// and an AABB obstacle (projected on XZ). Uses player radius for a soft collision boundary.
bool CollidePlayerWithBox(float px, float pz, float radius, Box box) {
    // Find closest point on the box to the player's center (xz only, treat y as irrelevant for walls)
    float closestX = clamp(px, box.minx, box.maxx);
    float closestZ = clamp(pz, box.minz, box.maxz);
    // Distance from player to this point
    float dx = px - closestX;
    float dz = pz - closestZ;
    float distSq = dx*dx + dz*dz;
    return distSq < (radius * radius);
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
            // Ensure bullet spawns above the floor
            {
                float floorY = GetFloorHeight(bullets[i].x, bullets[i].z);
                if (bullets[i].y <= floorY + 0.05f) {
                    bullets[i].y = floorY + 0.05f;
                }
            }
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

// Resolves collision between a player and a box by pushing the player out
void ResolvePlayerCollision(Player *p, Box box) {
    // Compute the closest point on the box (projected to the XZ plane)
    float closestX = clamp(p->x, box.minx, box.maxx);
    float closestZ = clamp(p->z, box.minz, box.maxz);
    
    // Compute vector from the closest point to the player
    float dx = p->x - closestX;
    float dz = p->z - closestZ;
    
    // Compute distance (in the XZ plane)
    float dist = sqrtf(dx * dx + dz * dz);
    
    // If the distance is less than the player's radius, there's penetration.
    if (dist < PLAYER_RADIUS) {
        float penetration = PLAYER_RADIUS - dist;
        // Avoid division by zero; if inside exactly, choose an arbitrary direction.
        if (dist == 0.0f) {
            dx = 1.0f;
            dz = 0.0f;
            dist = 1.0f;
        }
        // Push the player out along the collision normal
        p->x += (dx / dist) * penetration*2;
        p->z += (dz / dist) * penetration*2;
    }
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
        
        // Check against all collision obstacles
        for (int o = 0; o < numCollisionObstacles; ++o) {
            Box *b = &collisionObstacles[o];
            // Only collide if player is below the top of the obstacle (allow jumping over)
            if (p->y < b->maxy) {
                if (CollidePlayerWithBox(p->x, p->z, PLAYER_RADIUS, *b)) {
                    ResolvePlayerCollision(p, *b);
                    collisionFound = true;
                }
            }
        }
        // If no collisions were found, exit early
        if (!collisionFound)
            break;
    }
}

TTF_Font* font = NULL;
GLuint scoreTexture = 3;

// Initialize SDL_ttf and load a font
void InitFont() {
    TTF_Init();
    font = TTF_OpenFont("./leadcoat.ttf", 24);  // Use an appropriate font and size
    if (!font) {
        fprintf(stderr, "Failed to load font: %s\n", TTF_GetError());
    }
}

// Render a score string to an OpenGL texture
void CreateScoreTexture(const char* text, GLuint *texID) {
    SDL_Color white = {255, 255, 255};
    SDL_Surface* surface = TTF_RenderText_Blended(font, text, white);
    if (!surface) {
        fprintf(stderr, "TTF_RenderText_Blended error: %s\n", TTF_GetError());
        return;
    }
    
    // Convert the surface to a known 32-bit RGBA format.
    SDL_Surface* formatted = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGBA32, 0);
    SDL_FreeSurface(surface);  // free the original surface
    if (!formatted) {
        fprintf(stderr, "SDL_ConvertSurfaceFormat error: %s\n", SDL_GetError());
        return;
    }
    
    // Generate and bind texture.
    glGenTextures(1, texID);
    glBindTexture(GL_TEXTURE_2D, *texID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // Ensure rows are tightly packed.
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    // Upload texture data.
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, formatted->w, formatted->h,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, formatted->pixels);
    
    // Optionally, you can store the dimensions in global variables so you know how large your texture is.
    // For example: texWidth = formatted->w; texHeight = formatted->h;
    
    SDL_FreeSurface(formatted);
}


void RenderScoreTexture(GLuint texID, int texWidth, int texHeight) {
    // Set up orthographic projection
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT/2);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    // Enable blending so the alpha channel is used
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texID);
    glColor3f(1, 1, 1);
    
    glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2f(10, 10);
        glTexCoord2f(1, 0); glVertex2f(10 + texWidth, 10);
        glTexCoord2f(1, 1); glVertex2f(10 + texWidth, 10 + texHeight);
        glTexCoord2f(0, 1); glVertex2f(10, 10 + texHeight);
    glEnd();
    
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
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


    // Initialization code (before the main loop)
    InitFont();  // Initialize SDL_ttf and load the font
    // Create the initial score texture, for example with score "Score: 0"
    CreateScoreTexture("Score: 0", &scoreTexture);

    // Setup initial game state
    ResetGame();

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
                    FireBullet(0);
                } 
                if (event.key.keysym.scancode == P2_SHOOT && !event.key.repeat) {
                    FireBullet(1);
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
                float floorY = GetFloorHeight(players[i].x, players[i].z);
                float supportHeight = floorY + BASE_EYE_HEIGHT;
                for (int o = 0; o < numCollisionObstacles; ++o) {
                    Box *b = &collisionObstacles[o];
                    if (players[i].x >= b->minx - PLAYER_RADIUS && players[i].x <= b->maxx + PLAYER_RADIUS &&
                        players[i].z >= b->minz - PLAYER_RADIUS && players[i].z <= b->maxz + PLAYER_RADIUS) {
                        supportHeight = fmaxf(supportHeight, b->maxy);
                    }
                }
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

        // Update bullets
        for (int b = 0; b < MAX_BULLETS; ++b) {
            if (!bullets[b].active) continue;
            // Decrease lifetime and deactivate if expired
            bullets[b].life -= dt;
            if (bullets[b].life <= 0.0f) {
                bullets[b].active = false;
                continue;
            }
            // Move bullet
            bullets[b].x += bullets[b].vx * dt;
            bullets[b].y += bullets[b].vy * dt;
            bullets[b].z += bullets[b].vz * dt;
            // Check lifetime or out of bounds (e.g., beyond 50 units or out of arena)
            if (bullets[b].x < -10 || bullets[b].x > 10 || bullets[b].z < -10 || bullets[b].z > 10) {
                bullets[b].active = false;
                continue;
            }
            // Check vertical bounds
            if (bullets[b].y < -BULLET_WORLD_BOUND || bullets[b].y > BULLET_WORLD_BOUND) {
                bullets[b].active = false;
                continue;
            }

            if (!bullets[b].active) continue; // skip further checks if deactivated
            // Check player collisions
            for (int pi = 0; pi < 2; ++pi) {
                if (pi == bullets[b].owner) continue; // skip the shooter
                // Distance from bullet to player
                float dx = bullets[b].x - players[pi].x;
                float dy = bullets[b].y - players[pi].y;
                float dz = bullets[b].z - players[pi].z;
                float distSq = dx*dx + dy*dy + dz*dz;
                if (distSq < (PLAYER_RADIUS * PLAYER_RADIUS)) {
                    // Bullet hit player pi
                    printf("Player %d was hit by Player %d!\n", pi+1, bullets[b].owner+1);


                    // Award a point to the shooter
                    scores[bullets[b].owner] += 1;
                    // For prototype, reset the hit player position (or reduce health)
                    players[pi].x = randomInRange(-9.0f, 9.0f);
                    players[pi].y = BASE_EYE_HEIGHT;
                    players[pi].z = randomInRange(-9.0f, 9.0f);
                    players[pi].vy = 0.0f;
                    players[pi].onGround = true;
                    // Remove bullet
                    bullets[b].active = false;
                }
            }
        }

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
            
            // Draw bullets (small colored cubes)
            glDisable(GL_TEXTURE_2D);
            for (int b = 0; b < MAX_BULLETS; ++b) {
                if (!bullets[b].active) continue;
                // Use white color for bullets
                glColor3f(1.0f, 1.0f, 0.2f);
                // Draw as a small cube (or point)
                DrawColoredCube(bullets[b].x, bullets[b].y, bullets[b].z, 0.1f);
            }

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
            // You can call CreateScoreTexture() to update the texture when the score changes.
            // (For a production game, update only when needed.)
            CreateScoreTexture(scoreStr, &scoreTexture);
            // Render the score texture overlay onto the viewport
            RenderScoreTexture(scoreTexture, /*texWidth*/ 200, /*texHeight*/ 50);
        } // end for each viewport

        SDL_GL_SwapWindow(window);
    } // end main loop

    // Cleanup
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
