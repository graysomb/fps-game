#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
// Forward declarations for procedural floor
void GenerateFloor();
float GetFloorHeight(float x, float z);


// Constants for game mechanics
#define SCREEN_WIDTH  800
#define SCREEN_HEIGHT 600
#define MOVE_SPEED    5.0f    // units per second
#define TURN_SPEED    90.0f   // degrees per second
#define JUMP_SPEED    10.0f    // initial upward velocity for jumps
#define GRAVITY       9.8f    // gravity acceleration (units/s^2)
#define PLAYER_RADIUS 0.5f    // radius for player collision (approximate)
#define BASE_EYE_HEIGHT 1.0f  // eye height when standing on ground
#define WALL_HEIGHT       2.0f  // height of BSP walls and corridors
#define ACCELERATION   40.0f  // horizontal acceleration (units/s^2)
#define FRICTION       20.0f  // ground friction deceleration (units/s^2)

//constants for map
#define NUM_OBSTACLES 2 
// Procedural floor parameters
#define FLOOR_SIZE           20.0f     // extent of floor in X/Z (+/- FLOOR_SIZE)
#define FLOOR_RES            64        // grid resolution per axis
#define FLOOR_NOISE_SCALE    0.03f      // noise frequency scale
#define FLOOR_HEIGHT_SCALE   3.0f      // noise amplitude scale

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

typedef struct {
    float minx, miny, minz;
    float maxx, maxy, maxz;
} Box; // axis-aligned box for obstacles


// Parameters for the BSP generation
#define MIN_ROOM_SIZE 3.0f
#define MAX_BSP_LEVEL 4

// BSP node structure representing a rectangular partition
typedef struct BSPNode {
    float x, z, width, depth; // defines a rectangle on the XZ plane
    struct BSPNode* left;
    struct BSPNode* right;
    Box room;     // room carved inside this partition (if leaf)
    bool isLeaf;
} BSPNode;

//generate wall texture
// A simple integer noise function that returns a value in [-1, 1]
float noise(int x, int y) {
    int n = x + y * 57;
    n = (n << 13) ^ n;
    return ( 1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f );
}

// Interpolated (smoothed) noise
float smoothNoise(float x, float y) {
    int xi = (int)x;
    int yi = (int)y;
    float fracX = x - xi;
    float fracY = y - yi;
    
    float n00 = noise(xi, yi);
    float n10 = noise(xi+1, yi);
    float n01 = noise(xi, yi+1);
    float n11 = noise(xi+1, yi+1);
    
    float i1 = n00 * (1 - fracX) + n10 * fracX;
    float i2 = n01 * (1 - fracX) + n11 * fracX;
    return i1 * (1 - fracY) + i2 * fracY;
}

// Fractal Brownian Motion using multiple octaves of smooth noise
float fbm(float x, float y) {
    float total = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float persistence = 0.5f;
    int octaves = 4;  // Adjust the number of octaves as desired
    for (int i = 0; i < octaves; i++) {
        total += smoothNoise(x * frequency, y * frequency) * amplitude;
        amplitude *= persistence;
        frequency *= 2.0f;
    }
    return total;
}

// Set a random position within a given range
float randomInRange(float min, float max) {
    return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}


// Create a new BSP node
BSPNode* CreateBSPNode(float x, float z, float width, float depth) {
    BSPNode* node = (BSPNode*)malloc(sizeof(BSPNode));
    node->x = x;
    node->z = z;
    node->width = width;
    node->depth = depth;
    node->left = node->right = NULL;
    node->isLeaf = true;
    return node;
}

// Recursively split a node into two subregions
void SplitBSPNode(BSPNode* node, int level) {
    if (level <= 0 || (node->width < MIN_ROOM_SIZE * 2 && node->depth < MIN_ROOM_SIZE * 2))
        return;

    bool splitHorizontally;
    // Choose split orientation based on aspect ratio or randomly
    if (node->width / node->depth >= 1.25f) {
        splitHorizontally = false;
    } else if (node->depth / node->width >= 1.25f) {
        splitHorizontally = true;
    } else {
        splitHorizontally = (rand() % 2) == 0;
    }
    
    if (splitHorizontally && node->depth >= MIN_ROOM_SIZE * 2) {
        // Split along Z (vertical division)
        float split = MIN_ROOM_SIZE + ((float)rand() / RAND_MAX) * (node->depth - 2 * MIN_ROOM_SIZE);
        node->left = CreateBSPNode(node->x, node->z, node->width, split);
        node->right = CreateBSPNode(node->x, node->z + split, node->width, node->depth - split);
    } else if (!splitHorizontally && node->width >= MIN_ROOM_SIZE * 2) {
        // Split along X (horizontal division)
        float split = MIN_ROOM_SIZE + ((float)rand() / RAND_MAX) * (node->width - 2 * MIN_ROOM_SIZE);
        node->left = CreateBSPNode(node->x, node->z, split, node->depth);
        node->right = CreateBSPNode(node->x + split, node->z, node->width - split, node->depth);
    }
    
    if (node->left && node->right) {
        node->isLeaf = false;
        SplitBSPNode(node->left, level - 1);
        SplitBSPNode(node->right, level - 1);
    }
}

// In each leaf node, carve out a room with some margins
void CreateRooms(BSPNode* node) {
    if (node->isLeaf) {
        float margin = 1.0f;
        // Choose a random position within the node (leaving a margin)
        float roomX = node->x + margin + ((float)rand() / RAND_MAX) * (node->width - 2 * margin);
        float roomZ = node->z + margin + ((float)rand() / RAND_MAX) * (node->depth - 2 * margin);
        // Choose room size (ensure a minimum size and keep within the partition)
        float roomW = fmax(2.0f, (node->width - (roomX - node->x)) * 0.8f);
        float roomD = fmax(2.0f, (node->depth - (roomZ - node->z)) * 0.8f);
        node->room.minx = roomX;
        node->room.minz = roomZ;
        node->room.maxx = roomX + roomW;
        node->room.maxz = roomZ + roomD;
        // Set vertical extents based on procedural floor height
        {
            // sample floor height at room center
            float cx = (node->room.minx + node->room.maxx) * 0.5f;
            float cz = (node->room.minz + node->room.maxz) * 0.5f;
            float baseY = GetFloorHeight(cx, cz);
            node->room.miny = baseY;
            node->room.maxy = baseY + WALL_HEIGHT;
        }
    } else {
        if (node->left)
            CreateRooms(node->left);
        if (node->right)
            CreateRooms(node->right);
    }
}

// Structure for corridors connecting rooms
typedef struct {
    float minx, miny, minz;
    float maxx, maxy, maxz;
} Corridor;

// We’ll store corridors in a fixed array (adjust MAX_CORRIDORS as needed)
#define MAX_CORRIDORS 20
Corridor corridors[MAX_CORRIDORS];
int numCorridors = 0;

// Connect two rooms using an L-shaped corridor.
// The corridor is broken into two segments: one horizontal and one vertical.
void ConnectTwoRooms(Box a, Box b, Corridor* corridor1, Corridor* corridor2) {
    // Find centers of the two rooms
    float ax = (a.minx + a.maxx) / 2.0f;
    float az = (a.minz + a.maxz) / 2.0f;
    float bx = (b.minx + b.maxx) / 2.0f;
    float bz = (b.minz + b.maxz) / 2.0f;
    
    // Horizontal segment: from ax to bx at az
    corridor1->minx = fmin(ax, bx);
    corridor1->maxx = fmax(ax, bx);
    corridor1->minz = az - 0.5f; // corridor width of ~1 unit
    corridor1->maxz = az + 0.5f;
    // Vertical extents adapt to floor height at corridor center
    {
        float rx = (corridor1->minx + corridor1->maxx) * 0.5f;
        float rz = (corridor1->minz + corridor1->maxz) * 0.5f;
        float baseY = GetFloorHeight(rx, rz);
        corridor1->miny = baseY;
        corridor1->maxy = baseY + WALL_HEIGHT;
    }
    
    // Vertical segment: from az to bz at bx
    corridor2->minx = bx - 0.5f;
    corridor2->maxx = bx + 0.5f;
    corridor2->minz = fmin(az, bz);
    corridor2->maxz = fmax(az, bz);
    // Vertical extents adapt to floor height at corridor center
    {
        float rx = (corridor2->minx + corridor2->maxx) * 0.5f;
        float rz = (corridor2->minz + corridor2->maxz) * 0.5f;
        float baseY = GetFloorHeight(rx, rz);
        corridor2->miny = baseY;
        corridor2->maxy = baseY + WALL_HEIGHT;
    }
}

// Traverse the BSP tree and, for every internal node, connect the rooms
// of its two children.
void ConnectRooms(BSPNode* node) {
    if (!node || node->isLeaf)
        return;
    
    if (node->left && node->right) {
        // For simplicity, we take the room from each child (in a more robust solution,
        // you might want to get a representative room from each subtree)
        Box roomA = node->left->room;
        Box roomB = node->right->room;
        if (numCorridors + 2 <= MAX_CORRIDORS) {
            ConnectTwoRooms(roomA, roomB, &corridors[numCorridors], &corridors[numCorridors + 1]);
            numCorridors += 2;
        }
    }
    ConnectRooms(node->left);
    ConnectRooms(node->right);
}

// Free the BSP tree to avoid memory leaks
void FreeBSPTree(BSPNode* node) {
    if (!node)
        return;
    FreeBSPTree(node->left);
    FreeBSPTree(node->right);
    free(node);
}

// Call this once at startup to seed the random generator.
void SeedRandom() {
    srand((unsigned)time(NULL));
}

void GenerateRandomObstacles(Box obstacles[]) {
    int count = 0;
    // Define the arena bounds
    const float arenaMin = -9.0f;
    const float arenaMax = 9.0f;
    // For simplicity, divide the arena into a grid
    const int gridCells = 5;
    float cellSize = (arenaMax - arenaMin) / gridCells;
    
    for (int i = 0; i < gridCells && count < NUM_OBSTACLES; ++i) {
        for (int j = 0; j < gridCells && count < NUM_OBSTACLES; ++j) {
            // Use a probability to decide if an obstacle should be placed here
            if (rand() % 100 < 30) { // 30% chance to place an obstacle in this cell
                // Determine a random position within the cell
                float cellMinX = arenaMin + i * cellSize;
                float cellMinZ = arenaMin + j * cellSize;
                float offsetX = ((float)rand() / RAND_MAX) * (cellSize * 0.5f);
                float offsetZ = ((float)rand() / RAND_MAX) * (cellSize * 0.5f);
                
                // Random width and depth (ensuring the box stays within the cell, roughly)
                float width  = cellSize * 0.3f + ((float)rand() / RAND_MAX) * (cellSize * 0.2f);
                float depth  = cellSize * 0.3f + ((float)rand() / RAND_MAX) * (cellSize * 0.2f);
                
                obstacles[count].minx = cellMinX + offsetX;
                obstacles[count].minz = cellMinZ + offsetZ;
                obstacles[count].maxx = obstacles[count].minx + width;
                obstacles[count].maxz = obstacles[count].minz + depth;
                // Fixed vertical dimensions for now
                obstacles[count].miny = 0.0f;
                obstacles[count].maxy = 2.0f;
                
                count++;
            }
        }
    }
}


BSPNode* root = NULL;
// Global game state
Player players[2];
#define MAX_BULLETS 50
#define BULLET_LIFETIME 3.0f        // bullet time-to-live in seconds
#define BULLET_WORLD_BOUND 20.0f    // world bounds for bullet Y coordinate
Bullet bullets[MAX_BULLETS];
int numBullets = 0;  // or we can reuse bullets in a pool
// Define a few obstacles in the scene (e.g., 2 boxes)
Box obstacles[NUM_OBSTACLES];

// Texture ID
GLuint texChecker = 0;

// Initialize a checkerboard texture (64x64, 2-color)
void CreateCheckerTexture() {
    const int TEX_SIZE = 64;
    unsigned char image[TEX_SIZE * TEX_SIZE * 3];
    for(int j=0; j<TEX_SIZE; ++j) {
        for(int i=0; i<TEX_SIZE; ++i) {
            int color = (((i / 8) + (j / 8)) % 2) ? 255 : 50; // alternate colors
            // White (255) and dark gray (50) pattern
            image[(j*TEX_SIZE + i)*3 + 0] = color;
            image[(j*TEX_SIZE + i)*3 + 1] = color;
            image[(j*TEX_SIZE + i)*3 + 2] = color;
        }
    }
    glGenTextures(1, &texChecker);
    glBindTexture(GL_TEXTURE_2D, texChecker);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // Upload texture data (RGB format)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TEX_SIZE, TEX_SIZE, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, image);
}

// --- Create a Rocky Texture ---
// Texture ID
GLuint texRock = 1;
// Generate a procedural rock texture (128x128)
void CreateRockTexture() {
    const int TEX_SIZE = 128;
    unsigned char image[TEX_SIZE * TEX_SIZE * 3];
    float scale = 8.0f;  // Controls the "zoom" level of the noise
    
    // With 4 octaves, the maximum possible amplitude sum is ~1 + 0.5 + 0.25 + 0.125 = 1.875.
    // We'll normalize fbm() output (which is roughly in [-1.875, 1.875]) to [0,1].
    float amplitudeSum = 1.875f;
    
    // Loop over each pixel
    for (int j = 0; j < TEX_SIZE; ++j) {
        for (int i = 0; i < TEX_SIZE; ++i) {
            // Map pixel coordinates to noise space
            float x = ((float)i / TEX_SIZE) * scale;
            float y = ((float)j / TEX_SIZE) * scale;
            float value = fbm(x, y);
            // Normalize from roughly [-amplitudeSum, amplitudeSum] to [0,1]
            value = (value + amplitudeSum) / (2.0f * amplitudeSum);
            
            // Map the noise value to a rocky color.
            // Here we use a base rock color with slight variation.
            // You can experiment with these formulas to get different hues.
            int r = (int)(value * 80 + 70 + (rand() % 10));  // A dark, earthy red tone
            int g = (int)(value * 80 + 70 + (rand() % 10));  // Similar green value for grayish look
            int b = (int)(value * 80 + 90 + (rand() % 10));  // A touch more blue for a cool rock feel

            // Clamp to [0, 255]
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            
            image[(j * TEX_SIZE + i) * 3 + 0] = (unsigned char)r;
            image[(j * TEX_SIZE + i) * 3 + 1] = (unsigned char)g;
            image[(j * TEX_SIZE + i) * 3 + 2] = (unsigned char)b;
        }
    }
    
    // Create the OpenGL texture object
    glGenTextures(1, &texRock);
    glBindTexture(GL_TEXTURE_2D, texRock);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // Upload texture data (RGB format)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TEX_SIZE, TEX_SIZE, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, image);
    
    // You can store texRock in a global variable for later use in your render loop.
}


// --- Create a Grassy Texture ---
// Texture ID
GLuint texGrass = 2;
// Generate a procedural rock texture (128x128)
void CreateGrassTexture() {
    const int TEX_SIZE = 128;
    unsigned char image[TEX_SIZE * TEX_SIZE * 3];
    float scale = 8.0f;  // Controls the "zoom" level of the noise
    
    // With 4 octaves, the maximum possible amplitude sum is ~1 + 0.5 + 0.25 + 0.125 = 1.875.
    // We'll normalize fbm() output (which is roughly in [-1.875, 1.875]) to [0,1].
    float amplitudeSum = 1.875f;
    
    // Loop over each pixel
    for (int j = 0; j < TEX_SIZE; ++j) {
        for (int i = 0; i < TEX_SIZE; ++i) {
            // Map pixel coordinates to noise space
            float x = ((float)i / TEX_SIZE) * scale;
            float y = ((float)j / TEX_SIZE) * scale;
            float value = fbm(x, y);
            // Normalize from roughly [-amplitudeSum, amplitudeSum] to [0,1]
            value = (value + amplitudeSum) / (2.0f * amplitudeSum);
            
            // Map the noise value to a rocky color.
            // Here we use a base rock color with slight variation.
            // You can experiment with these formulas to get different hues.
            int r = (int)(value * 80 + 30 + (rand() % 10));  // A dark, earthy red tone
            int g = (int)(value * 80 + 70 + (rand() % 10));  // Similar green value for grayish look
            int b = (int)(value * 80 + 30 + (rand() % 10));  // A touch more blue for a cool rock feel

            // Clamp to [0, 255]
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            
            image[(j * TEX_SIZE + i) * 3 + 0] = (unsigned char)r;
            image[(j * TEX_SIZE + i) * 3 + 1] = (unsigned char)g;
            image[(j * TEX_SIZE + i) * 3 + 2] = (unsigned char)b;
        }
    }
    
    // Create the OpenGL texture object
    glGenTextures(1, &texGrass);
    glBindTexture(GL_TEXTURE_2D, texGrass);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // Upload texture data (RGB format)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TEX_SIZE, TEX_SIZE, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, image);
    
    // You can store texRock in a global variable for later use in your render loop.
}


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

// Render helper: draw a textured cube given a Box (assumes texture is bound and enabled)
void DrawTexturedBox(Box b) {
    // Each face as two triangles or a quad
    // Front face (z = maxz)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.minx, b.miny, b.maxz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.maxx, b.miny, b.maxz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.maxx, b.maxy, b.maxz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.minx, b.maxy, b.maxz);
    glEnd();
    // Back face (z = minz)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.maxx, b.miny, b.minz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.minx, b.miny, b.minz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.minx, b.maxy, b.minz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.maxx, b.maxy, b.minz);
    glEnd();
    // Left face (x = minx)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.minx, b.miny, b.minz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.minx, b.miny, b.maxz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.minx, b.maxy, b.maxz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.minx, b.maxy, b.minz);
    glEnd();
    // Right face (x = maxx)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.maxx, b.miny, b.maxz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.maxx, b.miny, b.minz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.maxx, b.maxy, b.minz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.maxx, b.maxy, b.maxz);
    glEnd();
    // Bottom face (y = miny)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.minx, b.miny, b.minz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.maxx, b.miny, b.minz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.maxx, b.miny, b.maxz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.minx, b.miny, b.maxz);
    glEnd();
    // Top face (y = maxy)
    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f); glVertex3f(b.minx, b.maxy, b.maxz);
      glTexCoord2f(1.0f, 0.0f); glVertex3f(b.maxx, b.maxy, b.maxz);
      glTexCoord2f(1.0f, 1.0f); glVertex3f(b.maxx, b.maxy, b.minz);
      glTexCoord2f(0.0f, 1.0f); glVertex3f(b.minx, b.maxy, b.minz);
    glEnd();
}

// Draw each room in the BSP tree (leaf nodes)
void DrawBSPRooms(BSPNode *node) {
    if (!node)
        return;
    if (node->isLeaf) {
        // Use a distinctive color (for example, gray) for rooms
        //glColor3f(0.7f, 0.7f, 0.7f);
        DrawTexturedBox(node->room);
    } else {
        DrawBSPRooms(node->left);
        DrawBSPRooms(node->right);
    }
}

// Draw all corridors stored in the corridors array
void DrawCorridors() {
    for (int i = 0; i < numCorridors; i++) {
        // Use a different color (for example, light gray) for corridors
        //glColor3f(0.9f, 0.9f, 0.9f);
        // Create a temporary Box from the Corridor fields
        Box corridorBox = {
            corridors[i].minx, corridors[i].miny, corridors[i].minz,
            corridors[i].maxx, corridors[i].maxy, corridors[i].maxz
        };
        DrawTexturedBox(corridorBox);
    }
}

#define WALL_THICKNESS 0.2f
// Height of BSP walls and corridors
#define WALL_HEIGHT    2.0f

// Add wall obstacles for a given room (leaf node)
void GenerateRoomWalls(Box room, Box obstacles[], int *count, int maxObstacles) {
    if (*count + 4 > maxObstacles) return; // Ensure we don't overflow

    // Left wall (along minx)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.minx + WALL_THICKNESS;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.maxz;
    // set vertical extents based on floor height at wall center
    {
        float wx = (obstacles[*count].minx + obstacles[*count].maxx) * 0.5f;
        float wz = (obstacles[*count].minz + obstacles[*count].maxz) * 0.5f;
        float baseY = GetFloorHeight(wx, wz);
        obstacles[*count].miny = baseY;
        obstacles[*count].maxy = baseY + WALL_HEIGHT;
    }
    (*count)++;

    // Right wall (along maxx)
    obstacles[*count].minx = room.maxx - WALL_THICKNESS;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.maxz;
    // set vertical extents based on floor height at wall center
    {
        float wx = (obstacles[*count].minx + obstacles[*count].maxx) * 0.5f;
        float wz = (obstacles[*count].minz + obstacles[*count].maxz) * 0.5f;
        float baseY = GetFloorHeight(wx, wz);
        obstacles[*count].miny = baseY;
        obstacles[*count].maxy = baseY + WALL_HEIGHT;
    }
    (*count)++;

    // Top wall (along maxz)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.maxz - WALL_THICKNESS;
    obstacles[*count].maxz = room.maxz;
    // set vertical extents based on floor height at wall center
    {
        float wx = (obstacles[*count].minx + obstacles[*count].maxx) * 0.5f;
        float wz = (obstacles[*count].minz + obstacles[*count].maxz) * 0.5f;
        float baseY = GetFloorHeight(wx, wz);
        obstacles[*count].miny = baseY;
        obstacles[*count].maxy = baseY + WALL_HEIGHT;
    }
    (*count)++;

    // Bottom wall (along minz)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.minz + WALL_THICKNESS;
    // set vertical extents based on floor height at wall center
    {
        float wx = (obstacles[*count].minx + obstacles[*count].maxx) * 0.5f;
        float wz = (obstacles[*count].minz + obstacles[*count].maxz) * 0.5f;
        float baseY = GetFloorHeight(wx, wz);
        obstacles[*count].miny = baseY;
        obstacles[*count].maxy = baseY + WALL_HEIGHT;
    }
    (*count)++;
}

#define MAX_COLLISION_OBSTACLES 100
Box collisionObstacles[MAX_COLLISION_OBSTACLES];
int numCollisionObstacles = 0;
// Health for each wall segment in collisionObstacles
#define WALL_MAX_HEALTH 3  // number of hits to destroy a wall segment
int wallHealth[MAX_COLLISION_OBSTACLES];

// Traverse the BSP tree and generate walls for each room (leaf node)
void CollectRoomWalls(BSPNode *node) {
    if (!node)
        return;
    if (node->isLeaf) {
        GenerateRoomWalls(node->room, collisionObstacles, &numCollisionObstacles, MAX_COLLISION_OBSTACLES);
    } else {
        CollectRoomWalls(node->left);
        CollectRoomWalls(node->right);
    }
}

// If you have corridors that should be collidable, you can add them similarly:
void CollectCorridorWalls() {
    for (int i = 0; i < numCorridors; i++) {
        // For corridors, you might treat the corridor itself as an obstacle.
        // Here we add the corridor box directly (or generate thin walls along its edges if desired).
        if (numCollisionObstacles < MAX_COLLISION_OBSTACLES) {
            collisionObstacles[numCollisionObstacles].minx = corridors[i].minx;
            collisionObstacles[numCollisionObstacles].maxx = corridors[i].maxx;
            collisionObstacles[numCollisionObstacles].minz = corridors[i].minz;
            collisionObstacles[numCollisionObstacles].maxz = corridors[i].maxz;
            collisionObstacles[numCollisionObstacles].miny = corridors[i].miny;
            collisionObstacles[numCollisionObstacles].maxy = corridors[i].maxy;
            numCollisionObstacles++;
        }
    }
}

// Floor height map array (size FLOOR_RES+1 by FLOOR_RES+1)
static float floorHeights[FLOOR_RES+1][FLOOR_RES+1];

// Generates a fractal noise-based heightmap for the floor
void GenerateFloor() {
    for (int iz = 0; iz <= FLOOR_RES; ++iz) {
        for (int ix = 0; ix <= FLOOR_RES; ++ix) {
            float x = -FLOOR_SIZE + (2.0f * FLOOR_SIZE) * ix / (float)FLOOR_RES;
            float z = -FLOOR_SIZE + (2.0f * FLOOR_SIZE) * iz / (float)FLOOR_RES;
            // Fractal Brownian Motion noise
            float h = fbm(x * FLOOR_NOISE_SCALE, z * FLOOR_NOISE_SCALE) * FLOOR_HEIGHT_SCALE;
            floorHeights[ix][iz] = h;
        }
    }
}

// Bilinearly sample the floor heightmap at world coordinates (x,z)
float GetFloorHeight(float x, float z) {
    float fx = (x + FLOOR_SIZE) * (FLOOR_RES / (2.0f * FLOOR_SIZE));
    float fz = (z + FLOOR_SIZE) * (FLOOR_RES / (2.0f * FLOOR_SIZE));
    int ix0 = (int)floorf(fx);
    int iz0 = (int)floorf(fz);
    float tx = fx - ix0;
    float tz = fz - iz0;
    if (ix0 < 0) ix0 = 0;
    if (iz0 < 0) iz0 = 0;
    if (ix0 >= FLOOR_RES) ix0 = FLOOR_RES - 1;
    if (iz0 >= FLOOR_RES) iz0 = FLOOR_RES - 1;
    int ix1 = ix0 + 1;
    int iz1 = iz0 + 1;
    float h00 = floorHeights[ix0][iz0];
    float h10 = floorHeights[ix1][iz0];
    float h01 = floorHeights[ix0][iz1];
    float h11 = floorHeights[ix1][iz1];
    float h0 = h00 * (1.0f - tx) + h10 * tx;
    float h1 = h01 * (1.0f - tx) + h11 * tx;
    return h0 * (1.0f - tz) + h1 * tz;
}
// Deform floor heightmap at given world position (hitX, hitZ) with crater effect
// Deform floor heightmap around (hitX, hitZ) within radius
void Deform(float hitX, float hitZ, float radius, float depth) {
    // Compute grid conversion factor (cells per unit)
    float invGrid = (FLOOR_RES / (2.0f * FLOOR_SIZE));
    // Convert world coords to grid indices
    float fx = (hitX + FLOOR_SIZE) * invGrid;
    float fz = (hitZ + FLOOR_SIZE) * invGrid;
    int ixCenter = (int)fx;
    int izCenter = (int)fz;
    // Compute affected radius in grid cells
    int radiusCells = (int)(radius * invGrid) + 1;
    // Clamp index ranges
    int ix0 = ixCenter - radiusCells; if (ix0 < 0) ix0 = 0;
    int iz0 = izCenter - radiusCells; if (iz0 < 0) iz0 = 0;
    int ix1 = ixCenter + radiusCells; if (ix1 > FLOOR_RES) ix1 = FLOOR_RES;
    int iz1 = izCenter + radiusCells; if (iz1 > FLOOR_RES) iz1 = FLOOR_RES;
    float radiusSq = radius * radius;
    // Only iterate over affected cells
    for (int iz = iz0; iz <= iz1; ++iz) {
        for (int ix = ix0; ix <= ix1; ++ix) {
            // Compute world position of cell center
            float x = -FLOOR_SIZE + (2.0f * FLOOR_SIZE) * ix / (float)FLOOR_RES;
            float z = -FLOOR_SIZE + (2.0f * FLOOR_SIZE) * iz / (float)FLOOR_RES;
            float dx = x - hitX;
            float dz = z - hitZ;
            float distSq = dx*dx + dz*dz;
            if (distSq <= radiusSq) {
                float falloff = 1.0f - (sqrtf(distSq) / radius);
                floorHeights[ix][iz] -= depth * falloff;
            }
        }
    }
}

// Draws the floor as a triangulated height mesh using the grass texture
void DrawProceduralFloor() {
    glBindTexture(GL_TEXTURE_2D, texGrass);
    glEnable(GL_TEXTURE_2D);
    glColor3f(1.0f, 1.0f, 1.0f);
    float step = (2.0f * FLOOR_SIZE) / FLOOR_RES;
    float texStep = 10.0f / FLOOR_RES; // reuse tileRepeat = 10
    for (int iz = 0; iz < FLOOR_RES; ++iz) {
        float z0 = -FLOOR_SIZE + step * iz;
        float z1 = z0 + step;
        float v0 = texStep * iz;
        float v1 = v0 + texStep;
        glBegin(GL_TRIANGLES);
        for (int ix = 0; ix < FLOOR_RES; ++ix) {
            float x0 = -FLOOR_SIZE + step * ix;
            float x1 = x0 + step;
            float u0 = texStep * ix;
            float u1 = u0 + texStep;
            float h00 = floorHeights[ix][iz];
            float h10 = floorHeights[ix+1][iz];
            float h11 = floorHeights[ix+1][iz+1];
            float h01 = floorHeights[ix][iz+1];
            // Triangle 1
            glTexCoord2f(u0, v0); glVertex3f(x0, h00, z0);
            glTexCoord2f(u1, v0); glVertex3f(x1, h10, z0);
            glTexCoord2f(u1, v1); glVertex3f(x1, h11, z1);
            // Triangle 2
            glTexCoord2f(u0, v0); glVertex3f(x0, h00, z0);
            glTexCoord2f(u1, v1); glVertex3f(x1, h11, z1);
            glTexCoord2f(u0, v1); glVertex3f(x0, h01, z1);
        }
        glEnd();
    }
    glDisable(GL_TEXTURE_2D);
}


// Reset game (players positions, etc.)
void ResetGame() {
    // Generate procedural floor heightmap before BSP and walls
    GenerateFloor();
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
    //GenerateRandomObstacles(obstacles);
    // gen map

    // Define your overall map area based on procedural floor size
    float mapX = -FLOOR_SIZE;
    float mapZ = -FLOOR_SIZE;
    float mapWidth = FLOOR_SIZE * 2.0f;
    float mapDepth = FLOOR_SIZE * 2.0f;
    root = CreateBSPNode(mapX, mapZ, mapWidth, mapDepth);

    // Recursively partition the map
    SplitBSPNode(root, MAX_BSP_LEVEL);

    // In each leaf, carve out a room
    CreateRooms(root);

    // Connect rooms between partitions with corridors
    numCorridors = 0;
    ConnectRooms(root);
    // Reset collision obstacles
    numCollisionObstacles = 0;
    // Generate walls from the rooms in the BSP tree
    CollectRoomWalls(root);
    // Optionally, add corridor obstacles as well
    CollectCorridorWalls();
    // Initialize health for each wall segment
    for (int w = 0; w < numCollisionObstacles; ++w) {
        wallHealth[w] = WALL_MAX_HEALTH;
    }

/*     // Define obstacle boxes (e.g., two blocks in the arena)
    obstacles[0].minx = -1.0f; obstacles[0].maxx = 1.0f;
    obstacles[0].miny = 0.0f;  obstacles[0].maxy = 2.0f;
    obstacles[0].minz = -4.0f; obstacles[0].maxz = -2.0f;
    // A wall or block in front of player1's side
    obstacles[1].minx = -0.5f; obstacles[1].maxx = 0.5f;
    obstacles[1].miny = 0.0f;  obstacles[1].maxy = 2.0f;
    obstacles[1].minz = 2.0f;  obstacles[1].maxz = 4.0f;
    // (You can add more obstacles or adjust positions/sizes as needed) */
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

    // Create textures (checkerboard)
    CreateCheckerTexture();
    CreateRockTexture();
    CreateGrassTexture();

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
            // Check floor collision and deform floor
            {
                float floorY = GetFloorHeight(bullets[b].x, bullets[b].z);
                if (bullets[b].y <= floorY) {
                    Deform(bullets[b].x, bullets[b].z, 0.6f, 0.3f);
                    bullets[b].active = false;
                    continue;
                }
            }
            // Check wall collisions and apply damage
            for (int o = 0; o < numCollisionObstacles; ++o) {
                if (wallHealth[o] <= 0) continue;
                Box *wall = &collisionObstacles[o];
                if (bullets[b].x >= wall->minx && bullets[b].x <= wall->maxx &&
                    bullets[b].y >= wall->miny && bullets[b].y <= wall->maxy &&
                    bullets[b].z >= wall->minz && bullets[b].z <= wall->maxz) {
                    // Bullet hit a wall: decrement health
                    wallHealth[o]--;
                    bullets[b].active = false;
                    break;
                }
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

            // Draw procedural floor mesh
            DrawProceduralFloor();
            // Draw wall segments from BSP (collision obstacles)
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, texRock);
            for (int w = 0; w < numCollisionObstacles; ++w) {
                if (wallHealth[w] <= 0) continue;
                DrawTexturedBox(collisionObstacles[w]);
            }
            glDisable(GL_TEXTURE_2D);

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
