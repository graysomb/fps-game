#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


// Constants for game mechanics
#define SCREEN_WIDTH  800
#define SCREEN_HEIGHT 600
#define MOVE_SPEED    5.0f    // units per second
#define TURN_SPEED    90.0f   // degrees per second
#define JUMP_SPEED    10.0f    // initial upward velocity for jumps
#define GRAVITY       9.8f    // gravity acceleration (units/s^2)
#define PLAYER_RADIUS 0.5f    // radius for player collision (approximate)
#define BASE_EYE_HEIGHT 1.0f  // eye height when standing on ground

//constants for map
#define NUM_OBSTACLES 2 

// Key mappings (SDL Scancodes for player actions)
const SDL_Scancode P1_FORWARD    = SDL_SCANCODE_W;
const SDL_Scancode P1_BACKWARD  = SDL_SCANCODE_S;
const SDL_Scancode P1_STRAFE_L  = SDL_SCANCODE_A;
const SDL_Scancode P1_STRAFE_R  = SDL_SCANCODE_D;
const SDL_Scancode P1_TURN_L    = SDL_SCANCODE_L;
const SDL_Scancode P1_TURN_R    = SDL_SCANCODE_J;
const SDL_Scancode P1_JUMP      = SDL_SCANCODE_SPACE;
const SDL_Scancode P1_SHOOT     = SDL_SCANCODE_LCTRL;

const SDL_Scancode P2_FORWARD   = SDL_SCANCODE_UP;
const SDL_Scancode P2_BACKWARD = SDL_SCANCODE_DOWN;
const SDL_Scancode P2_TURN_L   = SDL_SCANCODE_LEFT;
const SDL_Scancode P2_TURN_R   = SDL_SCANCODE_RIGHT;
const SDL_Scancode P2_STRAFE_L = SDL_SCANCODE_COMMA;   // ',' key
const SDL_Scancode P2_STRAFE_R = SDL_SCANCODE_PERIOD;  // '.' key
const SDL_Scancode P2_JUMP     = SDL_SCANCODE_RSHIFT;
const SDL_Scancode P2_SHOOT    = SDL_SCANCODE_RCTRL;

// Structure definitions
typedef struct {
    float x, y, z;
    float yaw;        // rotation around Y-axis in degrees
    float vy;         // vertical velocity
    bool  onGround;
} Player;

typedef struct {
    float x, y, z;
    float vx, vy, vz;
    bool active;
    int owner;  // which player fired (0 or 1)
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
        // Set vertical extents (you can adjust these as needed)
        node->room.miny = 0.0f;
        node->room.maxy = 2.0f;
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
    corridor1->miny = 0.0f;
    corridor1->maxy = 2.0f;
    
    // Vertical segment: from az to bz at bx
    corridor2->minx = bx - 0.5f;
    corridor2->maxx = bx + 0.5f;
    corridor2->minz = fmin(az, bz);
    corridor2->maxz = fmax(az, bz);
    corridor2->miny = 0.0f;
    corridor2->maxy = 2.0f;
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

// Add wall obstacles for a given room (leaf node)
void GenerateRoomWalls(Box room, Box obstacles[], int *count, int maxObstacles) {
    if (*count + 4 > maxObstacles) return; // Ensure we don't overflow

    // Left wall (along minx)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.minx + WALL_THICKNESS;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.maxz;
    obstacles[*count].miny = 0.0f;
    obstacles[*count].maxy = 2.0f;
    (*count)++;

    // Right wall (along maxx)
    obstacles[*count].minx = room.maxx - WALL_THICKNESS;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.maxz;
    obstacles[*count].miny = 0.0f;
    obstacles[*count].maxy = 2.0f;
    (*count)++;

    // Top wall (along maxz)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.maxz - WALL_THICKNESS;
    obstacles[*count].maxz = room.maxz;
    obstacles[*count].miny = 0.0f;
    obstacles[*count].maxy = 2.0f;
    (*count)++;

    // Bottom wall (along minz)
    obstacles[*count].minx = room.minx;
    obstacles[*count].maxx = room.maxx;
    obstacles[*count].minz = room.minz;
    obstacles[*count].maxz = room.minz + WALL_THICKNESS;
    obstacles[*count].miny = 0.0f;
    obstacles[*count].maxy = 2.0f;
    (*count)++;
}

#define MAX_COLLISION_OBSTACLES 100
Box collisionObstacles[MAX_COLLISION_OBSTACLES];
int numCollisionObstacles = 0;

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


// Reset game (players positions, etc.)
void ResetGame() {
    // Initialize players at different positions
    players[0].x = -3.0f; players[0].y = BASE_EYE_HEIGHT; players[0].z = 0.0f;
    players[0].yaw = 0.0f; players[0].vy = 0.0f; players[0].onGround = true;
    players[1].x = 3.0f;  players[1].y = BASE_EYE_HEIGHT; players[1].z = 0.0f;
    players[1].yaw = 180.0f; players[1].vy = 0.0f; players[1].onGround = true;
    // Clear bullets
    numBullets = 0;
    for(int i=0; i<MAX_BULLETS; ++i) bullets[i].active = false;
    //GenerateRandomObstacles(obstacles);
    // gen map

    // Define your overall map area (for example, an arena from -9 to +9 in X and Z)
    float mapX = -9.0f, mapZ = -9.0f;
    float mapWidth = 18.0f, mapDepth = 18.0f;
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
            float yawRad = players[playerIndex].yaw * M_PI / 180.0f;
            float dirx = sinf(-yawRad);
            float dirz = -cosf(yawRad);
            bullets[i].x = players[playerIndex].x + dirx * 0.8f; // offset forward a bit
            bullets[i].y = players[playerIndex].y; 
            bullets[i].z = players[playerIndex].z + dirz * 0.8f;
            bullets[i].vx = dirx * 20.0f; // bullet speed ~20 units/s
            bullets[i].vy = 0.0f;
            bullets[i].vz = dirz * 20.0f;
            bullets[i].active = true;
            bullets[i].owner = playerIndex;
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
        
        // Check against world bounds first (optional)
        if (p->x < -9.0f + PLAYER_RADIUS) { p->x = -9.0f + PLAYER_RADIUS; collisionFound = true; }
        else if (p->x > 9.0f - PLAYER_RADIUS) { p->x = 9.0f - PLAYER_RADIUS; collisionFound = true; }
        if (p->z < -9.0f + PLAYER_RADIUS) { p->z = -9.0f + PLAYER_RADIUS; collisionFound = true; }
        else if (p->z > 9.0f - PLAYER_RADIUS) { p->z = 9.0f - PLAYER_RADIUS; collisionFound = true; }
        
        // Check against all collision obstacles
        for (int o = 0; o < numCollisionObstacles; ++o) {
            if (CollidePlayerWithBox(p->x, p->z, PLAYER_RADIUS, collisionObstacles[o])) {
                ResolvePlayerCollision(p, collisionObstacles[o]);
                collisionFound = true;
            }
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

    // Create textures (checkerboard)
    CreateCheckerTexture();
    CreateRockTexture();

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
            // Movement forward/back/strafe
            float yawRad = players[i].yaw * (M_PI / 180.0f);
            float forwardX = sinf(-yawRad);
            float forwardZ = -cosf(yawRad);
            float rightX   = cosf(yawRad);
            float rightZ   = sinf(yawRad);
            float moveDX = 0.0f;
            float moveDZ = 0.0f;
            if ((i == 0 && keystate[P1_FORWARD]) || (i == 1 && keystate[P2_FORWARD])) {
                moveDX += forwardX * MOVE_SPEED * dt;
                moveDZ += forwardZ * MOVE_SPEED * dt;
            }
            if ((i == 0 && keystate[P1_BACKWARD]) || (i == 1 && keystate[P2_BACKWARD])) {
                moveDX -= forwardX * MOVE_SPEED * dt;
                moveDZ -= forwardZ * MOVE_SPEED * dt;
            }
            if ((i == 0 && keystate[P1_STRAFE_R]) || (i == 1 && keystate[P2_STRAFE_R])) {
                moveDX += rightX * MOVE_SPEED * dt;
                moveDZ += rightZ * MOVE_SPEED * dt;
            }
            if ((i == 0 && keystate[P1_STRAFE_L]) || (i == 1 && keystate[P2_STRAFE_L])) {
                moveDX -= rightX * MOVE_SPEED * dt;
                moveDZ -= rightZ * MOVE_SPEED * dt;
            }

            // Apply movement with collision check:
            if (moveDX != 0.0f || moveDZ != 0.0f) {
                /* // Check against obstacles by attempting movement in X and Z separately (simplified collision response)
                float newX = players[i].x + moveDX;
                float newZ = players[i].z + moveDZ;
                bool collisionX = false;
                bool collisionZ = false;
                // World bounds check (keep players within +/-9 in X,Z)
                if (newX < -9.0f || newX > 9.0f) collisionX = true;
                if (newZ < -9.0f || newZ > 9.0f) collisionZ = true;
                // Obstacle collision
                for (int o = 0; o < numCollisionObstacles; ++o) {
                    // Check intended X movement
                    if (!collisionX) {
                        if (CollidePlayerWithBox(newX, players[i].z, PLAYER_RADIUS, collisionObstacles[o])) {
                            collisionX = true;
                        }
                    }
                    // Check intended Z movement
                    if (!collisionZ) {
                        if (CollidePlayerWithBox(players[i].x, newZ, PLAYER_RADIUS, collisionObstacles[o])) {
                            collisionZ = true;
                        }
                    }
                }
                // Apply movement if no collision
                if (!collisionX) players[i].x = newX;
                if (!collisionZ) players[i].z = newZ; */

                // Compute intended new position based on input
                float newX = players[i].x + moveDX;
                float newZ = players[i].z + moveDZ;

                // Update position (even if this results in overlap)
                players[i].x = newX;
                players[i].z = newZ;

                // Now resolve any collisions (iteratively push out of obstacles)
                ResolvePlayerCollisions(&players[i]);
            }

            // Apply gravity and jumping physics for player i
            if (!players[i].onGround) {
                players[i].vy -= GRAVITY * dt;
            }
            players[i].y += players[i].vy * dt;
            // Ground collision
            if (players[i].y <= BASE_EYE_HEIGHT) {
                players[i].y = BASE_EYE_HEIGHT;
                players[i].vy = 0.0f;
                players[i].onGround = true;
            }
        } // end player loop

        // Update bullets
        for (int b = 0; b < MAX_BULLETS; ++b) {
            if (!bullets[b].active) continue;
            // Move bullet
            bullets[b].x += bullets[b].vx * dt;
            bullets[b].y += bullets[b].vy * dt;
            bullets[b].z += bullets[b].vz * dt;
            // Check lifetime or out of bounds (e.g., beyond 50 units or out of arena)
            if (bullets[b].x < -10 || bullets[b].x > 10 || bullets[b].z < -10 || bullets[b].z > 10) {
                bullets[b].active = false;
                continue;
            }
            // Check obstacle collisions
            for (int o = 0; o < NUM_OBSTACLES; ++o) {
                if (bullets[b].x >= obstacles[o].minx && bullets[b].x <= obstacles[o].maxx &&
                    bullets[b].y >= obstacles[o].miny && bullets[b].y <= obstacles[o].maxy &&
                    bullets[b].z >= obstacles[o].minz && bullets[b].z <= obstacles[o].maxz) {
                    // Bullet hit an obstacle
                    bullets[b].active = false;
                    continue;
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
                    // For prototype, reset the hit player position (or reduce health)
                    players[pi].x = (pi == 0 ? -3.0f : 3.0f);
                    players[pi].y = BASE_EYE_HEIGHT;
                    players[pi].z = 0.0f;
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
            // Camera rotations (note: no pitch implemented, only yaw)
            glRotatef(-players[i].yaw, 0.0f, 1.0f, 0.0f);
            glTranslatef(-players[i].x, -players[i].y, -players[i].z);

            // Draw floor (textured)
            glBindTexture(GL_TEXTURE_2D, texChecker);
            glEnable(GL_TEXTURE_2D);
            glColor3f(1,1,1);
            float floorSize = 10.0f;
            float tileRepeat = 10.0f;
            glBegin(GL_QUADS);
              glTexCoord2f(0.0f, 0.0f); glVertex3f(-floorSize, 0.0f, -floorSize);
              glTexCoord2f(tileRepeat, 0.0f); glVertex3f( floorSize, 0.0f, -floorSize);
              glTexCoord2f(tileRepeat, tileRepeat); glVertex3f( floorSize, 0.0f,  floorSize);
              glTexCoord2f(0.0f, tileRepeat); glVertex3f(-floorSize, 0.0f,  floorSize);
            glEnd();

            // Draw obstacles (textured cubes)
            glBindTexture(GL_TEXTURE_2D, texRock);
            for (int o = 0; o < NUM_OBSTACLES; ++o) {
                DrawTexturedBox(obstacles[o]);
            }


            // Disable texture if you want plain colored BSP elements
            //glDisable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, texRock);

            // Draw BSP generated rooms
            DrawBSPRooms(root);

            // Draw corridors connecting the rooms
            DrawCorridors();

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
        } // end for each viewport

        SDL_GL_SwapWindow(window);
    } // end main loop

    // Cleanup
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
