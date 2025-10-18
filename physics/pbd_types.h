#ifndef PBD_TYPES_H
#define PBD_TYPES_H

#include <stdint.h>
#include "raylib.h"

typedef struct {
    Vector3 pos;
    Vector3 prev_pos;
    Vector3 predicted_pos;
    Vector3 vel;
    float inv_mass;
    float radius;
    float friction;
    uint32_t flags;
} PbdParticle;

typedef struct {
    int indices[8];
    int is_boundary;
} PbdVoxelConstraint;

typedef struct {
    int voxelA;
    int voxelB;
    int faceIndex; // axis-aligned face id (0..5)
    float strain_min;
    float strain_max;
    int active;
} PbdFaceConstraint;

typedef struct {
    int numParticles;
    int numCells;
    Vector3 gridMin;
    Vector3 gridMax;
    int gridRes[3];
    Vector3 cellSize;
    int maxCellsPerElement;
    Vector4 wallConstraints[6];
    Vector4 gravity;
    float dt;
    float compliance;
    float omega_collision;
    float time;
} PbdSystemData;

#endif // PBD_TYPES_H
