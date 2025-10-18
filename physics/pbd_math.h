#ifndef PBD_MATH_H
#define PBD_MATH_H

#include <math.h>
#include "raylib.h"

static inline float pbd_clampf(float v, float lo, float hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static inline float pbd_mixf(float a, float b, float t) {
    return a * (1.0f - t) + b * t;
}

static inline Vector3 pbd_v3_add(Vector3 a, Vector3 b) {
    return (Vector3){ a.x + b.x, a.y + b.y, a.z + b.z };
}

static inline Vector3 pbd_v3_sub(Vector3 a, Vector3 b) {
    return (Vector3){ a.x - b.x, a.y - b.y, a.z - b.z };
}

static inline Vector3 pbd_v3_scale(Vector3 v, float s) {
    return (Vector3){ v.x * s, v.y * s, v.z * s };
}

static inline float pbd_v3_dot(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vector3 pbd_v3_cross(Vector3 a, Vector3 b) {
    return (Vector3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static inline float pbd_v3_length(Vector3 v) {
    return sqrtf(pbd_v3_dot(v, v));
}

static inline Vector3 pbd_v3_normalize(Vector3 v) {
    float len = pbd_v3_length(v);
    if (len <= 0.0f) {
        return (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    return pbd_v3_scale(v, 1.0f / len);
}

static inline Vector3 pbd_v3_project(Vector3 onto, Vector3 vec, float eps) {
    float denom = pbd_v3_dot(onto, onto);
    if (denom < eps) {
        return (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    float scale = pbd_v3_dot(onto, vec) / denom;
    return pbd_v3_scale(onto, scale);
}

static inline Vector3 pbd_v3_limit(Vector3 v, float max_len) {
    float len = pbd_v3_length(v);
    if (len > max_len && max_len > 0.0f) {
        return pbd_v3_scale(v, max_len / len);
    }
    return v;
}

#endif // PBD_MATH_H
