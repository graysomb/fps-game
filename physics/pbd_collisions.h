#ifndef PBD_COLLISIONS_H
#define PBD_COLLISIONS_H

#include <stdbool.h>
#include "pbd_math.h"

typedef struct {
    Vector3 normal;
    float magnitude;
} PbdFrictionConstraint;

static inline Vector3 pbd_wall_collision(Vector3 pos, float radius, float inv_mass, Vector4 plane, float omega, float eps, PbdFrictionConstraint *friction) {
    float d = pbd_v3_dot(pos, (Vector3){ plane.x, plane.y, plane.z }) - plane.w;
    float pen = radius - d;
    if (pen < 0.0f || inv_mass == 0.0f) return (Vector3){0};
    Vector3 normal = { plane.x, plane.y, plane.z };
    Vector3 dx = pbd_v3_scale(normal, omega * pen);
    if (friction && pen > friction->magnitude) {
        friction->normal = normal;
        friction->magnitude = pen;
    }
    return dx;
}

static inline void pbd_resolve_particle_pair(Vector3 *posA, float invMassA, float radiusA,
                                             Vector3 *posB, float invMassB, float radiusB,
                                             float omega, float eps)
{
    if (invMassA == 0.0f && invMassB == 0.0f) return;
    Vector3 diff = pbd_v3_sub(*posA, *posB);
    float dist = pbd_v3_length(diff);
    float target = radiusA + radiusB;
    if (dist >= target || dist <= eps) return;

    Vector3 normal = (dist > eps) ? pbd_v3_scale(diff, 1.0f / dist) : (Vector3){1.0f, 0.0f, 0.0f};
    float penetration = target - dist;
    float inv_sum = invMassA + invMassB;
    if (inv_sum <= 0.0f) return;
    float scale = omega * 0.5f * penetration / inv_sum;

    if (invMassA > 0.0f) {
        *posA = pbd_v3_add(*posA, pbd_v3_scale(normal, scale * invMassA));
    }
    if (invMassB > 0.0f) {
        *posB = pbd_v3_sub(*posB, pbd_v3_scale(normal, scale * invMassB));
    }
}

static inline Vector3 pbd_projectile_collision(Vector3 pos, float invMass, float radius, Vector4 projectile, float projectileInvMass, float omega, float eps) {
    float rj = projectile.w;
    if (rj <= 0.0f) return (Vector3){0};
    Vector3 diff = pbd_v3_sub(pos, (Vector3){ projectile.x, projectile.y, projectile.z });
    float dist = pbd_v3_length(diff) + eps;
    float pen = radius + rj - dist;
    if (pen <= 0.0f) return (Vector3){0};
    Vector3 normal = pbd_v3_scale(diff, 1.0f / dist);
    float inv_sum = invMass + projectileInvMass + eps;
    float scale = omega * 0.5f * pen / inv_sum;
    if (invMass > 0.0f) {
        return pbd_v3_scale(normal, scale * invMass);
    }
    return (Vector3){0};
}

static inline Vector3 pbd_limit_velocity(Vector3 vel, float radius, float dt) {
    float max_len = (dt > 0.0f) ? radius / dt : 0.0f;
    return pbd_v3_limit(vel, max_len);
}

#endif // PBD_COLLISIONS_H
