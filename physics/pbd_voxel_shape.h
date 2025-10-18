#ifndef PBD_VOXEL_SHAPE_H
#define PBD_VOXEL_SHAPE_H

#include <stdbool.h>
#include "pbd_math.h"

typedef struct {
    Vector3 *positions;       // length 8, in/out predicted positions
    const float *inv_mass;    // length 8
    float radius;             // particle radius
    float alpha;              // uAlpha
    float alpha_len;          // uAlphaLen
    int iterations;           // iteration count (typically 3)
    float eps;                // numerical epsilon
} PbdVoxelShapeParams;

static inline int pbd_index_of_min(Vector3 v) {
    if (v.x <= v.y) {
        if (v.x <= v.z) return 0;
    } else {
        if (v.y <= v.z) return 1;
    }
    return 2;
}

static inline bool pbd_voxel_shape_project(const PbdVoxelShapeParams *params) {
    if (!params || !params->positions || !params->inv_mass) {
        return false;
    }

    float sum_w = 0.0f;
    for (int i = 0; i < 8; ++i) {
        sum_w += params->inv_mass[i];
        if (params->inv_mass[i] < 0.0f) {
            return false; // deleted particle guard
        }
    }

    if (sum_w == 0.0f) {
        return false; // rigid voxel
    }

    Vector3 p[8];
    for (int i = 0; i < 8; ++i) {
        p[i] = params->positions[i];
    }

    const float eps = (params->eps > 0.0f) ? params->eps : 1e-12f;
    const float radius = params->radius;
    const int iterations = (params->iterations > 0) ? params->iterations : 3;

    for (int iter = 0; iter < iterations; ++iter) {
        Vector3 dp[3];
        dp[0] = pbd_v3_add(pbd_v3_sub(p[1], p[0]),
                pbd_v3_add(pbd_v3_sub(p[3], p[2]),
                pbd_v3_add(pbd_v3_sub(p[5], p[4]), pbd_v3_sub(p[7], p[6]))));
        dp[1] = pbd_v3_add(pbd_v3_sub(p[2], p[0]),
                pbd_v3_add(pbd_v3_sub(p[3], p[1]),
                pbd_v3_add(pbd_v3_sub(p[6], p[4]), pbd_v3_sub(p[7], p[5]))));
        dp[2] = pbd_v3_add(pbd_v3_sub(p[4], p[0]),
                pbd_v3_add(pbd_v3_sub(p[5], p[1]),
                pbd_v3_add(pbd_v3_sub(p[6], p[2]), pbd_v3_sub(p[7], p[3]))));

        Vector3 cen = {0};
        for (int i = 0; i < 8; ++i) {
            cen = pbd_v3_add(cen, p[i]);
        }
        cen = pbd_v3_scale(cen, 0.125f);

        Vector3 u0 = pbd_v3_sub(dp[0], pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp[1], dp[0], eps), pbd_v3_project(dp[2], dp[0], eps)), params->alpha));
        Vector3 u1 = pbd_v3_sub(dp[1], pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp[0], dp[1], eps), pbd_v3_project(dp[2], dp[1], eps)), params->alpha));
        Vector3 u2 = pbd_v3_sub(dp[2], pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp[0], dp[2], eps), pbd_v3_project(dp[1], dp[2], eps)), params->alpha));

        Vector3 lenu = { pbd_v3_length(u0) + eps, pbd_v3_length(u1) + eps, pbd_v3_length(u2) + eps };

        float V = pbd_v3_dot(pbd_v3_cross(u0, u1), u2);
        if (V < 0.0f) {
            int mix_ix = pbd_index_of_min(lenu);
            if (mix_ix == 0) u0 = pbd_v3_scale(u0, -1.0f);
            if (mix_ix == 1) u1 = pbd_v3_scale(u1, -1.0f);
            if (mix_ix == 2) u2 = pbd_v3_scale(u2, -1.0f);
        }

        Vector3 lenp = { pbd_v3_length(dp[0]) + eps, pbd_v3_length(dp[1]) + eps, pbd_v3_length(dp[2]) + eps };
        float r_v = powf(radius * radius * radius / (lenp.x * lenp.y * lenp.z), 0.3333f);
        u0 = pbd_v3_scale(u0, pbd_mixf(radius, lenp.x * r_v, params->alpha_len) / lenu.x);
        u1 = pbd_v3_scale(u1, pbd_mixf(radius, lenp.y * r_v, params->alpha_len) / lenu.y);
        u2 = pbd_v3_scale(u2, pbd_mixf(radius, lenp.z * r_v, params->alpha_len) / lenu.z);

        Vector3 new_p[8];
        new_p[0] = pbd_v3_sub(pbd_v3_sub(pbd_v3_sub(cen, u0), u1), u2);
        new_p[1] = pbd_v3_sub(pbd_v3_sub(pbd_v3_add(cen, u0), u1), u2);
        new_p[2] = pbd_v3_sub(pbd_v3_add(pbd_v3_sub(cen, u0), u1), u2);
        new_p[3] = pbd_v3_sub(pbd_v3_add(pbd_v3_add(cen, u0), u1), u2);
        new_p[4] = pbd_v3_add(pbd_v3_sub(pbd_v3_sub(cen, u0), u1), u2);
        new_p[5] = pbd_v3_add(pbd_v3_sub(pbd_v3_add(cen, u0), u1), u2);
        new_p[6] = pbd_v3_add(pbd_v3_add(pbd_v3_sub(cen, u0), u1), u2);
        new_p[7] = pbd_v3_add(pbd_v3_add(pbd_v3_add(cen, u0), u1), u2);

        for (int i = 0; i < 8; ++i) {
            if (params->inv_mass[i] > 0.0f) {
                p[i] = new_p[i];
            }
        }
    }

    for (int i = 0; i < 8; ++i) {
        params->positions[i] = p[i];
    }

    return true;
}

#endif // PBD_VOXEL_SHAPE_H
