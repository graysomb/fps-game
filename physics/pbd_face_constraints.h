#ifndef PBD_FACE_CONSTRAINTS_H
#define PBD_FACE_CONSTRAINTS_H

#include <stdbool.h>
#include <stdint.h>
#include "pbd_math.h"

typedef struct {
    Vector3 *face0_positions;      // length 4, in/out
    Vector3 *face1_positions;      // length 4, in/out
    const float *face0_inv_mass;   // length 4
    const float *face1_inv_mass;   // length 4
    float radius;
    float strain_limit_max;        // c.mStrainLimit[0]
    float strain_limit_min;        // c.mStrainLimit[1]
    int face_index;                // which face on voxel0 (0..5)
    float alpha;
    float alpha_len;
    bool enable_fracture;
    uint32_t seed;
} PbdFaceConstraintParams;

static inline uint32_t pbd_face_hash_step(uint32_t x) {
    x ^= x >> 17;
    x *= 0xed5ad4bbU;
    x ^= x >> 11;
    x *= 0xac4c1b51U;
    x ^= x >> 15;
    x *= 0x31848babU;
    x ^= x >> 14;
    return x;
}

static inline Vector3 pbd_face_rand(uint32_t seed) {
    uint32_t hx = pbd_face_hash_step(seed + 0u);
    uint32_t hy = pbd_face_hash_step(seed + 1u);
    uint32_t hz = pbd_face_hash_step(seed + 2u);
    const float inv = 1.0f / 4294967295.0f;
    return (Vector3){ hx * inv, hy * inv, hz * inv };
}

static inline Vector3 pbd_face_corner_offset(Vector3 d0, Vector3 d1, Vector3 d2, int idx, bool positive_side) {
    Vector3 offset = pbd_v3_add((idx & 1) ? d0 : pbd_v3_scale(d0, -1.0f),
                                (idx & 2) ? d1 : pbd_v3_scale(d1, -1.0f));
    offset = pbd_v3_add(offset, positive_side ? d2 : pbd_v3_scale(d2, -1.0f));
    return offset;
}

static inline bool pbd_face_constraint_project(const PbdFaceConstraintParams *params, bool *out_broken) {
    if (out_broken) *out_broken = false;
    if (!params || !params->face0_positions || !params->face1_positions) {
        return false;
    }

    Vector3 face0[4];
    Vector3 face1[4];
    float w0[4];
    float w1[4];
    float sum0 = 0.0f;
    float sum1 = 0.0f;

    for (int i = 0; i < 4; ++i) {
        face0[i] = params->face0_positions[i];
        face1[i] = params->face1_positions[i];
        w0[i] = params->face0_inv_mass ? params->face0_inv_mass[i] : 0.0f;
        w1[i] = params->face1_inv_mass ? params->face1_inv_mass[i] : 0.0f;
        sum0 += w0[i];
        sum1 += w1[i];
    }

    const float radius = params->radius;
    const float alpha = params->alpha;
    const float alpha_len = params->alpha_len;
    const int face_ix = params->face_index;

    if (face_ix < 0 || face_ix >= 6) {
        return false;
    }

    if (params->enable_fracture) {
        Vector3 rand = pbd_face_rand(params->seed);
        for (int i = 0; i < 4; ++i) {
            Vector3 u = pbd_v3_sub(face1[i], face0[i]);
            float L = pbd_v3_length(u);
            if (radius <= 0.0f) continue;
            float strain = (L - 2.0f * radius) / (2.0f * radius);
            if (strain > params->strain_limit_max || strain < params->strain_limit_min) {
                if (out_broken) *out_broken = true;
                if (sum0 > 0.0f && w0[i] > 0.0f) {
                    face0[i] = pbd_v3_sub(face0[i], pbd_v3_scale(pbd_v3_add(u, pbd_v3_scale(rand, 0.01f)), 0.001f));
                }
                if (sum1 > 0.0f && w1[i] > 0.0f) {
                    face1[i] = pbd_v3_add(face1[i], pbd_v3_scale(pbd_v3_add(u, pbd_v3_scale(rand, 0.01f)), 0.001f));
                }
                goto write_back;
            }
        }
    }

    Vector3 dp0 = {0}, dp1 = {0}, dp2 = {0};
    Vector3 cen = {0};

    if (sum0 == 0.0f || sum1 == 0.0f) {
        for (int i = 0; i < 4; ++i) {
            cen = pbd_v3_add(cen, face0[i]);
            cen = pbd_v3_add(cen, face1[i]);
        }
        cen = pbd_v3_scale(cen, 0.125f);

        Vector3 u1 = {0};
        Vector3 u2 = {0};
        if (sum0 == 0.0f) {
            u1 = pbd_v3_add(pbd_v3_sub(face0[1], face0[0]), pbd_v3_sub(face0[3], face0[2]));
            u2 = pbd_v3_add(pbd_v3_sub(face0[2], face0[0]), pbd_v3_sub(face0[3], face0[1]));
        } else {
            u1 = pbd_v3_add(pbd_v3_sub(face1[1], face1[0]), pbd_v3_sub(face1[3], face1[2]));
            u2 = pbd_v3_add(pbd_v3_sub(face1[2], face1[0]), pbd_v3_sub(face1[3], face1[1]));
        }

        Vector3 dir = pbd_v3_cross(u1, u2);
        if (face_ix == 3) dir = pbd_v3_scale(dir, -1.0f);

        dp0 = pbd_v3_scale(pbd_v3_normalize(u1), radius);
        dp1 = pbd_v3_scale(pbd_v3_normalize(u2), radius);
        dp2 = pbd_v3_scale(pbd_v3_normalize(dir), radius);

        if (sum0 > 0.0f) {
            for (int i = 0; i < 4; ++i) {
                if (w0[i] == 0.0f) continue;
                face0[i] = pbd_v3_add(cen, pbd_face_corner_offset(dp0, dp1, dp2, i, false));
            }
        }
        if (sum1 > 0.0f) {
            for (int i = 0; i < 4; ++i) {
                if (w1[i] == 0.0f) continue;
                face1[i] = pbd_v3_add(cen, pbd_face_corner_offset(dp0, dp1, dp2, i, true));
            }
        }
    } else {
        for (int con_it = 0; con_it < 3; ++con_it) {
            dp0 = pbd_v3_add(
                pbd_v3_add(pbd_v3_sub(face1[0], face0[0]), pbd_v3_sub(face1[1], face0[1])),
                pbd_v3_add(pbd_v3_sub(face1[2], face0[2]), pbd_v3_sub(face1[3], face0[3])));

            dp1 = pbd_v3_add(
                pbd_v3_add(pbd_v3_sub(face0[1], face0[0]), pbd_v3_sub(face0[3], face0[2])),
                pbd_v3_add(pbd_v3_sub(face1[1], face1[0]), pbd_v3_sub(face1[3], face1[2])));

            dp2 = pbd_v3_add(
                pbd_v3_add(pbd_v3_sub(face0[2], face0[0]), pbd_v3_sub(face0[3], face0[1])),
                pbd_v3_add(pbd_v3_sub(face1[2], face1[0]), pbd_v3_sub(face1[3], face1[1])));

            cen = (Vector3){0,0,0};
            for (int i = 0; i < 4; ++i) {
                cen = pbd_v3_add(cen, face0[i]);
                cen = pbd_v3_add(cen, face1[i]);
            }
            cen = pbd_v3_scale(cen, 0.125f);

            Vector3 u0 = pbd_v3_sub(dp0, pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp1, dp0, 1e-12f), pbd_v3_project(dp2, dp0, 1e-12f)), alpha));
            Vector3 u1 = pbd_v3_sub(dp1, pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp0, dp1, 1e-12f), pbd_v3_project(dp2, dp1, 1e-12f)), alpha));
            Vector3 u2 = pbd_v3_sub(dp2, pbd_v3_scale(pbd_v3_add(pbd_v3_project(dp0, dp2, 1e-12f), pbd_v3_project(dp1, dp2, 1e-12f)), alpha));

            float V = pbd_v3_dot(pbd_v3_cross(u0, u1), u2);
            if (face_ix == 3) V = -V;
            if (V < 0.0f) {
                if (out_broken) *out_broken = true;
                goto write_back;
            }

            Vector3 lenu = { pbd_v3_length(u0) + 1e-12f, pbd_v3_length(u1) + 1e-12f, pbd_v3_length(u2) + 1e-12f };
            Vector3 lenp = { pbd_v3_length(dp0) + 1e-12f, pbd_v3_length(dp1) + 1e-12f, pbd_v3_length(dp2) + 1e-12f };
            float r_v = powf(radius * radius * radius / (lenp.x * lenp.y * lenp.z), 0.3333f);
            dp0 = pbd_v3_scale(u0, pbd_mixf(radius, lenp.x * r_v, alpha_len) / lenu.x);
            dp1 = pbd_v3_scale(u1, pbd_mixf(radius, lenp.y * r_v, alpha_len) / lenu.y);
            dp2 = pbd_v3_scale(u2, pbd_mixf(radius, lenp.z * r_v, alpha_len) / lenu.z);

            for (int i = 0; i < 4; ++i) {
                if (w0[i] > 0.0f) face0[i] = pbd_v3_add(cen, pbd_face_corner_offset(dp0, dp1, dp2, i, false));
                if (w1[i] > 0.0f) face1[i] = pbd_v3_add(cen, pbd_face_corner_offset(dp0, dp1, dp2, i, true));
            }
        }
    }

write_back:
    for (int i = 0; i < 4; ++i) {
        if (!params->face0_inv_mass || params->face0_inv_mass[i] > 0.0f) {
            params->face0_positions[i] = face0[i];
        }
        if (!params->face1_inv_mass || params->face1_inv_mass[i] > 0.0f) {
            params->face1_positions[i] = face1[i];
        }
    }

    return !(out_broken && *out_broken);
}

#endif // PBD_FACE_CONSTRAINTS_H
