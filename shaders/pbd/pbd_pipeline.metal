#include <metal_stdlib>
using namespace metal;

struct ParticleState {
    float4 pos_radius;
    float4 prev_inv_mass;
    float4 predicted_base_inv_mass;
    float4 velocity;
};

struct VoxelState {
    uint4 particle_0_3;
    uint4 particle_4_7;
    float4 pos_rest_edge;
    float4 velocity_rest_volume;
    int4 flags;
    float4 bounds_min;
    float4 bounds_max;
    uint4 lifecycle;
};

struct StaticCollider { float4 center, bounds_min, bounds_max; };

struct FluidState {
    float4 density_lambda;
    float4 delta;
};

struct FluidNeighborData {
    uint count;
    uint neighbors[64];
};

struct GpuUniforms {
    int particle_count, sim_count, voxel_count, static_collider_count;
    int hash_size, static_hash_size, active_players, break_damp_frames;
    int mode;
    int vgs_color, fluid_count, integer_padding_2;
    float dt, voxel_size, floor_size, gravity;
    float velocity_damping, sor, collision_relaxation, vgs_alpha;
    float vgs_beta, vgs_epsilon, strain_threshold, shear_threshold;
    float tether_spring, tether_damping, rest_grid_step, padding;
    float4 players[4];
    float4 tether_targets[4];
};

static_assert(sizeof(ParticleState) == 64, "ParticleState layout mismatch");
static_assert(sizeof(VoxelState) == 128, "VoxelState layout mismatch");
static_assert(sizeof(StaticCollider) == 48, "StaticCollider layout mismatch");
static_assert(sizeof(GpuUniforms) == 240, "GpuUniforms layout mismatch");

constant int MODE_RESET = 0;
constant int MODE_INTEGRATE = 1;
constant int MODE_HASH_CLEAR = 2;
constant int MODE_HASH_BUILD = 3;
constant int MODE_SCENE_COLLISIONS = 4;
constant int MODE_PAIR_COLLISIONS = 5;
constant int MODE_APPLY = 6;
constant int MODE_VGS = 7;
constant int MODE_STATIC_COLLISIONS = 8;
constant int MODE_BREAK_MASK = 9;
constant int MODE_FINALIZE_PARTICLES = 10;
constant int MODE_FINALIZE_VOXELS = 11;
constant int MODE_SPLIT_BREAKS = 12;
constant int MODE_PREPARE_INDIRECT = 13;
constant int MODE_WAKE_GATHER = 14;
constant int MODE_WAKE_APPLY = 15;
constant int MODE_TOPOLOGY_REBUILD_SERIAL = 16;
constant int MODE_PBF_BUILD_NEIGHBORS = 17;
constant int MODE_PBF_LAMBDA = 18;
constant int MODE_PBF_DELTA = 19;
constant int MODE_PBF_APPLY = 20;
constant int MODE_PBF_STATIC_COLLISIONS = 21;
constant int MODE_PBF_DYNAMIC_SOLID_COLLISIONS = 22;
constant int MODE_PBF_DYNAMIC_SOLID_FINAL = 23;
constant int MODE_PBF_VISCOSITY = 24;
constant int MODE_PBF_APPLY_VISCOSITY = 25;
constant int CONTROL_FLAG_OVERFLOW = 1;
constant int CONTROL_FLAG_TOPOLOGY_DIRTY = 2;
constant int CONTROL_FLAG_BREAK_OCCURRED = 4;
constant uint PBF_MAX_NEIGHBORS = 64u;
constant float PARTICLE_MATERIAL_FLUID = 1.0f;

inline bool isFluid(ParticleState p) {
    return abs(p.velocity.w - PARTICLE_MATERIAL_FLUID) < 0.25f;
}

inline uint fluidId(uint gid, device uint *simId, constant GpuUniforms &u) {
    return simId[uint(u.sim_count) + gid];
}

constant int FACE_CORNERS[24] = {
    1,3,5,7, 0,2,4,6, 2,3,6,7,
    0,1,4,5, 4,5,6,7, 0,1,2,3
};

inline int controlLoad(device const int *control, int component) {
    return control[component];
}

inline uint voxelParticle(VoxelState v, int corner) {
    return corner < 4 ? v.particle_0_3[corner] : v.particle_4_7[corner - 4];
}

inline uint hashCoord(int3 c, int size) {
    uint h = uint(c.x) * 73856093u;
    h ^= uint(c.y) * 19349663u;
    h ^= uint(c.z) * 83492791u;
    return h & uint(size - 1);
}

inline void atomicAddFloat(device atomic_uint *value, float addend) {
    if (addend == 0.0f) return;
    atomic_fetch_add_explicit(reinterpret_cast<device atomic<float> *>(value), addend, memory_order_relaxed);
}

inline void accumulate(device atomic_uint *correction, device const int *control,
                       uint id, float3 delta, float weight) {
    if (weight <= 0.0f || id >= uint(controlLoad(control, 0))) return;
    device atomic_uint *base = correction + id * 4u;
    atomicAddFloat(base + 0, delta.x * weight);
    atomicAddFloat(base + 1, delta.y * weight);
    atomicAddFloat(base + 2, delta.z * weight);
    atomicAddFloat(base + 3, weight);
}

inline float3 projectOnto(float3 onto, float3 value, constant GpuUniforms &u) {
    float denom = dot(onto, onto);
    return denom < u.vgs_epsilon ? float3(0.0f) : onto * (dot(onto, value) / denom);
}

inline void resetCorrections(uint gid, device atomic_uint *correction,
                             device uint *simId, device const int *refcount,
                             device const int *control) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (refcount[id] <= 0) return;
    device float *base = reinterpret_cast<device float *>(correction) + id * 4u;
    base[0] = 0.0f;
    base[1] = 0.0f;
    base[2] = 0.0f;
    base[3] = 0.0f;
}

inline void integrateParticle(uint gid, device ParticleState *particle,
                              device int *tetherOwner, device uint *simId,
                              device const int *refcount, device const int *control,
                              constant GpuUniforms &u) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (refcount[id] <= 0) return;
    ParticleState p = particle[id];
    p.prev_inv_mass.xyz = p.pos_radius.xyz;
    p.predicted_base_inv_mass.xyz = p.pos_radius.xyz;
    if (p.prev_inv_mass.w > 0.0f) {
        p.velocity.xyz *= u.velocity_damping;
        p.predicted_base_inv_mass.xyz += p.velocity.xyz * u.dt;
        p.predicted_base_inv_mass.y -= u.gravity * u.dt * u.dt;
        int player = tetherOwner[id];
        if (player >= 0 && player < 4) {
            float scale = u.tether_targets[player].w;
            float3 accel = u.tether_targets[player].xyz * u.tether_spring * scale;
            p.predicted_base_inv_mass.xyz += accel * u.dt * u.dt;
            p.velocity.xyz += accel * u.dt;
            float damp = clamp(u.tether_damping * scale * u.dt, 0.0f, 0.9f);
            float3 disp = p.predicted_base_inv_mass.xyz - p.prev_inv_mass.xyz;
            p.predicted_base_inv_mass.xyz = p.prev_inv_mass.xyz + disp * (1.0f - damp);
            p.velocity.xyz *= 1.0f - damp;
        }
    }
    particle[id] = p;
}

inline void clearHash(uint gid, device atomic_int *hashHead, constant GpuUniforms &u) {
    if (gid < uint(u.hash_size)) atomic_store_explicit(&hashHead[gid], -1, memory_order_relaxed);
}

inline void buildHash(uint gid, device ParticleState *particle, device int4 *cell,
                      device uint *collisionId, device atomic_uint *collisionControl,
                      device atomic_int *hashHead,
                      device int *hashNext, device const int *refcount,
                      device const int *control, constant GpuUniforms &u) {
    if (gid >= atomic_load_explicit(&collisionControl[0], memory_order_relaxed)) return;
    uint id = collisionId[gid];
    if (refcount[id] <= 0) return;
    int3 c = int3(floor(particle[id].predicted_base_inv_mass.xyz / u.voxel_size));
    cell[id].xyz = c;
    uint h = hashCoord(c, u.hash_size);
    hashNext[gid] = atomic_exchange_explicit(&hashHead[h], int(gid), memory_order_relaxed);
}

inline void pairCollisions(uint gid, device ParticleState *particle,
                           device atomic_uint *correction, device int4 *cell,
                           device uint *collisionId, device atomic_uint *collisionControl,
                           device atomic_int *hashHead,
                           device int *hashNext, device int *collisionMeta,
                           device const int *refcount, device const int *control,
                           constant GpuUniforms &u) {
    int collisionCount = int(atomic_load_explicit(&collisionControl[0], memory_order_relaxed));
    if (gid >= uint(collisionCount)) return;
    uint aid = collisionId[gid];
    if (refcount[aid] <= 0) return;
    ParticleState a = particle[aid];
    if (isFluid(a)) return;
    float wa = a.prev_inv_mass.w;
    if (wa <= 0.0f) return;
    constexpr float eps = 1e-6f;

    float3 scenePos = a.predicted_base_inv_mass.xyz;
    float sceneRadius = a.pos_radius.w;
    for (int playerIndex = 0; playerIndex < u.active_players && playerIndex < 4; ++playerIndex) {
        if (u.players[playerIndex].w < 0.0f) continue;
        float halfSize = u.players[playerIndex].w;
        float3 nearest = clamp(scenePos, u.players[playerIndex].xyz - float3(halfSize),
                                         u.players[playerIndex].xyz + float3(halfSize));
        float3 sceneDelta = scenePos - nearest;
        float sceneDistSq = dot(sceneDelta, sceneDelta);
        if (sceneDistSq < sceneRadius * sceneRadius) {
            float sceneDist = sqrt(max(sceneDistSq, eps));
            float3 normal = sceneDist > eps ? -sceneDelta / sceneDist : float3(0.0f, 1.0f, 0.0f);
            scenePos += normal * (sceneRadius - sceneDist);
        }
    }
    accumulate(correction, control, aid, scenePos - a.predicted_base_inv_mass.xyz, 1.0f);

    int3 ac = cell[aid].xyz;
    for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        if (dz < 0 || (dz == 0 && dy < 0) || (dz == 0 && dy == 0 && dx < 0)) continue;
        int3 nc = ac + int3(dx, dy, dz);
        int item = atomic_load_explicit(&hashHead[hashCoord(nc, u.hash_size)], memory_order_relaxed);
        int traversed = 0;
        while (item >= 0 && item < collisionCount && traversed < collisionCount) {
            ++traversed;
            bool sameCell = dx == 0 && dy == 0 && dz == 0;
            if (!sameCell || uint(item) > gid) {
                uint bid = collisionId[item];
                if (refcount[bid] <= 0 || isFluid(particle[bid])) {
                    item = hashNext[item]; continue;
                }
                if (all(cell[bid].xyz == nc)) {
                    ParticleState b = particle[bid];
                    device int *metaA = collisionMeta + aid * 4u;
                    device int *metaB = collisionMeta + bid * 4u;
                    if (metaA[3] >= 0 && metaA[3] == metaB[3]) {
                        float3 restDelta = float3(metaA[0]-metaB[0], metaA[1]-metaB[1], metaA[2]-metaB[2]) * u.rest_grid_step;
                        float restTarget = a.pos_radius.w + b.pos_radius.w + 1e-5f;
                        if (dot(restDelta, restDelta) <= restTarget * restTarget) {
                            item = hashNext[item]; continue;
                        }
                    }
                    float wb = b.prev_inv_mass.w;
                    float wsum = wa + wb;
                    float3 delta = a.predicted_base_inv_mass.xyz - b.predicted_base_inv_mass.xyz;
                    float distSq = dot(delta, delta);
                    float target = a.pos_radius.w + b.pos_radius.w;
                    if (wsum > 0.0f && distSq < target * target) {
                        float dist = sqrt(max(distSq, eps));
                        float penetration = target - dist;
                        if (penetration > 0.0f) {
                            float3 normal = dist > eps ? delta / dist : float3(1.0f, 0.0f, 0.0f);
                            float scale = u.collision_relaxation * (0.5f * penetration) / wsum;
                            int maxTimer = max(cell[aid].w, cell[bid].w);
                            float damp = u.break_damp_frames > 0
                                ? clamp(1.0f - float(maxTimer) / float(u.break_damp_frames), 0.0f, 1.0f) : 1.0f;
                            if (wa > 0.0f) accumulate(correction, control, aid, normal * (scale * wa * damp), 1.0f);
                            if (wb > 0.0f) accumulate(correction, control, bid, normal * (-scale * wb * damp), 1.0f);
                        }
                    }
                }
            }
            item = hashNext[item];
        }
        if (item >= 0 && item < collisionCount)
            atomic_fetch_add_explicit(collisionControl + 3, 1u, memory_order_relaxed);
    }
}

inline void applyCorrections(uint gid, device ParticleState *particle,
                             device atomic_uint *correction, device uint *simId,
                             device const int *refcount, device const int *control,
                             constant GpuUniforms &u) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (refcount[id] <= 0) return;
    device float *base = reinterpret_cast<device float *>(correction) + id * 4u;
    float weight = base[3];
    if (weight > 0.0f) {
        float3 sum(base[0], base[1], base[2]);
        particle[id].predicted_base_inv_mass.xyz += sum * (u.sor / weight);
    }
    base[0] = 0.0f;
    base[1] = 0.0f;
    base[2] = 0.0f;
    base[3] = 0.0f;
}

inline void solveVgs(uint gid, device ParticleState *particle,
                     device atomic_uint *correction, device VoxelState *voxel,
                     device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.voxel_count)) return;
    VoxelState v = voxel[gid];
    if (v.flags.x == 0 || v.flags.y != 0 || v.flags.z != 0 || v.flags.w == 0) return;
    float3 p[8], original[8];
    float applyWeight[8];
    bool dynamicParticle = false;
    for (int i = 0; i < 8; ++i) {
        uint id = voxelParticle(v, i);
        if (id >= uint(controlLoad(control, 0))) return;
        p[i] = particle[id].predicted_base_inv_mass.xyz;
        original[i] = p[i];
        float w = particle[id].prev_inv_mass.w;
        applyWeight[i] = w > 0.0f ? w : 0.0f;
        dynamicParticle = dynamicParticle || w > 0.0f;
    }
    if (!dynamicParticle) return;
    float restEdge = v.pos_rest_edge.w, restVolume = v.velocity_rest_volume.w;
    for (int iteration = 0; iteration < 3; ++iteration) {
        float3 centerPos(0.0f);
        for (int i = 0; i < 8; ++i) centerPos += p[i];
        centerPos *= 0.125f;
        float3 v0=((p[1]-p[0])+(p[3]-p[2])+(p[5]-p[4])+(p[7]-p[6]))*0.25f;
        float3 v1=((p[2]-p[0])+(p[3]-p[1])+(p[6]-p[4])+(p[7]-p[5]))*0.25f;
        float3 v2=((p[4]-p[0])+(p[5]-p[1])+(p[6]-p[2])+(p[7]-p[3]))*0.25f;
        float3 u0=v0-u.vgs_alpha*(projectOnto(v1,v0,u)+projectOnto(v2,v0,u));
        float3 u1=v1-u.vgs_alpha*(projectOnto(v2,v1,u)+projectOnto(v0,v1,u));
        float3 u2=v2-u.vgs_alpha*(projectOnto(v0,v2,u)+projectOnto(v1,v2,u));
        float target0=mix(restEdge,length(v0),u.vgs_beta);
        float target1=mix(restEdge,length(v1),u.vgs_beta);
        float target2=mix(restEdge,length(v2),u.vgs_beta);
        if(length(u0)>u.vgs_epsilon)u0*=target0/length(u0);
        if(length(u1)>u.vgs_epsilon)u1*=target1/length(u1);
        if(length(u2)>u.vgs_epsilon)u2*=target2/length(u2);
        float volume=dot(cross(u0,u1),u2);
        if(fabs(volume)>u.vgs_epsilon){float scale=restVolume/volume;float root=pow(fabs(scale),1.0f/3.0f)*(scale<0.0f?-1.0f:1.0f);u0*=0.5f*root;u1*=0.5f*root;u2*=0.5f*root;}
        p[0]=centerPos-u0-u1-u2;p[1]=centerPos+u0-u1-u2;p[2]=centerPos-u0+u1-u2;p[3]=centerPos+u0+u1-u2;
        p[4]=centerPos-u0-u1+u2;p[5]=centerPos+u0-u1+u2;p[6]=centerPos-u0+u1+u2;p[7]=centerPos+u0+u1+u2;
    }
    for (int i = 0; i < 8; ++i) {
        if (applyWeight[i] > 0.0f) {
            accumulate(correction, control, voxelParticle(v, i), p[i] - original[i], applyWeight[i]);
        }
    }
}

inline int findStaticCell(int3 coord, device int4 *staticCell, constant GpuUniforms &u) {
    uint h=hashCoord(coord,u.static_hash_size);
    for(int probe=0;probe<u.static_hash_size;++probe){int4 entry=staticCell[h];if(entry.w==-1)return -1;if(all(entry.xyz==coord))return entry.w;h=(h+1u)&uint(u.static_hash_size-1);}return -1;
}

inline float3 pushOutOfBox(float3 pos,float radius,StaticCollider box,constant GpuUniforms &u){
    float3 bmin=box.bounds_min.xyz,bmax=box.bounds_max.xyz;bool inside=all(pos>=bmin)&&all(pos<=bmax);constexpr float eps=1e-6f;
    if(!inside){float3 closest=clamp(pos,bmin,bmax),delta=pos-closest;float distSq=dot(delta,delta);if(distSq>=radius*radius)return pos;float dist=sqrt(max(distSq,eps));float3 normal=dist>eps?delta/dist:float3(0,1,0);return pos+normal*max(radius-dist,0.0f);}
    float topGap=bmax.y-pos.y,bottomGap=pos.y-bmin.y,verticalBias=u.voxel_size*0.25f;bool preferTop=pos.y>=box.center.y&&topGap<=radius+verticalBias;bool preferBottom=pos.y<box.center.y&&bottomGap<=radius+verticalBias;
    if(preferTop||preferBottom){bool top=preferTop&&(!preferBottom||topGap<=bottomGap);float3 normal=top?float3(0,1,0):float3(0,-1,0);float gap=top?topGap:bottomGap;float penetration=radius-gap;if(penetration<0.0f)penetration=radius;return pos+normal*penetration;}
    float3 low=pos-bmin,high=bmax-pos;float distances[6]={low.x,high.x,low.y,high.y,low.z,high.z};int best=0;for(int i=1;i<6;++i)if(distances[i]<distances[best])best=i;float3 normal(0.0f);if(best==0)normal.x=-1;else if(best==1)normal.x=1;else if(best==2)normal.y=-1;else if(best==3)normal.y=1;else if(best==4)normal.z=-1;else normal.z=1;return pos+normal*(radius+distances[best]);
}

inline float3 pushOutOfPatch(float3 pos,float radius,StaticCollider patch){float3 closest=clamp(pos,patch.bounds_min.xyz,patch.bounds_max.xyz),normal=patch.center.xyz,delta=pos-closest;float signedDistance=dot(delta,normal);bool projectedInside;if(normal.x!=0.0f)projectedInside=pos.y>=patch.bounds_min.y&&pos.y<=patch.bounds_max.y&&pos.z>=patch.bounds_min.z&&pos.z<=patch.bounds_max.z;else if(normal.y!=0.0f)projectedInside=pos.x>=patch.bounds_min.x&&pos.x<=patch.bounds_max.x&&pos.z>=patch.bounds_min.z&&pos.z<=patch.bounds_max.z;else projectedInside=pos.x>=patch.bounds_min.x&&pos.x<=patch.bounds_max.x&&pos.y>=patch.bounds_min.y&&pos.y<=patch.bounds_max.y;if(projectedInside&&signedDistance<0.0f&&signedDistance>=-radius)return pos+normal*(radius-signedDistance);if(signedDistance<0.0f)return pos;float distanceSq=dot(delta,delta);if(distanceSq>=radius*radius)return pos;float distance=sqrt(max(distanceSq,1e-6f));float3 direction=distance>1e-6f?delta/distance:normal;return pos+direction*(radius-distance);}

inline void staticCollisionParticleDirect(uint id, device ParticleState *particle,
                                         device int4 *staticCell,
                                         device StaticCollider *staticCollider,
                                         device atomic_uint *collisionControl,
                                         constant GpuUniforms &u) {
    ParticleState p = particle[id];
    if (p.prev_inv_mass.w <= 0.0f) return;
    float radius = p.pos_radius.w;
    float3 pos = p.predicted_base_inv_mass.xyz;
    float terrainLimit = u.floor_size - radius, floorLimit = max(0.0f, 0.5f * u.voxel_size - radius);
    bool floorContact = pos.y < floorLimit;
    pos.y = max(pos.y, floorLimit);
    if (floorContact) {
        p.prev_inv_mass.xz = pos.xz - (pos.xz - p.prev_inv_mass.xz) * 0.05f;
        p.prev_inv_mass.y = pos.y;
    }
    pos.xz = clamp(pos.xz, float2(-terrainLimit), float2(terrainLimit));
    constexpr float eps = 1e-6f;
    for (int i = 0; i < u.active_players && i < 4; ++i) {
        if (u.players[i].w < 0.0f) continue;
        float halfSize = u.players[i].w;
        float3 nearest = clamp(pos, u.players[i].xyz - float3(halfSize), u.players[i].xyz + float3(halfSize));
        float3 delta = pos - nearest;
        float distSq = dot(delta, delta);
        if (distSq < radius * radius) {
            float dist = sqrt(max(distSq, eps));
            float3 normal = dist > eps ? delta / dist : float3(0, 1, 0);
            pos += normal * (radius - dist);
        }
    }
    bool surfaceMode = atomic_load_explicit(&collisionControl[4], memory_order_relaxed) != 0u;
    int3 center = int3(floor(pos / u.voxel_size));
    int seen[128]; int seenCount = 0;
    for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        int item = findStaticCell(center + int3(dx, dy, dz), staticCell, u);
        if (!surfaceMode && item <= -2) {
            int colliderId = -item - 2;
            if (colliderId >= 0 && colliderId < u.static_collider_count)
                pos = pushOutOfBox(pos, 0.25f * u.voxel_size, staticCollider[colliderId], u);
            continue;
        }
        while (item >= 0 && item < u.static_collider_count) {
            StaticCollider patch = staticCollider[item];
            int patchId = as_type<int>(patch.bounds_min.w);
            bool duplicate = false;
            for (int s = 0; s < seenCount; ++s) duplicate = duplicate || seen[s] == patchId;
            if (!duplicate && seenCount < 128) {
                seen[seenCount++] = patchId;
                pos = pushOutOfPatch(pos, 0.25f * u.voxel_size, patch);
            }
            item = as_type<int>(patch.center.w);
        }
    }
    center = int3(floor(pos / u.voxel_size));
    int recovery = findStaticCell(center, staticCell, u);
    if (recovery <= -2) {
        int colliderId = -recovery - 2;
        if (colliderId >= 0 && colliderId < u.static_collider_count)
            pos = pushOutOfBox(pos, 0.25f * u.voxel_size, staticCollider[colliderId], u);
    }
    p.predicted_base_inv_mass.xyz = pos;
    particle[id] = p;
}

inline void staticCollisions(uint gid, device ParticleState *particle, device uint *collisionId,
                            device atomic_uint *collisionControl, device int4 *staticCell,
                            device StaticCollider *staticCollider, device const int *refcount,
                            device const int *control, constant GpuUniforms &u) {
    if (gid >= atomic_load_explicit(&collisionControl[0], memory_order_relaxed)) return;
    uint id = collisionId[gid];
    if (refcount[id] <= 0 || isFluid(particle[id])) return;
    staticCollisionParticleDirect(id, particle, staticCell, staticCollider, collisionControl, u);
}

inline int topologyNeighbor(device int4 *topology,int voxelId,int face);

inline void gatherBreakMask(uint gid, device ParticleState *particle,
                            device VoxelState *voxel, device int4 *topology,
                            device const int *refcount,
                            device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.voxel_count)) return;
    VoxelState v = voxel[gid];
    if (v.flags.x == 0 || v.flags.y != 0 || v.flags.z != 0 || v.pos_rest_edge.w <= 0.0f || (v.lifecycle.w & 0x10u) == 0u) {
        voxel[gid] = v;
        return;
    }
    float3 p[8];
    for (int i = 0; i < 8; ++i) {
        uint id = voxelParticle(v, i);
        if (id >= uint(controlLoad(control, 0))) return;
        p[i] = particle[id].predicted_base_inv_mass.xyz;
    }
    float3 axis0 = ((p[1]-p[0])+(p[3]-p[2])+(p[5]-p[4])+(p[7]-p[6]))*0.25f;
    float3 axis1 = ((p[2]-p[0])+(p[3]-p[1])+(p[6]-p[4])+(p[7]-p[5]))*0.25f;
    float3 axis2 = ((p[4]-p[0])+(p[5]-p[1])+(p[6]-p[2])+(p[7]-p[3]))*0.25f;
    float3 lengths(length(axis0), length(axis1), length(axis2));
    float3 strain = fabs(lengths-float3(v.pos_rest_edge.w))/v.pos_rest_edge.w;
    uint glued = v.lifecycle.z, mask = 0u;
    bool exceeded = false;
    if (strain.x > u.strain_threshold) { exceeded=true; mask|=glued&3u; }
    if (strain.y > u.strain_threshold) { exceeded=true; mask|=glued&12u; }
    if (strain.z > u.strain_threshold) { exceeded=true; mask|=glued&48u; }
    float3 n0 = lengths.x > u.vgs_epsilon ? axis0/lengths.x : float3(0);
    float3 n1 = lengths.y > u.vgs_epsilon ? axis1/lengths.y : float3(0);
    float3 n2 = lengths.z > u.vgs_epsilon ? axis2/lengths.z : float3(0);
    if (fabs(dot(n0,n1)) > u.shear_threshold) { exceeded=true; mask|=glued&48u; }
    if (fabs(dot(n0,n2)) > u.shear_threshold) { exceeded=true; mask|=glued&12u; }
    if (fabs(dot(n1,n2)) > u.shear_threshold) { exceeded=true; mask|=glued&3u; }
    constexpr float hingeCosLimit = 0.9396926208f;
    for (int face = 0; face < 6; ++face) {
        if ((glued & (1u << face)) == 0u) continue;
        int sharedCornerCount = 0;
        for (int c = 0; c < 4; ++c) {
            uint sharedId = voxelParticle(v, FACE_CORNERS[face*4+c]);
            if (sharedId < uint(controlLoad(control, 0)) &&
                refcount[sharedId] > 1)
                ++sharedCornerCount;
        }
        if (sharedCornerCount > 1) continue;
        int neighborId = topologyNeighbor(topology, int(gid), face);
        if (neighborId < 0 || neighborId >= u.voxel_count) continue;
        VoxelState neighbor = voxel[neighborId];
        if (neighbor.flags.x == 0 || neighbor.flags.y != 0 || neighbor.flags.z != 0) continue;
        float3 q[8];
        for (int i = 0; i < 8; ++i) {
            uint id = voxelParticle(neighbor, i);
            if (id >= uint(controlLoad(control, 0))) return;
            q[i] = particle[id].predicted_base_inv_mass.xyz;
        }
        float3 m0=normalize(((q[1]-q[0])+(q[3]-q[2])+(q[5]-q[4])+(q[7]-q[6]))*0.25f);
        float3 m1=normalize(((q[2]-q[0])+(q[3]-q[1])+(q[6]-q[4])+(q[7]-q[5]))*0.25f);
        float3 m2=normalize(((q[4]-q[0])+(q[5]-q[1])+(q[6]-q[2])+(q[7]-q[3]))*0.25f);
        if (dot(n0,m0)<hingeCosLimit || dot(n1,m1)<hingeCosLimit || dot(n2,m2)<hingeCosLimit) {
            exceeded = true;
            mask |= 1u << face;
        }
    }
    if (exceeded) {
        v.flags.x = 0; // Deactivate VGS constraint
        v.lifecycle.y = 1u; // wake_source = true
        v.lifecycle.z = 0u; // Clear glued faces
        v.lifecycle.w &= ~0x8u; // No longer a fluid/solid collision volume
    }
    voxel[gid] = v;
}

inline void splitBrokenFaces(uint gid,device ParticleState *particle,device atomic_uint *correction,device int4 *cell,device uint *simId,device uint *collisionId,device atomic_int *collisionMember,device atomic_uint *collisionControl,device VoxelState *voxel,device int *tetherOwner,device int *collisionMeta,device const int *refcount,device const int *control,device int *cloneParent,constant GpuUniforms &u){
    // No-op in VGS-as-glue model: fracture never splits or clones particles.
    (void)gid; (void)particle; (void)correction; (void)cell; (void)simId;
    (void)collisionId; (void)collisionMember; (void)collisionControl;
    (void)voxel; (void)tetherOwner; (void)collisionMeta; (void)refcount;
    (void)control; (void)cloneParent; (void)u;
}

inline int topologyNeighbor(device int4 *topology,int voxelId,int face){return face<4?topology[voxelId*2][face]:topology[voxelId*2+1][face-4];}
inline int topologyRoot(device int4 *topology,int voxelId,constant GpuUniforms &u){int root=voxelId;for(int i=0;i<32;++i){int parent=topology[root*2+1].z;if(parent==root||parent<0||parent>=u.voxel_count)break;root=parent;}return root;}
inline void topologyUnion(device VoxelState *voxel,device int4 *topology,int gid,constant GpuUniforms &u){VoxelState v=voxel[gid];if((v.lifecycle.z&0x80000000u)==0u)return;for(int face=0;face<6;++face){if((v.lifecycle.z&(1u<<face))==0u)continue;int neighbor=topologyNeighbor(topology,gid,face);if(neighbor<0||neighbor>=u.voxel_count||(voxel[neighbor].lifecycle.z&0x80000000u)==0u)continue;int a=topologyRoot(topology,gid,u),b=topologyRoot(topology,neighbor,u);if(a!=b)topology[max(a,b)*2+1].z=min(a,b);}}
inline void assignGroup(device int *collisionMeta,uint id,int group){uint offset=id*4u+3u;int old=collisionMeta[offset];if(old==group||old==-2)return;collisionMeta[offset]=old==-1?group:-2;}
inline void topologyRebuild(uint gid,device ParticleState *particle,device uint *simId,device VoxelState *voxel,device int *collisionMeta,device const int *refcount,device const int *control,device int4 *topology,constant GpuUniforms &u){
    (void)gid; (void)particle; (void)simId; (void)voxel; (void)collisionMeta; (void)refcount; (void)control; (void)topology; (void)u;
}

inline void wakeGather(uint gid,device VoxelState *voxel,device int4 *topology,constant GpuUniforms &u){if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];if(v.flags.x==0||v.flags.y!=0)return;bool wake=v.lifecycle.y!=0u;if(!wake)for(int face=0;face<6;++face){int neighbor=topologyNeighbor(topology,int(gid),face);if(neighbor>=0&&neighbor<u.voxel_count&&voxel[neighbor].flags.y==0&&voxel[neighbor].lifecycle.y!=0u){wake=true;break;}}int timer=max(int(v.bounds_min.w)-1,0);if(wake)timer=30;v.bounds_min.w=float(timer);voxel[gid]=v;}
inline void wakeApply(uint gid,device VoxelState *voxel,constant GpuUniforms &u){(void)gid; (void)voxel; (void)u;}
inline void prepareIndirect(uint gid,device const int *control,device atomic_uint *collisionControl,device uint *dispatchArgs){if(gid!=0u)return;dispatchArgs[0]=(uint(max(controlLoad(control,1),0))+127u)/128u;dispatchArgs[1]=1u;dispatchArgs[2]=1u;dispatchArgs[3]=0u;dispatchArgs[4]=(atomic_load_explicit(&collisionControl[0], memory_order_relaxed)+127u)/128u;dispatchArgs[5]=1u;dispatchArgs[6]=1u;dispatchArgs[7]=0u;}
inline void finalizeParticle(uint gid,device ParticleState *particle,device int4 *cell,device uint *simId,device const int *refcount,device const int *control,constant GpuUniforms &u){if(gid>=uint(controlLoad(control,1)))return;uint id=simId[gid];if(refcount[id]<=0)return;ParticleState p=particle[id];float3 delta=p.predicted_base_inv_mass.xyz-p.prev_inv_mass.xyz;p.velocity.xyz=p.prev_inv_mass.w>0.0f&&u.dt>0.0f?delta/u.dt:float3(0);p.pos_radius.xyz=p.predicted_base_inv_mass.xyz;if(cell[id].w>0)cell[id].w--;particle[id]=p;}
inline void finalizeVoxel(uint gid,device ParticleState *particle,device VoxelState *voxel,device const int *control,constant GpuUniforms &u){if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];if(v.flags.x==0||v.flags.y!=0||v.flags.z!=0)return;float3 center(0),previous(0);for(int i=0;i<8;++i){uint id=voxelParticle(v,i);if(id>=uint(controlLoad(control,0)))return;center+=particle[id].predicted_base_inv_mass.xyz;previous+=particle[id].prev_inv_mass.xyz;}center*=0.125f;previous*=0.125f;if(v.lifecycle.x==0u&&u.dt>0.0f)v.velocity_rest_volume.xyz=(center-previous)/u.dt;else if(v.lifecycle.x>0u)v.lifecycle.x--;v.pos_rest_edge.xyz=center;voxel[gid]=v;}

inline float pbfPoly6(float distanceSq, float h, float poly6Coeff) {
    float hSq = h * h;
    if (distanceSq < 0.0f || distanceSq >= hSq) return 0.0f;
    float term = hSq - distanceSq;
    return poly6Coeff * term * term * term;
}

inline float3 pbfSpikyGradient(float3 separation, float h, float spikyCoeff) {
    float distanceSq = dot(separation, separation);
    if (distanceSq <= 1e-12f || distanceSq >= h * h) return float3(0.0f);
    float distance = sqrt(distanceSq);
    float term = h - distance;
    return separation * (spikyCoeff * term * term / distance);
}

inline void pbfBuildNeighbors(uint gid, device ParticleState *particle, device int4 *cell,
                              device uint *simId, device uint *collisionId,
                              device atomic_int *hashHead, device int *hashNext,
                              device FluidNeighborData *fluidNeighbors,
                              device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    int3 currentCell = cell[id].xyz;
    uint count = 0u;
    float h = u.voxel_size;
    float hSq = h * h;
    for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        int3 neighborCell = currentCell + int3(dx, dy, dz);
        int item = atomic_load_explicit(&hashHead[hashCoord(neighborCell, u.hash_size)], memory_order_relaxed);
        int traversed = 0;
        while (item >= 0 && traversed < u.particle_count) {
            ++traversed;
            uint candidate = collisionId[item];
            if (isFluid(particle[candidate]) && all(cell[candidate].xyz == neighborCell)) {
                float3 sep = particle[id].predicted_base_inv_mass.xyz - particle[candidate].predicted_base_inv_mass.xyz;
                if (dot(sep, sep) < hSq) {
                    if (count < PBF_MAX_NEIGHBORS) {
                        fluidNeighbors[gid].neighbors[count] = candidate;
                    }
                    count++;
                }
            }
            item = hashNext[item];
        }
    }
    fluidNeighbors[gid].count = count;
}

inline void pbfLambda(uint gid, device ParticleState *particle, device int4 *cell,
                      device uint *simId, device uint *collisionId,
                      device atomic_int *hashHead, device int *hashNext,
                      device FluidState *fluid, device FluidNeighborData *fluidNeighbors,
                      device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    ParticleState current = particle[id];
    float h = u.voxel_size;
    float spacing = u.voxel_size * 0.5f;
    float restDensity = 1000.0f;
    float particleMass = restDensity * spacing * spacing * spacing;
    float h2 = h * h, h6 = h2 * h2 * h2, h9 = h6 * h2 * h;
    float poly6Coeff = 315.0f / (64.0f * 3.14159265358979323846f * h9);
    float spikyCoeff = -45.0f / (3.14159265358979323846f * h6);
    float volume = particleMass / restDensity;
    float lambdaEpsilon = 1e-6f;

    float density = 0.0f;
    float3 gradientI = float3(0.0f);
    float gradientSumSq = 0.0f;
    uint cachedCount = fluidNeighbors[gid].count;
    if (cachedCount <= PBF_MAX_NEIGHBORS) {
        for (uint k = 0; k < cachedCount; ++k) {
            uint neighbor = fluidNeighbors[gid].neighbors[k];
            float3 sep = current.predicted_base_inv_mass.xyz - particle[neighbor].predicted_base_inv_mass.xyz;
            density += particleMass * pbfPoly6(dot(sep, sep), h, poly6Coeff);
            if (neighbor != id) {
                float3 gradientJ = -volume * pbfSpikyGradient(sep, h, spikyCoeff);
                gradientSumSq += dot(gradientJ, gradientJ);
                gradientI -= gradientJ;
            }
        }
    } else {
        int3 currentCell = cell[id].xyz;
        for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            int3 neighborCell = currentCell + int3(dx, dy, dz);
            int item = atomic_load_explicit(&hashHead[hashCoord(neighborCell, u.hash_size)], memory_order_relaxed);
            int traversed = 0;
            while (item >= 0 && traversed < u.particle_count) {
                ++traversed;
                uint neighbor = collisionId[item];
                if (isFluid(particle[neighbor]) && all(cell[neighbor].xyz == neighborCell)) {
                    float3 sep = current.predicted_base_inv_mass.xyz - particle[neighbor].predicted_base_inv_mass.xyz;
                    float distSq = dot(sep, sep);
                    if (distSq < h * h) {
                        density += particleMass * pbfPoly6(distSq, h, poly6Coeff);
                        if (neighbor != id) {
                            float3 gradientJ = -volume * pbfSpikyGradient(sep, h, spikyCoeff);
                            gradientSumSq += dot(gradientJ, gradientJ);
                            gradientI -= gradientJ;
                        }
                    }
                }
                item = hashNext[item];
            }
        }
    }
    gradientSumSq += dot(gradientI, gradientI);
    float constraint = density / restDensity - 1.0f;
    fluid[id].density_lambda.x = density;
    fluid[id].density_lambda.y = -constraint / (gradientSumSq + lambdaEpsilon);
}

inline void pbfDelta(uint gid, device ParticleState *particle, device int4 *cell,
                     device uint *simId, device uint *collisionId,
                     device atomic_int *hashHead, device int *hashNext,
                     device FluidState *fluid, device FluidNeighborData *fluidNeighbors,
                     device int4 *staticCell, device StaticCollider *staticCollider,
                     device atomic_uint *collisionControl,
                     device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    ParticleState current = particle[id];
    float h = u.voxel_size;
    float spacing = u.voxel_size * 0.5f;
    float restDensity = 1000.0f;
    float particleMass = restDensity * spacing * spacing * spacing;
    float h2 = h * h, h6 = h2 * h2 * h2, h9 = h6 * h2 * h;
    float poly6Coeff = 315.0f / (64.0f * 3.14159265358979323846f * h9);
    float spikyCoeff = -45.0f / (3.14159265358979323846f * h6);
    float volume = particleMass / restDensity;
    float scorrK = 0.2f;
    float scorrDeltaQ = h * 0.2f;
    float referenceKernel = pbfPoly6(scorrDeltaQ * scorrDeltaQ, h, poly6Coeff);

    float3 delta = float3(0.0f);
    uint cachedCount = fluidNeighbors[gid].count;
    if (cachedCount <= PBF_MAX_NEIGHBORS) {
        for (uint k = 0; k < cachedCount; ++k) {
            uint neighbor = fluidNeighbors[gid].neighbors[k];
            if (neighbor == id) continue;
            float3 sep = current.predicted_base_inv_mass.xyz - particle[neighbor].predicted_base_inv_mass.xyz;
            float wVal = pbfPoly6(dot(sep, sep), h, poly6Coeff);
            if (wVal <= 0.0f) continue;
            float ratio = (referenceKernel > 0.0f) ? wVal / referenceKernel : 0.0f;
            float ratioSq = ratio * ratio;
            float scorr = -scorrK * (0.25f * u.voxel_size * u.voxel_size) * (ratioSq * ratioSq);
            float scale = volume * (fluid[id].density_lambda.y + fluid[neighbor].density_lambda.y + scorr);
            delta += scale * pbfSpikyGradient(sep, h, spikyCoeff);
        }
    } else {
        int3 currentCell = cell[id].xyz;
        for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            int3 neighborCell = currentCell + int3(dx, dy, dz);
            int item = atomic_load_explicit(&hashHead[hashCoord(neighborCell, u.hash_size)], memory_order_relaxed);
            int traversed = 0;
            while (item >= 0 && traversed < u.particle_count) {
                ++traversed;
                uint neighbor = collisionId[item];
                if (neighbor != id && isFluid(particle[neighbor]) && all(cell[neighbor].xyz == neighborCell)) {
                    float3 sep = current.predicted_base_inv_mass.xyz - particle[neighbor].predicted_base_inv_mass.xyz;
                    float wVal = pbfPoly6(dot(sep, sep), h, poly6Coeff);
                    if (wVal > 0.0f) {
                        float ratio = (referenceKernel > 0.0f) ? wVal / referenceKernel : 0.0f;
                        float ratioSq = ratio * ratio;
                        float scorr = -scorrK * (0.25f * u.voxel_size * u.voxel_size) * (ratioSq * ratioSq);
                        float scale = volume * (fluid[id].density_lambda.y + fluid[neighbor].density_lambda.y + scorr);
                        delta += scale * pbfSpikyGradient(sep, h, spikyCoeff);
                    }
                }
                item = hashNext[item];
            }
        }
    }
    float deltaLength = length(delta);
    float deltaLimit = 0.0125f * u.voxel_size;
    if (deltaLength > deltaLimit) delta *= deltaLimit / deltaLength;
    fluid[id].delta = float4(0.0f);
    particle[id].predicted_base_inv_mass.xyz += delta;
    staticCollisionParticleDirect(id, particle, staticCell, staticCollider, collisionControl, u);
}

inline void pbfApply(uint gid, device ParticleState *particle,
                     device FluidState *fluid, device uint *simId,
                     constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    particle[id].predicted_base_inv_mass.xyz += fluid[id].delta.xyz;
    fluid[id].delta = float4(0.0f);
}

inline void pbfStaticCollisions(uint gid, device ParticleState *particle,
                                device int4 *staticCell,
                                device StaticCollider *staticCollider,
                                device atomic_uint *collisionControl,
                                device uint *simId, device const int *control,
                                constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    staticCollisionParticleDirect(id, particle, staticCell, staticCollider, collisionControl, u);
}

inline void pbfApplyAndStaticCollisions(uint gid, device ParticleState *particle,
                                       device FluidState *fluid,
                                       device int4 *staticCell,
                                       device StaticCollider *staticCollider,
                                       device atomic_uint *collisionControl,
                                       device uint *simId, device const int *control,
                                       constant GpuUniforms &u) {
    pbfApply(gid, particle, fluid, simId, u);
    pbfStaticCollisions(gid, particle, staticCell, staticCollider, collisionControl, simId, control, u);
}

inline void pbfDynamicSolidCollisions(uint gid, bool react,
                                      device ParticleState *particle,
                                      device atomic_uint *correction,
                                      device VoxelState *voxel,
                                      device uint *simId,
                                      device const int *control,
                                      constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    float radius = particle[id].pos_radius.w;
    float3 position = particle[id].predicted_base_inv_mass.xyz;
    float3 previous = particle[id].prev_inv_mass.xyz;
    const int4 faceCorner[6] = {
        int4(0, 2, 6, 4), int4(1, 5, 7, 3),
        int4(0, 4, 5, 1), int4(2, 3, 7, 6),
        int4(0, 1, 3, 2), int4(4, 6, 7, 5)
    };

    for (int voxelId = 0; voxelId < u.voxel_count; ++voxelId) {
        VoxelState solid = voxel[voxelId];
        if (solid.flags.y != 0 || solid.flags.z != 0 ||
            (solid.lifecycle.w & 0x8u) == 0u) continue;

        float maxReach = u.voxel_size * 1.5f + radius;
        float3 toPos = fabs(position - solid.pos_rest_edge.xyz);
        float3 toPrev = fabs(previous - solid.pos_rest_edge.xyz);
        if ((toPos.x > maxReach && toPrev.x > maxReach) ||
            (toPos.z > maxReach && toPrev.z > maxReach) ||
            (toPos.y > maxReach + u.voxel_size && toPrev.y > maxReach + u.voxel_size)) {
            continue;
        }

        float3 corner[8];
        float3 boundsMin(1e30f);
        float3 boundsMax(-1e30f);
        bool valid = true;
        for (int i = 0; i < 8; ++i) {
            uint cornerId = voxelParticle(solid, i);
            if (cornerId >= uint(controlLoad(control, 0))) { valid = false; break; }
            corner[i] = particle[cornerId].predicted_base_inv_mass.xyz;
            boundsMin = min(boundsMin, corner[i]);
            boundsMax = max(boundsMax, corner[i]);
        }
        if (!valid) continue;
        float3 expandedMin = boundsMin - float3(radius);
        float3 expandedMax = boundsMax + float3(radius);
        uint glueMask = solid.lifecycle.z;
        bool recoverBelow = !(glueMask & (1u << 2)) &&
                            position.x >= expandedMin.x && position.x <= expandedMax.x &&
                            position.z >= expandedMin.z && position.z <= expandedMax.z &&
                            position.y < expandedMin.y &&
                            position.y >= expandedMin.y - max(2.0f * radius, u.voxel_size);
        if (!recoverBelow && (
            (position.x < boundsMin.x - radius && previous.x < boundsMin.x - radius) ||
            (position.y < boundsMin.y - radius && previous.y < boundsMin.y - radius) ||
            (position.z < boundsMin.z - radius && previous.z < boundsMin.z - radius) ||
            (position.x > boundsMax.x + radius && previous.x > boundsMax.x + radius) ||
            (position.y > boundsMax.y + radius && previous.y > boundsMax.y + radius) ||
            (position.z > boundsMax.z + radius && previous.z > boundsMax.z + radius))) continue;

        bool currentInside = all(position >= expandedMin) && all(position <= expandedMax);
        bool previousInside = all(previous >= expandedMin) && all(previous <= expandedMax);
        int nearestFace = -1;
        float entryT = 0.0f;
        float exitT = 1.0f;
        float3 movement = position - previous;
        bool segmentHits = true;
        for (int axis = 0; axis < 3; ++axis) {
            float origin = previous[axis];
            float direction = movement[axis];
            if (fabs(direction) <= 1e-8f) {
                if (origin < expandedMin[axis] || origin > expandedMax[axis]) {
                    segmentHits = false;
                    break;
                }
                continue;
            }
            float tMin = (expandedMin[axis] - origin) / direction;
            float tMax = (expandedMax[axis] - origin) / direction;
            int minFace = axis * 2;
            int maxFace = minFace + 1;
            if (tMin > tMax) {
                float swapT = tMin; tMin = tMax; tMax = swapT;
                int swapFace = minFace; minFace = maxFace; maxFace = swapFace;
            }
            if (tMin > entryT) {
                entryT = tMin;
                nearestFace = minFace;
            }
            exitT = min(exitT, tMax);
            if (entryT > exitT) {
                segmentHits = false;
                break;
            }
        }
        if (!recoverBelow && !currentInside &&
            (!segmentHits || previousInside || entryT < 0.0f || entryT > 1.0f)) {
            continue;
        }
        if (recoverBelow) {
            nearestFace = 3;
        } else if (nearestFace < 0 || previousInside) {
            float3 sideSample = previousInside ? previous : position;
            float distances[6] = {
                fabs(sideSample.x - expandedMin.x), fabs(expandedMax.x - sideSample.x),
                fabs(sideSample.y - expandedMin.y), fabs(expandedMax.y - sideSample.y),
                fabs(sideSample.z - expandedMin.z), fabs(expandedMax.z - sideSample.z)
            };
            float nearest = 1e30f;
            for (int face = 0; face < 6; ++face) {
                if ((glueMask & (1u << uint(face))) != 0) continue;
                if (face == 2 && (previous.y >= (expandedMin.y + expandedMax.y) * 0.5f || movement.y <= 0.0f)) {
                    continue;
                }
                if (face == 3 && (previous.y < (expandedMin.y + expandedMax.y) * 0.5f && movement.y > 0.0f)) {
                    continue;
                }
                if (distances[face] < nearest) {
                    nearest = distances[face];
                    nearestFace = face;
                }
            }
            if (nearestFace < 0) {
                nearestFace = (movement.y <= 0.0f && !(glueMask & (1u << 3))) ? 3 : 2;
            }
        }
        if (nearestFace < 0) continue;

        float3 correctionDelta(0.0f);
        if (nearestFace == 0) correctionDelta.x = expandedMin.x - position.x - 1e-4f;
        else if (nearestFace == 1) correctionDelta.x = expandedMax.x - position.x + 1e-4f;
        else if (nearestFace == 2) correctionDelta.y = expandedMin.y - position.y - 1e-4f;
        else if (nearestFace == 3) correctionDelta.y = expandedMax.y - position.y + 1e-4f;
        else if (nearestFace == 4) correctionDelta.z = expandedMin.z - position.z - 1e-4f;
        else correctionDelta.z = expandedMax.z - position.z + 1e-4f;

        position += correctionDelta;
        if (nearestFace <= 1) previous.x = position.x;
        else if (nearestFace <= 3) previous.y = position.y;
        else previous.z = position.z;

        if (react) {
            int4 fc = faceCorner[nearestFace];
            float3 cornerReaction = -correctionDelta * (0.02f * 0.25f);
            float maxReaction = 0.05f * u.voxel_size;
            cornerReaction = clamp(cornerReaction, float3(-maxReaction), float3(maxReaction));
            for (int i = 0; i < 4; ++i) {
                uint cornerId = voxelParticle(solid, fc[i]);
                if (cornerId < uint(controlLoad(control, 0)) && particle[cornerId].prev_inv_mass.w > 0.0f) {
                    accumulate(correction, control, cornerId, cornerReaction, 1.0f);
                }
            }
        }
    }
    particle[id].predicted_base_inv_mass.xyz = position;
    particle[id].prev_inv_mass.xyz = previous;
}

inline void pbfViscosity(uint gid, device ParticleState *particle, device int4 *cell,
                         device uint *simId, device uint *collisionId,
                         device atomic_int *hashHead, device int *hashNext,
                         device FluidState *fluid, device FluidNeighborData *fluidNeighbors,
                         device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    ParticleState current = particle[id];
    float h = u.voxel_size;
    float spacing = u.voxel_size * 0.5f;
    float restDensity = 1000.0f;
    float particleMass = restDensity * spacing * spacing * spacing;
    float h2 = h * h, h6 = h2 * h2 * h2, h9 = h6 * h2 * h;
    float poly6Coeff = 315.0f / (64.0f * 3.14159265358979323846f * h9);
    float viscosity = 0.15f;

    float3 velocityDelta = float3(0.0f);
    uint cachedCount = fluidNeighbors[gid].count;
    if (cachedCount <= PBF_MAX_NEIGHBORS) {
        for (uint k = 0; k < cachedCount; ++k) {
            uint neighbor = fluidNeighbors[gid].neighbors[k];
            if (neighbor == id) continue;
            float3 sep = current.pos_radius.xyz - particle[neighbor].pos_radius.xyz;
            float wVal = pbfPoly6(dot(sep, sep), h, poly6Coeff);
            float neighborDensity = max(fluid[neighbor].density_lambda.x, restDensity * 0.1f);
            float scale = viscosity * particleMass * wVal / neighborDensity;
            velocityDelta += (particle[neighbor].velocity.xyz - current.velocity.xyz) * scale;
        }
    } else {
        int3 currentCell = cell[id].xyz;
        for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            int3 neighborCell = currentCell + int3(dx, dy, dz);
            int item = atomic_load_explicit(&hashHead[hashCoord(neighborCell, u.hash_size)], memory_order_relaxed);
            int traversed = 0;
            while (item >= 0 && traversed < u.particle_count) {
                ++traversed;
                uint neighbor = collisionId[item];
                if (neighbor != id && isFluid(particle[neighbor]) && all(cell[neighbor].xyz == neighborCell)) {
                    float3 sep = current.pos_radius.xyz - particle[neighbor].pos_radius.xyz;
                    float wVal = pbfPoly6(dot(sep, sep), h, poly6Coeff);
                    float neighborDensity = max(fluid[neighbor].density_lambda.x, restDensity * 0.1f);
                    float scale = viscosity * particleMass * wVal / neighborDensity;
                    velocityDelta += (particle[neighbor].velocity.xyz - current.velocity.xyz) * scale;
                }
                item = hashNext[item];
            }
        }
    }
    fluid[id].delta = float4(velocityDelta, 0.0f);
}

inline void pbfApplyViscosity(uint gid, device ParticleState *particle,
                              device FluidState *fluid, device uint *simId,
                              device const int *control, constant GpuUniforms &u) {
    if (gid >= uint(u.fluid_count)) return;
    uint id = fluidId(gid, simId, u);
    particle[id].velocity.xyz += fluid[id].delta.xyz;
    fluid[id].delta = float4(0.0f);
}

kernel void pbd_pipeline(
    device ParticleState *particle [[buffer(0)]],
    device atomic_uint *correction [[buffer(1)]],
    device int4 *cell [[buffer(2)]],
    device uint *simId [[buffer(3)]],
    device VoxelState *voxel [[buffer(4)]],
    device atomic_int *hashHead [[buffer(5)]],
    device int *hashNext [[buffer(6)]],
    device int4 *staticCell [[buffer(7)]],
    device int *tetherOwner [[buffer(8)]],
    device int *collisionMeta [[buffer(9)]],
    device StaticCollider *staticCollider [[buffer(10)]],
    device const int *refcount [[buffer(11)]],
    device const int *control [[buffer(12)]],
    device int *cloneParent [[buffer(13)]],
    device uint *dispatchArgs [[buffer(14)]],
    device int4 *topology [[buffer(15)]],
    device FluidState *fluid [[buffer(16)]],
    device FluidNeighborData *fluidNeighbors [[buffer(17)]],
    constant GpuUniforms &u [[buffer(18)]],
    uint gid [[thread_position_in_grid]]) {
    device uint *collisionId=simId+uint(controlLoad(control,2));
    device atomic_int *collisionMember=reinterpret_cast<device atomic_int *>(cloneParent);
    device atomic_uint *collisionControl=reinterpret_cast<device atomic_uint *>(dispatchArgs+8);
    switch(u.mode){
        case MODE_RESET:resetCorrections(gid,correction,simId,refcount,control);break;
        case MODE_INTEGRATE:integrateParticle(gid,particle,tetherOwner,simId,refcount,control,u);break;
        case MODE_HASH_CLEAR:clearHash(gid,hashHead,u);break;
        case MODE_HASH_BUILD:buildHash(gid,particle,cell,collisionId,collisionControl,hashHead,hashNext,refcount,control,u);break;
        case MODE_PAIR_COLLISIONS:pairCollisions(gid,particle,correction,cell,collisionId,collisionControl,hashHead,hashNext,collisionMeta,refcount,control,u);break;
        case MODE_APPLY:applyCorrections(gid,particle,correction,simId,refcount,control,u);break;
        case MODE_VGS:solveVgs(gid,particle,correction,voxel,control,u);break;
        case MODE_STATIC_COLLISIONS:staticCollisions(gid,particle,collisionId,collisionControl,staticCell,staticCollider,refcount,control,u);break;
        case MODE_BREAK_MASK:gatherBreakMask(gid,particle,voxel,topology,refcount,control,u);break;
        case MODE_FINALIZE_PARTICLES:finalizeParticle(gid,particle,cell,simId,refcount,control,u);break;
        case MODE_FINALIZE_VOXELS:finalizeVoxel(gid,particle,voxel,control,u);break;
        case MODE_SPLIT_BREAKS:splitBrokenFaces(gid,particle,correction,cell,simId,collisionId,collisionMember,collisionControl,voxel,tetherOwner,collisionMeta,refcount,control,cloneParent,u);break;
        case MODE_PREPARE_INDIRECT:prepareIndirect(gid,control,collisionControl,dispatchArgs);break;
        case MODE_WAKE_GATHER:wakeGather(gid,voxel,topology,u);break;
        case MODE_WAKE_APPLY:wakeApply(gid,voxel,u);break;
        case MODE_TOPOLOGY_REBUILD_SERIAL:break;
        case MODE_PBF_BUILD_NEIGHBORS:pbfBuildNeighbors(gid,particle,cell,simId,collisionId,hashHead,hashNext,fluidNeighbors,control,u);break;
        case MODE_PBF_LAMBDA:pbfLambda(gid,particle,cell,simId,collisionId,hashHead,hashNext,fluid,fluidNeighbors,control,u);break;
        case MODE_PBF_DELTA:pbfDelta(gid,particle,cell,simId,collisionId,hashHead,hashNext,fluid,fluidNeighbors,staticCell,staticCollider,collisionControl,control,u);break;
        case MODE_PBF_APPLY:pbfApplyAndStaticCollisions(gid,particle,fluid,staticCell,staticCollider,collisionControl,simId,control,u);break;
        case MODE_PBF_STATIC_COLLISIONS:pbfStaticCollisions(gid,particle,staticCell,staticCollider,collisionControl,simId,control,u);break;
        case MODE_PBF_DYNAMIC_SOLID_COLLISIONS:pbfDynamicSolidCollisions(gid,true,particle,correction,voxel,simId,control,u);break;
        case MODE_PBF_DYNAMIC_SOLID_FINAL:pbfDynamicSolidCollisions(gid,false,particle,correction,voxel,simId,control,u);break;
        case MODE_PBF_VISCOSITY:pbfViscosity(gid,particle,cell,simId,collisionId,hashHead,hashNext,fluid,fluidNeighbors,control,u);break;
        case MODE_PBF_APPLY_VISCOSITY:pbfApplyViscosity(gid,particle,fluid,simId,control,u);break;
    }
}
