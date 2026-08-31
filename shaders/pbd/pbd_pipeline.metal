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

struct GpuUniforms {
    int particle_count, sim_count, voxel_count, static_collider_count;
    int hash_size, static_hash_size, active_players, break_damp_frames;
    int mode;
    int integer_padding_0, integer_padding_1, integer_padding_2;
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
constant int CONTROL_FLAG_OVERFLOW = 1;
constant int CONTROL_FLAG_TOPOLOGY_DIRTY = 2;

constant int FACE_CORNERS[24] = {
    1,3,5,7, 0,2,4,6, 2,3,6,7,
    0,1,4,5, 4,5,6,7, 0,1,2,3
};

inline int controlLoad(device atomic_int *control, int component) {
    return atomic_load_explicit(&control[component], memory_order_relaxed);
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
    uint expected = atomic_load_explicit(value, memory_order_relaxed);
    for (;;) {
        uint desired = as_type<uint>(as_type<float>(expected) + addend);
        if (atomic_compare_exchange_weak_explicit(value, &expected, desired,
                                                  memory_order_relaxed,
                                                  memory_order_relaxed)) return;
    }
}

inline void accumulate(device atomic_uint *correction, device atomic_int *control,
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
                             device uint *simId, device atomic_int *refcount,
                             device atomic_int *control) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (atomic_load_explicit(&refcount[id], memory_order_relaxed) <= 0) return;
    for (int c = 0; c < 4; ++c)
        atomic_store_explicit(&correction[id * 4u + uint(c)], 0u, memory_order_relaxed);
}

inline void integrateParticle(uint gid, device ParticleState *particle,
                              device int *tetherOwner, device uint *simId,
                              device atomic_int *refcount, device atomic_int *control,
                              constant GpuUniforms &u) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (atomic_load_explicit(&refcount[id], memory_order_relaxed) <= 0) return;
    ParticleState p = particle[id];
    p.prev_inv_mass.xyz = p.pos_radius.xyz;
    p.predicted_base_inv_mass.xyz = p.pos_radius.xyz;
    if (p.prev_inv_mass.w > 0.0f) {
        p.velocity.xyz *= u.velocity_damping;
        p.predicted_base_inv_mass.xyz += p.velocity.xyz * u.dt;
        p.predicted_base_inv_mass.y -= u.gravity * u.dt * u.dt;
        int player = tetherOwner[id];
        if (player >= 0 && player < 4) {
            float3 accel = (u.tether_targets[player].xyz - p.predicted_base_inv_mass.xyz) * u.tether_spring;
            p.predicted_base_inv_mass.xyz += accel * u.dt * u.dt;
            p.velocity.xyz += accel * u.dt;
            p.velocity.xyz *= 1.0f - clamp(u.tether_damping * u.dt, 0.0f, 0.9f);
        }
    }
    particle[id] = p;
}

inline void clearHash(uint gid, device atomic_int *hashHead, constant GpuUniforms &u) {
    if (gid < uint(u.hash_size)) atomic_store_explicit(&hashHead[gid], -1, memory_order_relaxed);
}

inline void buildHash(uint gid, device ParticleState *particle, device int4 *cell,
                      device uint *simId, device atomic_int *hashHead,
                      device int *hashNext, device atomic_int *refcount,
                      device atomic_int *control, constant GpuUniforms &u) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (atomic_load_explicit(&refcount[id], memory_order_relaxed) <= 0) return;
    int3 c = int3(floor(particle[id].predicted_base_inv_mass.xyz / u.voxel_size));
    cell[id].xyz = c;
    uint h = hashCoord(c, u.hash_size);
    hashNext[gid] = atomic_exchange_explicit(&hashHead[h], int(gid), memory_order_relaxed);
}

inline void pairCollisions(uint gid, device ParticleState *particle,
                           device atomic_uint *correction, device int4 *cell,
                           device uint *simId, device atomic_int *hashHead,
                           device int *hashNext, device int *collisionMeta,
                           device atomic_int *refcount, device atomic_int *control,
                           constant GpuUniforms &u) {
    int simCount = controlLoad(control, 1);
    if (gid >= uint(simCount)) return;
    uint aid = simId[gid];
    if (atomic_load_explicit(&refcount[aid], memory_order_relaxed) <= 0) return;
    ParticleState a = particle[aid];
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
        while (item >= 0 && item < simCount) {
            bool sameCell = dx == 0 && dy == 0 && dz == 0;
            if (!sameCell || uint(item) > gid) {
                uint bid = simId[item];
                if (atomic_load_explicit(&refcount[bid], memory_order_relaxed) <= 0) {
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
    }
}

inline void applyCorrections(uint gid, device ParticleState *particle,
                             device atomic_uint *correction, device uint *simId,
                             device atomic_int *refcount, device atomic_int *control,
                             constant GpuUniforms &u) {
    if (gid >= uint(controlLoad(control, 1))) return;
    uint id = simId[gid];
    if (atomic_load_explicit(&refcount[id], memory_order_relaxed) <= 0) return;
    device atomic_uint *base = correction + id * 4u;
    float weight = as_type<float>(atomic_load_explicit(base + 3, memory_order_relaxed));
    if (weight > 0.0f) {
        float3 sum(as_type<float>(atomic_load_explicit(base + 0, memory_order_relaxed)),
                   as_type<float>(atomic_load_explicit(base + 1, memory_order_relaxed)),
                   as_type<float>(atomic_load_explicit(base + 2, memory_order_relaxed)));
        particle[id].predicted_base_inv_mass.xyz += sum * (u.sor / weight);
    }
    for (int c = 0; c < 4; ++c) atomic_store_explicit(base + c, 0u, memory_order_relaxed);
}

inline void solveVgs(uint gid, device ParticleState *particle,
                     device atomic_uint *correction, device VoxelState *voxel,
                     device atomic_int *control, constant GpuUniforms &u) {
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
    for(int i=0;i<8;++i)if(applyWeight[i]>0.0f)accumulate(correction,control,voxelParticle(v,i),p[i]-original[i],applyWeight[i]);
}

inline int findStaticCell(int3 coord, device int4 *staticCell, constant GpuUniforms &u) {
    uint h=hashCoord(coord,u.static_hash_size);
    for(int probe=0;probe<u.static_hash_size;++probe){int4 entry=staticCell[h];if(entry.w<0)return -1;if(all(entry.xyz==coord))return entry.w;h=(h+1u)&uint(u.static_hash_size-1);}return -1;
}

inline float3 pushOutOfBox(float3 pos,float radius,StaticCollider box,constant GpuUniforms &u){
    float3 bmin=box.bounds_min.xyz,bmax=box.bounds_max.xyz;bool inside=all(pos>=bmin)&&all(pos<=bmax);constexpr float eps=1e-6f;
    if(!inside){float3 closest=clamp(pos,bmin,bmax),delta=pos-closest;float distSq=dot(delta,delta);if(distSq>=radius*radius)return pos;float dist=sqrt(max(distSq,eps));float3 normal=dist>eps?delta/dist:float3(0,1,0);return pos+normal*max(radius-dist,0.0f);}
    float topGap=bmax.y-pos.y,bottomGap=pos.y-bmin.y,verticalBias=u.voxel_size*0.25f;bool preferTop=pos.y>=box.center.y&&topGap<=radius+verticalBias;bool preferBottom=pos.y<box.center.y&&bottomGap<=radius+verticalBias;
    if(preferTop||preferBottom){bool top=preferTop&&(!preferBottom||topGap<=bottomGap);float3 normal=top?float3(0,1,0):float3(0,-1,0);float gap=top?topGap:bottomGap;float penetration=radius-gap;if(penetration<0.0f)penetration=radius;return pos+normal*penetration;}
    float3 low=pos-bmin,high=bmax-pos;float distances[6]={low.x,high.x,low.y,high.y,low.z,high.z};int best=0;for(int i=1;i<6;++i)if(distances[i]<distances[best])best=i;float3 normal(0.0f);if(best==0)normal.x=-1;else if(best==1)normal.x=1;else if(best==2)normal.y=-1;else if(best==3)normal.y=1;else if(best==4)normal.z=-1;else normal.z=1;return pos+normal*(radius+distances[best]);
}

inline void staticCollisions(uint gid,device ParticleState *particle,device uint *simId,device int4 *staticCell,device StaticCollider *staticCollider,device atomic_int *refcount,device atomic_int *control,constant GpuUniforms &u){
    if(gid>=uint(controlLoad(control,1)))return;uint id=simId[gid];if(atomic_load_explicit(&refcount[id],memory_order_relaxed)<=0)return;ParticleState p=particle[id];if(p.prev_inv_mass.w<=0.0f)return;float radius=p.pos_radius.w;float3 pos=p.predicted_base_inv_mass.xyz;float terrainLimit=u.floor_size-radius,floorLimit=max(0.0f,0.5f*u.voxel_size-radius);pos.y=max(pos.y,floorLimit);pos.xz=clamp(pos.xz,float2(-terrainLimit),float2(terrainLimit));constexpr float eps=1e-6f;
    for(int i=0;i<u.active_players&&i<4;++i){if(u.players[i].w<0.0f)continue;float halfSize=u.players[i].w;float3 nearest=clamp(pos,u.players[i].xyz-float3(halfSize),u.players[i].xyz+float3(halfSize));float3 delta=pos-nearest;float distSq=dot(delta,delta);if(distSq<radius*radius){float dist=sqrt(max(distSq,eps));float3 normal=dist>eps?delta/dist:float3(0,1,0);pos+=normal*(radius-dist);}}
    int3 center=int3(floor(pos/u.voxel_size));int seen[27];int seenCount=0;for(int z=-1;z<=1;++z)for(int y=-1;y<=1;++y)for(int x=-1;x<=1;++x){int voxelId=findStaticCell(center+int3(x,y,z),staticCell,u);if(voxelId<0||voxelId>=u.static_collider_count)continue;bool duplicate=false;for(int s=0;s<seenCount;++s)duplicate=duplicate||seen[s]==voxelId;if(duplicate)continue;seen[seenCount++]=voxelId;pos=pushOutOfBox(pos,0.25f*u.voxel_size,staticCollider[voxelId],u);}particle[id].predicted_base_inv_mass.xyz=pos;
}

inline void gatherBreakMask(uint gid,device ParticleState *particle,device VoxelState *voxel,device atomic_int *control,constant GpuUniforms &u){
    if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];v.lifecycle.w=0u;if(v.flags.x==0||v.flags.y!=0||v.flags.z!=0||v.pos_rest_edge.w<=0.0f){voxel[gid]=v;return;}float3 p[8];for(int i=0;i<8;++i){uint id=voxelParticle(v,i);if(id>=uint(controlLoad(control,0)))return;p[i]=particle[id].predicted_base_inv_mass.xyz;}float3 axis0=((p[1]-p[0])+(p[3]-p[2])+(p[5]-p[4])+(p[7]-p[6]))*0.25f;float3 axis1=((p[2]-p[0])+(p[3]-p[1])+(p[6]-p[4])+(p[7]-p[5]))*0.25f;float3 axis2=((p[4]-p[0])+(p[5]-p[1])+(p[6]-p[2])+(p[7]-p[3]))*0.25f;float3 lengths(length(axis0),length(axis1),length(axis2));float3 strain=fabs(lengths-float3(v.pos_rest_edge.w))/v.pos_rest_edge.w;uint glued=v.lifecycle.z,mask=0u;bool exceeded=false;if(strain.x>u.strain_threshold){exceeded=true;mask|=glued&3u;}if(strain.y>u.strain_threshold){exceeded=true;mask|=glued&12u;}if(strain.z>u.strain_threshold){exceeded=true;mask|=glued&48u;}float3 n0=lengths.x>u.vgs_epsilon?axis0/lengths.x:float3(0);float3 n1=lengths.y>u.vgs_epsilon?axis1/lengths.y:float3(0);float3 n2=lengths.z>u.vgs_epsilon?axis2/lengths.z:float3(0);if(fabs(dot(n0,n1))>u.shear_threshold){exceeded=true;mask|=glued&48u;}if(fabs(dot(n0,n2))>u.shear_threshold){exceeded=true;mask|=glued&12u;}if(fabs(dot(n1,n2))>u.shear_threshold){exceeded=true;mask|=glued&3u;}v.lifecycle.w=mask;if(exceeded)v.lifecycle.y=1u;voxel[gid]=v;
}

inline uint allocateClone(uint parent,device ParticleState *particle,device atomic_uint *correction,device int4 *cell,device uint *simId,device int *tetherOwner,device int *collisionMeta,device atomic_int *refcount,device atomic_int *control,device int *cloneParent,constant GpuUniforms &u){
    int slot=-1;for(;;){int current=controlLoad(control,0);if(current>=controlLoad(control,2)){atomic_store_explicit(&control[3],CONTROL_FLAG_OVERFLOW,memory_order_relaxed);return uint(-1);}int expected=current;if(atomic_compare_exchange_weak_explicit(&control[0],&expected,current+1,memory_order_relaxed,memory_order_relaxed)){slot=current;break;}}uint id=uint(slot);particle[id]=particle[parent];cell[id]=cell[parent];cell[id].w=u.break_damp_frames;for(int c=0;c<4;++c)collisionMeta[id*4u+uint(c)]=collisionMeta[parent*4u+uint(c)];tetherOwner[id]=tetherOwner[parent];for(int c=0;c<4;++c)atomic_store_explicit(&correction[id*4u+uint(c)],0u,memory_order_relaxed);atomic_store_explicit(&refcount[id],1,memory_order_relaxed);cloneParent[id]=int(parent);if(particle[id].predicted_base_inv_mass.w>0.0f){int simSlot=atomic_fetch_add_explicit(&control[1],1,memory_order_relaxed);if(simSlot<controlLoad(control,2))simId[simSlot]=id;else atomic_store_explicit(&control[3],CONTROL_FLAG_OVERFLOW,memory_order_relaxed);}return id;
}

inline void setVoxelParticle(thread VoxelState &v,int corner,uint id){if(corner<4)v.particle_0_3[corner]=id;else v.particle_4_7[corner-4]=id;}

inline void detachCorner(thread VoxelState &v,int corner,device ParticleState *particle,device atomic_uint *correction,device int4 *cell,device uint *simId,device int *tetherOwner,device int *collisionMeta,device atomic_int *refcount,device atomic_int *control,device int *cloneParent,constant GpuUniforms &u){
    uint oldId=voxelParticle(v,corner);if(oldId>=uint(controlLoad(control,0)))return;for(;;){int refs=atomic_load_explicit(&refcount[oldId],memory_order_relaxed);if(refs<=1)return;int expected=refs;if(!atomic_compare_exchange_weak_explicit(&refcount[oldId],&expected,refs-1,memory_order_relaxed,memory_order_relaxed))continue;uint newId=allocateClone(oldId,particle,correction,cell,simId,tetherOwner,collisionMeta,refcount,control,cloneParent,u);if(newId==uint(-1)){atomic_fetch_add_explicit(&refcount[oldId],1,memory_order_relaxed);return;}cell[oldId].w=u.break_damp_frames;setVoxelParticle(v,corner,newId);return;}
}

inline void splitBrokenFaces(uint gid,device ParticleState *particle,device atomic_uint *correction,device int4 *cell,device uint *simId,device VoxelState *voxel,device int *tetherOwner,device int *collisionMeta,device atomic_int *refcount,device atomic_int *control,device int *cloneParent,constant GpuUniforms &u){
    if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];uint mask=v.lifecycle.w&63u;if(mask==0u)return;atomic_fetch_or_explicit(&control[3],CONTROL_FLAG_TOPOLOGY_DIRTY,memory_order_relaxed);for(int face=0;face<6;++face){if((mask&(1u<<face))==0u)continue;for(int c=0;c<4;++c)detachCorner(v,FACE_CORNERS[face*4+c],particle,correction,cell,simId,tetherOwner,collisionMeta,refcount,control,cloneParent,u);}v.lifecycle.z&=~mask;v.lifecycle.w=0u;voxel[gid]=v;
}

inline int topologyNeighbor(device int4 *topology,int voxelId,int face){return face<4?topology[voxelId*2][face]:topology[voxelId*2+1][face-4];}
inline int topologyRoot(device int4 *topology,int voxelId,constant GpuUniforms &u){int root=voxelId;for(int i=0;i<32;++i){int parent=topology[root*2+1].z;if(parent==root||parent<0||parent>=u.voxel_count)break;root=parent;}return root;}
inline void topologyUnion(device VoxelState *voxel,device int4 *topology,int gid,constant GpuUniforms &u){VoxelState v=voxel[gid];if((v.lifecycle.z&0x80000000u)==0u)return;for(int face=0;face<6;++face){if((v.lifecycle.z&(1u<<face))==0u)continue;int neighbor=topologyNeighbor(topology,gid,face);if(neighbor<0||neighbor>=u.voxel_count||(voxel[neighbor].lifecycle.z&0x80000000u)==0u)continue;int a=topologyRoot(topology,gid,u),b=topologyRoot(topology,neighbor,u);if(a!=b)topology[max(a,b)*2+1].z=min(a,b);}}
inline void assignGroup(device int *collisionMeta,uint id,int group){uint offset=id*4u+3u;int old=collisionMeta[offset];if(old==group||old==-2)return;collisionMeta[offset]=old==-1?group:-2;}
inline void topologyRebuild(uint gid,device ParticleState *particle,device uint *simId,device VoxelState *voxel,device int *collisionMeta,device atomic_int *refcount,device atomic_int *control,device int4 *topology,constant GpuUniforms &u){
    if(gid!=0u||(controlLoad(control,3)&CONTROL_FLAG_TOPOLOGY_DIRTY)==0)return;for(int i=0;i<u.voxel_count;++i)topology[i*2+1].z=i;for(int i=0;i<u.voxel_count;++i)topologyUnion(voxel,topology,i,u);for(int pass=0;pass<18;++pass)for(int i=0;i<u.voxel_count;++i)topology[i*2+1].z=topologyRoot(topology,i,u);int simCount=controlLoad(control,1);for(int i=0;i<simCount;++i){uint id=simId[i];if(atomic_load_explicit(&refcount[id],memory_order_relaxed)>0)collisionMeta[id*4u+3u]=-1;}for(int i=0;i<u.voxel_count;++i){VoxelState v=voxel[i];if((v.lifecycle.z&0x80000000u)==0u||v.flags.y!=0||v.flags.z!=0)continue;int group=topologyRoot(topology,i,u);for(int c=0;c<8;++c){uint id=voxelParticle(v,c);if(id<uint(controlLoad(control,0))&&atomic_load_explicit(&refcount[id],memory_order_relaxed)>0&&particle[id].predicted_base_inv_mass.w>0.0f)assignGroup(collisionMeta,id,group);}}for(int i=0;i<simCount;++i){uint id=simId[i];if(atomic_load_explicit(&refcount[id],memory_order_relaxed)>0&&collisionMeta[id*4u+3u]==-2)collisionMeta[id*4u+3u]=-1;}atomic_fetch_and_explicit(&control[3],~CONTROL_FLAG_TOPOLOGY_DIRTY,memory_order_relaxed);
}

inline void wakeGather(uint gid,device VoxelState *voxel,device int4 *topology,constant GpuUniforms &u){if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];if(v.flags.x==0||v.flags.y!=0)return;bool wake=v.lifecycle.y!=0u;if(!wake)for(int face=0;face<6;++face){int neighbor=topologyNeighbor(topology,int(gid),face);if(neighbor>=0&&neighbor<u.voxel_count&&voxel[neighbor].flags.y==0&&voxel[neighbor].lifecycle.y!=0u){wake=true;break;}}int timer=max(int(v.bounds_min.w)-1,0);if(wake)timer=30;v.bounds_max.w=float(timer);voxel[gid]=v;}
inline void wakeApply(uint gid,device VoxelState *voxel,constant GpuUniforms &u){if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];v.bounds_min.w=v.bounds_max.w;voxel[gid]=v;}
inline void prepareIndirect(uint gid,device atomic_int *control,device uint *dispatchArgs){if(gid!=0u)return;dispatchArgs[0]=(uint(max(controlLoad(control,1),0))+127u)/128u;dispatchArgs[1]=1u;dispatchArgs[2]=1u;}
inline void finalizeParticle(uint gid,device ParticleState *particle,device int4 *cell,device uint *simId,device atomic_int *refcount,device atomic_int *control,constant GpuUniforms &u){if(gid>=uint(controlLoad(control,1)))return;uint id=simId[gid];if(atomic_load_explicit(&refcount[id],memory_order_relaxed)<=0)return;ParticleState p=particle[id];float3 delta=p.predicted_base_inv_mass.xyz-p.prev_inv_mass.xyz;p.velocity.xyz=p.prev_inv_mass.w>0.0f&&u.dt>0.0f?delta/u.dt:float3(0);p.pos_radius.xyz=p.predicted_base_inv_mass.xyz;if(cell[id].w>0)cell[id].w--;particle[id]=p;}
inline void finalizeVoxel(uint gid,device ParticleState *particle,device VoxelState *voxel,device atomic_int *control,constant GpuUniforms &u){if(gid>=uint(u.voxel_count))return;VoxelState v=voxel[gid];if(v.flags.x==0||v.flags.y!=0||v.flags.z!=0)return;float3 center(0),previous(0);for(int i=0;i<8;++i){uint id=voxelParticle(v,i);if(id>=uint(controlLoad(control,0)))return;center+=particle[id].predicted_base_inv_mass.xyz;previous+=particle[id].prev_inv_mass.xyz;}center*=0.125f;previous*=0.125f;if(v.lifecycle.x==0u&&u.dt>0.0f)v.velocity_rest_volume.xyz=(center-previous)/u.dt;else if(v.lifecycle.x>0u)v.lifecycle.x--;v.pos_rest_edge.xyz=center;voxel[gid]=v;}

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
    device atomic_int *refcount [[buffer(11)]],
    device atomic_int *control [[buffer(12)]],
    device int *cloneParent [[buffer(13)]],
    device uint *dispatchArgs [[buffer(14)]],
    device int4 *topology [[buffer(15)]],
    constant GpuUniforms &u [[buffer(16)]],
    uint gid [[thread_position_in_grid]]) {
    switch(u.mode){
        case MODE_RESET:resetCorrections(gid,correction,simId,refcount,control);break;
        case MODE_INTEGRATE:integrateParticle(gid,particle,tetherOwner,simId,refcount,control,u);break;
        case MODE_HASH_CLEAR:clearHash(gid,hashHead,u);break;
        case MODE_HASH_BUILD:buildHash(gid,particle,cell,simId,hashHead,hashNext,refcount,control,u);break;
        case MODE_PAIR_COLLISIONS:pairCollisions(gid,particle,correction,cell,simId,hashHead,hashNext,collisionMeta,refcount,control,u);break;
        case MODE_APPLY:applyCorrections(gid,particle,correction,simId,refcount,control,u);break;
        case MODE_VGS:solveVgs(gid,particle,correction,voxel,control,u);break;
        case MODE_STATIC_COLLISIONS:staticCollisions(gid,particle,simId,staticCell,staticCollider,refcount,control,u);break;
        case MODE_BREAK_MASK:gatherBreakMask(gid,particle,voxel,control,u);break;
        case MODE_FINALIZE_PARTICLES:finalizeParticle(gid,particle,cell,simId,refcount,control,u);break;
        case MODE_FINALIZE_VOXELS:finalizeVoxel(gid,particle,voxel,control,u);break;
        case MODE_SPLIT_BREAKS:splitBrokenFaces(gid,particle,correction,cell,simId,voxel,tetherOwner,collisionMeta,refcount,control,cloneParent,u);break;
        case MODE_PREPARE_INDIRECT:prepareIndirect(gid,control,dispatchArgs);break;
        case MODE_WAKE_GATHER:wakeGather(gid,voxel,topology,u);break;
        case MODE_WAKE_APPLY:wakeApply(gid,voxel,u);break;
        case MODE_TOPOLOGY_REBUILD_SERIAL:topologyRebuild(gid,particle,simId,voxel,collisionMeta,refcount,control,topology,u);break;
    }
}
