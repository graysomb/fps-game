/* Appended to pbd_pipeline.metal and mesh.metal at runtime. The VGS numerical
 * projection below uses the existing solver equations; output is private. */
struct AParticleInput { int xyz[3];uint arena;uint incident[8];float inverse_mass,base_inverse_mass;uint count,simulated; };
struct AControl {uint source[8];float weight[8];uint independent,attachment,begin,end;};
struct AShape {uint corner[8];uint source,level,active,reserved;int origin[3];uint padding;};
struct APhysicsParams {uint mode,count,n,p,extent,bit,src,dst,aux;};
static_assert(sizeof(AParticleInput)==64,"particle input ABI");
static_assert(sizeof(AControl)==80,"control ABI");
static_assert(sizeof(AShape)==64,"shape ABI");
inline float a_weight(int3 point,int3 origin,uint level,uint corner){float3 uv=float3(point-origin)/float(1<<level);return ((corner&1)?uv.x:1-uv.x)*((corner&2)?uv.y:1-uv.y)*((corner&4)?uv.z:1-uv.z);}
kernel void adaptive_physics(
    device ParticleState *particle [[buffer(0)]],device VoxelState *voxel [[buffer(1)]],
    device const AParticleInput *info [[buffer(2)]],device AControl *controls [[buffer(3)]],
    device const ALeaf *leaves [[buffer(4)]],device const ALeaf *fine [[buffer(5)]],
    device uint *h [[buffer(6)]],device AShape *shapes [[buffer(7)]],
    device uint2 *refs [[buffer(8)]],device uint2 *other [[buffer(9)]],
    device float4 *delta [[buffer(10)]],device uint *work [[buffer(11)]],
    device uint *owner [[buffer(12)]],constant APhysicsParams &a [[buffer(13)]],
    constant GpuUniforms &u [[buffer(14)]],device uint *sim [[buffer(15)]],
    device int *pcontrol [[buffer(16)]],device const uint *arena_ids [[buffer(17)]],
    device const uint2 *arenas [[buffer(18)]],device float4 *momentum [[buffer(19)]],device const uint *slot_order [[buffer(20)]],
    device const uint4 *ref_ranges [[buffer(21)]],device const uint *fine_arena [[buffer(22)]],uint gid [[thread_position_in_grid]],
    uint tid [[thread_index_in_threadgroup]],uint group [[threadgroup_position_in_grid]]) {
    threadgroup uint scan[128];
    if(h[2])return;
    if((a.mode<30||(a.mode>=40&&a.mode<=42))&&!h[9])return;
    if(a.mode==3){uint v=gid<a.count?work[a.src+gid]:0;scan[tid]=v;threadgroup_barrier(mem_flags::mem_threadgroup);
        for(uint step=1;step<128;step*=2){uint t=tid>=step?scan[tid-step]:0;threadgroup_barrier(mem_flags::mem_threadgroup);scan[tid]+=t;threadgroup_barrier(mem_flags::mem_threadgroup);}
        if(gid<a.count)work[a.dst+gid]=scan[tid]-v;if(tid==127)work[a.aux+group]=scan[tid];return;}
    if(gid>=a.count)return;
    uint leaf_count=h[0],shape_slots=8*a.n;
    switch(a.mode){
    case 45: {
        uint broken=0;
        if(gid<leaf_count&&a_needs_refine(leaves[gid],fine,a.n,reinterpret_cast<device const uint *>(voxel)))broken=1;
        work[gid]=broken;break;
    }
    case 47: if(!gid)h[9]=1;break;
    case 46: if(!gid){uint total=work[a.dst+a.n-1]+work[a.n-1];h[16]=0;h[17]=total?a.bit:0;h[10]+=h[17];h[9]=uint(total!=0);}break;
    case 42: if(gid<a.p&&momentum[8*info[gid].arena+7].x>0&&particle[gid].prev_inv_mass.w>0){uint arena=info[gid].arena;float3 dv=momentum[arena*8+4].xyz,omega=momentum[arena*8+5].xyz,center=momentum[arena*8+6].xyz;particle[gid].velocity.xyz+=dv+cross(omega,particle[gid].pos_radius.xyz-center);}break;
    case 0: if(gid<a.n&&momentum[8*fine_arena[fine[gid].source]+7].x>0){uint leaf=alookup(leaves,leaf_count,fine[gid].arena,fine[gid].key);owner[fine[gid].source]=leaves[leaf].source;}break;
    case 1: if(gid<leaf_count&&momentum[8*fine_arena[leaves[gid].source]+7].x>0){ALeaf c=leaves[gid];int3 p=aorigin(c.key);AShape shape;shape.level=c.level;shape.source=fine[alookup(fine,a.n,c.arena,c.key)].source;shape.active=1;shape.reserved=0;shape.origin[0]=p.x;shape.origin[1]=p.y;shape.origin[2]=p.z;shape.padding=0;
        for(uint k=0;k<8;k++){int3 q=p+int3(k&1,(k>>1)&1,(k>>2)&1)*((1<<c.level)-1);uint i=alookup(fine,a.n,c.arena,akey(q));shape.corner[k]=voxelParticle(voxel[fine[i].source],k);}
        shapes[c.source]=shape;}break;
    case 2: if(gid<a.p&&momentum[8*info[gid].arena+7].x>0){AParticleInput pi=info[gid];int3 p(pi.xyz[0],pi.xyz[1],pi.xyz[2]);AControl c;uint coarse=0xffffffffu;uint highest=0;bool independent=pi.count==0;
        for(uint j=0;j<pi.count;j++){uint leaf=owner[pi.incident[j]/8];AShape l=shapes[leaf];bool used=false;
            for(uint k=0;k<8;k++)used=used||shapes[leaf].corner[k]==gid;
            if(used)independent=true;
            else if(coarse==0xffffffffu||l.level>highest){coarse=leaf;highest=l.level;}}
        c.independent=uint(independent);c.attachment=uint(independent&&coarse!=0xffffffffu);c.begin=c.end=0;
        for(uint k=0;k<8;k++){c.source[k]=gid;c.weight[k]=k==0?1.0f:0.0f;}
        if(coarse!=0xffffffffu)for(uint k=0;k<8;k++){c.source[k]=shapes[coarse].corner[k];c.weight[k]=a_weight(p,int3(shapes[coarse].origin[0],shapes[coarse].origin[1],shapes[coarse].origin[2]),shapes[coarse].level,k);}
        controls[gid]=c;}break;
    case 4: work[a.dst+gid]+=work[a.src+gid/128];break;
    case 5: {uint slot=slot_order[gid],arena=slot<shape_slots?fine_arena[slot/8]:info[(slot-shape_slots)/9].arena;
        if(momentum[8*arena+7].x<=0)break;uint dest=0xffffffffu;
        if(slot<shape_slots){uint shape=slot/8;if(shapes[shape].active)dest=shapes[shape].corner[slot%8];}
        else {uint p=(slot-shape_slots)/9,k=(slot-shape_slots)%9;AControl c=controls[p];if(!k){if(c.independent)dest=p;}else if(c.attachment||!c.independent){if(c.weight[k-1]>0)dest=c.source[k-1];}}
        refs[gid]=uint2(dest,slot);break;}
    case 6: {uint slot=slot_order[gid],arena=slot<shape_slots?fine_arena[slot/8]:info[(slot-shape_slots)/9].arena;
        work[gid]=momentum[8*arena+7].x>0?uint(!((refs[gid].x>>a.bit)&1u)):0;break;}
    case 7: {uint slot=slot_order[gid],arena=slot<shape_slots?fine_arena[slot/8]:info[(slot-shape_slots)/9].arena;
        if(momentum[8*arena+7].x<=0)break;uint2 range=ref_ranges[arena].xy;uint base=work[a.dst+range.x];
        uint zeros=work[a.dst+range.x+range.y-1]+work[range.x+range.y-1]-base,rank=work[a.dst+gid]-base;
        other[range.x+(((refs[gid].x>>a.bit)&1u)?zeros+gid-range.x-rank:rank)]=refs[gid];break;}
    case 8: if(gid<a.p&&momentum[8*info[gid].arena+7].x>0){uint2 range=ref_ranges[info[gid].arena].xy;uint lo=range.x,hi=range.x+range.y;while(lo<hi){uint m=lo+(hi-lo)/2;if(refs[m].x<gid)lo=m+1;else hi=m;}controls[gid].begin=lo;
        hi=range.x+range.y;while(lo<hi){uint m=lo+(hi-lo)/2;if(refs[m].x<=gid)lo=m+1;else hi=m;}controls[gid].end=lo;}break;
    case 9: if(gid<a.p&&momentum[8*info[gid].arena+7].x>0){AControl c=controls[gid];float mass=0;
        for(uint i=c.begin;i<c.end;i++){uint slot=refs[i].y;if(slot<shape_slots)continue;uint p=(slot-shape_slots)/9,k=(slot-shape_slots)%9;float inv=info[p].inverse_mass;if(inv<=0)continue;
            if(!k&&controls[p].independent)mass+=1/inv;
            else if(k&&!controls[p].independent)mass+=controls[p].weight[k-1]/inv;}
        particle[gid].prev_inv_mass.w=c.independent&&info[gid].inverse_mass>0&&mass>0?1/mass:0;
        particle[gid].predicted_base_inv_mass.w=particle[gid].prev_inv_mass.w;
        work[gid]=uint(c.independent&&info[gid].simulated); }break;
    case 10: if(gid<a.p&&work[gid])sim[work[a.dst+gid]]=gid;
        if(!gid)pcontrol[1]=int(work[a.dst+a.p-1]+work[a.p-1]);break;
    case 12: {uint slot=slot_order[gid],arena=slot<shape_slots?fine_arena[slot/8]:info[(slot-shape_slots)/9].arena;if(momentum[8*arena+7].x>0)other[gid]=refs[gid];break;}
    case 13: if(gid<a.p)work[gid]=uint(controls[gid].independent&&info[gid].simulated);break;
    case 48: {uint raw_arena=ref_ranges[gid].z,lo=0,hi=leaf_count;while(lo<hi){uint m=lo+(hi-lo)/2;if(leaves[m].arena<raw_arena)lo=m+1;else hi=m;}uint first=lo;
        hi=leaf_count;while(lo<hi){uint m=lo+(hi-lo)/2;if(leaves[m].arena<=raw_arena)lo=m+1;else hi=m;}
        uint count=lo>first?work[a.dst+lo-1]+work[lo-1]-work[a.dst+first]:0;momentum[8*gid+7].x=count?1.0f:0.0f;break;}
    case 49: momentum[8*gid+7].x=1;break;
    case 50: if(gid<a.n&&momentum[8*fine_arena[gid]+7].x>0)shapes[gid].active=0;break;
    case 11: if(!gid)h[9]=0;break;
    case 30: if(gid<leaf_count){AShape sh=shapes[leaves[gid].source];VoxelState v=voxel[sh.source];float3 p[8],original[8];float weight[8];bool dynamic=false;
        for(uint k=0;k<8;k++){uint id=sh.corner[k];p[k]=original[k]=particle[id].predicted_base_inv_mass.xyz;weight[k]=particle[id].prev_inv_mass.w;dynamic=dynamic||weight[k]>0;delta[8*sh.source+k]=float4(0);}
        if(!dynamic||v.flags.x==0||v.flags.y!=0||v.flags.z!=0||v.flags.w==0)break;
        float edge=v.pos_rest_edge.w*float(1<<sh.level),volume=edge*edge*edge;
        for(int it=0;it<3;it++){float3 center(0);for(int k=0;k<8;k++)center+=p[k];center*=0.125f;
            float3 v0=((p[1]-p[0])+(p[3]-p[2])+(p[5]-p[4])+(p[7]-p[6]))*0.25f;
            float3 v1=((p[2]-p[0])+(p[3]-p[1])+(p[6]-p[4])+(p[7]-p[5]))*0.25f;
            float3 v2=((p[4]-p[0])+(p[5]-p[1])+(p[6]-p[2])+(p[7]-p[3]))*0.25f;
            float3 u0=v0-u.vgs_alpha*(projectOnto(v1,v0,u)+projectOnto(v2,v0,u));
            float3 u1=v1-u.vgs_alpha*(projectOnto(v2,v1,u)+projectOnto(v0,v1,u));
            float3 u2=v2-u.vgs_alpha*(projectOnto(v0,v2,u)+projectOnto(v1,v2,u));
            if(length(u0)>u.vgs_epsilon)u0*=mix(edge,length(v0),u.vgs_beta)/length(u0);
            if(length(u1)>u.vgs_epsilon)u1*=mix(edge,length(v1),u.vgs_beta)/length(u1);
            if(length(u2)>u.vgs_epsilon)u2*=mix(edge,length(v2),u.vgs_beta)/length(u2);
            float vol=dot(cross(u0,u1),u2);if(fabs(vol)>u.vgs_epsilon){float scale=volume/vol,root=pow(fabs(scale),1.0f/3.0f)*(scale<0?-1.0f:1.0f);u0*=0.5f*root;u1*=0.5f*root;u2*=0.5f*root;}
            for(uint k=0;k<8;k++)p[k]=center+((k&1)?u0:-u0)+((k&2)?u1:-u1)+((k&4)?u2:-u2);
        }
        for(uint k=0;k<8;k++)if(weight[k]>0)delta[8*sh.source+k]=float4((p[k]-original[k])*weight[k],weight[k]);
        }break;
    case 31: if(gid<a.p){AControl c=controls[gid];uint base=shape_slots+9*gid;for(uint k=0;k<9;k++)delta[base+k]=float4(0);
        if(!c.attachment)break;float wp=particle[gid].prev_inv_mass.w,denom=wp;float3 target(0);
        for(uint k=0;k<8;k++){float w=c.weight[k];target+=w*particle[c.source[k]].predicted_base_inv_mass.xyz;denom+=w*w*particle[c.source[k]].prev_inv_mass.w;}
        if(denom<=1e-12f)break;float3 error=particle[gid].predicted_base_inv_mass.xyz-target;
        if(wp>0)delta[base]=float4(-error*wp/denom,1);
        for(uint k=0;k<8;k++){float factor=c.weight[k]*particle[c.source[k]].prev_inv_mass.w/denom;if(factor>0)delta[base+1+k]=float4(error*factor,1);}
        }break;
    case 34: case 35: if(gid<a.p&&controls[gid].independent){AControl c=controls[gid];float4 sum(0);
        for(uint i=c.begin;i<c.end;i++){uint slot=refs[i].y;if((a.mode==34)==(slot<shape_slots))sum+=delta[slot];}
        if(sum.w>0&&particle[gid].prev_inv_mass.w>0)particle[gid].predicted_base_inv_mass.xyz+=sum.xyz/sum.w*u.sor;
        }break;
    case 32: if(gid<a.p&&!controls[gid].independent){AControl c=controls[gid];ParticleState p=particle[gid];float3 pos(0),prev(0),pred(0),vel(0);
        for(uint k=0;k<8;k++){ParticleState s=particle[c.source[k]];float w=c.weight[k];pos+=w*s.pos_radius.xyz;prev+=w*s.prev_inv_mass.xyz;pred+=w*s.predicted_base_inv_mass.xyz;vel+=w*s.velocity.xyz;}
        p.pos_radius.xyz=pos;p.prev_inv_mass.xyz=prev;p.predicted_base_inv_mass.xyz=pred;p.velocity.xyz=vel;particle[gid]=p;}break;
    }
}

/* Topology-only reduction gets its own pipeline: reserving 8 KiB of shared
 * memory in every VGS/gather workgroup would reduce steady-state occupancy. */
kernel void adaptive_momentum(device ParticleState *particle [[buffer(0)]],
    device uint *h [[buffer(6)]],constant APhysicsParams &a [[buffer(13)]],
    device const uint *arena_ids [[buffer(17)]],device const uint2 *arenas [[buffer(18)]],
    device float4 *momentum [[buffer(19)]],uint tid [[thread_index_in_threadgroup]],
    uint group [[threadgroup_position_in_grid]]) {
    threadgroup float sums[16][128];
    if(h[2]||!h[9]||momentum[8*group+7].x<=0)return;
    if(a.mode==40||a.mode==41){
        float values[16]={0};uint2 range=arenas[group];
        for(uint j=tid;j<range.y;j+=128){uint id=arena_ids[range.x+j];ParticleState p=particle[id];
            float inv=p.prev_inv_mass.w;if(inv<=0)continue;float mass=1/inv;float3 x=p.pos_radius.xyz,v=p.velocity.xyz,linear=mass*v,angular=cross(x,linear);
            values[0]+=linear.x;values[1]+=linear.y;values[2]+=linear.z;values[3]+=angular.x;values[4]+=angular.y;values[5]+=angular.z;
            values[6]+=mass;values[7]+=mass*x.x;values[8]+=mass*x.y;values[9]+=mass*x.z;
            values[10]+=mass*(x.y*x.y+x.z*x.z);values[11]+=mass*(x.x*x.x+x.z*x.z);values[12]+=mass*(x.x*x.x+x.y*x.y);
            values[13]-=mass*x.x*x.y;values[14]-=mass*x.x*x.z;values[15]-=mass*x.y*x.z;
        }
        for(uint c=0;c<16;c++)sums[c][tid]=values[c];threadgroup_barrier(mem_flags::mem_threadgroup);
        for(uint step=64;step;step/=2){if(tid<step)for(uint c=0;c<16;c++)sums[c][tid]+=sums[c][tid+step];threadgroup_barrier(mem_flags::mem_threadgroup);}
        if(!tid){float3 linear(sums[0][0],sums[1][0],sums[2][0]),angular(sums[3][0],sums[4][0],sums[5][0]);float mass=sums[6][0];
            if(a.mode==40){momentum[group*8]=float4(linear,mass);momentum[group*8+1]=float4(angular,0);}
            else {float3 center=mass>0?float3(sums[7][0],sums[8][0],sums[9][0])/mass:float3(0);
                float3 dp=momentum[group*8].xyz-linear,dv=mass>0?dp/mass:float3(0);
                float3 dl=momentum[group*8+1].xyz-angular-cross(center,dp);
                float xx=sums[10][0]-mass*(center.y*center.y+center.z*center.z),yy=sums[11][0]-mass*(center.x*center.x+center.z*center.z),zz=sums[12][0]-mass*(center.x*center.x+center.y*center.y);
                float xy=sums[13][0]+mass*center.x*center.y,xz=sums[14][0]+mass*center.x*center.z,yz=sums[15][0]+mass*center.y*center.z;
                float3 r0(xx,xy,xz),r1(xy,yy,yz),r2(xz,yz,zz);float determinant=dot(r0,cross(r1,r2));
                float3 omega=fabs(determinant)>1e-12f?(dl.x*cross(r1,r2)+dl.y*cross(r2,r0)+dl.z*cross(r0,r1))/determinant:float3(0);
                momentum[group*8+4]=float4(dv,0);momentum[group*8+5]=float4(omega,0);momentum[group*8+6]=float4(center,0);
            }}return;
    }
}
