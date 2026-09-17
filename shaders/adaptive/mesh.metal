#include <metal_stdlib>
using namespace metal;
struct ALeaf { ulong key; uint arena, material, level, flags, source, reserved; };
struct AMeshParams { uint mode, count, bit, level, src, dst, aux, surface; };
static_assert(sizeof(ALeaf)==32,"leaf ABI");
inline int3 aorigin(ulong key) {
    uint3 p(0);
    for(uint b=0;b<21;b++)for(uint d=0;d<3;d++)p[d]|=uint((key>>(3*b+d))&1ul)<<b;
    return int3(p)-int3(1<<20);
}
inline ulong akey(int3 p) {
    uint3 q=uint3(p+int3(1<<20));ulong k=0;
    for(uint b=0;b<21;b++)for(uint d=0;d<3;d++)k|=ulong((q[d]>>b)&1u)<<(3*b+d);
    return k;
}
inline uint alookup(device const ALeaf *a,uint n,uint arena,ulong k) {
    uint lo=0,hi=n;
    while(lo<hi){uint m=lo+(hi-lo)/2;bool le=a[m].arena<arena||(a[m].arena==arena&&a[m].key<=k);if(le)lo=m+1;else hi=m;}
    if(!lo)return n;uint i=lo-1;
    return a[i].arena==arena&&k-a[i].key<(1ul<<(3*a[i].level))?i:n;
}
inline bool abit(ALeaf a,uint bit) {return bit<63?((a.key>>bit)&1ul):((a.arena>>(bit-63))&1u);}
inline bool a_needs_refine(ALeaf c,device const ALeaf *fine,uint n,device const uint *voxels) {
    if(!c.level)return false;
    uint first=alookup(fine,n,c.arena,c.key);ulong end=c.key+(1ul<<(3*c.level));
    for(uint i=first;i<n&&fine[i].arena==c.arena&&fine[i].key<end;i++)if(voxels[fine[i].source*32+16]==0)return true;
    int3 origin=aorigin(c.key);int side=1<<c.level;
    for(uint f=0;f<6;f++)for(int v=0;v<side;v++)for(int u=0;u<side;u++){
        uint d=f/2;int3 p=origin;p[d]+=(f&1)?side:-1;p[(d+1)%3]+=u;p[(d+2)%3]+=v;
        if(any(p<int3(-(1<<20)))||any(p>=int3(1<<20)))continue;
        uint neighbor=alookup(fine,n,c.arena,akey(p));if(neighbor<n&&voxels[fine[neighbor].source*32+16]==0)return true;
    }
    return false;
}
kernel void adaptive_mesh(device ALeaf *aa [[buffer(0)]],device ALeaf *bb [[buffer(1)]],
    device uint *r [[buffer(2)]],device uint *s [[buffer(3)]],device uint *w [[buffer(4)]],
    device uint *h [[buffer(5)]],constant AMeshParams &u [[buffer(6)]],
    device const uint *voxels [[buffer(7)]],device ALeaf *fine [[buffer(8)]],
    uint gid [[thread_position_in_grid]],uint tid [[thread_index_in_threadgroup]],
    uint group [[threadgroup_position_in_grid]]) {
    threadgroup uint scan[128];
    device ALeaf *a=h[8]?bb:aa; device ALeaf *b=h[8]?aa:bb;
    if(u.mode==20){if(!gid){h[6]=h[2]?0:1;h[5]=0;}return;}
    if(h[2])return;
    if(!h[6])return;
    uint n=h[0],cap=h[7];
    if(u.mode==2) {
        uint v=gid<u.count?w[u.src+gid]:0;scan[tid]=v;threadgroup_barrier(mem_flags::mem_threadgroup);
        for(uint step=1;step<128;step*=2){uint x=tid>=step?scan[tid-step]:0;threadgroup_barrier(mem_flags::mem_threadgroup);scan[tid]+=x;threadgroup_barrier(mem_flags::mem_threadgroup);}
        if(gid<u.count)w[u.dst+gid]=scan[tid]-v;
        if(tid==127)w[u.aux+group]=scan[tid];return;
    }
    if(gid>=u.count)return;
    switch(u.mode) {
    case 1: w[gid]=gid<n?uint(!abit(a[gid],u.bit)):0;break;
    case 3: w[u.dst+gid]+=w[u.src+gid/128];break;
    case 4: if(gid<n){uint zeros=w[u.dst+u.count-1]+w[u.count-1];uint rank=w[u.dst+gid];b[abit(a[gid],u.bit)?zeros+gid-rank:rank]=a[gid];}break;
    case 5: {
        uint flag=0;
        if(gid<n) {
            ALeaf c=a[gid];int3 p=aorigin(c.key);
            if(u.surface)for(uint f=0;f<6;f++){int3 q=p;q[f/2]+=(f&1)?1:-1;if(any(q<int3(-(1<<20)))||any(q>=int3(1<<20))||alookup(a,n,c.arena,akey(q))==n)flag=1;}
        }
        w[gid]=flag;break;
    }
    case 6: if(gid<n)a[gid].flags|=w[gid];break;
    case 7: {
        uint count=gid<n?1:0;
        if(gid<n&&a[gid].level==u.level){uint child=uint((a[gid].key>>(3*u.level))&7ul);
            if(gid>=child){uint first=gid-child;ulong span=1ul<<(3*u.level);bool valid=first+8<=n&&!(a[first].key&(8*span-1));
                for(uint j=0;valid&&j<8;j++)valid=a[first+j].level==u.level&&a[first+j].arena==a[first].arena&&a[first+j].material==a[first].material&&!(a[first+j].flags&1u)&&a[first+j].key==a[first].key+j*span;
                if(valid&&child)count=0;}}
        w[gid]=count;break;
    }
    case 8: if(gid<n&&w[gid]){ALeaf c=a[gid];if(gid+1<n&&!w[gid+1])c.level++;b[w[u.dst+gid]]=c;}break;
    case 9: if(!gid){uint nn=w[u.dst+u.aux-1]+w[u.aux-1];h[5]=uint(nn!=n);h[0]=nn;h[3]+=h[5];h[9]|=h[5];h[6]=h[5];h[8]^=1;}break;
    case 10: {
        uint cell=gid/6,f=gid%6;uint req=0xffffffffu;
        if(cell<n){ALeaf c=a[cell];int3 p=aorigin(c.key);p[f/2]+=(f&1)?(1<<c.level):-1;
            if(all(p>=int3(-(1<<20)))&&all(p<int3(1<<20))){uint other=alookup(a,n,c.arena,akey(p));if(other<n&&a[other].level>c.level+1)req=other;}}
        r[gid]=req;w[gid]=uint(req!=0xffffffffu);break;
    }
    case 11: if(w[gid])s[w[u.dst+gid]]=r[gid];if(!gid)h[1]=w[u.dst+u.count-1]+w[u.count-1];break;
    case 12: w[gid]=gid<h[1]?uint(!((r[gid]>>u.bit)&1u)):0;break;
    case 13: if(gid<h[1]){uint zeros=w[u.dst+u.count-1]+w[u.count-1],rank=w[u.dst+gid];s[((r[gid]>>u.bit)&1u)?zeros+gid-rank:rank]=r[gid];}break;
    case 14: w[gid]=uint(gid<h[1]&&(!gid||r[gid]!=r[gid-1]));break;
    case 15: if(w[gid])s[w[u.dst+gid]]=r[gid];if(!gid)h[1]=w[u.dst+u.count-1]+w[u.count-1];break;
    case 16: {
        uint count=gid<n?1:0;
        if(gid<n){uint lo=0,hi=h[1];while(lo<hi){uint m=lo+(hi-lo)/2;if(r[m]<gid)lo=m+1;else hi=m;}if(lo<h[1]&&r[lo]==gid)count=8;}
        w[gid]=count;break;
    }
    case 17: if(gid<n){uint parts=w[gid],dest=w[u.dst+gid];if(dest+parts<=cap){ALeaf c=a[gid];if(parts==8)c.level--;for(uint j=0;j<parts;j++){b[dest+j]=c;b[dest+j].key+=ulong(j)<<(3*c.level);uint source=alookup(fine,cap,c.arena,b[dest+j].key);b[dest+j].source=fine[source].source;}}}break;
    case 18: w[gid]=uint(gid<n&&gid>0&&a[gid].arena==a[gid-1].arena&&a[gid].key==a[gid-1].key);break;
    case 19: if(!gid){h[2]=w[u.dst+u.aux-1]+w[u.aux-1];h[6]=uint(h[2]==0);}break;
    case 27: if(gid<n)fine[gid]=a[gid];if(!gid)h[9]=1;break;
    case 30: {
        uint count=gid<n?1:0;
        if(gid<n&&a_needs_refine(a[gid],fine,cap,voxels))count=8;
        w[gid]=count;break;
    }
    case 31: {
        uint count=gid<n?1:0;
        if(gid<n&&a[gid].level>1){ALeaf c=a[gid];int3 origin=aorigin(c.key);int side=1<<c.level;
            for(uint f=0;f<6&&count!=8;f++)for(int v=0;v<side&&count!=8;v++)for(int u=0;u<side;u++){
                uint d=f/2;int3 p=origin;p[d]+=(f&1)?side:-1;p[(d+1)%3]+=u;p[(d+2)%3]+=v;
                if(any(p<int3(-(1<<20)))||any(p>=int3(1<<20)))continue;
                uint other=alookup(a,n,c.arena,akey(p));if(other<n&&a[other].level+1<c.level){count=8;break;}}}
        w[gid]=count;break;
    }
    case 23: if(gid<n&&h[8])aa[gid]=bb[gid];break;
    case 24: if(!gid)h[8]=0;break;
    case 25: {
        uint bad=0;
        if(gid<n){ALeaf c=a[gid];int3 p=aorigin(c.key);
            for(uint f=0;f<6;f++){int3 q=p;q[f/2]+=(f&1)?(1<<c.level):-1;
                if(all(q>=int3(-(1<<20)))&&all(q<int3(1<<20))){uint other=alookup(a,n,c.arena,akey(q));if(other<n&&a[other].level>c.level+1)bad=1;}}}
        w[gid]=bad;break;
    }
    case 21: if(!gid){uint nn=w[u.dst+u.aux-1]+w[u.aux-1];if(nn>cap){h[2]=2;h[6]=0;}}break;
    }
}
