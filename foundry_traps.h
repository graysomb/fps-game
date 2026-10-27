/* Arena-only gameplay state; stable voxel tags survive compaction and GPU readback. */
#define FOUNDRY_TRAPS 32
#define FOUNDRY_TRAP_PAYLOAD 2048
typedef enum { FOUNDRY_DROP, FOUNDRY_SWITCH, FOUNDRY_TIP, FOUNDRY_BAY, FOUNDRY_PILLAR } FoundryTrapKind;
typedef struct {
    FoundryTrapKind kind;
    int x0,x1,y0,y1,z0,z1,rx,ry,rz,owner;
    bool active,pending,fourPosts;
    unsigned activation;
    uint32_t victims;
    Color color;
    Voxel *original;
    int originalCount;
} FoundryTrap;
static FoundryTrap foundryTraps[FOUNDRY_TRAPS];
static int foundryTrapCount;
static int foundryFragileReferenceCells,foundryFragileCells;
static int foundryFragileReferenceBoxes[32][6],foundryFragileReferenceCount;
static int foundry_bay_index(const FoundryTrap *t,int x,int y,int z) {
    if(x<t->x0 || x>t->x1 || y<t->y0 || y>t->y1 || z<t->z0 || z>t->z1)return -1;
    int width=t->x1-t->x0+1,height=t->y1-t->y0+1;
    int key=((z-t->z0)*height+y-t->y0)*width+x-t->x0;
    int lo=0,hi=t->originalCount-1;
    while(lo<=hi){int mid=(lo+hi)/2;const Voxel *v=&t->original[mid];
        int other=((v->gz-t->z0)*height+v->gy-t->y0)*width+v->gx-t->x0;
        if(key==other)return mid;if(key<other)hi=mid-1;else lo=mid+1;
    }
    return -1;
}
static bool foundry_bay_crack(const FoundryTrap *t,int x,int y,int z) {
    int k=foundry_bay_index(t,x,y,z);
    return k>=0 && ColorToInt(t->original[k].color)==ColorToInt(ORANGE);
}
static bool foundry_fractured(const FoundryTrap *t) { return t->kind==FOUNDRY_BAY || t->kind==FOUNDRY_PILLAR; }
static bool foundry_trap_contains(const FoundryTrap *t,int x,int y,int z) {
    if(foundry_fractured(t))return y>0 && foundry_bay_index(t,x,y,z)>=0 && !foundry_bay_crack(t,x,y,z);
    return x>=t->x0 && x<=t->x1 && y>=t->y0 && y<=t->y1 && z>=t->z0 && z<=t->z1;
}
static void foundry_trap_define(FoundryTrapKind kind,int x0,int x1,int y0,int y1,int z0,int z1,int rx,int ry,int rz) {
    if(foundryTrapCount>=FOUNDRY_TRAPS)return;
    int cells=(x1-x0+1)*(y1-y0+1)*(z1-z0+1);
    if(cells<1 || cells>FOUNDRY_TRAP_PAYLOAD)return;
    Voxel *original=calloc(cells+32,sizeof(Voxel));if(!original)return;
    FoundryTrap *t=&foundryTraps[foundryTrapCount++];
    *t=(FoundryTrap){.kind=kind,.x0=x0,.x1=x1,.y0=y0,.y1=y1,.z0=z0,.z1=z1,.rx=rx,.ry=ry,.rz=rz,.owner=-1,.color={127,173,178,255},.original=original};
    add_static_box_at_grid(x0,x1,y0,y1,z0,z1,t->color);
    add_static_voxel_at_grid(rx,ry,rz,ORANGE,0);
    for(int z=z0;z<=z1;z++)for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++) {
        int id=table_get_static_only(x,y,z);
        if(id>=0 && t->originalCount<FOUNDRY_TRAP_PAYLOAD)t->original[t->originalCount++]=voxels[id];
    }
    int release=table_get_static_only(rx,ry,rz);
    if(release>=0)t->original[t->originalCount++]=voxels[release];
}

static bool foundry_trap_trigger(const FoundryTrap *t,int x,int y,int z) {
    if(foundry_fractured(t))return foundry_bay_crack(t,x,y,z);
    if(t->fourPosts)return y==t->ry && (x==t->x0 || x==t->x1) && (z==t->z0 || z==t->z1);
    return x==t->rx && y==t->ry && z==t->rz;
}
static void foundry_trap_posts(FoundryTrap *t) {
    t->fourPosts=true;
    for(int x=t->x0;x<=t->x1;x+=t->x1-t->x0)for(int z=t->z0;z<=t->z1;z+=t->z1-t->z0)
        for(int y=0;y<t->y0;y++) {
            if(x==t->rx && y==t->ry && z==t->rz)continue;
            add_static_voxel_at_grid(x,y,z,y>=t->ry?ORANGE:(Color){82,102,117,255},0);
            if(y>=t->ry){int id=table_get_static_only(x,y,z);if(id>=0)t->original[t->originalCount++]=voxels[id];}
        }
}

/* Midpoint displacement produces repeatable jagged fractures on an unfolded
 * wall/roof/wall surface. Branch endpoints share their parent's exact cell. */
static unsigned foundry_crack_random(unsigned *seed) {
    *seed=*seed*1664525u+1013904223u;return *seed;
}
static void foundry_crack_path(bool mask[64][64],int ax,int ay,int bx,int by,unsigned *seed) {
    int dx=abs(bx-ax),dy=abs(by-ay);
    if(dx+dy<=2){
        mask[ay][ax]=true;
        while(ax!=bx){ax+=bx>ax?1:-1;mask[ay][ax]=true;}
        while(ay!=by){ay+=by>ay?1:-1;mask[ay][ax]=true;}
        return;
    }
    int mx=(ax+bx)/2,my=(ay+by)/2;
    int jitter=(int)(foundry_crack_random(seed)%3)-1;
    if(dx>dy)my+=jitter;else mx+=jitter;
    mx=mx<1?1:mx>62?62:mx;my=my<1?1:my>62?62:my;
    /* Dominant-axis midpoint guarantees termination. */
    foundry_crack_path(mask,ax,ay,mx,my,seed);
    foundry_crack_path(mask,mx,my,bx,by,seed);
}
/* Register only existing authored cells: open bay interiors stay empty. */
static void foundry_bay_trap(int x,int z) {
    if(foundryTrapCount>=FOUNDRY_TRAPS)return;
    FoundryTrap t={.kind=FOUNDRY_BAY,.x0=(int)floorf((x-6)*1.875f),.x1=(int)ceilf((x+7)*1.875f)-1,
        .z0=(int)floorf((z-5)*1.875f),.z1=(int)ceilf((z+6)*1.875f)-1,.y0=0,.y1=11,.ry=2,.owner=-1};
    int count=0;
    for(int gz=t.z0;gz<=t.z1;gz++)for(int gy=0;gy<=11;gy++)for(int gx=t.x0;gx<=t.x1;gx++)if(table_get_static_only(gx,gy,gz)>=0)count++;
    if(count>FOUNDRY_TRAP_PAYLOAD){TraceLog(LOG_ERROR,"Foundry bay exceeds trap capacity");return;}
    t.original=calloc(count,sizeof(Voxel));if(!t.original)return;
    int width=t.x1-t.x0+1;
    int back=(int)floorf((z+4)*1.875f),depth=back-t.z0;
    int height=22+depth;
    bool cracks[64][64]={0};
    unsigned seed=0x9e3779b9u ^ (unsigned)x*73856093u ^ (unsigned)z*19349663u;
    int rootX=width/2,rootY=11+depth/2;
    foundry_crack_path(cracks,0,rootY,rootX,rootY,&seed);
    foundry_crack_path(cracks,rootX,rootY,width-1,rootY+2,&seed);
    for(int branch=0;branch<3;branch++){
        int u=3+branch*(width-6)/2;
        foundry_crack_path(cracks,rootX,rootY,u,11,&seed);
        foundry_crack_path(cracks,u,11,u,1,&seed);
        foundry_crack_path(cracks,rootX,rootY,u,11+depth,&seed);
        foundry_crack_path(cracks,u,11+depth,u,height-1,&seed);
        /* Secondary branches divide walls without deleting a whole band. */
        if(branch!=1){
            int edge=branch==0?0:width-1;
            foundry_crack_path(cracks,u,11,edge,4,&seed);
            foundry_crack_path(cracks,u,11+depth,edge,height-4,&seed);
        }
    }
    t.rx=t.x0+3;t.ry=1;t.rz=t.z0;
    for(int gz=t.z0;gz<=t.z1;gz++)for(int gy=0;gy<=11;gy++)for(int gx=t.x0;gx<=t.x1;gx++){
        int id=table_get_static_only(gx,gy,gz);if(id<0)continue;
        int v=gy==11?11+(gz-t.z0<depth?gz-t.z0:depth):
            gz<(t.z0+back)/2?gy:height-gy;
        if(gy>0 && cracks[v][gx-t.x0]){add_static_voxel_at_grid(gx,gy,gz,ORANGE,0);id=table_get_static_only(gx,gy,gz);}
        t.original[t.originalCount++]=voxels[id];
    }
    foundryTraps[foundryTrapCount++]=t;
}

/* Four metres square, ten metres tall; fracture sheets continue through the
 * solid interior so visible cracks also separate the physical chunks. */
static void foundry_corner_pillar(int sx,int sz) {
    if(foundryTrapCount>=FOUNDRY_TRAPS)return;
    int cx=sx*70,cz=sz*70;
    FoundryTrap t={.kind=FOUNDRY_PILLAR,.x0=cx-4,.x1=cx+3,.y0=0,.y1=19,
        .z0=cz-4,.z1=cz+3,.owner=-1};
    t.original=calloc(8*8*20,sizeof(Voxel));if(!t.original)return;
    unsigned seed=0xa137u^(unsigned)cx*73856093u^(unsigned)cz*19349663u;
    bool cracks[64][64]={0};
    for(int y=4;y<=16;y+=6)foundry_crack_path(cracks,0,y,7,y+1,&seed);
    foundry_crack_path(cracks,3,4,4,19,&seed);
    bool marked=false;
    for(int z=t.z0;z<=t.z1;z++)for(int y=0;y<=19;y++)for(int x=t.x0;x<=t.x1;x++){
        bool crack=y>0 && cracks[y][x-t.x0];
        Color color=crack?ORANGE:(Color){127,173,178,255};
        add_static_voxel_at_grid(x,y,z,color,0);
        int id=table_get_static_only(x,y,z);
        if(id>=0)t.original[t.originalCount++]=voxels[id];
        if(crack && !marked){t.rx=x;t.ry=y;t.rz=z;marked=true;}
    }
    foundryTraps[foundryTrapCount++]=t;
}
