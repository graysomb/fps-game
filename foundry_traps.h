/* Arena-only gameplay state; stable voxel tags survive compaction and GPU readback. */
#define FOUNDRY_TRAPS 32
#define FOUNDRY_TRAP_PAYLOAD 256
typedef enum { FOUNDRY_DROP, FOUNDRY_SWITCH, FOUNDRY_TIP } FoundryTrapKind;
typedef struct {
    FoundryTrapKind kind;
    int x0,x1,y0,y1,z0,z1,rx,ry,rz,owner;
    bool active,pending,fourPosts;
    unsigned activation;
    uint32_t victims;
    Color color;
    Voxel original[FOUNDRY_TRAP_PAYLOAD+32];
    int originalCount;
} FoundryTrap;
static FoundryTrap foundryTraps[FOUNDRY_TRAPS];
static int foundryTrapCount;
static int foundryFragileReferenceCells,foundryFragileCells;
static int foundryFragileReferenceBoxes[32][6],foundryFragileReferenceCount;
static bool foundry_trap_contains(const FoundryTrap *t,int x,int y,int z) {
    return x>=t->x0 && x<=t->x1 && y>=t->y0 && y<=t->y1 && z>=t->z0 && z<=t->z1;
}
static void foundry_trap_define(FoundryTrapKind kind,int x0,int x1,int y0,int y1,int z0,int z1,int rx,int ry,int rz) {
    if(foundryTrapCount>=FOUNDRY_TRAPS)return;
    FoundryTrap *t=&foundryTraps[foundryTrapCount++];
    *t=(FoundryTrap){.kind=kind,.x0=x0,.x1=x1,.y0=y0,.y1=y1,.z0=z0,.z1=z1,.rx=rx,.ry=ry,.rz=rz,.owner=-1,.color={127,173,178,255}};
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
