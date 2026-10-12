/* Whole rendering chunks are scheduling units; only due missing cells are restored. */
#define REGEN_CELL_COUNT (CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE)
#define REGEN_PREVIEW_FRAMES 36
#define REGEN_MAX_CANDIDATES 2048

typedef struct {
    bool present;
    Color color;
    int type, owner;
} RegenerationCell;
typedef struct { int x, y, z, oldest; } RegenerationChunk;
static struct {
    bool active, meshDirty;
    int cx, cy, cz, capturedFrame, startedFrame, count, nextTry, receiveTotal;
    uint32_t serial;
    double remoteStart;
    float elapsed;
    RegenerationCell cells[REGEN_CELL_COUNT];
} regeneration;
static RegenerationCell regenerationScratch[REGEN_CELL_COUNT];
static uint64_t regenerationCommits;
static double regenerationLastCommitMs;
static void regeneration_release_mesh(void);

static int regen_min(int a,int b){return a<b?a:b;}
static int regen_max(int a,int b){return a>b?a:b;}
static int regeneration_cell_index(int x, int y, int z) {
    return (z * CHUNK_SIZE + y) * CHUNK_SIZE + x;
}
static void regeneration_bounds(const Voxel *v, int *x0,int *x1,int *y0,int *y1,int *z0,int *z1) {
    bool original=v->orig_min_gx<=v->orig_max_gx && v->orig_min_gy<=v->orig_max_gy && v->orig_min_gz<=v->orig_max_gz;
    *x0=original?v->orig_min_gx:v->rest_min_gx; *x1=original?v->orig_max_gx:v->rest_max_gx;
    *y0=original?v->orig_min_gy:v->rest_min_gy; *y1=original?v->orig_max_gy:v->rest_max_gy;
    *z0=original?v->orig_min_gz:v->rest_min_gz; *z1=original?v->orig_max_gz:v->rest_max_gz;
}
static int regeneration_collect(int cx,int cy,int cz,int cutoff,RegenerationCell *cells) {
    memset(cells,0,sizeof(RegenerationCell)*REGEN_CELL_COUNT);
    int count=0,ox=cx*CHUNK_SIZE,oy=cy*CHUNK_SIZE,oz=cz*CHUNK_SIZE;
    for(int n=0;n<recycleQueueCount;n++) {
        const Voxel *v=&recycleQueue[(recycleQueueHead+n)%MAX_VOXELS];
        if(cutoff-v->lifeFrames<RECYCLE_STATIC_RESTORE_DELAY)continue;
        int x0,x1,y0,y1,z0,z1;regeneration_bounds(v,&x0,&x1,&y0,&y1,&z0,&z1);
        x0=regen_max(x0,ox);x1=regen_min(x1,ox+CHUNK_MASK);
        y0=regen_max(y0,oy);y1=regen_min(y1,oy+CHUNK_MASK);
        z0=regen_max(z0,oz);z1=regen_min(z1,oz+CHUNK_MASK);
        for(int z=z0;z<=z1;z++)for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++) {
            if(table_get_static_only(x,y,z)>=0)continue; /* Never overwrite construction. */
            RegenerationCell *cell=&cells[regeneration_cell_index(x-ox,y-oy,z-oz)];
            if(cell->present)continue;
            cell->present=true;cell->color=v->color;cell->type=v->type;cell->owner=v->owner;count++;
        }
    }
    return count;
}
static bool regeneration_intersects(BoundingBox bounds,int cx,int cy,int cz,const RegenerationCell *cells) {
    int ox=cx*CHUNK_SIZE,oy=cy*CHUNK_SIZE,oz=cz*CHUNK_SIZE;
    int x0=regen_max(0,(int)floorf(bounds.min.x/VOXEL_SIZE)-ox),x1=regen_min(CHUNK_MASK,(int)floorf((bounds.max.x-.001f)/VOXEL_SIZE)-ox);
    int y0=regen_max(0,(int)floorf(bounds.min.y/VOXEL_SIZE)-oy),y1=regen_min(CHUNK_MASK,(int)floorf((bounds.max.y-.001f)/VOXEL_SIZE)-oy);
    int z0=regen_max(0,(int)floorf(bounds.min.z/VOXEL_SIZE)-oz),z1=regen_min(CHUNK_MASK,(int)floorf((bounds.max.z-.001f)/VOXEL_SIZE)-oz);
    for(int z=z0;z<=z1;z++)for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++)
        if(cells[regeneration_cell_index(x,y,z)].present)return true;
    return false;
}
static bool regeneration_blocked(int cx,int cy,int cz,const RegenerationCell *cells) {
    for(int i=0;i<activePlayers;i++) {
        Player *p=&players[i];if(p->respawn_timer>0)continue;
        BoundingBox box={{p->pos.x-PLAYER_RADIUS,p->pos.y-BASE_EYE_HEIGHT,p->pos.z-PLAYER_RADIUS},
                        {p->pos.x+PLAYER_RADIUS,p->pos.y+.25f,p->pos.z+PLAYER_RADIUS}};
        if(regeneration_intersects(box,cx,cy,cz,cells))return true;
    }
    for(int i=0;i<voxel_count;i++)if(voxels[i].simulate) {
        VoxelWorldBounds b;voxel_world_bounds(&voxels[i],&b);
        if(regeneration_intersects((BoundingBox){{b.minx,b.miny,b.minz},{b.maxx,b.maxy,b.maxz}},cx,cy,cz,cells))return true;
    }
    return false;
}
static void regeneration_reset(void) {
    regeneration_release_mesh();
    uint32_t serial=regeneration.serial;
    memset(&regeneration,0,sizeof(regeneration));regeneration.serial=serial;
}
static void regeneration_cancel(void) {
    regeneration_net_send(-1,false);
    regeneration_reset();regeneration.nextTry=recycleFrameCounter+1;
}
static int regeneration_compare(const void *a,const void *b) {
    const RegenerationChunk *x=a,*y=b;
    if(x->oldest!=y->oldest)return x->oldest<y->oldest?-1:1;
    if(x->y!=y->y)return x->y<y->y?-1:1;
    if(x->z!=y->z)return x->z<y->z?-1:1;
    return (x->x>y->x)-(x->x<y->x);
}
static void regeneration_step(float dt) {
    if(netTransport.role==NET_ROLE_CLIENT)return;
    if(regeneration.active) {
        regeneration.elapsed+=fmaxf(0,dt);
        if(regeneration.elapsed+.00001f<.6f)return;
        int count=regeneration_collect(regeneration.cx,regeneration.cy,regeneration.cz,regeneration.capturedFrame,regenerationScratch);
        if(count!=regeneration.count || memcmp(regenerationScratch,regeneration.cells,sizeof(regenerationScratch)) ||
           regeneration_blocked(regeneration.cx,regeneration.cy,regeneration.cz,regeneration.cells) ||
           count>MAX_VOXELS-voxel_count || count>(MAX_PARTICLES-particle_pool_count+free_particle_count)/8) {
            regeneration_cancel();return;
        }
        double start=GetTime();int before=voxel_count;
        int ox=regeneration.cx*CHUNK_SIZE,oy=regeneration.cy*CHUNK_SIZE,oz=regeneration.cz*CHUNK_SIZE;
        staticEditBatch=true;
        bool complete=true;
        for(int n=0;n<REGEN_CELL_COUNT;n++)if(regeneration.cells[n].present) {
            int x=ox+n%CHUNK_SIZE,y=oy+(n/CHUNK_SIZE)%CHUNK_SIZE,z=oz+n/(CHUNK_SIZE*CHUNK_SIZE);
            RegenerationCell *c=&regeneration.cells[n];
            int idx=addVoxel((x+.5f)*VOXEL_SIZE,(y+.5f)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE,true,false,c->color,c->type);
            if(idx<0){complete=false;break;}
            voxels[idx].owner=c->owner;
            voxels[idx].activationCooldownFrames=STATIC_REBUILD_ACTIVATION_COOLDOWN_FRAMES;
        }
        /* Capacity is preflighted; retain rollback for unexpected allocation failure. */
        if(!complete)while(voxel_count>before)remove_voxel_index(voxel_count-1);
        staticEditBatch=false;
        mark_static_hash_dirty();staticBeliefsForceFullRefresh=true;
        mark_chunk_dirty_for_grid(ox,oy,oz);
        mark_chunks_dirty_for_region(ox,ox+CHUNK_MASK,oy,oy+CHUNK_MASK,oz,oz+CHUNK_MASK);
        refresh_surfaces_in_region(ox,ox+CHUNK_MASK,oy,oy+CHUNK_MASK,oz,oz+CHUNK_MASK);
        if(complete) {
            recycle_queue_prune_fully_restored_in_region(ox,ox+CHUNK_MASK,oy,oy+CHUNK_MASK,oz,oz+CHUNK_MASK);
            regenerationCommits++;
        }
        regenerationLastCommitMs=(GetTime()-start)*1000;
        regeneration_cancel();return;
    }
    if(!recycleQueueCount || recycleFrameCounter<regeneration.nextTry)return;
    regeneration.nextTry=recycleFrameCounter+15;
    static RegenerationChunk candidates[REGEN_MAX_CANDIDATES];int count=0;
    for(int n=0;n<recycleQueueCount;n++) {
        const Voxel *v=&recycleQueue[(recycleQueueHead+n)%MAX_VOXELS];
        if(recycleFrameCounter-v->lifeFrames<RECYCLE_STATIC_RESTORE_DELAY)continue;
        int x0,x1,y0,y1,z0,z1;regeneration_bounds(v,&x0,&x1,&y0,&y1,&z0,&z1);
        for(int z=grid_to_chunk_coord(z0);z<=grid_to_chunk_coord(z1);z++)
        for(int y=grid_to_chunk_coord(y0);y<=grid_to_chunk_coord(y1);y++)
        for(int x=grid_to_chunk_coord(x0);x<=grid_to_chunk_coord(x1);x++) {
            int k=0;while(k<count && (candidates[k].x!=x || candidates[k].y!=y || candidates[k].z!=z))k++;
            if(k==count){if(count==REGEN_MAX_CANDIDATES)continue;candidates[count++]=(RegenerationChunk){x,y,z,v->lifeFrames};}
            else if(v->lifeFrames<candidates[k].oldest)candidates[k].oldest=v->lifeFrames;
        }
    }
    qsort(candidates,count,sizeof(*candidates),regeneration_compare);
    for(int k=0;k<count;k++) {
        RegenerationChunk c=candidates[k];
        int cells=regeneration_collect(c.x,c.y,c.z,recycleFrameCounter,regeneration.cells);
        if(!cells) {
            recycle_queue_prune_fully_restored_in_region(c.x*CHUNK_SIZE,c.x*CHUNK_SIZE+CHUNK_MASK,c.y*CHUNK_SIZE,c.y*CHUNK_SIZE+CHUNK_MASK,c.z*CHUNK_SIZE,c.z*CHUNK_SIZE+CHUNK_MASK);
            continue;
        }
        if(cells>MAX_VOXELS-voxel_count || cells>(MAX_PARTICLES-particle_pool_count+free_particle_count)/8 || regeneration_blocked(c.x,c.y,c.z,regeneration.cells))continue;
        regeneration.active=regeneration.meshDirty=true;
        regeneration.cx=c.x;regeneration.cy=c.y;regeneration.cz=c.z;regeneration.count=cells;
        regeneration.capturedFrame=regeneration.startedFrame=recycleFrameCounter;
        if(!++regeneration.serial)++regeneration.serial;
        regeneration_net_send(-1,true);return;
    }
}
