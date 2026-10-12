#define ACTIVATION_PROFILE
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

typedef struct {float belief,ground;uint8_t mask,neighbors;} SupportSample;
static SupportSample samples[MAX_VOXELS];
static void check_support(void) {
    staticBeliefsForceFullRefresh=true;refresh_static_voxel_beliefs();
    for(int i=0;i<voxel_count;i++)samples[i]=(SupportSample){voxels[i].freezeBelief,voxels[i].groundSupport,voxels[i].supportMask,voxels[i].neighborSupport};
    recompute_static_freeze_beliefs_path_length();
    for(int i=0;i<voxel_count;i++) {
        Voxel *v=&voxels[i];SupportSample a=samples[i];
        if(a.belief!=v->freezeBelief || a.ground!=v->groundSupport || a.mask!=v->supportMask || a.neighbors!=v->neighborSupport){fprintf(stderr,"Support mismatch at %d (%d,%d,%d) belief=%g/%g mask=%u/%u\n",i,v->gx,v->gy,v->gz,a.belief,v->freezeBelief,a.mask,v->supportMask);abort();}
    }
}
static uint64_t face_key(int x,int y,int z,int face){uint64_t h=(uint32_t)x;h=h*1099511628211ull+(uint32_t)y;h=h*1099511628211ull+(uint32_t)z;return h*7+face;}
static void check_surfaces(void) {
    assert(ensure_static_surface_cache());uint64_t expected=0,actual=0;int ec=0,ac=0;
    for(int v=0;v<voxel_count;v++)if(!voxels[v].simulate && !voxel_is_fluid(&voxels[v])) {
        Voxel *p=&voxels[v];
        for(int z=p->min_gz;z<=p->max_gz;z++)for(int y=p->min_gy;y<=p->max_gy;y++)for(int x=p->min_gx;x<=p->max_gx;x++)if(table_get_static_only(x,y,z)==v)
        for(int f=0;f<6;f++)if(!is_static_solid_cell(x+face_offsets[f][0],y+face_offsets[f][1],z+face_offsets[f][2])){expected+=face_key(x,y,z,f);ec++;}
    }
    for(int p=0;p<patchCount;p++)for(int j=0;j<patches[p].dj;j++)for(int i=0;i<patches[p].di;i++) {
        Patch *v=&patches[p];int a=v->i0+i,b=v->j0+j;
        int x=v->plane==2?v->layer:a,y=v->plane==0?b:v->plane==1?v->layer:a,z=v->plane==0?v->layer:b;
        int f=(v->plane==2?0:v->plane==1?2:4)+!v->positive;
        assert(is_static_solid_cell(x,y,z));assert(!is_static_solid_cell(x+face_offsets[f][0],y+face_offsets[f][1],z+face_offsets[f][2]));actual+=face_key(x,y,z,f);ac++;
    }
    assert(ac==ec && actual==expected);
}
static void check_gpu(void) {
    if(physicsBackend.active!=PHYSICS_BACKEND_GPU_GL43)return;
    assert(gpu_rebuild_static_cells());
    bool cells=!gpuPhysics.static_cell_buffer.id || gpuPhysics.static_cell_buffer.capacity<(size_t)HASH_SIZE*sizeof(GpuStaticCell);
    bool colliders=!gpuPhysics.static_collider_buffer.id || gpuPhysics.static_collider_buffer.capacity<(size_t)gpuPhysics.static_collider_count*sizeof(GpuStaticCollider);
    assert(gpu_buffer_ensure(&gpuPhysics.static_cell_buffer,(size_t)HASH_SIZE*sizeof(GpuStaticCell)));
    assert(gpu_buffer_ensure(&gpuPhysics.static_collider_buffer,(size_t)gpuPhysics.static_collider_count*sizeof(GpuStaticCollider)));
    gpu_upload_static_pages(cells || !gpuStaticIncremental,colliders || !gpuStaticIncremental);
    GpuStaticCell *read=malloc((size_t)HASH_SIZE*sizeof(*read));assert(read);
    gpu_buffer_read(&gpuPhysics.static_cell_buffer,read,(size_t)HASH_SIZE*sizeof(*read),0);
    assert(!memcmp(read,gpuPhysics.static_stage,(size_t)HASH_SIZE*sizeof(*read)));free(read);
    for(int c=0;c<staticCacheChunkCount;c++) {
        GpuCacheChunk *chunk=&gpuCacheChunks[c];size_t size=(size_t)chunk->count*sizeof(GpuStaticCollider);
        GpuStaticCollider *data=malloc(size?size:1);assert(data);gpu_buffer_read(&gpuPhysics.static_collider_buffer,data,size,(size_t)chunk->start*sizeof(*data));
        assert(!size || !memcmp(data,gpuPhysics.static_collider_stage+chunk->start,size));free(data);
    }
    for(int c=0;c<staticCacheChunkCount;c++)for(int k=0;k<staticCacheChunks[c].cellCount;k++) {
        StaticSurfaceCell *cell=&staticSurfaceCells[staticCacheChunks[c].cellStart+k];unsigned h=gpu_static_hash_coord(cell->x,cell->y,cell->z);
        while(gpuPhysics.static_stage[h].value[3]!=-1){GpuStaticCell v=gpuPhysics.static_stage[h];if(v.value[0]==cell->x && v.value[1]==cell->y && v.value[2]==cell->z)break;h=(h+1)&(HASH_SIZE-1);}
        int id=gpuPhysics.static_stage[h].value[3];assert(id!=-1);
        if(table_get_static_only(cell->x,cell->y,cell->z)>=0){assert(id<=-2);continue;}
        int count=0;while(id>=0){assert(id<gpuPhysics.static_collider_count);GpuStaticCollider *v=&gpuPhysics.static_collider_stage[id];int p;memcpy(&p,&v->bounds_min[3],sizeof(p));bool found=false;for(int r=0;r<cell->count;r++)found|=staticSurfacePatchRefs[cell->start+r]==p;assert(found);count++;memcpy(&id,&v->center[3],sizeof(id));assert(count<=cell->count);}assert(count==cell->count);
    }
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Static cache regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    activePlayers=0;clear_world_voxels();
    for(int x=-34;x<=34;x++)for(int z=-1;z<=1;z++)add_static_voxel_at_grid(x,0,z,GRAY,0);
    for(int y=0;y<8;y++){add_static_voxel_at_grid(-32,y,0,GRAY,0);add_static_voxel_at_grid(128,y,0,BLUE,0);}
    for(int x=-32;x<=-20;x++)add_static_voxel_at_grid(x,7,0,GRAY,0);
    check_support();check_surfaces();check_gpu();
    int distant=cache_chunk_find(4,0,0,false);assert(distant>=0);uint64_t revision=staticCacheChunks[distant].revision;
    uint64_t visited=supportVisitedNodes,bytes=gpuStaticBytesUploaded;
    remove_voxel_index(table_get_static_only(-32,4,0));check_support();check_surfaces();check_gpu();
    assert(supportVisitedNodes-visited<(uint64_t)voxel_count);assert(staticCacheChunks[distant].revision==revision);
    if(physicsBackend.active==PHYSICS_BACKEND_GPU_GL43)assert(gpuStaticBytesUploaded-bytes<(uint64_t)HASH_SIZE*sizeof(GpuStaticCell));
    // Alternate changes across seams, reuse freed slots, split and rejoin support.
    int arenaHigh=0,patchHigh=0,cellHigh=0,refHigh=0;
    for(int n=0;n<48;n++) {
        int x=n%2?-33:32,y=1+n%6;
        int old=table_get_static_only(x,y,0);if(old>=0)remove_voxel_index(old);else add_static_voxel_at_grid(x,y,0,GRAY,0);
        check_support();check_surfaces();check_gpu();
        if(n==23){arenaHigh=gpuColliderArena.high;patchHigh=staticPatchArena.high;cellHigh=staticCellArena.high;refHigh=staticRefArena.high;}
    }
    assert(gpuColliderArena.high==arenaHigh && staticPatchArena.high==patchHigh && staticCellArena.high==cellHigh && staticRefArena.high==refHigh);
    add_static_voxel_at_grid(-32,4,0,GRAY,0);check_support();check_surfaces();check_gpu();
    // Restored cells and static-to-dynamic activation use the same invalidation.
    int i=table_get_static_only(-31,0,0);assert(recycle_queue_push(&voxels[i]));remove_voxel_index(i);
    check_support();check_surfaces();check_gpu();
    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;for(int f=0;f<REGEN_PREVIEW_FRAMES+2;f++)recycle_dead_voxels(1.0f/60);
    check_support();check_surfaces();check_gpu();
    assert(activate_static_voxel_for_tether(table_get_static_only(-24,7,0),0,ACTIVATION_TETHER_BELIEF));check_support();check_surfaces();check_gpu();
    clear_world_voxels();
    for(int x=29;x<=35;x++)for(int z=-2;z<=2;z++)add_static_voxel_at_grid(x,0,z,GRAY,0);
    check_support();check_surfaces();check_gpu();
    // Chunk seams must behave as flat surfaces rather than produce edge catches.
    for(int n=-2;n<=2;n++){Particle p={0};p.predicted_pos=(Vector3){16+n*.025f,.52f,.25f};collide_particle_with_static_surface(&p,.125f);assert(isfinite(p.predicted_pos.y));assert(fabsf(p.predicted_pos.x-(16+n*.025f))<.001f);assert(fabsf(p.predicted_pos.y-.625f)<.001f);}
    // Same grid bounds but changed physical ground contact and collider size.
    int coarse=add_static_voxel_at_grid(96,1,0,GRAY,0);check_support();check_surfaces();check_gpu();
    voxel_table_unregister(&voxels[coarse]);voxels[coarse].rest_edge=1.5f;voxel_table_register(&voxels[coarse],coarse);mark_static_hash_dirty();
    check_support();check_surfaces();check_gpu();
    // Exercise backward-shift deletion through a deliberately colliding probe chain.
    if(physicsBackend.active==PHYSICS_BACKEND_GPU_GL43) {
        assert(gpu_cache_insert(0,123,456,-2));assert(gpu_cache_insert(HASH_SIZE,123,456,-2));
        gpu_cache_remove((CacheKey){0,123,456});unsigned h=gpu_static_hash_coord(HASH_SIZE,123,456);
        assert(gpuPhysics.static_stage[h].value[0]==HASH_SIZE && gpuPhysics.static_stage[h].value[3]==-2);
        gpu_cache_remove((CacheKey){HASH_SIZE,123,456});
    }
    // Capacity/failure recovery publishes a complete full cache, then returns
    // to incremental operation on the next edit.
    staticCacheOverflow=true;remove_voxel_index(table_get_static_only(29,0,-2));
    check_surfaces();assert(!staticCacheReady);check_gpu();
    add_static_voxel_at_grid(29,0,-2,GRAY,0);check_support();check_surfaces();assert(staticCacheReady);check_gpu();
    uint64_t epoch=staticCacheEpoch;
    add_static_voxel_at_grid(35,1,0,GRAY,0);mark_static_hash_dirty();
    check_support();check_surfaces();assert(staticCacheEpoch!=epoch);check_gpu();
    clear_world_voxels();check_support();check_surfaces();check_gpu();
    shutdown_pbd_thread_pool();gpu_physics_shutdown();assert(!supportNodes && !patches && !staticSurfaceHash);CloseWindow();puts("Static support/collision cache regression passed");return 0;
}
