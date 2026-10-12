#define GPU_STATIC_PAGE 256
#define GPU_STATIC_PAGES (HASH_SIZE/GPU_STATIC_PAGE)
typedef struct {int x,y,z;} CacheKey;
typedef struct {
    uint64_t revision;
    int start,capacity,count,keyCount,keyCapacity;
    CacheKey *keys;
    bool upload;
} GpuCacheChunk;
static GpuCacheChunk gpuCacheChunks[STATIC_CACHE_CHUNKS];
static CacheArena gpuColliderArena;
static uint64_t gpuCacheEpoch;
static bool gpuStaticPages[GPU_STATIC_PAGES],gpuStaticIncremental;
static CacheKey *gpuOldKeys;static int gpuOldKeyCapacity;
static uint32_t gpuKeyStamp[HASH_SIZE],gpuKeySerial;
static uint64_t gpuStaticBytesUploaded;
static double gpuStaticRebuildMs;
static void gpu_static_cache_reset(void) {
    for(int c=0;c<STATIC_CACHE_CHUNKS;c++)free(gpuCacheChunks[c].keys);
    memset(gpuCacheChunks,0,sizeof(gpuCacheChunks));gpuColliderArena=(CacheArena){0};
    free(gpuOldKeys);gpuOldKeys=NULL;gpuOldKeyCapacity=0;
    memset(gpuKeyStamp,0,sizeof(gpuKeyStamp));gpuKeySerial=0;
    gpuCacheEpoch=0;gpuStaticIncremental=false;memset(gpuStaticPages,0,sizeof(gpuStaticPages));
}
static void gpu_static_write(int index,GpuStaticCell value) {
    if(memcmp(&gpuPhysics.static_stage[index],&value,sizeof(value))){gpuPhysics.static_stage[index]=value;gpuStaticPages[index/GPU_STATIC_PAGE]=true;}
}
static bool gpu_cache_insert(int x,int y,int z,int value) {
    unsigned h=gpu_static_hash_coord(x,y,z);
    for(int n=0;n<HASH_SIZE;n++,h=(h+1)&(HASH_SIZE-1)) {
        GpuStaticCell *c=&gpuPhysics.static_stage[h];
        if(c->value[3]==-1 || (c->value[0]==x && c->value[1]==y && c->value[2]==z)){
            gpu_static_write(h,(GpuStaticCell){{x,y,z,value}});gpuKeyStamp[h]=gpuKeySerial;return true;
        }
    }return false;
}
static void gpu_cache_remove(CacheKey key) {
    unsigned mask=HASH_SIZE-1,h=gpu_static_hash_coord(key.x,key.y,key.z);bool found=false;
    for(int n=0;n<HASH_SIZE;n++,h=(h+1)&mask){GpuStaticCell *c=&gpuPhysics.static_stage[h];if(c->value[3]==-1)return;if(c->value[0]==key.x && c->value[1]==key.y && c->value[2]==key.z){found=true;break;}}
    if(!found)return;
    unsigned hole=h;
    for(unsigned j=(h+1)&mask;j!=h && gpuPhysics.static_stage[j].value[3]!=-1;j=(j+1)&mask) {
        GpuStaticCell v=gpuPhysics.static_stage[j];unsigned home=gpu_static_hash_coord(v.value[0],v.value[1],v.value[2]);
        if(((j-home)&mask)>=((hole-home)&mask)){gpu_static_write(hole,v);gpuKeyStamp[hole]=gpuKeyStamp[j];hole=j;}
    }
    gpu_static_write(hole,(GpuStaticCell){{0,0,0,-1}});gpuKeyStamp[hole]=0;
}
static bool gpu_cache_key(GpuCacheChunk *chunk,int x,int y,int z,int value) {
    if(!cache_array((void**)&chunk->keys,&chunk->keyCapacity,chunk->keyCount+1,sizeof(CacheKey)))return false;
    chunk->keys[chunk->keyCount++]=(CacheKey){x,y,z};return gpu_cache_insert(x,y,z,value);
}
static bool gpu_build_static_chunks(void) {
    if(!gpuPhysics.static_stage){gpuPhysics.static_stage=malloc((size_t)HASH_SIZE*sizeof(GpuStaticCell));if(!gpuPhysics.static_stage)return false;}
    if(gpuCacheEpoch!=staticCacheEpoch) {
        gpu_static_cache_reset();
        for(int i=0;i<HASH_SIZE;i++)gpuPhysics.static_stage[i]=(GpuStaticCell){{0,0,0,-1}};
        memset(gpuStaticPages,1,sizeof(gpuStaticPages));gpuCacheEpoch=staticCacheEpoch;
    }
    for(int c=0;c<staticCacheChunkCount;c++) {
        StaticCacheChunk *src=&staticCacheChunks[c];GpuCacheChunk *dst=&gpuCacheChunks[c];
        if(dst->revision==src->revision)continue;
        int oldKeyCount=dst->keyCount;
        if(!cache_array((void**)&gpuOldKeys,&gpuOldKeyCapacity,oldKeyCount,sizeof(CacheKey)))return false;
        if(oldKeyCount)memcpy(gpuOldKeys,dst->keys,(size_t)oldKeyCount*sizeof(CacheKey));
        dst->keyCount=0;if(++gpuKeySerial==0){memset(gpuKeyStamp,0,sizeof(gpuKeyStamp));gpuKeySerial=1;}
        int ox=src->x*CHUNK_SIZE,oy=src->y*CHUNK_SIZE,oz=src->z*CHUNK_SIZE;
        int boxes=0;
        for(int i=0;i<voxel_count;i++){Voxel *v=&voxels[i];if(!v->simulate && v->max_gx>=ox && v->min_gx<ox+CHUNK_SIZE && v->max_gy>=oy && v->min_gy<oy+CHUNK_SIZE && v->max_gz>=oz && v->min_gz<oz+CHUNK_SIZE)boxes++;}
        if(!cache_reserve(&gpuColliderArena,&dst->start,&dst->capacity,src->refCount+boxes))return false;
        if((size_t)gpuColliderArena.high>gpuPhysics.static_collider_capacity){
            size_t capacity=gpuPhysics.static_collider_capacity?gpuPhysics.static_collider_capacity:256;
            while(capacity<(size_t)gpuColliderArena.high)capacity*=2;
            if(!gpu_resize_array((void**)&gpuPhysics.static_collider_stage,capacity,sizeof(GpuStaticCollider)))return false;
            gpuPhysics.static_collider_capacity=capacity;
        }
        dst->count=0;
        for(int k=0;k<src->cellCount;k++) {
            StaticSurfaceCell *cell=&staticSurfaceCells[src->cellStart+k];int head=-1;
            for(int r=0;r<cell->count;r++) {
                int p=staticSurfacePatchRefs[cell->start+r],id=dst->start+dst->count++;
                gpu_fill_static_patch(&gpuPhysics.static_collider_stage[id],&patches[p],p,head);head=id;
            }
            if(table_get_static_only(cell->x,cell->y,cell->z)<0 && !gpu_cache_key(dst,cell->x,cell->y,cell->z,head))return false;
        }
        for(int i=0;i<voxel_count;i++) {
            Voxel *v=&voxels[i];if(v->simulate || v->max_gx<ox || v->min_gx>=ox+CHUNK_SIZE || v->max_gy<oy || v->min_gy>=oy+CHUNK_SIZE || v->max_gz<oz || v->min_gz>=oz+CHUNK_SIZE)continue;
            int id=dst->start+dst->count++;GpuStaticCollider *box=&gpuPhysics.static_collider_stage[id];VoxelWorldBounds b;voxel_world_bounds(v,&b);
            *box=(GpuStaticCollider){{v->pos.x,v->pos.y,v->pos.z,0},{b.minx,b.miny,b.minz,0},{b.maxx,b.maxy,b.maxz,0}};
            int x0=v->min_gx>ox?v->min_gx:ox,x1=v->max_gx<ox+CHUNK_MASK?v->max_gx:ox+CHUNK_MASK;
            int y0=v->min_gy>oy?v->min_gy:oy,y1=v->max_gy<oy+CHUNK_MASK?v->max_gy:oy+CHUNK_MASK;
            int z0=v->min_gz>oz?v->min_gz:oz,z1=v->max_gz<oz+CHUNK_MASK?v->max_gz:oz+CHUNK_MASK;
            for(int z=z0;z<=z1;z++)for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++)if(!gpu_cache_key(dst,x,y,z,-id-2))return false;
        }
        for(int k=0;k<oldKeyCount;k++) {
            CacheKey key=gpuOldKeys[k];unsigned h=gpu_static_hash_coord(key.x,key.y,key.z);
            for(int probe=0;probe<HASH_SIZE && gpuPhysics.static_stage[h].value[3]!=-1;probe++,h=(h+1)&(HASH_SIZE-1)) {
                GpuStaticCell *entry=&gpuPhysics.static_stage[h];
                if(entry->value[0]==key.x && entry->value[1]==key.y && entry->value[2]==key.z){if(gpuKeyStamp[h]!=gpuKeySerial)gpu_cache_remove(key);break;}
            }
        }
        dst->revision=src->revision;dst->upload=true;
    }
    gpuPhysics.static_collider_count=gpuColliderArena.high;gpuStaticIncremental=true;return true;
}
static size_t gpu_upload_static_pages(bool fullCells,bool fullColliders) {
    size_t bytes=0;
    if(fullCells){size_t size=(size_t)HASH_SIZE*sizeof(GpuStaticCell);gpu_buffer_update(&gpuPhysics.static_cell_buffer,gpuPhysics.static_stage,size,0);bytes+=size;}
    else for(int p=0;p<GPU_STATIC_PAGES;) {
        if(!gpuStaticPages[p]){p++;continue;}int first=p;while(p<GPU_STATIC_PAGES && gpuStaticPages[p])p++;
        size_t offset=(size_t)first*GPU_STATIC_PAGE*sizeof(GpuStaticCell),size=(size_t)(p-first)*GPU_STATIC_PAGE*sizeof(GpuStaticCell);
        gpu_buffer_update(&gpuPhysics.static_cell_buffer,(char*)gpuPhysics.static_stage+offset,size,offset);bytes+=size;
    }
    if(fullColliders){size_t size=(size_t)gpuPhysics.static_collider_count*sizeof(GpuStaticCollider);gpu_buffer_update(&gpuPhysics.static_collider_buffer,gpuPhysics.static_collider_stage,size,0);bytes+=size;}
    else for(int c=0;c<staticCacheChunkCount;c++)if(gpuCacheChunks[c].upload){GpuCacheChunk *chunk=&gpuCacheChunks[c];size_t size=(size_t)chunk->count*sizeof(GpuStaticCollider),offset=(size_t)chunk->start*sizeof(GpuStaticCollider);gpu_buffer_update(&gpuPhysics.static_collider_buffer,gpuPhysics.static_collider_stage+chunk->start,size,offset);bytes+=size;}
    memset(gpuStaticPages,0,sizeof(gpuStaticPages));for(int c=0;c<staticCacheChunkCount;c++)gpuCacheChunks[c].upload=false;
    gpuStaticBytesUploaded+=bytes;return bytes;
}
