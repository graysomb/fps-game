/* Persistent ranges shared by CPU patches/cells and GPU collider staging. */
#define STATIC_CACHE_CHUNKS 2048
#define STATIC_CACHE_SPANS 4096
typedef struct {int start,count;} CacheSpan;
typedef struct {CacheSpan free[STATIC_CACHE_SPANS];int count,high;} CacheArena;
static bool cache_reserve(CacheArena *a,int *start,int *capacity,int need) {
    if(need<=*capacity)return true;
    int size=16;while(size<need){if(size>INT_MAX/2)return false;size*=2;}
    int slot=-1;
    for(int i=0;i<a->count;i++)if(a->free[i].count>=size){slot=i;break;}
    int next;
    if(slot>=0){next=a->free[slot].start;a->free[slot].start+=size;a->free[slot].count-=size;if(!a->free[slot].count)a->free[slot]=a->free[--a->count];}
    else {if(a->high>INT_MAX-size)return false;next=a->high;a->high+=size;}
    if(*capacity){
        if(a->count>=STATIC_CACHE_SPANS)return false;
        a->free[a->count++]=(CacheSpan){*start,*capacity};
        // Coalesce ranges so repeated growth does not accumulate tiny holes.
        for(int i=0;i<a->count;i++)for(int j=i+1;j<a->count;) {
            CacheSpan *x=&a->free[i],*y=&a->free[j];
            if(x->start+x->count==y->start || y->start+y->count==x->start){if(y->start<x->start)x->start=y->start;x->count+=y->count;a->free[j]=a->free[--a->count];j=i+1;}else j++;
        }
    }
    *start=next;*capacity=size;return true;
}
static bool cache_array(void **p,int *capacity,int need,size_t size) {
    if(need<=*capacity)return true;
    int cap=*capacity?*capacity:16;while(cap<need){if(cap>INT_MAX/2)return false;cap*=2;}
    void *next=realloc(*p,(size_t)cap*size);if(!next)return false;
    memset((char*)next+(size_t)*capacity*size,0,(size_t)(cap-*capacity)*size);*p=next;*capacity=cap;return true;
}
typedef struct {
    int x,y,z;
    int patchStart,patchCapacity,patchCount;
    int cellStart,cellCapacity,cellCount,refStart,refCapacity,refCount;
    bool patchesDirty,dataDirty;
    uint64_t revision;
} StaticCacheChunk;
static StaticCacheChunk staticCacheChunks[STATIC_CACHE_CHUNKS];
static int staticCacheChunkCount;
static CacheArena staticPatchArena,staticCellArena,staticRefArena;
static bool staticCacheReady, staticCacheOverflow,staticCacheSuppressEdits,staticCacheForceRebuild;
static uint64_t staticCacheMutationSerial,staticCacheMarkedSerial;
static uint64_t staticCacheEpoch=1,staticCacheRevision;
static uint64_t staticCacheChunksBuilt,staticCachePatchesBuilt;
static double staticCacheBuildMs;
static StaticSurfaceFace *cacheFaces;static int cacheFaceCapacity;
static StaticSurfaceRef *cacheRefs;static int cacheRefCapacity;
static Patch *cachePatchScratch;static int cachePatchCapacity;
static int cache_chunk_find(int x,int y,int z,bool create) {
    for(int i=0;i<staticCacheChunkCount;i++)if(staticCacheChunks[i].x==x && staticCacheChunks[i].y==y && staticCacheChunks[i].z==z)return i;
    if(!create)return -1;
    if(staticCacheChunkCount==STATIC_CACHE_CHUNKS){staticCacheOverflow=true;return -1;}
    int i=staticCacheChunkCount++;staticCacheChunks[i]=(StaticCacheChunk){.x=x,.y=y,.z=z,.patchesDirty=true,.dataDirty=true};return i;
}
static void static_cache_dirty_bounds_local(int minx,int maxx,int miny,int maxy,int minz,int maxz) {
    if(!staticCacheReady || staticCacheSuppressEdits)return;
    ++staticCacheMutationSerial;
    for(int z=grid_to_chunk_coord(minz-1);z<=grid_to_chunk_coord(maxz+1);z++)
    for(int y=grid_to_chunk_coord(miny-1);y<=grid_to_chunk_coord(maxy+1);y++)
    for(int x=grid_to_chunk_coord(minx-1);x<=grid_to_chunk_coord(maxx+1);x++){
        int c=cache_chunk_find(x,y,z,true);if(c>=0)staticCacheChunks[c].patchesDirty=true;
    }
}
static void static_cache_dirty_bounds(int minx,int maxx,int miny,int maxy,int minz,int maxz) {
    static_cache_dirty_bounds_local(minx,maxx,miny,maxy,minz,maxz);
    if(!worldPeriodic || !staticCacheReady || staticCacheSuppressEdits)return;
    int half=(int)lroundf(FLOOR_SIZE/VOXEL_SIZE),period=2*half;
    // Only seam-touching edits have counterpart chunks. Avoid creating empty remote chunks.
    if(minx<=-half)static_cache_dirty_bounds_local(minx+period,maxx+period,miny,maxy,minz,maxz);
    if(maxx>=half-1)static_cache_dirty_bounds_local(minx-period,maxx-period,miny,maxy,minz,maxz);
    if(minz<=-half)static_cache_dirty_bounds_local(minx,maxx,miny,maxy,minz+period,maxz+period);
    if(maxz>=half-1)static_cache_dirty_bounds_local(minx,maxx,miny,maxy,minz-period,maxz-period);
}
static void static_cache_reset(void) {
    staticCacheReady=false;staticCacheOverflow=false;staticCacheForceRebuild=false;staticCacheChunkCount=0;
    staticCacheMutationSerial=staticCacheMarkedSerial=0;
    staticPatchArena=(CacheArena){0};staticCellArena=(CacheArena){0};staticRefArena=(CacheArena){0};
    free(patches);patches=NULL;patchCount=patchCapacity=0;
    free(staticSurfaceCells);staticSurfaceCells=NULL;staticSurfaceCellCount=staticSurfaceCellCapacity=0;
    free(staticSurfacePatchRefs);staticSurfacePatchRefs=NULL;staticSurfacePatchRefCount=staticSurfacePatchRefCapacity=0;
    free(staticSurfaceHash);staticSurfaceHash=NULL;staticSurfaceHashSize=0;
    free(cacheFaces);cacheFaces=NULL;cacheFaceCapacity=0;
    free(cacheRefs);cacheRefs=NULL;cacheRefCapacity=0;
    free(cachePatchScratch);cachePatchScratch=NULL;cachePatchCapacity=0;
    staticSurfaceOccupiedCellCount=0;staticSurfaceGeneration=0;++staticCacheEpoch;
}
