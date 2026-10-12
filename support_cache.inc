/* Static support topology. Node slots survive voxel-array compaction; identity
 * lookup is rebuilt once per refresh, not once per edited voxel. */
typedef struct {
    uint64_t identity;
    int index, component, count;
    int bounds[6], priorBounds[6];
    int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
    uint8_t mask;
    bool alive, seen, dirty, affected, edited, grounded;
} SupportNode;
static SupportNode *supportNodes;
static int supportNodeCount;
static int supportByIndex[MAX_VOXELS], supportIdentityHash[HASH_SIZE];
static int supportParents[MAX_VOXELS];
static bool supportComponentsDirty[MAX_VOXELS];
static bool supportCacheReady;
static uint64_t supportVisitedNodes, supportNeighborBuilds;

static void support_cache_reset(void) {
    free(supportNodes); supportNodes=NULL; supportNodeCount=0;supportCacheReady=false;
}
static unsigned support_identity_hash(uint64_t id) {return (unsigned)hashVoxelKey(id);}
static int support_find(uint64_t id) {
    unsigned h=support_identity_hash(id);
    for(int n=0;n<HASH_SIZE;n++,h=(h+1)&(HASH_SIZE-1)) {
        int s=supportIdentityHash[h];if(s<0)return -1;
        if(supportNodes[s].identity==id)return s;
    }return -1;
}
static int support_root(int s) {
    while(supportParents[s]!=s){supportParents[s]=supportParents[supportParents[s]];s=supportParents[s];}return s;
}
static void support_dirty_boundary(const int *b) {
    for(int axis=0;axis<3;axis++)for(int side=0;side<2;side++) {
        int u=(axis+1)%3,v=(axis+2)%3;
        for(int a=b[u*2];a<=b[u*2+1];a++)for(int c=b[v*2];c<=b[v*2+1];c++){
            int xyz[3];xyz[axis]=b[axis*2+side]+(side?1:-1);xyz[u]=a;xyz[v]=c;
            int index=table_get_static_only(xyz[0],xyz[1],xyz[2]);
            if(index>=0){int n=supportByIndex[index];if(n>=0){supportNodes[n].dirty=true;supportComponentsDirty[supportNodes[n].component]=true;}}
        }
    }
}
static bool support_cache_prepare(void) {
    if(!supportNodes){supportNodes=calloc(MAX_VOXELS,sizeof(*supportNodes));if(!supportNodes)return false;}
    memset(supportIdentityHash,0xff,sizeof(supportIdentityHash));
    memset(supportByIndex,0xff,sizeof(supportByIndex));
    memset(supportComponentsDirty,0,sizeof(supportComponentsDirty));
    for(int s=0;s<supportNodeCount;s++) {
        SupportNode *n=&supportNodes[s];n->seen=false;n->dirty=false;n->affected=false;n->edited=false;
        if(!n->alive)continue;
        unsigned h=support_identity_hash(n->identity);while(supportIdentityHash[h]>=0)h=(h+1)&(HASH_SIZE-1);
        supportIdentityHash[h]=s;
    }
    int freeCursor=0;
    for(int i=0;i<voxel_count;i++) {
        Voxel *v=&voxels[i];if(v->simulate)continue;
        int s=support_find(v->identity);bool added=s<0;
        if(added){
            while(freeCursor<supportNodeCount && supportNodes[freeCursor].alive)freeCursor++;
            s=freeCursor;
            if(s>=MAX_VOXELS){support_cache_reset();return support_cache_prepare();}
            if(s==supportNodeCount)supportNodeCount++;
            supportNodes[s]=(SupportNode){.identity=v->identity,.component=s,.alive=true};
        }
        SupportNode *n=&supportNodes[s];int bounds[6];
        voxel_grid_bounds(v,&bounds[0],&bounds[1],&bounds[2],&bounds[3],&bounds[4],&bounds[5]);
        n->index=i;n->seen=true;supportByIndex[i]=s;
        VoxelWorldBounds worldBounds;voxel_world_bounds(v,&worldBounds);bool grounded=worldBounds.miny<=GRID_EPSILON;
        if(added || grounded!=n->grounded || memcmp(bounds,n->bounds,sizeof(bounds))){
            n->dirty=n->edited=true;supportComponentsDirty[n->component]=true;
            memcpy(n->priorBounds,added?bounds:n->bounds,sizeof(bounds));memcpy(n->bounds,bounds,sizeof(bounds));n->grounded=grounded;
        }
    }
    // Removed nodes dirty their old component and cached neighbors before slots
    // can be reused on a later transaction.
    for(int s=0;s<supportNodeCount;s++) {
        SupportNode *n=&supportNodes[s];if(!n->alive || n->seen)continue;
        supportComponentsDirty[n->component]=true;
        for(int k=0;k<n->count;k++)supportNodes[n->neighbors[k]].dirty=true;
        memcpy(n->priorBounds,n->bounds,sizeof(n->bounds));n->edited=true;n->alive=false;
    }
    for(int n=0;n<supportNodeCount;n++)if(supportNodes[n].edited){support_dirty_boundary(supportNodes[n].priorBounds);if(supportNodes[n].alive)support_dirty_boundary(supportNodes[n].bounds);}
    // Discover the new boundary before rebuilding its reciprocal connections.
    for(int s=0;s<supportNodeCount;s++)if(supportNodes[s].alive && supportNodes[s].dirty) {
        int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
        int count=gather_static_face_neighbors(&voxels[supportNodes[s].index],neighbors,MAX_STATIC_COLLISION_NEIGHBORS);
        for(int k=0;k<count;k++){int t=supportByIndex[neighbors[k]];if(t>=0){supportNodes[t].affected=true;supportComponentsDirty[supportNodes[t].component]=true;}}
    }
    for(int s=0;s<supportNodeCount;s++) {
        SupportNode *n=&supportNodes[s];if(!n->alive)continue;
        if(n->dirty || n->affected || !supportCacheReady){
            int neighbors[MAX_STATIC_COLLISION_NEIGHBORS];
            n->count=gather_static_face_neighbors(&voxels[n->index],neighbors,MAX_STATIC_COLLISION_NEIGHBORS);
            for(int k=0;k<n->count;k++)n->neighbors[k]=supportByIndex[neighbors[k]];
            n->mask=compute_static_support_mask(&voxels[n->index]);supportNeighborBuilds++;
        }
        n->affected=!supportCacheReady || n->dirty || n->affected || supportComponentsDirty[n->component];
        supportParents[s]=s;
    }
    for(int s=0;s<supportNodeCount;s++)if(supportNodes[s].alive)
        for(int k=0;k<supportNodes[s].count;k++){int t=supportNodes[s].neighbors[k];if(t>=0 && supportNodes[t].alive){int a=support_root(s),b=support_root(t);supportParents[a]=b;}}
    memset(supportComponentsDirty,0,sizeof(supportComponentsDirty));
    for(int s=0;s<supportNodeCount;s++)if(supportNodes[s].alive && supportNodes[s].affected)supportComponentsDirty[support_root(s)]=true;
    for(int s=0;s<supportNodeCount;s++)if(supportNodes[s].alive){supportNodes[s].component=support_root(s);supportNodes[s].affected=supportComponentsDirty[supportNodes[s].component];}
    supportCacheReady=true;return true;
}
static int support_cached_neighbors(int index,int *out) {
    SupportNode *n=&supportNodes[supportByIndex[index]];
    for(int k=0;k<n->count;k++)out[k]=supportNodes[n->neighbors[k]].index;
    return n->count;
}
