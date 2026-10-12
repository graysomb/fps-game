static void cache_cell_insert(int index) {
    StaticSurfaceCell *c=&staticSurfaceCells[index];unsigned mask=(unsigned)staticSurfaceHashSize-1;
    unsigned h=static_surface_hash_coord(c->x,c->y,c->z,staticSurfaceHashSize);
    while(staticSurfaceHash[h]>=0)h=(h+1)&mask;
    staticSurfaceHash[h]=index;
}
static void cache_cell_remove(const StaticSurfaceCell *c) {
    unsigned mask=(unsigned)staticSurfaceHashSize-1,h=static_surface_hash_coord(c->x,c->y,c->z,staticSurfaceHashSize);
    while(staticSurfaceHash[h]>=0){StaticSurfaceCell *v=&staticSurfaceCells[staticSurfaceHash[h]];if(v->x==c->x && v->y==c->y && v->z==c->z)break;h=(h+1)&mask;}
    if(staticSurfaceHash[h]<0)return;
    unsigned hole=h;
    for(unsigned j=(h+1)&mask;staticSurfaceHash[j]>=0;j=(j+1)&mask){
        StaticSurfaceCell *v=&staticSurfaceCells[staticSurfaceHash[j]];
        unsigned home=static_surface_hash_coord(v->x,v->y,v->z,staticSurfaceHashSize);
        if(((j-home)&mask)>=((hole-home)&mask)){staticSurfaceHash[hole]=staticSurfaceHash[j];hole=j;}
    }
    staticSurfaceHash[hole]=-1;
}
static bool cache_build_patches(StaticCacheChunk *chunk) {
    int face_count=0,face_capacity=cacheFaceCapacity;StaticSurfaceFace *faces=cacheFaces;
    int ox=chunk->x*CHUNK_SIZE,oy=chunk->y*CHUNK_SIZE,oz=chunk->z*CHUNK_SIZE;
    Patch *saved=patches;int savedCount=patchCount,savedCapacity=patchCapacity;
    patches=cachePatchScratch;patchCount=0;patchCapacity=cachePatchCapacity;
    for(int v=0;v<voxel_count;v++) {
        Voxel *voxel=&voxels[v];if(voxel->simulate || voxel_is_fluid(voxel))continue;
        int x0=voxel->min_gx>ox?voxel->min_gx:ox,x1=voxel->max_gx<ox+CHUNK_MASK?voxel->max_gx:ox+CHUNK_MASK;
        int y0=voxel->min_gy>oy?voxel->min_gy:oy,y1=voxel->max_gy<oy+CHUNK_MASK?voxel->max_gy:oy+CHUNK_MASK;
        int z0=voxel->min_gz>oz?voxel->min_gz:oz,z1=voxel->max_gz<oz+CHUNK_MASK?voxel->max_gz:oz+CHUNK_MASK;
        for(int z=z0;z<=z1;z++)for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++) {
            if(table_get_static_only(x,y,z)!=v)continue;
            if(!is_static_solid_cell(x+1,y,z) && !append_static_surface_face(&faces,&face_count,&face_capacity,2,true,x,y,z))goto fail;
            if(!is_static_solid_cell(x-1,y,z) && !append_static_surface_face(&faces,&face_count,&face_capacity,2,false,x,y,z))goto fail;
            if(!is_static_solid_cell(x,y+1,z) && !append_static_surface_face(&faces,&face_count,&face_capacity,1,true,y,x,z))goto fail;
            if(!is_static_solid_cell(x,y-1,z) && !append_static_surface_face(&faces,&face_count,&face_capacity,1,false,y,x,z))goto fail;
            if(!is_static_solid_cell(x,y,z+1) && !append_static_surface_face(&faces,&face_count,&face_capacity,0,true,z,x,y))goto fail;
            if(!is_static_solid_cell(x,y,z-1) && !append_static_surface_face(&faces,&face_count,&face_capacity,0,false,z,x,y))goto fail;
        }
    }
    qsort(faces, (size_t)face_count, sizeof(*faces), compare_static_surface_faces);
    for (int index = 0; index < face_count; ++index) {
        StaticSurfaceFace *face = &faces[index];
        if (face->used) continue;
        int width = 1;
        for (;;) {
            int found = find_static_surface_face(faces, face_count, face->plane, face->positive,
                                                 face->layer, face->i + width, face->j);
            if (found < 0 || faces[found].used) break;
            ++width;
        }
        int height = 1;
        for (;;) {
            bool complete = true;
            for (int di = 0; di < width; ++di) {
                int found = find_static_surface_face(faces, face_count, face->plane, face->positive,
                                                     face->layer, face->i + di, face->j + height);
                if (found < 0 || faces[found].used) { complete = false; break; }
            }
            if (!complete) break;
            ++height;
        }
        for (int dj = 0; dj < height; ++dj) {
            for (int di = 0; di < width; ++di) {
                int found = find_static_surface_face(faces, face_count, face->plane, face->positive,
                                                     face->layer, face->i + di, face->j + dj);
                if (found >= 0) faces[found].used = true;
            }
        }
        if (!emit_patch(face->plane, face->layer, face->i, face->j,
                        width, height, face->positive)) goto fail;
    }

    cacheFaces=faces;cacheFaceCapacity=face_capacity;
    cachePatchScratch=patches;cachePatchCapacity=patchCapacity;int count=patchCount;
    patches=saved;patchCount=savedCount;patchCapacity=savedCapacity;
    if(chunk->patchCapacity)memset(patches+chunk->patchStart,0,(size_t)chunk->patchCapacity*sizeof(Patch));
    if(!cache_reserve(&staticPatchArena,&chunk->patchStart,&chunk->patchCapacity,count) ||
       !cache_array((void**)&patches,&patchCapacity,staticPatchArena.high,sizeof(Patch)))return false;
    if(count)memcpy(patches+chunk->patchStart,cachePatchScratch,(size_t)count*sizeof(Patch));
    chunk->patchCount=count;patchCount=staticPatchArena.high;
    ++staticCacheChunksBuilt;staticCachePatchesBuilt+=count;return true;
fail:
    cacheFaces=faces;cacheFaceCapacity=face_capacity;cachePatchScratch=patches;cachePatchCapacity=patchCapacity;
    patches=saved;patchCount=savedCount;patchCapacity=savedCapacity;return false;
}
static bool cache_build_cells(StaticCacheChunk *chunk) {
    int count=0;
    for(int c=0;c<staticCacheChunkCount;c++) {
        StaticCacheChunk *owner=&staticCacheChunks[c];
        if(abs(owner->x-chunk->x)>1 || abs(owner->y-chunk->y)>1 || abs(owner->z-chunk->z)>1)continue;
        for(int p=owner->patchStart;p<owner->patchStart+owner->patchCount;p++) {
            Patch *pt=&patches[p];
            for(int j=0;j<pt->dj;j++)for(int i=0;i<pt->di;i++) {
                int a=pt->i0+i,b=pt->j0+j,layer=pt->layer+(pt->positive?1:-1);
                int x=pt->plane==2?layer:a,y=pt->plane==0?b:pt->plane==1?layer:a,z=pt->plane==0?layer:b;
                if(grid_to_chunk_coord(x)!=chunk->x || grid_to_chunk_coord(y)!=chunk->y || grid_to_chunk_coord(z)!=chunk->z)continue;
                if(!append_static_surface_ref(&cacheRefs,&count,&cacheRefCapacity,x,y,z,p))return false;
            }
        }
    }
    qsort(cacheRefs,(size_t)count,sizeof(*cacheRefs),compare_static_surface_refs);
    for(int k=0;k<chunk->cellCount;k++)cache_cell_remove(&staticSurfaceCells[chunk->cellStart+k]);
    if(chunk->cellCapacity)memset(staticSurfaceCells+chunk->cellStart,0,(size_t)chunk->cellCapacity*sizeof(StaticSurfaceCell));
    if(!cache_reserve(&staticCellArena,&chunk->cellStart,&chunk->cellCapacity,count) ||
       !cache_reserve(&staticRefArena,&chunk->refStart,&chunk->refCapacity,count) ||
       !cache_array((void**)&staticSurfaceCells,&staticSurfaceCellCapacity,staticCellArena.high,sizeof(StaticSurfaceCell)) ||
       !cache_array((void**)&staticSurfacePatchRefs,&staticSurfacePatchRefCapacity,staticRefArena.high,sizeof(int)))return false;
    chunk->refCount=count;chunk->cellCount=0;
    for(int k=0;k<count;) {
        int begin=k,x=cacheRefs[k].x,y=cacheRefs[k].y,z=cacheRefs[k].z;
        do{staticSurfacePatchRefs[chunk->refStart+k]=cacheRefs[k].patch;k++;}while(k<count && cacheRefs[k].x==x && cacheRefs[k].y==y && cacheRefs[k].z==z);
        int slot=chunk->cellStart+chunk->cellCount++;
        staticSurfaceCells[slot]=(StaticSurfaceCell){x,y,z,chunk->refStart+begin,k-begin};
    }
    int live=0;for(int c=0;c<staticCacheChunkCount;c++)live+=staticCacheChunks[c].cellCount;
    if(live>staticSurfaceHashSize/2) {
        int size=staticSurfaceHashSize;while(size<live*2)size*=2;
        int *next=realloc(staticSurfaceHash,(size_t)size*sizeof(int));if(!next)return false;
        staticSurfaceHash=next;staticSurfaceHashSize=size;memset(next,0xff,(size_t)size*sizeof(int));
        for(int c=0;c<staticCacheChunkCount;c++)for(int k=0;k<staticCacheChunks[c].cellCount;k++)cache_cell_insert(staticCacheChunks[c].cellStart+k);
    }else for(int k=0;k<chunk->cellCount;k++)cache_cell_insert(chunk->cellStart+k);
    staticSurfaceCellCount=staticCellArena.high;staticSurfacePatchRefCount=staticRefArena.high;
    chunk->revision=++staticCacheRevision;chunk->dataDirty=false;return true;
}
static bool cache_rebuild_chunks(void) {
    if(!staticCacheReady) {
        static_cache_reset();staticSurfaceHashSize=HASH_SIZE;
        staticSurfaceHash=malloc((size_t)staticSurfaceHashSize*sizeof(int));if(!staticSurfaceHash)return false;
        memset(staticSurfaceHash,0xff,(size_t)staticSurfaceHashSize*sizeof(int));
        for(int i=0;i<voxel_count;i++)if(!voxels[i].simulate) {
            Voxel *v=&voxels[i];
            for(int z=grid_to_chunk_coord(v->min_gz);z<=grid_to_chunk_coord(v->max_gz);z++)
            for(int y=grid_to_chunk_coord(v->min_gy);y<=grid_to_chunk_coord(v->max_gy);y++)
            for(int x=grid_to_chunk_coord(v->min_gx);x<=grid_to_chunk_coord(v->max_gx);x++)if(cache_chunk_find(x,y,z,true)<0)return false;
        }
    }
    if(staticCacheOverflow)return false;
    int owners=staticCacheChunkCount;
    for(int c=0;c<owners;c++)if(staticCacheChunks[c].patchesDirty) {
        StaticCacheChunk *chunk=&staticCacheChunks[c];if(!cache_build_patches(chunk))return false;
        chunk->patchesDirty=false;
        // External-face references can land in the neighboring air chunk.
        for(int z=-1;z<=1;z++)for(int y=-1;y<=1;y++)for(int x=-1;x<=1;x++) {
            if(abs(x)+abs(y)+abs(z)>1)continue;
            int old=staticCacheChunkCount,n=cache_chunk_find(chunk->x+x,chunk->y+y,chunk->z+z,true);if(n<0)return false;
            if(n==old)staticCacheChunks[n].patchesDirty=false;
            staticCacheChunks[n].dataDirty=true;
        }
    }
    for(int c=0;c<staticCacheChunkCount;c++)if(staticCacheChunks[c].dataDirty && !cache_build_cells(&staticCacheChunks[c]))return false;
    staticSurfaceOccupiedCellCount=0;
    for(int i=0;i<voxel_count;i++)if(!voxels[i].simulate && !voxel_is_fluid(&voxels[i])) {
        Voxel *v=&voxels[i];staticSurfaceOccupiedCellCount+=(v->max_gx-v->min_gx+1)*(v->max_gy-v->min_gy+1)*(v->max_gz-v->min_gz+1);
    }
    staticCacheReady=true;staticSurfaceGeneration=staticHashGeneration;return true;
}
static bool rebuild_static_surface_cache(void) {
    double started=GetTime();
    bool result;
    if(getenv("FPS_COLLISION_CACHE_FULL")){static_cache_reset();result=rebuild_static_surface_cache_legacy();}
    else {
        bool dirty=false;for(int c=0;c<staticCacheChunkCount;c++)dirty|=staticCacheChunks[c].patchesDirty;
        // Unknown bulk edits (load/reset tools) still invalidate the full cache.
        if(staticCacheReady && (!dirty || staticCacheForceRebuild))static_cache_reset();
        result=cache_rebuild_chunks();
        if(!result){static_cache_reset();result=rebuild_static_surface_cache_legacy();}
    }
    staticCacheBuildMs+=(GetTime()-started)*1000;return result;
}
