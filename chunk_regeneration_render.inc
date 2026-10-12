/* Preview meshes are cached once per job; only the fill uniform changes. */
#define REGEN_MESH_PATCHES 8192
#define REGEN_MAX_MESHES 24
static Mesh regenerationMeshes[REGEN_MAX_MESHES];
static int regenerationMeshCount;
static Material regenerationMaterial;
static bool regenerationMaterialReady,regenerationCustomShader;
static int regenerationFillLoc;
static void regeneration_release_mesh(void) {
    for(int i=0;i<regenerationMeshCount;i++)UnloadMesh(regenerationMeshes[i]);
    regenerationMeshCount=0;
}
static void regeneration_shutdown(void) {
    regeneration_reset();
    if(regenerationMaterialReady) {
        if(regenerationCustomShader)UnloadShader(regenerationMaterial.shader);
        regenerationMaterial.shader=(Shader){0};UnloadMaterial(regenerationMaterial);
    }
    regenerationMaterialReady=regenerationCustomShader=false;
}
static bool regeneration_has(int x,int y,int z) {
    return x>=0 && x<CHUNK_SIZE && y>=0 && y<CHUNK_SIZE && z>=0 && z<CHUNK_SIZE &&
           regeneration.cells[regeneration_cell_index(x,y,z)].present;
}
static bool regeneration_build_mesh(void) {
    regeneration_release_mesh();
    uint32_t faces[6][CHUNK_SIZE][CHUNK_SIZE]={0};
    for(int z=0;z<CHUNK_SIZE;z++)for(int y=0;y<CHUNK_SIZE;y++)for(int x=0;x<CHUNK_SIZE;x++)if(regeneration_has(x,y,z)) {
        if(!regeneration_has(x,y,z-1))faces[0][z][y]|=1u<<x;
        if(!regeneration_has(x,y,z+1))faces[1][z][y]|=1u<<x;
        if(!regeneration_has(x,y-1,z))faces[2][y][z]|=1u<<x;
        if(!regeneration_has(x,y+1,z))faces[3][y][z]|=1u<<x;
        if(!regeneration_has(x-1,y,z))faces[4][x][z]|=1u<<y;
        if(!regeneration_has(x+1,y,z))faces[5][x][z]|=1u<<y;
    }
    Patch *buffer=malloc(REGEN_MESH_PATCHES*sizeof(*buffer));if(!buffer)return false;
    int count=0,ox=regeneration.cx*CHUNK_SIZE,oy=regeneration.cy*CHUNK_SIZE,oz=regeneration.cz*CHUNK_SIZE;
    for(int face=0;face<6;face++)for(int layer=0;layer<CHUNK_SIZE;layer++)for(int j=0;j<CHUNK_SIZE;j++) {
        uint32_t *rows=faces[face][layer];
        while(rows[j]) {
            int i=0,width=1;while(!(rows[j]&(1u<<i)))i++;
            while(i+width<CHUNK_SIZE && (rows[j]&(1u<<(i+width))))width++;
            uint32_t run=(UINT32_MAX>>(CHUNK_SIZE-width))<<i;int height=1;
            while(j+height<CHUNK_SIZE && (rows[j+height]&run)==run)height++;
            for(int row=j;row<j+height;row++)rows[row]&=~run;
            int plane=face/2;
            int cell=plane==0?regeneration_cell_index(i,j,layer):plane==1?regeneration_cell_index(i,layer,j):regeneration_cell_index(layer,i,j);
            Color color=regeneration.cells[cell].color;color.a=105;
            buffer[count++]=(Patch){.plane=plane,.positive=(face&1)!=0,
                .layer=layer+(plane==0?oz:plane==1?oy:ox),.i0=i+(plane==2?oy:ox),.j0=j+(plane==0?oy:oz),
                .di=width,.dj=height,.col=color};
            if(count==REGEN_MESH_PATCHES){regenerationMeshes[regenerationMeshCount++]=build_mesh_from_patches(buffer,count);count=0;}
        }
    }
    if(count)regenerationMeshes[regenerationMeshCount++]=build_mesh_from_patches(buffer,count);
    free(buffer);regeneration.meshDirty=false;return true;
}
static void regeneration_draw(void) {
    if(!regeneration.active)return;
    if(regeneration.meshDirty && !regeneration_build_mesh())return;
    if(!regenerationMaterialReady) {
        regenerationMaterial=LoadMaterialDefault();
        Shader shader={0};
        if(!getenv("FPS_REGEN_SHADER_FALLBACK"))shader=LoadShader("shaders/regeneration.vert","shaders/regeneration.frag");
        if(shader.id && shader.id!=rlGetShaderIdDefault()) {
            regenerationMaterial.shader=shader;regenerationCustomShader=true;
            regenerationFillLoc=GetShaderLocation(shader,"fillHeight");
        }
        regenerationMaterialReady=true;
    }
    float progress=netTransport.role==NET_ROLE_CLIENT?(float)((GetTime()-regeneration.remoteStart)/.6):
        regeneration.elapsed/.6f;
    progress=Clamp(progress,0,1);
    float low=INFINITY,high=-INFINITY;
    for(int n=0;n<REGEN_CELL_COUNT;n++)if(regeneration.cells[n].present) {
        float y=(regeneration.cy*CHUNK_SIZE+(n/CHUNK_SIZE)%CHUNK_SIZE)*VOXEL_SIZE;
        low=fminf(low,y);high=fmaxf(high,y+VOXEL_SIZE);
    }
    float fill=low+(high-low)*progress;
    if(regenerationCustomShader)SetShaderValue(regenerationMaterial.shader,regenerationFillLoc,&fill,SHADER_UNIFORM_FLOAT);
    rlDrawRenderBatchActive();rlDisableDepthMask();rlDisableBackfaceCulling();
    BeginBlendMode(BLEND_ALPHA);
    for(int i=0;i<regenerationMeshCount;i++)DrawMesh(regenerationMeshes[i],regenerationMaterial,MatrixIdentity());
    EndBlendMode();rlDrawRenderBatchActive();rlEnableDepthMask();rlEnableBackfaceCulling();
}
