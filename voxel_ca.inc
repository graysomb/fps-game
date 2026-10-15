/* Shared, purely cosmetic CA. Update only outside a viewport/3D render pass. */
#define VOXEL_CA_SIZE 256
#define VOXEL_CA_STEP (1.0/12.0)
static RenderTexture2D voxelCaTargets[2];
static Shader voxelCaStepShader;
static bool voxelCaReady,voxelCaAttempted;
static int voxelCaCurrent,voxelCaStepLoc,voxelCaSourceLoc;
static unsigned voxelCaGeneration;
static double voxelCaLastTime,voxelCaAccumulator;
static float voxelCaBlend=1;
typedef struct {unsigned id;int enabled,blend;} VoxelCaBinding;
static VoxelCaBinding voxelCaBindings[2];
static uint32_t voxel_ca_hash(uint32_t x) {
    x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);
}
static void voxel_ca_shutdown(void) {
    for(int i=0;i<2;i++){if(voxelCaTargets[i].id)UnloadRenderTexture(voxelCaTargets[i]);voxelCaTargets[i]=(RenderTexture2D){0};}
    if(voxelCaStepShader.id && voxelCaStepShader.id!=rlGetShaderIdDefault())UnloadShader(voxelCaStepShader);
    voxelCaStepShader=(Shader){0};voxelCaReady=voxelCaAttempted=false;
    memset(voxelCaBindings,0,sizeof(voxelCaBindings));
    voxelCaCurrent=0;voxelCaGeneration=0;voxelCaAccumulator=0;voxelCaBlend=1;
}
static bool voxel_ca_init(void) {
    if(voxelCaAttempted || !voxelCaEnabled)return voxelCaReady;
    voxelCaAttempted=true;
    if(getenv("FPS_VOXEL_CA_FALLBACK"))return false;
    voxelCaStepShader=LoadShader(NULL,"shaders/voxel_ca_step.frag");
    if(!voxelCaStepShader.id || voxelCaStepShader.id==rlGetShaderIdDefault())goto failed;
    for(int i=0;i<2;i++) {
        voxelCaTargets[i]=LoadRenderTexture(VOXEL_CA_SIZE,VOXEL_CA_SIZE);
        if(!IsRenderTextureValid(voxelCaTargets[i]))goto failed;
        SetTextureFilter(voxelCaTargets[i].texture,TEXTURE_FILTER_POINT);
    }
    {
        Image seed=GenImageColor(VOXEL_CA_SIZE,VOXEL_CA_SIZE,BLACK);
        if(!seed.data)goto failed;
        Color *cells=seed.data;
        for(unsigned i=0;i<VOXEL_CA_SIZE*VOXEL_CA_SIZE;i++) {
            unsigned char v=voxel_ca_hash(i+0x1327u)%100<30?255:0;
            cells[i]=(Color){v,v,0,255};
        }
        for(int i=0;i<2;i++)UpdateTexture(voxelCaTargets[i].texture,seed.data);
        UnloadImage(seed);
    }
    voxelCaStepLoc=GetShaderLocation(voxelCaStepShader,"generation");
    voxelCaSourceLoc=GetShaderLocation(voxelCaStepShader,"previousState");
    if(voxelCaStepLoc<0 || voxelCaSourceLoc<0)goto failed;
    voxelCaLastTime=GetTime();voxelCaReady=true;return true;
failed:
    voxel_ca_shutdown();voxelCaAttempted=true;
    TraceLog(LOG_WARNING,"Voxel CA unavailable; using plain voxel surfaces");return false;
}
static void voxel_ca_step(void) {
    if(!voxelCaReady)return;
    int next=1-voxelCaCurrent,generation=(int)voxelCaGeneration;
    BeginTextureMode(voxelCaTargets[next]);
    BeginShaderMode(voxelCaStepShader);
    SetShaderValue(voxelCaStepShader,voxelCaStepLoc,&generation,SHADER_UNIFORM_INT);
    SetShaderValueTexture(voxelCaStepShader,voxelCaSourceLoc,voxelCaTargets[voxelCaCurrent].texture);
    DrawRectangle(0,0,VOXEL_CA_SIZE,VOXEL_CA_SIZE,WHITE);
    EndShaderMode();EndTextureMode();
    voxelCaCurrent=next;voxelCaGeneration++;
}
static void voxel_ca_update(double now) {
    if(!voxelCaReady || !voxelCaEnabled){voxelCaLastTime=now;return;}
    double dt=now-voxelCaLastTime;voxelCaLastTime=now;
    if(dt<0)dt=0;
    // Discard excess catch-up work after a stall, rather than building a backlog.
    voxelCaAccumulator+=fmin(dt,2*VOXEL_CA_STEP);
    for(int n=0;n<2 && voxelCaAccumulator>=VOXEL_CA_STEP;n++) {
        voxel_ca_step();voxelCaAccumulator-=VOXEL_CA_STEP;
    }
    voxelCaBlend=(float)(voxelCaAccumulator/VOXEL_CA_STEP);
}
static Texture2D voxel_ca_texture(void) {
    return voxelCaReady && voxelCaEnabled?voxelCaTargets[voxelCaCurrent].texture:
        (Texture2D){.id=rlGetTextureIdDefault(),.width=1,.height=1,.mipmaps=1,.format=PIXELFORMAT_UNCOMPRESSED_R8G8B8A8};
}
static void voxel_ca_bind(Shader shader) {
    VoxelCaBinding *b=NULL;
    for(int i=0;i<2;i++)if(voxelCaBindings[i].id==shader.id){b=&voxelCaBindings[i];break;}
    if(!b) {
        b=&voxelCaBindings[voxelCaBindings[0].id?1:0];
        *b=(VoxelCaBinding){shader.id,GetShaderLocation(shader,"caEnabled"),GetShaderLocation(shader,"caBlend")};
    }
    if(b->enabled<0 || b->blend<0)voxelCaEnabled=false;
    int enabled=voxelCaReady && voxelCaEnabled;
    SetShaderValue(shader,b->enabled,&enabled,SHADER_UNIFORM_INT);
    SetShaderValue(shader,b->blend,&voxelCaBlend,SHADER_UNIFORM_FLOAT);
}
