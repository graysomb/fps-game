/* Shared post-process and bounded periodic camera images. No simulation copies. */
static Shader periodicFogShader,nightSkyShader;
static Texture2D generated_space_texture(void);
static RenderTexture2D periodicFogScratch,periodicIntroScene;
static const Color periodicSky={151,170,181,255};
static Camera3D periodicViewCamera;
static bool periodicViewActive;
static float periodicAspect=16.f/9.f;
static BoundingBox periodicTileBounds;
static void periodic_tile_bounds(void) {
    periodicTileBounds=(BoundingBox){{-FLOOR_SIZE-2,-2,-FLOOR_SIZE-2},{FLOOR_SIZE+2,2*FLOOR_SIZE,FLOOR_SIZE+2}};
    ensure_dynamic_voxels();
    for(int d=0;d<dynamic_voxel_count;d++)for(int c=0;c<8;c++) {
        Particle *p=voxels[dynamic_voxels[d]].particles[c];if(!p)continue;
        periodicTileBounds.min=Vector3Min(periodicTileBounds.min,p->pos);
        periodicTileBounds.max=Vector3Max(periodicTileBounds.max,p->pos);
    }
}
static bool periodic_visible_box(Camera3D camera,Vector3 center,Vector3 extent) {
    Vector3 f=Vector3Normalize(Vector3Subtract(camera.target,camera.position));
    Vector3 r=Vector3Normalize(Vector3CrossProduct(f,camera.up)),u=Vector3CrossProduct(r,f);
    float vertical=tanf(camera.fovy*DEG2RAD*.5f);
    // Use the actual viewport, including split-screen render targets.
    float aspect=periodicAspect;
    Vector3 planes[5]={f,Vector3Add(Vector3Scale(f,vertical*aspect),r),Vector3Subtract(Vector3Scale(f,vertical*aspect),r),Vector3Add(Vector3Scale(f,vertical),u),Vector3Subtract(Vector3Scale(f,vertical),u)};
    Vector3 d=Vector3Subtract(center,camera.position);
    for(int i=0;i<5;i++)if(Vector3DotProduct(planes[i],d)+fabsf(planes[i].x)*extent.x+fabsf(planes[i].y)*extent.y+fabsf(planes[i].z)*extent.z<0)return false;
    return true;
}
static bool periodic_visible(Vector3 center,float radius) {
    if(!worldPeriodic || !periodicViewActive)return true;
    return Vector3DistanceSqr(center,periodicViewCamera.position)<(56+radius)*(56+radius) &&
        periodic_visible_box(periodicViewCamera,center,(Vector3){radius,radius,radius});
}
static RenderTexture2D world_render_target(int w,int h) {
    RenderTexture2D t={0};t.id=rlLoadFramebuffer();
    t.texture=(Texture2D){rlLoadTexture(NULL,w,h,PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,1),w,h,1,PIXELFORMAT_UNCOMPRESSED_R8G8B8A8};
    t.depth=(Texture2D){rlLoadTextureDepth(w,h,false),w,h,1,19};
    rlFramebufferAttach(t.id,t.texture.id,RL_ATTACHMENT_COLOR_CHANNEL0,RL_ATTACHMENT_TEXTURE2D,0);
    rlFramebufferAttach(t.id,t.depth.id,RL_ATTACHMENT_DEPTH,RL_ATTACHMENT_TEXTURE2D,0);
    if(!rlFramebufferComplete(t.id)){UnloadRenderTexture(t);t=LoadRenderTexture(w,h);t.depth.mipmaps=0;}
    return t;
}
static void periodic_visuals_shutdown(void) {
    if(periodicIntroScene.id)UnloadRenderTexture(periodicIntroScene);
    periodicIntroScene=(RenderTexture2D){0};
    if(periodicFogScratch.id)UnloadRenderTexture(periodicFogScratch);
    if(periodicFogShader.id && periodicFogShader.id!=rlGetShaderIdDefault())UnloadShader(periodicFogShader);
    if(nightSkyShader.id && nightSkyShader.id!=rlGetShaderIdDefault())UnloadShader(nightSkyShader);
    nightSkyShader=(Shader){0};
    periodicFogScratch=(RenderTexture2D){0};periodicFogShader=(Shader){0};
}
static void periodic_fullscreen(Texture2D t,int w,int h) {
    DrawTexturePro(t,(Rectangle){0,0,(float)t.width,-(float)t.height},(Rectangle){0,0,(float)w,(float)h},(Vector2){0},0,WHITE);
}
// Called outside 3D mode, with scene bound. Leaves scene bound and depth intact.
static void periodic_fog(RenderTexture2D scene,Camera3D cam,bool emission,Texture2D depth) {
    periodicViewActive=false;
    if(!worldPeriodic || !depth.mipmaps)return;
    if(!periodicFogShader.id)periodicFogShader=LoadShader(NULL,"shaders/world_fog.frag");
    if(periodicFogShader.id==rlGetShaderIdDefault())return;
    int w=scene.texture.width,h=scene.texture.height;
    if(periodicFogScratch.texture.width!=w || periodicFogScratch.texture.height!=h) {
        if(periodicFogScratch.id)UnloadRenderTexture(periodicFogScratch);
        periodicFogScratch=LoadRenderTexture(w,h);
    }
    EndTextureMode();BeginTextureMode(periodicFogScratch);ClearBackground(BLACK);
    float aspect=(float)w/h,tanHalf=tanf(cam.fovy*DEG2RAD*.5f);
    Vector4 settings={aspect,tanHalf,(float)rlGetCullDistanceNear(),(float)rlGetCullDistanceFar()};
    Vector3 color=emission?(Vector3){0}:(Vector3){periodicSky.r/255.f,periodicSky.g/255.f,periodicSky.b/255.f};
    SetShaderValue(periodicFogShader,GetShaderLocation(periodicFogShader,"settings"),&settings,SHADER_UNIFORM_VEC4);
    SetShaderValue(periodicFogShader,GetShaderLocation(periodicFogShader,"fogColor"),&color,SHADER_UNIFORM_VEC3);
    BeginShaderMode(periodicFogShader);
    SetShaderValueTexture(periodicFogShader,GetShaderLocation(periodicFogShader,"worldDepth"),depth);
    periodic_fullscreen(scene.texture,w,h);EndShaderMode();EndTextureMode();
    BeginTextureMode(scene);periodic_fullscreen(periodicFogScratch.texture,w,h);rlDrawRenderBatchActive();
}
static bool periodic_camera(Camera3D camera,int image,Camera3D *result) {
    if(!worldPeriodic){if(image)return false;*result=camera;return true;}
    if(image==0)periodic_tile_bounds();
    float period=2*FLOOR_SIZE;
    Vector3 origin=Vector3Subtract(camera.position,world_wrap(camera.position));
    Vector3 shift={origin.x+(image%3-1)*period,0,origin.z+(image/3-1)*period};
    *result=camera;result->position=Vector3Subtract(camera.position,shift);result->target=Vector3Subtract(camera.target,shift);
    Vector3 nearest=Vector3Clamp(result->position,periodicTileBounds.min,periodicTileBounds.max);
    Vector3 center=Vector3Scale(Vector3Add(periodicTileBounds.min,periodicTileBounds.max),.5f);
    Vector3 extent=Vector3Scale(Vector3Subtract(periodicTileBounds.max,periodicTileBounds.min),.5f);
    bool visible=Vector3DistanceSqr(nearest,result->position)<56*56 && periodic_visible_box(*result,center,extent);
    if(visible){periodicViewCamera=*result;periodicViewActive=true;}
    return visible;
}
// Sample once per rendered frame; every split-screen view shares one sky event.
static double nightSkyTime;
static float shooting_star_random(unsigned slot,unsigned salt) {
    uint32_t x=slot*747796405u+salt*2891336453u+277803737u;
    x=((x>>((x>>28u)+4u))^x)*277803737u;x=(x>>22u)^x;
    return (x&0xffffffu)/16777216.f;
}
static void shooting_star_sample(double time,Vector4 *trail,float *glow) {
    *trail=(Vector4){0};*glow=0;
    if(time<0)return;
    unsigned slot=(unsigned)(time/5.0);
    float age=(float)(time-slot*5.0)-(2.7f+.7f*shooting_star_random(slot,1));
    if(age<0 || age>1.4f)return;
    float angle=shooting_star_random(slot,2)*2*PI;
    float radius=.08f+.12f*shooting_star_random(slot,3);
    float direction=shooting_star_random(slot,4)*2*PI;
    float dx=cosf(direction),dy=sinf(direction),travel=age*.28f;
    float x=cosf(angle)*radius+dx*(travel-.196f),y=sinf(angle)*radius+dy*(travel-.196f);
    float length=fminf(.14f,travel);
    *trail=(Vector4){x,y,x-dx*length,y-dy*length};
    *glow=clampf(age/.08f,0,1)*clampf((1.4f-age)/.35f,0,1);
}
static void periodic_skyline(Camera3D camera) {
    // One sky per viewport, before periodic world images. Direction only: no
    // parallax, tile repetition or discontinuity when the player wraps.
    if(!nightSkyShader.id)nightSkyShader=LoadShader(NULL,"shaders/night_sky.frag");
    if(nightSkyShader.id!=rlGetShaderIdDefault() && nightSkyShader.id) {
        int w=rlGetFramebufferWidth(),h=rlGetFramebufferHeight();
        Vector3 forward=Vector3Normalize(Vector3Subtract(camera.target,camera.position));
        Vector3 right=Vector3Normalize(Vector3CrossProduct(forward,camera.up));
        Vector3 up=Vector3CrossProduct(right,forward);
        Vector4 viewport={(float)w,(float)h,tanf(camera.fovy*DEG2RAD*.5f),(float)w/h};
        Vector3 haze={periodicSky.r/255.f,periodicSky.g/255.f,periodicSky.b/255.f};
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewForward"),&forward,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewRight"),&right,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewUp"),&up,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewport"),&viewport,SHADER_UNIFORM_VEC4);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"hazeColor"),&haze,SHADER_UNIFORM_VEC3);
        Vector4 trail;float glow;shooting_star_sample(nightSkyTime,&trail,&glow);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"meteorTrail"),&trail,SHADER_UNIFORM_VEC4);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"meteorGlow"),&glow,SHADER_UNIFORM_FLOAT);
        Texture2D space=generated_space_texture();
        rlDisableDepthMask();BeginShaderMode(nightSkyShader);
        SetShaderValueTexture(nightSkyShader,GetShaderLocation(nightSkyShader,"spaceTexture"),space);
        DrawRectangle(0,0,w,h,WHITE);EndShaderMode();
        rlDrawRenderBatchActive();rlEnableDepthMask();
    }
    if(!worldPeriodic)return;
    // Sky geometry follows translation but retains world orientation, like a skybox.
    BeginMode3D(camera);rlDisableDepthMask();
    for(int n=0;n<12;n++) {
        float angle=n*(2*PI/12)+.13f,r=180,height=20+(n*17%43);
        Vector3 p=Vector3Add(camera.position,(Vector3){sinf(angle)*r,height*.5f-8,cosf(angle)*r});
        Color c={137+(n%3)*3,155+(n%3)*3,166+(n%3)*3,255};
        DrawCube(p,10+(n%3)*5,height,12,c);
        if(n%3==0) {
            DrawCube(Vector3Add(p,(Vector3){9,height*.25f,0}),30,3,5,c);
            DrawCube(Vector3Add(p,(Vector3){-3,height*.5f+5,0}),3,10,3,c);
        }
    }
    rlDrawRenderBatchActive();rlEnableDepthMask();EndMode3D();
}
