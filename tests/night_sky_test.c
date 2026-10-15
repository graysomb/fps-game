#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static Image sky_image(RenderTexture2D target,Camera3D cam) {
    BeginTextureMode(target);ClearBackground(BLACK);periodic_skyline(cam);EndTextureMode();
    Image im=LoadImageFromTexture(target.texture);ImageFlipVertical(&im);return im;
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Night sky regression");
    worldPeriodic=false;
    Camera3D camera={{0,3,0},{0,23,-3},{0,1,0},65,CAMERA_PERSPECTIVE};
    for(int small=0;small<2;small++) {
        RenderTexture2D target=LoadRenderTexture(small?480:960,small?270:540);
        Image a=sky_image(target,camera);assert(nightSkyShader.id && nightSkyShader.id!=rlGetShaderIdDefault());assert(worldSpaceTexture.id && generated_space_texture().id==worldSpaceTexture.id);
        Camera3D shifted=camera;shifted.position.x+=80;shifted.target.x+=80;
        Image b=sky_image(target,shifted);Color *pa=LoadImageColors(a),*pb=LoadImageColors(b);int bright=0;
        for(int i=0;i<a.width*a.height;i++){assert(abs(pa[i].r-pb[i].r)<=1 && abs(pa[i].b-pb[i].b)<=1);bright+=pa[i].b>100 && i%a.width>a.width/4 && i%a.width<a.width*3/4 && i/a.width>a.height/4 && i/a.width<a.height*3/4;}
        assert(bright>25);fprintf(stderr,"Sky %dx%d: %d star pixels; translation invariant\n",a.width,a.height,bright);
        UnloadImageColors(pa);UnloadImageColors(pb);UnloadImage(a);UnloadImage(b);
        Camera3D horizon={{0,3,0},{0,3,-20},{0,1,0},65,CAMERA_PERSPECTIVE};
        Image fog=sky_image(target,horizon);Color *pixels=LoadImageColors(fog);
        for(int x=0;x<fog.width;x++) {
            Color c=pixels[(fog.height/2)*fog.width+x];
            assert(abs(c.r-periodicSky.r)<=1 && abs(c.g-periodicSky.g)<=1 && abs(c.b-periodicSky.b)<=1);
        }
        UnloadImageColors(pixels);UnloadImage(fog);UnloadRenderTexture(target);
    }
    voxel_ca_init();character_renderer_init();weaponRenderingReady=weapon_renderer_init();
    for(int map=0;map<2;map++) {
        gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;
        currentWorldType=WORLD_TYPE_CROSSCUT;multiplayerActivePlayers=firefightHumans=1;ResetGame();assert(worldPeriodic);
        RenderTexture2D target=world_render_target(960,540);periodicAspect=960.f/540;
        for(int view=0;view<2;view++) {
        Camera3D cam={{12,3,15},{0,3,-8},{0,1,0},65,CAMERA_PERSPECTIVE};
        if(view)cam.target=(Vector3){12,23,12};
        voxel_render_frame++;BeginTextureMode(target);ClearBackground(periodicSky);periodic_skyline(cam);
        for(int n=0;n<9;n++){Camera3D c;if(!periodic_camera(cam,n,&c))continue;BeginMode3D(c);draw_world_surfaces();DrawVoxels(c);draw_pickups(c);EndMode3D();}
        periodic_fog(target,cam,false,target.depth);EndTextureMode();
        assert(periodicFogShader.id!=rlGetShaderIdDefault());
        Image im=LoadImageFromTexture(target.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/night-sky-map%d-view%d.png",map,view));UnloadImage(im);
        }
        UnloadRenderTexture(target);
    }
    clear_all_chunks();clear_world_voxels();shutdown_world_visuals();assert(!nightSkyShader.id && !worldSpaceTexture.id);
    character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();puts("Night sky regression passed");
}
