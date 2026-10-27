#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,500,"Drop slab trigger preview");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    activePlayers=0;worldPeriodic=false;worldFloorHalfSize=40;
    gameMode=GAME_MODE_FIREFIGHT;firefightHumans=1;firefightPractice=false;firefightFoundry=true;
    clear_all_chunks();clear_world_voxels();
    for(unsigned seed=1;seed<10000;seed++){
        SetRandomSeed(seed);int b=GetRandomValue(0,3),p=GetRandomValue(0,3);
        if(b==3 && p==3){SetRandomSeed(seed);break;}
    }
    buildFoundryHoldout();voxel_ca_init();
    bool bay=getenv("FPS_BAY_PREVIEW")!=NULL;
    bool pillar=getenv("FPS_PILLAR_PREVIEW")!=NULL;
    FoundryTrap *trap=&foundryTraps[pillar?3:bay?2:1];
    Vector3 trigger={(trap->x1+.5f)*VOXEL_SIZE,(trap->ry+.5f)*VOXEL_SIZE,(trap->z1+.5f)*VOXEL_SIZE};
    Camera3D camera={{29,6.5f,12},{21.5f,1.8f,0},{0,1,0},55,CAMERA_PERSPECTIVE};
    if(bay){trigger=(Vector3){(trap->rx+.5f)*VOXEL_SIZE,(trap->ry+.5f)*VOXEL_SIZE,(trap->rz+.5f)*VOXEL_SIZE};camera=(Camera3D){{35,10,8},{21.5f,2.5f,21},{0,1,0},55,CAMERA_PERSPECTIVE};}
    if(pillar){trigger=(Vector3){(trap->rx+.5f)*VOXEL_SIZE,(trap->ry+.5f)*VOXEL_SIZE,(trap->rz+.5f)*VOXEL_SIZE};camera=(Camera3D){{46,13,21},{35,4,35},{0,1,0},55,CAMERA_PERSPECTIVE};}
    float approach=(bay||pillar)?-3:3;
    RenderTexture2D target=LoadRenderTexture(800,500);
    for(int frame=0;frame<100;frame++) {
        if(frame==24){
            Vector3 origin=v_add(trigger,(Vector3){0,0,approach});
            int shot=addVoxel(origin.x,origin.y,origin.z,false,true,WHITE,0);
            assert(shot>=0);voxels[shot].isBullet=true;voxels[shot].owner=-1;voxels[shot].vel=(Vector3){0,0,(bay||pillar)?30:-30};
        }
        if(frame>=24)for(int step=0;step<3;step++){
            update_projectiles(1.f/60);foundry_traps_prepare();simulate_voxel_pbd_steps(1.f/60,1);
        }
        voxel_ca_step();voxel_render_frame++;
        BeginTextureMode(target);ClearBackground((Color){30,37,46,255});BeginMode3D(camera);
        DrawPlane((Vector3){0,-.02f,0},(Vector2){80,80},(Color){55,61,68,255});DrawVoxels(camera);foundry_traps_draw();
        if(frame>=24 && frame<28)DrawLine3D(v_add(trigger,(Vector3){0,0,approach}),trigger,YELLOW);
        EndMode3D();DrawText(trap->active?(pillar?"CORNER PILLAR COLLAPSE":bay?"FACTORY BAY COLLAPSE":"ALL FOUR SUPPORTS RELEASED"):"SHOOT ANY ORANGE TRIGGER",24,22,24,trap->active?ORANGE:RAYWHITE);EndTextureMode();
        Image im=LoadImageFromTexture(target.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/%s/frame-%03d.png",pillar?"pillar-collapse":bay?"bay-collapse":"drop-slab",frame));UnloadImage(im);
    }
    assert(trap->active);UnloadRenderTexture(target);voxel_ca_shutdown();clear_all_chunks();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    return 0;
}
