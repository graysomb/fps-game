#define TOWER_PREVIEW
#include "builder_towers_test.c"
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,500,"Builder tower preview");init_pbd_thread_pool();assert(initialize_physics_backend());
    weaponRenderingReady=weapon_renderer_init();character_renderer_init();voxel_ca_init();fixture();
    Camera3D camera={{17,12,28},{2,3,14},{0,1,0},50,CAMERA_PERSPECTIVE};RenderTexture2D target=LoadRenderTexture(800,500);
    EnemyStructure *t=NULL;int output=0;bool repair=false,jump=false;const char *label="BUILDER CONSTRUCTS A WATCHTOWER";
    for(int f=0;f<1200;f++){
        t=built();
        if(f==720 && t){for(int v=voxel_count-1;v>=0;v--)if(voxels[v].towerTag==t->tag && voxels[v].gy==4)remove_voxel_index(v);repair=true;label="NEARBY BUILDER REPAIRS DAMAGE";}
        if(f==900 && t){players[0].pos=tower_point(t,(Vector3){5,1,7.5f});players[0].matter=0;players[0].isExposed=true;ai_select(5,0,players[0].pos,AI_ALERT_SEEN);jump=true;label="EXPOSED PLAYER: RIFLEMAN JUMPS DOWN";}
        step(1.f/60);
        if(tower_stationed(5) && !repair)label="RIFLEMAN CLIMBS AND HOLDS THE TOWER";
        if(f%6)continue;
        weaponFrameSampled=true;weaponFrameTime=(float)GetTime();for(int i=0;i<activePlayers;i++)sample_character(i,weaponFrameTime);
        voxel_render_frame++;voxel_ca_step();drawViewPlayerIndex=-1;
        BeginTextureMode(target);ClearBackground((Color){30,37,46,255});BeginMode3D(camera);
        DrawPlane((Vector3){0,-.02f,10},(Vector2){80,80},(Color){55,61,68,255});DrawVoxels(camera);draw_players();EndMode3D();
        DrawText(label,18,18,20,SKYBLUE);EndTextureMode();
        Image im=LoadImageFromTexture(target.texture);ImageFlipVertical(&im);assert(ExportImage(im,TextFormat("artifacts/builder-tower/frame-%03d.png",output++)));UnloadImage(im);
    }
    assert(repair&&jump);UnloadRenderTexture(target);enemy_cleanup(true);voxel_ca_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
