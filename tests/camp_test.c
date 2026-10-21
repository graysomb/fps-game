#define FIREFIGHT_TESTING
#include "raylib.h"
static bool testKey(int key){(void)key;return false;}
#define IsKeyPressed testKey
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(int humans){
    netTransport.role=NET_ROLE_OFFLINE;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;enemyLab=false;
    multiplayerActivePlayers=humans;ResetGame();gameState=GAME_STATE_PLAYING;firefightEndless=true;
    for(int i=0;i<MAX_PLAYERS;i++)playerInput[i]=INPUT_TYPE_KEYBOARD;
}
static void timing(void){
    fixture(2);ff_wave_complete();assert(camp.active && firefightIntermissionTimer==30);
    for(int i=0;i<2;i++)prepReady[i]=true;
    update_firefight_logic(4.9f);assert(firefightWaveState==FIREFIGHT_STATE_INTERMISSION && camp.active);
    prepReady[1]=false;update_firefight_logic(.2f);assert(firefightWaveState==FIREFIGHT_STATE_INTERMISSION);
    prepReady[1]=true;update_firefight_logic(.01f);assert(firefightWave==2 && !camp.active);
    ff_wave_complete();update_firefight_logic(29.9f);assert(firefightWaveState==FIREFIGHT_STATE_INTERMISSION);
    update_firefight_logic(.2f);assert(firefightWave==3 && !camp.active);
    for(int n=0;n<3;n++)teamUpgradeLevels[n]=teamUpgradeCaps[n];
    ff_wave_complete();update_firefight_logic(1);assert(camp.active && firefightWaveState==FIREFIGHT_STATE_INTERMISSION);
    for(int i=0;i<2;i++)prepReady[i]=true;update_firefight_logic(4.1f);
    assert(firefightWave==4 && !camp.active && ff_upgrades_maxed());
}
static void placement(void){
    for(int humans=1;humans<=4;humans++){
        fixture(humans);int count=voxel_count;ff_wave_complete();assert(camp.active);
        for(int i=0;i<humans;i++){
            assert(v_isfinite(players[i].pos));assert(spawn_position_clear_ex(players[i].pos,spawn_clear_half(i),0,i));
            Vector3 d=world_delta(camp.pos,players[i].pos);d.y=0;assert(v_dot(v_norm(d),player_forward(&players[i]))>.999f);
        }
        Vector3 old=camp.pos;players[0].pos.x+=5;camp_update(.1f);assert(Vector3Distance(old,camp.pos)==0);
        (void)count;camp_end();
    }
    fixture(4);clear_world_voxels();clear_pickups();init_static_hash();rebuild_voxel_hash();
    // Obstruct the preferred camp, forcing a deterministic safe alternative.
    add_static_box_at_grid(-2,2,0,6,-2,2,GRAY);init_static_hash();rebuild_voxel_hash();camp_begin();
    assert(camp.active && (fabsf(camp.pos.x)>1 || fabsf(camp.pos.z)>1));
    // Periodic center must stay near players straddling the world seam.
    clear_world_voxels();init_static_hash();rebuild_voxel_hash();firefightFoundry=false;worldPeriodic=true;firefightHumans=2;
    players[0].pos=(Vector3){FLOOR_SIZE-1,1,0};players[1].pos=(Vector3){-FLOOR_SIZE+1,1,0};camp_begin();
    assert(camp.active && fabsf(camp.pos.x)>FLOOR_SIZE-3);
    Vector3 previous=players[0].pos;float size=FLOOR_SIZE;FLOOR_SIZE=1;worldPeriodic=false;camp_begin();
    assert(!camp.active && Vector3Distance(previous,players[0].pos)==0);FLOOR_SIZE=size;
    enemyLab=true;camp_begin();assert(!camp.active);enemyLab=false;firefightPractice=true;camp_begin();assert(!camp.active);
}
static void captures(void){
    fixture(4);ff_wave_complete();assert(camp.active);camp_init();assert(campShader.id && campShader.id!=rlGetShaderIdDefault());
    gameState=GAME_STATE_PAUSED;float age=camp.age;camp_update(2);assert(camp.active && camp.age==age);gameState=GAME_STATE_PLAYING;
    Model models[3]={campSolid,campFlames,campSmoke};for(int m=0;m<3;m++)for(int n=0;n<models[m].meshes[0].vertexCount*3;n++)assert(isfinite(models[m].meshes[0].vertices[n]));
    unsigned int vao=campSolid.meshes[0].vaoId;int count=voxel_count;
    for(int n=0;n<40;n++){camp_begin();camp_update(.1f);camp_init();assert(campSolid.meshes[0].vaoId==vao && voxel_count==count);}
    Camera3D cam={v_add(camp.pos,(Vector3){2.4f,2.3f,3.4f}),v_add(camp.pos,(Vector3){0,.9f,0}),{0,1,0},42,CAMERA_PERSPECTIVE};
    RenderTexture2D rt=LoadRenderTexture(960,720);
    for(int n=0;n<3;n++){
        camp.age=n*1.3f;
        BeginTextureMode(rt);ClearBackground((Color){150,158,169,255});BeginMode3D(cam);
        DrawPlane(camp.pos,(Vector2){200,200},(Color){130,140,151,255});camp_draw(cam,false);EndMode3D();EndTextureMode();
        Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/camp-detail-%d.png",n));UnloadImage(im);
    }
    Camera3D blocked={v_add(camp.pos,(Vector3){0,1.3f,4}),v_add(camp.pos,(Vector3){0,1.3f,0}),{0,1,0},42,CAMERA_PERSPECTIVE};
    Image covered[2];for(int n=0;n<2;n++){
        BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(blocked);
        DrawCube(v_add(camp.pos,(Vector3){0,1.3f,2}),4,4,.2f,RED);if(n)camp_draw(blocked,false);
        EndMode3D();EndTextureMode();covered[n]=LoadImageFromTexture(rt.texture);
    }
    Color *a=LoadImageColors(covered[0]),*b=LoadImageColors(covered[1]);assert(!memcmp(a,b,960*720*sizeof(Color)));
    UnloadImageColors(a);UnloadImageColors(b);UnloadImage(covered[0]);UnloadImage(covered[1]);UnloadRenderTexture(rt);
    weaponRenderingReady=weapon_renderer_init();character_renderer_init();RenderTexture2D screens[4]={0};int views=0,w=0,h=0;
    for(int bloom=0;bloom<2;bloom++){
        weaponBloomEnabled=bloom!=0;
        for(int mode=0;mode<2;mode++){
            camp.active=mode!=0;
            for(int warm=0;warm<10;warm++)render_gameplay_view(screens,&views,&w,&h,false);
            double start=GetTime();for(int n=0;n<35;n++)render_gameplay_view(screens,&views,&w,&h,false);
            fprintf(stderr,"Camp four views bloom=%d camp=%d %.3f ms/frame\n",bloom,mode,(GetTime()-start)*1000/35);
            if(mode){Image im=LoadImageFromScreen();ExportImage(im,TextFormat("artifacts/camp-four-bloom-%d.png",bloom));UnloadImage(im);}
        }
    }
    for(int i=0;i<views;i++)UnloadRenderTexture(screens[i]);
    memset(screens,0,sizeof(screens));views=w=h=0;firefightHumans=multiplayerActivePlayers=1;
    for(int n=0;n<3;n++)render_gameplay_view(screens,&views,&w,&h,false);
    Image full=LoadImageFromScreen();ExportImage(full,"artifacts/camp-full.png");UnloadImage(full);UnloadRenderTexture(screens[0]);
    // Exercise the real fixed-function fallback without a live custom shader.
    Shader saved=campShader;Shader fallback={.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};
    campShader=(Shader){0};campSolid.materials[0].shader=campFlames.materials[0].shader=campSmoke.materials[0].shader=fallback;
    for(int n=0;n<3;n++){BeginDrawing();ClearBackground(BLACK);BeginMode3D(cam);camp_draw(cam,false);EndMode3D();EndDrawing();}
    Image fallbackShot=LoadImageFromScreen();ExportImage(fallbackShot,"artifacts/camp-fallback.png");UnloadImage(fallbackShot);campShader=saved;
    camp_shutdown();assert(!campResources && !camp.active);
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Camp validation");
    init_pbd_thread_pool();assert(initialize_physics_backend());timing();placement();captures();
    for(int i=0;i<4;i++)weapon_bloom_unload(&weaponBlooms[i]);shutdown_world_visuals();weapon_renderer_shutdown();character_renderer_shutdown();
    clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Camp checks passed");return 0;
}
