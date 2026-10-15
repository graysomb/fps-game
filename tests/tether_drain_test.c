#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void setup(void){
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;
    ResetGame();enemy_cleanup(true);clear_world_voxels();init_static_hash();worldPeriodic=false;
    gameState=GAME_STATE_PLAYING;activePlayers=2;enemyLab=false;
    for(int i=0;i<2;i++){players[i].pos=(Vector3){0,2,-i*5};players[i].yaw=players[i].pitch=0;
        players[i].respawn_timer=players[i].invuln_timer=0;players[i].matter=100;players[i].kd_ratio=1;players[i].isExposed=false;
        playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;}
    players[0].matter=0;players[0].isExposed=true;
}
static unsigned capture(float time,const char *path,bool bloom,bool covered){
    RenderTexture2D image=LoadRenderTexture(800,600);WeaponBloom effect={0};
    Camera3D camera={.position={3,3,4},.target={0,2,-3},.up={0,1,0},.fovy=55};
    Vector3 start={0,2,0},end={0,2,-5};
    BeginTextureMode(image);ClearBackground((Color){16,22,30,255});BeginMode3D(camera);
    DrawPlane((Vector3){0,0,0},(Vector2){20,20},DARKGRAY);DrawCube(end,.7f,1.2f,.7f,GRAY);
    if(covered)DrawCube((Vector3){0,2,2},20,10,.5f,GRAY);
    weapon_draw_drain(start,end,time);EndMode3D();EndTextureMode();
    if(bloom){assert(weapon_bloom_resize(&effect,800,600));weapon_bloom_begin(&effect,image);BeginMode3D(camera);weapon_draw_drain(start,end,time);EndMode3D();weapon_bloom_composite(&effect,image);EndTextureMode();}
    Image shot=LoadImageFromTexture(image.texture);ImageFlipVertical(&shot);assert(ExportImage(shot,path));
    Color *pixels=LoadImageColors(shot);unsigned hash=2166136261u;int colored=0;
    for(int n=0;n<shot.width*shot.height;n++){hash=(hash^pixels[n].g)*16777619u;colored+=pixels[n].g>pixels[n].r+40;}
    assert(covered?colored==0:colored>100);UnloadImageColors(pixels);UnloadImage(shot);
    weapon_bloom_unload(&effect);UnloadRenderTexture(image);
    return hash;
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,600,"Tether drain checks");init_pbd_thread_pool();assert(initialize_physics_backend());
    for(int hz=30;hz<=120;hz*=2){setup();start_tether(0);assert(players[0].drainTarget==2 && !players[0].tetherHolding);
        Vector3 pos=players[1].pos,velocity=players[1].vel;
        for(int frame=0;frame<hz;frame++)update_tether_drains(1.f/hz);
        assert(fabsf(players[1].matter-(100-2*player_bullet_damage(0)))<.01f);
        assert(players[0].matter==0 && Vector3Equals(pos,players[1].pos) && Vector3Equals(velocity,players[1].vel));
        update_tether_drains(100);assert(players[1].matter==0 && players[1].isExposed && players[1].respawn_timer==0 && !players[0].drainTarget);
    }
    setup();start_tether(0);players[1].pos.x=2;update_tether_drains(.1f);assert(players[0].drainTarget==2);
    Vector3 endpoint;assert(player_tether_visual_target(0,&endpoint) && Vector3Equals(endpoint,players[1].pos));
    release_tether(0);assert(!players[0].drainTarget);
    setup();start_tether(0);
    uint8_t data[1024];NetWriter writer;NetReader reader;NetPlayerWireState state;int slot;
    net_writer_init(&writer,data,sizeof(data));net_write_player_state(&writer,0);
    net_reader_init(&reader,data,writer.length);assert(net_read_player_state(&reader,&slot,&state));
    assert(slot==0 && (state.visual.flags&NET_PLAYER_VISUAL_DRAIN) && (state.visual.flags&NET_PLAYER_VISUAL_TETHER));
    assert(state.visual.tether_z==-5);
    setup();start_tether(0);players[1].pos.z=-9;update_tether_drains(.1f);assert(!players[0].drainTarget);
    setup();start_tether(0);players[0].yaw=180;update_tether_drains(.1f);assert(!players[0].drainTarget);
    setup();start_tether(0);add_static_box_at_grid(-3,3,0,8,-5,-5,GRAY);init_static_hash();update_tether_drains(.1f);assert(!players[0].drainTarget && players[1].matter==100);
    assert(!start_tether_drain(0));
    setup();gameMode=GAME_MODE_FIREFIGHT;combatantTeam[0]=combatantTeam[1]=0;assert(!start_tether_drain(0));
    combatantTeam[1]=1;assert(start_tether_drain(0));players[1].invuln_timer=1;update_tether_drains(.5f);assert(players[1].matter==100);
    players[1].invuln_timer=0;netTransport.role=NET_ROLE_CLIENT;update_tether_drains(.5f);assert(players[1].matter==100);netTransport.role=NET_ROLE_OFFLINE;
    ai_reset(1);assert(!players[0].drainTarget);
    setup();players[1].pos.z=-30;int v=addVoxel(0,2,-3,false,true,GRAY,0);rebuild_voxel_hash();start_tether(0);assert(players[0].tetherHolding && !players[0].drainTarget);release_tether(0);
    setup();players[0].matter=100;players[0].isExposed=false;players[1].yaw=180;
    start_tether(0);start_tether(1);update_tether_drains(.1f);
    assert(players[0].drainTarget==2 && players[1].drainTarget==1 && players[0].matter<100 && players[1].matter<100); // No parry.
    enemy_lab_start(ENEMY_TYPE_FOUNDRY);EnemySpecialState *boss=&enemySpecial[4];
    players[0].pos=v_add(players[4].pos,(Vector3){0,0,-4});players[0].yaw=180;
    assert(!start_tether_drain(0));
    for(int n=boss->equipmentCount-1;n>=0;n--)if(boss->equipment[n].group<2)enemy_detach(enemy_voxel(boss->equipment[n].id));
    enemy_update_equipment(4);rebuild_voxel_hash();assert(start_tether_drain(0));
    update_tether_drains(100);assert(boss->coreHealth==1 && players[4].respawn_timer==0 && !players[0].drainTarget);
    assert(weapon_renderer_init());assert(weapon_drain_shader_ready());
    MakeDirectory("artifacts");unsigned first=capture(0,"artifacts/tether-drain.png",false,false);
    assert(first!=capture(.15f,"artifacts/tether-drain-animated.png",false,false));capture(.15f,"artifacts/tether-drain-bloom.png",true,false);
    capture(0,"artifacts/tether-drain-occluded.png",true,true);
    weapon_renderer_shutdown();assert(!weapon_drain_shader_ready());
    capture(0,"artifacts/tether-drain-fallback.png",false,false);
    shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Tether drain checks passed");return 0;
}
