#define FIREFIGHT_TESTING
#include "raylib.h"
/* Exercise the actual couch shell with deterministic virtual devices. */
static bool testKeys[512],testPads[4],testButtons[4][32];
static bool test_key_pressed(int key) {return key>=0&&key<512&&testKeys[key];}
static bool test_pad_available(int pad) {return pad>=0&&pad<4&&testPads[pad];}
static bool test_button_pressed(int pad,int button) {return test_pad_available(pad)&&button>=0&&button<32&&testButtons[pad][button];}
static float test_pad_axis(int pad,int axis) {(void)pad;(void)axis;return 0;}
#define IsKeyPressed test_key_pressed
#define IsGamepadAvailable test_pad_available
#define IsGamepadButtonPressed test_button_pressed
#define GetGamepadAxisMovement test_pad_axis
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static int compare_ms(const void *a,const void *b) {
    double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);
}

static void test_down(int i) {
    players[i].invuln_timer=0;
    players[i].isExposed=true;
    kill_player(i,4,false,false);
    assert(combatantLife[i]==LIFE_DOWNED);
}

int main(int argc,char **argv) {
    if(!parse_physics_arguments(argc,argv))return 2;
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1920,1080,"Firefight integration");
    assert(IsWindowReady());
    init_pbd_thread_pool();assert(initialize_physics_backend());
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;
    multiplayerActivePlayers=4;firefightPractice=false;
    ResetGame();gameState=GAME_STATE_PLAYING;
    assert(activePlayers==16 && firefightHumans==4);
    for(int i=0;i<4;i++) assert(combatantLife[i]==LIFE_ALIVE);
    for(int i=4;i<16;i++) assert(players[i].respawn_timer>0);
    start_firefight_wave(5);
    for(int n=0;n<30;n++)update_firefight_logic(2.1f);
    int alive=0,giants=0;
    for(int i=4;i<16;i++) if(players[i].respawn_timer<=0) {alive++;giants+=players[i].enemyType==ENEMY_TYPE_GOLIATH;}
    assert(alive==12 && giants<=2);
    /* A real 16-combatant solver dispatch catches GPU uniform/tether bounds. */
    simulate_voxel_pbd_steps(1.0f/120,2);
    players[0].pos=(Vector3){0,4,0};players[1].pos=(Vector3){1,4,0};
    test_down(1);assert(firefightTeamLives==3);
    ffTestReviveHeld[0]=true;players[0].last_damage_time=-1000;
    ff_update_coop(1);assert(reviveProgress[1]>0);
    ffTestReviveHeld[0]=false;ff_update_coop(.1f);assert(reviveProgress[1]==0);
    players[0].matter=25;ffTestReviveHeld[0]=true;ff_update_coop(.1f);assert(reviveProgress[1]==0);
    players[0].matter=100;players[0].pos.x=10;ff_update_coop(.1f);assert(reviveProgress[1]==0);
    players[0].pos.x=0;players[0].last_damage_time=(float)GetTime();ff_update_coop(.1f);assert(reviveProgress[1]==0);
    players[0].last_damage_time=-1000;
    ffTestReviveHeld[0]=true;
    ff_update_coop(3.1f);
    assert(combatantLife[1]==LIFE_ALIVE && players[1].matter==35);
    assert(players[0].matter==75 && reviveCount[0]==1);
    ffTestReviveHeld[0]=false;
    test_down(1);ff_update_coop(16);assert(combatantLife[1]==LIFE_WAITING && firefightTeamLives==2);
    ff_update_coop(5.1f);assert(combatantLife[1]==LIFE_ALIVE);
    firefightTeamLives=0;test_down(1);ff_update_coop(16);assert(combatantLife[1]==LIFE_SPECTATING);
    firefightWave=2;teamUpgrade=UPGRADE_THROW;ff_wave_complete();
    assert(combatantLife[1]==LIFE_ALIVE && players[1].matter==MATTER_MAX_DEFAULT);
    assert(firefightTeamLives==1 && teamUpgrade==UPGRADE_NONE);
    for(int i=0;i<4;i++){prepVotes[i]=UPGRADE_HARVEST;prepReady[i]=true;}
    update_firefight_logic(5.1f);assert(firefightWave==3 && teamUpgrade==UPGRADE_HARVEST);
    assert(ff_harvest_amount(0)==15 && ff_harvest_amount(4)==MATTER_MELEE_HARVEST);
    firefightWave=5;firefightEndless=false;ff_wave_complete();assert(firefightWaveState==FIREFIGHT_STATE_VICTORY);
    firefightEndless=true;ff_wave_complete();assert(firefightWaveState==FIREFIGHT_STATE_INTERMISSION);
    /* Four views render without indexing enemy slots as cameras. */
    RenderTexture2D screens[4]={0};int count=0,w=0,h=0;
    gameState=GAME_STATE_PLAYING;
    render_gameplay_view(screens,&count,&w,&h,false);
    assert(count==4);
    Image view=LoadImageFromTexture(screens[0].texture);
    ImageFlipVertical(&view);
    ExportImage(view,"artifacts/firefight-quarter-view.png");UnloadImage(view);
    for(int humans=1;humans<=4;humans++) {
        firefightHumans=humans;
        render_gameplay_view(screens,&count,&w,&h,false);assert(count==humans);
    }
    firefightHumans=4;
    /* Refill the arena and measure actual simulation plus four-view rendering. */
    start_firefight_wave(5);
    for(int n=0;n<30;n++)update_firefight_logic(2.1f);
    double samples[30];
    for(int n=0;n<35;n++) {
        double before=GetTime();
        simulate_voxel_pbd_steps(1.0f/120,2);
        render_gameplay_view(screens,&count,&w,&h,false);
        if(n>=5)samples[n-5]=(GetTime()-before)*1000;
    }
    qsort(samples,30,sizeof(double),compare_ms);
    printf("Four-view firefight (%s): p50 %.2f ms, p95 %.2f ms\n",physics_backend_name(physicsBackend.active),samples[15],samples[28]);
    for(int i=0;i<count;i++)UnloadRenderTexture(screens[i]);
    /* Team wipe with no reserves must terminate, not leave downed bodies forever. */
    firefightWaveState=FIREFIGHT_STATE_ACTIVE;firefightTeamLives=0;
    for(int i=0;i<4;i++)test_down(i);
    ff_update_coop(.1f);assert(firefightWaveState==FIREFIGHT_STATE_DEFEAT);
    /* Solo uses reserves directly and cannot wait for an impossible revive. */
    multiplayerActivePlayers=1;firefightEndless=false;ResetGame();gameState=GAME_STATE_PLAYING;
    players[0].invuln_timer=0;kill_player(0,4,false,false);
    assert(combatantLife[0]==LIFE_WAITING && firefightTeamLives==2);
    ff_update_coop(5.1f);assert(combatantLife[0]==LIFE_ALIVE);
    firefightPractice=true;ResetGame();gameState=GAME_STATE_PLAYING;
    practiceStep=1;update_firefight_logic(2.1f);assert(firefightEnemiesSpawned==1);
    players[4].invuln_timer=0;apply_matter_damage(4,0,100);assert(practiceStep==2);
    kill_player(4,0,true,false);assert(practiceStep==3);
    update_firefight_logic(2.1f);players[4].invuln_timer=0;
    apply_matter_damage(4,0,100);kill_player(4,0,true,true);
    assert(practiceStep==4 && firefightWaveState==FIREFIGHT_STATE_VICTORY);
    firefightPractice=false;
    /* Controller-only menu -> four-seat lobby -> game -> disconnect pause -> rematch. */
    gameState=GAME_STATE_MENU;couchSelection=0;couchLobby=false;
    for(int i=0;i<4;i++)testPads[i]=true;
    testButtons[0][GAMEPAD_BUTTON_MIDDLE_RIGHT]=true;couch_frame();memset(testButtons,0,sizeof(testButtons));
    testButtons[0][GAMEPAD_BUTTON_RIGHT_FACE_DOWN]=true;couch_frame();memset(testButtons,0,sizeof(testButtons));
    assert(couchLobby && multiplayerActivePlayers==1);
    for(int i=1;i<4;i++)testButtons[i][GAMEPAD_BUTTON_MIDDLE_RIGHT]=true;
    couch_frame();memset(testButtons,0,sizeof(testButtons));assert(multiplayerActivePlayers==4);
    for(int i=0;i<4;i++)testButtons[i][GAMEPAD_BUTTON_RIGHT_FACE_DOWN]=true;
    couch_frame();memset(testButtons,0,sizeof(testButtons));assert(gameState==GAME_STATE_PLAYING);
    testPads[2]=false;couch_frame();assert(gameState==GAME_STATE_PAUSED);
    testPads[2]=true;testButtons[0][GAMEPAD_BUTTON_RIGHT_FACE_DOWN]=true;
    couch_frame();memset(testButtons,0,sizeof(testButtons));assert(gameState==GAME_STATE_PLAYING);
    gameState=GAME_STATE_GAMEOVER;couchSelection=0;
    testButtons[0][GAMEPAD_BUTTON_RIGHT_FACE_DOWN]=true;couch_frame();memset(testButtons,0,sizeof(testButtons));
    assert(gameState==GAME_STATE_PLAYING && multiplayerActivePlayers==4 && firefightWave==1);
    /* Consume all authored queues through the real death handler. */
    multiplayerActivePlayers=1;ResetGame();gameState=GAME_STATE_PLAYING;
    for(int wave=1;wave<=5;wave++) {
        int watchdog=200;
        while(firefightWaveState==FIREFIGHT_STATE_ACTIVE && watchdog-->0) {
            update_firefight_logic(2.1f);
            for(int i=4;i<16;i++) if(players[i].respawn_timer<=0) kill_player(i,0,true,true);
        }
        assert(watchdog>0 && firefightWaveEnemiesRemaining==0);
        if(wave<5) {
            assert(firefightWaveState==FIREFIGHT_STATE_INTERMISSION);
            prepReady[0]=true;update_firefight_logic(5.1f);
            assert(firefightWave==wave+1);
        }
    }
    assert(firefightWaveState==FIREFIGHT_STATE_VICTORY);
    assert(physicsBackend.requested==PHYSICS_BACKEND_CPU_MT || physics_backend_is_gpu(physicsBackend.active));
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;ResetGame();assert(activePlayers==2);
    currentWorldType=WORLD_TYPE_CROSSCUT;useCustomMap=false;
    for(int seats=2;seats<=4;seats++) {
        multiplayerActivePlayers=seats;ResetGame();
        assert(activePlayers==seats && pickups[0].type==PICKUP_GOLD_TETHER);
        for(int i=0;i<seats;i++) {
            assert(spawn_position_clear_ex(players[i].pos,spawn_clear_half(i),2,i));
            for(int j=0;j<i;j++) assert(is_view_occluded_by_voxels(players[i].pos,players[j].pos));
        }
    }
    int originalCount=voxel_count;
    /* Cutting one isolated slab's feet must activate its debris locally. */
    remove_static_voxels_in_region_recycle(-3,2,0,3,5,5);
    assert(activate_static_voxels_near_region(-3,2,0,3,5,5,0));
    bool slabDynamic=false;
    for(int v=0;v<voxel_count;v++) if(voxels[v].simulate && fabsf(voxels[v].pos.x)<1.6f && voxels[v].pos.z>=2 && voxels[v].pos.z<4) slabDynamic=true;
    assert(slabDynamic);ResetGame();
    for(int phase=0;phase<2;phase++) {
        if(phase) {
            for(int v=voxel_count-1;v>=0;v--) if(fabsf(voxels[v].pos.x)<9 && fabsf(voxels[v].pos.z)<9) remove_voxel_index(v);
            init_static_hash();rebuild_all_voxel_surfaces();meshDirty=true;
        }
        RenderTexture2D views[MAX_PLAYERS]={0};int viewCount=0,vw=0,vh=0;
        double timings[30];
        for(int n=0;n<35;n++) {
            double start=GetTime();simulate_voxel_pbd_steps(1.0f/120,2);
            render_gameplay_view(views,&viewCount,&vw,&vh,false);
            if(n>=5)timings[n-5]=(GetTime()-start)*1000;
        }
        qsort(timings,30,sizeof(double),compare_ms);
        printf("Crosscut four-view %s (%s): p50 %.2f ms, p95 %.2f ms\n",phase?"crossing removed":"intact",physics_backend_name(physicsBackend.active),timings[15],timings[28]);
        if(!phase) {
            Image eye=LoadImageFromTexture(views[0].texture);ImageFlipVertical(&eye);
            ExportImage(eye,"artifacts/crosscut-player-view.png");UnloadImage(eye);
        }
        for(int i=0;i<viewCount;i++)UnloadRenderTexture(views[i]);
    }
    ResetGame();
    RenderTexture2D mapShot=LoadRenderTexture(1000,1000);
    Camera3D mapCamera={.position={29,38,29},.target={0,0,0},.up={0,1,0},.fovy=48,.projection=CAMERA_PERSPECTIVE};
    BeginTextureMode(mapShot);ClearBackground((Color){30,35,43,255});BeginMode3D(mapCamera);
    DrawPlane((Vector3){0,0,0},(Vector2){40,40},(Color){74,77,72,255});
    for(int v=0;v<voxel_count;v++) DrawCube(voxels[v].pos,VOXEL_SIZE,VOXEL_SIZE,VOXEL_SIZE,voxels[v].color);
    EndMode3D();EndTextureMode();
    Image mapImage=LoadImageFromTexture(mapShot.texture);ImageFlipVertical(&mapImage);
    ExportImage(mapImage,"artifacts/crosscut-overview.png");UnloadImage(mapImage);UnloadRenderTexture(mapShot);
    /* Remove the crossing and verify all authored ground supplies/spawns survive. */
    for(int v=voxel_count-1;v>=0;v--) if(fabsf(voxels[v].pos.x)<9 && fabsf(voxels[v].pos.z)<9) remove_voxel_index(v);
    init_static_hash();
    for(int i=0;i<4;i++) {
        Vector3 p=pick_player_spawn(i);
        assert(spawn_position_clear_ex(p,spawn_clear_half(i),3,i));
    }
    ResetGame();assert(voxel_count==originalCount && currentWorldType==WORLD_TYPE_CROSSCUT);
    gameState=GAME_STATE_MENU;couchLobby=true;
    testButtons[0][GAMEPAD_BUTTON_LEFT_FACE_RIGHT]=true;couch_frame();memset(testButtons,0,sizeof(testButtons));
    assert(currentWorldType==WORLD_TYPE_GREEK_TEMPLE && !couchReady[0]);
    testButtons[0][GAMEPAD_BUTTON_LEFT_FACE_LEFT]=true;couch_frame();memset(testButtons,0,sizeof(testButtons));
    assert(currentWorldType==WORLD_TYPE_CROSSCUT);
    puts("Crosscut: 2/3/4 seats, screened spawns, destruction fallback, rematch passed");
    shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Firefight integration passed");return 0;
}
