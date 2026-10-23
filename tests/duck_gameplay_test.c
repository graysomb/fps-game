#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(bool firefight) {
    netTransport.role=NET_ROLE_OFFLINE;gameMode=firefight?GAME_MODE_FIREFIGHT:GAME_MODE_DEATHMATCH;
    currentWorldType=WORLD_TYPE_CROSSCUT;firefightFoundry=true;firefightPractice=false;enemyLab=false;creativeModeActive=false;useCustomMap=false;
    multiplayerActivePlayers=2;ResetGame();gameState=GAME_STATE_PLAYING;
    for(int i=0;i<activePlayers;i++)players[i].pos=(Vector3){0,8,0};
    ducks_tick(.01f);assert(ducksInitialized);
}
static void placements(void) {
    for(int mode=0;mode<2;mode++) {
        fixture(mode);
        for(int n=0;n<4;n++) {
            if(!mapDucks[n].alive){fprintf(stderr,"Blocked duck map=%d index=%d at %.2f %.2f\n",mode,n,mapDucks[n].pos.x,mapDucks[n].pos.z);
                BoundingBox b=duck_bounds(mapDucks[n].pos,.06f);for(int z=floorf(b.min.z/VOXEL_SIZE);z<=floorf(b.max.z/VOXEL_SIZE);z++)for(int x=floorf(b.min.x/VOXEL_SIZE);x<=floorf(b.max.x/VOXEL_SIZE);x++)for(int y=0;y<2;y++){int v=table_get(x,y,z);if(v>=0)fprintf(stderr," occupied cell %d %d %d voxel center %.2f %.2f %.2f\n",x,y,z,voxels[v].pos.x,voxels[v].pos.y,voxels[v].pos.z);}}
            assert(mapDucks[n].alive);
            bool crossed=false;
            Vector3 previous=duck_route(&mapDucks[n],0);
            for(int step=1;step<=128;step++) {
                Vector3 at=duck_route(&mapDucks[n],step*2*PI/128);
                assert(duck_clear(at,true));
                assert(fmaxf(fabsf(at.x),fabsf(at.z))>worldFloorHalfSize-5);
                assert(v_length(world_delta(at,previous))<.16f);
                crossed|=v_length(v_sub(at,previous))>worldFloorHalfSize;
                previous=at;
            }
            assert(!worldPeriodic || crossed);
        }
        int count=voxel_count;for(int frame=0;frame<120;frame++)ducks_tick(1.f/60);
        assert(voxel_count==count);for(int n=0;n<4;n++)assert(mapDucks[n].alive);
    }
    puts("Duck placement: both maps, four clear perimeter routes, continuous seam crossings and no physics allocations passed");
}
static void empty_fixture(void) {
    fixture(false);clear_world_voxels();init_static_hash();rebuild_voxel_hash();
    worldPeriodic=false;activePlayers=1;players[0].pos=(Vector3){0,1,4};players[0].yaw=0;players[0].pitch=-10;
    playerInput[0]=INPUT_TYPE_KEYBOARD;players[0].matter=100;players[0].matterMax=100;players[0].isExposed=false;combatantLife[0]=LIFE_ALIVE;
    for(int n=0;n<4;n++)mapDucks[n]=(MapDuck){.home={20,0,20},.pos={20,0,20},.alive=false,.respawn=90};
    mapDucks[0]=(MapDuck){.home={-DUCK_ROUTE_RADIUS,0,0},.pos={0,0,0},.alive=true};
}
static Ray shot(void){return (Ray){{0,.3f,3},{0,0,-1}};}
static void revive_duck(void){mapDucks[0].alive=true;mapDucks[0].respawn=0;}
static void hits(void) {
    empty_fixture();int count=voxel_count,score=firefightScore;
    assert(ducks_ray(shot(),5,0,-1,false));assert(!mapDucks[0].alive && mapDucks[0].respawn==90 && mapDucks[0].bursts==1);
    assert(!ducks_ray(shot(),5,0,-1,false));assert(mapDucks[0].bursts==1);
    ducks_tick(89.9f);assert(!mapDucks[0].alive);ducks_tick(.11f);assert(mapDucks[0].alive);
    assert(voxel_count==count && firefightScore==score);
    gameState=GAME_STATE_PAUSED;duck_pop(0);ducks_tick(100);assert(mapDucks[0].respawn==90);gameState=GAME_STATE_PLAYING;
    // Cover blocks all ordinary contacts and splash; rail beams deliberately penetrate.
    revive_duck();add_static_box_at_grid(-2,2,0,3,2,2,GRAY);init_static_hash();
    assert(!ducks_ray(shot(),5,0,-1,false));ducks_blast((Vector3){0,.3f,2},4);assert(mapDucks[0].alive);
    assert(ducks_ray(shot(),5,0,-1,true));
    clear_world_voxels();init_static_hash();revive_duck();ducks_blast((Vector3){0,.3f,2},4);assert(!mapDucks[0].alive);
    revive_duck();ducks_segment((Vector3){0,.3f,3},(Vector3){0,.3f,-3},.3f,-1);assert(!mapDucks[0].alive);
    revive_duck();players[0].pos=(Vector3){0,1,0};ducks_tick(.01f);assert(!mapDucks[0].alive);
    ducks_tick(91);assert(!mapDucks[0].alive);players[0].pos=(Vector3){0,1,4};ducks_tick(1.1f);assert(mapDucks[0].alive);
    // Wrapped traces and blasts see exactly one copy.
    worldPeriodic=true;worldFloorHalfSize=32;mapDucks[0].pos=(Vector3){-31,0,0};revive_duck();
    assert(ducks_ray((Ray){{31,.3f,0},{1,0,0}},4,0,-1,false));
    puts("Duck hits: rays, swept debris, blast occlusion, contact, pause, safe 90-second respawn, periodic seam passed");
}
static void weapons(void) {
    empty_fixture();players[0].pos=(Vector3){0,.3f,3};players[0].pitch=0;
    const int guns[]={EQ_GLOCK,EQ_DMR,EQ_SHOTGUN,EQ_CUTTER,EQ_RAIL,EQ_ROCKET,EQ_GLOB};
    for(int k=0;k<7;k++) {
        revive_duck();eq_reset_all();players[0].matter=100;players[0].isExposed=false;
        equipment[0].gun[0]=guns[k];equipment[0].selected=0;
        if(guns[k]==EQ_CUTTER)eq_command(0,&(NetInputCommand){.held=NET_INPUT_FIRE},.1f);
        else if(guns[k]==EQ_RAIL){eq_command(0,&(NetInputCommand){.held=NET_INPUT_FIRE},.65f);eq_command(0,&(NetInputCommand){0},.01f);}
        else {assert(eq_fire(0));for(int frame=0;frame<20;frame++)eq_tick(.01f);}
        assert(!mapDucks[0].alive);
    }
    revive_duck();clarity_explode((Vector3){0,.3f,1},0);assert(!mapDucks[0].alive);
    revive_duck();melee_hit_players(0,(Vector3){0,.3f,1},(Vector3){0,.3f,0},(Vector3){0,0,-1});assert(!mapDucks[0].alive);
    revive_duck();players[0].tetherHolding=false;players[0].drainTarget=0;start_tether(0);assert(!mapDucks[0].alive);
    revive_duck();eq_reset_all();players[0].matter=100;players[0].isExposed=false;players[0].last_shot_time=-100;
    FireVoxel(0);update_projectiles(.1f);assert(!mapDucks[0].alive);
    revive_duck();equipment[0].gun[0]=EQ_BLADE;players[0].last_melee_time=-100;players[0].meleeSwingActive=false;
    assert(perform_sword(0));for(int frame=0;frame<60;frame++)eq_sword_tick(1.f/60);assert(!mapDucks[0].alive);
    puts("Duck actual weapon paths: default, Glock, DMR, shotgun, cutter, railgun, rocket, Glob, grenade, melee, sword and tether passed");
}
static void network(void) {
    empty_fixture();duck_pop(0);duckClock=12.5f;
    uint8_t data[256];NetWriter w;net_writer_init(&w,data,sizeof(data));ducks_write(&w);assert(!w.failed);
    MapDuck expected[4];memcpy(expected,mapDucks,sizeof(expected));ducks_reset();netTransport.role=NET_ROLE_CLIENT;
    NetReader r;net_reader_init(&r,data,w.length);assert(ducks_receive(NET_MSG_DUCKS,&r,100));
    assert(ducksInitialized && !mapDucks[0].alive && mapDucks[0].respawn==90 && mapDucks[0].bursts==1 && duckClock==12.5f);
    assert(!ducks_ray(shot(),5,0,-1,false));ducks_tick(91);assert(!mapDucks[0].alive);
    net_reader_init(&r,data,w.length);ducks_receive(NET_MSG_DUCKS,&r,99);assert(mapDucks[0].respawn==0);
    net_reader_init(&r,data,w.length-1);ducks_receive(NET_MSG_DUCKS,&r,101);assert(duckReceiveTick==100);
    netTransport.role=NET_ROLE_OFFLINE;memcpy(mapDucks,expected,sizeof(expected));
    puts("Duck network: late state, host authority, old and truncated snapshot rejection passed");
}
static void lan(bool host) {
    fixture(false);assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27946));net_set_host_local_players(1);netLobbyStarted=true;duck_pop(0);}
    else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27946));}
    double start=GetTime(),joined=-1;bool sawInitial=false,sawHit=false,sawRespawn=false,popped=false,restored=false;
    while(GetTime()-start<9) {
        net_transport_pump(&netTransport,net_on_receive,net_on_connect,NULL);
        bool connected=host?netTransport.peers[1].welcomed:netWorldReady;
        if(connected){if(joined<0)joined=GetTime();float age=GetTime()-joined;
            if(host){
                if(age>1 && !popped){duck_pop(1);popped=true;}
                if(age>3 && !restored){mapDucks[0].respawn=mapDucks[1].respawn=.01f;restored=true;}
                ducks_tick(1.f/60);net_host_tick();
            }else {
                sawInitial|=ducksInitialized && !mapDucks[0].alive && mapDucks[0].bursts==1;
                sawHit|=ducksInitialized && !mapDucks[1].alive && mapDucks[1].bursts==1;
                sawRespawn|=sawInitial && sawHit && mapDucks[0].alive && mapDucks[1].alive;
                if(age>5)break;
            }
            if(host && age>6)break;
        }
        WaitTime(1.f/60);
    }
    if(host)assert(popped && restored && mapDucks[0].alive && mapDucks[1].alive);
    else assert(sawInitial && sawHit && sawRespawn);
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();
    puts(host?"Mallard LAN host passed":"Mallard LAN client: late join, feather event and respawn passed");
}
static void captures(void) {
    fixture(true);duck_renderer_init(NULL);weaponRenderingReady=weapon_renderer_init();character_renderer_init();
    RenderTexture2D rt=LoadRenderTexture(960,640);MapDuck *d=&mapDucks[0];
    Camera3D camera={v_add(d->pos,(Vector3){2,1.5f,-2.8f}),v_add(d->pos,(Vector3){0,.35f,0}),{0,1,0},44,CAMERA_PERSPECTIVE};
    for(int frame=0;frame<50;frame++) {
        if(frame==15)duck_pop(0);
        if(frame>15)d->respawn=90-(frame-15)/20.f;
        BeginTextureMode(rt);ClearBackground((Color){64,79,90,255});BeginMode3D(camera);
        draw_world_surfaces();DrawVoxels(camera);ducks_draw();EndMode3D();EndTextureMode();
        Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/duck/map-%03d.png",frame));UnloadImage(im);
    }
    UnloadRenderTexture(rt);duck_renderer_shutdown();weapon_renderer_shutdown();character_renderer_shutdown();
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,640,"Map ducks");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    const char *role=getenv("FPS_DUCK_LAN");
    if(role)lan(!strcmp(role,"host"));else {placements();hits();weapons();network();captures();}
    fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Map mallards passed");return 0;
}
