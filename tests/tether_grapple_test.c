#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];

static void fixture(void){
    enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
    memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));
    gameMode=GAME_MODE_DEATHMATCH;gameState=GAME_STATE_PLAYING;activePlayers=2;
    enemyLab=enemyLabPaused=creativeModeActive=false;worldPeriodic=false;worldFloorHalfSize=40;netTransport.role=NET_ROLE_OFFLINE;
    for(int i=0;i<MAX_COMBATANTS;i++){
        playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;
        players[i].matter=players[i].matterMax=100;players[i].pos=(Vector3){i*20,BASE_EYE_HEIGHT,0};players[i].onGround=true;
    }
    eq_reset_all();clarity_reset();
}
static int anchor(Vector3 pos,bool fixed){
    int v=addVoxel(pos.x,pos.y,pos.z,false,true,GRAY,0);assert(v>=0);
    if(fixed){voxels[v].fixed=true;for(int c=0;c<8;c++)voxels[v].particles[c]->inv_mass=0;}
    init_static_hash();rebuild_voxel_hash();acquire_tether_hit(0,v,(Vector3){0,0,1});assert(players[0].tetherHolding);return v;
}
static void step(float dt,bool physics){
    NetInputCommand command={.held=NET_INPUT_TETHER};net_apply_movement_command(&players[0],&command,dt);
    net_integrate_player(&players[0],dt,true);
    if(physics){prepare_tether_forces();simulate_voxel_pbd_steps(dt,1);rebuild_voxel_hash();}
    assert(v_isfinite(players[0].pos)&&v_isfinite(players[0].vel));
}
static void behavior(void){
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    Vector3 initial=players[0].pos;
    for(int view=0;view<4;view++)prepare_tether_forces();
    assert(v_length(players[0].vel)==0 && v_length(v_sub(players[0].pos,initial))==0);
    tether_grapple_pull(0,1);assert(v_length(players[0].vel)<=TETHER_GRAPPLE_MAX_APPROACH+.001f);
    players[0].vel=(Vector3){0};
    for(int n=0;n<30;n++)step(1.f/60,false);
    assert(players[0].pos.z< -1 && players[0].matter==100);
    players[0].vel=(Vector3){0,0,-4};HandleKeyboardInput(0,.01f);assert(players[0].vel.z==-4);
    HandleGamepadInput(0,.01f);assert(players[0].vel.z==-4);
    release_tether(0);assert(!tether_grapple_active(0));
    Vector3 velocity=players[0].vel;tether_grapple_pull(0,.1f);assert(v_length(v_sub(velocity,players[0].vel))==0);
    fixture();anchor((Vector3){0,6,-3},true);for(int n=0;n<30;n++)step(1.f/60,false);
    assert(players[0].pos.y>BASE_EYE_HEIGHT+.5f && !players[0].onGround);
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-1},true);tether_grapple_pull(0,.1f);assert(v_length(players[0].vel)==0);
    fixture();int id=anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    players[0].goldTetherHolding=true;assert(!tether_grapple_active(0));tether_grapple_pull(0,.1f);assert(v_length(players[0].vel)==0);
    players[0].goldTetherHolding=false;players[0].drainTarget=2;assert(!tether_grapple_active(0));players[0].drainTarget=0;
    combatantLife[0]=LIFE_DOWNED;assert(!tether_grapple_active(0));combatantLife[0]=LIFE_ALIVE;
    players[0].respawn_timer=1;assert(!tether_grapple_active(0));players[0].respawn_timer=0;
    playerInput[0]=INPUT_TYPE_BOT_HARD;assert(!tether_grapple_active(0));playerInput[0]=INPUT_TYPE_KEYBOARD;
    remove_voxel_index(id);assert(!tether_grapple_active(0));
    fixture();players[0].goldTetherCharged=true;players[0].matter=0;anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    assert(!players[0].goldTetherHolding && tether_grapple_active(0));
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    add_static_box_at_grid(-4,4,0,10,-5,-5,GRAY);init_static_hash();rebuild_voxel_hash();
    for(int n=0;n<120;n++)step(1.f/60,false);assert(players[0].pos.z>=-1.751f);
    fixture();worldPeriodic=true;players[0].pos.x=39;anchor((Vector3){-37,BASE_EYE_HEIGHT,0},true);
    tether_grapple_pull(0,.1f);assert(players[0].vel.x>0 && fabsf(players[0].vel.z)<players[0].vel.x*.3f);
    puts("Grapple controls, upward lift, exclusions, lifecycle, collisions and periodic distance passed");
}
static void moving_voxels(void){
    fixture();int id=anchor((Vector3){0,BASE_EYE_HEIGHT,-5},false);uint64_t identity=voxels[id].identity;
    Vector3 start=voxels[id].pos;
    for(int n=0;n<12;n++)step(1.f/120,true);
    assert(resolve_tether_voxel(&players[0]) && players[0].tetherVoxelIdentity==identity);
    assert(players[0].pos.z<-.02f && voxels[players[0].tetherVoxel].pos.z>start.z+.02f);
    // Supported terrain is acquired through the unchanged activation path.
    fixture();add_static_box_at_grid(-2,2,0,15,-12,-10,GRAY);init_static_hash();rebuild_voxel_hash();start_tether(0);
    assert(players[0].tetherHolding && tether_grapple_active(0));for(int n=0;n<12;n++)step(1.f/120,true);
    assert(players[0].pos.z<-.02f);
    puts("Actual voxel physics: both endpoints move; supported structure activation passed");
}
static void rates_and_prediction(void){
    Vector3 positions[3];int rates[]={30,60,144};
    for(int r=0;r<3;r++){fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
        for(int n=0;n<rates[r]/2;n++)step(1.f/rates[r],false);positions[r]=players[0].pos;}
    assert(v_length(v_sub(positions[0],positions[2]))<.25f);
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);Vector3 target,vel;assert(tether_grapple_anchor(0,&target,&vel));
    tether_grapple_pull(0,.01f);Vector3 hostVelocity=players[0].vel;players[0].vel=(Vector3){0};
    netTransport.role=NET_ROLE_CLIENT;players[0].grappleInputHeld=true;tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER,target,1);
    tether_grapple_pull(0,.01f);assert(v_length(v_sub(hostVelocity,players[0].vel))<.0001f);
    NetInputCommand release={0};net_apply_movement_command(&players[0],&release,.01f);assert(!tether_grapple_active(0));
    players[0].grappleInputHeld=true;tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER|NET_PLAYER_VISUAL_GOLD_TETHER,target,2);assert(!tether_grapple_active(0));
    tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER|NET_PLAYER_VISUAL_DRAIN,target,3);assert(!tether_grapple_active(0));
    tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER,target,4);
    tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER,v_add(target,(Vector3){.1f,0,0}),10);
    assert(fabsf(players[0].netGrappleAnchorVelocity.x-1)<.001f);
    eq_reset(0);assert(!tether_grapple_active(0));netTransport.role=NET_ROLE_OFFLINE;
    puts("30/60/144 Hz movement, client prediction, release and anchor velocity passed");
}
#ifndef TETHER_GRAPPLE_LAN
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Grapple tests");init_pbd_thread_pool();assert(initialize_physics_backend());
    behavior();moving_voxels();rates_and_prediction();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Default tether grapple passed");return 0;
}


#endif
