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
static float accumulator;
static void step(float dt,bool physics){
    NetInputCommand command={.held=players[0].tetherHolding?NET_INPUT_TETHER:0};net_apply_movement_command(&players[0],&command,dt);
    net_integrate_player(&players[0],dt,true);
    if(physics){accumulator+=dt;while(accumulator>=1.f/60-1e-7f){prepare_tether_forces();simulate_voxel_pbd_steps(1.f/60,1);rebuild_voxel_hash();accumulator-=1.f/60;}}
    assert(v_isfinite(players[0].pos)&&v_isfinite(players[0].vel));
}
static Vector3 momentum(void){Vector3 sum=v_mul(players[0].vel,TETHER_GRAPPLE_PLAYER_MASS);sum=v_add(sum,v_mul(players[1].vel,TETHER_GRAPPLE_PLAYER_MASS));
    int stamp=++particle_sync_stamp;
    for(int v=0;v<voxel_count;v++)for(int k=0;k<8;k++){Particle *p=voxels[v].particles[k];if(p&&p->inv_mass>0&&p->sync_stamp!=stamp){p->sync_stamp=stamp;sum=v_add(sum,v_mul(p->vel,fmaxf(1,p->glue_count)/p->inv_mass));}}return sum;
}
static void conservation(void){
    for(int size=1;size<=16;size*=4){
        fixture();int first=-1;
        for(int n=0;n<size;n++){int v=addVoxel((n%4)*.5f,3+(n/4)*.5f,-5,false,true,GRAY,0);glue_neighbor_faces_for_voxel(v);if(first<0)first=v;}
        rebuild_voxel_hash();acquire_tether_hit(0,first,(Vector3){0,0,1});
        Vector3 before=momentum();weak_tether_prepare(.01f);Vector3 impulse=weakTethers[0].impulse;assert(v_length(impulse)>0);
        weak_tether_integrate(.01f);Vector3 after=momentum();assert(v_length(v_sub(after,before))<.001f);
        if(size==1){Vector3 cube={0};for(int k=0;k<8;k++)cube=v_add(cube,v_mul(voxels[first].particles[k]->vel,.125f));float ratio=v_length(cube)/v_length(players[0].vel);printf("Single cube/player speed change ratio %.3f (expected 16)\n",ratio);assert(fabsf(ratio-16)<.001f);}
    }
    fixture();int v=anchor((Vector3){0,3,-5},false);players[1].pos=(Vector3){0,3,-10};acquire_tether_hit(1,v,(Vector3){0,0,1});
    Vector3 before=momentum();weak_tether_prepare(.01f);weak_tether_integrate(.01f);assert(v_length(v_sub(momentum(),before))<.001f);
    // Relative damping removes energy, rather than damping toward world zero.
    fixture();v=anchor((Vector3){0,3,-5},false);WeakTetherCoupling c=weakTethers[0];assert(weak_tether_sample(0,&c));
    players[0].pos=v_sub(c.target,v_sub(player_tether_hold_position(&players[0]),players[0].pos));players[0].vel=(Vector3){2,0,0};
    for(int k=0;k<8;k++)voxels[v].particles[k]->vel=(Vector3){-2,0,0};
    weak_tether_prepare(.001f);assert(weakTethers[0].impulse.x<0);
    Vector3 j=weak_tether_impulse((Vector3){1,0,0},(Vector3){0},.125f,.01f);
    Vector3 moving=weak_tether_impulse((Vector3){1,0,0},v_sub((Vector3){12,0,0},(Vector3){12,0,0}),.125f,.01f);assert(Vector3Distance(j,moving)<1e-5f);
    puts("Isolated momentum, connected/shared particles, two-player pulls and relative damping passed");
}
static void behavior(void){
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);Vector3 initial=players[0].pos;
    for(int view=0;view<4;view++)prepare_tether_forces();assert(v_length(players[0].vel)==0&&Vector3Distance(initial,players[0].pos)==0);
    accumulator=0;for(int n=0;n<60;n++)step(1.f/60,true);
    printf("Anchored grounded grapple travelled %.3fm\n",-players[0].pos.z);assert(players[0].pos.z< -1 && players[0].matter==100);
    fixture();anchor((Vector3){0,6,-3},true);accumulator=0;for(int n=0;n<30;n++)step(1.f/60,true);assert(players[0].pos.y>BASE_EYE_HEIGHT+.5f&&!players[0].onGround&&players[0].grappleAirborne);
    Vector3 velocity=players[0].vel;Vector3 zero={0};release_tether_velocity(0,&zero);assert(!tether_grapple_active(0));
    NetInputCommand none={0};net_apply_movement_command(&players[0],&none,.01f);assert(fabsf(players[0].vel.x-velocity.x)<.001f&&fabsf(players[0].vel.z-velocity.z)<.001f);
    for(int n=0;n<600;n++){net_apply_movement_command(&players[0],&none,1.f/60);net_integrate_player(&players[0],1.f/60,false);}
    assert(players[0].onGround&&!players[0].grappleAirborne&&fabsf(players[0].vel.x)+fabsf(players[0].vel.z)<.001f);
    fixture();int id=anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    players[0].goldTetherHolding=true;assert(!tether_grapple_active(0));players[0].goldTetherHolding=false;players[0].drainTarget=2;assert(!tether_grapple_active(0));players[0].drainTarget=0;
    combatantLife[0]=LIFE_DOWNED;assert(!tether_grapple_active(0));combatantLife[0]=LIFE_ALIVE;
    players[0].respawn_timer=1;assert(!tether_grapple_active(0));players[0].respawn_timer=0;
    playerInput[0]=INPUT_TYPE_BOT_HARD;assert(!tether_grapple_active(0));playerInput[0]=INPUT_TYPE_KEYBOARD;remove_voxel_index(id);assert(!tether_grapple_active(0));
    fixture();players[0].goldTetherCharged=true;players[0].matter=0;players[0].isExposed=true;anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);
    assert(!players[0].goldTetherHolding&&tether_grapple_active(0));weak_tether_prepare(.01f);assert(v_length(players[0].vel)>0&&players[0].matter==0);
    fixture();worldPeriodic=true;players[0].pos.x=39;anchor((Vector3){-37,BASE_EYE_HEIGHT,0},true);weak_tether_prepare(.1f);assert(players[0].vel.x>0);
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);add_static_box_at_grid(-4,4,0,10,-5,-5,GRAY);init_static_hash();rebuild_voxel_hash();accumulator=0;for(int n=0;n<120;n++)step(1.f/60,true);assert(players[0].pos.z>=-1.751f);
    puts("Anchored motion, airborne carry, landing brakes, exclusions, lifecycle and collisions passed");
}
static void moving_voxels(void){
    fixture();int id=anchor((Vector3){0,BASE_EYE_HEIGHT,-5},false);uint64_t identity=voxels[id].identity;Vector3 start=voxels[id].pos;accumulator=0;
    for(int n=0;n<12;n++)step(1.f/120,true);assert(resolve_tether_voxel(&players[0])&&players[0].tetherVoxelIdentity==identity);assert(voxels[players[0].tetherVoxel].pos.z>start.z+.02f);
    assert(fabsf(players[0].pos.z)<.1f);players[0].vel=(Vector3){0,0,-4};NetInputCommand none={.held=NET_INPUT_TETHER};net_apply_movement_command(&players[0],&none,.02f);assert(v_length(players[0].vel)==0);
    for(int n=0;n<180;n++)step(1.f/60,true);Vector3 settled=players[0].pos;
    for(int n=0;n<60;n++)step(1.f/60,true);printf("Loose cube hold drift over last second %.4fm\n",Vector3Distance(settled,players[0].pos));assert(Vector3Distance(settled,players[0].pos)<.15f);
    fixture();add_static_box_at_grid(-2,2,0,15,-12,-10,GRAY);init_static_hash();rebuild_voxel_hash();start_tether(0);assert(players[0].tetherHolding&&tether_grapple_active(0));accumulator=0;for(int n=0;n<12;n++)step(1.f/120,true);
    puts("Loose cube stays lightweight, ground braking, settled drift and terrain activation passed");
}
static void rates_and_prediction(void){
    Vector3 positions[3];int rates[]={30,60,144};
    for(int r=0;r<3;r++){fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);accumulator=0;for(int n=0;n<rates[r];n++)step(1.f/rates[r],true);positions[r]=players[0].pos;}
    printf("30 vs 144 FPS grapple difference %.3fm\n",Vector3Distance(positions[0],positions[2]));assert(Vector3Distance(positions[0],positions[2])<.5f);
    fixture();anchor((Vector3){0,BASE_EYE_HEIGHT,-6},true);WeakTetherCoupling c=weakTethers[0];assert(weak_tether_sample(0,&c));weak_tether_prepare(.01f);Vector3 hostVelocity=players[0].vel;players[0].vel=(Vector3){0};
    netTransport.role=NET_ROLE_CLIENT;players[0].grappleInputHeld=true;tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER,c.target,c.velocity,0,1);tether_grapple_pull(0,.01f);assert(Vector3Distance(hostVelocity,players[0].vel)<.0001f);
    NetInputCommand release={0};net_apply_movement_command(&players[0],&release,.01f);assert(!tether_grapple_active(0));
    players[0].grappleInputHeld=true;tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER|NET_PLAYER_VISUAL_GOLD_TETHER,c.target,(Vector3){0},0,2);assert(!tether_grapple_active(0));
    tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER|NET_PLAYER_VISUAL_DRAIN,c.target,(Vector3){0},0,3);assert(!tether_grapple_active(0));
    tether_grapple_receive(0,NET_PLAYER_VISUAL_TETHER,c.target,(Vector3){1,0,0},.125f,4);assert(players[0].netGrappleAnchorVelocity.x==1&&players[0].netGrappleInverseMass==.125f);
    eq_reset(0);assert(!tether_grapple_active(0)&&!players[0].grappleAirborne);netTransport.role=NET_ROLE_OFFLINE;
    puts("30/60/144 Hz, shared client prediction, explicit target velocity/mass and resets passed");
}
static void physical_momentum(void){
    assert(GRAVITY==0 && VELOCITY_DAMPING==1);
    for(int owners=1;owners<=2;owners++)for(int size=1;size<=16;size*=4){
        if(owners==2&&size!=4)continue;
        fixture();int first=-1;for(int n=0;n<size;n++){int v=addVoxel((n%4)*.5f,6+(n/4)*.5f,-5,false,true,GRAY,0);glue_neighbor_faces_for_voxel(v);if(first<0)first=v;}
        rebuild_voxel_hash();acquire_tether_hit(0,first,(Vector3){0,0,1});
        if(owners==2){players[1].pos=(Vector3){0,6,-10};acquire_tether_hit(1,first,(Vector3){0,0,1});}
        prepare_tether_forces();Vector3 before=momentum();simulate_voxel_pbd_steps(1.f/60,1);Vector3 after=momentum();
        float exchanged=v_length(weakTethers[0].impulse)+v_length(weakTethers[1].impulse),error=Vector3Distance(before,after)/fmaxf(exchanged,.001f);
        fprintf(stderr,"Physical momentum cells=%d players=%d backend=%s error=%.5f%%\n",size,owners,physics_backend_name(physicsBackend.active),100*error);assert(error<.01f);
        if(size==1){float ratio=v_length(voxels[first].vel)/v_length(players[0].vel);assert(fabsf(ratio-16)<.1f);}
    }
}
#ifndef TETHER_GRAPPLE_LAN
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Grapple tests");init_pbd_thread_pool();assert(initialize_physics_backend());
#ifdef TETHER_COUPLING_ISOLATED
    physical_momentum();
#else
    conservation();behavior();moving_voxels();rates_and_prediction();
#endif
    clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Default tether grapple passed");return 0;
}


#endif
