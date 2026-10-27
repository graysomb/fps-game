#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void arena(int difficulty) {
    enemyLab=false;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=false;
    currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=2;ResetGame();
    enemy_cleanup(true);clear_world_voxels();init_static_hash();
    gameState=GAME_STATE_PLAYING;activePlayers=6;
    for(int i=0;i<activePlayers;i++){players[i].respawn_timer=9999;ai_reset(i);}
    for(int i=0;i<2;i++){
        players[i].respawn_timer=0;players[i].pos=(Vector3){i*5,2,-20};
        playerInput[i]=INPUT_TYPE_KEYBOARD;combatantTeam[i]=0;
    }
    enemy_lab_spawn(4,ENEMY_TYPE_STANDARD,(Vector3){0,2,0});
    playerInput[4]=difficulty;players[4].yaw=0;players[4].pitch=0;
    players[4].matter=10000;players[4].invuln_timer=0;
    players[1].respawn_timer=9999;ai_reset(4);
}
static void swarmer_checks(void) {
    arena(INPUT_TYPE_BOT_EASY);worldPeriodic=false;players[0].pos=(Vector3){0,2,-5};
    activePlayers=7;
    enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,(Vector3){0,2,0});
    enemy_lab_spawn(5,ENEMY_TYPE_SWARMER,(Vector3){12,2,0});
    enemy_lab_spawn(6,ENEMY_TYPE_SWARMER,(Vector3){29,2,0});
    players[4].yaw=0;players[5].yaw=players[6].yaw=180;
    ai_tick(4,.1f);
    assert(aiCombat[4].visible && aiCombat[4].target==0);
    assert(aiCombat[5].target==0 && aiCombat[5].source==AI_ALERT_HEARD);
    assert(!aiCombat[5].visible && aiCombat[6].target==-1);
    ai_swarmer_rally(5);assert(aiCombat[6].target==-1);
    Vector3 known=aiCombat[5].known;players[0].pos.x=3;
    ai_tick(4,.1f);assert(v_length(v_sub(known,aiCombat[5].known))<.001f);
    ai_tick(4,1);assert(v_length(v_sub(players[0].pos,aiCombat[5].known))<.001f);
    aiCombat[4].visible=false;aiCombat[4].rallyCooldown=0;
    players[0].pos.x=7;ai_swarmer_rally(4);assert(aiCombat[5].known.x==3);
    arena(INPUT_TYPE_BOT_EASY);firefightWave=2;firefightWaveEnemiesRemaining=999;
    players[0].matter=100;players[0].isExposed=false;players[0].invuln_timer=0;
    for(int hit=1;hit<=5;hit++) {
        enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,players[0].pos);
        players[4].invuln_timer=0;players[4].contactDamageTimer=0;
        update_enemy_contact_damage(.1f);
        assert(players[4].respawn_timer>0);
        if(hit<=4){assert(players[0].matter==100-25*hit);assert(players[0].respawn_timer<=0);}
    }
    assert(combatantLife[0]==LIFE_DOWNED);
    arena(INPUT_TYPE_BOT_EASY);enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,players[0].pos);
    players[0].invuln_timer=1;players[0].matter=100;players[0].isExposed=false;
    update_enemy_contact_damage(.1f);assert(players[0].matter==100 && players[4].respawn_timer<=0);
    netTransport.role=NET_ROLE_CLIENT;players[0].invuln_timer=0;
    update_enemy_contact_damage(.1f);assert(players[0].matter==100);netTransport.role=NET_ROLE_OFFLINE;
    fprintf(stderr,"Swarmer 25 damage, exposed finish, invulnerability, authority and local rally passed\n");
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);
    SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(640,480,"AI combat regression");init_pbd_thread_pool();assert(initialize_physics_backend());
    swarmer_checks();
    int levels[]={INPUT_TYPE_BOT_EASY,INPUT_TYPE_BOT_MEDIUM,INPUT_TYPE_BOT_HARD};
    int rates[]={30,60,120};
    for(int d=0;d<3;d++)for(int r=0;r<3;r++){
        arena(levels[d]);double last=-100,first=-1;int shots=0;
        for(int f=0;f<rates[r]*8;f++){
            float dt=1.f/rates[r];ai_tick(4,dt);ai_aim(4,aiCombat[4].known,dt);
            if(ai_ready_to_fire(4)){
                assert(FireVoxel(4));double now=aiCombat[4].clock;
                assert(now-last>=1./3.-.00001);
                if(shots && shots%(d==0?2:3)==0)assert(now-last>=ai_pause(4)-.00001);
                if(first<0)first=now;last=now;shots++;
                assert(!FireVoxel(4));
            }
        }
        assert(first>=ai_reaction(4) && shots>=6 && shots<20);
        fprintf(stderr,"AI difficulty %d @ %d Hz: %d shots, first %.3f s\n",d,rates[r],shots,first);
        ai_tick(4,2);ai_aim(4,aiCombat[4].known,.01f);
        if(ai_ready_to_fire(4)){assert(FireVoxel(4));assert(!FireVoxel(4));}
    }
    arena(INPUT_TYPE_BOT_MEDIUM);players[0].pos.z=12;
    assert(!check_ai_awareness_of_player(4,0,NULL));
    ai_alert(4,0,players[0].pos,AI_ALERT_HIT);ai_tick(4,.01f);
    assert(aiCombat[4].target==0 && aiCombat[4].visible);
    players[1].respawn_timer=0;players[1].pos=(Vector3){0,2,-5};
    ai_tick(4,.1f);assert(aiCombat[4].target==0);
    ai_tick(4,1.5f);assert(aiCombat[4].target==1);
    ai_alert(4,0,players[0].pos,AI_ALERT_HIT);assert(aiCombat[4].target==0);
    ai_alert(4,1,players[1].pos,AI_ALERT_HIT);assert(aiCombat[4].target==0);
    arena(INPUT_TYPE_BOT_MEDIUM);ai_tick(4,.1f);
    Vector3 known=aiCombat[4].known;players[0].pos=(Vector3){0,2,-100};
    ai_tick(4,1);assert(!aiCombat[4].visible && !ai_ready_to_fire(4));
    assert(v_length(v_sub(known,aiCombat[4].known))<.001f);
    ai_tick(4,3.1f);assert(aiCombat[4].target==-1);
    arena(INPUT_TYPE_BOT_MEDIUM);players[0].pos=(Vector3){0,2,20};
    ai_gunfire(0,players[0].pos);assert(aiCombat[4].target==0);
    ai_reset(4);players[0].pos.z=25;ai_gunfire(0,players[0].pos);assert(aiCombat[4].target==-1);
    int bullet=addVoxel(1,2,-5,false,true,WHITE,0);
    voxels[bullet].isBullet=true;voxels[bullet].owner=0;voxels[bullet].shotOrigin=players[0].pos;
    ai_near_miss(bullet,(Vector3){1,2,-5},(Vector3){1,2,5});assert(aiCombat[4].target==0);
    ai_reset(4);ai_near_miss(bullet,(Vector3){2,2,-5},(Vector3){2,2,5});assert(aiCombat[4].target==-1);
    netTransport.role=NET_ROLE_CLIENT;ai_alert(4,0,players[0].pos,AI_ALERT_HIT);assert(aiCombat[4].target==-1);netTransport.role=NET_ROLE_OFFLINE;
    arena(INPUT_TYPE_BOT_HARD);ai_tick(4,.1f);ai_tick(4,1);
    ai_aim(4,aiCombat[4].known,.01f);assert(fabsf(aiCombat[4].error-1)<.001f);
    players[4].vel.x=3;ai_aim(4,aiCombat[4].known,.01f);assert(fabsf(aiCombat[4].error-2.5f)<.001f);
    ai_alert(4,0,players[0].pos,AI_ALERT_HIT);assert(aiCombat[4].tracking==0);
    ai_reset(0);assert(aiCombat[4].target==-1);
    arena(INPUT_TYPE_BOT_MEDIUM);
    add_static_box_at_grid(-8,8,0,8,-10,-10,GRAY);init_static_hash();
    ai_gunfire(0,players[0].pos);assert(aiCombat[4].target==-1);
    players[0].pos.z=-7;ai_gunfire(0,players[0].pos);assert(aiCombat[4].target==0);
    ai_tick(4,.1f);assert(!aiCombat[4].visible && !ai_ready_to_fire(4));
    ai_reset(4);bullet=addVoxel(0,2,-10,false,true,WHITE,0);
    voxels[bullet].isBullet=true;voxels[bullet].owner=0;voxels[bullet].shotOrigin=players[0].pos;
    ai_near_miss(bullet,(Vector3){0,2,-10},(Vector3){0,2,3});assert(aiCombat[4].target==-1);
    arena(INPUT_TYPE_BOT_MEDIUM);enemy_lab_spawn(4,ENEMY_TYPE_LOBBER,(Vector3){0,2,0});
    ai_alert(4,0,players[0].pos,AI_ALERT_HIT);enemySpecial[4].windup=.5f;
    players[1].respawn_timer=0;aiCombat[4].overrideDelay=0;
    ai_alert(4,1,players[1].pos,AI_ALERT_HIT);assert(aiCombat[4].target==0 && aiCombat[4].pendingTarget==1);
    enemySpecial[4].windup=0;ai_tick(4,.01f);assert(aiCombat[4].target==1);
    int hits[2]={0},attempts[2]={0};
    for(int moving=0;moving<2;moving++){
        arena(INPUT_TYPE_BOT_HARD);players[0].pos.z=-30;
        for(int frame=0;frame<1800;frame++){
            float phase=fmodf(frame/60.f,4);players[0].pos.x=moving?(phase<2?-6+phase*6:18-phase*6):0;
            players[0].vel.x=moving?(phase<2?6:-6):0;
            ai_tick(4,1.f/60);ai_aim(4,aiCombat[4].known,1.f/60);
            if(ai_ready_to_fire(4)){
                Vector3 direction=player_forward(&players[4]);
                float flight=(players[0].pos.z-players[4].pos.z)/(direction.z*60);
                Vector3 impact=v_add(players[4].pos,v_mul(direction,60*flight));
                Vector3 target=v_add(players[0].pos,v_mul(players[0].vel,flight));
                attempts[moving]++;hits[moving]+=fabsf(impact.x-target.x)<PLAYER_SIZE*.5f && fabsf(impact.y-target.y)<.5f;
                ai_shot_fired(4);
            }
        }
    }
    fprintf(stderr,"30m aim trial: stationary %d/%d, strafe %d/%d (ballistic projection)\n",hits[0],attempts[0],hits[1],attempts[1]);
    assert(hits[0]>hits[1] && attempts[0]>20 && attempts[1]>20);
    for(int fixture=8;fixture<15;fixture++){enemy_lab_fixture(fixture);assert(enemyLabFixture==fixture);enemyLabPaused=false;enemy_lab_ai_step(.1f);UpdateBot(4,.1f);}
    fprintf(stderr,"AI combat regression PASS\n");shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}

