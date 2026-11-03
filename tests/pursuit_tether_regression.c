#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fresh(void) {
    enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
    memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(holograms,0,sizeof(holograms));
    memset(fpPairs,0,sizeof(fpPairs));memset(fpCarrierType,0,sizeof(fpCarrierType));
    gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;firefightHumans=1;activePlayers=MAX_COMBATANTS;
    firefightPractice=enemyLab=enemyLabPaused=creativeModeActive=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;
    worldPeriodic=false;worldFloorHalfSize=40;netTransport.role=NET_ROLE_OFFLINE;
    for(int i=0;i<MAX_COMBATANTS;i++) {
        playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i?1:0;
        players[i].respawn_timer=i?9999:0;players[i].matter=players[i].matterMax=100;
        players[i].pos=(Vector3){0,BASE_EYE_HEIGHT,-8};players[i].onGround=true;players[i].invuln_timer=1000;ai_reset(i);
    }
    eq_reset_all();clarity_reset();
    init_firefight_bot_entity(4,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_MEDIUM);
    players[4].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[4].yaw=0;
}
static void wall(void) {add_static_box_at_grid(-12,12,0,8,-6,-4,GRAY);init_static_hash();rebuild_voxel_hash();}
static void held(void) {
    for(int gold=0;gold<2;gold++) {
        fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[0].goldTetherCharged=gold;
        int v=addVoxel(0,2,-4,false,true,GRAY,0);assert(v>=0);
        int neighbor=addVoxel(VOXEL_SIZE,2,-4,false,true,GRAY,0);assert(neighbor>=0);
        init_static_hash();rebuild_voxel_hash();rebuild_glue_constraints();
        acquire_tether_hit(0,v,(Vector3){0,0,1});assert(players[0].tetherHolding);prepare_tether_forces();
        uint64_t id=voxels[v].identity,neighborId=voxels[neighbor].identity;
        voxels[v].lifeFrames=RECYCLE_DYNAMIC_MAX_FRAMES+1;
        voxels[neighbor].lifeFrames=RECYCLE_DYNAMIC_MAX_FRAMES+1;
        voxels[neighbor].restorationQueued=true;voxels[neighbor].restorationDeadlineFrame=recycleFrameCounter+1;
        voxels[v].restorationQueued=true;voxels[v].restorationDeadlineFrame=recycleFrameCounter+1;
        for(int f=0;f<120;f++)recycle_dead_voxels(1.f/60);
        assert(players[0].tetherHolding && voxels[players[0].tetherVoxel].identity==id);
        assert(voxels[players[0].tetherVoxel].simulate);
        bool connected=false;for(int n=0;n<voxel_count;n++)connected|=voxels[n].identity==neighborId;assert(connected);
        release_tether(0);recycle_dead_voxels(1.f/60);
        bool found=false;for(int n=0;n<voxel_count;n++)found|=voxels[n].identity==id;assert(!found);
    }
    puts("Held voxel survives lifetime/deadline; released voxel recycles: PASS");
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};
    int old=addVoxel(8.25f,2.25f,-4.25f,false,true,GRAY,0);
    int held=addVoxel(.25f,2.25f,-4.25f,false,true,GRAY,0);
    init_static_hash();rebuild_voxel_hash();acquire_tether_hit(0,held,(Vector3){0,0,1});prepare_tether_forces();
    uint64_t identity=voxels[held].identity;
    voxels[old].lifeFrames=RECYCLE_DYNAMIC_MAX_FRAMES+1;
    voxels[held].restorationQueued=true;voxels[held].restorationDeadlineFrame=recycleFrameCounter+1;
    recycle_dead_voxels(1.f/60);
    assert(players[0].tetherHolding && enemy_voxel(identity)>=0 && voxels[enemy_voxel(identity)].simulate);
    puts("Held protection survives other-island restoration/compaction: PASS");
}
static void memory(void) {
    for(int holo=0;holo<2;holo++) {
        fresh();wall();int target=holo?MAX_COMBATANTS:0;
        if(holo){holograms[0].life=20;holograms[0].body=players[0];players[0].respawn_timer=9999;}
        ai_select(4,target,(Vector3){0,BASE_EYE_HEIGHT,-8},AI_ALERT_SEEN);
        for(int n=0;n<68;n++)ai_tick(4,.1f);
        assert(aiCombat[4].target==target && !aiCombat[4].visible);
        for(int n=0;n<4;n++)ai_tick(4,.1f);
        assert(aiCombat[4].target==-1);
    }
    puts("Player and hologram pursuit lasts seven seconds: PASS");
}
static void routes(void) {
    fresh();wall();ai_select(4,0,players[0].pos,AI_ALERT_SEEN);
    float maxSide=0;double start=GetTime();
    for(int n=0;n<360;n++) {
        UpdateBot(4,1.f/60);net_integrate_player(&players[4],1.f/60,true);
        maxSide=fmaxf(maxSide,fabsf(players[4].pos.x));assert(v_isfinite(players[4].pos));
    }
    printf("Wall route: position %.2f %.2f, detour %.2f, simulation %.2fms\n",players[4].pos.x,players[4].pos.z,maxSide,(GetTime()-start)*1000);
    assert(maxSide>6 && aiCombat[4].visible);
    // Direct pursuit reaches the goal, not just the first corner.
    for(int n=0;n<360;n++) {
        Vector3 delta=world_delta(players[0].pos,players[4].pos);delta.y=0;
        if(v_length(delta)<.4f)break;
        bot_apply_move(&players[4],&botStates[4],delta,MOVE_SPEED,1.f/60);
        net_integrate_player(&players[4],1.f/60,true);
    }
    assert(v_length(world_delta(players[0].pos,players[4].pos))<.5f);
    puts("Rifleman routes around wall and resumes sight; route reaches destination: PASS");
    fresh();worldPeriodic=true;players[4].pos=(Vector3){39,BASE_EYE_HEIGHT,0};
    Vector3 goal={-36,BASE_EYE_HEIGHT,0};
    add_static_box_at_grid(79,79,0,8,-4,4,GRAY);init_static_hash();rebuild_voxel_hash();
    for(int n=0;n<240;n++){
        Vector3 delta=world_delta(goal,players[4].pos);delta.y=0;if(v_length(delta)<.4f)break;
        bot_apply_move(&players[4],&botStates[4],delta,MOVE_SPEED,1.f/60);net_integrate_player(&players[4],1.f/60,true);
    }
    assert(v_length(world_delta(goal,players[4].pos))<.5f);
    puts("Route crosses periodic seam around obstruction: PASS");
}
static void arrival(void) {
    fresh();wall();Hologram *h=&holograms[0];h->body=players[0];h->life=9;h->health=50;h->travelRemaining=.1f;h->direction=(Vector3){1,0,0};
    uint32_t before=clarityNextId;hologram_tick(.1f);
    assert(h->arrivalAnnounced && clarityNextId==before+1);
    assert(aiCombat[4].target==MAX_COMBATANTS && aiCombat[4].memory==AI_PURSUIT_SECONDS);
    for(int n=0;n<10;n++)hologram_tick(.1f);
    assert(clarityNextId==before+1);
    puts("Arrival squeaks once and attracts enemies through cover: PASS");
    fresh();wall();init_firefight_bot_entity(5,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_MEDIUM);
    players[5].pos=(Vector3){1,BASE_EYE_HEIGHT,0};
    fpCarrierType[4]=PICKUP_GLOB;fpCarrierType[5]=PICKUP_GLOCK;
    equipment[4].gun[0]=EQ_GLOB;equipment[5].gun[0]=EQ_GLOCK;fp_pair_link(4,5);
    h=&holograms[0];h->body=players[0];h->life=10;h->health=50;
    hologram_arrived(0);fp_pairs_tick(.1f);fp_update_carrier(4,.1f);
    assert(fpPairs[4].target==MAX_COMBATANTS && fpPairs[5].target==MAX_COMBATANTS);
    assert(aiCombat[4].target==MAX_COMBATANTS && v_length(players[4].vel)>0);
    puts("Paired carriers pursue heard decoy inside normal firing range: PASS");
}
static void contact_and_sound(void) {
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[4].pos=(Vector3){0,BASE_EYE_HEIGHT,-4};
    int v=addVoxel(0,BASE_EYE_HEIGHT,-4,false,true,GRAY,0);init_static_hash();rebuild_voxel_hash();
    acquire_tether_hit(0,v,(Vector3){0,0,1});prepare_tether_forces();uint64_t id=voxels[v].identity;
    players[4].invuln_timer=0;players[4].matter=players[4].matterMax=1000;
    handle_pbd_projectile_hits();assert(players[4].matter<1000);assert(players[0].tetherHolding && enemy_voxel(id)>=0);
    float matter=players[4].matter;handle_pbd_projectile_hits();assert(players[4].matter==matter);
    puts("Held contact preserves voxel and cannot damage every frame: PASS");
    fresh();players[0].pos=(Vector3){-20,BASE_EYE_HEIGHT,-20};
    EnemyStructure tower={.tower=true,.tag=999,.origin={15,0,15}};
    v=addVoxel(-20,2,-23,false,true,GRAY,0);voxels[v].towerTag=tower.tag;
    init_static_hash();rebuild_voxel_hash();acquire_tether_hit(0,v,(Vector3){0,0,1});prepare_tether_forces();
    id=voxels[v].identity;assert(!tower_place(&tower));assert(players[0].tetherHolding && enemy_voxel(id)>=0);
    release_tether(0);assert(tower_place(&tower));
    puts("Tower repair waits for held piece and resumes after release: PASS");
    Wave wave=make_hologram_squeak();assert(wave.data && wave.frameCount==10584);
    short *samples=wave.data;int peak=0;for(unsigned n=0;n<wave.frameCount;n++)peak=abs(samples[n])>peak?abs(samples[n]):peak;
    assert(peak>1000 && peak<32000 && samples[0]==0);UnloadWave(wave);
    uint8_t packet[64];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
    uint32_t event=clarityLastEvent+1;net_write_u32(&w,event);net_write_u8(&w,4);net_write_u8(&w,0);net_write_u8(&w,0);
    for(int n=0;n<4;n++)net_write_f32(&w,0);
    NetReader r;net_reader_init(&r,packet,w.length);assert(clarity_net_receive(NET_MSG_GAME_EVENT,&r,0));assert(clarityLastEvent==event);
    net_reader_init(&r,packet,w.length);clarity_net_receive(NET_MSG_GAME_EVENT,&r,0);assert(clarityLastEvent==event);
    puts("Squeak waveform and reliable event decoding/deduplication: PASS");
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Pursuit regressions");setvbuf(stdout,NULL,_IONBF,0);
    held();memory();routes();arrival();contact_and_sound();enemy_cleanup(true);clear_world_voxels();CloseWindow();return 0;
}

