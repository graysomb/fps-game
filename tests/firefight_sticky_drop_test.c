#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(int humans){
    netTransport.role=NET_ROLE_OFFLINE;gameMode=GAME_MODE_FIREFIGHT;
    firefightPractice=enemyLab=equipmentLab=false;firefightFoundry=true;firefightEndless=true;firefightCheckpoints=false;
    useCustomMap=false;multiplayerActivePlayers=humans;
    for(int i=0;i<humans;i++)multiplayerPlayerInput[i]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;
    firefightWaveEnemiesRemaining=firefightWaveEnemiesTotal=100;
}
static size_t collect_carrier(int type,int actor){
    size_t n=fpCount;fpPendingCount=1;fpPending[0]=fp_rule(type);
    init_firefight_bot_entity(actor,(EnemyType)fpRules[fpPending[0]].chassis,INPUT_TYPE_BOT_EASY);
    fp_spawn(actor,fpRules[fpPending[0]].chassis);assert(fpCount==n+1);
    players[actor].pos=players[0].pos;players[actor].invuln_timer=0;
    kill_player(actor,0,false,false);assert(fpItems[n].owner==-1);
    assert(fp_collect(0,n,false));assert(fpItems[n].owner==0);return n;
}
static void solo(bool loaded){
    fixture(1);size_t sticky=collect_carrier(PICKUP_ADHESIVE,4);
    size_t gun=0,build=0;
    if(loaded){gun=collect_carrier(PICKUP_GLOCK,5);build=collect_carrier(PICKUP_REINFORCED,6);}
    assert(eq_upgrade(0,EQ_ADHESIVE));
    grenadeInventory[0]=2;assert(clarity_throw(0));grenadeCooldown[0]=0;assert(clarity_throw(0));
    assert(!grenadeInventory[0] && eq_upgrade(0,EQ_ADHESIVE) && fpItems[sticky].owner==0);
    players[0].pos=(Vector3){0,3,3}; // Die away from the carriers' dropped physical armor.
    players[0].invuln_timer=0;Vector3 death=players[0].pos;kill_player(0,7,false,false);
    assert(combatantLife[0]==LIFE_WAITING && !eq_upgrade(0,EQ_ADHESIVE));
    assert(fpItems[sticky].owner==-1);
    assert(fpCount==(loaded?3u:1u)); // kill_player resolves death twice; never duplicate loot.
    if(loaded)assert(fpItems[gun].owner==-1 && fpItems[build].owner==-1);
    if(loaded){
        assert(v_length(world_delta(fpItems[sticky].pos,fpItems[gun].pos))>.5f);
        assert(v_length(world_delta(fpItems[sticky].pos,fpItems[build].pos))>.5f);
        assert(v_length(world_delta(fpItems[gun].pos,fpItems[build].pos))>.5f);
    }else{Vector3 offset=world_delta(fpItems[sticky].pos,death);assert(fabsf(offset.x)+fabsf(offset.z)<.001f);}
    assert(v_length(world_delta(fpItems[sticky].pos,death))<1.4f);
    fp_tick(.01f);assert(fpItems[sticky].owner==-1);
    ff_update_coop(PLAYER_RESPAWN_TIME+.1f);assert(combatantLife[0]==LIFE_ALIVE);
    players[0].pos=v_add(fpItems[sticky].pos,(Vector3){0,.2f,0});
    fp_tick(.01f);
    assert(fpItems[sticky].owner==0 && eq_upgrade(0,EQ_ADHESIVE));
    puts(loaded?"Solo: sticky drops separately from guns/build upgrades and can be recovered":"Solo: sticky-only death drop remains centered and can be recovered after using both grenades");
}
static void coop(void){
    fixture(2);size_t sticky=collect_carrier(PICKUP_ADHESIVE,4);
    players[0].invuln_timer=0;kill_player(0,7,false,false);
    assert(combatantLife[0]==LIFE_DOWNED && fpItems[sticky].owner==0 && eq_upgrade(0,EQ_ADHESIVE));
    ff_restore_human(0,35);assert(fpItems[sticky].owner==0 && eq_upgrade(0,EQ_ADHESIVE));
    players[0].invuln_timer=players[1].invuln_timer=0;kill_player(0,7,false,false);kill_player(1,7,false,false);
    ff_update_coop(.01f);assert(fpItems[sticky].owner==-1 && !eq_upgrade(0,EQ_ADHESIVE));
    puts("Co-op: retains sticky while downed/revived, drops it on resolved team-wipe death passed");
}
static void checkpoint(void){
    fixture(1);size_t sticky=collect_carrier(PICKUP_ADHESIVE,4);uint32_t id=fpItems[sticky].id;
    firefightCheckpoints=true;start_firefight_wave(5);fp_restore();
    assert(fpItems[sticky].id==id && fpItems[sticky].owner==0 && eq_upgrade(0,EQ_ADHESIVE));
    players[0].pos=(Vector3){0,3,3};players[0].invuln_timer=0;kill_player(0,7,false,false);
    assert(fpItems[sticky].id==id && fpItems[sticky].owner==-1 && !eq_upgrade(0,EQ_ADHESIVE));
    puts("Sticky ownership and death drops survive wave-checkpoint restoration");
}
int main(void){
    SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Sticky death-drop regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());solo(false);solo(true);coop();checkpoint();
    fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
