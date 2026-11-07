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
static void equipment_grab(EnemyType type,bool gold){
    enemy_lab_start(type);update_goliath_armor_positions(4);
    players[0].pos=v_add(players[4].pos,(Vector3){0,0,-4});
    players[0].yaw=180;players[0].pitch=0;players[0].goldTetherCharged=gold;
    rebuild_voxel_hash();float matter=players[4].matter;
    start_tether(0);
    assert(players[0].tetherHolding && !players[0].drainTarget);
    assert(get_goliath_owner_of_voxel(players[0].tetherVoxel)<0);
    update_tether_drains(1);assert(players[4].matter==matter);release_tether(0);
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,600,"Physical tether checks");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    for(int mode=0;mode<2;mode++)for(int gold=0;gold<2;gold++){
        setup();gameMode=mode?GAME_MODE_FIREFIGHT:GAME_MODE_DEATHMATCH;
        combatantTeam[0]=0;combatantTeam[1]=1;players[0].goldTetherCharged=gold;
        start_tether(0);assert(!players[0].drainTarget && !players[0].tetherHolding);
        players[0].drainTarget=2; // Stale locks cannot continue dealing damage.
        update_tether_drains(10);assert(!players[0].drainTarget && players[1].matter==100);
        assert(players[0].matter==0);
        uint8_t data[1024];NetWriter writer;NetReader reader;NetPlayerWireState state;int slot;
        net_writer_init(&writer,data,sizeof(data));net_write_player_state(&writer,0);
        net_reader_init(&reader,data,writer.length);assert(net_read_player_state(&reader,&slot,&state));
        assert(!(state.visual.flags&NET_PLAYER_VISUAL_DRAIN));
    }
    setup();int v=addVoxel(0,2,-3,false,true,GRAY,0);assert(v>=0);rebuild_voxel_hash();
    start_tether(0);assert(players[0].tetherHolding && !players[0].drainTarget);release_tether(0);
    for(int gold=0;gold<2;gold++){
        equipment_grab(ENEMY_TYPE_GOLIATH,gold);
        equipment_grab(ENEMY_TYPE_SHIELDBEARER,gold);
    }
    shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Physical tether / no actor drain checks passed");return 0;
}
