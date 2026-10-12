#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void arena(void) {
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;enemyLab=false;
    multiplayerActivePlayers=1;multiplayerPlayerInput[0]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;
    for(int i=4;i<MAX_COMBATANTS;i++)players[i].respawn_timer=9999;
    players[0].pos=(Vector3){0,3,0};players[0].invuln_timer=1000;
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,720,"Firefight playtest checks");SetTargetFPS(0);
    init_pbd_thread_pool();assert(initialize_physics_backend());
    arena();assert(FLOOR_SIZE==40);
    int scrap=0;for(int v=0;v<voxel_count;v++)scrap+=voxels[v].owner==ENEMY_SCRAP_OWNER;
    fprintf(stderr,"Arena: %d voxels, %d supply voxels\n",voxel_count,scrap);assert(scrap==64);
    int before=voxel_count;enemy_stock_scrap(false);assert(voxel_count==before);
    for(int k=4;k<8;k++) {
        assert(pickups[k].active && pickups[k].type==(k>=6?PICKUP_GOLD_TETHER:PICKUP_DYNAMIC_SHOT));
        players[0].pos=pickups[k].pos;players[0].goldTetherCharged=players[0].dynamicShotActive=false;
        update_pickups(.01f);assert(!pickups[k].active);
        assert(k>=6?players[0].goldTetherCharged:players[0].dynamicShotActive);
        players[0].pos=(Vector3){0,3,0};update_pickups(31);assert(pickups[k].active);
        pickups[k].active=false;pickups[k].respawnTimer=20;
    }
    ff_wave_complete();for(int k=4;k<8;k++)assert(pickups[k].active);
    // Real terrain, normal AI/movement/tether/physics: no forced ability state.
    EnemyType types[]={ENEMY_TYPE_BUILDER,ENEMY_TYPE_LOBBER};
    for(int n=0;n<2;n++) {
        arena();players[4].enemyType=types[n];Vector3 spawn;assert(ff_spawn_position(4,&spawn));
        enemy_lab_spawn(4,types[n],spawn);
        assert(spawn_position_clear_ex(players[4].pos,spawn_clear_half(4),0,4));
        for(int frame=0;frame<2400 && !enemySpecial[4].actions;frame++) {
            enemy_autonomous_step(1.0f/60);
            if(frame%600==599)fprintf(stderr,"%s elapsed=%d state=%d held=%d selected=%llu pos=%.1f,%.1f,%.1f reason=%s\n",enemy_names[types[n]],frame/60,enemySpecial[4].state,players[4].tetherHolding,(unsigned long long)enemySpecial[4].selected,players[4].pos.x,players[4].pos.y,players[4].pos.z,enemySpecial[4].reason);
        }
        EnemySpecialState *s=&enemySpecial[4];
        fprintf(stderr,"Arena %s: actions=%d pos=(%.1f %.1f %.1f) reason=%s\n",enemy_names[types[n]],s->actions,players[4].pos.x,players[4].pos.y,players[4].pos.z,s->reason);
        assert(s->actions>0);
    }
    // Untouched supplies survive the normal ten-second recycle lifetime and remain usable.
    enemy_cleanup(true);clear_world_voxels();enemy_scrap_bundle((Vector3){0,0,0});
    uint64_t supply=voxels[0].identity;
    prepare_tether_forces();simulate_voxel_pbd_steps(1.0f/120,2);
    assert(voxels[0].sleeping && voxels[0].vgs_active);
    assert(build_glue_cluster_indices(0,glueClusterIndices)==8);
    for(int v=0;v<voxel_count;v++)voxels[v].lifeFrames=RECYCLE_DYNAMIC_MAX_FRAMES+10;
    recycle_dead_voxels();assert(enemy_voxel(supply)>=0 && voxels[enemy_voxel(supply)].simulate);
    assert(wake_sleeping_cluster_in_place(enemy_voxel(supply)));
    // A charge dislodges a substantial region and accelerates surviving debris.
    clear_world_voxels();for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;
    add_static_box_at_grid(-6,6,0,7,-2,2,GRAY);init_static_hash();
    enemy_blast((EnemyCharge){.owner=4,.pos={0,1,0}});
    int moving=0;for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && v_length(voxels[v].vel)>.1f) {
        assert(v_dot(voxels[v].vel,v_sub(voxels[v].pos,(Vector3){0,1,0}))>-.01f);moving++;
    }
    assert(moving>8);prepare_tether_forces();simulate_voxel_pbd_steps(1.0f/120,2);
    for(int p=0;p<sim_particle_count;p++)assert(v_isfinite(sim_particles[p]->pos));
    if(physics_backend_is_gpu(physicsBackend.requested))assert(physicsBackend.active==physicsBackend.requested);
    enemy_cleanup(true);shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Firefight playtest checks passed");return 0;
}
