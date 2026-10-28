#define TOWER_PREVIEW
#include "builder_towers_test.c"
static void setup(void){
 enemyLab=false;gameState=GAME_STATE_PLAYING;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;firefightHumans=1;activePlayers=MAX_COMBATANTS;firefightEndless=false;ResetGame();start_firefight_wave(12);
 players[0].invuln_timer=10000;update_firefight_logic(3);
 assert(goActive&&goCombatEnabled&&!goCombatDead);assert(firefightEnemiesSpawned==1);assert(ffRemaining[ENEMY_TYPE_FOUNDRY]==0);
 assert(players[4].enemyType==ENEMY_TYPE_STANDARD&&go_combat_actor(4));
 assert(goCount>0 && go_move_clear(goFight.position));
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Firefight golem regression");init_pbd_thread_pool();assert(initialize_physics_backend());
 setup();Vector3 initial=goFight.position;
 players[0].pos=v_add(initial,(Vector3){initial.x>0?-9:9,1,0});
 for(int f=0;f<1800;f++){go_combat_prepare(1.f/60);golem_physics_tick(1.f/60);go_combat_commit(1.f/60);}
 fprintf(stderr,"Golem cells=%d bullets=%d rockets=%d movement=%.2f\n",goCount,goArsenal.bullets,goArsenal.rockets,Vector3Distance(initial,goFight.position));
 assert(goArsenal.bullets>0&&goArsenal.rockets>0);
 int remaining=firefightWaveEnemiesRemaining,kills=players[0].kills;
 eq_damage(0,4,1000,true,DEATH_RAIL,players[0].pos);
 assert(goCombatDead&&players[0].kills==kills+1&&firefightWaveEnemiesRemaining==remaining-1);
 go_combat_damage_from(1000,false,0);assert(players[0].kills==kills+1&&firefightWaveEnemiesRemaining==remaining-1);
 for(int f=0;f<100 && firefightWaveState==FIREFIGHT_STATE_ACTIVE;f++){
  update_firefight_logic(1);
  assert(go_combat_actor(4)&&players[4].respawn_timer>0);
  for(int i=5;i<activePlayers;i++)if(players[i].respawn_timer<=0){players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,0,true,false);}
 }
 assert(firefightWaveState==FIREFIGHT_STATE_VICTORY&&firefightWaveEnemiesRemaining==0);
 setup();start_firefight_wave(1);assert(!goActive&&!goCombatEnabled);for(int v=0;v<voxel_count;v++)assert(!voxels[v].golemOwned);
 enemy_cleanup(true);clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Firefight golem spawn, attacks, kill credit, reserves, victory and reset passed");return 0;
}
