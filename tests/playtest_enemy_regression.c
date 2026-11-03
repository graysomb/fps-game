#define FIREFIGHT_TESTING
#include "raylib.h"
static double testAdvance;
static double enemy_test_time(void){return GetTime()+testAdvance;}
#define GetTime enemy_test_time
#define main fps_application_main
#include "../fps_ray.c"
#undef main
bool ffTestReviveHeld[4];
static void fresh(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();init_static_hash();
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;firefightHumans=1;activePlayers=MAX_COMBATANTS;
 firefightPractice=enemyLab=equipmentLab=creativeModeActive=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;worldPeriodic=false;worldFloorHalfSize=40;netTransport.role=NET_ROLE_OFFLINE;
 fp_new_run();eq_reset_all();clarity_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){players[i].respawn_timer=9999;combatantLife[i]=LIFE_UNUSED;playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_EASY;combatantTeam[i]=i<4?0:1;}
 players[0].pos=(Vector3){-.9f,1,-12};players[0].respawn_timer=0;players[0].matter=100;combatantLife[0]=LIFE_ALIVE;
 init_firefight_bot_entity(4,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_EASY);players[4].pos=(Vector3){-.9f,1,0};players[4].yaw=0;players[4].pitch=0;
 init_firefight_bot_entity(5,ENEMY_TYPE_SHIELDBEARER,INPUT_TYPE_BOT_EASY);players[5].pos=(Vector3){0,1,-3};players[5].yaw=0;players[5].pitch=0;enemy_update_equipment(5);rebuild_voxel_hash();
 ai_select(4,0,players[0].pos,AI_ALERT_SEEN);aiCombat[4].reaction=0;
}

static int failures;
static void check(bool ok,const char *name){printf("%s: %s\n",name,ok?"PASS":"REPRODUCED");failures+=!ok;}
static void behaviors(void){
 fresh();players[5].respawn_timer=9999;enemy_release(5);
 fpCarrierType[4]=PICKUP_GLOB;equipment[4].gun[0]=EQ_GLOB;players[4].matter=0;players[4].isExposed=true;
 players[0].respawn_timer=9999;ai_reset(4);add_static_box_at_grid(-2,-1,2,3,-3,-3,GRAY);init_static_hash();rebuild_voxel_hash();
 for(int n=0;n<10;n++)UpdateBot(4,.1f);
 check(players[4].meleeSwingActive,"Carrier mines with no visible enemy");
 fresh();players[0].pos=(Vector3){0,1,-20};players[4].pos=(Vector3){-1,1,0};enemy_release(5);init_firefight_bot_entity(5,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_EASY);players[5].pos=(Vector3){1,1,0};
 fpCarrierType[4]=PICKUP_GLOB;fpCarrierType[5]=PICKUP_GLOCK;equipment[4].gun[0]=EQ_GLOB;equipment[5].gun[0]=EQ_GLOCK;
 clear_world_voxels();init_static_hash();rebuild_voxel_hash();fp_pair_link(4,5);fp_pairs_tick(.1f);holograms[0].life=10;holograms[0].body=players[0];holograms[0].body.pos=(Vector3){0,1,-7};
 for(int n=0;n<30;n++)fp_pairs_tick(.1f);
 check(fpPairs[4].target==MAX_COMBATANTS,"Pair retargets nearby visible hologram");holograms[0].body.pos.z=-20;fp_pairs_tick(.1f);fp_pair_move(4,.1f);check(players[4].vel.z<0,"Pair follows moving hologram");holograms[0].life=0;
 fresh();players[4].respawn_timer=9999;players[0].pos=(Vector3){0,1,-4.9f};players[0].invuln_timer=0;ai_select(5,0,players[0].pos,AI_ALERT_SEEN);
 for(int n=0;n<20;n++)enemy_update(5,.1f);
 check(players[0].matter<100,"Shield bearer melee across own shield");
 for(int gun=EQ_GLOB;gun<=EQ_DMR;gun++){
  fresh();players[0].pos=(Vector3){-.9f,1,-9};Ray ray={players[0].pos,{0,0,1}};ex_reset();ex_fire(0,gun,ray);if(gun==EQ_GLOB)ex_tick(.2f);
  int hits=0;for(int g=0;g<4;g++)hits+=SHIELD_SECTION_HITS-enemySpecial[5].shieldHits[g];
  check(hits>0,gun==EQ_GLOB?"Glob damages shield":gun==EQ_GLOCK?"Glock damages shield":"DMR damages shield");
 }
}
static void harvesting(void){
 for(int gun=EQ_GLOB;gun<=EQ_GLOCK;gun++){
  fresh();enemy_cleanup(true);clear_world_voxels();players[5].respawn_timer=9999;players[0].respawn_timer=9999;ai_reset(4);
  fpCarrierType[4]=gun==EQ_GLOB?PICKUP_GLOB:PICKUP_GLOCK;equipment[4].gun[0]=gun;players[4].matter=0;players[4].isExposed=true;players[4].last_melee_time=-100;
  for(int n=0;n<4;n++)add_static_box_at_grid(-2+n*10,2+n*10,0,6,-7,-5,GRAY);init_static_hash();rebuild_voxel_hash();
  int f=0;for(;f<3600 && players[4].matter<player_max_matter(&players[4]);f++){
   testAdvance+=1.0/60;UpdateBot(4,1.f/60);update_melee_swings();net_integrate_player(&players[4],1.f/60,true);rebuild_voxel_hash();
  }
  printf("Carrier %d harvested %.1f Matter in %.2fs\n",gun,players[4].matter,f/60.f);
  check(players[4].matter==player_max_matter(&players[4]),"Carrier actually harvests to full");
 }
}
static void profile(void){
 enemy_cleanup(true);gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightHumans=1;ResetGame();start_firefight_wave(12);update_firefight_logic(3);assert(goActive);
 for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;
 for(int f=0;f<30;f++)simulate_voxel_pbd_steps(1.f/60,1);
 goRegen=true;go_release_part(GOLEM_ARM_R);goClock+=goRepairDelay-.1f;go_repair_prepare();
 check(goCells[goMount[GOLEM_ARM_R]].index<0,"Golem waits four seconds before regrowth");goClock=10;
 double sum=0,peak=0;int before=voxel_count;
 for(int n=0;n<40;n++){double start=GetTime();go_repair_prepare();double ms=(GetTime()-start)*1000;sum+=ms;peak=fmax(peak,ms);}
 printf("Foundry repair check: world=%d cells=%d mean=%.3fms peak=%.3fms\n",before,goCount,sum/40,peak);
 // The combat stance can intersect nearby cover: profile successful publication
 // separately after clearing only those non-golem obstacles in the work area.
 for(int v=voxel_count-1;v>=0;v--)if(!voxels[v].golemOwned){
  Vector3 b[8];go_voxel_points(v,b);bool overlap=false;
  for(int c=0;c<goCount&&!overlap;c++)if(goCells[c].index<0&&go_repair_group(goCells[c].part)==GOLEM_ARM_R){Vector3 a[8];go_cell_points(c,a);overlap=go_boxes_overlap(a,b);}
  if(overlap)remove_voxel_index(v);
 }
 for(int c=0;c<goCount;c++)goCells[c].index=go_find(goCells[c].id,goCells[c].index);
 init_static_hash();rebuild_voxel_hash();
 sum=peak=0;double commit=0;
 for(int n=0;n<120;n++){
  bool missing=goCells[goMount[GOLEM_ARM_R]].index<0;goClock+=1.f/60;
  double start=GetTime();go_repair_prepare();double ms=(GetTime()-start)*1000;sum+=ms;peak=fmax(peak,ms);
  if(missing&&goCells[goMount[GOLEM_ARM_R]].index>=0)commit=ms;
 }
 printf("Foundry repair timeline: mean=%.3fms peak=%.3fms commit=%.3fms\n",sum/120,peak,commit);
 check(goCells[goMount[GOLEM_ARM_R]].index>=0,"Foundry arm rebuild completes");
}
int main(int argc,char **argv){SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Enemy playtest regressions");setvbuf(stdout,NULL,_IONBF,0);assert(parse_physics_arguments(argc,argv));init_pbd_thread_pool();assert(initialize_physics_backend());behaviors();harvesting();profile();golem_physics_reset();enemy_cleanup(true);clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return failures!=0;}
