#define FIREFIGHT_TESTING
#include "raylib.h"
static double testClock=10;
static double test_time(void){return testClock;}
#define GetTime test_time
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fresh(float matter){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;firefightHumans=1;activePlayers=MAX_COMBATANTS;
 firefightPractice=enemyLab=equipmentLab=creativeModeActive=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;worldPeriodic=false;worldFloorHalfSize=40;netTransport.role=NET_ROLE_OFFLINE;
 fp_new_run();eq_reset_all();clarity_reset();ex_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){players[i].respawn_timer=9999;combatantLife[i]=LIFE_UNUSED;playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_MEDIUM;combatantTeam[i]=i<4?0:1;}
 players[0].pos=(Vector3){0,1,-15};players[0].respawn_timer=0;players[0].matter=100;players[0].invuln_timer=1000;combatantLife[0]=LIFE_ALIVE;
 init_firefight_bot_entity(4,ENEMY_TYPE_LOBBER,INPUT_TYPE_BOT_MEDIUM);players[4].pos=(Vector3){0,1,0};players[4].yaw=players[4].pitch=0;players[4].onGround=true;players[4].matter=matter;players[4].isExposed=matter<=0;
 players[4].last_melee_time=-100;ai_select(4,0,players[0].pos,AI_ALERT_SEEN);aiCombat[4].reaction=0;
}
static void terrain(void){
 for(int n=0;n<4;n++)add_static_box_at_grid(-2+n*10,2+n*10,0,6,-7,-5,GRAY);
 init_static_hash();rebuild_voxel_hash();
}
static void tick(void){testClock+=1.0/60;UpdateBot(4,1.f/60);update_melee_swings();net_integrate_player(&players[4],1.f/60,true);rebuild_voxel_hash();}
static void mining(void){
 for(int lab=0;lab<2;lab++)for(int empty=0;empty<2;empty++){
  fresh(empty?0:52);terrain();enemyLab=lab;enemyLabPaused=false;
  float initial=players[4].matter;int frames=0;
  while((players[4].matter<120 || enemySpecial[4].harvesting || players[4].meleeSwingActive) && frames++<3600)tick();
  fprintf(stderr,"Lobber lab%d exposed%d: %.1f -> %.1f in %.2fs, mining%d, target%llu, reason=%s pos %.1f %.1f %.1f\n",lab,empty,initial,players[4].matter,frames/60.f,enemySpecial[4].harvesting,(unsigned long long)enemySpecial[4].harvestTarget,enemySpecial[4].reason,players[4].pos.x,players[4].pos.y,players[4].pos.z);
  assert(frames<3600 && players[4].matter==120 && !players[4].isExposed && !enemySpecial[4].harvesting);
  assert(voxel_count>0);
 }
 fresh(80);terrain();UpdateBot(4,.1f);assert(!enemySpecial[4].harvesting);
 fresh(10);terrain();netTransport.role=NET_ROLE_CLIENT;UpdateBot(4,1);assert(players[4].matter==10 && !enemySpecial[4].harvesting);netTransport.role=NET_ROLE_OFFLINE;
 enemyLab=true;enemyLabPaused=true;UpdateBot(4,1);assert(!enemySpecial[4].harvesting);enemyLabPaused=false;UpdateBot(4,.1f);assert(enemySpecial[4].harvesting);
 fresh(0);UpdateBot(4,.1f);assert(enemySpecial[4].harvesting && !enemySpecial[4].harvestTarget && !players[4].tetherHolding && players[4].matter==0);
 // Compaction and missing targets are reacquired without retaining stale indices.
 fresh(0);terrain();UpdateBot(4,.1f);int v=enemy_voxel(enemySpecial[4].harvestTarget);assert(v>=0);remove_voxel_index(v);init_static_hash();rebuild_voxel_hash();for(int f=0;f<80;f++)tick();assert(v_isfinite(players[4].pos));
 // Low Matter interrupts an attack wind-up and safely releases carried rubble.
 fresh(10);int cube=addVoxel(0,1,-2,false,true,GRAY,0);players[4].goldTetherCharged=true;acquire_tether_hit(4,cube,(Vector3){0,1,0});assert(players[4].tetherHolding);
 enemySpecial[4].windup=1;enemySpecial[4].state=3;terrain();UpdateBot(4,.1f);assert(!players[4].tetherHolding && enemySpecial[4].windup==0 && enemySpecial[4].state==0);
 puts("Lobber recovery: real melee harvesting to full in Firefight/lab, zero Matter, interruption, lifecycle and authority passed");
}
static void drops(void){
 fresh(120);fp_spawn(4,ENEMY_TYPE_LOBBER);assert(fpCount==1 && !fp_carrier(4) && fpItems[0].type==PICKUP_GOLD_TETHER && fpItems[0].owner==4);
 players[4].isExposed=true;players[4].matter=0;players[4].invuln_timer=0;kill_player(4,0,false,false);assert(fpItems[0].owner==-1);
 assert(fp_collect(0,0,false) && players[0].goldTetherCharged && fpBannerWrite==1);
 init_firefight_bot_entity(5,ENEMY_TYPE_LOBBER,INPUT_TYPE_BOT_MEDIUM);fp_spawn(5,ENEMY_TYPE_LOBBER);assert(fpCount==2 && !fp_carrier(5));
 fp_death(5,players[5].pos);assert(!fp_collect(0,1,false) && fpItems[1].owner==-1);
 players[0].death_pos=players[0].pos;ff_resolve_death(0);assert(fpItems[0].owner==-1 && !players[0].goldTetherCharged);
 puts("Lobber drops: every native lobber, collection, duplicates, unlock and owner death passed");
}
static void cover(void){
 fresh(120);assert(player_max_matter(&players[4])==120);
 add_static_box_at_grid(-8,8,0,6,-12,-10,GRAY);init_static_hash();rebuild_voxel_hash();
 Vector3 move;assert(lobber_cover(4,.1f,&move));
 assert(sapper_hidden(4,enemySpecial[4].coverPosition,players[0].pos,NULL));
 assert(v_isfinite(move));
 clear_world_voxels();init_static_hash();rebuild_voxel_hash();
 assert(!lobber_cover(4,.1f,&move));
 puts("Lobber Matter capacity and cover selection/invalidation passed");
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Lobber harvest tests");
 init_pbd_thread_pool();assert(initialize_physics_backend());mining();drops();cover();enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Lobber harvest passed");return 0;
}
