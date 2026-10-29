#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
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
static void formations(void){
 fresh();shield_escort_tick();assert(shieldEscort[4]==6);assert(shield_escort_clear_shot(4));
 shield_escort_move(4,.1f);assert(players[4].vel.z<0);
 players[4].pos=(Vector3){-.9f,1,-1.6f};shield_escort_move(4,.1f);assert(players[4].vel.x==0 && players[4].vel.z==0);
 players[5].yaw=90;shield_escort_move(4,.1f);assert(v_length(players[4].vel)>0);
 for(int i=6;i<=7;i++){init_firefight_bot_entity(i,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_EASY);players[i].pos=(Vector3){2,1,0};}
 shield_escort_tick();assert(shieldEscort[6]==6 && !shieldEscort[7]);
 players[4].pos.x=0;players[4].yaw=0;assert(!shield_escort_clear_shot(4));
 players[5].respawn_timer=9;shield_escort_tick();assert(!shieldEscort[4] && !shieldEscort[6]);
 players[5].respawn_timer=0;shield_escort_tick();assert(shieldEscort[4]);ai_reset(5);assert(!shieldEscort[4]);
 shield_escort_tick();enemy_release(5);shield_escort_tick();assert(!shieldEscort[4]);
 puts("Formation: follow, stop, turn, two-rifle limit, ally body, death/reset/shield loss passed");
}
static void shots(void){
 fresh();VoxelHit hit;Ray ray={players[4].pos,{0,0,-1}};
 assert(first_voxel_hit_detailed(ray,12,-1,false,&hit));assert(clarity_shield_voxel(hit.id));
 assert(!rifle_shield_ray(ray,12,-1,4,&hit));assert(!rifle_view_blocked(4,players[0].pos));
 assert(rifle_shield_ray(ray,12,-1,0,&hit)); // human cannot phase through it
 int bullet=addVoxel(-.9f,1,-.8f,false,true,WHITE,0);voxels[bullet].isBullet=true;voxels[bullet].owner=4;voxels[bullet].vel=(Vector3){0,0,-60};
 rebuild_voxel_hash();update_projectiles(.2f);assert(players[0].matter<100);
 for(int g=0;g<4;g++)assert(enemySpecial[5].shieldHits[g]==SHIELD_SECTION_HITS);
 fresh();bullet=addVoxel(-.9f,1,-9,false,true,WHITE,0);voxels[bullet].isBullet=true;voxels[bullet].owner=0;voxels[bullet].vel=(Vector3){0,0,60};
 rebuild_voxel_hash();update_projectiles(.1f);
 int hits=0;for(int g=0;g<4;g++)hits+=SHIELD_SECTION_HITS-enemySpecial[5].shieldHits[g];assert(hits==1);
 fresh();add_static_box_at_grid(-4,4,0,5,-15,-14,GRAY);init_static_hash();rebuild_voxel_hash();
 assert(rifle_shield_ray(ray,12,-1,4,&hit));assert(!voxels[hit.id].simulate);assert(rifle_view_blocked(4,players[0].pos));
  fresh();assert(ex_fire(4,EQ_GLOCK,ray));assert(players[0].matter<100);
 for(int g=0;g<4;g++)assert(enemySpecial[5].shieldHits[g]==SHIELD_SECTION_HITS);
 fresh();ex_reset();assert(ex_fire(4,EQ_GLOB,ray));ex_tick(.21f);assert(players[0].matter<100);
 for(int g=0;g<4;g++)assert(enemySpecial[5].shieldHits[g]==SHIELD_SECTION_HITS);
 fresh();combatantTeam[5]=0;assert(rifle_shield_ray(ray,12,-1,4,&hit));
 puts("Shots: actual friendly projectile passage, player interception, health, walls and hostile shields passed");
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Shield escort tests");
 init_pbd_thread_pool();assert(initialize_physics_backend());formations();shots();enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Shield escort tests passed");return 0;
}