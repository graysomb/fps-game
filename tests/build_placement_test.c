#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(goliathStates,0,sizeof(goliathStates));
 worldPeriodic=false;worldFloorHalfSize=40;activePlayers=6;firefightHumans=2;gameMode=GAME_MODE_DEATHMATCH;gameState=GAME_STATE_PLAYING;
 enemyLab=equipmentLab=creativeModeActive=firefightPractice=useCustomMap=false;netTransport.role=NET_ROLE_OFFLINE;fp_new_run();eq_reset_all();clarity_reset();combat_ui_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){players[i].matter=players[i].matterMax=100;players[i].pos=(Vector3){20+i*3,BASE_EYE_HEIGHT,20};players[i].last_build_time=players[i].last_shot_time=-100;playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_MEDIUM;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;ai_reset(i);}
 players[0].pos=(Vector3){.25f,BASE_EYE_HEIGHT,.25f};
}

int main(void){
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Building placement");
 for(int upgrade=0;upgrade<2;upgrade++){
  fixture();equipment[0].upgrades=upgrade?EQ_BUILD:0;players[0].pos.y=5;players[0].pitch=0;
  Vector3 at[EQ_PANEL_CUBES];assert(eq_panel_positions(0,at));
  assert(fabsf(v_length(world_delta(at[eq_build_count(0)/2],players[0].pos))-(upgrade?16.0f:8.0f))<VOXEL_SIZE);
  assert(eq_build(0));
  assert(voxel_count==(upgrade?125:27));assert(players[0].matter==90);
  puts("Unsupported midair build: PASS");
  fixture();equipment[0].upgrades=upgrade?EQ_BUILD:0;players[0].pos.y=5;players[0].pitch=-10;
  assert(eq_panel_positions(0,at));assert(v_length(world_delta(at[eq_build_count(0)/2],players[0].pos))>8);
  assert(eq_build(0));puts("Distant ground build: PASS");
  fixture();equipment[0].upgrades=upgrade?EQ_BUILD:0;players[0].pos.y=5;
  add_static_voxel_at_grid(0,10,-40,GRAY,0);rebuild_voxel_hash();
  assert(eq_panel_positions(0,at));assert(v_length(world_delta(at[eq_build_count(0)/2],players[0].pos))>8);assert(eq_build(0));
  puts("Distant surface build: PASS");
  fixture();equipment[0].upgrades=upgrade?EQ_BUILD:0;players[0].pos.y=5;
  assert(eq_panel_positions(0,at));Vector3 edge=at[0];
  int existing=add_static_voxel_at_grid(floorf(edge.x/VOXEL_SIZE),floorf(edge.y/VOXEL_SIZE),floorf(edge.z/VOXEL_SIZE),RED,0);
  voxels[existing].owner=77;
  Vector3 loose=at[1];int debris=addVoxel(loose.x,loose.y,loose.z,false,true,GRAY,0);assert(debris>=0);rebuild_voxel_hash();
  assert(eq_panel_positions(0,at));assert(eq_build(0));assert(voxels[existing].owner==77);assert(voxels[existing].color.r==RED.r);
  assert(voxel_count==eq_build_count(0)+1);assert(players[0].matter==90);
  puts("Terrain and dynamic debris overlap: PASS");
  fixture();equipment[0].upgrades=upgrade?EQ_BUILD:0;players[0].pos.y=5;
  assert(eq_panel_positions(0,at));players[1].pos=at[0];float lift[MAX_COMBATANTS];assert(eq_build_lifts(eq_build_count(0),at,lift));
  Vector3 ceiling=v_add(players[1].pos,(Vector3){0,lift[1],0});
  add_static_voxel_at_grid(floorf(ceiling.x/VOXEL_SIZE),floorf(ceiling.y/VOXEL_SIZE),floorf(ceiling.z/VOXEL_SIZE),GRAY,0);rebuild_voxel_hash();
  assert(player_aabb_hits_world_ex(ceiling,player_collision_half(&players[1]),1));
  assert(eq_build(0));assert(players[0].matter==90);puts("Blocked headroom no longer vetoes build: PASS");
 }
 fixture();CloseWindow();return 0;
}
