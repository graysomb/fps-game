#define FIREFIGHT_TESTING
#include "raylib.h"
static bool keys[512],buttons[32],pad;
static bool test_pressed(int key){return key>=0&&key<512&&keys[key];}
static bool connected(int i){return pad&&i==0;}
static bool test_button(int i,int b){return connected(i)&&b>=0&&b<32&&buttons[b];}
#define IsKeyPressed test_pressed
#define IsGamepadAvailable connected
#define IsGamepadButtonPressed test_button
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();worldPeriodic=false;
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(goliathStates,0,sizeof(goliathStates));
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;activePlayers=MAX_COMBATANTS;firefightHumans=2;firefightPractice=enemyLab=equipmentLab=creativeModeActive=false;firefightCheckpoints=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;
 fp_new_run();eq_reset_all();clarity_reset();combat_ui_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_EASY;combatantLife[i]=i<2?LIFE_ALIVE:LIFE_UNUSED;combatantTeam[i]=i<4?0:1;players[i].respawn_timer=i<2?0:9999;players[i].matter=players[i].matterMax=100;players[i].pos=(Vector3){i*3,1,0};players[i].last_shot_time=players[i].last_melee_time=players[i].last_build_time=-100;ai_reset(i);}
}
static size_t carrier(int type,int slot){int r=fp_rule(type);assert(r>=0);fpPendingCount=1;fpPending[0]=r;players[slot].respawn_timer=0;combatantLife[slot]=LIFE_ALIVE;players[slot].enemyType=fpRules[r].chassis;fp_spawn(slot,fpRules[r].chassis);assert(fp_carrier(slot));return fpCount-1;}
static void schedule(void){
 for(int humans=1;humans<=4;humans++){fixture();firefightHumans=humans;
  for(int wave=1;wave<=36;wave++){
   int before[10],sum=0;for(int t=0;t<10;t++)sum+=(before[t]=ffRemaining[t]=ff_count(wave,humans,t));
   fp_prepare_wave(wave);int after=0;for(int t=0;t<10;t++){after+=ffRemaining[t];assert(ffRemaining[t]>=0);if(t>=2)assert(ffRemaining[t]>=before[t]);}assert(sum==after);
   int introduced=0;for(int r=0;r<FP_RULE_COUNT;r++)if(fpRules[r].wave==wave){int found=0;for(int n=0;n<fpPendingCount;n++)found+=fpPending[n]==r;assert(found==1);introduced++;}
   assert(fpPendingCount==introduced+(wave>=3?(humans+1)/2:0));for(int n=0;n<fpPendingCount;n++)assert(fpRules[fpPending[n]].wave<=wave);
  }
 }
 puts("Schedules: waves 1-36, 1-4 players, guaranteed introductions and preserved quotas passed");
}
static void opening_rounds(void){
 fixture();firefightFoundry=false;firefightHumans=1;start_firefight_wave(1);
 for(int round=1;round<=3;round++){
  int count=6<<(round-1);assert(firefightRound==round && firefightWaveEnemiesTotal==count);
  for(int step=0;step<500 && firefightEnemiesSpawned<count;step++)update_firefight_logic(.1f);
  int living=0;for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++)living+=players[i].respawn_timer<=0;
  assert(living==count && firefightEnemiesSpawned==count);
  for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++)if(players[i].respawn_timer<=0){players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,0,false,false);}
  if(round<3){assert(firefightWaveState==FIREFIGHT_STATE_ACTIVE && ffSpawnTimer==5);update_firefight_logic(1);assert(firefightEnemiesSpawned==0);}
 }
 assert(firefightWaveState==FIREFIGHT_STATE_LOOT);
 start_firefight_wave(1);assert(firefightRound==1 && firefightWaveEnemiesTotal==6);
 start_firefight_wave(2);assert(firefightRound==1);ff_wave_complete();assert(firefightRound==2 && firefightWaveState==FIREFIGHT_STATE_ACTIVE);
 puts("Opening rounds: 6/12/24 simultaneous swarmers, pauses, loot and reset passed");
}
static void second_wave_rounds(void){
 for(int humans=1;humans<=4;humans++){
  fixture();firefightHumans=humans;firefightFoundry=false;firefightTeamLives=3;start_firefight_wave(2);
  for(int round=1;round<=3;round++){
   int total=0;for(int t=0;t<FF_TYPE_COUNT;t++)total+=ff_round_count(2,round,humans,t);
   assert(firefightRound==round && firefightWaveEnemiesTotal==total);
   int watchdog=1000;
   while(firefightWaveState==FIREFIGHT_STATE_ACTIVE && firefightRound==round && watchdog-->0){
    update_firefight_logic(.5f);
    for(int i=4;i<MAX_COMBATANTS;i++)if(players[i].respawn_timer<=0){
     players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,0,false,false);
    }
   }
   assert(watchdog>0);
   if(round<3){assert(firefightWaveState==FIREFIGHT_STATE_ACTIVE && ffSpawnTimer==5 && firefightTeamLives==3);update_firefight_logic(1);assert(firefightEnemiesSpawned==0);}
  }
  assert(firefightWaveState==FIREFIGHT_STATE_LOOT && firefightTeamLives==4);
  int glob=0,build=0;for(size_t n=0;n<fpCount;n++){glob+=fpItems[n].type==PICKUP_GLOB;build+=fpItems[n].type==PICKUP_REINFORCED;}
  assert(glob==1 && build==1);
 }
 puts("Wave 2: three escalating rounds, 1-4 players, pauses, one reward and no duplicate unlock carriers passed");
}
static void ownership(void){
 fixture();size_t glob=carrier(PICKUP_GLOB,4);uint32_t id=fpItems[glob].id;fp_death(4,(Vector3){0,1,0});assert(fpItems[glob].owner==-1&&fpItems[glob].id==id);
 assert(fp_collect(0,glob,false)&&eq_gun(0)==EQ_GLOB);assert(fpUnlocked==(1u<<fp_rule(PICKUP_GLOB)));int banners=fpBannerWrite;assert(!fp_collect(1,glob,false));
 size_t duplicate=carrier(PICKUP_GLOB,5);fp_death(5,players[0].pos);assert(!fp_collect(0,duplicate,true)&&fpItems[duplicate].owner==-1);
 size_t rifle=carrier(PICKUP_DMR,6);fp_death(6,players[0].pos);assert(!fp_collect(0,rifle,false));assert(fp_collect(0,rifle,true));assert(fpItems[glob].owner==-1&&fpItems[rifle].owner==0&&fpCount==3);
 assert(fp_collect(1,glob,false)&&fpItems[glob].id==id&&fpBannerWrite==banners+1);
 size_t build=carrier(PICKUP_REINFORCED,7);fp_death(7,players[0].pos);assert(fp_collect(0,build,false));
 size_t gold=carrier(PICKUP_GOLD_TETHER,8);fp_death(8,players[0].pos);assert(fp_collect(0,gold,false));
 players[0].invuln_timer=0;players[0].matter=0;players[0].isExposed=true;kill_player(0,4,false,false);assert(combatantLife[0]==LIFE_DOWNED&&fpItems[rifle].owner==0&&eq_upgrade(0,EQ_BUILD));
 ff_restore_human(0,35);assert(fpItems[gold].owner==0&&players[0].goldTetherCharged&&eq_gun(0)==EQ_DMR);
 players[0].death_pos=players[0].pos;ff_resolve_death(0);assert(fpItems[rifle].owner==-1&&fpItems[build].owner==-1&&fpItems[gold].owner==-1&&!players[0].goldTetherCharged&&eq_gun(0)==EQ_DEFAULT);
 ff_restore_human(0,100);assert(grenadeInventory[0]==0);
 size_t grenade=carrier(PICKUP_GRENADE,9);assert(grenadeInventory[9]==2);fp_death(9,players[0].pos);assert(fp_collect(0,grenade,false)&&grenadeInventory[0]==1&&fpItems[grenade].owner==-2);
 firefightWave=2;firefightEndless=true;ff_wave_complete();assert(grenadeInventory[0]==1&&grenadeInventory[1]==0&&fpCount==6);
 puts("Conservation, replacement, duplicates, down/revive, death and consumable drops passed");
}
static void enemy_weapons(void){
 Player preview={0};assert(player_weapon_melee_sample(&preview,0).sword_direction==0&&player_weapon_kind(&preview)==WEAPON_STANDARD);
 fixture();players[1].pos.x=30;players[0].pos=(Vector3){0,2,-3};carrier(PICKUP_RAILGUN,15);players[15].pos=(Vector3){0,2,0};players[15].yaw=players[15].pitch=0;
 NetInputCommand held={.held=NET_INPUT_FIRE},released={0};eq_command(15,&held,.65f);eq_command(15,&released,.01f);assert(players[0].matter==50&&players[15].matter==85);eq_tick(2);eq_command(15,&held,.65f);eq_command(15,&released,.01f);assert(combatantLife[0]==LIFE_DOWNED);
 fixture();players[1].pos.x=30;carrier(PICKUP_ROCKET,15);players[0].pos=(Vector3){0,2,-3};eq_rocket_blast((EquipmentRocket){.owner=15,.pos={0,2,-3}});assert(combatantLife[0]==LIFE_DOWNED);
 fixture();players[1].pos.x=30;carrier(PICKUP_GLOB,15);players[15].pos=(Vector3){0,2,0};players[0].pos=(Vector3){0,2,-3};assert(eq_fire(15));ex_tick(.1f);assert(players[0].matter<100&&ex_suppressed(0));
 fixture();carrier(PICKUP_ADHESIVE,15);players[15].pos=(Vector3){0,2,0};assert(clarity_throw(15));int thrown=0;for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id&&grenades[n].owner==15){assert(grenades[n].adhesive);thrown++;}assert(thrown==1&&grenadeInventory[15]==1);
 fixture();carrier(PICKUP_IGLOO,15);players[0].pos.x=players[1].pos.x=players[15].pos.x=20;assert(ex_igloo((GrenadeState){.owner=15,.pos={.25f,.15f,.25f},.launchOrigin={0,1,-4}}));assert(voxel_count>100);
 fixture();carrier(PICKUP_REINFORCED,15);players[15].pos=(Vector3){0,1,0};players[15].pitch=-20;players[0].pos.x=players[1].pos.x=20;assert(eq_build(15)&&voxel_count==125&&players[15].matter==90);
 puts("Enemy slot 15: shared railgun/rocket lethality, Glob suppression, sticky throw, igloo and building passed");
}
static void autonomous_carriers(void){
 for(int r=0;r<FP_RULE_COUNT;r++){
  fixture();int type=fpRules[r].type;carrier(type,15);players[15].pos=(Vector3){0,1,0};players[15].yaw=players[15].pitch=0;players[0].pos=(Vector3){0,1,type==PICKUP_SWORD?-3:-8};players[1].pos=(Vector3){30,1,30};players[0].invuln_timer=100;
  if(type==PICKUP_GOLD_TETHER){add_static_box_at_grid(3,4,0,1,-4,-3,GRAY);init_static_hash();rebuild_voxel_hash();}
  ai_select(15,0,players[0].pos,AI_ALERT_SEEN);aiCombat[15].reaction=0;aiCombat[15].wasVisible=true;
  bool action=false;for(int n=0;n<100;n++){
   UpdateBot(15,1.f/60);action|=players[15].matter<100||grenadeInventory[15]<2&&(type==PICKUP_GRENADE||type==PICKUP_ADHESIVE)||ex_grenade_count(15,GRENADE_GRAVITY)<2&&type==PICKUP_BLACK_HOLE||ex_grenade_count(15,GRENADE_IGLOO)<2&&type==PICKUP_IGLOO;
   if(action)break;eq_tick(1.f/60);
  }
  fprintf(stderr,"carrier %s action=%d matter=%.1f\n",pickup_type_label(type),action,players[15].matter);assert(action);
 }
 puts("All 15 carrier types autonomously use their equipment passed");
}
static void checkpoints(void){
 fixture();firefightFoundry=true;multiplayerActivePlayers=2;firefightCheckpoints=true;ResetGame();gameState=GAME_STATE_PLAYING;
 assert(fpCheckpoint.valid&&fpCheckpoint.wave==1);for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)assert(pickups[n].type==PICKUP_HEALTH||pickups[n].type==PICKUP_AMMO||pickups[n].type==PICKUP_GRENADE);assert(!grenadeInventory[0]);
 size_t glob=carrier(PICKUP_GLOB,4);fp_death(4,players[0].pos);assert(fp_collect(0,glob,false));size_t ground=carrier(PICKUP_REINFORCED,5);fp_death(5,players[1].pos);
 grenadeInventory[0]=1;firefightScore=123;players[0].kills=7;teamUpgradeLevels[1]=2;start_firefight_wave(3);size_t count=fpCount;uint32_t next=fpNextId,unlocks=fpUnlocked;int voxelsBefore=voxel_count;
 for(int attempt=0;attempt<3;attempt++){
  size_t rocket=carrier(PICKUP_ROCKET,6);fp_death(6,players[0].pos);assert(fp_collect(0,rocket,true));assert(fp_collect(0,ground,false));firefightTeamLives=0;firefightScore=999;grenadeInventory[0]=0;teamUpgradeLevels[1]=8;
  remove_voxel_index(voxel_count-1);for(int i=0;i<2;i++){players[i].invuln_timer=0;players[i].isExposed=true;players[i].matter=0;kill_player(i,8,false,false);}ff_update_coop(.01f);
  assert(firefightWaveState==FIREFIGHT_STATE_RESTART&&fpRestartTimer==3&&gameState==GAME_STATE_PLAYING);fp_tick(3.1f);
  assert(firefightWave==3&&firefightWaveState==FIREFIGHT_STATE_ACTIVE&&firefightScore==123&&players[0].kills==7&&teamUpgradeLevels[1]==2);
  assert(!memcmp(ffRemaining,fpCheckpoint.remaining,sizeof(ffRemaining))&&!memcmp(fpPending,fpCheckpoint.pending,sizeof(fpPending)));
  assert(fpCount==count&&fpNextId==next&&fpUnlocked==unlocks&&fpItems[glob].owner==0&&fpItems[ground].owner==-1&&eq_gun(0)==EQ_GLOB&&!eq_upgrade(0,EQ_BUILD));
  assert(grenadeInventory[0]==1&&players[0].matter==100&&firefightTeamLives==3&&voxel_count==voxelsBefore);
 }
 firefightCheckpoints=false;firefightTeamLives=0;for(int i=0;i<2;i++){players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,8,false,false);}ff_update_coop(.1f);assert(gameState==GAME_STATE_GAMEOVER&&firefightWaveState==FIREFIGHT_STATE_DEFEAT);
 puts("Checkpoint retries: fresh terrain, equipment/loot rollback, scores, lives, upgrades and option-off defeat passed");
}
static void lobby(void){
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_MENU;firefightPractice=false;firefightCheckpoints=false;couch_begin_lobby(false);keys[KEY_C]=true;couch_frame();memset(keys,0,sizeof(keys));assert(firefightCheckpoints&&!couchReady[0]);
 pad=true;buttons[GAMEPAD_BUTTON_RIGHT_FACE_LEFT]=true;couch_frame();memset(buttons,0,sizeof(buttons));assert(!firefightCheckpoints);pad=false;
 RenderTexture2D rt=LoadRenderTexture(960,540);BeginTextureMode(rt);couch_title("FIREFIGHT LOBBY");couch_row(0,"Twelve-wave Firefight",false);couch_row(1,"Arena: Foundry Holdout",false);couch_row(2,"Checkpoints: OFF | X / C: toggle",false);couch_row(3,"P1: Keyboard | Not ready",false);EndTextureMode();Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,".build/playtest-run/firefight-checkpoints.png");UnloadImage(im);UnloadRenderTexture(rt);
 puts("Keyboard/controller checkpoint lobby toggle passed");
}
int main(int argc,char **argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Firefight progression tests");init_pbd_thread_pool();assert(initialize_physics_backend());schedule();opening_rounds();second_wave_rounds();ownership();enemy_weapons();autonomous_carriers();checkpoints();lobby();fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Firefight progression passed");return 0;}


