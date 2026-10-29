#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();worldPeriodic=false;
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;activePlayers=MAX_COMBATANTS;firefightHumans=2;
 firefightPractice=enemyLab=equipmentLab=creativeModeActive=false;firefightCheckpoints=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;
 netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;fp_new_run();eq_reset_all();clarity_reset();ex_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_EASY;combatantLife[i]=i<2?LIFE_ALIVE:LIFE_UNUSED;combatantTeam[i]=i<4?0:1;players[i].respawn_timer=i<2?0:9999;players[i].matter=players[i].matterMax=100;players[i].pos=(Vector3){i*3,1,0};ai_reset(i);}
}
static void recovery(void){
 for(int fps=10;fps<=120;fps+=10){
  fixture();assert(eq_pickup(0,PICKUP_SHIELD_REGEN,false));players[0].matter=0;players[0].isExposed=true;
  for(int n=0;n<fps*4;n++)eq_shield_tick(1.f/fps);assert(players[0].matter==0);
  for(int n=0;n<fps*2;n++)eq_shield_tick(1.f/fps);assert(fabsf(players[0].matter-10)<.01f && !players[0].isExposed);
  for(int n=0;n<fps*10;n++)eq_shield_tick(1.f/fps);assert(fabsf(players[0].matter-50)<.001f);
 }
 fixture();eq_pickup(0,PICKUP_SHIELD_REGEN,false);players[0].matter=10;eq_shield_tick(5.25f);assert(fabsf(players[0].matter-12.5f)<.001f);
 players[0].matter=80;eq_shield_tick(2);assert(players[0].matter==80 && !equipment[0].shieldRegenerating);
 players[0].matterMax=200;eq_shield_tick(3);assert(players[0].matter==100);
 eq_spend(0,20);eq_shield_tick(1);assert(players[0].matter==90); // spending is not damage
 apply_matter_damage(0,4,10);assert(equipment[0].shieldDelay==5 && !equipment[0].shieldRegenerating);
 eq_shield_tick(4.9f);assert(players[0].matter==80);eq_shield_tick(.2f);assert(fabsf(players[0].matter-81)<.001f);
 equipment[0].shieldDelay=0;players[0].invuln_timer=2;apply_matter_damage(0,4,10);assert(equipment[0].shieldDelay==0 && fabsf(players[0].matter-81)<.001f);
 players[0].invuln_timer=0;apply_matter_damage(0,1,10);assert(equipment[0].shieldDelay==0);
 players[4].respawn_timer=0;players[4].matter=100;apply_matter_damage(4,0,10);assert(equipment[0].shieldDelay==0);
 tetherDrainDamage=true;apply_matter_damage(0,4,1);tetherDrainDamage=false;assert(equipment[0].shieldDelay==5);
 equipment[0].shieldDelay=0;exSuppression[0]=5;float matter=players[0].matter;eq_shield_tick(1);assert(players[0].matter==matter && !equipment[0].shieldRegenerating);
 exSuppression[0]=0;eq_shield_tick(1);assert(players[0].matter>matter);
 gameState=GAME_STATE_MENU;matter=players[0].matter;eq_shield_tick(10);assert(players[0].matter==matter);
 gameState=GAME_STATE_PLAYING;combatantLife[0]=LIFE_DOWNED;eq_shield_tick(10);assert(players[0].matter==matter);
 ff_restore_human(0,35);assert(equipment[0].shieldDelay==5 && eq_upgrade(0,EQ_SHIELD));eq_shield_tick(4);assert(players[0].matter==35);
 netTransport.role=NET_ROLE_CLIENT;eq_shield_tick(10);assert(players[0].matter==35 && equipment[0].shieldDelay==1);netTransport.role=NET_ROLE_OFFLINE;
 assert(!eq_pickup(4,PICKUP_SHIELD_REGEN,false));
 puts("Recovery: delay, rate, cap, exposure, damage filters, spending, suppression, revive and client authority passed");
}
static void drops(void){
 fixture();
 for(int i=4;i<7;i++){
  init_firefight_bot_entity(i,ENEMY_TYPE_SHIELDBEARER,INPUT_TYPE_BOT_MEDIUM);fp_spawn(i,ENEMY_TYPE_SHIELDBEARER);
  assert(!fp_carrier(i) && !(equipment[i].upgrades&EQ_SHIELD));
  players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,0,false,false);
 }
 assert(fpCount==3);for(int n=0;n<3;n++)assert(fpItems[n].type==PICKUP_SHIELD_REGEN && fpItems[n].owner==-1);
 uint32_t id=fpItems[0].id;assert(fp_collect(0,0,false));assert(!fp_collect(0,1,false));assert(fpItems[1].owner==-1);
 assert(fp_collect(1,1,false));assert(fpBannerWrite==1);assert(eq_gun(0)==EQ_DEFAULT && equipment[0].gun[1]==EQ_EMPTY);
 players[0].invuln_timer=0;players[0].matter=0;players[0].isExposed=true;kill_player(0,4,false,false);
 assert(combatantLife[0]==LIFE_DOWNED && fpItems[0].owner==0 && eq_upgrade(0,EQ_SHIELD));
 ff_restore_human(0,35);assert(eq_upgrade(0,EQ_SHIELD));
 players[0].death_pos=players[0].pos;ff_resolve_death(0);assert(!eq_upgrade(0,EQ_SHIELD) && fpItems[0].owner==-1 && fpItems[0].id==id);
 ff_restore_human(0,100);assert(fp_collect(0,0,false));
 firefightCheckpoints=true;firefightWave=4;fp_capture();assert(fpCheckpoint.equipment[0].upgrades&EQ_SHIELD);
 eq_reset(0);fp_drop(&fpItems[0],players[0].pos);fp_restore();assert(eq_upgrade(0,EQ_SHIELD) && fpItems[0].owner==0 && fpItems[0].id==id && equipment[0].shieldDelay==5);
 fp_new_run();eq_reset_all();assert(fpCount==0 && !eq_upgrade(0,EQ_SHIELD));
 puts("Drops: every bearer, native AI, duplicates, one celebration, downs/revives/death, checkpoint and reset passed");
}
static void packet(float delay,int regen,int flags,uint32_t tick){
 uint8_t data[2048];NetWriter w;net_writer_init(&w,data,sizeof(data));net_write_u8(&w,0);
 for(int i=0;i<4;i++){
  net_write_u8(&w,EQ_DEFAULT);net_write_u8(&w,EQ_EMPTY);net_write_u8(&w,0);net_write_u8(&w,i?0:flags);net_write_u8(&w,0);
  for(int n=0;n<3;n++)net_write_f32(&w,0);net_write_f32(&w,i?0:delay);net_write_u8(&w,i?0:regen);
  for(int n=0;n<EQ_GUN_COUNT;n++)net_write_f32(&w,0);eq_write_vec(&w,(Vector3){0});net_write_u8(&w,0);net_write_u32(&w,0);
 }
 NetReader r;net_reader_init(&r,data,w.length);assert(eq_receive(NET_MSG_EQUIPMENT,&r,tick));
}
static void networking(void){
 fixture();netTransport.role=NET_ROLE_CLIENT;packet(3,0,EQ_SHIELD,10);assert(eq_upgrade(0,EQ_SHIELD) && equipment[0].shieldDelay==3);
 packet(0,1,EQ_SHIELD,11);assert(equipment[0].shieldRegenerating);
 packet(NAN,1,EQ_SHIELD,12);assert(eqPartTick[0]==11);
 packet(0,2,EQ_SHIELD,13);assert(eqPartTick[0]==11);
 packet(0,1,0,14);assert(eqPartTick[0]==11);
 packet(0,0,128,15);assert(eqPartTick[0]==11);
 packet(0,0,0,16);assert(!eq_upgrade(0,EQ_SHIELD));netTransport.role=NET_ROLE_OFFLINE;
 puts("Network: upgrade/delay/recovery round trip, reset and malformed-state rejection passed");
}
static void captures(void){
 fixture();weapon_renderer_init();character_renderer_init();eq_pickup(0,PICKUP_SHIELD_REGEN,false);players[0].matter=20;equipment[0].shieldDelay=0;eq_shield_tick(.1f);
 drawViewPlayerIndex=0;weaponFrameSampled=true;weaponFrameTime=1;
 for(int small=0;small<2;small++)for(int aim=0;aim<2;aim++){
  int w=small?480:960,h=small?270:540;RenderTexture2D rt=LoadRenderTexture(w,h);
  Camera3D camera={{3,2.3f,6},{0,1,0},{0,1,0},aim?30:60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){25,32,43,255});BeginMode3D(camera);
  DrawPlane((Vector3){0,0,0},(Vector2){20,20},(Color){60,67,75,255});
  CharacterPose pose=character_pose((Vector3){-.7f,BASE_EYE_HEIGHT,0},145,0,0,true,false,BLUE,1);pose.effect=SKYBLUE;pose.effectAmount=.2f;character_draw(pose,false);
  eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_COLOR);eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_GLASS);
  EndMode3D();draw_hud_bars(0,&players[0],w,h);eq_hud(0,w,h);DrawCircle(w/2,h/2,2,WHITE);EndTextureMode();
  Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat(".build/shield-regen-%d-%d.png",small,aim));UnloadImage(im);UnloadRenderTexture(rt);
 }
 RenderTexture2D split=LoadRenderTexture(960,540),tile=LoadRenderTexture(480,270);
 for(int i=0;i<4;i++){
  combatantLife[i]=LIFE_ALIVE;players[i].respawn_timer=0;players[i].matter=20+i*8;eq_pickup(i,PICKUP_SHIELD_REGEN,false);
  equipment[i].shieldDelay=i==1?3:0;equipment[i].shieldRegenerating=i!=1;
  Camera3D c={{0,1.8f,4},{0,1,0},{0,1,0},60,CAMERA_PERSPECTIVE};
  drawViewPlayerIndex=i;
  BeginTextureMode(tile);ClearBackground((Color){25,32,43,255});BeginMode3D(c);
  DrawPlane((Vector3){0,0,0},(Vector2){20,20},(Color){60,67,75,255});
  CharacterPose pose=character_pose((Vector3){-.7f,BASE_EYE_HEIGHT,0},135+i*20,0,0,true,false,player_palette_color(i),1);
  if(i!=1){pose.effect=SKYBLUE;pose.effectAmount=.2f;}character_draw(pose,false);
  eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_COLOR);eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_GLASS);
  if(i==3)DrawCube((Vector3){.8f,1.2f,.7f},1,1.2f,.3f,GRAY);
  EndMode3D();draw_hud_bars(i,&players[i],480,270);eq_hud(i,480,270);EndTextureMode();
  BeginTextureMode(split);if(i==0)ClearBackground(BLACK);DrawTextureRec(tile.texture,(Rectangle){0,0,480,-270},(Vector2){(i%2)*480,(i/2)*270},WHITE);EndTextureMode();
 }
 Image four=LoadImageFromTexture(split.texture);ImageFlipVertical(&four);ExportImage(four,".build/shield-regen-four-player.png");UnloadImage(four);UnloadRenderTexture(tile);UnloadRenderTexture(split);
 character_renderer_shutdown();weapon_renderer_shutdown();weaponFrameSampled=false;
}
int main(int argc,char**argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Shield regeneration tests");
 init_pbd_thread_pool();assert(initialize_physics_backend());recovery();drops();networking();captures();
 enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Shield regeneration passed");return 0;
}