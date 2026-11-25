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
static void combat_feel(void){
 fixture();
 Player *p=&players[0];
 for(int role=0;role<3;role++){
  netTransport.role=role;
  p->onGround=true;p->vel=(Vector3){0};
  NetInputCommand jump={.pressed=NET_INPUT_JUMP,.move_y=32767};
  net_apply_movement_command(p,&jump,.016f);
  assert(!p->onGround && p->vel.z<0 && p->vel.y==JUMP_SPEED);
  Vector3 takeoff=p->vel;
  NetInputCommand reverse={.move_y=-32767,.move_x=32767,.mouse_look=1,.aim_yaw=90};
  net_apply_movement_command(p,&reverse,.1f);
  net_apply_movement_command(p,&(NetInputCommand){0},.1f);
  assert(Vector3Equals(takeoff,p->vel) && p->yaw==90);
  p->onGround=false;p->vel=(Vector3){0};net_apply_movement_command(p,&reverse,.02f);
  assert(Vector3Length(p->vel)==0);
  p->onGround=true;net_apply_movement_command(p,&reverse,.02f);assert(Vector3Length(p->vel)>0);
  p->yaw=0;
 }
 fixture();p->onGround=true;p->last_melee_time=-10;
 players[1].pos=v_add(p->pos,(Vector3){0,0,-2.61f});assert(melee_lunge_target(0)<0);
 players[1].pos.z+=.02f;
 players[2].pos=v_add(p->pos,(Vector3){.2f,0,-2});
 assert(melee_lunge_target(0)==1);
 gameMode=GAME_MODE_FIREFIGHT;combatantTeam[1]=combatantTeam[0];
 assert(melee_lunge_target(0)==2);gameMode=GAME_MODE_DEATHMATCH;
 perform_melee(0);float z=p->pos.z;
 update_melee_lunge(0,.05f);assert(fabsf(p->pos.z-z+.3f)<.002f);
 update_melee_lunge(0,.10f);assert(fabsf(p->pos.z-z+.6f)<.002f);
 update_melee_lunge(0,.12f);assert(fabsf(p->pos.z-z+.6f)<.002f);
 p->meleeSwingStartTime=(float)GetTime()-MELEE_ANIM_DURATION_SECONDS*.52f;
 update_melee_swings();float afterHit=players[1].matter;
 assert(afterHit<100 && p->meleeSwingHitApplied);
 update_melee_swings();assert(players[1].matter==afterHit);
 add_static_box_at_grid(0,0,0,3,-3,-3,spawnGlassColor);init_static_hash();
 assert(melee_lunge_target(0)<0);
 fixture();p->pos.y+=.5f;p->onGround=true;p->last_melee_time=-10;
 add_static_voxel_at_grid(0,0,0,GRAY,0);init_static_hash();
 players[1].pos=v_add(p->pos,(Vector3){0,0,-2});
 perform_melee(0);update_melee_lunge(0,.1f);
 assert(p->pos.z> -CHARACTER_HALF_WIDTH && p->pos.z<.25f);
 fixture();p->onGround=true;p->last_melee_time=-10;
 players[1].pos=v_add(p->pos,(Vector3){0,0,-.85f});perform_melee(0);
 z=p->pos.z;update_melee_lunge(0,.1f);assert(p->pos.z>z-.09f);
 fixture();p->kills=1;p->debrisKills=1;assert(player_points(p)==2);
 for(int sticky=0;sticky<2;sticky++){
  fixture();grenadeInventory[0]=1;
  if(sticky){equipment[0].upgrades|=EQ_ADHESIVE;exSelectedGrenade[0]=GRENADE_STICKY;}
  assert(clarity_throw(0));GrenadeState *g=&grenades[CLARITY_GRENADES-1];
  assert(g->id && g->adhesive==sticky && fabsf(g->vel.z+20)<.001f);
  g->pos.y=100;
  for(int n=0;n<60;n++)clarity_grenade_step(g,1.f/120);
  assert(g->fuse==-1);
  g->pos.y=.15f;g->vel.y=-5;
  for(int n=0;n<41;n++)clarity_grenade_step(g,1.f/120);
  assert(g->fuse>0 && g->attached==sticky);
  clarity_grenade_step(g,.009f);assert(g->fuse==0);
 }
 puts("Movement momentum, melee targeting and stops, smush score, grenade launch and contact fuse passed");
}
static void slime_capture(void){
 fixture();players[0].matter=0;players[0].isExposed=true;
 const float remaining[]={5,2.5f,.05f,0};
 RenderTexture2D rt=LoadRenderTexture(360,72);
 Image strip=GenImageColor(360,4*50,DARKGRAY);
 for(int n=0;n<4;n++){
  exSuppression[0]=remaining[n];
  BeginTextureMode(rt);ClearBackground(DARKGRAY);draw_hud_bars(0,&players[0],360,72);EndTextureMode();
  Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);
  int green=0;
  for(int y=26;y<72;y++)for(int x=0;x<360;x++){
   Color c=GetImageColor(image,x,y);green+=c.g>130 && c.g>c.r+35 && c.g>c.b+70;
  }
  assert((green>0)==(remaining[n]>0));
  players[1].matter=0;players[1].isExposed=true;
  BeginTextureMode(rt);ClearBackground(DARKGRAY);draw_hud_bars(1,&players[1],360,72);EndTextureMode();
  Image other=LoadImageFromTexture(rt.texture);
  for(int y=0;y<46;y++)for(int x=0;x<360;x++){
   Color c=GetImageColor(other,x,y);assert(!(c.g>130 && c.g>c.r+35 && c.g>c.b+70));
  }
  UnloadImage(other);
  ImageDraw(&strip,image,(Rectangle){0,22,360,50},(Rectangle){0,n*50,360,50},WHITE);
  UnloadImage(image);
 }
 ExportImage(strip,".build/playtest-glob-bar.png");UnloadImage(strip);UnloadRenderTexture(rt);
}
static void combat(void){
 for(int target=1;target<=4;target+=3){fixture();players[target].pos=(Vector3){.25f,BASE_EYE_HEIGHT,-2.75f};players[target].matter=10;
  assert(eq_pickup(0,PICKUP_GLOB,false));equipment[0].switchTime=0;
  for(int n=0;n<4;n++){equipment[0].cooldown[EQ_GLOB]=0;assert(eq_fire(0));ex_tick(.1f);assert(players[target].respawn_timer==0&&combatantLife[target]==LIFE_ALIVE);}
  assert(players[target].matter==0&&players[target].isExposed&&exSuppression[target]>4.8f);
 }
 fixture();players[4].enemyType=ENEMY_TYPE_FOUNDRY;enemySpecial[4].coreHealth=players[4].matter=5;enemySpecial[4].coreMax=100;enemySpecial[4].coreOpen=true;
 eq_damage(0,4,18,false,DEATH_GLOB,players[0].pos);assert(players[4].respawn_timer==0&&players[4].matter>0&&players[4].matter<5);
 fixture();players[1].pos=(Vector3){.25f,BASE_EYE_HEIGHT,-2.75f};assert(eq_pickup(0,PICKUP_GLOCK,false));equipment[0].switchTime=0;assert(eq_fire(0));assert(players[1].matter==92);
 equipment[0].cooldown[EQ_GLOCK]=.25f;NetInputCommand swap={.pressed=NET_INPUT_SWITCH};eq_command(0,&swap,.01f);eq_tick(.12f);assert(equipment[0].switchTime==0);eq_command(0,&swap,.01f);assert(equipment[0].cooldown[EQ_GLOCK]>.12f);
 fixture();gameMode=GAME_MODE_FIREFIGHT;combatantTeam[4]=combatantTeam[5]=1;players[4].enemyType=ENEMY_TYPE_ANCHOR;players[4].pos=(Vector3){3,2,0};players[5].pos=(Vector3){5,2,0};players[5].matter=20;enemySpecial[4].equipmentCount=1;
 feedback_tick(1);assert(fabsf(players[5].matter-36)<.01f);exSuppression[5]=5;feedback_tick(1);assert(players[5].matter==36);
 puts("Nonlethal globs, Glock tuning, switch cooldowns, Anchor healing passed");
}
static void construction(void){
 for(int upgraded=0;upgraded<2;upgraded++)for(int target=0;target<=4;target++){
  fixture();players[0].pitch=target==0?-89:-40;equipment[0].upgrades=upgraded?EQ_BUILD:0;
  Vector3 points[EQ_PANEL_CUBES];assert(eq_panel_positions(0,points));if(target)players[target].pos=points[eq_build_count(0)/2];
  float old=players[target].pos.y;assert(eq_build(0));assert(players[target].pos.y>old);assert(players[0].matter==90);
  assert(!player_aabb_hits_world_ex(players[target].pos,player_collision_half(&players[target]),target));
 }
 puts("Both builds lift self, teammates, and enemies onto the completed block");
}
static void refills(void){
 fixture();grenadeInventory[0]=0;exGrenades[0][0]=exGrenades[0][1]=0;
 pickups[0]=(Pickup){.type=PICKUP_GRENADE,.active=true,.pos=players[0].pos};update_pickups(0);
 assert(!pickups[0].active&&pickups[0].respawnTimer==30);for(int k=0;k<3;k++)assert(ex_grenade_count(0,k)==2);
 pickups[0].active=true;update_pickups(0);assert(pickups[0].active);
 exGrenades[0][0]=0;update_pickups(0);assert(!pickups[0].active&&exGrenades[0][0]==2);
 fixture();gameMode=GAME_MODE_FIREFIGHT;grenadeInventory[0]=0;assert(!ex_refill_grenades(0));
 fp_queue_unlock(PICKUP_GRENADE);assert(ex_refill_grenades(0));assert(grenadeInventory[0]==2&&!exGrenades[0][0]&&!exGrenades[0][1]);
 fp_queue_unlock(PICKUP_BLACK_HOLE);fp_queue_unlock(PICKUP_IGLOO);assert(ex_refill_grenades(0));assert(exGrenades[0][0]==2&&exGrenades[0][1]==2);
 puts("Refills fill all available types, preserve unlock tiers, and respawn in 30s");
}
static void attraction(void){
 fixture();players[4].pos=(Vector3){3,BASE_EYE_HEIGHT,0};players[4].enemyType=ENEMY_TYPE_SHIELDBEARER;
 exFields[0]=(GravityField){.id=1,.owner=0,.pos={0,BASE_EYE_HEIGHT,0},.age=1};
 int shield=addVoxel(2.5f,BASE_EYE_HEIGHT,0,false,true,GRAY,0);enemySpecial[4].equipmentCount=1;enemySpecial[4].equipment[0].id=voxels[shield].identity;
 rebuild_voxel_hash();assert(ex_in_field(4));ex_actor_pull(4,.1f);assert(players[4].vel.x<0);
 remove_voxel_index(shield);rebuild_voxel_hash();exBotPull[4]=exBotApplied[4]=(Vector3){0};
 for(int f=0;f<120;f++){players[4].vel=(Vector3){MOVE_SPEED,0,0};exBotApplied[4]=(Vector3){0};net_integrate_player(&players[4],1.f/60,false);}
 printf("Enemy running against field: x %.3f from 3.0\n",players[4].pos.x);assert(players[4].pos.x<2);
 fixture();players[4].pos=(Vector3){3,BASE_EYE_HEIGHT,0};exFields[0]=(GravityField){.id=1,.owner=0,.pos={0,BASE_EYE_HEIGHT,0},.age=1};add_static_box_at_grid(3,3,0,8,-2,2,GRAY);init_static_hash();rebuild_voxel_hash();ex_actor_pull(4,.1f);assert(v_length(players[4].vel)==0);
 puts("Armor does not block attraction, running AI moves inward, static cover still blocks");
}
static void visuals(void){
 fixture();assert(weapon_renderer_init());drawViewPlayerIndex=0;players[0].pos=(Vector3){0,2,0};
 Camera3D camera={players[0].pos,{0,2,-1},{0,1,0},60,CAMERA_PERSPECTIVE};RenderTexture2D rt=LoadRenderTexture(640,480);
 Image captures[2];for(int n=0;n<2;n++){
  if(n)ex_splat((Vector3){0,2,-.1f},0,-1,false);
  BeginTextureMode(rt);ClearBackground(DARKGRAY);BeginMode3D(camera);ex_draw(WEAPON_COLOR);ex_draw(WEAPON_EMISSION);EndMode3D();EndTextureMode();captures[n]=LoadImageFromTexture(rt.texture);ImageFormat(&captures[n],PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
 }
 assert(memcmp(captures[0].data,captures[1].data,640*480*4)==0);for(int n=0;n<2;n++)UnloadImage(captures[n]);
 fixture();gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;drawViewPlayerIndex=0;players[0].pos=(Vector3){0,1.2f,3};assert(fp_reserve(1));fpCount=1;fpItems[0]=(FfItem){.id=1,.type=PICKUP_GOLD_TETHER,.owner=-1,.pos={0,1,0}};
 camera=(Camera3D){players[0].pos,{0,1,0},{0,1,0},45,CAMERA_PERSPECTIVE};weaponFrameSampled=true;weaponFrameTime=0;
 BeginTextureMode(rt);ClearBackground((Color){20,26,35,255});BeginMode3D(camera);fp_draw(WEAPON_COLOR);EndMode3D();EndTextureMode();Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);ExportImage(image,"artifacts/playtest-tether-pickup.png");UnloadImage(image);
 UnloadRenderTexture(rt);weapon_renderer_shutdown();puts("Local glob attachment leaves color/emission capture unchanged; tether pickup captured");
}
int main(void){SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,480,"Playtest fixes");combat_feel();slime_capture();combat();construction();refills();attraction();visuals();clear_world_voxels();CloseWindow();puts("Playtest fixes passed");return 0;}
