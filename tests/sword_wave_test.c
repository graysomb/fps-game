#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void setup(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(goliathStates,0,sizeof(goliathStates));
 gameMode=GAME_MODE_DEATHMATCH;gameState=GAME_STATE_PLAYING;activePlayers=2;worldPeriodic=false;worldFloorHalfSize=40;
 enemyLab=creativeModeActive=false;netTransport.role=NET_ROLE_OFFLINE;
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;players[i].matter=players[i].matterMax=100;players[i].kd_ratio=1;players[i].pos=(Vector3){i*30.f,3,0};players[i].last_melee_time=-100;}
 eq_reset_all();clarity_reset();combat_ui_reset();eq_pickup(0,PICKUP_SWORD,false);equipment[0].switchTime=0;
}
static void advance(float seconds,int fps){int steps=(int)roundf(seconds*fps);for(int n=0;n<steps;n++)eq_tick(1.f/fps);}
static void behavior(void){
 for(int fps=30;fps<=120;fps*=2){
  setup();players[1].pos=(Vector3){0,3,-5};perform_sword(0);assert(swordDirection[0]==0);
  uint32_t seq=swordWaves[0].sequence;perform_sword(0);assert(swordWaves[0].sequence==seq && swordNextDirection[0]==1);
  advance(.2f,fps);assert(players[1].respawn_timer==0);advance(.4f,fps);assert(players[1].respawn_timer>0);
  assert(swordWaves[0].hit[1] && swordWaves[0].cuts==0);
  players[0].last_melee_time=-100;perform_sword(0);assert(swordDirection[0]==1);advance(.6f,fps);
  players[0].last_melee_time=-100;perform_sword(0);assert(swordDirection[0]==0);
 }
 setup();players[1].pos=(Vector3){0,3,-7.5f};perform_sword(0);advance(.6f,60);assert(players[1].respawn_timer==0);
 setup();players[1].pos=(Vector3){0,3,-4};players[1].invuln_timer=2;perform_sword(0);advance(.6f,60);assert(players[1].matter==100);
 setup();gameMode=GAME_MODE_FIREFIGHT;firefightHumans=2;combatantTeam[1]=0;players[1].pos=(Vector3){0,3,-4};perform_sword(0);advance(.6f,60);assert(players[1].matter==100);
 setup();perform_sword(0);eq_tick(.1f);players[0].respawn_timer=9;eq_tick(.1f);assert(!swordWaves[0].active);
 setup();players[1].pos=(Vector3){0,3,-5};perform_sword(0);eq_tick(.2f);Vector3 origin=swordWaves[0].origin;
 players[0].pos.x=10;players[0].yaw=180;kill_player(0,-1,false,false);advance(.4f,60);
 assert(players[1].respawn_timer>0 && v_length(v_sub(origin,swordWaves[0].origin))==0);
 setup();worldPeriodic=true;worldFloorHalfSize=16;players[0].pos.z=-14;players[1].pos=(Vector3){0,3,14};perform_sword(0);advance(.6f,60);assert(players[1].respawn_timer>0);
 for(int dir=0;dir<2;dir++){
  setup();swordNextDirection[0]=dir;perform_sword(0);
  for(int n=0;n<4;n++){float age=.15f+n*.05f;WeaponMeleePose m={.progress=age/.6f,.sword_direction=dir};
   WeaponPose p=weapon_pose(WEAPON_SWORD,(Vector3){0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},true,1.77f,&m,0);
   SwordSwingSample shared=weapon_sword_sample(age,dir);Vector3 d=shared.direction;assert(v_dot(v_norm(v_sub(p.muzzle,p.grip_heel)),d)>.98f);
   Vector3 edge=v_norm(v_sub(Vector3Transform((Vector3){-1,0,0},p.transform),Vector3Transform((Vector3){0},p.transform)));
   assert(v_dot(edge,shared.edge_direction)>.999f);}
 }
}
static void swing_cost(void){
 int balances[]={0,24,25,100};
 for(int n=0;n<4;n++){
  setup();players[0].matter=balances[n];players[0].isExposed=balances[n]==0;
  perform_sword(0);bool affordable=balances[n]>=25;
  assert(players[0].meleeSwingActive==affordable && swordWaves[0].active==affordable);
  assert(players[0].matter==balances[n]-(affordable?25:0));
  assert(swordNextDirection[0]==(affordable?1:0));float matter=players[0].matter;
  perform_sword(0);assert(players[0].matter==matter);assert(!eq_fire(0));assert(players[0].matter==matter);
  if(affordable){assert(players[0].isExposed==(balances[n]==25));advance(.6f,60);assert(!swordWaves[0].active);}
 }
 setup();assert(eq_fire(0));assert(players[0].matter==75);advance(.6f,60);perform_sword(0);assert(players[0].matter==75); // Real-time cooldown rejects another attack.
 setup();netTransport.role=NET_ROLE_CLIENT;perform_sword(0);assert(players[0].matter==100 && swordWaves[0].active);netTransport.role=NET_ROLE_OFFLINE;
}
static void ordinary_melee_with_sword(void){
 float harvested[2];int remaining[2];
 for(int blade=0;blade<2;blade++){
  setup();equipment[0].gun[equipment[0].selected]=blade?EQ_BLADE:EQ_DEFAULT;players[0].matter=0;players[0].isExposed=true;
  add_static_box_at_grid(-3,3,2,9,-3,-3,GRAY);init_static_hash();rebuild_voxel_hash();
  perform_melee(0);assert(players[0].meleeSwingActive && !players[0].swordSwing && !swordWaves[0].active && players[0].matter==0);
  assert(swordNextDirection[0]==0 && player_weapon_kind(&players[0])==WEAPON_STANDARD);
  players[0].meleeSwingStartTime=(float)GetTime()-MELEE_ANIM_DURATION_SECONDS*WEAPON_MELEE_PEAK;
  update_melee_swings();harvested[blade]=players[0].matter;remaining[blade]=voxel_count;
  assert(players[0].meleeSwingHitApplied && !swordWaves[0].active);
 }
 assert(harvested[0]>0 && harvested[0]==harvested[1] && remaining[0]==remaining[1]);
 setup();perform_melee(0);assert(!eq_fire(0) && players[0].matter==100); // No overlap with a normal punch.
}
static void terrain(void){
 int counts[3]={0};int trial=0;
 for(int fps=30;fps<=120;fps*=2){
  setup();add_static_box_at_grid(-7,7,0,12,-8,-8,GRAY);init_static_hash();rebuild_voxel_hash();
  int before=voxel_count;perform_sword(0);advance(.6f,fps);counts[trial++]=swordWaves[0].cuts;assert(voxel_count==before);int active=0;for(int v=0;v<voxel_count;v++)active+=voxels[v].simulate;assert(active>=swordWaves[0].cuts);assert(counts[trial-1]>0 && swordWaves[0].cuts<=96);
  int first=voxel_count;players[0].last_melee_time=-100;perform_sword(0);advance(.6f,fps);assert(voxel_count==first && swordWaves[0].cuts>0);
 }
 assert(counts[0]==counts[1] && counts[1]==counts[2]);
 setup();players[1].pos=(Vector3){0,3,-5};add_static_box_at_grid(-10,10,0,15,-8,-4,GRAY);init_static_hash();rebuild_voxel_hash();
 perform_sword(0);eq_tick(.15f);swordWaves[0].cuts=96;advance(.45f,60);assert(players[1].matter==100);
 int stopped=0;for(int n=0;n<SWORD_LANES;n++)stopped+=swordWaves[0].stop[n]<6;assert(stopped>0);
 setup();players[1].pos.x=30;int v=add_static_voxel_at_grid(0,5,-5,GRAY,0);init_static_hash();rebuild_voxel_hash();
 SwordWave wave={0};assert(sword_cut_hit(0,&wave,(VoxelHit){.id=v,.x=0,.y=5,.z=-5}));wave.cuts=1;
 v=enemy_voxel(wave.cutIds[0]);assert(v>=0 && voxels[v].simulate && voxels[v].glueEligible && voxels[v].cutDetached);
 float height=voxels[v].pos.y;for(int n=0;n<12;n++)simulate_voxel_pbd_steps(1.f/60,1);
 v=enemy_voxel(wave.cutIds[0]);assert(v>=0 && voxels[v].pos.y<height && v_isfinite(voxels[v].pos));
 acquire_tether_hit(0,v,(Vector3){0,1,0});assert(players[0].tetherHolding);release_tether(0);
 setup();v=add_static_voxel_at_grid(0,5,-5,GRAY,0);init_static_hash();int count=voxel_count,saved=particle_pool_count,freeSaved=free_particle_count;particle_pool_count=MAX_PARTICLES;free_particle_count=0;
 assert(!sword_cut_hit(0,&wave,(VoxelHit){.id=v,.x=0,.y=5,.z=-5}) && voxel_count==count && !voxels[v].simulate);particle_pool_count=saved;free_particle_count=freeSaved;
 fprintf(stderr,"Sword cut counts at 30/60/120 FPS: %d %d %d\n",counts[0],counts[1],counts[2]);
}
static void lifecycle_and_protection(void){
 setup();perform_sword(0);advance(.6f,60);NetInputCommand change={.pressed=NET_INPUT_SWITCH};
 eq_command(0,&change,.01f);assert(eq_gun(0)==EQ_DEFAULT && swordNextDirection[0]==1);
 eq_command(0,&change,.01f);assert(eq_gun(0)==EQ_BLADE && swordNextDirection[0]==1);
 gameMode=GAME_MODE_FIREFIGHT;firefightHumans=2;firefightPractice=false;kill_player(0,-1,false,false);
 assert(combatantLife[0]==LIFE_DOWNED && swordNextDirection[0]==1);ff_restore_human(0,35);assert(swordNextDirection[0]==1);
 eq_reset_all();assert(!swordNextDirection[0] && !swordWaves[0].active);
 setup();players[1].pos=(Vector3){0,3,-2};perform_sword(0);advance(.4f,60);assert(swordWaves[0].hit[1]);
 players[1].respawn_timer=0;players[1].matter=100;players[1].isExposed=false;advance(.2f,60);assert(players[1].matter==100);
 setup();players[1].pos=(Vector3){0,3,-4};players[1].enemyType=ENEMY_TYPE_GOLIATH;
 int armor=addVoxel(0,3,-3,false,true,GRAY,0);goliathStates[1].count=1;goliathStates[1].voxels[0].voxelIdentity=voxels[armor].identity;rebuild_voxel_hash();
 perform_sword(0);eq_tick(.2f);assert(players[1].matter==100);advance(.4f,60);assert(goliathStates[1].count==0 && players[1].respawn_timer>0);
 setup();activePlayers=5;playerInput[4]=INPUT_TYPE_BOT_HARD;enemy_lab_spawn(4,ENEMY_TYPE_FOUNDRY,(Vector3){0,3,-4});
 float health=enemySpecial[4].coreHealth;perform_sword(0);advance(.6f,60);assert(enemySpecial[4].coreHealth==health);
}
static void stress(void){
 setup();activePlayers=4;for(int i=0;i<4;i++){players[i].pos=(Vector3){10.f*i,3,0};players[i].invuln_timer=999;equipment[i].gun[0]=EQ_BLADE;equipment[i].selected=0;
 add_static_box_at_grid(i*20-7,i*20+7,0,12,-8,-8,GRAY);}
 init_static_hash();rebuild_voxel_hash();double first=0,last=0,worst=0;
 for(int frame=0;frame<480;frame++){
  double start=GetTime();if(frame%48==0)for(int i=0;i<4;i++){players[i].last_melee_time=-100;players[i].matter=100;perform_sword(i);}
  eq_tick(1.f/60);int active=0;for(int i=0;i<4;i++){active+=swordWaves[i].active;assert(swordWaves[i].cuts<=96 && v_isfinite(swordWaves[i].origin));}
  assert(active<=4);double ms=1000*(GetTime()-start);worst=fmax(worst,ms);if(frame<60)first+=ms;if(frame>=420)last+=ms;
 }
 eq_reset_all();for(int i=0;i<4;i++)assert(!swordWaves[i].active);
 fprintf(stderr,"Four-wave stress first %.3f last %.3f worst %.3f ms\n",first/60,last/60,worst);
}
static void preview(void){
 setup();activePlayers=1;weaponRenderingReady=weapon_renderer_init();character_renderer_init();weaponFrameSampled=false;
 add_static_box_at_grid(-7,7,0,12,-8,-8,GRAY);init_static_hash();rebuild_voxel_hash();
 RenderTexture2D rt=LoadRenderTexture(960,540);weaponViewAspect=960.f/540;
 for(int frame=0;frame<48;frame++){
  if(frame==0 || frame==24){players[0].last_melee_time=-100;perform_sword(0);}
  eq_tick(1.f/30);drawViewPlayerIndex=0;
  Camera3D cam={players[0].pos,v_add(players[0].pos,(Vector3){0,0,-1}),{0,1,0},60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){21,29,43,255});BeginMode3D(cam);DrawPlane((Vector3){0},(Vector2){60,60},DARKGRAY);
  for(int v=0;v<voxel_count;v++){DrawCube(voxels[v].pos,.49f,.49f,.49f,GRAY);DrawCubeWires(voxels[v].pos,.5f,.5f,.5f,DARKGRAY);}
  eq_world();EndMode3D();weapon_clear_depth();BeginMode3D(cam);
  WeaponPose pose=player_weapon_pose(0,true);weapon_draw(pose,&weaponVisuals[0],cam.position,GetTime(),false,WEAPON_COLOR);character_draw_paws(pose,BLUE,false);EndMode3D();
  DrawText(frame<24?"RIGHT TO LEFT":"LEFT TO RIGHT",18,18,20,SKYBLUE);EndTextureMode();
  Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/sword-wave-%03d.png",frame));UnloadImage(im);
 }
 assert(eqSlashShader.id && eqSlashShader.id!=rlGetShaderIdDefault());
 for(int view=0;view<3;view++){
  swordWaves[0]=(SwordWave){0};players[0].matter=100;eq_sword(0);players[0].meleeSwingActive=true;eq_tick(.27f);
  bool third=view==1;drawViewPlayerIndex=third?-1:0;weaponViewAspect=view==0?480.f/270:960.f/540;
  RenderTexture2D shot=LoadRenderTexture(view==0?480:960,view==0?270:540);
  Camera3D camera={third?(Vector3){5,5,4}:players[0].pos,third?(Vector3){0,3,-2}:v_add(players[0].pos,(Vector3){0,0,-1}),{0,1,0},60,CAMERA_PERSPECTIVE};
  if(view==2){eq_beam_shutdown();eqSlashTried=true;}
  BeginTextureMode(shot);ClearBackground((Color){21,29,43,255});BeginMode3D(camera);
  DrawPlane((Vector3){0},(Vector2){60,60},DARKGRAY);eq_world();
  if(third){drawWeaponModels=true;draw_rat_world(0,false);weapon_draw(player_weapon_pose(0,false),&weaponVisuals[0],camera.position,GetTime(),false,WEAPON_COLOR);}
  else{weapon_clear_depth();WeaponPose p=player_weapon_pose(0,true);weapon_draw(p,&weaponVisuals[0],camera.position,GetTime(),false,WEAPON_COLOR);character_draw_paws(p,BLUE,false);}
  EndMode3D();EndTextureMode();Image im=LoadImageFromTexture(shot.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/sword-wave-view-%d.png",view));UnloadImage(im);UnloadRenderTexture(shot);
 }

 UnloadRenderTexture(rt);eq_beam_shutdown();
 setup();activePlayers=multiplayerActivePlayers=4;drawWeaponModels=true;weaponBloomEnabled=true;
 for(int i=0;i<4;i++){players[i].pos=(Vector3){i*10.f,3,0};equipment[i].gun[0]=EQ_BLADE;equipment[i].gun[1]=EQ_EMPTY;equipment[i].selected=0;equipment[i].switchTime=0;perform_sword(i);}
 eq_tick(.27f);RenderTexture2D views[4]={0};int count=0,width=0,height=0;render_gameplay_view(views,&count,&width,&height,false);assert(count==4);
 Image four=GenImageColor(960,540,BLACK);
 for(int i=0;i<count;i++){Image tile=LoadImageFromTexture(views[i].texture);ImageFlipVertical(&tile);ImageDraw(&four,tile,(Rectangle){0,0,tile.width,tile.height},(Rectangle){(i%2)*480,(i/2)*270,480,270},WHITE);UnloadImage(tile);}
 ExportImage(four,"artifacts/sword-wave-four.png");UnloadImage(four);
 for(int i=0;i<count;i++){UnloadRenderTexture(views[i]);weapon_bloom_unload(&weaponBlooms[i]);}
 eq_beam_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Sword wave regression");
 init_pbd_thread_pool();assert(initialize_physics_backend());behavior();swing_cost();ordinary_melee_with_sword();terrain();lifecycle_and_protection();stress();preview();
 clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Sword wave tests passed");return 0;
}
