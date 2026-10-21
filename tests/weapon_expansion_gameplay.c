#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();worldPeriodic=false;
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(goliathStates,0,sizeof(goliathStates));
 gameMode=GAME_MODE_DEATHMATCH;gameState=GAME_STATE_PLAYING;activePlayers=2;enemyLab=false;creativeModeActive=false;netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;players[i].matter=players[i].matterMax=100;players[i].kd_ratio=1;players[i].pos=(Vector3){.25f,2.25f,-3*i+.25f};players[i].last_shot_time=players[i].last_build_time=-100;}
 eq_reset_all();clarity_reset();combat_ui_reset();
}
static void gun(int type){assert(eq_pickup(0,type,true));equipment[0].switchTime=0;}
static void guns(void){
 fixture();gun(PICKUP_DMR);for(int n=0;n<4;n++){equipment[0].cooldown[EQ_DMR]=0;assert(eq_fire(0));if(n<3)assert(players[1].matter==75-25*n&&players[1].respawn_timer==0);}fprintf(stderr,"dmr dead=%f matter=%f\n",players[1].respawn_timer,players[0].matter);assert(players[1].respawn_timer>0);
 fixture();gun(PICKUP_GLOCK);players[1].matter=15;assert(eq_fire(0));assert(players[1].isExposed&&players[1].respawn_timer==0);equipment[0].cooldown[EQ_GLOCK]=0;assert(eq_fire(0)&&players[1].respawn_timer>0);
 fixture();gun(PICKUP_GLOB);assert(eq_fire(0));ex_tick(.1f);assert(players[1].matter==82&&exSuppression[1]==5);add_player_matter(&players[1],100);assert(players[1].matter==82);ex_tick(4);equipment[0].cooldown[EQ_GLOB]=0;assert(eq_fire(0));ex_tick(.1f);assert(exSuppression[1]==5&&players[1].matter==64);ex_tick(5);add_player_matter(&players[1],10);assert(players[1].matter==74);
 fixture();gun(PICKUP_GLOB);for(int n=0;n<EX_GLOBS;n++)exGlobs[n].id=n+1;assert(!eq_fire(0)&&players[0].matter==100);
 fixture();gun(PICKUP_GLOB);add_static_box_at_grid(0,0,4,4,-3,-3,GRAY);init_static_hash();rebuild_voxel_hash();assert(eq_fire(0));ex_tick(.1f);assert(players[1].matter==100&&!exSuppression[1]);int splats=0;for(int n=0;n<EX_SPLATS;n++)splats+=exSplats[n].id!=0;assert(splats==1);
 fixture();gun(PICKUP_DMR);add_static_box_at_grid(0,0,4,4,-3,-3,GRAY);init_static_hash();rebuild_voxel_hash();assert(eq_fire(0));assert(players[1].matter==100&&!occupied(0,4,-3));
 fixture();gun(PICKUP_DMR);players[1].invuln_timer=1;assert(eq_fire(0)&&players[1].matter==100);
 fixture();gun(PICKUP_GLOB);gameMode=GAME_MODE_FIREFIGHT;combatantTeam[1]=combatantTeam[0];assert(eq_fire(0));ex_tick(.1f);assert(players[1].matter==100&&!exSuppression[1]);
 fixture();gun(PICKUP_GLOCK);players[0].matter=1;assert(eq_fire(0)&&players[0].isExposed&&eq_gun(0)==EQ_GLOCK);
 puts("guns/suppression/costs/occlusion passed");
}
static void grenades_test(void){
 fixture();assert(ex_grenade_pickup(0,GRENADE_GRAVITY));assert(ex_grenade_pickup(0,GRENADE_IGLOO));ex_cycle(0);assert(exSelectedGrenade[0]==1);assert(clarity_throw(0));assert(!ex_throw_available(0,1));assert(exSelectedGrenade[0]==2);grenadeCooldown[0]=0;assert(clarity_throw(0));assert(exSelectedGrenade[0]==0);
 GrenadeState g={.kind=GRENADE_GRAVITY,.owner=0,.pos={10,.2f,10},.vel={0,-1,0},.fuse=-1};clarity_grenade_step(&g,.1f);assert(g.fuse>1.3f);for(int n=0;n<15;n++)clarity_grenade_step(&g,.1f);assert(g.fuse==0);ex_grenade_deploy(g);assert(exFields[0].id&&players[0].matter==100);
 fixture();exGrenades[0][1]=1;gameMode=GAME_MODE_FIREFIGHT;firefightHumans=2;firefightPractice=false;kill_player(0,-1,false,false);assert(combatantLife[0]==LIFE_DOWNED&&exGrenades[0][1]==1);exSuppression[0]=4;ff_restore_human(0,35);assert(exGrenades[0][1]==1&&!exSuppression[0]);ex_reset_player(0);assert(!exGrenades[0][1]);
 puts("grenade selection/reservation/fuse/lifecycle passed");
}
static void gravity(void){
 fixture();players[0].pos=(Vector3){2,2,0};players[1].pos=(Vector3){-2,2,0};GravityField f={.id=1,.owner=0,.pos={0,2,0},.age=1};exFields[0]=f;ex_actor_pull(0,.1f);ex_actor_pull(1,.1f);assert(players[0].vel.x<0&&players[1].vel.x>0);assert(players[0].matter==100&&players[1].matter==100);
 players[1].vel=(Vector3){0};add_static_box_at_grid(-3,-3,0,8,-2,2,GRAY);init_static_hash();rebuild_voxel_hash();ex_actor_pull(1,.1f);assert(v_length(players[1].vel)==0);
 int v=addVoxel(2,2,1,false,true,GRAY,0);uint64_t id=voxels[v].identity;assert(freeze_dynamic_cluster_in_place(&v,1));assert(voxels[v].sleeping);ex_gravity_pull(&f,.1f);assert(voxels[v].vel.x<0&&!voxels[v].sleeping&&voxels[v].particles[0]->inv_mass>0&&voxels[v].explosionUntil>clarityClock);float x=voxels[v].pos.x;simulate_voxel_pbd_steps(1.f/60,1);v=enemy_voxel(id);assert(v>=0&&voxels[v].pos.x<x);
 fixture();add_static_box_at_grid(-7,7,1,7,-7,7,GRAY);init_static_hash();rebuild_voxel_hash();f=(GravityField){.id=1,.owner=0,.pos={0,2,0},.age=1};
 for(int frame=0;frame<120;frame++){int budget=32,old=f.activated;ex_gravity_activate(&f,&budget);assert(f.activated-old<=32&&f.activated<=256);}
 assert(f.activated==256);printf("gravity activated %d units, bounded per frame\n",f.activated);
}
static void shelter(void){
 fixture();players[0].pos=(Vector3){.25f,BASE_EYE_HEIGHT,.25f};players[1].pos.x=20;
 GrenadeState g={.owner=0,.kind=GRENADE_IGLOO,.pos={.25f,.15f,.25f},.launchOrigin={.25f,1,-4},.incarnation=exIncarnation[0]};
 assert(ex_igloo(g));int count=voxel_count;assert(count>150&&count<IGLOO_MAX_CELLS);for(int v=0;v<voxel_count;v++)assert(voxels[v].reinforced&&voxels[v].persistentConstruction&&voxels[v].owner==0);
 for(int v=0;v<voxel_count;v++)voxels[v].lifeFrames=RECYCLE_OWNED_STATIC_MAX_FRAMES+20;recycle_dead_voxels(1.f/60);assert(voxel_count==count&&recycleQueueCount==0);
 assert(!occupied(0,6,0)&&!occupied(0,0,-4));assert(!ex_igloo(g)&&voxel_count==count);
 // Activate the entire shelter and check that its reinforced shell stays standing.
 assert(activate_all_static_voxels(0));for(int n=0;n<1800;n++){simulate_voxel_pbd_steps(1.f/60,1);recycle_dead_voxels(1.f/60);deactivate_sleeping_voxels();}
 int high=0;for(int v=0;v<voxel_count;v++){assert(v_isfinite(voxels[v].pos));if(voxels[v].pos.y>2.5f)high++;}fprintf(stderr,"igloo cells=%d high-after-30s=%d\n",voxel_count,high);assert(high>=8);
 fixture();g.incarnation=exIncarnation[0];g.pos.y=8;ex_grenade_deploy(g);assert(exGrenades[0][1]==1&&voxel_count==0&&exFailedBuild[0]>0);
 ex_reset_player(0);ex_grenade_deploy(g);assert(exGrenades[0][1]==0);
 puts("igloo atomic placement/refund/30s activated stability passed");
}
static void maps(void){for(int map=0;map<2;map++){fixture();currentWorldType=map?WORLD_TYPE_CROSSCUT:WORLD_TYPE_TEST;gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=!map;firefightPractice=false;useCustomMap=false;multiplayerActivePlayers=2;firefightHumans=2;ResetGame();int counts[PICKUP_TYPE_COUNT]={0};for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active){counts[pickups[n].type]++;Vector3 at;assert(clarity_supply_position(pickups[n].pos,0,&at));assert(fabsf(at.y-pickups[n].pos.y)<.01f);}for(int t=PICKUP_GLOB;t<PICKUP_TYPE_COUNT;t++){fprintf(stderr,"map%d type%d count%d\n",map,t,counts[t]);assert(counts[t]==(map?4:0));}}}
static void captures(void){
 fixture();weaponRenderingReady=weapon_renderer_init();assert(weaponRenderingReady);assert(character_renderer_init());RenderTexture2D rt=LoadRenderTexture(1280,720);int guns[]={EQ_GLOB,EQ_GLOCK,EQ_DMR};
 for(int k=0;k<3;k++){equipment[0].gun[0]=guns[k];equipment[0].selected=0;weaponVisuals[0].shot_time=-1000;players[0].pos=(Vector3){0,1,0};Camera3D camera={players[0].pos,{0,1,-1},{0,1,0},60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){22,28,38,255});BeginMode3D(camera);DrawGrid(30,1);WeaponPose pose=player_weapon_pose(0,true);weapon_draw(pose,&weaponVisuals[0],camera.position,exClock,false,WEAPON_COLOR);character_draw_paws(pose,player_palette_color(0),false);weapon_draw(pose,&weaponVisuals[0],camera.position,exClock,false,WEAPON_GLASS);EndMode3D();eq_hud(0,1280,720);EndTextureMode();Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);ExportImage(image,TextFormat("artifacts/expansion-gameplay-%s.png",eqGuns[guns[k]].name));UnloadImage(image);
 }
 fixture();players[0].pos=(Vector3){.25f,BASE_EYE_HEIGHT,-3.25f};players[1].pos.x=20;
 assert(ex_igloo((GrenadeState){.owner=0,.kind=GRENADE_IGLOO,.pos={.25f,.15f,.25f},.launchOrigin={0,1,-4}}));exReveals[0].age=.6f;
 rebuild_all_voxel_surfaces();init_static_hash();Camera3D exterior={{-6,4,-7},{0,1,0},{0,1,0},48,CAMERA_PERSPECTIVE};drawViewPlayerIndex=-1;sample_character(0,0);
 BeginTextureMode(rt);ClearBackground((Color){24,32,45,255});BeginMode3D(exterior);DrawPlane((Vector3){0},(Vector2){50,50},DARKGRAY);DrawVoxels(exterior);ex_draw(WEAPON_COLOR);draw_rat_world(0,false);EndMode3D();EndTextureMode();Image shelterImage=LoadImageFromTexture(rt.texture);ImageFlipVertical(&shelterImage);ExportImage(shelterImage,"artifacts/expansion-gameplay-Igloo.png");UnloadImage(shelterImage);
 UnloadRenderTexture(rt);weapon_renderer_shutdown();character_renderer_shutdown();
}
static void edge_cases(void){
 fixture();gun(PICKUP_GLOB);players[1].pos.x=20;players[0].matter=50;exSuppression[0]=5;pickups[0]=(Pickup){.active=true,.type=PICKUP_AMMO,.pos=players[0].pos};update_pickups(0);assert(pickups[0].active&&players[0].matter==50);
 for(int n=0;n<5;n++)ex_splat(players[1].pos,1,-1,false);int count=0;for(int n=0;n<EX_SPLATS;n++)count+=exSplats[n].id!=0;assert(count==3);players[1].pos.x+=2;ex_tick(.1f);for(int n=0;n<EX_SPLATS;n++)if(exSplats[n].id)assert(exSplats[n].pos.x==players[1].pos.x);
 fixture();gun(PICKUP_GLOB);int armor=addVoxel(.25f,2.25f,-1,false,true,GRAY,0);goliathStates[1].count=1;goliathStates[1].voxels[0].voxelIdentity=voxels[armor].identity;rebuild_voxel_hash();assert(eq_fire(0));ex_tick(.1f);assert(!exSuppression[1]&&players[1].matter==100&&voxel_count==1);
 fixture();gun(PICKUP_DMR);worldPeriodic=true;worldFloorHalfSize=20;players[1].pos.z=10;assert(eq_fire(0));assert(players[1].matter==75);worldPeriodic=false;
 fixture();exGrenades[0][0]=2;exSelectedGrenade[0]=1;for(int n=0;n<EX_FIELDS;n++)exFields[n]=(GravityField){.id=n+1,.owner=n};assert(!clarity_throw(0)&&exGrenades[0][0]==2);
 fixture();players[1].pos.x=20;GrenadeState g={.owner=0,.kind=GRENADE_IGLOO,.pos={.25f,.15f,.25f}};int saved=particle_pool_count;particle_pool_count=MAX_PARTICLES;assert(!ex_igloo(g)&&voxel_count==0);particle_pool_count=saved;
 fixture();netTransport.role=NET_ROLE_CLIENT;uint8_t bytes[8];NetWriter w;net_writer_init(&w,bytes,sizeof(bytes));net_write_u8(&w,255);NetReader r;net_reader_init(&r,bytes,w.length);assert(ex_receive(NET_MSG_EXPANSION,&r,10));netTransport.role=NET_ROLE_OFFLINE;
 fixture();playerInput[1]=INPUT_TYPE_BOT_MEDIUM;players[1].pos=(Vector3){2,2,0};exFields[0]=(GravityField){.id=1,.owner=0,.pos={0,2,0},.age=1};
 for(int n=0;n<60;n++){players[1].vel=(Vector3){0};ex_actor_pull(1,1.f/60);}assert(exBotPull[1].x < -3);
 puts("periodic rays, armor, pickups, attachments, allocation failures, AI pull passed");
}
#ifndef glFinish
__declspec(dllimport) void __stdcall glFinish(void);
#endif
static int compare_ms(const void *a,const void *b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static void stress(void){
 weaponRenderingReady=weapon_renderer_init();assert(weaponRenderingReady);assert(character_renderer_init());
 fixture();gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;useCustomMap=false;multiplayerActivePlayers=firefightHumans=4;currentWorldType=WORLD_TYPE_TEST;ResetGame();gameState=GAME_STATE_PLAYING;
 for(int i=4;i<MAX_COMBATANTS;i++)players[i].respawn_timer=9999;
 for(int i=0;i<4;i++){players[i].invuln_timer=999;players[i].pitch=-12;equipment[i].gun[0]=i%2?EQ_DMR:EQ_GLOB;exGrenades[i][0]=exGrenades[i][1]=1;exSelectedGrenade[i]=i%2+1;}
 firefightWaveState=FIREFIGHT_STATE_ACTIVE;firefightWaveBannerTimer=0;RenderTexture2D views[4]={0};int count=0,w=0,h=0;drawWeaponModels=true;
 for(int phase=0;phase<3;phase++){
  if(phase==1)for(int i=0;i<4;i++)exFields[i]=(GravityField){.id=exNextId++,.owner=i,.pos=crosscut_rotate((Vector3){20,2,14},i)};
  if(phase==2)for(int frame=0;frame<1200;frame++){eq_tick(1.f/60);simulate_voxel_pbd_steps(1.f/60,1);recycle_dead_voxels(1.f/60);deactivate_sleeping_voxels();rebuild_static_hash_if_dirty();}
  double times[360],simTime=0,drawTime=0;int frames=phase==1?360:120;
  for(int frame=0;frame<frames;frame++){double start=GetTime();
   if(phase==1 && frame<240)for(int i=0;i<4;i++){players[i].matter=100;players[i].isExposed=false;if(frame%10==0)eq_fire(i);}
   double firedAt=GetTime();clarity_tick(1.f/60);eq_tick(1.f/60);double tickedAt=GetTime();simulate_voxel_pbd_steps(1.f/60,1);double physicsAt=GetTime();recycle_dead_voxels(1.f/60);rebuild_static_hash_if_dirty();simTime+=GetTime()-start;double renderStart=GetTime();render_gameplay_view(views,&count,&w,&h,false);rlDrawRenderBatchActive();glFinish();drawTime+=GetTime()-renderStart;times[frame]=(GetTime()-start)*1000;if(phase==1 && times[frame]>45)fprintf(stderr,"SLOW f%d fire%.2f tick%.2f physics%.2f recycle%.2f render%.2f\n",frame,(firedAt-start)*1000,(tickedAt-firedAt)*1000,(physicsAt-tickedAt)*1000,(renderStart-physicsAt)*1000,(GetTime()-renderStart)*1000);
   int globs=0,fields=0;for(int n=0;n<EX_GLOBS;n++)globs+=exGlobs[n].id!=0;for(int n=0;n<EX_FIELDS;n++)fields+=exFields[n].id!=0;assert(globs<=128&&fields<=4);for(int v=0;v<voxel_count;v++)assert(v_isfinite(voxels[v].pos));
   if(phase==1&&frame==90)for(int i=0;i<4;i++){Image im=LoadImageFromTexture(views[i].texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/expansion-four-%d.png",i));UnloadImage(im);}
  }
  fprintf(stderr,"EXPANSION phase%d mean sim=%.3f draw=%.3f ms\n",phase,simTime*1000/frames,drawTime*1000/frames);
  qsort(times+20,frames-20,sizeof(double),compare_ms);fprintf(stderr,"EXPANSION STRESS phase%d median=%.3f p95=%.3f worst=%.3f ms voxels=%d\n",phase,times[20+(frames-20)/2],times[20+(frames-20)*95/100],times[frames-1],voxel_count);
 }
 for(int n=0;n<EX_FIELDS;n++)assert(!exFields[n].id);for(int n=0;n<EX_GLOBS;n++)assert(!exGlobs[n].id);
 for(int i=0;i<count;i++){UnloadRenderTexture(views[i]);weapon_bloom_unload(&weaponBlooms[i]);}weapon_renderer_shutdown();character_renderer_shutdown();
}
int main(int argc,char **argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Expansion gameplay");init_pbd_thread_pool();assert(initialize_physics_backend());if(getenv("FPS_EXPANSION_ONLY_STRESS"))stress();else {guns();grenades_test();gravity();shelter();maps();edge_cases();captures();if(getenv("FPS_EXPANSION_STRESS"))stress();}clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Expansion gameplay passed");return 0;}



