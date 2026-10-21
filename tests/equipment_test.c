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
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;players[i].matter=players[i].matterMax=100;players[i].kd_ratio=1;players[i].pos=(Vector3){0,2,-3*i};players[i].last_shot_time=-100;players[i].last_build_time=-100;}
 eq_reset_all();clarity_reset();combat_ui_reset();
}
static void inventory(void){fixture();assert(eq_gun(0)==EQ_DEFAULT);assert(eq_pickup(0,PICKUP_SHOTGUN,false));assert(eq_gun(0)==EQ_SHOTGUN);assert(!eq_pickup(0,PICKUP_SHOTGUN,true));assert(!eq_pickup(0,PICKUP_RAILGUN,false));equipment[0].cooldown[EQ_SHOTGUN]=.8f;NetInputCommand c={.pressed=NET_INPUT_SWITCH};eq_command(0,&c,.01f);assert(eq_gun(0)==EQ_DEFAULT);eq_command(0,&c,.01f);assert(equipment[0].cooldown[EQ_SHOTGUN]==.8f);assert(eq_pickup(0,PICKUP_RAILGUN,true));players[0].matter=0;players[0].isExposed=true;eq_tick(.3f);assert(eq_gun(0)==EQ_RAIL);kill_player(0,1,false,false);assert(eq_gun(0)==EQ_DEFAULT);}
static void attacks(void){fixture();eq_pickup(0,PICKUP_RAILGUN,false);equipment[0].switchTime=0;NetInputCommand hold={.held=NET_INPUT_FIRE},release={0};eq_command(0,&hold,.65f);eq_command(0,&release,.01f);assert(players[1].matter==50 && players[1].respawn_timer==0 && players[0].matter==85);eq_tick(2);eq_command(0,&hold,.65f);eq_command(0,&release,.01f);assert(players[1].respawn_timer>0);
 fixture();eq_rocket_blast((EquipmentRocket){.owner=0,.pos={0,2,-3}});assert(players[1].respawn_timer>0);
 fixture();players[1].pos.z=-1;eq_pickup(0,PICKUP_SHOTGUN,false);equipment[0].switchTime=0;players[1].matter=5;assert(eq_fire(0));assert(players[1].isExposed && players[1].respawn_timer==0);assert(players[0].matter==94);assert(!eq_fire(0));
 fixture();eq_pickup(0,PICKUP_ROCKET,false);equipment[0].switchTime=0;for(int n=0;n<32;n++)eqRockets[n].active=true;assert(!eq_fire(0) && players[0].matter==100);memset(eqRockets,0,sizeof(eqRockets));players[0].matter=25;assert(eq_fire(0)&&players[0].isExposed&&eq_gun(0)==EQ_ROCKET);
 fixture();add_static_box_at_grid(-2,2,0,6,-4,-4,GRAY);init_static_hash();rebuild_voxel_hash();eq_pickup(0,PICKUP_CUTTER,false);equipment[0].switchTime=0;eq_command(0,&hold,.01f);assert(players[1].matter==100);eq_rocket_blast((EquipmentRocket){.owner=0,.pos={0,2,-1}});assert(players[1].matter==100);
 fixture();players[1].invuln_timer=1;eq_damage(0,1,100,true,DEATH_ROCKET,players[0].pos);assert(players[1].matter==100);gameMode=GAME_MODE_FIREFIGHT;combatantTeam[1]=0;players[1].invuln_timer=0;eq_damage(0,1,100,true,DEATH_RAIL,players[0].pos);assert(players[1].matter==100);
}
static void panels(void){fixture();eq_pickup(0,PICKUP_REINFORCED,false);players[0].pitch=-40;players[1].pos.x=20;Vector3 points[EQ_PANEL_CUBES];assert(eq_panel_positions(0,points));float cost=3*ff_build_cost(0);assert(eq_build(0));assert(voxel_count==EQ_PANEL_CUBES && players[0].matter==100-cost);for(int i=0;i<EQ_PANEL_CUBES;i++)assert(voxels[i].reinforced);float m=players[0].matter;assert(!eq_build(0)&&players[0].matter==m);}
static void sticky(void){fixture();GrenadeState g={.id=1,.owner=0,.pos={0,.2f,0},.vel={0,-1,0},.fuse=-1,.adhesive=true};clarity_grenade_step(&g,.1f);assert(g.attached && g.fuse>1.3f);for(int n=0;n<13;n++)clarity_grenade_step(&g,.1f);assert(g.fuse>0);clarity_grenade_step(&g,.11f);assert(g.fuse==0);
 int v=addVoxel(0,2,0,false,true,GRAY,0);g=(GrenadeState){.id=2,.attached=true,.adhesive=true,.support=voxels[v].identity,.attachment={.5f,.5f,.5f},.fuse=.4f};clarity_grenade_step(&g,.01f);Vector3 old=g.pos;for(int p=0;p<8;p++)voxels[v].particles[p]->pos.x+=2;clarity_grenade_step(&g,.01f);assert(fabsf(g.pos.x-old.x-2)<.001f);remove_voxel_index(v);clarity_grenade_step(&g,.01f);assert(!g.attached && g.fuse<.4f && g.fuse>0);}
static void supplies_regenerate(void){fixture();gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;useCustomMap=false;pickups[0]=(Pickup){.active=true,.type=PICKUP_RAILGUN,.pos={5,.8f,5}};add_static_box_at_grid(9,11,0,3,9,11,GRAY);init_static_hash();clarity_supply_update();assert(pickups[0].pos.y>2.7f);clear_world_voxels();init_static_hash();clarity_supply_update();assert(fabsf(pickups[0].pos.y-.8f)<.01f);}
static void maps(void){for(int map=0;map<2;map++){fixture();currentWorldType=map?WORLD_TYPE_CROSSCUT:WORLD_TYPE_TEST;gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=!map;firefightPractice=false;useCustomMap=false;multiplayerActivePlayers=2;firefightHumans=2;ResetGame();int counts[PICKUP_TYPE_COUNT]={0};for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active){counts[pickups[n].type]++;assert(v_isfinite(pickups[n].pos));}for(int t=PICKUP_SHOTGUN;t<=PICKUP_REINFORCED;t++){fprintf(stderr,"map%d type%d count%d\n",map,t,counts[t]);assert(counts[t]==(map?4:2));}}}
static void lifecycle(void){
 fixture();gameMode=GAME_MODE_FIREFIGHT;firefightHumans=2;combatantTeam[1]=0;firefightPractice=false;eq_pickup(0,PICKUP_ROCKET,false);eq_pickup(0,PICKUP_ADHESIVE,false);players[0].isExposed=true;kill_player(0,-1,false,false);assert(combatantLife[0]==LIFE_DOWNED&&eq_gun(0)==EQ_ROCKET&&eq_upgrade(0,EQ_ADHESIVE));ff_restore_human(0,35);assert(eq_gun(0)==EQ_ROCKET&&eq_upgrade(0,EQ_ADHESIVE));firefightWave=1;firefightEndless=true;firefightFoundry=false;ff_wave_complete();assert(eq_gun(0)==EQ_ROCKET);eq_reset_all();assert(eq_gun(0)==EQ_DEFAULT&&!equipment[0].upgrades);
}
static void sword(void){fixture();eq_pickup(0,PICKUP_SWORD,false);players[1].pos.z=-2;eq_sword(0);assert(players[1].respawn_timer>0);}
static void shotgun_structures(void){
 fixture();players[1].pos.z=-5;
 add_static_box_at_grid(-2,2,2,6,-6,-6,GRAY);init_static_hash();rebuild_voxel_hash();
 int loose=addVoxel(0,2,-6,false,true,BLUE,0);uint64_t hidden=voxels[loose].identity;rebuild_voxel_hash();
 eq_pickup(0,PICKUP_SHOTGUN,false);equipment[0].switchTime=0;int before=voxel_count;assert(eq_fire(0));
 int moving=0;float momentum=0;
 for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && voxels[v].identity!=hidden){
   moving++;assert(v_isfinite(voxels[v].vel) && voxels[v].vel.z<0 && voxels[v].activator==0);
   momentum+=-voxels[v].vel.z*voxel_mass_units(&voxels[v]);
 }
 assert(moving>0 && voxel_count==before && moving<=VOXEL_ACTIVATION_UNIT_BUDGET);
 assert(momentum>0 && momentum<=eqShotgunImpulse+.01f);
 assert(players[1].matter==100 && players[0].matter==94);
 assert(v_length(voxels[enemy_voxel(hidden)].vel)==0); // No impulse through the wall.
 float beforeZ=0;for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && voxels[v].identity!=hidden)beforeZ+=voxels[v].pos.z;
 for(int step=0;step<10;step++)simulate_voxel_pbd_steps(1.f/60,1);
 float afterZ=0;int afterCount=0;for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && voxels[v].identity!=hidden){afterZ+=voxels[v].pos.z;afterCount++;assert(v_isfinite(voxels[v].pos));}
 assert(afterCount==moving && afterZ<beforeZ-.01f);
 // A wide loose target receives a single volley's bounded impulse and wakes.
 fixture();players[1].pos.x=20;int target=addVoxel(0,2,-1,false,true,GRAY,0);voxels[target].sleeping=true;
 uint64_t id=voxels[target].identity;rebuild_voxel_hash();eq_pickup(0,PICKUP_SHOTGUN,false);equipment[0].switchTime=0;
 assert(eq_fire(0));target=enemy_voxel(id);assert(target>=0 && !voxels[target].sleeping && voxels[target].vel.z<0);
 assert(v_length(voxels[target].vel)<=eqShotgunImpulse+.01f);
 // Ordinary pellets still stop at attached armor, without ripping it free.
 fixture();players[1].enemyType=ENEMY_TYPE_GOLIATH;players[1].pos.z=-20;
 target=addVoxel(0,2,-1,false,true,GRAY,0);id=voxels[target].identity;
 goliathStates[1].count=1;goliathStates[1].voxels[0].voxelIdentity=id;rebuild_voxel_hash();
 eq_pickup(0,PICKUP_SHOTGUN,false);equipment[0].switchTime=0;assert(eq_fire(0));
 assert(goliathStates[1].count==1 && v_length(voxels[enemy_voxel(id)].vel)==0 && players[1].matter==100);
}
static void protection(void){
 fixture();players[1].enemyType=ENEMY_TYPE_GOLIATH;int armor=addVoxel(0,2,-1.5f,false,true,GRAY,0);goliathStates[1].count=1;goliathStates[1].voxels[0].voxelIdentity=voxels[armor].identity;rebuild_voxel_hash();eq_pickup(0,PICKUP_RAILGUN,false);equipment[0].switchTime=0;NetInputCommand hold={.held=NET_INPUT_FIRE},release={0};eq_command(0,&hold,.65f);eq_command(0,&release,.01f);assert(players[1].matter==50&&goliathStates[1].count==0);
 fixture();activePlayers=5;playerInput[4]=INPUT_TYPE_BOT_HARD;enemy_lab_spawn(4,ENEMY_TYPE_FOUNDRY,(Vector3){0,2,-5});float health=enemySpecial[4].coreHealth;eq_damage(0,4,100,true,DEATH_ROCKET,players[0].pos);assert(enemySpecial[4].coreHealth==health);enemySpecial[4].coreOpen=2;eq_damage(0,4,50,true,DEATH_RAIL,players[0].pos);assert(enemySpecial[4].coreHealth==health-50);
 fixture();worldPeriodic=true;worldFloorHalfSize=32;players[1].pos.z=10;eq_pickup(0,PICKUP_RAILGUN,false);equipment[0].switchTime=0;eq_command(0,&hold,.65f);eq_command(0,&release,.01f);assert(players[1].matter==50); // Behind locally, ahead in the adjacent periodic image; hit only once.
}
static void stress(void){
 fixture();gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;useCustomMap=false;multiplayerActivePlayers=4;firefightHumans=4;currentWorldType=WORLD_TYPE_TEST;ResetGame();gameState=GAME_STATE_PLAYING;activePlayers=4;for(int i=0;i<4;i++){players[i].pitch=-15;players[i].invuln_timer=999;equipment[i].gun[1]=EQ_ROCKET;equipment[i].selected=1;equipment[i].upgrades=7;}
 double worst=0,first=0,last=0;int frames=600;
 for(int frame=0;frame<frames;frame++){double start=GetTime();for(int i=0;i<4;i++){players[i].matter=100;players[i].isExposed=false;if(frame%90==0){equipment[i].selected=1;eq_fire(i);}if(frame%90>30&&frame%90<60){equipment[i].gun[0]=EQ_CUTTER;equipment[i].selected=0;NetInputCommand c={.held=NET_INPUT_FIRE};eq_command(i,&c,1.f/60);}if(frame%90==60){players[i].last_build_time=-100;float pitch=players[i].pitch;players[i].pitch=-45;eq_build(i);players[i].pitch=pitch;}}eq_tick(1.f/60);simulate_voxel_pbd_steps(1.f/60,1);int rockets=0;for(int n=0;n<32;n++)rockets+=eqRockets[n].active;assert(rockets<=32);for(int n=0;n<voxel_count;n++)assert(v_isfinite(voxels[n].pos));double ms=(GetTime()-start)*1000;worst=fmax(worst,ms);if(frame<60)first+=ms;if(frame>=frames-60)last+=ms;}
 for(int n=0;n<100;n++)eq_tick(.1f);for(int n=0;n<32;n++)assert(!eqRockets[n].active);fprintf(stderr,"Four-player Foundry rocket/cutter/panel stress first %.3f last %.3f worst %.3f ms, voxels %d\n",first/60,last/60,worst,voxel_count);
 firefightWaveState=FIREFIGHT_STATE_ACTIVE;firefightWaveBannerTimer=0;RenderTexture2D views[4]={0};int count=0,w=0,h=0;multiplayerActivePlayers=4;drawWeaponModels=true;render_gameplay_view(views,&count,&w,&h,false);assert(count==4);for(int i=0;i<count;i++){Image im=LoadImageFromTexture(views[i].texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/equipment-four-%d.png",i));UnloadImage(im);UnloadRenderTexture(views[i]);}
}
static void captures(void){fixture();weapon_renderer_init();character_renderer_init();for(int size=0;size<2;size++){int w=size?480:960,h=size?270:540;RenderTexture2D rt=LoadRenderTexture(w,h);weaponViewAspect=(float)w/h;drawViewPlayerIndex=0;weaponFrameSampled=false;for(int gun=0;gun<9;gun++){equipment[0].gun[0]=gun<=6?gun:0;equipment[0].upgrades=gun==8?EQ_BUILD:0;equipment[0].buildTime=gun==8?.4f:0;players[0].meleeSwingActive=gun==6;players[0].meleeSwingStartTime=GetTime()-.5f*MELEE_ANIM_DURATION_SECONDS;weaponVisuals[0].charge=.8f;Camera3D c={players[0].pos,v_add(players[0].pos,(Vector3){0,0,-1}),{0,1,0},60,CAMERA_PERSPECTIVE};BeginTextureMode(rt);ClearBackground((Color){30,40,55,255});BeginMode3D(c);WeaponPose p=player_weapon_pose(0,true);if(gun==7)p.kind=WEAPON_ADHESIVE_GRENADE;weapon_draw(p,&weaponVisuals[0],c.position,GetTime(),false,WEAPON_COLOR);weapon_draw(p,&weaponVisuals[0],c.position,GetTime(),false,WEAPON_GLASS);character_draw_paws(p,BLUE,false);EndMode3D();eq_hud(0,w,h);DrawCircle(w/2,h/2,2,WHITE);EndTextureMode();Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/equipment-%d-%d.png",gun,size));UnloadImage(im);}UnloadRenderTexture(rt);}character_renderer_shutdown();weapon_renderer_shutdown();}
static void revised_weapons(void){
 fixture();players[1].pos.z=-6;add_static_box_at_grid(-2,2,1,6,-8,-4,GRAY);init_static_hash();rebuild_voxel_hash();
 eq_pickup(0,PICKUP_RAILGUN,false);equipment[0].switchTime=0;
 NetInputCommand held={.held=NET_INPUT_FIRE},release={0};eq_command(0,&held,.65f);eq_command(0,&release,.01f);
 assert(players[1].matter==50);VoxelHit hit;assert(!first_voxel_hit_detailed((Ray){players[0].pos,{0,0,-1}},10,-1,false,&hit));
 assert(eqTraces[0].life>0 && v_length(v_sub(eqTraces[0].start,players[0].pos))>.2f);
 fixture();players[1].pos.x=20;add_static_box_at_grid(-1,1,1,6,-20,-2,GRAY);init_static_hash();rebuild_voxel_hash();
 eq_pickup(0,PICKUP_CUTTER,false);equipment[0].switchTime=0;eq_command(0,&held,.25f);assert(fabsf(players[0].matter-99)<.001f);
 assert(eqCutterUnitSeconds==1.f/12);
 fixture();eq_pickup(0,PICKUP_SWORD,false);equipment[0].switchTime=0;assert(eq_gun(0)==EQ_BLADE && !(equipment[0].upgrades&EQ_SWORD));
 players[0].matter=0;players[0].isExposed=true;players[0].last_melee_time=-100;assert(eq_fire(0));assert(players[0].meleeSwingActive);
 fixture();players[1].pos.x=20;add_static_box_at_grid(-5,5,0,8,-7,-3,GRAY);init_static_hash();rebuild_voxel_hash();int before=voxel_count;
 eq_sword(0);assert(voxel_count<before);for(int v=0;v<voxel_count;v++)assert(v_isfinite(voxels[v].pos));
}
static void beam_captures(void){
 fixture();weapon_renderer_init();weaponFrameSampled=false;drawViewPlayerIndex=0;
 for(int small=0;small<2;small++)for(int gun=EQ_CUTTER;gun<=EQ_RAIL;gun++){
  int w=small?480:960,h=small?270:540;weaponViewAspect=(float)w/h;
  equipment[0].gun[0]=gun;equipment[0].held=gun==EQ_CUTTER;equipment[0].beamEnd=(Vector3){0,2,-12};
  memset(eqTraces,0,sizeof(eqTraces));if(gun==EQ_RAIL)eq_trace(player_weapon_pose(0,false).muzzle,(Vector3){0,2,-100},gun,0);
  RenderTexture2D rt=LoadRenderTexture(w,h);Camera3D cam={players[0].pos,{0,2,-20},{0,1,0},60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){24,32,45,255});BeginMode3D(cam);
  DrawPlane((Vector3){0},(Vector2){60,60},DARKGRAY);eq_world();weapon_draw(player_weapon_pose(0,true),&weaponVisuals[0],cam.position,GetTime(),false,WEAPON_COLOR);
  EndMode3D();EndTextureMode();assert(eqBeamShader.id && eqBeamShader.id!=rlGetShaderIdDefault());
  Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/revised-beam-%d-%d.png",gun,small));UnloadImage(im);UnloadRenderTexture(rt);
 }
 eq_beam_shutdown();weapon_renderer_shutdown();
}
int main(int argc,char**argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Equipment tests");init_pbd_thread_pool();assert(initialize_physics_backend());inventory();revised_weapons();attacks();panels();sticky();lifecycle();sword();shotgun_structures();protection();supplies_regenerate();maps();stress();captures();beam_captures();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Equipment tests passed");return 0;}
