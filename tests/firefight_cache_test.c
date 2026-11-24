#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(int humans){
    netTransport.role=NET_ROLE_OFFLINE;gameMode=GAME_MODE_FIREFIGHT;
    firefightPractice=enemyLab=equipmentLab=false;firefightFoundry=true;firefightEndless=true;firefightCheckpoints=false;
    useCustomMap=false;multiplayerActivePlayers=humans;
    for(int i=0;i<humans;i++)multiplayerPlayerInput[i]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;
    firefightWaveEnemiesRemaining=firefightWaveEnemiesTotal=100;
}
static void move_away(void){for(int i=0;i<firefightHumans;i++)players[i].pos=(Vector3){25,3,25};}
static size_t ground_item(int type,Vector3 pos){assert(fp_reserve(fpCount+1));size_t n=fpCount++;fpItems[n]=(FfItem){.id=fpNextId++,.type=type,.owner=-1,.slot=-1,.pos=pos};return n;}
static void timers(void){
 fixture(1);move_away();
 for(int r=0;r<FP_RULE_COUNT;r++)ground_item(fpRules[r].type,(Vector3){0,3,0});
 fp_tick(59);for(size_t n=0;n<fpCount;n++)assert(fpItems[n].owner==-1);
 fp_tick(1);for(size_t n=0;n<fpCount;n++)assert(fpItems[n].owner==-2);
 size_t n=ground_item(PICKUP_GLOCK,players[0].pos);assert(fp_collect(0,n,false));fp_tick(100);assert(fpItems[n].owner==0);
 fp_drop(&fpItems[n],(Vector3){0,3,0});assert(fpItems[n].groundAge==0);fp_tick(59);assert(fpItems[n].owner==-1);fp_tick(1);assert(fpItems[n].owner==-2);
 puts("Drop timers: all equipment types, held items, fresh death/replacement timers passed");
}
static void duplicates(void){
 fixture(1);size_t original=ground_item(PICKUP_GLOCK,players[0].pos);assert(fp_collect(0,original,false));
 for(int slot=0;slot<2;slot++){
  if(slot==1){int g=equipment[0].gun[0];equipment[0].gun[0]=equipment[0].gun[1];equipment[0].gun[1]=g;}
  players[0].matter=40;size_t n=ground_item(PICKUP_GLOCK,players[0].pos);fp_tick(.01f);assert(fpItems[n].owner==-2 && players[0].matter==65 && fpItems[original].owner==0);
  fp_tick(.01f);assert(players[0].matter==65);
 }
 players[0].matter=player_max_matter(&players[0])-5;ground_item(PICKUP_GLOCK,players[0].pos);fp_tick(.01f);assert(players[0].matter==player_max_matter(&players[0]));
 players[0].matter=0;players[0].isExposed=true;ground_item(PICKUP_GLOCK,players[0].pos);fp_tick(.01f);assert(players[0].matter==25 && !players[0].isExposed);
 puts("Duplicate guns: both slots, full loadout, single consumption, matter cap and exposed recovery passed");
}
static void cache(void){
 fixture(2);assert(fpCacheReady && !fpUnlocked);players[0].pos=fpCachePositions[fp_rule(PICKUP_GLOCK)];assert(!fp_cache_collect(0));
 size_t n=ground_item(PICKUP_GLOCK,players[0].pos);assert(fp_collect(0,n,false));int rule=fp_rule(PICKUP_GLOCK);assert(fpUnlocked&(1u<<rule));
 players[1].pos=fpCachePositions[rule];assert(fp_cache_near(1)==rule);size_t count=fpCount;
 fp_tick(.01f);assert(fpCount==count);assert(fp_cache_collect(1));assert(eq_gun(1)==EQ_GLOCK || equipment[1].gun[0]==EQ_GLOCK || equipment[1].gun[1]==EQ_GLOCK);
 assert(!fp_cache_collect(1)); // No automatic duplicate conversion from the cache.
 fp_queue_unlock(PICKUP_GLOB);players[1].pos=fpCachePositions[fp_rule(PICKUP_GLOB)];assert(fp_cache_collect(1));
 fp_queue_unlock(PICKUP_ADHESIVE);players[1].pos=fpCachePositions[fp_rule(PICKUP_ADHESIVE)];assert(!fp_cache_collect(1));fp_tick(.01f);assert(eq_upgrade(1,EQ_ADHESIVE));assert(!fp_cache_collect(1));
 move_away();fp_tick(120);assert(fpCacheReady && (fpUnlocked&(1u<<rule)));
 fp_capture();uint32_t unlocked=fpUnlocked;fp_restore();assert(fpUnlocked==unlocked && fpCacheReady);
 netTransport.role=NET_ROLE_CLIENT;players[0].pos=fpCachePositions[rule];assert(!fp_cache_collect(0));netTransport.role=NET_ROLE_OFFLINE;
 fp_new_run();assert(!fpUnlocked && !fpCacheReady);
 puts("Cache: collected-only unlocks, shared team supplies, manual pickup/replacement, upgrades, checkpoint/reset and authority passed");
}
static void supplies_and_glass(void) {
 fixture(2);
 int rule=fp_rule(PICKUP_GRENADE);
 for(int i=0;i<2;i++)players[i].pos=fpCachePositions[rule];
 fp_tick(.01f);assert(!grenadeInventory[0] && !grenadeInventory[1]);
 fp_queue_unlock(PICKUP_GRENADE);fp_tick(.01f);
 assert(grenadeInventory[0]==1 && grenadeInventory[1]==1);
 grenadeInventory[0]=0;fp_tick(.01f);assert(!grenadeInventory[0]);
 move_away();fp_tick(.01f);players[0].pos=fpCachePositions[rule];fp_tick(.01f);
 assert(grenadeInventory[0]==1);
 int glass=0;bool tested=false;
 for(int v=0;v<voxel_count;v++)if(spawn_glass(v)){
  glass++;
  for(int axis=0;axis<2 && !tested;axis++){
   Vector3 dir=axis?(Vector3){0,0,1}:(Vector3){1,0,0};
   Vector3 a=v_sub(voxels[v].pos,dir),b=v_add(voxels[v].pos,dir);
   VoxelHit hit;
   if(!first_voxel_hit_detailed((Ray){a,dir},2,-1,true,&hit)||hit.id!=v)continue;
   if(!first_voxel_hit_detailed((Ray){b,v_mul(dir,-1)},2,-1,true,&hit)||hit.id!=v)continue;
   assert(!is_view_occluded_by_voxels(a,b));
   players[4].pos=a;assert(!rifle_view_blocked(4,b));
   assert(player_aabb_hits_world_ex(voxels[v].pos,.25f,4));
   tested=true;
  }
 }
 assert(glass>0 && tested);
 // On a flat arena, cross each opening using the same body sweeps as enemy navigation.
 clear_world_voxels();init_static_hash();players[0].pos=(Vector3){0,1,0};firefightHumans=1;fp_cache_place();
 Vector3 low=fpCacheOrigin,high=fpCacheOrigin;
 for(int v=0;v<voxel_count;v++)if(spawn_glass(v)){low=Vector3Min(low,voxels[v].pos);high=Vector3Max(high,voxels[v].pos);}
 Vector3 middle=v_mul(v_add(low,high),.5f);middle.y=1;
 players[4].enemyType=ENEMY_TYPE_STANDARD;players[4].respawn_timer=0;
 for(int axis=0;axis<2;axis++)for(int side=0;side<2;side++){
  Vector3 edge=middle,dir=axis?(Vector3){0,0,1}:(Vector3){1,0,0};
  if(axis)edge.z=side?high.z:low.z;else edge.x=side?high.x:low.x;
  assert(ai_nav_clear(4,v_sub(edge,dir),v_add(edge,dir)));
  assert(!first_voxel_hit_detailed((Ray){v_sub(edge,dir),dir},2,-1,true,NULL));
 }
 Vector3 wall={low.x,1,low.z+1};
 players[4].pos=v_add(wall,(Vector3){-1.5f,0,0});Vector3 goal=v_add(wall,(Vector3){2,0,0});
 assert(!ai_nav_clear(4,players[4].pos,goal));
 ai_nav_search(4,goal);assert(aiNavigation[4].count>0);
 assert(Vector3Distance(aiNavigation[4].points[aiNavigation[4].count-1],goal)<1);
 puts("Walk-up supplies, independent entry latches, glass collision and sight passed");
}
static void controller_cache(void){
 fixture(1);debug_console_exec("wave 4");clear_pickups();
 int rule=fp_rule(PICKUP_GLOCK);players[0].pos=v_add(fpCachePositions[rule],(Vector3){0,.2f,0});
 playerInput[0]=INPUT_TYPE_GAMEPAD;aimStates[0]=(AimState){.tapArmed=true};NetInputCommand c={0};
 assert(fp_near(0,false)<0 && fp_cache_near(0)==rule);
 // Merge regression: remappable combined interaction must favor revival over the cache.
 firefightHumans=2;combatantLife[1]=LIFE_DOWNED;players[1].death_pos=players[0].pos;
 aim_interact(0,true,.26f,&c);assert(c.held&NET_INPUT_REVIVE);assert(!(c.pressed&NET_INPUT_EQUIP));
 firefightHumans=1;combatantLife[1]=LIFE_UNUSED;aimStates[0]=(AimState){.tapArmed=true};c=(NetInputCommand){0};
 aim_interact(0,true,.26f,&c);assert(c.pressed&NET_INPUT_EQUIP);assert(!(c.held&NET_INPUT_REVIVE));
 eq_command(0,&c,.016f);assert(equipment[0].gun[0]==EQ_GLOCK || equipment[0].gun[1]==EQ_GLOCK);
 size_t count=fpCount;c=(NetInputCommand){0};aim_interact(0,true,.3f,&c);assert(!(c.pressed&NET_INPUT_EQUIP) && fpCount==count);
 aim_interact(0,false,.01f,&c);
 // LAN clients must also generate the interaction for the authoritative host.
 eq_reset(0);aimStates[0].tapArmed=true;netTransport.role=NET_ROLE_CLIENT;c=(NetInputCommand){0};aim_interact(0,true,.26f,&c);assert(c.pressed&NET_INPUT_EQUIP);
 netTransport.role=NET_ROLE_OFFLINE;
 puts("Controller cache: hold X generates one pickup, equips weapon, and works for LAN input passed");
}
static void console_unlocks(void){
 fixture(1);start_firefight_wave(5);assert(!fpUnlocked); // Normal progression still requires collection.
 debug_console_exec("wave 4");
 uint32_t early=(1u<<fp_rule(PICKUP_GLOB))|(1u<<fp_rule(PICKUP_REINFORCED))|(1u<<fp_rule(PICKUP_GLOCK))|(1u<<fp_rule(PICKUP_GRENADE));
 assert(firefightWave==4 && fpUnlocked==early && fpCheckpoint.unlocked==early);
 players[0].pos=fpCachePositions[fp_rule(PICKUP_GLOCK)];assert(fp_cache_collect(0));
 debug_console_exec("nextwave");assert(firefightWave==5 && (fpUnlocked&(1u<<fp_rule(PICKUP_ADHESIVE))) && (fpUnlocked&(1u<<fp_rule(PICKUP_SHIELD_REGEN))));
 assert(!(fpUnlocked&(1u<<fp_rule(PICKUP_SHOTGUN))));
 uint32_t unlocked=fpUnlocked;int banners=fpBannerWrite;debug_console_exec("wave 5");assert(fpUnlocked==unlocked && fpBannerWrite==banners);
 fp_restore();assert(fpUnlocked==unlocked);
 debug_console_exec("wave 1");assert(fpUnlocked==unlocked); // Going back never removes earned unlocks.
 debug_console_exec("wave 12");assert(fpUnlocked==(((1u<<FP_RULE_COUNT)-1)&~(1u<<fp_rule(PICKUP_IGLOO))));
 netTransport.role=NET_ROLE_CLIENT;debug_console_exec("wave 5");assert(firefightWave==12);netTransport.role=NET_ROLE_OFFLINE;
 fixture(1);debug_console_exec("wave 1");assert(!fpUnlocked);
 puts("Console skips: earlier-wave unlocks, cache claims, nextwave, checkpoint persistence, idempotence and host authority passed");
}
static void glass_visible_at_rest(void) {
 clear_world_voxels();clear_all_chunks();init_static_hash();
 int id=add_static_voxel_at_grid(0,2,0,spawnGlassColor,0);
 assert(id>=0);init_static_hash();
 RenderTexture2D rt=LoadRenderTexture(80,80);
 Vector3 center=voxels[id].pos;
 const Vector3 views[]={{2,0,0},{-2,0,0},{0,0,2},{0,0,-2},{0,2,0},{0,-2,0}};
 for(int i=0;i<6;i++){
  Camera3D cam={.position=v_add(center,views[i]),.target=center,.up=i<4?(Vector3){0,1,0}:(Vector3){0,0,1},.fovy=45};
  BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(cam);
  DrawVoxels(cam);draw_spawn_glass();
  EndMode3D();EndTextureMode();
  Image image=LoadImageFromTexture(rt.texture);Color pixel=GetImageColor(image,40,40);
  assert(pixel.b>60 && pixel.b<150 && pixel.b>pixel.r+40);
  UnloadImage(image);
 }
 UnloadRenderTexture(rt);
 puts("Unhit glass transmits the background and shows its tint from all six sides");
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Firefight cache regression");init_pbd_thread_pool();assert(initialize_physics_backend());
 timers();duplicates();cache();supplies_and_glass();glass_visible_at_rest();console_unlocks();controller_cache();
 fixture(1);assert(weapon_renderer_init());fpUnlocked=(1u<<FP_RULE_COUNT)-1;drawViewPlayerIndex=0;
 Camera3D camera={.position=v_add(fpCacheOrigin,(Vector3){13,14,20}),.target=v_add(fpCacheOrigin,(Vector3){0,0,7}),.up={0,1,0},.fovy=55};
 RenderTexture2D preview=LoadRenderTexture(960,540);BeginTextureMode(preview);ClearBackground((Color){25,29,35,255});BeginMode3D(camera);DrawPlane((Vector3){0},(Vector2){80,80},DARKGRAY);fp_draw(WEAPON_COLOR);fp_draw(WEAPON_GLASS);EndMode3D();DrawText("FIREFIGHT SPAWN CACHE",20,20,24,SKYBLUE);EndTextureMode();Image previewImage=LoadImageFromTexture(preview.texture);ImageFlipVertical(&previewImage);ExportImage(previewImage,".build/firefight-cache-preview.png");UnloadImage(previewImage);UnloadRenderTexture(preview);
 enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Firefight cache checks passed");return 0;
}
