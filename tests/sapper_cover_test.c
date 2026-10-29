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
 players[0].pos=(Vector3){0,1,-12};players[0].respawn_timer=0;players[0].matter=100;combatantLife[0]=LIFE_ALIVE;
 init_firefight_bot_entity(4,ENEMY_TYPE_SAPPER,INPUT_TYPE_BOT_EASY);players[4].pos=(Vector3){0,1,0};players[4].yaw=0;enemySpecial[4].cooldown=0;
 ai_select(4,0,players[0].pos,AI_ALERT_SEEN);aiCombat[4].reaction=0;
}
static void arc_and_cover(void){
 fresh();add_static_box_at_grid(-6,6,0,3,-7,-6,GRAY);init_static_hash();rebuild_voxel_hash();
 Vector3 velocity;assert(sapper_arc(players[4].pos,players[0].pos,&velocity));assert(velocity.y>0);
 assert(sapper_find_cover(4));assert(sapper_hidden(4,enemySpecial[4].buildOrigin,players[0].pos,NULL));
 assert(!player_aabb_hits_world_ex(enemySpecial[4].buildOrigin,player_collision_half(&players[4]),4));
 players[4].pos=enemySpecial[4].buildOrigin;
 for(int n=0;n<150 && !enemySpecial[4].actions;n++){UpdateBot(4,1.f/60);net_integrate_player(&players[4],1.f/60,true);}
 assert(enemySpecial[4].actions==1);int index=-1;for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id)index=n;assert(index>=0);
 GrenadeState g=grenades[index];bool bounced=false;float flight=0;
 while(!bounced && flight<4){clarity_grenade_step(&g,1.f/120);flight+=1.f/120;bounced=g.fuse>=0;}
 fprintf(stderr,"Sapper lob: first bounce %.2fs, miss %.2fm, vertical launch %.2f\n",flight,v_length((Vector3){g.pos.x-players[0].pos.x,0,g.pos.z-players[0].pos.z}),grenades[index].vel.y);
 assert(bounced && v_length((Vector3){g.pos.x-players[0].pos.x,0,g.pos.z-players[0].pos.z})<1);
 assert(flight+g.fuse>=1.49f);
 // A new ceiling blocks the launch; pool and action count remain unchanged.
 memset(grenades,0,sizeof(grenades));add_static_box_at_grid(-12,12,5,6,-12,4,GRAY);init_static_hash();
 assert(!sapper_arc(players[4].pos,players[0].pos,&velocity));
 clear_world_voxels();init_static_hash();enemySpecial[4].searchTimer=1;UpdateBot(4,.01f);assert(!enemySpecial[4].selected);
 puts("Cover selection, real hidden throw, ballistic landing, fuse and obstruction checks passed");
}
static void supplies_and_safety(void){
 fresh();Vector3 v;assert(sapper_arc(players[4].pos,players[0].pos,&v));
 grenadeInventory[4]=2;equipment[4].upgrades|=EQ_ADHESIVE;fpCarrierType[4]=PICKUP_ADHESIVE;
 assert(sapper_throw(4,v));assert(grenadeInventory[4]==1);
 int g=-1;for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id)g=n;assert(g>=0 && grenades[g].adhesive);
 memset(grenades,0,sizeof(grenades));exSelectedGrenade[4]=GRENADE_GRAVITY;exGrenades[4][GRENADE_GRAVITY-1]=1;
 assert(sapper_throw(4,v));assert(!exGrenades[4][GRENADE_GRAVITY-1]);
 memset(grenades,0,sizeof(grenades));grenadeInventory[4]=0;exSelectedGrenade[4]=GRENADE_REGULAR;
 assert(sapper_throw(4,v));assert(grenadeInventory[4]==0);
 for(int n=0;n<CLARITY_GRENADES;n++)grenades[n].id=n+1;
 grenadeInventory[4]=2;assert(!sapper_throw(4,v) && grenadeInventory[4]==2);
 memset(grenades,0,sizeof(grenades));players[0].pos.z=-3;assert(!sapper_arc(players[4].pos,players[0].pos,&v));
 ai_select(4,0,players[0].pos,AI_ALERT_SEEN);UpdateBot(4,.1f);assert(players[4].vel.z>0 && !enemySpecial[4].windup);
 netTransport.role=NET_ROLE_CLIENT;int actions=enemySpecial[4].actions;UpdateBot(4,1);assert(actions==enemySpecial[4].actions);
 netTransport.role=NET_ROLE_OFFLINE;players[4].isExposed=true;UpdateBot(4,.1f);assert(!enemySpecial[4].windup);
 puts("Carrier supplies, innate bombs, allocation failure, close retreat, exposed and client safety passed");
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Sapper cover tests");
 init_pbd_thread_pool();assert(initialize_physics_backend());arc_and_cover();supplies_and_safety();enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Sapper cover tests passed");return 0;
}
