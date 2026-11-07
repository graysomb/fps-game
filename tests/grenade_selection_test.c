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

static void selection(void){
 fixture(1);grenadeInventory[0]=2;ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_REGULAR);
 assert(eq_pickup(0,PICKUP_ADHESIVE,false));ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_STICKY);
 ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_REGULAR);
 assert(clarity_throw(0));assert(grenades[CLARITY_GRENADES-1].kind==GRENADE_REGULAR && !grenades[CLARITY_GRENADES-1].adhesive && grenadeInventory[0]==1);
 grenadeCooldown[0]=0;ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_STICKY);assert(clarity_throw(0));
 assert(grenades[CLARITY_GRENADES-2].kind==GRENADE_REGULAR && grenades[CLARITY_GRENADES-2].adhesive && grenadeInventory[0]==0);
 grenadeCooldown[0]=0;assert(!clarity_throw(0));
 exGrenades[0][0]=exGrenades[0][1]=1;ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_GRAVITY);
 ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_IGLOO);
 grenadeInventory[0]=2;ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_REGULAR);ex_cycle(0);assert(exSelectedGrenade[0]==GRENADE_STICKY);
 fp_capture();fp_restore();assert(exSelectedGrenade[0]==GRENADE_STICKY && grenadeInventory[0]==2);
 equipment[0].upgrades&=~EQ_ADHESIVE;assert(ex_grenade_count(0,GRENADE_STICKY)==0);ex_cycle(0);assert(exSelectedGrenade[0]!=GRENADE_STICKY);
 fp_queue_unlock(PICKUP_GRENADE);grenadeInventory[0]=0;assert(ex_refill_grenades(0));assert(grenadeInventory[0]==2);
 ex_reset_player(0);assert(exSelectedGrenade[0]==GRENADE_REGULAR);
 puts("Regular/sticky cycle, shared ammo, throw flags, unavailable modes, refill and checkpoint passed");
}
static void network(void){
 fixture(1);uint8_t packet[1024];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_u8(&w,0);
 for(int i=0;i<4;i++){net_write_u8(&w,i==0?GRENADE_STICKY:GRENADE_REGULAR);net_write_u8(&w,0);net_write_u8(&w,0);net_write_f32(&w,0);}
 for(int i=0;i<MAX_COMBATANTS;i++)net_write_f32(&w,0);
 for(int i=0;i<EX_FIELDS+8;i++)net_write_u32(&w,0);
 assert(!w.failed);NetReader r;net_reader_init(&r,packet,w.length);netTransport.role=NET_ROLE_CLIENT;
 assert(ex_receive(NET_MSG_EXPANSION,&r,1) && !r.failed && exSelectedGrenade[0]==GRENADE_STICKY);
 netTransport.role=NET_ROLE_OFFLINE;puts("LAN snapshot accepts sticky selection");
}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Grenade selection regression");init_pbd_thread_pool();assert(initialize_physics_backend());
 selection();network();enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
