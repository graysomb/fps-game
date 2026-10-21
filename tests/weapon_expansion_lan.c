#define FIREFIGHT_TESTING
#include "raylib.h"
static bool focused(void){return true;}
#define IsWindowFocused focused
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static bool joined,sawGlob,sawDmr,sawGravity,sawIgloo,sawSuppression,sawCycle,sawReset,sawLastPickup;
static void receive(NetTransport *t,int peer,uint8_t channel,const uint8_t *data,size_t length,void *user){
 assert(length<=FPS_NET_MAX_PACKET);NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);assert(net_read_header(&r,&h));net_on_receive(t,peer,channel,data,length,user);
 if(t->role==NET_ROLE_HOST && h.type==NET_MSG_HELLO && netPeerPlayerMask[peer]){
  players[1].pos=(Vector3){0,2,0};players[1].respawn_timer=0;players[1].invuln_timer=0;players[1].matter=100;players[1].isExposed=false;eq_reset(1);equipment[1].gun[1]=EQ_GLOB;equipment[1].selected=1;
  grenadeInventory[1]=0;exGrenades[1][0]=exGrenades[1][1]=2;exSelectedGrenade[1]=GRENADE_GRAVITY;exSuppression[1]=3;
  pickups[MAX_PICKUPS-1]=(Pickup){.type=PICKUP_DMR,.active=true,.pos={30,.8f,30}};joined=true;
 }
 if(t->role==NET_ROLE_CLIENT){
  sawLastPickup|=pickups[MAX_PICKUPS-1].active&&pickups[MAX_PICKUPS-1].type==PICKUP_DMR;
  sawSuppression|=exSuppression[1]>0;sawCycle|=exSelectedGrenade[1]==GRENADE_IGLOO;
  for(int n=0;n<EX_GLOBS;n++)sawGlob|=exGlobs[n].id!=0;
  for(int n=0;n<EX_FIELDS;n++)sawGravity|=exFields[n].id!=0;
  for(int n=0;n<8;n++)sawIgloo|=exReveals[n].id!=0;
  for(int n=0;n<32;n++)sawDmr|=eqTraces[n].life>0&&eqTraces[n].gun==EQ_DMR;
  sawReset|=sawGravity&&equipment[1].gun[1]==EQ_EMPTY&&exGrenades[1][0]==0&&exGrenades[1][1]==0;
 }
}
int main(int argc,char **argv){
 bool host=argc>1&&!strcmp(argv[1],"host");SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Expansion LAN");gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();gameState=GAME_STATE_PLAYING;assert(net_transport_global_init());
 if(host){assert(net_transport_host(&netTransport,27940));net_set_host_local_players(1);netLobbyStarted=true;players[0].pos=(Vector3){20,2,20};}else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27940));}
 double start=GetTime(),begin=-1,previous=start;uint32_t seq=0;bool cycled=false,igloo=false,gravity=false,switched=false,fired=false,reset=false,given=false;
 while(GetTime()-start<12){net_transport_pump(&netTransport,receive,net_on_connect,NULL);double now=GetTime();float dt=now-previous;previous=now;
  if((host&&joined)||(!host&&netWorldReady)){if(begin<0)begin=now;float age=now-begin;
   if(host){clarity_tick(dt);eq_tick(dt);net_host_tick();if(age>3.5&&!given){equipment[1].gun[0]=EQ_DMR;given=true;}if(age>8.5&&!reset){assert(exGrenades[1][1]==1&&exGrenades[1][0]==1);assert(weaponVisuals[1].shot_sequence>1);kill_player(1,0,true,false);reset=true;}}
   else{NetInputCommand c={.sequence=++seq,.mouse_look=true,.aim_yaw=0,.aim_pitch=-20};if(age>.3&&!cycled){c.pressed|=NET_INPUT_GRENADE_CYCLE;cycled=true;}if(age>.8&&!igloo){c.pressed|=NET_INPUT_GRENADE;igloo=true;}if(age>1.4&&age<2.5)c.held|=NET_INPUT_FIRE;
    if(age>3&&!gravity){c.pressed|=NET_INPUT_GRENADE_CYCLE|NET_INPUT_GRENADE;gravity=true;}if(age>4&&!switched){c.pressed|=NET_INPUT_SWITCH;switched=true;}if(age>4.5&&!fired){c.pressed|=NET_INPUT_FIRE;fired=true;}
    uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);net_write_input(&w,&c);net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,c.pressed!=0);if(age>10)break;
   }
  }WaitTime(1./60);
 }
 if(host){assert(joined&&reset);puts("Expansion LAN host passed");}else{fprintf(stderr,"glob%d dmr%d gravity%d igloo%d suppression%d cycle%d reset%d pickup%d\n",sawGlob,sawDmr,sawGravity,sawIgloo,sawSuppression,sawCycle,sawReset,sawLastPickup);assert(sawGlob&&sawDmr&&sawGravity&&sawIgloo&&sawSuppression&&sawCycle&&sawReset&&sawLastPickup);puts("Expansion LAN client passed");}
 net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();return 0;
}
