#define FIREFIGHT_TESTING
#include "raylib.h"
static bool focused(void){return true;}
#define IsWindowFocused focused
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static bool joined,sawItems,sawGun,sawCharge,sawCarrier,sawRestart,sawRestore,sawReady,sawRejected,sawLoot;
static uint32_t initialSession;
static size_t savedCount;
static void receive(NetTransport*t,int peer,uint8_t channel,const uint8_t*data,size_t length,void*u){
 NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);if(!net_read_header(&r,&h))return;
 net_on_receive(t,peer,channel,data,length,u);
 if(t->role==NET_ROLE_HOST&&h.type==NET_MSG_HELLO&&netPeerPlayerMask[peer]){
  net_host_start_match();gameState=GAME_STATE_PLAYING;firefightWave=6;players[0].pos=(Vector3){20,1,20};players[1].pos=(Vector3){-36,1,0};players[1].invuln_timer=0;
  assert(fp_reserve(122));fpCount=122;for(size_t n=0;n<fpCount;n++)fpItems[n]=(FfItem){.id=fpNextId++,.type=n==1?PICKUP_RAILGUN:PICKUP_GLOB,.owner=-1,.slot=-1,.pos=n<2?players[1].pos:(Vector3){25,1,25}};++fpRevision;
  players[15].pos=(Vector3){-10,1,-10};players[15].respawn_timer=0;players[15].matter=100;players[15].enemyType=ENEMY_TYPE_STANDARD;combatantLife[15]=LIFE_ALIVE;fpCarrierType[15]=PICKUP_RAILGUN;equipment[15].gun[0]=EQ_RAIL;equipment[15].charge=.4f;equipment[15].held=true;
  initialSession=netTransport.session_id;joined=true;
 }
 if(t->role==NET_ROLE_CLIENT){
  
  sawItems|=fpCount==122;sawGun|=equipment[1].gun[equipment[1].selected]==EQ_RAIL;sawCharge|=equipment[1].charge>.3f;
  sawCarrier|=fpCarrierType[15]==PICKUP_RAILGUN&&players[15].respawn_timer==0&&equipment[15].charge>.3f;
  sawRestart|=firefightWaveState==FIREFIGHT_STATE_RESTART;
  sawLoot|=firefightWaveState==FIREFIGHT_STATE_LOOT&&firefightIntermissionTimer>0&&firefightIntermissionTimer<=FF_LOOT_SECONDS;
  if(!initialSession&&netWorldReady)initialSession=t->session_id;
  sawRestore|=sawRestart&&netWorldReady&&t->session_id!=initialSession&&firefightWave==6&&fpCount==122&&eq_gun(1)==EQ_RAIL;
 }
}
int main(int argc,char**argv){bool late=argc>1&&!strcmp(argv[1],"late-host");bool host=late||(argc>1&&!strcmp(argv[1],"host"));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Firefight progression LAN");
 gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;firefightCheckpoints=true;multiplayerActivePlayers=1;ResetGame();assert(net_transport_global_init());
 if(host){assert(net_transport_host(&netTransport,27942));net_set_host_local_players(1);netLobbyStarted=late;gameState=late?GAME_STATE_PLAYING:GAME_STATE_LOBBY;if(late)activePlayers=MAX_COMBATANTS;}
 else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27942));}
 double start=GetTime(),at=-1;uint32_t seq=0;bool replaced=false,captured=false,failed=false,restored=false,vote=false;uint32_t original=0;
 while(GetTime()-start<30){net_transport_pump(&netTransport,receive,net_on_connect,NULL);double now=GetTime();
  if((host&&joined)||(!host&&netWorldReady)){if(at<0)at=now;float age=now-at;
   if(host){
    eq_tick(1.f/60);net_host_tick();
    if(age>3&&eq_gun(1)==EQ_RAIL&&!captured){fp_capture();savedCount=fpCount;assert(savedCount==122);captured=true;original=netTransport.session_id;}
    if(captured&&age>5.5f&&equipment[1].charge>=eqRailCharge&&!failed){assert(fp_reserve(fpCount+1));fpItems[fpCount++]=(FfItem){.id=fpNextId++,.type=PICKUP_ROCKET,.owner=-1,.slot=-1,.pos={0,1,0}};++fpRevision;firefightScore=999;firefightTeamLives=0;for(int i=0;i<2;i++){combatantLife[i]=LIFE_SPECTATING;players[i].respawn_timer=9999;}ff_update_coop(.01f);assert(fpRestartTimer>0);failed=true;}
    if(failed&&netTransport.session_id!=original&&!restored){assert(fpCount==savedCount&&firefightScore!=999&&eq_gun(1)==EQ_RAIL);restored=true;
     // A previous-session interaction must not enter authoritative input state.
     uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,original,0);net_write_u8(&w,1);NetInputCommand bad={.sequence=99999,.pressed=NET_INPUT_EQUIP};net_write_input(&w,&bad);net_on_receive(&netTransport,1,FPS_NET_CHANNEL_INPUT,packet,w.length,NULL);assert(!netHasPendingInput[1]);sawRejected=true;
    }
    if(age>9&&restored&&!vote){ff_wave_complete();vote=true;}
    if(age>11&&firefightWaveState==FIREFIGHT_STATE_LOOT){assert(!prepReady[1]);update_firefight_logic(FF_LOOT_SECONDS);}
    if(vote&&prepReady[1])sawReady=true;
   }else{
    NetInputCommand c={.sequence=++seq};if(age>1&&eq_gun(1)==EQ_GLOB&&!replaced){c.pressed|=NET_INPUT_EQUIP;replaced=true;}
    if(age>3.2&&!sawRestart)c.held|=NET_INPUT_FIRE;
    if(firefightWaveState==FIREFIGHT_STATE_LOOT)c.pressed|=NET_INPUT_PREP_RIGHT|NET_INPUT_PREP_READY;
    if(firefightWaveState==FIREFIGHT_STATE_INTERMISSION&&!vote){c.pressed|=NET_INPUT_PREP_RIGHT|NET_INPUT_PREP_READY;vote=true;}
    uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);net_write_input(&w,&c);net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,c.pressed!=0);
   }
   if(host?sawReady:sawRestore&&vote){if(age>12)break;}
  }
  WaitTime(1.f/60);
 }
 if(host){fprintf(stderr,"host joined%d capture%d failed%d restored%d rejected%d ready%d\n",joined,captured,failed,restored,sawRejected,sawReady);assert(joined&&captured&&failed&&restored&&sawRejected&&sawReady);puts("Firefight LAN host passed");}
 else {fprintf(stderr,"client items%d gun%d charge%d carrier%d restart%d restore%d loot%d\n",sawItems,sawGun,sawCharge,sawCarrier,sawRestart,sawRestore,sawLoot);assert(sawItems&&sawGun&&sawCharge&&sawCarrier&&sawRestart&&sawRestore&&sawLoot&&firefightCheckpoints);puts("Firefight LAN client passed");}
 net_transport_shutdown(&netTransport);net_transport_global_shutdown();fp_shutdown();CloseWindow();return 0;
}
