#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
bool ffTestReviveHeld[4];
static bool joined;
static void receive(NetTransport*t,int peer,uint8_t channel,const uint8_t*data,size_t length,void*u){
 NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);if(!net_read_header(&r,&h))return;
 net_on_receive(t,peer,channel,data,length,u);
 if(t->role==NET_ROLE_HOST&&h.type==NET_MSG_HELLO&&netPeerPlayerMask[peer]){net_host_start_match();gameState=GAME_STATE_PLAYING;joined=true;}
}
int main(int argc,char**argv){
 bool host=argc>1&&!strcmp(argv[1],"host"),victory=false,watched=false,restarted=false;
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Victory LAN regression");SetExitKey(KEY_NULL);
 character_renderer_init();weapon_renderer_init();duck_renderer_init(NULL);
 gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;multiplayerActivePlayers=1;ResetGame();assert(net_transport_global_init());
 if(host){assert(net_transport_host(&netTransport,27943));net_set_host_local_players(1);netLobbyStarted=false;gameState=GAME_STATE_LOBBY;}
 else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27943));}
 double start=GetTime(),at=-1,last=start;
 while(GetTime()-start<35){
  net_transport_pump(&netTransport,receive,net_on_connect,NULL);double now=GetTime();float dt=(float)(now-last);last=now;
  if((host&&joined)||(!host&&netWorldReady)){
   if(at<0)at=now;
   if(host&&now-at>1&&!victory){firefightWave=FF_WAVE_COUNT;ff_wave_complete();assert(gameState==GAME_STATE_GAMEOVER);victory=true;}
   if(host)net_host_tick();
   ffend_frame(dt);
   if(ffend.phase==FFEND_RESULTS&&!watched){assert(ffend_start_movie());watched=true;}
   if(host&&now-at>8&&watched&&!restarted){ResetGame();gameState=GAME_STATE_PLAYING;fp_network_reset();restarted=true;}
   if(!host&&watched&&gameState==GAME_STATE_PLAYING&&ffend.phase==FFEND_NONE){assert(!town.ready);restarted=true;break;}
   if(host&&restarted&&now-at>22)break;
  }
  WaitTime(.016);
 }
 fprintf(stderr,"host=%d watched=%d restart=%d state=%d phase=%d ready=%d connected=%d generation=%u wave=%d waveState=%d\n",host,watched,restarted,gameState,ffend.phase,netWorldReady,netTransport.connected,fpGeneration,firefightWave,firefightWaveState);assert(watched&&restarted);printf("Victory LAN %s: playback and authoritative restart passed\n",host?"host":"client");
 ffend_reset();net_transport_shutdown(&netTransport);net_transport_global_shutdown();camp_shutdown();duck_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();return 0;
}
