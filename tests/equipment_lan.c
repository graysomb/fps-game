#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static bool joinedSetup,sawRocket,sawCharge,sawUpgrade,sawSwitch,sawReset,sawRail,sawSticky;
static void receive(NetTransport*t,int peer,uint8_t channel,const uint8_t*data,size_t length,void*u){
 NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);if(!net_read_header(&r,&h))return;net_on_receive(t,peer,channel,data,length,u);
 if(t->role==NET_ROLE_HOST&&h.type==NET_MSG_HELLO&&netPeerPlayerMask[peer]){players[1].pos=(Vector3){0,2,0};players[1].respawn_timer=0;players[1].matter=100;players[1].isExposed=false;players[1].invuln_timer=0;eq_reset(1);equipment[1].gun[1]=EQ_ROCKET;equipment[1].upgrades=7;pickups[0]=(Pickup){.active=true,.type=PICKUP_RAILGUN,.pos={0,2,0}};grenades[0]=(GrenadeState){.id=77,.owner=1,.pos={3,.15f,0},.fuse=1,.adhesive=true,.attached=true};joinedSetup=true;}
 if(t->role==NET_ROLE_CLIENT){sawSticky|=grenades[0].adhesive&&grenades[0].attached;sawUpgrade|=equipment[1].upgrades==7;sawSwitch|=eq_gun(1)==EQ_ROCKET;sawCharge|=equipment[1].charge>.3f;for(int n=0;n<32;n++){sawRocket|=eqRockets[n].active;sawRail|=eqTraces[n].life>0&&eqTraces[n].gun==EQ_RAIL;}sawReset|=sawCharge&&equipment[1].gun[1]==EQ_EMPTY&&equipment[1].upgrades==0;}
}
int main(int argc,char**argv){bool host=argc>1&&!strcmp(argv[1],"host");SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Equipment LAN");gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();gameState=GAME_STATE_PLAYING;assert(net_transport_global_init());if(host){assert(net_transport_host(&netTransport,27938));net_set_host_local_players(1);netLobbyStarted=true;players[0].pos=(Vector3){20,2,20};}else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27938));}
 double start=GetTime(),joined=-1,prev=start;uint32_t seq=0;bool switched=false,fired=false,equipped=false,railGiven=false,reset=false;
 while(GetTime()-start<8){net_transport_pump(&netTransport,receive,net_on_connect,NULL);double now=GetTime();float dt=now-prev;prev=now;if((host&&joinedSetup)||(!host&&netWorldReady)){if(joined<0)joined=now;float age=now-joined;
 if(host){eq_tick(dt);net_host_tick();if(age>2&&!railGiven){assert(eq_gun(1)==EQ_RAIL && !pickups[0].active && pickups[0].respawnTimer==60);railGiven=true;}if(age>4.5&&!reset){eq_reset(1);reset=true;}}
 else{NetInputCommand c={.sequence=++seq,.mouse_look=true,.aim_yaw=0,.aim_pitch=20};if(age>.3&&!switched){c.pressed=NET_INPUT_SWITCH;switched=true;}if(age>.9&&!fired){c.pressed=NET_INPUT_FIRE;fired=true;}if(age>1.4&&!equipped){c.pressed=NET_INPUT_EQUIP;equipped=true;}if(age>2.3&&age<3.4)c.held=NET_INPUT_FIRE;uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);net_write_input(&w,&c);net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,false);if(age>6)break;}}
 WaitTime(1.0/60);}
 if(host){assert(joinedSetup&&railGiven&&reset);puts("Equipment LAN host passed");}else{fprintf(stderr,"rocket%d charge%d upgrade%d switch%d reset%d rail%d\n",sawRocket,sawCharge,sawUpgrade,sawSwitch,sawReset,sawRail);assert(sawRocket&&sawCharge&&sawUpgrade&&sawSwitch&&sawReset&&sawRail&&sawSticky);puts("Equipment LAN client passed");}net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();return 0;}
