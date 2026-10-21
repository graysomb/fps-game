#define TETHER_GRAPPLE_LAN
#include "tether_grapple_test.c"
static bool joined,sawPull,sawRelease,predictedPull;
static void receive(NetTransport *t,int peer,uint8_t channel,const uint8_t *data,size_t length,void *user){
    NetReader reader;NetPacketHeader header;net_reader_init(&reader,data,length);
    if(!net_read_header(&reader,&header))return;
    net_on_receive(t,peer,channel,data,length,user);
    if(t->role==NET_ROLE_HOST && header.type==NET_MSG_HELLO && netPeerPlayerMask[peer]){
        players[1].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[1].vel=(Vector3){0};
        players[1].respawn_timer=0;players[1].onGround=true;joined=true;
    }
    if(t->role==NET_ROLE_CLIENT){
        sawPull|=players[1].netTetherVisualActive && players[1].pos.z<-.5f;
        sawRelease|=sawPull && !players[1].netTetherVisualActive;
    }
}
int main(int argc,char **argv){
    bool host=argc>1&&!strcmp(argv[1],"host");SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Grapple LAN");
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;ResetGame();fixture();
    int v=addVoxel(0,BASE_EYE_HEIGHT,-6,true,true,GRAY,0);assert(v>=0);init_static_hash();rebuild_voxel_hash();
    players[0].pos.x=20;players[1].pos=(Vector3){0,BASE_EYE_HEIGHT,0};
    assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27941));net_set_host_local_players(1);netLobbyStarted=true;}
    else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27941));}
    double start=GetTime(),begin=-1,next=start;uint32_t sequence=0;bool hostPulled=false,hostReleased=false;
    while(GetTime()-start<6){
        net_transport_pump(&netTransport,receive,net_on_connect,NULL);double now=GetTime();
        if((host&&joined)||(!host&&netWorldReady)){
            if(begin<0)begin=now;float age=now-begin;
            if(host){net_host_tick();hostPulled|=players[1].tetherHolding && players[1].pos.z<-.5f;
                hostReleased|=hostPulled && !players[1].tetherHolding;
                if(age>3.5f)break;
            }else if(now>=next){
                next=now+1./60;NetInputCommand command={.sequence=++sequence,.mouse_look=true};
                if(age<1.1f)command.held=NET_INPUT_TETHER;
                if(age>=1.1f&&age<1.15f)command.pressed=NET_INPUT_TETHER_CANCEL;
                uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
                net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);net_write_input(&w,&command);
                net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,false);
                net_apply_movement_command(&players[1],&command,1.f/60);
                Vector3 before=players[1].vel;tether_grapple_pull(1,1.f/60);
                predictedPull|=players[1].vel.z<before.z-.001f;
                if(!command.held)assert(!tether_grapple_active(1));
                if(age>3)break;
            }
        }
        WaitTime(.001);
    }
    if(host){assert(joined && hostPulled && hostReleased);puts("Grapple LAN host passed");}
    else{fprintf(stderr,"pull%d release%d prediction%d\n",sawPull,sawRelease,predictedPull);assert(sawPull&&sawRelease&&predictedPull);puts("Grapple LAN client passed");}
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();return 0;
}
