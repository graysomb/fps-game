#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static bool joined,sawSpin[2],sawBeam[2];
static void receive(NetTransport *t,int peer,uint8_t channel,const uint8_t *data,size_t length,void *user){
    NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);
    assert(net_read_header(&r,&h));net_on_receive(t,peer,channel,data,length,user);
    if(t->role==NET_ROLE_HOST && h.type==NET_MSG_HELLO && netPeerPlayerMask[peer] && !joined){
        players[1].pos=(Vector3){0,BASE_EYE_HEIGHT,0};
        // Below the Versus kill-based capacity caps, so only charging changes this value.
        players[1].matter=45;
        players[1].respawn_timer=players[1].invuln_timer=0;eq_reset(1);equipment[1].gun[0]=EQ_RAIL;joined=true;
    }
    if(t->role!=NET_ROLE_CLIENT)return;
    int trial=players[1].kills;
    if(trial<2 && rail_spinning(1) && fabsf(players[1].yaw)>45)sawSpin[trial]=true;
    if(trial>0 && trial<=2)for(int n=0;n<32;n++)if(eqTraces[n].life>0 && eqTraces[n].owner==1 && eqTraces[n].gun==EQ_RAIL)sawBeam[trial-1]=true;
    assert(!rail_spinning(0));
}
int main(int argc,char **argv){
    bool host=argc==2&&!strcmp(argv[1],"host");
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"360 rail LAN");SetLoggingEnabled(false);
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;ResetGame();clear_world_voxels();clear_pickups();init_static_hash();
    gameState=GAME_STATE_PLAYING;assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27958));net_set_host_local_players(1);netLobbyStarted=true;
        players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,-8};players[0].yaw=42;players[0].invuln_timer=0;
    }else{netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27958));}
    double start=GetTime(),begin=-1,previous=start;uint32_t sequence=0;bool revived=false;
    while(GetTime()-start<10){
        net_transport_pump(&netTransport,receive,net_on_connect,NULL);
        double now=GetTime();float dt=now-previous;previous=now;
        if((host&&joined)||(!host&&netWorldReady)){
            if(begin<0)begin=now;
            float age=now-begin;
            if(host){
                if(players[1].kills==1 && !revived){
                    players[0].respawn_timer=0;players[0].invuln_timer=0;players[0].matter=100;players[0].isExposed=false;
                    players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,-8};revived=true;
                }
                if(revived&&rail_spinning(1))players[0].pos.x+=dt*4;
                eq_tick(dt);net_host_tick();
                assert(players[0].yaw==42);
            }else{
                eq_tick(dt);
                NetInputCommand c={.sequence=++sequence,.mouse_look=1};
                if((age>.25f&&age<1.25f)||(age>3.6f&&age<4.6f))c.held=NET_INPUT_FIRE;
                if(rail_spinning(1))c.held=NET_INPUT_RAIL_SPIN;
                uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
                net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,1);net_write_input(&w,&c);
                net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,false);
            }
            if(age>5.8f)break;
        }
        WaitTime(1.f/120);
    }
    fprintf(stderr,"role=%d kills=%d matter=%.2f spins=%d,%d beams=%d,%d\n",netTransport.role,players[1].kills,players[1].matter,sawSpin[0],sawSpin[1],sawBeam[0],sawBeam[1]);
    assert(players[1].kills==2 && fabsf(players[1].matter-15)<.01f);
    if(host)assert(players[0].death_pos.x>1);
    else assert(sawSpin[0]&&sawSpin[1]&&sawBeam[0]&&sawBeam[1]);
    puts(host?"Host confirmed stationary and moving hits":"Client observed both spins, hits, and beams");
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();
}
