#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static bool joined,preview,solid,pulse,removed;
static void receive(NetTransport *transport,int peer,uint8_t channel,const uint8_t *data,size_t length,void *user){
    NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);if(!net_read_header(&r,&h))return;
    net_on_receive(transport,peer,channel,data,length,user);
    if(transport->role==NET_ROLE_HOST && h.type==NET_MSG_HELLO && netPeerPlayerMask[peer]){
        net_host_start_match();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();
        for(int i=4;i<MAX_COMBATANTS;i++){players[i].respawn_timer=9999;combatantLife[i]=LIFE_UNUSED;}
        players[0].pos=(Vector3){-20,1,-20};players[1].pos=(Vector3){-15,1,-20};
        assert(golem_physics_start_at((Vector3){0,.04f,10},0));go_combat_start();joined=true;
    }
    if(transport->role==NET_ROLE_CLIENT){
        preview|=goRemoteActive&&!goRemoteDead;solid|=goRemoteActive&&netVoxelProxyCount>100;pulse|=goRemoteActive&&goRemoteShield<150;removed|=preview&&!goRemoteActive;
    }
}
int main(int argc,char **argv){
    bool host=argc>1&&!strcmp(argv[1],"host");SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Golem LAN regression");
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=false;firefightPractice=false;currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=1;ResetGame();assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27948));net_set_host_local_players(1);netLobbyStarted=false;gameState=GAME_STATE_LOBBY;}
    else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27948));}
    double start=GetTime(),at=-1;bool built=false,repaired=false,cleaned=false;
    while(GetTime()-start<20){
        net_transport_pump(&netTransport,receive,net_on_connect,NULL);
        if(host&&joined){
            if(at<0)at=GetTime();double age=GetTime()-at;
            if(age>1&&!built){go_combat_damage_from(25,true,0);built=true;}
            if(age>2.5&&!repaired){go_combat_damage_from(1000,false,0);repaired=true;}
            if(age>4&&!cleaned){enemy_cleanup(true);cleaned=true;}
            net_host_tick();if(age>6)break;
        }
        if(!host&&removed)break;
        WaitTime(1.f/60);
    }
    if(host){assert(joined&&built&&repaired&&cleaned);puts("Golem LAN host passed");}
    else {fprintf(stderr,"Golem LAN preview=%d solid=%d repair=%d cleanup=%d\n",preview,solid,pulse,removed);assert(preview&&solid&&pulse&&removed);puts("Golem LAN client passed");}
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();enemy_cleanup(true);clear_world_voxels();CloseWindow();return 0;
}
