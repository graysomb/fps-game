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
        enemy_lab_spawn(4,ENEMY_TYPE_BUILDER,(Vector3){-1,1,6.5f});
        enemyStructures[0]=(EnemyStructure){.tower=true,.owner=4,.generation=enemySpecial[4].generation,.origin={0},.quarter=0,.occupant=-1,.tag=towerNextTag++,.construction=1};joined=true;
    }
    if(transport->role==NET_ROLE_CLIENT){
        EnemyStructure *t=&enemyStructures[0];preview|=t->tower&&!t->built&&t->construction>0;solid|=t->tower&&t->built&&voxel_count>=200;pulse|=t->tower&&t->pulse>0;
        removed|=solid&&!t->tower&&voxel_count==0;
    }
}
int main(int argc,char **argv){
    bool host=argc>1&&!strcmp(argv[1],"host");SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Tower LAN regression");
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=false;firefightPractice=false;currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=1;ResetGame();assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27947));net_set_host_local_players(1);netLobbyStarted=false;gameState=GAME_STATE_LOBBY;}
    else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27947));}
    double start=GetTime(),at=-1;bool built=false,repaired=false,cleaned=false;
    while(GetTime()-start<20){
        net_transport_pump(&netTransport,receive,net_on_connect,NULL);
        if(host&&joined){
            if(at<0)at=GetTime();double age=GetTime()-at;
            if(age>1&&!built){assert(tower_place(&enemyStructures[0]));built=true;}
            if(age>2.5&&!repaired){int id=enemy_voxel(enemyStructures[0].ids[0]);assert(id>=0);remove_voxel_index(id);assert(tower_place(&enemyStructures[0]));repaired=true;}
            if(age>4&&!cleaned){enemy_cleanup(true);cleaned=true;}
            net_host_tick();if(age>6)break;
        }
        if(!host&&removed)break;
        WaitTime(1.f/60);
    }
    if(host){assert(joined&&built&&repaired&&cleaned);puts("Tower LAN host passed");}
    else {fprintf(stderr,"Tower LAN preview=%d solid=%d repair=%d cleanup=%d\n",preview,solid,pulse,removed);assert(preview&&solid&&pulse&&removed);puts("Tower LAN client passed");}
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();enemy_cleanup(true);clear_world_voxels();CloseWindow();return 0;
}
