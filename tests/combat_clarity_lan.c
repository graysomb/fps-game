/* Real ENet host/client regression. Run host first, then client on loopback. */
#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static int explosions,deaths,snapshots;
static bool sawLateGrenade;
static void receive_test(NetTransport *t,int peer,uint8_t channel,const uint8_t *data,size_t length,void *user) {
    NetReader r;NetPacketHeader h;net_reader_init(&r,data,length);
    if(!net_read_header(&r,&h))return;
    /* Drop whole application snapshots; reorder/partial cases have unit coverage. */
    if(t->role==NET_ROLE_CLIENT && h.type==NET_MSG_GRENADES && ++snapshots%7==0)return;
    uint32_t before=clarityLastEvent;
    net_on_receive(t,peer,channel,data,length,user);
    if(t->role==NET_ROLE_HOST && h.type==NET_MSG_HELLO && netPeerPlayerMask[peer]) {
        for(int i=1;i<4;i++)if(netPeerPlayerMask[peer]&(1u<<i)){grenadeInventory[i]=2;players[i].pos=(Vector3){20,BASE_EYE_HEIGHT,20};}
    }
    if(t->role==NET_ROLE_CLIENT && h.type==NET_MSG_GAME_EVENT) {
        if(clarityLastEvent>before){(void)net_read_u32(&r);int kind=net_read_u8(&r);if(kind==2)explosions++;else deaths++;}
        net_on_receive(t,peer,channel,data,length,user); /* Duplicate delivery must be harmless. */
    }
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id && grenades[n].owner==0)sawLateGrenade=true;
}
static void send_throw(void) {
    uint8_t data[128];NetWriter w;net_writer_init(&w,data,sizeof(data));
    net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);
    NetInputCommand command={.sequence=10,.pressed=NET_INPUT_GRENADE};net_write_input(&w,&command);
    assert(net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,data,w.length,true));
}
int main(int argc,char **argv) {
    bool host=argc>1 && !strcmp(argv[1],"host");
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"LAN clarity test");
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();
    netRequestedCreative=false;gameState=GAME_STATE_PLAYING;worldFloorHalfSize=32;
    assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27937));net_set_host_local_players(1);netLobbyStarted=true;players[0].pos=(Vector3){-20,BASE_EYE_HEIGHT,-20};players[0].respawn_timer=0;assert(clarity_throw(0));feedback_spawn_body((PendingRemains){.pos={0,1,0},.edge=1.25f,.color=RED});}
    else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27937));}
    double start=GetTime(),previous=start,joined=-1;bool sent=false,repeated=false,deathSent=false;
    while(GetTime()-start<9) {
        net_transport_pump(&netTransport,receive_test,net_on_connect,NULL);
        double now=GetTime();float dt=(float)(now-previous);previous=now;
        if(host) {
            if(netTransport.peers[1].welcomed){if(joined<0)joined=now;clarity_tick(dt);net_host_tick();if(!deathSent && now-joined>4){clarity_death(0,1);deathSent=true;}}
        } else {
            if(netWorldReady){if(joined<0)joined=now;clarity_tick(dt);if(!sent && grenadeInventory[netLocalPlayerSlots[0]]==2){send_throw();send_throw();sent=true;}if(sent && !repeated && now-joined>.8){send_throw();repeated=true;}if(now-joined>6)break;}
        }
        WaitTime(1.0/60);
    }
    if(host){assert(joined>=0 && deathSent);assert(grenadeInventory[1]==1);for(int n=0;n<32;n++)assert(!grenades[n].id);puts("LAN host: one throw despite duplicates; ordnance cleaned up");}
    else {assert(netWorldReady && sawLateGrenade && sent && repeated);assert(explosions==2 && deaths==1);assert(netVoxelProxyCount==1 && voxel_count==0);assert(fabsf(v_length(netVoxelProxies[0].x_axis)-1.25f)<.01f);assert(grenadeInventory[netLocalPlayerSlots[0]]==1);for(int n=0;n<32;n++)assert(!grenades[n].id);puts("LAN client: late join, inventory, snapshot loss and unique events passed");}
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();return 0;
}
