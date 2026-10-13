/* Real two-process Crosscut decoration replication, without new message types. */
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
int main(int argc,char **argv) {
    bool host=argc>1 && !strcmp(argv[1],"host");
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Setting LAN regression");
    gameMode=GAME_MODE_DEATHMATCH;currentWorldType=host?WORLD_TYPE_CROSSCUT:WORLD_TYPE_TEST;
    multiplayerActivePlayers=2;ResetGame();netRequestedCreative=false;gameState=GAME_STATE_PLAYING;
    assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27938));net_set_host_local_players(1);netLobbyStarted=true;}
    else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27938));}
    int chosen=-1,phase=0;bool sawHidden=false,sawRestored=false;Color color=GRAY;int cell[3]={0};
    double started=GetTime(),joined=-1;
    while(GetTime()-started<8) {
        net_transport_pump(&netTransport,net_on_receive,net_on_connect,NULL);
        double now=GetTime();
        if(host && netTransport.peers[1].welcomed) {
            if(joined<0)joined=now;
            if(now-joined>1 && phase==0) {
                setting_update();for(int d=0;d<settingDetailCount;d++)if(settingDetails[d].visible){chosen=d;break;}assert(chosen>=0);
                memcpy(cell,settingDetails[chosen].cells[2],sizeof(cell));int idx=table_get_static_only(cell[0],cell[1],cell[2]);assert(idx>=0);color=voxels[idx].color;
                remove_voxel_index(idx);net_send_world_to(1);phase=1;
            }
            if(now-joined>3 && phase==1){add_static_voxel_at_grid(cell[0],cell[1],cell[2],color,0);net_send_world_to(1);phase=2;}
        }
        if(!host && netWorldReady) {
            assert(settingMap==2);setting_update();
            if(chosen<0){for(int d=0;d<settingDetailCount;d++)if(settingDetails[d].visible){chosen=d;break;}assert(chosen>=0);}
            if(!settingDetails[chosen].visible)sawHidden=true;
            if(sawHidden && settingDetails[chosen].visible){sawRestored=true;break;}
        }
        WaitTime(.01);
    }
    if(host)assert(phase==2);else assert(sawHidden && sawRestored);
    puts(host?"Setting LAN host passed":"Setting LAN client: identity, destruction, restoration passed");
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();CloseWindow();return 0;
}
