/* Real ENet host/client regression. Run host first, then client on loopback. */
#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static int explosions,deaths,snapshots;
static bool sawLateGrenade,sawRegeneration,sawSeamCrossing,sawHologram,sawDamageCue,sawRevive,sawCharge;
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
    if(t->role==NET_ROLE_CLIENT && holograms[1].life>0 && holograms[1].health==50)sawHologram=true;
    if(t->role==NET_ROLE_CLIENT){
        for(int n=0;n<8;n++)if(damageCues[1][n].until>0)sawDamageCue=true;
        if(reviveProgress[0]==1.5f && reviveOwner[0]==1)sawRevive=true;
        if(enemyCharges[40].active && enemyCharges[40].owner==4)sawCharge=true;
    }
    if(t->role==NET_ROLE_CLIENT && regeneration.active)sawRegeneration=true;
    if(t->role==NET_ROLE_CLIENT && players[0].pos.x<-30)sawSeamCrossing=true;
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id && grenades[n].owner==0)sawLateGrenade=true;
}
static void send_throw(void) {
    uint8_t data[128];NetWriter w;net_writer_init(&w,data,sizeof(data));
    net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&w,netLocalPlayerSlots[0]);
    NetInputCommand command={.sequence=10,.pressed=NET_INPUT_GRENADE|NET_INPUT_HOLOGRAM};net_write_input(&w,&command);
    assert(net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,data,w.length,true));
}
int main(int argc,char **argv) {
    bool host=argc>1 && !strcmp(argv[1],"host");
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"LAN clarity test");
    gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();
    netRequestedCreative=false;gameState=GAME_STATE_PLAYING;worldFloorHalfSize=32;worldPeriodic=host;
    assert(net_transport_global_init());
    if(host){assert(net_transport_host(&netTransport,27937));net_set_host_local_players(1);netLobbyStarted=true;players[0].pos=(Vector3){-20,BASE_EYE_HEIGHT,-20};players[0].respawn_timer=0;assert(clarity_throw(0));feedback_spawn_body((PendingRemains){.pos={0,1,0},.edge=1.25f,.color=RED});}
    else {netRequestedLocalPlayers=1;net_prepare_client_local_players(1);assert(net_transport_connect(&netTransport,"127.0.0.1",27937));}
    assert(voxel_ca_init());
    double start=GetTime(),previous=start,joined=-1;bool sent=false,repeated=false,deathSent=false,regenQueued=false,regenResent=false,crossed=false;
    while(GetTime()-start<9) {
        net_transport_pump(&netTransport,receive_test,net_on_connect,NULL);
        double now=GetTime();float dt=(float)(now-previous);previous=now;
        if(host) {
            if(joined>=0 && now-joined>2 && !crossed) {
                players[0].pos=(Vector3){31.9f,BASE_EYE_HEIGHT,0};players[0].vel=(Vector3){8,0,0};
                net_integrate_player(&players[0],.1f,true);players[0].vel=(Vector3){0};
                assert(players[0].pos.x<-30);crossed=true;
            }
            if(netTransport.peers[1].welcomed){if(joined<0)joined=now;clarity_tick(dt);hologram_tick(dt);net_host_tick();if(!deathSent && now-joined>4){
                clarity_death(0,1);combat_ui_damage(1,0);reviveProgress[0]=1.5f;reviveOwner[0]=1;
                enemyCharges[40]=(EnemyCharge){.active=true,.owner=4,.id=12345,.pos={2,1,2},.fuse=4};
                combat_ui_snapshot(-1);deathSent=true;}}
            if(joined>=0 && now-joined>1) {
                if(!regenQueued) {
                    int idx=add_static_voxel_at_grid(8,0,0,GRAY,0);assert(idx>=0);
                    assert(recycle_queue_push(&voxels[idx]));remove_voxel_index(idx);
                    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;regenQueued=true;
                }
                recycle_dead_voxels(dt);
                if(regeneration.active && !regenResent && now-joined>1.25){net_send_world_to(1);regenResent=true;}
            }
        } else {
            if(netWorldReady){assert(worldPeriodic);if(joined<0)joined=now;clarity_tick(dt);if(!sent && grenadeInventory[netLocalPlayerSlots[0]]==2){send_throw();send_throw();sent=true;}if(sent && !repeated && now-joined>.8){send_throw();repeated=true;}if(now-joined>6)break;}
        }
        voxel_ca_update(now);voxel_render_frame++;
        Camera3D camera={{3,4,8},{0,1,0},{0,1,0},60,CAMERA_PERSPECTIVE};
        BeginDrawing();ClearBackground(BLACK);BeginMode3D(camera);DrawVoxels(camera);EndMode3D();EndDrawing();
        assert(voxelCaReady);
        WaitTime(1.0/60);
    }
    if(host){assert(joined>=0 && deathSent && regenResent && recycleQueueCount==0);assert(grenadeInventory[1]==1);assert(holograms[1].generation==1);for(int n=0;n<32;n++)assert(!grenades[n].id);puts("LAN host: one throw despite duplicates; ordnance cleaned up");}
    else {assert(netWorldReady && sawLateGrenade && sent && repeated && sawHologram && sawDamageCue && sawRevive && sawCharge);assert(explosions==2 && deaths==1);assert(netVoxelProxyCount==1 && voxel_count==1 && table_get_static_only(8,0,0)>=0 && sawRegeneration && !regeneration.active);assert(fabsf(v_length(netVoxelProxies[0].x_axis)-1.25f)<.01f);assert(grenadeInventory[netLocalPlayerSlots[0]]==1);for(int n=0;n<32;n++)assert(!grenades[n].id);puts("LAN client: late join, inventory, snapshot loss and unique events passed");}
    regeneration_shutdown();
    assert(host?crossed:sawSeamCrossing);
    net_transport_shutdown(&netTransport);net_transport_global_shutdown();shutdown_world_visuals();CloseWindow();return 0;
}
