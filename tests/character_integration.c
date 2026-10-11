#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

static void shot(RenderTexture2D rt,const char *name) {
    char path[256];snprintf(path,sizeof(path),"artifacts/rat/%s.png",name);
    Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);assert(ExportImage(image,path));UnloadImage(image);
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Rat gameplay integration");
    assert(IsWindowReady());SetTargetFPS(0);
    init_static_hash();clear_world_voxels();clear_pickups();
    weaponRenderingReady=weapon_renderer_init();assert(character_renderer_init());
    gameMode=GAME_MODE_DEATHMATCH;activePlayers=4;
    for(int i=0;i<4;i++) {
        players[i]=(Player){.pos={i==0?0:(i-2)*.8f,1,i==0?0:-2},.matter=100,.matterMax=100,
            .onGround=true,.yaw=i==0?0:180,.tetherVoxel=-1,.last_melee_time=-1000};
        playerInput[i]=INPUT_TYPE_KEYBOARD;weapon_visual_reset(&weaponVisuals[i]);
        assert(player_uses_rat(i));
    }
    playerInput[3]=INPUT_TYPE_BOT_EASY;assert(!player_uses_rat(3));
    assert(!player_uses_rat(4));
    for(int role=NET_ROLE_HOST;role<=NET_ROLE_CLIENT;role++) {
        netTransport.role=role;netPlayerPresent[3]=true;assert(player_uses_rat(3));
        netPlayerPresent[3]=false;assert(!player_uses_rat(3));
    }
    netTransport.role=NET_ROLE_OFFLINE;
    RenderTexture2D screens[4]={0};int count=0,w=0,h=0;
    for(int n=1;n<=4;n++) {
        activePlayers=n;render_gameplay_view(screens,&count,&w,&h,false);
        char name[64];snprintf(name,sizeof(name),"gameplay-%d-view",n);shot(screens[0],name);
    }
    activePlayers=4;
    double started=GetTime();
    for(int frame=0;frame<60;frame++)render_gameplay_view(screens,&count,&w,&h,false);
    printf("Rat four-view render average %.2f ms\n",(GetTime()-started)*1000/60);
    float matter=players[0].matter;Vector3 position=players[0].pos;
    players[0].dynamicShotActive=true;players[0].meleeSwingActive=true;
    players[0].meleeSwingStartTime=(float)GetTime()-.4f*MELEE_ANIM_DURATION_SECONDS;
    render_gameplay_view(screens,&count,&w,&h,false);shot(screens[0],"gameplay-launcher-melee");
    assert(players[0].matter==matter && Vector3Distance(players[0].pos,position)==0);
    players[0].meleeSwingActive=false;
    players[1].isExposed=true;players[2].invuln_timer=2;
    render_gameplay_view(screens,&count,&w,&h,false);shot(screens[0],"gameplay-status");
    gameMode=GAME_MODE_FIREFIGHT;firefightHumans=3;
    combatantLife[1]=LIFE_DOWNED;players[1].death_pos=players[1].pos;players[1].respawn_timer=10;
    render_gameplay_view(screens,&count,&w,&h,false);assert(characterSamples[1].downed);shot(screens[0],"gameplay-downed");
    gameMode=GAME_MODE_DEATHMATCH;players[1].respawn_timer=0;combatantLife[1]=LIFE_ALIVE;
    weaponBloomEnabled=false;render_gameplay_view(screens,&count,&w,&h,false);shot(screens[0],"gameplay-no-bloom");
    // Replay an actual host snapshot into client state; remote keyboard settings are irrelevant.
    uint8_t packet[2048];NetWriter writer;net_writer_init(&writer,packet,sizeof(packet));
    netTransport.role=NET_ROLE_HOST;netTransport.session_id=77;
    net_write_header(&writer,NET_MSG_PLAYER_STATE,0,77,200);net_write_u8(&writer,1);net_write_player_state(&writer,2);
    netTransport.role=NET_ROLE_CLIENT;playerInput[2]=INPUT_TYPE_BOT_HARD;netPlayerPresent[2]=false;
    NetTransport host={0};
    assert(net_transport_host(&host,28891));
    assert(net_transport_connect(&netTransport,"127.0.0.1",28891));
    netTransport.session_id=77;
    double deadline=GetTime()+5;
    while((!netTransport.connected || !host.peers[1].connected) && GetTime()<deadline) {
        net_transport_pump(&host,NULL,NULL,NULL);net_transport_pump(&netTransport,NULL,NULL,NULL);WaitTime(.001);
    }
    assert(netTransport.connected && host.peers[1].connected);
    assert(net_transport_send(&host,1,FPS_NET_CHANNEL_SNAPSHOT,packet,writer.length,true));
    while(!netPlayerPresent[2] && GetTime()<deadline) {
        net_transport_pump(&host,NULL,NULL,NULL);net_transport_pump(&netTransport,net_on_receive,NULL,NULL);WaitTime(.001);
    }
    assert(netPlayerPresent[2] && player_uses_rat(2));
    netLocalPlayerCount=1;netLocalPlayerSlots[0]=0;
    render_gameplay_view(screens,&count,&w,&h,false);shot(screens[0],"gameplay-client");
    for(int i=0;i<count;i++)UnloadRenderTexture(screens[i]);
    for(int i=0;i<4;i++)weapon_bloom_unload(&weaponBlooms[i]);
    net_transport_shutdown(&netTransport);net_transport_shutdown(&host);
    character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();
    puts("Rat integration passed");return 0;
}
