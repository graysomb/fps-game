/* Compile this instead of fps_ray.c to exercise the actual gameplay renderer,
 * FireVoxel, tether resolver, split-screen allocation, and LAN snapshots. */
#define main fps_game_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

static void capture_view(RenderTexture2D screen,const char *out,const char *name) {
    char path[1024]; snprintf(path,sizeof(path),"%s/%s.png",out,name);
    Image image=LoadImageFromTexture(screen.texture); ImageFlipVertical(&image);
    assert(ExportImage(image,path)); UnloadImage(image);
}

int main(int argc,char **argv) {
    const char *out=argc>1 ? argv[1] : "artifacts/weapon/gameplay";
    SetLoggingEnabled(false); SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN); InitWindow(1280,720,"Weapon integration test");
    assert(IsWindowReady()); SetRandomSeed(1);
    clear_world_voxels(); clear_pickups(); init_static_hash();
    weaponRenderingReady=weapon_renderer_init(); assert(weaponRenderingReady);
    gameMode=GAME_MODE_DEATHMATCH; activePlayers=4;
    for (int i=0;i<4;++i) {
        players[i]=(Player){.pos={i==0 ? 0 : (i-2)*1.5f,2.0f,i==0 ? 0 : -4},.matter=100,.matterMax=100,
                            .last_shot_time=-1000,.tetherVoxel=-1,.last_melee_time=-1000};
        weapon_visual_reset(&weaponVisuals[i]); playerInput[i]=INPUT_TYPE_KEYBOARD;
    }
    // Feedback occurs only after accepted shots, never on cooldown/resource failure.
    FireVoxel(0); assert(weaponVisuals[0].shot_sequence==1);
    int bullet_count=voxel_count; float matter=players[0].matter;
    FireVoxel(0); assert(weaponVisuals[0].shot_sequence==1 && voxel_count==bullet_count);
    assert(players[0].matter==matter);
    players[0].last_shot_time=-1000; players[0].matter=0;
    FireVoxel(0); assert(weaponVisuals[0].shot_sequence==1 && voxel_count==bullet_count);
    players[0].matter=100;
    Vector3 bullet_start=v_add(players[0].pos,v_mul(player_forward(&players[0]),.8f));
    assert(v_length(v_sub(voxels[0].pos,bullet_start))<.0001f);
    assert(voxels[0].isBullet);
    // Round-trip the same player snapshot used by the game.
    uint8_t packet[1024]; NetWriter writer; NetReader reader; NetPlayerWireState state;
    net_writer_init(&writer,packet,sizeof(packet)); net_write_player_state(&writer,0);
    net_reader_init(&reader,packet,writer.length); int slot=-1;
    assert(net_read_player_state(&reader,&slot,&state)); assert(slot==0);
    assert(state.visual.shot_sequence==1 && state.visual.shot_age_ms<120);
    // Feed an authoritative packet through the production client receive path.
    net_writer_init(&writer,packet,sizeof(packet));
    net_write_header(&writer,NET_MSG_PLAYER_STATE,0,77,100);
    net_write_u8(&writer,1); net_write_player_state(&writer,0);
    netTransport.role=NET_ROLE_CLIENT; netTransport.session_id=77;
    netLocalPlayerCount=0; weapon_visual_reset(&weaponVisuals[0]);
    net_on_receive(&netTransport,0,FPS_NET_CHANNEL_SNAPSHOT,packet,writer.length,NULL);
    assert(weaponVisuals[0].shot_sequence==1);
    float received_time=weaponVisuals[0].shot_time;
    net_on_receive(&netTransport,0,FPS_NET_CHANNEL_SNAPSHOT,packet,writer.length,NULL);
    assert(weaponVisuals[0].shot_time==received_time); // duplicate never restarts recoil
    netTransport.role=NET_ROLE_OFFLINE;
    clear_world_voxels(); weaponVisuals[0].shot_time=-1000;
    // Valid host tether target; its solver hold point must remain unchanged.
    int held=addVoxel(.55f,2.2f,-1.5f,false,true,BLUE,0); assert(held>=0);
    players[0].tetherHolding=true; players[0].tetherVoxel=held;
    players[0].tetherVoxelIdentity=voxels[held].identity; players[0].tetherFace=-1;
    Vector3 hold_before=player_tether_hold_position(&players[0]);
    weaponVisuals[0].claw_open=1;
    RenderTexture2D screens[MAX_PLAYERS]={0}; int count=0,w=0,h=0;
    for (int views=1;views<=4;++views) {
        activePlayers=views;
        render_gameplay_view(screens,&count,&w,&h,false);
        char name[64]; snprintf(name,sizeof(name),"%d-player-view",views);
        capture_view(screens[0],out,name);
        assert(count==views);
        assert(w==GetScreenWidth()/(views>1 ? 2 : 1));
        assert(h==GetScreenHeight()/(views>2 ? 2 : 1));
        for (int i=0;i<views;++i) assert(weaponBlooms[i].width==w && weaponBlooms[i].height==h);
    }
    assert(v_length(v_sub(hold_before,player_tether_hold_position(&players[0])))<.0001f);
    weaponBloomEnabled=false; render_gameplay_view(screens,&count,&w,&h,false);
    for (int i=0;i<4;++i) assert(weaponBlooms[i].emission.id==0);
    weaponBloomEnabled=true;
    players[0].meleeSwingActive=true; players[0].meleeSwingStartTime=(float)GetTime()-.12f;
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee");
    players[0].meleeSwingActive=false; players[0].respawn_timer=2;
    render_gameplay_view(screens,&count,&w,&h,false);
    assert(!player_has_gun(0)); capture_view(screens[0],out,"respawn");
    players[0].respawn_timer=0;
    render_gameplay_view(screens,&count,&w,&h,true); capture_view(screens[0],out,"creative");
    assert(!drawWeaponModels); for (int i=0;i<4;++i) assert(weaponBlooms[i].emission.id==0);
    // Client tether uses its replicated target and gold state, with no local voxel ID.
    netTransport.role=NET_ROLE_CLIENT; netLocalPlayerCount=2;
    netLocalPlayerSlots[0]=0; netLocalPlayerSlots[1]=1;
    players[0].netTetherVisualActive=true; players[0].netGoldTetherVisualActive=true;
    players[0].netTetherVisualTarget=(Vector3){.55f,2.2f,-1.5f};
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"client-gold-tether");
    assert(count==2 && weaponBlooms[2].emission.id==0 && weaponBlooms[3].emission.id==0);
    Vector3 target; assert(player_tether_visual_target(0,&target));
    assert(v_length(v_sub(target,players[0].netTetherVisualTarget))<.0001f);
    netTransport.role=NET_ROLE_OFFLINE;
    for (int i=0;i<count;++i) UnloadRenderTexture(screens[i]);
    for (int i=0;i<MAX_PLAYERS;++i) weapon_bloom_unload(&weaponBlooms[i]);
    weapon_renderer_shutdown(); shutdown_shard_rendering(); shutdown_world_visuals(); clear_all_chunks();
    if (greedyMesh.vertices) UnloadMesh(greedyMesh);
    if (voxelInstanceVboId) rlUnloadVertexBuffer(voxelInstanceVboId);
    if (greedyMaterialInit) UnloadMaterial(greedyMaterial);
    CloseWindow(); puts("Weapon gameplay integration passed"); return 0;
}
