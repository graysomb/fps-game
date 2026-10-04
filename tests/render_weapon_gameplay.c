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

static WeaponPose pose_at(Player p,float progress,bool first_person,float aspect) {
    p.meleeSwingActive=true; p.meleeSwingStartTime=10-progress*MELEE_ANIM_DURATION_SECONDS;
    WeaponMeleePose melee=player_weapon_melee_sample(&p,10);
    Vector3 forward,right,up; melee_swing_basis(&p,&forward,&right,&up);
    return weapon_pose(p.pos,forward,right,up,first_person,aspect,&melee,0);
}

static void test_melee_alignment(void) {
    const float aspects[]={16.0f/9,8.0f/9};
    const float pitches[]={0,89,-89};
    for (int view=0;view<2;++view) for (int a=0;a<2;++a) for (int p=0;p<3;++p) {
        Player player=players[0]; player.pitch=pitches[p]; player.yaw=37;
        for (int phase=29;phase<=61;++phase) {
            float progress=phase/100.0f;
            float fraction=melee_reach_fraction(progress);
            Vector3 end=v_add(melee_swing_origin(&player,fraction),
                             v_mul(melee_swing_direction(&player,fraction),MELEE_RANGE*fraction));
            WeaponPose pose=pose_at(player,progress,view==0,aspects[a]);
            assert(v_isfinite(pose.muzzle) && v_isfinite(pose.grip_heel));
            assert(v_length(v_sub(pose.grip_heel,end))<.0001f);
        }
    }
}

static void set_melee_phase(float progress) {
    players[0].meleeSwingStartTime=(float)GetTime()-progress*MELEE_ANIM_DURATION_SECONDS;
}

static void capture_melee_sequence(const char *out) {
    // Deterministic 60 fps pose samples on a plain backdrop, without world/HUD noise.
    RenderTexture2D scene=LoadRenderTexture(960,540); WeaponBloom bloom={0};
    assert(weapon_bloom_resize(&bloom,960,540));
    Camera3D camera={.position=players[0].pos,.target=v_add(players[0].pos,player_forward(&players[0])),
        .up={0,1,0},.fovy=60,.projection=CAMERA_PERSPECTIVE};
    WeaponVisual visual; weapon_visual_reset(&visual);
    bool top_clipped=false;
    for (int frame=0;frame<=20;++frame) {
        WeaponPose pose=pose_at(players[0],frame/20.0f,true,16.0f/9);
        BeginTextureMode(scene); ClearBackground((Color){29,55,94,255});
        BeginMode3D(camera); weapon_draw(pose,&visual,camera.position,1,false,false,false);
        EndMode3D(); EndTextureMode();
        weapon_bloom_begin(&bloom,scene);
        BeginMode3D(camera); weapon_draw(pose,&visual,camera.position,1,false,false,true);
        EndMode3D(); weapon_bloom_composite(&bloom,scene);
        DrawLine(474,270,486,270,WHITE); DrawLine(480,264,480,276,WHITE); EndTextureMode();
        char name[64]; snprintf(name,sizeof(name),"swing-frame-%02d",frame); capture_view(scene,out,name);
        // No silver geometry should be clipped by the top during the swing.
        Image image=LoadImageFromTexture(scene.texture); Color *pixels=LoadImageColors(image);
        for (int x=0;x<image.width;++x) {
            Color c=pixels[(image.height-1)*image.width+x];
            if (c.r>80 && c.g>80 && c.b>80) {
                fprintf(stderr,"Melee top clipping at phase %.2f\n",frame/20.0f);
                top_clipped=true; break;
            }
        }
        UnloadImageColors(pixels); UnloadImage(image);
    }
    assert(!top_clipped);
    weapon_bloom_unload(&bloom); UnloadRenderTexture(scene);
}

static void test_melee_hits(void) {
    Player saved[MAX_PLAYERS]; memcpy(saved,players,sizeof(players));
    activePlayers=2;
    players[1].pos=v_add(players[0].pos,(Vector3){0,0,-1.5f});
    perform_melee(0); assert(players[0].meleeSwingActive);
    float started=players[0].meleeSwingStartTime;
    perform_melee(0); assert(players[0].meleeSwingStartTime==started);
    set_melee_phase(.15f); update_melee_swings(); assert(players[1].matter==100);
    WeaponPose miss=pose_at(players[0],.52f,true,16.0f/9);
    set_melee_phase(.52f); update_melee_swings();
    assert(players[0].meleeSwingHitApplied);
    assert(fabsf(players[1].matter-(100-MATTER_MELEE_DAMAGE))<.0001f);
    assert(fabsf(players[1].vel.z+MELEE_KNOCKBACK_SPEED)<.0001f);
    assert(fabsf(players[1].vel.y-MELEE_UPWARD_BOOST)<.0001f);
    WeaponPose hit=pose_at(players[0],.52f,true,16.0f/9);
    assert(v_length(v_sub(miss.muzzle,hit.muzzle))<.0001f);
    update_melee_swings(); assert(fabsf(players[1].matter-(100-MATTER_MELEE_DAMAGE))<.0001f);
    set_melee_phase(1.01f); update_melee_swings(); assert(!players[0].meleeSwingActive);
    // A completed swing accepts a repeat; an empty swing leaves matter untouched.
    players[1].pos=(Vector3){10,2,-10}; players[0].last_melee_time=-1000;
    perform_melee(0); set_melee_phase(.52f); update_melee_swings();
    assert(!players[0].meleeSwingHitApplied);
    players[0].meleeSwingActive=false; players[0].last_melee_time=-1000;
    // Mine a single static cube on the production sweep, outside the overlap box.
    activePlayers=1; players[0].matter=40;
    float fraction=melee_reach_fraction(.52f);
    Vector3 target=v_add(melee_swing_origin(&players[0],fraction),
                        v_mul(melee_swing_direction(&players[0],fraction),1.5f));
    int cube=addVoxel(target.x,target.y,target.z,true,false,BLUE,0); assert(cube>=0);
    rebuild_static_hash_if_dirty();
    perform_melee(0); set_melee_phase(.52f); update_melee_swings();
    assert(players[0].meleeSwingHitApplied && players[0].matter==40+MATTER_MELEE_HARVEST);
    WeaponPose mining=pose_at(players[0],.52f,true,16.0f/9);
    assert(v_length(v_sub(miss.muzzle,mining.muzzle))<.0001f);
    update_melee_swings(); assert(players[0].matter==40+MATTER_MELEE_HARVEST);
    players[0].respawn_timer=1; update_melee_swings(); assert(!players[0].meleeSwingActive);
    clear_world_voxels(); memcpy(players,saved,sizeof(players)); activePlayers=4;
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
    test_melee_alignment(); test_melee_hits(); capture_melee_sequence(out);
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
    players[0].meleeSwingActive=true; set_melee_phase(.45f);
    uint8_t packet[1024]; NetWriter writer; NetReader reader; NetPlayerWireState state;
    net_writer_init(&writer,packet,sizeof(packet)); net_write_player_state(&writer,0);
    net_reader_init(&reader,packet,writer.length); int slot=-1;
    assert(net_read_player_state(&reader,&slot,&state)); assert(slot==0);
    assert(state.visual.shot_sequence==1 && state.visual.shot_age_ms<120);
    assert(state.visual.flags & NET_PLAYER_VISUAL_MELEE);
    assert(fabsf(state.visual.melee_progress/255.0f-.45f)<.01f);
    // Feed an authoritative packet through the production client receive path.
    net_writer_init(&writer,packet,sizeof(packet));
    net_write_header(&writer,NET_MSG_PLAYER_STATE,0,77,100);
    net_write_u8(&writer,1); net_write_player_state(&writer,0);
    netTransport.role=NET_ROLE_CLIENT; netTransport.session_id=77;
    netLocalPlayerCount=0; weapon_visual_reset(&weaponVisuals[0]);
    net_on_receive(&netTransport,0,FPS_NET_CHANNEL_SNAPSHOT,packet,writer.length,NULL);
    assert(weaponVisuals[0].shot_sequence==1);
    assert(players[0].meleeSwingActive);
    assert(fabsf(melee_anim_progress(&players[0],(float)GetTime())-.45f)<.02f);
    float received_time=weaponVisuals[0].shot_time;
    net_on_receive(&netTransport,0,FPS_NET_CHANNEL_SNAPSHOT,packet,writer.length,NULL);
    assert(weaponVisuals[0].shot_time==received_time); // duplicate never restarts recoil
    netTransport.role=NET_ROLE_OFFLINE;
    players[0].meleeSwingActive=false;
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
    const float phases[]={0,.14f,.28f,.52f,.62f,.82f,1};
    const char *names[]={"idle","windup","active-start","peak","active-end","recovery","recovered"};
    // Capture the full-size view and another player's world view of each phase.
    activePlayers=2; netTransport.role=NET_ROLE_CLIENT; netLocalPlayerCount=1;
    netLocalPlayerSlots[0]=0;
    players[0].netTetherVisualActive=true; players[0].netTetherVisualTarget=voxels[held].pos;
    for (int phase=0;phase<7;++phase) {
        players[0].meleeSwingActive=true; set_melee_phase(phases[phase]);
        render_gameplay_view(screens,&count,&w,&h,false);
        char name[64]; snprintf(name,sizeof(name),"melee-%s",names[phase]);
        capture_view(screens[0],out,name);
        assert(!weaponFrameSampled);
        // Cached data are exactly those used by both the color and emission passes.
        WeaponMeleePose sample=weaponMeleeSamples[0];
        if (melee_is_active_window(sample.progress)) {
            Vector3 forward,right,up; melee_swing_basis(&players[0],&forward,&right,&up);
            WeaponPose pose=weapon_pose(players[0].pos,forward,right,up,true,weaponViewAspect,&sample,0);
            assert(v_length(v_sub(pose.grip_heel,sample.strike_endpoint))<.0001f);
        }
        netLocalPlayerSlots[0]=1;
        players[0].tetherHolding=false; players[0].netTetherVisualActive=false;
        players[1].pos=(Vector3){-3,2,-1.5f}; players[1].yaw=-90;
        set_melee_phase(phases[phase]);
        render_gameplay_view(screens,&count,&w,&h,false);
        snprintf(name,sizeof(name),"melee-world-%s",names[phase]); capture_view(screens[0],out,name);
        netLocalPlayerSlots[0]=0;
        players[0].tetherHolding=true; players[0].netTetherVisualActive=true;
        players[0].netTetherVisualTarget=voxels[held].pos;
    }
    netTransport.role=NET_ROLE_OFFLINE; activePlayers=2;
    players[0].meleeSwingActive=true; set_melee_phase(.28f);
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee-narrow-active-start");
    set_melee_phase(.52f);
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee-narrow-peak");
    activePlayers=4;
    players[0].meleeSwingActive=true; set_melee_phase(.52f);
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee-quarter");
    players[0].pitch=89; set_melee_phase(.52f);
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee-pitch-up");
    players[0].pitch=-89; set_melee_phase(.52f);
    render_gameplay_view(screens,&count,&w,&h,false); capture_view(screens[0],out,"melee-pitch-down");
    players[0].pitch=0; players[0].last_shot_time=-1000;
    FireVoxel(0); assert(weaponVisuals[0].shot_sequence==2);
    set_melee_phase(.52f); render_gameplay_view(screens,&count,&w,&h,false);
    capture_view(screens[0],out,"melee-firing-tether");
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
