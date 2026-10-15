#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void empty_periodic(void) {
    clear_all_chunks();clear_world_voxels();activePlayers=0;
    worldFloorHalfSize=32;worldPeriodic=true;
}
static void topology_checks(void) {
    empty_periodic();
    assert(world_wrap_axis(32)==-32 && world_wrap_axis(-32)==-32);
    assert(world_wrap_axis(32+64*12)==-32 && world_grid(64)==-64);
    Vector3 d=world_delta((Vector3){-31,2,-31},(Vector3){31,1,31});
    assert(d.x==2 && d.y==1 && d.z==2);
    add_static_voxel_at_grid(-64,2,-64,GRAY,0);
    assert(table_get_static_only(64,2,64)>=0);
    VoxelHit hit;
    assert(first_voxel_hit_detailed((Ray){{31,1.25f,-31.75f},{1,0,0}},4,-1,true,&hit));
    assert(hit.x==-64 && hit.z==-64 && fabsf(hit.t-1)<.001f);
    Particle probe={0};probe.predicted_pos=(Vector3){31.95f,1.25f,32.25f};
    assert(push_particle_out_of_static(&voxels[hit.id],&probe,.125f));
    assert(probe.predicted_pos.x<31.95f);
    assert(ensure_static_surface_cache());
    probe.predicted_pos=(Vector3){31.95f,1.25f,32.25f};
    collide_particle_with_static_surface(&probe,.125f);
    assert(probe.predicted_pos.x<=31.876f);
    clear_all_chunks();clear_world_voxels();
    activePlayers=1;playerInput[0]=INPUT_TYPE_KEYBOARD;
    players[0].pos=(Vector3){31.9f,BASE_EYE_HEIGHT,31.9f};players[0].vel=(Vector3){6,0,6};players[0].onGround=true;
    net_integrate_player(&players[0],.1f,true);
    assert(players[0].pos.x<0 && players[0].pos.z<0 && players[0].vel.x>0);
    GrenadeState g={.pos={31.9f,5,0},.vel={10,0,0},.fuse=-1};
    clarity_grenade_step(&g,.1f);assert(g.pos.x<0 && g.fuse==-1);
    worldPeriodic=false;assert(world_grid(64)==64 && world_wrap_axis(33)==33);
    puts("Periodic coordinates, corners, seam ray, static contact, player and grenade passed");
}
static void physics_checks(void) {
    empty_periodic();
    int id=addVoxel(31.8f,3,0,false,true,GRAY,0);assert(id>=0);
    Voxel *v=&voxels[id];
    for(int c=0;c<8;c++)v->particles[c]->vel=(Vector3){12,0,0};
    v->vel=(Vector3){12,0,0};
    for(int f=0;f<24;f++)simulate_voxel_pbd_steps(1.f/60,1);
    assert(v->pos.x<0 && isfinite(v->pos.y));
    for(int c=1;c<8;c++)assert(Vector3Distance(v->particles[c]->pos,v->particles[0]->pos)<1.2f);
    assert(!voxel_outside_world_bounds(v));
    fprintf(stderr,"Periodic physics backend=%s position=%.3f,%.3f,%.3f\n",physics_backend_name(physicsBackend.active),v->pos.x,v->pos.y,v->pos.z);
    empty_periodic();
    int a=addVoxel(31.75f,3.25f,.25f,false,true,GRAY,0);
    int b=addVoxel(-31.75f,3.25f,.25f,false,true,GRAY,0);
    rebuild_voxel_hash();glue_neighbor_faces();
    assert(Vector3Distance(voxels[a].pos,voxels[b].pos)<.6f);
    for(int i=0;i<sim_particle_count;i++)sim_particles[i]->vel=(Vector3){-12,0,0};
    for(int f=0;f<24;f++)simulate_voxel_pbd_steps(1.f/60,1);
    assert(Vector3Distance(voxels[a].pos,voxels[b].pos)<.7f);
    for(int i=0;i<2;i++)for(int c=1;c<8;c++)assert(Vector3Distance(voxels[i].particles[c]->pos,voxels[i].particles[0]->pos)<1.2f);
    puts("Connected island assembled across seam remains intact");
}
static void mesh_checks(void) {
    empty_periodic();worldFloorHalfSize=40;
    add_static_voxel_at_grid(-80,2,0,GRAY,0);
    add_static_voxel_at_grid(79,2,0,GRAY,0);
    int count=0;
    for(int i=0;i<voxel_chunk_count;i++) {
        VoxelChunk *c=&voxel_chunks[i];rebuild_chunk(c);count+=c->voxel_count;
        for(int v=0;v<c->mesh.vertexCount;v++)assert(fabsf(c->mesh.vertices[v*3])<=40);
    }
    assert(count==2);
    for(int i=0;i<1000;i++) {
        Vector3 p=world_wrap((Vector3){39.9f+80*i,3,-39.9f-80*i});
        assert(fabsf(p.x-39.9f)<.01f && fabsf(p.z+39.9f)<.01f);
    }
    assert(voxel_count==2);
    puts("Non-aligned Foundry chunk ownership and repeated laps passed");
}
static void stress_checks(void) {
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;
    firefightHumans=multiplayerActivePlayers=1;ResetGame();activePlayers=0;
    /* Destructible blocks straddle the seam; later blasts hit real map cover. */
    add_static_box_at_grid(78,79,0,5,-4,4,GRAY);
    add_static_box_at_grid(-80,-79,0,5,-4,4,GRAY);
    int original=voxel_count;double worst=0;
    for(int f=0;f<180;f++) {
        double start=GetTime();
        if(f%30==0) {
            Vector3 at={39,1,0};
            if(f)for(int v=voxel_count/2;v<voxel_count;v++)if(!voxels[v].simulate && voxels[v].pos.y<3){at=voxels[v].pos;break;}
            clarity_explode(at,0);
        }
        simulate_voxel_pbd_steps(1.f/60,1);
        for(int i=0;i<voxel_count;i++)if(voxels[i].simulate) {
            assert(isfinite(voxels[i].pos.x) && isfinite(voxels[i].pos.y) && isfinite(voxels[i].pos.z));
            assert(fabsf(voxels[i].pos.x)<2*FLOOR_SIZE && fabsf(voxels[i].pos.z)<2*FLOOR_SIZE);
        }
        worst=fmax(worst,(GetTime()-start)*1000);
    }
    assert(voxel_count<original);
    fprintf(stderr,"Foundry seam blast stress: 6 blasts/180 steps, %d voxels, worst step %.2f ms\n",voxel_count,worst);
}
static void capture_maps(void) {
    voxel_ca_init();character_renderer_init();weaponRenderingReady=weapon_renderer_init();assert(weaponRenderingReady);
    for(int map=0;map<2;map++) {
        gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;
        currentWorldType=WORLD_TYPE_CROSSCUT;multiplayerActivePlayers=1;ResetGame();assert(worldPeriodic);
        RenderTexture2D target=world_render_target(960,540);
        Color *seamBefore=NULL;
        for(int view=0;view<5;view++) {
            Camera3D cam={{12,4,21},{0,3,0},{0,1,0},65,CAMERA_PERSPECTIVE};
            if(view==1){cam.position=(Vector3){FLOOR_SIZE-2,3,0};cam.target=(Vector3){FLOOR_SIZE+12,3,0};}
            if(view==2){cam.position=(Vector3){FLOOR_SIZE-2,3,FLOOR_SIZE-2};cam.target=(Vector3){FLOOR_SIZE+12,3,FLOOR_SIZE+12};}
            if(view>=3){cam.position=world_wrap((Vector3){FLOOR_SIZE+(view==3?-.0001f:.0001f),3,0});cam.target=v_add(cam.position,(Vector3){10,0,0});}
            voxel_render_frame++;BeginTextureMode(target);ClearBackground(periodicSky);periodic_skyline(cam);
            for(int image=0;image<9;image++) {
                Camera3D c;if(!periodic_camera(cam,image,&c))continue;
                BeginMode3D(c);draw_world_surfaces();DrawVoxels(c);draw_pickups(c);EndMode3D();
            }
            periodic_fog(target,cam,false,target.depth);EndTextureMode();
            Image image=LoadImageFromTexture(target.texture);ImageFlipVertical(&image);
            if(view==3)seamBefore=LoadImageColors(image);
            if(view==4) {
                Color *after=LoadImageColors(image);int changed=0;
                for(int p=0;p<960*540;p++)if(abs(after[p].r-seamBefore[p].r)>4 || abs(after[p].g-seamBefore[p].g)>4 || abs(after[p].b-seamBefore[p].b)>4)changed++;
                fprintf(stderr,"Map %d seam continuity: %.3f%% changed pixels\n",map,changed*100.0/(960*540));
                assert(changed<960*540/20);UnloadImageColors(after);UnloadImageColors(seamBefore);
            }
            ExportImage(image,TextFormat("artifacts/periodic-map%d-view%d.png",map,view));UnloadImage(image);
        }
        UnloadRenderTexture(target);
    }
}
static int compare_time(const void *a,const void *b) {
    double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);
}
static void profile_views(void) {
    for(int map=0;map<2;map++)for(int count=1;count<=4;count+=3) {
        gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightHumans=multiplayerActivePlayers=count;
        currentWorldType=WORLD_TYPE_CROSSCUT;firefightFoundry=true;ResetGame();gameState=GAME_STATE_PLAYING;
        firefightWaveState=FIREFIGHT_STATE_ACTIVE;firefightWaveBannerTimer=0;
        for(int p=0;p<count;p++){playerInput[p]=INPUT_TYPE_KEYBOARD;players[p].pos=(Vector3){FLOOR_SIZE-2,BASE_EYE_HEIGHT,p*.8f};players[p].yaw=-90;}
        RenderTexture2D targets[4]={0};int views=0,w=0,h=0;
        for(int periodic=0;periodic<2;periodic++) {
            worldPeriodic=periodic;double time[60];
            for(int f=0;f<70;f++) {
                double start=GetTime();voxel_ca_update(start);render_gameplay_view(targets,&views,&w,&h,false);
                if(f>=10)time[f-10]=(GetTime()-start)*1000;
            }
            qsort(time,60,sizeof(double),compare_time);
            fprintf(stderr,"map=%d views=%d periodic=%d render p50=%.3f p95=%.3f ms\n",map,count,periodic,time[30],time[57]);
        }
        Image image=GenImageColor(960,540,BLACK);
        for(int p=0;p<views;p++) {
            Image part=LoadImageFromTexture(targets[p].texture);ImageFlipVertical(&part);
            int x,y,pw,ph;get_viewport(p,views,&x,&y,&pw,&ph);
            ImageDraw(&image,part,(Rectangle){0,0,part.width,part.height},(Rectangle){x,y,pw,ph},WHITE);UnloadImage(part);
        }
        ExportImage(image,TextFormat("artifacts/periodic-map%d-players%d.png",map,count));UnloadImage(image);
        for(int p=0;p<views;p++)UnloadRenderTexture(targets[p]);
    }
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Periodic world regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    topology_checks();physics_checks();mesh_checks();stress_checks();capture_maps();profile_views();
    clear_all_chunks();clear_world_voxels();shutdown_world_visuals();
    assert(!periodicFogScratch.id && !periodicIntroScene.id && !periodicFogShader.id);
    for(int i=0;i<MAX_PLAYERS;i++)weapon_bloom_unload(&weaponBlooms[i]);
    character_renderer_shutdown();weapon_renderer_shutdown();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Periodic world regression passed");return 0;
}
