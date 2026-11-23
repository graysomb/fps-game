/* Shared by deterministic backend checks and the GIF frame exporter. */
static void trap_lab_run_preview(const char *name) {
    int which=0;for(int i=1;i<=LAB_FOREST;i++)if(!strcmp(name,trapLabNames[i]))which=i;
    assert(which);
    const char *directory=getenv("FPS_TRAP_FRAMES");
    const char *titles[]={"","Rolling logs","Swinging log","The broken floor","Topple tower","Tower forest"};
    trap_lab_start(which);
    Vector3 target=which>=LAB_TOWER?(Vector3){-.25f,.25f,.25f}:(Vector3){-2.75f,1.25f,-2.25f};
    Ray shot={players[0].pos,v_norm(v_sub(target,players[0].pos))};
    if(which!=LAB_PIT){
        assert(ex_fire(0,EQ_GLOCK,(Ray){players[0].pos,{0,1,0}}));
        assert(!lab.triggered && (!foundryTrapCount || !foundryTraps[0].pending));
        Vector3 mid=Vector3Lerp(shot.position,target,.5f);
        int wall=add_static_voxel_at_grid(floorf(mid.x*2),floorf(mid.y*2),floorf(mid.z*2),GRAY,0);
        assert(wall>=0 && ex_fire(0,EQ_GLOCK,shot));
        assert(!lab.triggered && (!foundryTrapCount || !foundryTraps[0].pending));
        remove_voxel_index(wall);eq_tick(.25f);
    }
    int total=which==LAB_PIT?555:which==LAB_FOREST?360:which==LAB_LOGS?360:255;
    Camera3D camera={.position={11,7,9},.target={0,1,0},.up={0,1,0},.fovy=47};
    if(which==LAB_LOGS)camera=(Camera3D){{17,10,12},{0,1.1f,-.2f},{0,1,0},47};
    if(which==LAB_PIT)camera=(Camera3D){{15,9,7},{0,1.55f,1.6f},{0,1,0},45};
    if(which==LAB_TOWER)camera=(Camera3D){{11,7,11},{.5f,2.9f,2.3f},{0,1,0},46};
    if(which==LAB_FOREST)camera=(Camera3D){{21,13,22},{5.5f,1.9f,4.4f},{0,1,0},47};
    RenderTexture2D rt={0};if(directory){rt=LoadRenderTexture(800,500);character_renderer_init();}
    bool floorOpened=false,heldOpen=false;
    for(int frame=0;frame<total;frame++){
        if(frame==42 && which!=LAB_PIT)assert(ex_fire(0,EQ_GLOCK,shot));
        if(which==LAB_PIT){
            Player *p=&players[0];
            p->vel.z=(frame>=30 && p->pos.z<-.2f) || (frame>=350 && p->pos.z<8.6f)?3.2f:0;
            net_integrate_player(p,1.f/60,true);
        }
        players[0].matter_flash_timer=fmaxf(0,players[0].matter_flash_timer-1.f/60);
        eq_tick(1.f/60);
        simulate_voxel_pbd_steps(1.f/60,1);
        if(which==LAB_PIT){floorOpened|=lab.open;heldOpen|=lab.open && lab.elapsed>4 && players[0].pos.z<1;}
        if(which==LAB_SWING)assert(fabsf(Vector3Distance(lab.logs[0].pos,swingPivot)-4.6f)<.001f);
        if(!directory || frame%3)continue;
        voxel_render_frame++;drawViewPlayerIndex=-1;
        BeginTextureMode(rt);ClearBackground((Color){57,74,83,255});BeginMode3D(camera);
        DrawPlane((Vector3){0,-.02f,0},(Vector2){80,80},(Color){55,68,69,255});
        for(int j=-16;j<=22;j+=2){DrawLine3D((Vector3){-16,0,j},(Vector3){22,0,j},(Color){78,92,90,255});DrawLine3D((Vector3){j,0,-16},(Vector3){j,0,22},(Color){78,92,90,255});}
        DrawVoxels(camera);draw_players();
        if(which!=LAB_PIT && frame>=42 && frame<51)DrawCylinderEx(shot.position,target,.025f,.025f,6,YELLOW);
        EndMode3D();
        DrawRectangle(0,0,800,43,(Color){22,34,40,245});DrawText(titles[which],17,11,22,RAYWHITE);
        const char *state=which==LAB_PIT?(players[0].pos.z>8?"Back on the upper lane":lab.landed?"Survived   /   exit ramp":lab.open?"Floor open":lab.triggered?"Panel shaking":"Step onto the cracked panel"):
          which==LAB_FOREST?TextFormat("%d of 9 released",foundryTraps[0].activation+foundryTraps[1].activation+foundryTraps[2].activation):
          frame<42?"Shoot the orange release":players[4].respawn_timer>0?"Impact   /   side target safe":"Released";
        DrawText(state,783-MeasureText(state,17),15,17,(Color){246,194,112,255});
        if(which==LAB_PIT)DrawText(TextFormat("Matter  %.0f",players[0].matter),17,466,19,RAYWHITE);
        EndTextureMode();Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);
        assert(ExportImage(im,TextFormat("%s/frame-%03d.png",directory,frame/3)));UnloadImage(im);
    }
    fprintf(stderr,"%s: target %.1f safe %.1f matter %.1f pos %.2f %.2f %.2f opened %d held %d credit %d\n",name,players[4].respawn_timer,players[5].respawn_timer,players[0].matter,players[0].pos.x,players[0].pos.y,players[0].pos.z,floorOpened,heldOpen,players[0].debrisKills);
    if(which==LAB_PIT)assert(floorOpened && heldOpen && !lab.open && players[0].matter==75 && players[0].pos.y>=3.95f && players[0].pos.z>=8);
    else assert(players[4].respawn_timer>0 && players[5].respawn_timer<=0 && players[0].debrisKills==1);
    if(which==LAB_LOGS)for(int i=0;i<3;i++)assert(lab.logs[i].angle>5 && lab.logs[i].pos.z>5);
    if(which>=LAB_TOWER)for(int i=0;i<foundryTrapCount;i++){
        assert(foundryTraps[i].activation==(i<3?1:0));
        if(i<3)assert(lab.towers[i].broken && foundryTraps[i].owner==0);
    }
    trap_lab_start(which);assert(!lab.triggered && !lab.open && players[4].respawn_timer<=0);
    if(directory)UnloadRenderTexture(rt);
    fprintf(stderr,"%s PASS (%s)\n",name,physics_backend_name(physicsBackend.active));
}
