/* Interactive fixtures and smoke checks use the same specialist update functions. */
static WorldType enemyLabPreviousWorld;
static bool enemyLabPreviousFoundry,enemyLabPreviousCustom;
static void enemy_lab_exit(void) {
    enemy_cleanup(true);enemyLab=false;enemyLabAll=false;enemyLabPaused=false;
    currentWorldType=enemyLabPreviousWorld;firefightFoundry=enemyLabPreviousFoundry;useCustomMap=enemyLabPreviousCustom;
    gameState=GAME_STATE_MENU;couchLobby=false;
}
static int enemy_parse_case(const char *name) {
    if(!name || !strcmp(name,"all")) return ENEMY_TYPE_SCRAPPER;
    if(!strcmp(name,"mixed")) return ENEMY_TYPE_COUNT;
    for(int i=ENEMY_TYPE_SCRAPPER;i<ENEMY_TYPE_COUNT;i++) if(!strcmp(name,enemy_names[i])) return i;
    return -1;
}
static void enemy_lab_spawn(int slot,EnemyType type,Vector3 position) {
    init_firefight_bot_entity(slot,type,INPUT_TYPE_BOT_MEDIUM);
    players[slot].pos=players[slot].death_pos=position;players[slot].yaw=0;players[slot].invuln_timer=0;
    if(enemy_is_special(slot)) enemy_update_equipment(slot);
}
static void enemy_lab_start(int which) {
    if(!enemyLab) {enemyLabPreviousWorld=currentWorldType;enemyLabPreviousFoundry=firefightFoundry;enemyLabPreviousCustom=useCustomMap;enemyLabGod=true;enemyLabPaused=false;enemyLabStep=false;}
    enemyLab=true;enemyLabCase=which;enemyLabClock=0;
    gameMode=GAME_MODE_FIREFIGHT;firefightPractice=false;firefightFoundry=false;useCustomMap=false;
    currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=1;firefightHumans=1;
    ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();worldFloorHalfSize=32;
    SetRandomSeed(7);srand(7);activePlayers=MAX_COMBATANTS;
    for(int i=0;i<MAX_COMBATANTS;i++) {players[i].respawn_timer=9999;combatantLife[i]=LIFE_UNUSED;}
    players[0].pos=players[0].death_pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[0].respawn_timer=0;
    players[0].matter=players[0].matterMax=100;players[0].isExposed=false;players[0].yaw=180;
    playerInput[0]=multiplayerPlayerInput[0];combatantLife[0]=LIFE_ALIVE;combatantTeam[0]=0;
    firefightWaveEnemiesTotal=firefightWaveEnemiesRemaining=999;firefightEnemiesSpawned=999;
    firefightWaveBannerTimer=0;
    memset(ffRemaining,0,sizeof(ffRemaining));firefightWaveState=FIREFIGHT_STATE_ACTIVE;
    if(which==ENEMY_TYPE_COUNT) {
        for(int t=ENEMY_TYPE_SCRAPPER;t<ENEMY_TYPE_COUNT;t++) {
            int n=t-ENEMY_TYPE_SCRAPPER;
            enemy_lab_spawn(4+n,(EnemyType)t,(Vector3){(n%3-1)*8,BASE_EYE_HEIGHT,10+n/3*7});
        }
    } else enemy_lab_spawn(4,(EnemyType)which,(Vector3){0,BASE_EYE_HEIGHT,which==ENEMY_TYPE_LOBBER?16:10});
    if(which==ENEMY_TYPE_SAPPER || which==ENEMY_TYPE_COUNT) {
        add_static_box_at_grid(-3,3,0,5,10,10,GRAY);
        add_static_box_at_grid(-6,6,6,6,9,14,LIGHTGRAY);
    }
    if(which==ENEMY_TYPE_SCRAPPER || which==ENEMY_TYPE_FOUNDRY || which==ENEMY_TYPE_COUNT) {
        for(int n=0;n<12;n++) {
            int v=addVoxel(-3+(n%4)*.6f,.3f,8+(n/4)*.6f,false,true,BROWN,0);
            if(v>=0) {voxels[v].owner=voxels[v].activator=-1;set_voxel_velocity(&voxels[v],(Vector3){0});}
        }
    }
    if(which==ENEMY_TYPE_SHIELDBEARER || which==ENEMY_TYPE_ANCHOR || which==ENEMY_TYPE_BUILDER) {
        enemy_lab_spawn(5,ENEMY_TYPE_STANDARD,(Vector3){3,BASE_EYE_HEIGHT,11});
        if(which==ENEMY_TYPE_BUILDER) {
            add_static_box_at_grid(-4,4,0,3,-6,-3,GRAY);
            players[0].pos=(Vector3){0,BASE_EYE_HEIGHT+2,-2};
        }
    }
    init_static_hash();rebuild_voxel_hash();rebuild_glue_constraints();rebuild_all_voxel_surfaces();meshDirty=true;
    gameState=GAME_STATE_PLAYING;couchLobby=false;couchAdvanced=false;
}
static bool enemy_lab_controls(void) {
    enemyLabStep=false;
    if(IsKeyPressed(KEY_F3)||couch_button(0,GAMEPAD_BUTTON_LEFT_TRIGGER_1)) {
        int next=enemyLabCase+1;if(next>ENEMY_TYPE_COUNT) next=ENEMY_TYPE_SCRAPPER;
        enemy_lab_start(next);
    }
    if(IsKeyPressed(KEY_F4)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_UP)) enemy_lab_start(enemyLabCase);
    if(IsKeyPressed(KEY_F5)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_LEFT)) enemyLabPaused=!enemyLabPaused;
    if(IsKeyPressed(KEY_F6)||couch_button(0,GAMEPAD_BUTTON_RIGHT_TRIGGER_1)) {enemyLabPaused=true;enemyLabStep=true;}
    if(IsKeyPressed(KEY_F7)||couch_button(0,GAMEPAD_BUTTON_MIDDLE_LEFT)) enemyLabGod=!enemyLabGod;
    if(enemyLabAll && enemyLabClock>=30) {int n=enemyLabCase+1;enemy_lab_start(n>=ENEMY_TYPE_COUNT?ENEMY_TYPE_SCRAPPER:n);}
    if(IsKeyPressed(KEY_ESCAPE)) enemy_lab_exit();
    return true;
}
static void enemy_lab_hud(int w,int h) {
    int fs=w<900?14:18;
    DrawRectangle(8,70,w-16,130,Fade(BLACK,.8f));
    DrawText(TextFormat("ENEMY LAB: %s | %s | Invulnerable: %s",enemyLabCase==ENEMY_TYPE_COUNT?"mixed":enemy_names[enemyLabCase],enemyLabPaused?"PAUSED":"RUNNING",enemyLabGod?"ON":"OFF"),16,78,fs,GOLD);
    DrawText("F3/LB: next | F4/Top: reset | F5/Left face: pause",16,102,fs,RAYWHITE);
    DrawText("F6/RB: step | F7/Back button: invulnerability | Esc: exit",16,126,fs,RAYWHITE);
    if(enemyLabCase<ENEMY_TYPE_COUNT) DrawText(enemy_hints[enemyLabCase-ENEMY_TYPE_SCRAPPER],16,150,fs,RAYWHITE);
    EnemySpecialState *s=&enemySpecial[4];
    DrawText(TextFormat("State %d | Target P%d | Cooldown %.1f | Wind-up %.1f | Equipment %d",s->state,s->target+1,fmaxf(0,s->cooldown),s->windup,s->equipmentCount),16,174,fs,SKYBLUE);
    (void)h;
}
#define ENEMY_CHECK(condition,message) do {if(!(condition)) {fprintf(stderr,"enemy-smoke FAIL: %s\n",message);return false;}} while(0)
static void enemy_lab_capture(int type) {
    if(!DirectoryExists("artifacts")) MakeDirectory("artifacts");
    RenderTexture2D shot=LoadRenderTexture(960,720);
    Camera3D cam={.position={12,12,21},.target={0,1.5f,9},.up={0,1,0},.fovy=50,.projection=CAMERA_PERSPECTIVE};
    BeginTextureMode(shot);ClearBackground((Color){24,28,36,255});BeginMode3D(cam);
    DrawPlane((Vector3){0,0,8},(Vector2){30,30},(Color){65,70,76,255});
    for(int v=0;v<voxel_count;v++) DrawCube(voxels[v].pos,.5f,.5f,.5f,voxels[v].color);
    for(int i=4;i<activePlayers;i++) if(players[i].respawn_timer<=0) {
        DrawCube(players[i].pos,.5f,.5f,.5f,enemy_color(players[i].enemyType));
        DrawCubeWires(players[i].pos,.5f,.5f,.5f,WHITE);
    }
    enemy_draw();EndMode3D();DrawText(enemy_names[type],24,24,28,RAYWHITE);
    DrawText(enemy_hints[type-ENEMY_TYPE_SCRAPPER],24,660,18,GOLD);EndTextureMode();
    Image im=LoadImageFromTexture(shot.texture);ImageFlipVertical(&im);
    ExportImage(im,TextFormat("artifacts/enemy-%s.png",enemy_names[type]));UnloadImage(im);UnloadRenderTexture(shot);
}
static bool enemy_smoke_case(int type) {
    enemy_lab_start(type);enemyLabPaused=false;int i=4;
    enemy_lab_capture(type);
    EnemySpecialState *s=&enemySpecial[i];
    if(type==ENEMY_TYPE_SCRAPPER) {
        s->cooldown=0;enemy_update(i,.1f);ENEMY_CHECK(s->windup>0,"scrapper starts collection wind-up");
        enemy_update(i,1.1f);ENEMY_CHECK(s->equipmentCount==1,"scrapper consumes actual debris");
        int v=enemy_voxel(s->equipment[0].id);ENEMY_CHECK(v>=0 && !enemy_scrap_eligible(v,i),"equipment is ineligible scrap");
        enemy_detach(v);ENEMY_CHECK(!s->equipmentCount,"scrapper armor can be stripped");
        for(int n=0;n<20;n++) {s->cooldown=0;enemy_update(i,.01f);enemy_update(i,1.1f);}
        ENEMY_CHECK(s->equipmentCount<=12,"scrapper armor stays bounded");
    } else if(type==ENEMY_TYPE_SHIELDBEARER) {
        rebuild_voxel_hash();
        int front=first_voxel_hit((Ray){{0,BASE_EYE_HEIGHT,6},{0,0,1}},4,-1);
        ENEMY_CHECK(front>=0 && enemy_equipment_owner(front)==i,"shield intercepts frontal shot path");
        ENEMY_CHECK(first_voxel_hit((Ray){{3,BASE_EYE_HEIGHT,10},{-1,0,0}},2.5f,-1)<0,"shield leaves side vulnerable");
        float before=players[i].yaw;players[0].pos.x=10;enemy_update(i,.1f);
        ENEMY_CHECK(fabsf(angle_diff_deg(before,players[i].yaw))<=6.01f,"shield turn-rate cap");
        enemy_detach(front);ENEMY_CHECK(s->equipmentCount==19,"shield pieces detach");
    } else if(type==ENEMY_TYPE_SAPPER) {
        s->selected=enemy_sapper_target(i,0);int support=enemy_voxel(s->selected);
        ENEMY_CHECK(support>=0,"sapper finds structural support");
        players[i].pos=v_add(voxels[support].pos,(Vector3){0,0,1});
        ENEMY_CHECK(enemy_plant(i),"sapper plants charge");enemy_tick(1);
        ENEMY_CHECK(enemyCharges[i].active,"charge honors fuse");
        /* Reusing the old Sapper's actor slot must not turn its planted charge invulnerable. */
        init_firefight_bot_entity(i,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_EASY);
        int v=enemy_voxel(enemyCharges[i].id);ENEMY_CHECK(enemy_equipment_owner(v)==i,"charge identity survives actor reuse");
        detach_goliath_armor_voxel(v);enemy_tick(.1f);
        ENEMY_CHECK(!enemyCharges[i].active && enemy_voxel(s->selected)>=0,"shooting charge disarms without blast");
        enemy_lab_start(type);s=&enemySpecial[i];s->selected=enemy_sapper_target(i,0);
        support=enemy_voxel(s->selected);players[i].pos=v_add(voxels[support].pos,(Vector3){0,0,1});
        ENEMY_CHECK(enemy_plant(i),"charge can be replanted");
        remove_voxel_index(enemy_voxel(s->selected));enemy_tick(.1f);
        ENEMY_CHECK(!enemyCharges[i].active,"removing support disarms charge");
        enemy_lab_start(type);s=&enemySpecial[i];s->selected=enemy_sapper_target(i,0);
        support=enemy_voxel(s->selected);players[i].pos=v_add(voxels[support].pos,(Vector3){0,0,1});
        ENEMY_CHECK(enemy_plant(i),"explosive fixture");uint64_t supportId=s->selected;
        enemy_tick(4.1f);ENEMY_CHECK(!enemyCharges[i].active && enemy_voxel(supportId)<0,"detonation destroys local support");
    } else if(type==ENEMY_TYPE_LOBBER) {
        s->aim=players[0].pos;ENEMY_CHECK(enemy_launch(i),"lobber ballistic launch");
        ENEMY_CHECK(enemy_launch(i) && !enemy_launch(i),"two-projectile limit");
        int v=enemy_voxel(enemyLobs[MAX_COMBATANTS*2-1].id);ENEMY_CHECK(v>=0,"projectile exists");
        tetherTag[v]=1;enemy_tick(.1f);ENEMY_CHECK(voxels[v].owner==0,"catch transfers attribution");tetherTag[v]=0;
        Vector3 velocity;add_static_box_at_grid(-4,4,7,7,24,34,GRAY);init_static_hash();
        ENEMY_CHECK(!enemy_lob_velocity(v_add(players[i].pos,(Vector3){0,2.3f,0}),s->aim,&velocity),"ceiling blocks trajectory");
    } else if(type==ENEMY_TYPE_BUILDER) {
        s->buildOrigin=(Vector3){0,0,6};s->buildKind=1;
        ENEMY_CHECK(enemy_build(i),"builder makes traversable ramp");
        ENEMY_CHECK(!enemy_build(i),"builder cannot overlap existing structure");
        s->buildOrigin=(Vector3){10,0,6};s->buildKind=0;ENEMY_CHECK(enemy_build(i),"second structure");
        s->buildOrigin=(Vector3){-10,0,6};ENEMY_CHECK(!enemy_build(i),"two-structure cap");
        s->buildOrigin=players[0].pos;ENEMY_CHECK(!enemy_build_valid(i,s->buildOrigin,0),"no building through players");
        uint64_t block=enemyStructures[0].ids[0];int v=enemy_voxel(block);ENEMY_CHECK(v>=0,"construction has physical voxels");
        int recycledBefore=recycleQueueCount,gx=voxels[v].gx,gy=voxels[v].gy,gz=voxels[v].gz;
        remove_static_voxels_in_region_recycle(gx,gx,gy,gy,gz,gz);
        ENEMY_CHECK(enemy_voxel(block)<0 && recycleQueueCount==recycledBefore,"harvested construction never restores as permanent map geometry");
    } else if(type==ENEMY_TYPE_ANCHOR) {
        ENEMY_CHECK(enemy_knockback_scale(5)==.25f,"anchor braces nearby ally");
        players[5].pos.x=20;ENEMY_CHECK(enemy_knockback_scale(5)==1,"anchor range enforced");players[5].pos.x=3;
        add_static_box_at_grid(3,3,0,5,18,24,GRAY);init_static_hash();
        ENEMY_CHECK(enemy_knockback_scale(5)==1,"solid cover breaks anchor line of sight");
        remove_static_voxels_in_region(3,3,0,5,18,24);init_static_hash();
        int v=enemy_voxel(s->equipment[0].id);enemy_detach(v);
        ENEMY_CHECK(enemy_knockback_scale(5)==1 && s->disabled==8,"anchor block breaks link immediately");
    } else if(type==ENEMY_TYPE_FOUNDRY) {
        float matter=players[i].matter;apply_matter_damage(i,0,100);
        ENEMY_CHECK(players[i].matter==matter,"boss core starts protected");
        for(int n=s->equipmentCount-1;n>=0;n--) if(s->equipment[n].group<2) enemy_detach(enemy_voxel(s->equipment[n].id));
        enemy_update_equipment(i);ENEMY_CHECK(enemy_core_open(i),"two supports expose core");
        apply_matter_damage(i,0,20);ENEMY_CHECK(players[i].matter==matter-20,"exposed core takes damage");
        enemy_tick(10.1f);ENEMY_CHECK(!enemy_core_open(i),"partial exposure expires");
        for(int n=s->equipmentCount-1;n>=0;n--) if(s->equipment[n].group<4) enemy_detach(enemy_voxel(s->equipment[n].id));
        enemy_update_equipment(i);enemy_tick(11);ENEMY_CHECK(enemy_core_open(i),"all supports expose core permanently");
        for(int n=0;n<8;n++) {s->spawnTimer=0;enemy_update(i,.1f);s->windup=0;}
        ENEMY_CHECK(s->spawned==6,"boss reinforcement budget");
        s->cooldown=100;s->windup=0;
        for(int n=0;n<30;n++) {s->repair=0;enemy_update(i,.01f);}
        ENEMY_CHECK(s->repairs==12 && enemy_core_open(i),"boss repair budget never restores destroyed supports");
    }
    s->windup=1;s->state=3;players[i].isExposed=true;enemy_update(i,.01f);
    ENEMY_CHECK(s->windup==0,"exposure interrupts ability wind-up");
    enemy_cleanup(true);
    ENEMY_CHECK(enemySpecial[4].equipmentCount==0 && !enemyCharges[4].active,"cleanup clears equipment and charges");
    for(int n=0;n<ENEMY_STRUCTURE_MAX;n++) ENEMY_CHECK(!enemyStructures[n].count,"cleanup clears construction");
    fprintf(stderr,"enemy-smoke PASS: %s\n",enemy_names[type]);return true;
}
static bool enemy_run_smoke(const char *which) {
    SetTargetFPS(0);
    int type=enemy_parse_case(which);if(type<0) {fprintf(stderr,"Unknown enemy case: %s\n",which);return false;}
    if(!strcmp(which,"all") || type==ENEMY_TYPE_COUNT) {
        for(int t=ENEMY_TYPE_SCRAPPER;t<ENEMY_TYPE_COUNT;t++) if(!enemy_smoke_case(t)) return false;
    } else if(!enemy_smoke_case(type)) return false;
    /* Sustained mixed simulation exercises compaction, ownership, physics and render paths. */
    for(int scenario=0;scenario<2;scenario++) {
    enemy_lab_start(scenario?ENEMY_TYPE_COUNT:ENEMY_TYPE_FOUNDRY);double samples[30];
    RenderTexture2D views[MAX_PLAYERS]={0};int count=0,w=0,h=0;
    for(int i=1;i<4;i++) {players[i]=players[0];players[i].pos.x=i*2;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=0;playerInput[i]=INPUT_TYPE_KEYBOARD;}
    firefightHumans=multiplayerActivePlayers=4;
    for(int frame=0;frame<600;frame++) {
        float dt=1.0f/30;double start=GetTime();enemy_tick(dt);
        for(int i=4;i<activePlayers;i++) if(players[i].respawn_timer<=0) {UpdateBot(i,dt);net_integrate_player(&players[i],dt,true);}
        simulate_voxel_pbd_steps(1.0f/120,4);update_projectiles(dt);handle_pbd_projectile_hits();
        if(frame>=570) {render_gameplay_view(views,&count,&w,&h,false);samples[frame-570]=(GetTime()-start)*1000;}
        ENEMY_CHECK(voxel_count<MAX_VOXELS && activePlayers<=MAX_COMBATANTS,"bounded mixed simulation");
        for(int i=4;i<activePlayers;i++) {
            ENEMY_CHECK(enemySpecial[i].equipmentCount<=ENEMY_EQUIPMENT_MAX && enemySpecial[i].repairs<=12,"equipment and repair budgets");
            ENEMY_CHECK(isfinite(players[i].pos.x) && isfinite(players[i].pos.y),"finite actor positions");
        }
    }
    for(int i=0;i<30;i++) for(int j=i+1;j<30;j++) if(samples[j]<samples[i]) {double t=samples[i];samples[i]=samples[j];samples[j]=t;}
    fprintf(stderr,"enemy-smoke %s four-view (%s): p50 %.2f ms p95 %.2f ms\n",scenario?"mixed":"boss",physics_backend_name(physicsBackend.active),samples[15],samples[28]);
    if(count) {
        Image im=LoadImageFromTexture(views[0].texture);ImageFlipVertical(&im);
        ExportImage(im,scenario?"artifacts/enemy-lab-quarter.png":"artifacts/enemy-boss-quarter.png");UnloadImage(im);
    }
    for(int i=0;i<count;i++) UnloadRenderTexture(views[i]);
    }
    enemy_cleanup(true);enemyLab=false;return true;
}
#undef ENEMY_CHECK
