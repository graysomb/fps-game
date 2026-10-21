/* Interactive fixtures and smoke checks use the same specialist update functions. */
static WorldType enemyLabPreviousWorld;
static int enemyLabFixture=-1;
static float enemyLabAiTime,enemyLabAiNext;
static const char *enemyLabFixtureNames[]={"Grenade range: C / D-pad Up; full friendly fire", "Bomb avoidance: enemies flee a live Sapper charge", "Death feedback: exposed enemies, slot reuse", "Melee shield: remove one section per swing", "Lobber: elevated target and real structure flight", "Anchor: damage linked allies; exposed recovery after 2 seconds", "Footsteps: listen while enemies move", "Tether accuracy: throw the side-offset block at the target", "AI rear attack: fire at the bot facing away", "AI near miss: fire past the bot shoulder", "AI shield alert: shoot the shield from either side", "AI alternating attackers: two scripted shooters", "AI occlusion: attacker behind wall; investigate without tracking", "AI stationary target: observe bursts and recovery", "AI reversing target: observe tracking accuracy"};
static bool enemyLabPreviousFoundry,enemyLabPreviousCustom;
static void enemy_lab_exit(void) {
    enemy_cleanup(true);enemyLab=false;enemyLabAll=false;enemyLabPaused=false;
    currentWorldType=enemyLabPreviousWorld;firefightFoundry=enemyLabPreviousFoundry;useCustomMap=enemyLabPreviousCustom;
    gameState=GAME_STATE_MENU;couchLobby=false;
}
static int enemy_parse_case(const char *name) {
    if(!name || !strcmp(name,"all") || !strcmp(name,"scrapper")) return ENEMY_TYPE_GOLIATH;
    if(!strcmp(name,"goliath"))return ENEMY_TYPE_GOLIATH;
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
    if(which==ENEMY_TYPE_SCRAPPER)which=ENEMY_TYPE_GOLIATH;
    if(!enemyLab) {
        enemyLabPreviousWorld=currentWorldType;enemyLabPreviousFoundry=firefightFoundry;enemyLabPreviousCustom=useCustomMap;
        enemyLabGod=true;enemyLabPaused=false;enemyLabStep=false;
        /* The lab skips the lobby, so bind its human seat on entry. Resets and
         * case changes retain the seat, including explicit device reassignment. */
        int device=couch_menu_device();
        multiplayerPlayerInput[0]=device>=0?INPUT_TYPE_GAMEPAD:INPUT_TYPE_KEYBOARD;
        if(device>=0)localSeats[0].device=device;
    }
    enemyLabFixture=-1;
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
        for(int t=ENEMY_TYPE_GOLIATH;t<ENEMY_TYPE_COUNT;t++) {
            if(t==ENEMY_TYPE_SCRAPPER)continue;
            int n=t-ENEMY_TYPE_GOLIATH-(t>ENEMY_TYPE_SCRAPPER);
            enemy_lab_spawn(4+n,(EnemyType)t,(Vector3){(n%3-1)*8,BASE_EYE_HEIGHT,10+n/3*7});
        }
    } else enemy_lab_spawn(4,(EnemyType)which,(Vector3){0,BASE_EYE_HEIGHT,which==ENEMY_TYPE_LOBBER?16:10});
    if(which==ENEMY_TYPE_SAPPER || which==ENEMY_TYPE_COUNT) {
        add_static_box_at_grid(-3,3,0,5,10,10,GRAY);
        add_static_box_at_grid(-6,6,6,6,9,14,LIGHTGRAY);
    }
    if(which==ENEMY_TYPE_GOLIATH || which==ENEMY_TYPE_FOUNDRY || which==ENEMY_TYPE_COUNT) {
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
    if(which==ENEMY_TYPE_LOBBER || which==ENEMY_TYPE_COUNT) {
        enemy_scrap_bundle((Vector3){-3,0,12});enemy_scrap_bundle((Vector3){3,0,12});
    }
    init_static_hash();rebuild_voxel_hash();rebuild_glue_constraints();rebuild_all_voxel_surfaces();meshDirty=true;
    gameState=GAME_STATE_PLAYING;couchLobby=false;couchAdvanced=false;
}
static void enemy_lab_fixture(int which) {
    which=(which+15)%15;
    enemy_lab_start(which>=8?(which==10?ENEMY_TYPE_SHIELDBEARER:ENEMY_TYPE_STANDARD):which==3?ENEMY_TYPE_SHIELDBEARER:which==4?ENEMY_TYPE_LOBBER:which==5?ENEMY_TYPE_ANCHOR:which==6?ENEMY_TYPE_SWARMER:ENEMY_TYPE_SAPPER);
    enemyLabFixture=which;enemyLabAiTime=enemyLabAiNext=0;grenadeInventory[0]=2;
    if(which==0){enemy_cleanup(true);players[4].respawn_timer=9999;add_static_box_at_grid(-4,4,0,5,16,16,GRAY);enemy_lab_spawn(5,ENEMY_TYPE_STANDARD,(Vector3){3,BASE_EYE_HEIGHT,8});}
    if(which==1){int v=addVoxel(0,.3f,9,false,true,ORANGE,0);enemyCharges[4]=(EnemyCharge){.active=true,.owner=4,.id=voxels[v].identity,.pos={0,.3f,9},.fuse=4};for(int i=5;i<8;i++)enemy_lab_spawn(i,ENEMY_TYPE_STANDARD,(Vector3){(i-6)*2,BASE_EYE_HEIGHT,10});}
    if(which==2){enemy_cleanup(true);players[4].respawn_timer=9999;for(int i=5;i<8;i++){enemy_lab_spawn(i,ENEMY_TYPE_SWARMER,(Vector3){(i-6)*3,BASE_EYE_HEIGHT,6});players[i].matter=0;players[i].isExposed=true;}}
    if(which==3){players[4].pos=(Vector3){0,BASE_EYE_HEIGHT,2};players[4].yaw=0;enemy_update_equipment(4);}
    if(which==4){add_static_box_at_grid(-3,3,0,5,-3,3,GRAY);players[0].pos.y=BASE_EYE_HEIGHT+3;}
    if(which==5){for(int i=5;i<8;i++){enemy_lab_spawn(i,ENEMY_TYPE_STANDARD,(Vector3){(i-6)*2,BASE_EYE_HEIGHT,8});players[i].matter=i==5?0:20;players[i].isExposed=i==5;}}
    if(which==7){enemy_cleanup(true);players[4].respawn_timer=9999;add_static_box_at_grid(-2,2,0,8,40,40,ORANGE);int v=addVoxel(1,BASE_EYE_HEIGHT,3,false,true,GRAY,0);players[0].yaw=180;acquire_tether_hit(0,v,(Vector3){0,0,1});}
    if(which>=8) {
        for(int i=5;i<activePlayers;i++)players[i].respawn_timer=9999;
        players[4].pos=(Vector3){0,BASE_EYE_HEIGHT,12};players[4].yaw=which==8?180:0;
        players[0].pos=(Vector3){which==9?1.2f:0,BASE_EYE_HEIGHT,0};players[0].yaw=180;
        if(which==10)enemy_update_equipment(4);
        if(which==12)add_static_box_at_grid(-7,7,0,8,12,12,GRAY);
        if(which==11 || which>=13) {
            firefightHumans=which==11?3:2;
            for(int i=1;i<firefightHumans;i++) {
                players[i]=players[0];players[i].pos=(Vector3){i==1?-5:5,BASE_EYE_HEIGHT,0};
                playerInput[i]=INPUT_TYPE_KEYBOARD;combatantTeam[i]=0;combatantLife[i]=LIFE_ALIVE;
            }
            players[0].pos=(Vector3){25,BASE_EYE_HEIGHT,20};players[0].yaw=90;
        }
    }
    init_static_hash();rebuild_voxel_hash();meshDirty=true;
}
static void enemy_lab_ai_step(float dt) {
    if(!enemyLab || enemyLabFixture<8 || (enemyLabPaused && !enemyLabStep))return;
    enemyLabAiTime+=dt;
    if(enemyLabFixture==11 && enemyLabAiTime>=enemyLabAiNext) {
        int source=1+((int)(enemyLabAiNext/2)%2);enemyLabAiNext=enemyLabAiTime+2;
        players[4].matter=players[4].matterMax;players[4].isExposed=false;
        apply_matter_damage(4,source,1);
    }
    if(enemyLabFixture>=13) {
        players[1].pos=(Vector3){0,BASE_EYE_HEIGHT,-8};players[1].vel=(Vector3){0};
        if(enemyLabFixture==14) {
            float phase=fmodf(enemyLabAiTime,4);
            players[1].pos.x=phase<2?-6+phase*6:18-phase*6;
            players[1].vel.x=phase<2?6:-6;
        }
    }
}
static bool enemy_lab_controls(void) {
    enemyLabStep=false;
    if(equipmentLab){if(IsKeyPressed(KEY_F4))eq_lab();if(IsKeyPressed(KEY_F7))enemyLabGod=!enemyLabGod;if(IsKeyPressed(KEY_ESCAPE))enemy_lab_exit();return true;}
    if(IsKeyPressed(KEY_F3)||couch_button(0,GAMEPAD_BUTTON_LEFT_FACE_RIGHT)) {
        int next=enemyLabCase+1;if(next==ENEMY_TYPE_SCRAPPER)next++;if(next>ENEMY_TYPE_COUNT) next=ENEMY_TYPE_GOLIATH;
        enemy_lab_start(next);
    }
    if(IsKeyPressed(KEY_F4)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_UP)) (enemyLabFixture>=0?enemy_lab_fixture(enemyLabFixture):enemy_lab_start(enemyLabCase));
    if(IsKeyPressed(KEY_F8)||couch_button(0,GAMEPAD_BUTTON_LEFT_FACE_DOWN))enemy_lab_fixture(enemyLabFixture+1);
    if(IsKeyPressed(KEY_F5)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_LEFT)) enemyLabPaused=!enemyLabPaused;
    if(IsKeyPressed(KEY_F6)||couch_button(0,GAMEPAD_BUTTON_RIGHT_TRIGGER_1)) {enemyLabPaused=true;enemyLabStep=true;}
    if(IsKeyPressed(KEY_F7)||couch_button(0,GAMEPAD_BUTTON_MIDDLE_LEFT)) enemyLabGod=!enemyLabGod;
    if(enemyLabAll && enemyLabClock>=30) {int n=enemyLabCase+1;if(n==ENEMY_TYPE_SCRAPPER)n++;enemy_lab_start(n>=ENEMY_TYPE_COUNT?ENEMY_TYPE_GOLIATH:n);}
    if(IsKeyPressed(KEY_ESCAPE)) enemy_lab_exit();
    return true;
}
static void enemy_lab_hud(int w,int h) {
    int fs=w<900?14:18;
    if(equipmentLab){DrawText("WEAPON LAB | F4 reset | F7 god | Esc exit",16,75,fs,GOLD);DrawText("Cycle grenades: D-pad Down / Tab / ;  |  Console: give glob|glock|dmr|gravity|igloo",16,99,w<900?12:16,SKYBLUE);return;}
    DrawRectangle(8,70,w-16,232,Fade(BLACK,.8f));
    DrawText(TextFormat("ENEMY LAB: %s | %s | Invulnerable: %s",enemyLabCase==ENEMY_TYPE_COUNT?"mixed":enemy_names[enemyLabCase],enemyLabPaused?"PAUSED":"RUNNING",enemyLabGod?"ON":"OFF"),16,78,fs,GOLD);
    DrawText("F3/D-right: next | F4/Top: reset | F5/Left face: pause",16,102,fs,RAYWHITE);
    DrawText("F6/RB: step | F7/Back button: invulnerability | Esc: exit",16,126,fs,RAYWHITE);
    if(enemyLabCase<ENEMY_TYPE_COUNT) DrawText(enemyLabFixture>=0?enemyLabFixtureNames[enemyLabFixture]:enemy_description(enemyLabCase),16,150,fs,RAYWHITE);
    EnemySpecialState *s=&enemySpecial[4];
    DrawText(TextFormat("State %d | Target P%d | Cooldown %.1f | Wind-up %.1f | Equipment %d",s->state,s->target+1,fmaxf(0,s->cooldown),s->windup,s->equipmentCount),16,174,fs,SKYBLUE);
    DrawText(TextFormat("%s | work %d | actions %d | resource %llu",s->reason?s->reason:"seeking",s->work,s->actions,(unsigned long long)s->selected),16,200,fs,SKYBLUE);
    DrawText("F8 / D-pad Down: combat fixtures | C / D-pad Up: grenade",16,224,fs,ORANGE);
    AiCombat *a=&aiCombat[4];const char *alerts[]={"none","seen","heard","near miss","hit"};
    DrawText(TextFormat("AI %s | P%d | %s | memory %.1f | aim %.1f deg",alerts[a->source],a->target+1,a->visible?"visible":"remembered",a->memory,a->error),16,250,fs,LIME);
    DrawText(TextFormat("Burst %d | reaction %.2f | recovery %.2f",a->burstShots,a->reaction,a->recovery),16,276,fs,LIME);
    (void)h;
}
#define ENEMY_CHECK(condition,message) do {if(!(condition)) {fprintf(stderr,"enemy-smoke FAIL: %s\n",message);return false;}} while(0)
static void enemy_lab_capture(int type) {
    if(!DirectoryExists("artifacts")) MakeDirectory("artifacts");
    RenderTexture2D shot=LoadRenderTexture(960,720);
    Camera3D cam={.position={12,12,21},.target={0,1.5f,9},.up={0,1,0},.fovy=50,.projection=CAMERA_PERSPECTIVE};
    BeginTextureMode(shot);ClearBackground((Color){24,28,36,255});BeginMode3D(cam);
    DrawPlane((Vector3){0,0,8},(Vector2){30,30},(Color){65,70,76,255});
    for(int v=0;v<voxel_count;v++) if(!clarity_shield_voxel(v)) DrawCube(voxels[v].pos,.5f,.5f,.5f,voxels[v].color);
    for(int i=4;i<activePlayers;i++) if(players[i].respawn_timer<=0) {
        DrawCube(players[i].pos,.5f,.5f,.5f,enemy_color(players[i].enemyType));
        DrawCubeWires(players[i].pos,.5f,.5f,.5f,WHITE);
    }
    enemy_draw();clarity_draw();EndMode3D();DrawText(enemy_names[type],24,24,28,RAYWHITE);
    DrawText(enemy_description(type),24,660,18,GOLD);EndTextureMode();
    Image im=LoadImageFromTexture(shot.texture);ImageFlipVertical(&im);
    ExportImage(im,TextFormat("artifacts/enemy-%s.png",enemy_names[type]));UnloadImage(im);UnloadRenderTexture(shot);
}
static bool enemy_smoke_case(int type) {
    enemy_lab_start(type);enemyLabPaused=false;int i=4;
    enemy_lab_capture(type);
    EnemySpecialState *s=&enemySpecial[i];
    if(type==ENEMY_TYPE_GOLIATH) {
        ENEMY_CHECK(goliathStates[i].count==GOLIATH_MAX_ATTACHED,"Goliath has meaningful armor");
        int part=enemy_voxel(goliathStates[i].voxels[0].voxelIdentity);detach_goliath_armor_voxel(part);
        players[i].matter=100;ENEMY_CHECK(goliath_try_regen_armor(i),"Goliath rebuilds armor");
    } else if(type==ENEMY_TYPE_SHIELDBEARER) {
        rebuild_voxel_hash();
        int front=first_voxel_hit((Ray){{0,BASE_EYE_HEIGHT,6},{0,0,1}},4,-1);
        ENEMY_CHECK(front>=0 && enemy_equipment_owner(front)==i,"shield intercepts frontal shot path");
        ENEMY_CHECK(first_voxel_hit((Ray){{3,BASE_EYE_HEIGHT,10},{-1,0,0}},2.5f,-1)<0,"shield leaves side vulnerable");
        float before=players[i].yaw;players[0].pos.x=10;enemy_update(i,.1f);
        ENEMY_CHECK(fabsf(angle_diff_deg(before,players[i].yaw))<=12.01f,"shield turn-rate cap");
        enemy_detach(front);ENEMY_CHECK(s->equipmentCount==19,"shield pieces detach");
    } else if(type==ENEMY_TYPE_SAPPER) {
        s->selected=enemy_sapper_target(i,0);int support=enemy_voxel(s->selected);
        ENEMY_CHECK(support>=0,"sapper finds structural support");
        players[i].pos=v_add(voxels[support].pos,(Vector3){0,0,1});
        ENEMY_CHECK(enemy_plant(i),"sapper plants charge");enemy_tick(1);
        ENEMY_CHECK(enemyCharges[i].active,"charge honors fuse");
        /* Reusing the old Sapper's actor slot must not turn its planted charge invulnerable. */
        uint64_t chargeSupportId=s->selected;
        init_firefight_bot_entity(i,ENEMY_TYPE_STANDARD,INPUT_TYPE_BOT_EASY);
        int v=enemy_voxel(enemyCharges[i].id);ENEMY_CHECK(enemy_equipment_owner(v)==i,"charge identity survives actor reuse");
        detach_goliath_armor_voxel(v);enemy_tick(.1f);
        ENEMY_CHECK(!enemyCharges[i].active && enemy_voxel(chargeSupportId)>=0,"shooting charge disarms without blast");
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
        int v=-1;uint64_t selected=enemy_find_payload(i);v=enemy_voxel(selected);
        ENEMY_CHECK(v>=0,"lobber finds a real group");
        players[i].goldTetherCharged=true;acquire_tether_hit(i,v,(Vector3){0,1,0});
        ENEMY_CHECK(players[i].tetherHolding,"lobber holds group with gold tether");
        s->aim=players[0].pos;ENEMY_CHECK(enemy_launch(i),"lobber ballistic group launch");
        EnemyLob *lob=&enemyLobs[MAX_COMBATANTS*2-1];ENEMY_CHECK(lob->count==8,"real eight-voxel payload");
        v=enemy_voxel(lob->ids[0]);ENEMY_CHECK(v>=0,"payload exists");
        enemyLabGod=false;players[0].invuln_timer=0;players[0].matter=100;players[0].isExposed=false;
        ENEMY_CHECK(enemy_lob_hit(v,0) && enemy_lob_hit(enemy_voxel(lob->ids[1]),0),"payload members share impact routing");
        ENEMY_CHECK(players[0].matter==70,"one payload cannot multiply damage across its members");enemyLabGod=true;
        tetherTag[v]=1;enemy_tick(.1f);ENEMY_CHECK(voxels[v].owner==0 && !lob->id,"catch transfers whole payload attribution");tetherTag[v]=0;
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
        enemy_tick(10.1f);ENEMY_CHECK(enemy_core_open(i),"core opening is permanent");
        for(int n=s->equipmentCount-1;n>=0;n--) if(s->equipment[n].group<4) enemy_detach(enemy_voxel(s->equipment[n].id));
        enemy_update_equipment(i);enemy_tick(11);ENEMY_CHECK(enemy_core_open(i),"all supports expose core permanently");
        for(int n=0;n<8;n++) {s->spawnTimer=0;enemy_update(i,.1f);s->windup=0;}
        ENEMY_CHECK(s->spawned==6,"boss reinforcement budget");
        s->cooldown=100;s->windup=0;
        for(int n=0;n<30;n++) {s->repair=0;enemy_update(i,.01f);}
        ENEMY_CHECK(s->repairs<=12 && enemy_core_open(i),"boss repair budget never restores destroyed supports");
    }
    s->windup=1;s->state=3;players[i].isExposed=true;enemy_update(i,.01f);
    ENEMY_CHECK(type==ENEMY_TYPE_FOUNDRY? !players[i].isExposed : s->windup==0,"boss keeps attacking; exposure interrupts specialists");
    enemy_cleanup(true);
    ENEMY_CHECK(enemySpecial[4].equipmentCount==0 && !enemyCharges[4].active,"cleanup clears equipment and charges");
    for(int n=0;n<ENEMY_STRUCTURE_MAX;n++) ENEMY_CHECK(!enemyStructures[n].count,"cleanup clears construction");
    fprintf(stderr,"enemy-smoke PASS: %s\n",enemy_names[type]);return true;
}
static void enemy_autonomous_step(float dt) {
    clarity_tick(dt);enemy_tick(dt);
    for(int i=MAX_PLAYERS;i<activePlayers;i++)if(players[i].respawn_timer<=0) {
        UpdateBot(i,dt);net_integrate_player(&players[i],dt,true);
        if(enemy_is_special(i))enemy_update_equipment(i);
    }
    prepare_tether_forces();simulate_voxel_pbd_steps(1.0f/120,2);
    recycle_dead_voxels(1.0f/60);deactivate_sleeping_voxels();
    rebuild_voxel_hash();update_projectiles(dt);handle_pbd_projectile_hits();
}
static bool enemy_autonomous_checks(void) {
    const EnemyType types[]={ENEMY_TYPE_BUILDER,ENEMY_TYPE_LOBBER};
    for(int t=0;t<2;t++) {
        enemy_lab_start(types[t]);enemyLabPaused=false;
        for(int frame=0;frame<1800 && !enemySpecial[4].actions;frame++)enemy_autonomous_step(1.0f/60);
        if(!enemySpecial[4].actions)fprintf(stderr,"autonomous %s stuck: %s pos=(%.1f %.1f %.1f) selected=%llu held=%d\n",enemy_names[types[t]],enemySpecial[4].reason,players[4].pos.x,players[4].pos.y,players[4].pos.z,(unsigned long long)enemySpecial[4].selected,players[4].tetherHolding);
        ENEMY_CHECK(enemySpecial[4].actions>0,"autonomous enemy executes its role without forced wind-up");
        enemy_lab_capture(types[t]);
    }
    enemy_lab_start(ENEMY_TYPE_LOBBER);clear_world_voxels();init_static_hash();
    for(int f=0;f<180;f++)enemy_autonomous_step(1.0f/60);
    ENEMY_CHECK(!enemySpecial[4].actions && voxel_count==0,"empty arena never generates Lobber ammunition");
    enemy_lab_start(ENEMY_TYPE_LOBBER);
    int payload=enemy_voxel(enemy_find_payload(4));ENEMY_CHECK(payload>=0,"interruption fixture has ammunition");
    players[4].goldTetherCharged=true;acquire_tether_hit(4,payload,(Vector3){0,1,0});
    enemySpecial[4].state=3;enemySpecial[4].windup=1;players[4].isExposed=true;enemy_update(4,.1f);
    ENEMY_CHECK(!players[4].tetherHolding && !enemySpecial[4].windup && !voxels[payload].goliathProjectile,"interruption releases ordinary debris");
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);
    for(int q=0;q<4;q++) {
        float angle=q*PI/2;players[0].pos=v_add(players[4].pos,(Vector3){-sinf(angle)*8,0,-cosf(angle)*8});
        for(int frame=0;frame<120;frame++)enemy_update(4,1.0f/60);
        Vector3 plate=v_sub(enemy_local(4,(Vector3){0,0,-1}),players[4].pos);
        ENEMY_CHECK(v_dot(v_norm(plate),v_norm(v_sub(players[0].pos,players[4].pos)))>.99f,"shield faces each quadrant");
    }
    enemy_lab_start(ENEMY_TYPE_FOUNDRY);EnemySpecialState *boss=&enemySpecial[4];
    float health=boss->coreHealth;apply_matter_damage(4,0,999);ENEMY_CHECK(boss->coreHealth==health,"protected core rejects damage");
    for(int n=boss->equipmentCount-1;n>=0;n--)if(boss->equipment[n].group<2)enemy_detach(enemy_voxel(boss->equipment[n].id));
    enemy_update_equipment(4);enemy_tick(20);ENEMY_CHECK(enemy_core_open(4),"two supports permanently open core");
    players[0].pos=v_add(players[4].pos,(Vector3){0,0,-2});
    for(int f=0;f<120;f++)enemy_autonomous_step(1.0f/60);
    ENEMY_CHECK(boss->actions>0,"open boss retains a close-range attack");
    int kills=players[0].kills;apply_matter_damage(4,0,999);
    ENEMY_CHECK(players[4].respawn_timer>0 && players[0].kills==kills+1,"ranged damage kills core and awards one kill");
    apply_matter_damage(4,0,999);ENEMY_CHECK(players[0].kills==kills+1,"no duplicate boss death credit");
    fprintf(stderr,"enemy-smoke PASS: autonomous roles and boss lifecycle\n");return true;
}
static bool enemy_run_smoke(const char *which) {
    SetTargetFPS(0);
    int type=enemy_parse_case(which);if(type<0) {fprintf(stderr,"Unknown enemy case: %s\n",which);return false;}
    if(!strcmp(which,"all") || type==ENEMY_TYPE_COUNT) {
        for(int t=ENEMY_TYPE_GOLIATH;t<ENEMY_TYPE_COUNT;t++) if(t!=ENEMY_TYPE_SCRAPPER && !enemy_smoke_case(t)) return false;
    } else if(!enemy_smoke_case(type)) return false;
    if(!enemy_autonomous_checks())return false;
    /* Sustained mixed simulation exercises compaction, ownership, physics and render paths. */
    for(int scenario=0;scenario<2;scenario++) {
    enemy_lab_start(scenario?ENEMY_TYPE_COUNT:ENEMY_TYPE_FOUNDRY);double samples[30];
    RenderTexture2D views[MAX_PLAYERS]={0};int count=0,w=0,h=0;
    for(int i=1;i<4;i++) {players[i]=players[0];players[i].pos.x=i*2;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=0;playerInput[i]=INPUT_TYPE_KEYBOARD;}
    firefightHumans=multiplayerActivePlayers=4;
    for(int frame=0;frame<600;frame++) {
        float dt=1.0f/30;double start=GetTime();clarity_tick(dt);enemy_tick(dt);
        for(int i=4;i<activePlayers;i++) if(players[i].respawn_timer<=0) {UpdateBot(i,dt);net_integrate_player(&players[i],dt,true);}
        prepare_tether_forces();simulate_voxel_pbd_steps(1.0f/120,4);update_projectiles(dt);handle_pbd_projectile_hits();
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
