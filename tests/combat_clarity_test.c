#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void empty_arena(void){enemyLab=false;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=false;firefightPractice=false;useCustomMap=false;currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=2;ResetGame();enemy_cleanup(true);clear_world_voxels();init_static_hash();gameState=GAME_STATE_PLAYING;enemyLabGod=false;for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;for(int i=0;i<2;i++){players[i].respawn_timer=0;players[i].invuln_timer=0;players[i].pos=(Vector3){i*1.5f,1,0};players[i].matter=100;players[i].isExposed=false;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=0;}worldFloorHalfSize=40;}
static void check_crosscut_routes(void) {
    enum {SIDE=129};static unsigned char blocked[SIDE*SIDE],seen[SIDE*SIDE];static int queue[SIDE*SIDE];
    for(int z=0;z<SIDE;z++)for(int x=0;x<SIDE;x++)blocked[z*SIDE+x]=spawn_aabb_blocked_for((Vector3){(x-64)*.5f,BASE_EYE_HEIGHT+.05f,(z-64)*.5f},PLAYER_SIZE*.5f,0);
    for(int q=0;q<4;q++) {
        memset(seen,0,sizeof(seen));Vector3 spawn=crosscut_spawn_point(q*2);
        int start=((int)roundf(spawn.z*2)+64)*SIDE+(int)roundf(spawn.x*2)+64,head=0,tail=0;
        assert(!blocked[start]);queue[tail++]=start;seen[start]=1;
        while(head<tail){int cell=queue[head++],x=cell%SIDE,z=cell/SIDE;int dx[4]={1,-1,0,0},dz[4]={0,0,1,-1};
            for(int n=0;n<4;n++){int nx=x+dx[n],nz=z+dz[n];if(nx<0||nx>=SIDE||nz<0||nz>=SIDE)continue;int next=nz*SIDE+nx;if(!blocked[next]&&!seen[next]){seen[next]=1;queue[tail++]=next;}}}
        Vector3 at=pickups[q*3+1].pos;int end=((int)roundf(at.z*2)+64)*SIDE+(int)roundf(at.x*2)+64;
        assert(seen[end]);
    }
}
static void snapshot_packet(int part,uint32_t seq,bool broken){uint8_t data[1000];NetWriter w;net_writer_init(&w,data,sizeof(data));net_write_u32(&w,seq);net_write_u8(&w,part);for(int n=0;n<16;n++){net_write_u32(&w,n==0 && part==0?123:0);net_write_u8(&w,0);for(int c=0;c<6;c++)net_write_f32(&w,c==1?2:0);net_write_f32(&w,n==0?1:0);}NetReader r;net_reader_init(&r,data,broken?w.length-1:w.length);assert(clarity_net_receive(NET_MSG_GRENADES,&r,seq));}
int main(int argc,char **argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,720,"Combat clarity regression");init_pbd_thread_pool();assert(initialize_physics_backend());
    empty_arena();players[1].pos=(Vector3){10,1,0};grenadeInventory[0]=2;
    assert(clarity_throw(0));assert(grenadeInventory[0]==1);assert(!clarity_throw(0));clarity_tick(.5f);assert(clarity_throw(0));assert(grenadeInventory[0]==0);assert(!clarity_throw(0));
    for(int f=0;f<220;f++)clarity_tick(1.0f/60);for(int n=0;n<32;n++)assert(!grenades[n].id);assert(!debugDrawParticles);
    empty_arena();clarity_explode((Vector3){0,1,0},0);assert(players[0].matter==40 && players[1].matter<100);assert(!explosionDamageContext);
    players[1].isExposed=true;players[1].matter=0;clarity_explode((Vector3){0,1,0},0);assert(combatantLife[1]==LIFE_DOWNED);assert(players[0].kills==0);
    empty_arena();players[1].pos=(Vector3){3,1,0};add_static_box_at_grid(2,2,0,5,-3,3,GRAY);init_static_hash();clarity_explode((Vector3){0,1,0},0);assert(players[1].matter==100);
    empty_arena();enemy_lab_spawn(4,ENEMY_TYPE_STANDARD,(Vector3){2,1,0});players[4].invuln_timer=0;enemyCharges[5]=(EnemyCharge){.active=true,.pos={0,1,0},.fuse=3};assert(clarity_evade(4,.1f));assert(players[4].vel.x>0);enemyCharges[5].active=false;
    clarity_explode((Vector3){1,1,0},4);assert(players[4].matter<100);
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);enemyLabPaused=false;int v=enemy_voxel(enemySpecial[4].equipment[0].id);assert(clarity_shield_melee(v,0));assert(enemySpecial[4].equipmentCount==15);assert(players[4].matter==75);
    for(int n=0;n<3;n++){v=enemy_voxel(enemySpecial[4].equipment[0].id);assert(clarity_shield_melee(v,0));}assert(enemySpecial[4].equipmentCount==0);
    /* Shield sections absorb six projectiles, independent of melee and actor Matter. */
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);uint64_t plate=enemySpecial[4].equipment[0].id;
    for(int n=0;n<5;n++){assert(clarity_shield_projectile(enemy_voxel(plate)));assert(enemySpecial[4].equipmentCount==20);}
    assert(clarity_shield_projectile(enemy_voxel(plate)));assert(enemySpecial[4].equipmentCount==15 && players[4].matter==75);
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);players[4].pos=(Vector3){0,2,4};players[4].yaw=0;enemy_update_equipment(4);
    for(int shot=0;shot<6;shot++) {
        int bullet=addVoxel(.25f,2,0,false,true,WHITE,0);voxels[bullet].isBullet=true;voxels[bullet].owner=0;voxels[bullet].vel=(Vector3){0,0,60};
        rebuild_voxel_hash();update_projectiles(.1f);assert(players[4].matter==75);
    }
    assert(enemySpecial[4].equipmentCount==15);
    /* Continuous, nonstacking healing with a finishing window and disabled-link cutoff. */
    empty_arena();enemy_lab_spawn(4,ENEMY_TYPE_ANCHOR,(Vector3){0,2,0});enemy_lab_spawn(5,ENEMY_TYPE_STANDARD,(Vector3){3,2,0});
    enemy_lab_spawn(6,ENEMY_TYPE_ANCHOR,(Vector3){0,2,3});players[5].matter=20;
    for(int n=0;n<60;n++)feedback_tick(1.0f/60);assert(fabsf(players[5].matter-28)<.01f);
    players[5].matter=0;players[5].isExposed=true;for(int n=0;n<119;n++)feedback_tick(1.0f/60);assert(players[5].isExposed);
    feedback_tick(.04f);assert(!players[5].isExposed && players[5].matter>0);
    enemySpecial[4].disabled=enemySpecial[6].disabled=8;float hp=players[5].matter;feedback_tick(1);assert(players[5].matter==hp);
    enemySpecial[4].disabled=0;players[5].matter=99;feedback_tick(1);assert(players[5].matter==100);
    /* Remains are physical, unowned, bounded and protected while tethered. */
    empty_arena();enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,(Vector3){0,2,4});feedback_body(4);feedback_tick(.01f);
    int corpse=enemy_voxel(enemyRemains[0].id);assert(corpse>=0 && voxels[corpse].simulate && voxels[corpse].owner<0);
    uint64_t corpseId=voxels[corpse].identity;
    for(int f=0;f<180;f++){feedback_tick(1.0f/60);prepare_tether_forces();simulate_voxel_pbd_steps(1.0f/120,2);recycle_dead_voxels(1.0f/60);deactivate_sleeping_voxels();}
    corpse=enemy_voxel(corpseId);assert(corpse>=0 && voxels[corpse].simulate);assert(clarity_debris_hit(corpse,0));
    acquire_tether_hit(0,corpse,(Vector3){0,0,-1});assert(players[0].tetherHolding);
    clarityClock+=46;feedback_tick(.01f);assert(enemy_voxel(corpseId)>=0);release_tether_velocity(0,NULL);feedback_tick(.01f);assert(enemy_voxel(corpseId)<0);
    for(int n=0;n<40;n++){feedback_body(4);feedback_tick(.01f);}int bodies=0;for(int n=0;n<32;n++)bodies+=enemy_voxel(enemyRemains[n].id)>=0;assert(bodies==32);
    for(int view=0;view<3;view++) {
        int width=view==0?960:view==1?480:1280,height=view==0?720:360;
        Camera3D camera={.position={1,2,3},.target={1,2,-7},.up={0,1,0},.fovy=60};
        Ray ray=GetScreenToWorldRayEx((Vector2){width*.5f,height*.5f},camera,width,height);
        assert(v_dot(ray.direction,(Vector3){0,0,-1})>.99999f);
    }
    /* Side-offset throws converge on the camera target, without sideways inherited velocity. */
    for(int gold=0;gold<2;gold++)for(int yaw=0;yaw<360;yaw+=90) {
        empty_arena();players[1].respawn_timer=9999;players[0].pos=(Vector3){0,2,0};players[0].yaw=yaw;players[0].pitch=0;
        Vector3 forward=player_forward(&players[0]),right=v_cross(forward,(Vector3){0,1,0});
        Vector3 at=v_add(players[0].pos,v_add(v_mul(forward,3),v_mul(right,1)));
        int block=addVoxel(at.x,at.y,at.z,false,true,GRAY,0);players[0].goldTetherCharged=gold;acquire_tether_hit(0,block,(Vector3){0,1,0});
        assert(players[0].tetherHolding);voxels[block].vel=(Vector3){9,0,-7};release_tether_velocity(0,NULL);
        Vector3 expected=v_norm(v_sub(v_add(players[0].pos,v_mul(forward,100)),at));assert(v_dot(v_norm(voxels[block].vel),expected)>.99996f);
    }
    /* Flight checks use actual glued particles and both available solvers. */
    for(int scenario=0;scenario<3;scenario++) {
        empty_arena();for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;
        players[4].respawn_timer=0;players[4].enemyType=ENEMY_TYPE_LOBBER;players[4].pos=(Vector3){0,1,-4};enemy_init(4);
        enemy_scrap_bundle((Vector3){0,5,0});int root=0;players[4].goldTetherCharged=true;acquire_tether_hit(4,root,(Vector3){0,1,0});
        Vector3 target={.5f,scenario==2?4:1,scenario==0?10:20};enemySpecial[4].aim=target;assert(enemy_launch(4));
        float previousY=voxels[0].pos.y,best=100;bool crossed=false;
        for(int f=0;f<300;f++){prepare_tether_forces();simulate_voxel_pbd_steps(PBD_MAX_STEP_DT,1);Vector3 center={0};int count=0;for(int n=0;n<voxel_count;n++)if(voxels[n].simulate){center=v_add(center,voxels[n].pos);count++;}center=v_mul(center,1.0f/count);
            if(previousY>target.y && center.y<=target.y){best=v_length(v_sub(center,target));crossed=true;break;}previousY=center.y;}
        fprintf(stderr,"Lob flight %d miss %.3f\n",scenario,best);assert(crossed && best<1.5f);
    }
    empty_arena();players[4].enemyType=ENEMY_TYPE_STANDARD;players[4].respawn_timer=0;players[4].invuln_timer=0;players[4].pos=(Vector3){5,1,0};players[4].isExposed=true;combatantTeam[4]=1;kill_player(4,0,true,false);assert(killConfirm[0]>0);int effects=0;for(int n=0;n<64;n++)effects+=clarityEffects[n].id!=0;assert(effects==1);kill_player(4,0,true,false);int after=0;for(int n=0;n<64;n++)after+=clarityEffects[n].id!=0;assert(after==effects);
    clarity_net_reset();memset(grenades,0,sizeof(grenades));snapshot_packet(1,10,false);assert(!grenades[0].id);snapshot_packet(0,10,true);assert(!grenades[0].id);snapshot_packet(0,10,false);assert(grenades[0].id==123);snapshot_packet(0,9,false);assert(clarityReceiveSequence==10);
    empty_arena();
    GrenadeState bounce={.id=1,.owner=0,.pos={0,.16f,0},.vel={0,-5,0},.fuse=-1};
    clarity_grenade_step(&bounce,.02f);assert(bounce.vel.y>0 && fabsf(bounce.fuse-1.48f)<.0001f);
    /* Further floor contacts cannot extend the original deadline. */
    bounce.pos.y=.16f;bounce.vel.y=-5;clarity_grenade_step(&bounce,.02f);
    assert(fabsf(bounce.fuse-1.46f)<.0001f);
    for(int n=0;n<145;n++)clarity_grenade_step(&bounce,.01f);
    assert(bounce.fuse>0);clarity_grenade_step(&bounce,.011f);assert(bounce.fuse==0);
    /* No airburst at the old three-second deadline. */
    GrenadeState late={.id=2,.pos={0,1000,0},.fuse=-1};
    for(int n=0;n<400;n++)clarity_grenade_step(&late,.01f);
    assert(late.fuse==-1 && late.age==1.5f);
    late.pos.y=.16f;late.vel.y=-5;clarity_grenade_step(&late,.02f);
    assert(late.fuse>.48f && late.fuse<.5f);
    for(int n=0;n<47;n++)clarity_grenade_step(&late,.01f);
    assert(late.fuse>0);clarity_grenade_step(&late,.03f);assert(late.fuse==0);
    empty_arena();players[0].invuln_timer=2;clarity_explode(players[0].pos,1);assert(players[0].matter==100);
    empty_arena();grenadeInventory[0]=0;combatantLife[0]=LIFE_DOWNED;ff_restore_human(0,30);assert(grenadeInventory[0]==0);combatantLife[0]=LIFE_WAITING;ff_restore_human(0,100);assert(grenadeInventory[0]==1);
    empty_arena();enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,(Vector3){4,1,4});firefightWaveState=FIREFIGHT_STATE_ACTIVE;memset(ffRemaining,0,sizeof(ffRemaining));for(int n=0;n<81;n++)clarity_tick(.1f);assert(swarmHintActive[4]);firefightWaveState=FIREFIGHT_STATE_INTERMISSION;clarity_tick(.1f);assert(!swarmHintActive[4]);
    gameMode=GAME_MODE_DEATHMATCH;currentWorldType=WORLD_TYPE_CROSSCUT;firefightFoundry=false;enemyLab=false;ResetGame();int supplyCount=0;for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)supplyCount++;fprintf(stderr,"Crosscut supplies %d\n",supplyCount);assert(supplyCount==13);for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)assert(clarity_supply_clear(pickups[n].pos));
    for(int q=0;q<4;q++)for(int k=0;k<3;k++){Vector3 expected=crosscut_rotate(pickups[k].pos,q);assert(v_length(v_sub(pickups[q*3+k].pos,expected))<.01f);assert(pickups[q*3+k].pos.y<1.4f);}
    check_crosscut_routes();
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;enemyLab=false;ResetGame();int types[PICKUP_TYPE_COUNT]={0};for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)types[pickups[n].type]++;assert(types[PICKUP_GRENADE]==2 && types[PICKUP_HEALTH]==2 && types[PICKUP_GOLD_TETHER]==2 && types[PICKUP_DYNAMIC_SHOT]==2);
    ff_wave_complete();int validSupplies=0;for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active){assert(pickups[n].pos.y>=0);validSupplies++;}assert(validSupplies==10);
    character_renderer_init();gameState=GAME_STATE_MENU;BeginDrawing();couch_title("CONSTRUCTOR BRAWLER");clarity_menu();couch_row(0,"Local Multiplayer",true);couch_row(1,"Firefight",false);couch_row(2,"Enemy Lab",false);rlDrawRenderBatchActive();TakeScreenshot("artifacts/clarity-menu.png");EndDrawing();
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);enemy_lab_capture(ENEMY_TYPE_SHIELDBEARER);
    for(int fixture=0;fixture<8;fixture++){enemy_lab_fixture(fixture);assert(enemyLabFixture==fixture && grenadeInventory[0]==2);}
    enemy_lab_fixture(3);enemyLab=false;gameMode=GAME_MODE_DEATHMATCH;firefightHumans=multiplayerActivePlayers=4;
    for(int i=1;i<4;i++){players[i]=players[0];players[i].pos.x=i*2;playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;}
    RenderTexture2D views[4]={0};int viewCount=0,vw=0,vh=0;render_gameplay_view(views,&viewCount,&vw,&vh,false);assert(viewCount==4);
    Image quarter=LoadImageFromTexture(views[0].texture);ImageFlipVertical(&quarter);ExportImage(quarter,"artifacts/clarity-quarter.png");UnloadImage(quarter);for(int i=0;i<viewCount;i++)UnloadRenderTexture(views[i]);
    empty_arena();ClarityEffect blast={.id=777,.kind=2,.pos={0,1,0},.age=.18f};
    Camera3D blastCam={.position={8,5,8},.target={0,1,0},.up={0,1,0},.fovy=60};
    players[0].pos=blastCam.position;drawViewPlayerIndex=0;
    for(int fallback=0;fallback<2;fallback++) {
        if(fallback){feedback_shutdown();blastShaderAttempted=true;}
        BeginDrawing();ClearBackground((Color){24,28,32,255});BeginMode3D(blastCam);DrawPlane((Vector3){0},(Vector2){20,20},DARKGRAY);DrawCube((Vector3){-2,1,0},1,2,4,BLUE);feedback_blast(&blast);EndMode3D();rlDrawRenderBatchActive();TakeScreenshot(fallback?"artifacts/blast-fallback.png":"artifacts/blast-shader.png");EndDrawing();
    }
    feedback_shutdown();
    empty_arena();init_sfx();if(sfxReady){play_sfx_at(SFX_SHIELD,(Vector3){0,1,2},4,3,1);bool voice=false;for(int n=0;n<32;n++)voice|=clarityVoices[n].sound.frameCount>0;assert(voice);clarity_reset();shutdown_sfx();for(int n=0;n<32;n++)assert(!clarityVoices[n].sound.frameCount);}
    for(int n=0;n<3;n++){Wave step=make_footstep_wave(n);assert(step.frameCount==3969);ExportWave(step,TextFormat("artifacts/footstep-%d.wav",n+1));UnloadWave(step);}
    assert(enemy_parse_case("scrapper")==ENEMY_TYPE_GOLIATH);assert(!debugDrawParticles);
    enemy_cleanup(true);character_renderer_shutdown();clarity_menu_shutdown();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Combat clarity checks passed");return 0;}
