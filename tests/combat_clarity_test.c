#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void empty_arena(void){enemyLab=false;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=false;firefightPractice=false;useCustomMap=false;currentWorldType=WORLD_TYPE_TEST;multiplayerActivePlayers=2;ResetGame();enemy_cleanup(true);clear_world_voxels();init_static_hash();gameState=GAME_STATE_PLAYING;enemyLabGod=false;for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;for(int i=0;i<2;i++){players[i].respawn_timer=0;players[i].invuln_timer=0;players[i].pos=(Vector3){i*1.5f,1,0};players[i].matter=100;players[i].isExposed=false;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=0;}worldFloorHalfSize=40;}
static void snapshot_packet(int part,uint32_t seq,bool broken){uint8_t data[1000];NetWriter w;net_writer_init(&w,data,sizeof(data));net_write_u32(&w,seq);net_write_u8(&w,part);for(int n=0;n<16;n++){net_write_u32(&w,n==0 && part==0?123:0);net_write_u8(&w,0);for(int c=0;c<6;c++)net_write_f32(&w,c==1?2:0);net_write_f32(&w,n==0?2:0);}NetReader r;net_reader_init(&r,data,broken?w.length-1:w.length);assert(clarity_net_receive(NET_MSG_GRENADES,&r,seq));}
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
    empty_arena();GrenadeState bounce={.id=1,.owner=0,.pos={0,.16f,0},.vel={0,-5,0},.fuse=3};clarity_grenade_step(&bounce,.02f);assert(bounce.vel.y>0 && bounce.fuse>2.9f);
    empty_arena();players[0].invuln_timer=2;clarity_explode(players[0].pos,1);assert(players[0].matter==100);
    empty_arena();grenadeInventory[0]=0;combatantLife[0]=LIFE_DOWNED;ff_restore_human(0,30);assert(grenadeInventory[0]==0);combatantLife[0]=LIFE_WAITING;ff_restore_human(0,100);assert(grenadeInventory[0]==1);
    empty_arena();enemy_lab_spawn(4,ENEMY_TYPE_SWARMER,(Vector3){4,1,4});firefightWaveState=FIREFIGHT_STATE_ACTIVE;memset(ffRemaining,0,sizeof(ffRemaining));for(int n=0;n<81;n++)clarity_tick(.1f);assert(swarmHintActive[4]);firefightWaveState=FIREFIGHT_STATE_INTERMISSION;clarity_tick(.1f);assert(!swarmHintActive[4]);
    gameMode=GAME_MODE_DEATHMATCH;currentWorldType=WORLD_TYPE_CROSSCUT;firefightFoundry=false;enemyLab=false;ResetGame();int supplyCount=0;for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)supplyCount++;fprintf(stderr,"Crosscut supplies %d\n",supplyCount);assert(supplyCount==13);
    for(int q=0;q<4;q++)for(int k=0;k<3;k++){Vector3 expected=crosscut_rotate(pickups[k].pos,q);assert(v_length(v_sub(pickups[q*3+k].pos,expected))<.01f);assert(pickups[q*3+k].pos.y<1.4f);}
    gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;enemyLab=false;ResetGame();int types[PICKUP_TYPE_COUNT]={0};for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active)types[pickups[n].type]++;assert(types[PICKUP_GRENADE]==2 && types[PICKUP_HEALTH]==2 && types[PICKUP_GOLD_TETHER]==2 && types[PICKUP_DYNAMIC_SHOT]==2);
    ff_wave_complete();int validSupplies=0;for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active){assert(pickups[n].pos.y>=0);validSupplies++;}assert(validSupplies==10);
    character_renderer_init();gameState=GAME_STATE_MENU;BeginDrawing();couch_title("CONSTRUCTOR BRAWLER");clarity_menu();couch_row(0,"Local Multiplayer",true);couch_row(1,"Firefight",false);couch_row(2,"Enemy Lab",false);rlDrawRenderBatchActive();TakeScreenshot("artifacts/clarity-menu.png");EndDrawing();
    enemy_lab_start(ENEMY_TYPE_SHIELDBEARER);enemy_lab_capture(ENEMY_TYPE_SHIELDBEARER);
    for(int fixture=0;fixture<5;fixture++){enemy_lab_fixture(fixture);assert(enemyLabFixture==fixture && grenadeInventory[0]==2);}
    enemy_lab_fixture(3);enemyLab=false;gameMode=GAME_MODE_DEATHMATCH;firefightHumans=multiplayerActivePlayers=4;
    for(int i=1;i<4;i++){players[i]=players[0];players[i].pos.x=i*2;playerInput[i]=INPUT_TYPE_KEYBOARD;combatantLife[i]=LIFE_ALIVE;}
    RenderTexture2D views[4]={0};int viewCount=0,vw=0,vh=0;render_gameplay_view(views,&viewCount,&vw,&vh,false);assert(viewCount==4);
    Image quarter=LoadImageFromTexture(views[0].texture);ImageFlipVertical(&quarter);ExportImage(quarter,"artifacts/clarity-quarter.png");UnloadImage(quarter);for(int i=0;i<viewCount;i++)UnloadRenderTexture(views[i]);
    empty_arena();init_sfx();if(sfxReady){play_sfx_at(SFX_SHIELD,(Vector3){0,1,2},4,3,1);bool voice=false;for(int n=0;n<32;n++)voice|=clarityVoices[n].sound.frameCount>0;assert(voice);clarity_reset();shutdown_sfx();for(int n=0;n<32;n++)assert(!clarityVoices[n].sound.frameCount);}
    assert(enemy_parse_case("scrapper")==ENEMY_TYPE_GOLIATH);assert(!debugDrawParticles);
    enemy_cleanup(true);character_renderer_shutdown();clarity_menu_shutdown();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Combat clarity checks passed");return 0;}
