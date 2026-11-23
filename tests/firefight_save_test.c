#define FIREFIGHT_TESTING
#include "raylib.h"
static bool keys[512],buttons[32],pad;
static bool test_pressed(int key){return key>=0 && key<512 && keys[key];}
static bool test_pad(int i){return pad && i==0;}
static bool test_button(int i,int b){return test_pad(i) && b>=0 && b<32 && buttons[b];}
#define IsKeyPressed test_pressed
#define IsGamepadAvailable test_pad
#define IsGamepadButtonPressed test_button
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static const char *savePath=".build/firefight-save-test.sav";
static void input_clear(void){memset(keys,0,sizeof(keys));memset(buttons,0,sizeof(buttons));}
static void new_run(void){
    gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_MENU;netTransport.role=NET_ROLE_OFFLINE;
    firefightPractice=enemyLab=equipmentLab=false;firefightFoundry=false;firefightEndless=true;firefightCheckpoints=false;
    currentWorldType=WORLD_TYPE_TEST;useCustomMap=false;multiplayerActivePlayers=2;
    multiplayerPlayerInput[0]=multiplayerPlayerInput[1]=INPUT_TYPE_KEYBOARD;
    fp_start_new_game();gameState=GAME_STATE_PLAYING;
    assert(fpSaveEnabled && !fpSaveError && fpSaveWave==1);
}
static void prepare_saved_wave(void){
    new_run();
    assert(fp_reserve(3));fpCount=3;fpNextId=4;
    fpItems[0]=(FfItem){1,PICKUP_GLOB,0,0,{0,1,0}};
    fpItems[1]=(FfItem){2,PICKUP_REINFORCED,-1,-1,{3,1,0}};
    fpItems[2]=(FfItem){3,PICKUP_ADHESIVE,0,-1,{0,1,0}};
    equipment[0].gun[0]=EQ_GLOB;equipment[0].gun[1]=EQ_GLOCK;equipment[0].selected=1;equipment[0].upgrades=EQ_SHIELD|EQ_SWORD|EQ_ADHESIVE;
    players[0].kills=7;players[0].deaths=2;players[0].debrisKills=3;players[0].goldTetherCharged=true;
    grenadeInventory[0]=2;exGrenades[0][0]=1;exSelectedGrenade[0]=GRENADE_STICKY;reviveCount[0]=4;
    firefightScore=12345;firefightKillsTotal=39;firefightTeamLives=5;teamUpgradeLevels[1]=2;
    start_firefight_wave(7);assert(fpSaveWave==7 && !fpSaveError);
}
static void verify_saved_wave(void){
    fp_shutdown();fpSaveWave=-1;
    assert(fp_save_available() && fpSaveWave==7 && fpSaveHumans==2);
    gameMode=GAME_MODE_DEATHMATCH;firefightPractice=true;currentWorldType=WORLD_TYPE_BLOOD;
    multiplayerActivePlayers=1;firefightScore=0;firefightCheckpoints=true;
    assert(fp_save_continue());
    assert(gameState==GAME_STATE_PLAYING && firefightWave==7 && firefightRound==1 && !firefightCheckpoints && firefightEndless);
    assert(currentWorldType==WORLD_TYPE_TEST && !firefightFoundry && firefightHumans==2 && multiplayerActivePlayers==2);
    assert(firefightScore==12345 && firefightKillsTotal==39 && firefightTeamLives==5 && teamUpgradeLevels[1]==2);
    assert(players[0].shields==75 && players[0].health==25);
    assert(players[0].kills==7 && players[0].deaths==2 && players[0].debrisKills==3 && players[0].goldTetherCharged);
    assert(fabsf(players[0].kd_ratio-8.f/3)<.001f && reviveCount[0]==4);
    assert(eq_gun(0)==EQ_GLOCK && equipment[0].gun[0]==EQ_GLOB && eq_upgrade(0,EQ_SHIELD) && eq_upgrade(0,EQ_SWORD));
    assert(grenadeInventory[0]==2 && exGrenades[0][0]==1 && exSelectedGrenade[0]==GRENADE_STICKY);
    assert(fpCount==3 && fpItems[0].owner==0 && fpItems[1].owner==-1 && fpNextId==4);
    assert(fpItems[2].type==PICKUP_ADHESIVE && fpItems[2].owner==0 && eq_upgrade(0,EQ_ADHESIVE));
    assert(!memcmp(fpOpeningRoster,fpCheckpoint.remaining,sizeof(fpOpeningRoster)));
    assert(!memcmp(fpOpeningCarriers,fpCheckpoint.pending,fpPendingCount*sizeof(int)));
}
static void malformed_files(void){
    FpSavedRun s;FfItem *items;assert(fp_save_read(&s,&items));free(items);
    /* A schema change, truncation, checksum failure and invalid counts must all fail safely. */
    FILE *f=fopen(savePath,"wb");assert(f);fwrite("FPSFFSV1",1,8,f);fclose(f);
    assert(!fp_save_continue() && firefightWave==7 && gameState==GAME_STATE_PLAYING);
    fp_save_write();f=fopen(savePath,"r+b");assert(f);fseek(f,12,SEEK_SET);uint32_t version=999;fwrite(&version,4,1,f);fclose(f);
    assert(!fp_save_continue() && firefightWave==7);
    fp_save_write();f=fopen(savePath,"r+b");assert(f);fseek(f,-1,SEEK_END);fputc(0xff,f);fclose(f);
    assert(!fp_save_continue() && firefightWave==7);
    fp_save_write();s.count=FP_SAVE_MAX_ITEMS+1;uint32_t hash=fp_save_hash(2166136261u,&s,sizeof(s));
    f=fopen(savePath,"wb");assert(f);fwrite("FPSFFSV1",1,8,f);fwrite(&hash,4,1,f);fwrite(&s,sizeof(s),1,f);fclose(f);
    assert(!fp_save_continue() && firefightWave==7);
    fp_save_write();
    puts("Invalid saves rejected without changing the running match");
}
static void menu(void){
    gameState=GAME_STATE_MENU;couchAdvanced=couchLobby=couchSettings=false;
    couch_firefight_options();couch_frame();
    RenderTexture2D target=LoadRenderTexture(960,540);
    BeginTextureMode(target);couch_draw_firefight_options(true);EndTextureMode();
    Image preview=LoadImageFromTexture(target.texture);ImageFlipVertical(&preview);
    assert(ExportImage(preview,".build/firefight-continue.png"));UnloadImage(preview);UnloadRenderTexture(target);
    keys[KEY_ENTER]=true;couch_frame();input_clear();
    assert(gameState==GAME_STATE_PLAYING && firefightWave==7 && !couchFirefightMenu);
    gameState=GAME_STATE_PAUSED;couchSelection=2;keys[KEY_ENTER]=true;couch_frame();input_clear();
    assert(gameState==GAME_STATE_MENU && couchFirefightMenu);
    keys[KEY_DOWN]=true;couch_frame();input_clear();keys[KEY_ENTER]=true;couch_frame();input_clear();
    assert(couchLobby && !couchFirefightMenu && fpSaveWave==7); // choosing setup does not erase the save
    couchLobby=false;couch_firefight_options();pad=true;buttons[GAMEPAD_BUTTON_RIGHT_FACE_DOWN]=true;couch_frame();input_clear();pad=false;
    assert(gameState==GAME_STATE_PLAYING && firefightWave==7);
    puts("Keyboard/controller Continue, New Game setup and pause-return flow passed");
}
int main(int argc,char **argv){
#ifdef _WIN32
    assert(_putenv_s("FPS_FIREFIGHT_SAVE",savePath)==0);
#else
    assert(setenv("FPS_FIREFIGHT_SAVE",savePath,1)==0);
#endif
    SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Firefight save regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    if(argc>1 && !strcmp(argv[1],"--read")){verify_saved_wave();puts("Cross-process Continue passed");}
    else {
        remove(savePath);fpSaveWave=-1;assert(!fp_save_available());
        prepare_saved_wave();verify_saved_wave();malformed_files();menu();
        /* A failed replacement leaves the previous complete save available. */
        assert(MakeDirectory(".build/firefight-save-test.sav.tmp")==0);
        fp_save_write();assert(fpSaveError);
        FpSavedRun retained;FfItem *retainedItems;
        assert(fp_save_read(&retained,&retainedItems) && retained.wave==7);free(retainedItems);
#ifdef _WIN32
        assert(RemoveDirectoryA(".build/firefight-save-test.sav.tmp"));
#else
        assert(rmdir(".build/firefight-save-test.sav.tmp")==0);
#endif
        fp_save_write();assert(!fpSaveError);
        /* Menu previews, other modes, practice and LAN cannot overwrite the save. */
        gameState=GAME_STATE_MENU;ResetGame();assert(!fpSaveEnabled);start_firefight_wave(2);
        FpSavedRun s;FfItem *items;assert(fp_save_read(&s,&items) && s.wave==7);free(items);
        firefightPractice=true;fp_start_new_game();assert(!fpSaveEnabled);
        firefightPractice=false;netTransport.role=NET_ROLE_CLIENT;fp_start_new_game();assert(!fpSaveEnabled);netTransport.role=NET_ROLE_OFFLINE;
        assert(fp_save_read(&s,&items) && s.wave==7);free(items);
        assert(fp_save_continue());firefightCheckpoints=true;assert(fp_retry());fp_tick(3.1f);assert(fpSaveEnabled && firefightWave==7);
        start_firefight_wave(8);assert(fp_save_read(&s,&items) && s.wave==8);free(items);
        fp_start_new_game();assert(fp_save_read(&s,&items) && s.wave==1 && s.score==0 && !s.count);free(items);
        firefightEndless=false;firefightWave=FF_WAVE_COUNT;ff_wave_complete();assert(!FileExists(savePath) && !fpSaveEnabled);
        new_run();firefightTeamLives=0;for(int i=0;i<2;i++){combatantLife[i]=LIFE_SPECTATING;players[i].respawn_timer=9999;}ff_update_coop(.1f);
        assert(gameState==GAME_STATE_GAMEOVER && !FileExists(savePath));
        for(int humans=1;humans<=4;humans++) {
            new_run();multiplayerActivePlayers=humans;
            for(int i=0;i<humans;i++)multiplayerPlayerInput[i]=i<2?INPUT_TYPE_KEYBOARD:INPUT_TYPE_GAMEPAD;
            fp_start_new_game();start_firefight_wave(2);multiplayerActivePlayers=1;
            assert(fp_save_continue() && firefightHumans==humans && multiplayerActivePlayers==humans);
        }
        prepare_saved_wave(); /* Leave this for the separate --read process. */
        puts("Persistence, retries, autosave isolation, New Game and completed-run cleanup passed");
    }
    fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
