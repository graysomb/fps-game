#define FIREFIGHT_TESTING
#include "raylib.h"
static int pressed,held,pad=-1,padPressed,padHeld;
static bool key_press(int k){return k==pressed&&k!=0;}
static bool key_hold(int k){return k==held&&k!=0;}
static bool pad_available(int p){return p==pad&&pad>=0;}
static bool pad_press(int p,int b){return pad_available(p)&&b==padPressed&&b!=0;}
static bool pad_hold(int p,int b){return pad_available(p)&&b==padHeld&&b!=0;}
static float pad_axis(int p,int a){(void)p;(void)a;return 0;}
#define IsKeyPressed key_press
#define IsKeyDown key_hold
#define IsGamepadAvailable pad_available
#define IsGamepadButtonPressed pad_press
#define IsGamepadButtonDown pad_hold
#define GetGamepadAxisMovement pad_axis
#define main fps_application_main
#include "../fps_ray.c"
#undef main
bool ffTestReviveHeld[4];
static void frame(float dt){bool consumed=ffend_frame(dt);assert(consumed);pressed=padPressed=0;}
static void screenshot(const char *path){Image im=LoadImageFromScreen();bool saved=ExportImage(im,path);assert(saved);UnloadImage(im);}
static void win(void){
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_GAMEOVER;
 firefightWave=FF_WAVE_COUNT;firefightWaveState=FIREFIGHT_STATE_VICTORY;
 firefightPractice=firefightEndless=enemyLab=false;firefightHumans=4;firefightKillsTotal=123;firefightScore=456;
 for(int i=0;i<4;i++){players[i].kills=10+i;players[i].debrisKills=i;reviveCount[i]=2*i;}
}
int main(void){
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Victory flow regression");SetExitKey(KEY_NULL);
 InitAudioDevice();character_renderer_init();weapon_renderer_init();duck_renderer_init(NULL);
 BeginDrawing();ClearBackground(DARKBLUE);EndDrawing();
 clear_all_chunks();clear_world_voxels();add_static_box_at_grid(0,1,0,1,0,1,GRAY);init_static_hash();rebuild_voxel_hash();
 int count=voxel_count;Voxel before=voxels[0];int originalInputs=activePlayers;
 win();ffend_reset();ffend.runSeconds=1234;frame(.016f);assert(ffend.phase==FFEND_BANNER);
 for(int i=0;i<12;i++)frame(.25f);assert(ffend.phase==FFEND_RESULTS);frame(.016f);
 frame(.016f);assert(ffend.runSeconds==1234);screenshot("artifacts/victory-town/in-game-results.png");
 players[0].kills=999;assert(ffend.playerKills[0]==10);
 pressed=KEY_ENTER;frame(.016f);assert(ffend.phase==FFEND_EPILOGUE&&town.ready);
 assert(voxel_count==count&&memcmp(&before,&voxels[0],sizeof(before))==0&&activePlayers==originalInputs);
 for(int i=0;i<4;i++)assert(ColorToInt(townPlayerColors[i])==ColorToInt(ffend.colors[i]));
 frame(.016f);frame(.016f);screenshot("artifacts/victory-town/in-game-epilogue.png");
 held=KEY_ESCAPE;for(int i=0;i<6;i++)frame(.2f);held=0;
 assert(ffend.phase==FFEND_RESULTS&&!town.ready&&!ffend.scenes[0].id&&!IsMusicValid(ffend.music));
 assert(gameState==GAME_STATE_GAMEOVER&&firefightScore==456&&ffend.playerKills[0]==10);
 /* Controller replay, and a tap must not skip. */
 pad=0;padPressed=GAMEPAD_BUTTON_RIGHT_FACE_DOWN;frame(.016f);assert(ffend.phase==FFEND_EPILOGUE);
 padHeld=GAMEPAD_BUTTON_RIGHT_FACE_RIGHT;frame(.1f);padHeld=0;frame(.1f);assert(ffend.phase==FFEND_EPILOGUE&&ffend.skip==0);
 StopMusicStream(ffend.music);ffend.time=TOWN_MUSIC_SECONDS-.3f;
 for(int i=0;i<4;i++)frame(.2f);assert(ffend.phase==FFEND_RESULTS&&!town.ready);
 /* LAN clients watch locally; a new authoritative run cancels their playback. */
 pad=-1;couchMenuDevice=-1;netTransport.role=NET_ROLE_CLIENT;netTransport.connected=true;
 pressed=KEY_ENTER;frame(.016f);assert(ffend.phase==FFEND_EPILOGUE);
 ++fpGeneration;gameState=GAME_STATE_PLAYING;assert(!ffend_frame(.016f));assert(ffend.phase==FFEND_NONE&&!town.ready);
 win();frame(.016f);for(int i=0;i<12;i++)frame(.25f);
 ffend.selection=1;pressed=KEY_ENTER;frame(.016f);assert(gameState==GAME_STATE_GAMEOVER&&ffend.error);
 /* Authoritative packet refreshes must not restart the local movie. */
 ffend.selection=0;pressed=KEY_ENTER;frame(.016f);gameState=GAME_STATE_GAMEOVER;frame(.016f);assert(ffend.phase==FFEND_EPILOGUE);
 netTransport.connected=false;assert(!ffend_frame(.016f));assert(gameState==GAME_STATE_MENU&&!town.ready);
 netTransport.role=NET_ROLE_OFFLINE;ffend_reset();win();firefightPractice=true;assert(!ffend_frame(.016f));firefightPractice=false;firefightWaveState=FIREFIGHT_STATE_DEFEAT;assert(!ffend_frame(.016f));
 ffend_reset();clear_world_voxels();clear_all_chunks();camp_shutdown();duck_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();shutdown_world_visuals();CloseAudioDevice();CloseWindow();
 puts("Victory flow: results, replay, keyboard/controller skip, natural finish, frozen world, LAN restart/disconnect and exclusions passed");return 0;
}
