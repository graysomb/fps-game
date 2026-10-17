#define FIREFIGHT_TESTING
#include "raylib.h"
static bool focused=true,keys[512],pressed[512],buttons[8],clicks[8];
static Vector2 delta;
static int captures,releases;
static bool key_down(int k){return keys[k];}static bool key_pressed(int k){return pressed[k];}
static bool button_down(int k){return buttons[k];}static bool button_pressed(int k){return clicks[k];}
static bool focus(void){return focused;}static Vector2 motion(void){return delta;}
static void test_capture(void){captures++;}static void release(void){releases++;}
#define IsKeyDown key_down
#define IsKeyPressed key_pressed
#define IsMouseButtonDown button_down
#define IsMouseButtonPressed button_pressed
#define IsWindowFocused focus
#define GetMouseDelta motion
#define DisableCursor test_capture
#define EnableCursor release
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void setup(void){
    netTransport.role=NET_ROLE_OFFLINE;gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;
    ResetGame();enemy_cleanup(true);clear_world_voxels();init_static_hash();gameState=GAME_STATE_PLAYING;
    enemyLab=debugConsoleOpen=false;focused=true;mouseCaptured=false;mouseSlot=mouseSeat=-1;mouseTetherHeld=false;
    memset(keys,0,sizeof(keys));memset(pressed,0,sizeof(pressed));memset(buttons,0,sizeof(buttons));memset(clicks,0,sizeof(clicks));delta=(Vector2){0};
    for(int i=0;i<2;i++){playerInput[i]=INPUT_TYPE_KEYBOARD;players[i].respawn_timer=0;players[i].invuln_timer=0;players[i].pos=(Vector3){0,2,-i*5};players[i].yaw=players[i].pitch=0;players[i].matter=100;players[i].isExposed=false;localSeats[i].sensitivity=1;localSeats[i].invert=false;}
}
int main(void){
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,480,"Mouse tests");init_pbd_thread_pool();
    for(int hz=30;hz<=120;hz*=2){Player p={0};for(int f=0;f<hz;f++)mouse_angles(&p,(Vector2){600.f/hz,120.f/hz},0,0,1.f/hz,0);assert(fabsf(p.yaw+90)<.001f && fabsf(p.pitch+18)<.001f);}
    setup();delta=(Vector2){900,900};mouse_frame(.016f);assert(mouseSlot==0 && players[0].yaw==0);delta=(Vector2){10,5};mouse_frame(.016f);assert(players[0].yaw==-1.5f && players[1].yaw==0);
    float yaw=players[0].yaw;delta=(Vector2){0};players[0].yaw_vel=100;mouse_frame(.016f);HandleKeyboardInput(0,.016f);assert(players[0].yaw==yaw && players[0].yaw_vel==0);
    localSeats[0].sensitivity=2;localSeats[0].invert=true;mouse_angles(&players[0],(Vector2){0,1000},0,0,.01f,0);assert(players[0].pitch==89);
    setup();playerInput[0]=INPUT_TYPE_GAMEPAD;mouse_frame(.016f);assert(mouseSeat==1 && mouseSlot==1);delta.x=10;mouse_frame(.016f);assert(players[1].yaw==-1.5f && players[0].yaw==0);
    setup();buttons[MOUSE_BUTTON_LEFT]=clicks[MOUSE_BUTTON_LEFT]=true;mouse_frame(.016f);assert(!mouseButtonsArmed);mouse_frame(.016f);assert(voxel_count==0);buttons[0]=clicks[0]=false;mouse_frame(.016f);assert(mouseButtonsArmed);
    players[0].last_shot_time=-100;buttons[0]=clicks[0]=true;mouse_frame(.016f);assert(voxel_count==1);clicks[0]=false;mouse_frame(.016f);assert(voxel_count==1);
    setup();mouse_frame(.016f);mouse_frame(.016f);buttons[MOUSE_BUTTON_RIGHT]=true;mouse_frame(.016f);assert(players[0].drainTarget==2);keys[KEY_R]=true;buttons[MOUSE_BUTTON_RIGHT]=false;mouse_frame(.016f);assert(players[0].drainTarget==2);keys[KEY_R]=false;mouse_frame(.016f);assert(!players[0].drainTarget);
    buttons[MOUSE_BUTTON_RIGHT]=true;mouse_frame(.016f);focused=false;mouse_frame(.016f);assert(!mouseCaptured && !players[0].drainTarget && releases>0);
    focused=true;delta.x=500;mouse_frame(.016f);assert(players[0].yaw==0);gameState=GAME_STATE_PAUSED;mouse_sync();assert(!mouseCaptured);gameState=GAME_STATE_PLAYING;debugConsoleOpen=true;mouse_sync();assert(!mouseCaptured);
    setup();mouse_frame(.016f);mouse_frame(.016f);players[0].last_melee_time=-100;clicks[MOUSE_BUTTON_MIDDLE]=true;mouse_frame(.016f);assert(players[0].meleeSwingActive);
    clicks[MOUSE_BUTTON_MIDDLE]=false;players[0].respawn_timer=2;mouse_frame(.016f);assert(mouseSlot==0);
    playerInput[0]=INPUT_TYPE_GAMEPAD;delta.x=500;mouse_frame(.016f);assert(mouseSlot==1 && players[1].yaw==0);
    setup();gameState=GAME_STATE_CREATIVE;mouse_frame(.016f);mouse_frame(.016f);delta.x=10;clicks[0]=true;players[0].last_shot_time=-100;mouse_frame(.016f);assert(players[0].yaw==-1.5f && voxel_count==0);
    setup();netTransport.role=NET_ROLE_CLIENT;netLocalPlayerCount=2;netLocalPlayerSlots[0]=2;netLocalPlayerSlots[1]=3;mouse_frame(.016f);assert(mouseSlot==2 && mouseSeat==0);
    mouseFrameDelta=(Vector2){20,-10};players[2].yaw=0;players[2].pitch=0;players[2].respawn_timer=0;combatantLife[2]=LIFE_ALIVE;NetInputCommand c={0};mouse_command(&c,0,2,.016f);assert(c.mouse_look && c.aim_yaw==-3 && c.aim_pitch==1.5f);
    Player authoritative={0};c.look_x=32767;net_apply_movement_command(&authoritative,&c,.016f);assert(authoritative.yaw==-3);net_apply_movement_command(&authoritative,&c,.016f);assert(authoritative.yaw==-3);net_apply_creative_command(&authoritative,&c,.016f);assert(authoritative.pitch==1.5f);
    uint8_t data[128];NetWriter w;NetReader r;NetInputCommand decoded;net_writer_init(&w,data,sizeof(data));net_write_input(&w,&c);net_reader_init(&r,data,w.length);assert(net_read_input(&r,&decoded) && decoded.aim_yaw==c.aim_yaw);
    c.aim_yaw=NAN;net_writer_init(&w,data,sizeof(data));net_write_input(&w,&c);net_reader_init(&r,data,w.length);assert(!net_read_input(&r,&decoded));
    c.aim_yaw=-30;c.look_x=0;c.sequence=100;
    netTransport.role=NET_ROLE_HOST;netTransport.session_id=77;netPeerPlayerMask[1]=2;netLastProcessedInput[1]=0;netHasPendingInput[1]=false;
    for(int n=0;n<2;n++){
        c.sequence=99+n;c.aim_yaw=n?-30:-15;
        net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_INPUT,0,77,0);net_write_u8(&w,1);net_write_input(&w,&c);
        net_on_receive(&netTransport,1,FPS_NET_CHANNEL_INPUT,data,w.length,NULL);
    }
    assert(netPendingInput[1].sequence==100 && netPendingInput[1].aim_yaw==-30);
    netPredictionHead[1]=0;netPredictionCount[1]=2;
    netPredictionHistory[1][0]=(NetPredictedCommand){.command=c,.dt=.01f};
    c.sequence=101;c.aim_yaw=-45;netPredictionHistory[1][1]=(NetPredictedCommand){.command=c,.dt=.01f};
    NetPlayerWireState state={.acknowledged_input=100,.pos={0,2,0},.yaw=-30};
    net_reconcile_local(1,&state);assert(players[1].yaw==-45 && netPredictionCount[1]==1);
    setup();mouse_frame(.016f);mouse_frame(.016f);players[1].pos.z=-30;
    int voxel=addVoxel(0,2,-3,false,true,WHITE,0);rebuild_voxel_hash();buttons[MOUSE_BUTTON_RIGHT]=true;mouse_frame(.016f);assert(players[0].tetherHolding);
    focused=false;mouse_frame(.016f);assert(!players[0].tetherHolding && v_length(voxels[voxel].vel)<.001f);
    netTransport.role=NET_ROLE_OFFLINE;shutdown_pbd_thread_pool();CloseWindow();puts("Mouse input checks passed");return 0;
}

