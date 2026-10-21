#define FIREFIGHT_TESTING
#include "raylib.h"
static bool testFocused=true,padAvailable,keysDown[512],keysPressed[512],padDown[32],padPressed[32];
static bool t_key(int k){return keysDown[k];}static bool t_press(int k){return keysPressed[k];}
static bool t_pad(int d){(void)d;return padAvailable;}static bool t_down(int d,int b){(void)d;return padDown[b];}static bool t_button(int d,int b){(void)d;return padPressed[b];}
#define IsKeyDown t_key
#define IsKeyPressed t_press
#define IsGamepadAvailable t_pad
#define IsGamepadButtonDown t_down
#define IsGamepadButtonPressed t_button
static bool test_focus(void){return testFocused;}
#define IsWindowFocused test_focus
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 netTransport.role=NET_ROLE_OFFLINE;gameMode=GAME_MODE_DEATHMATCH;multiplayerActivePlayers=2;currentWorldType=WORLD_TYPE_TEST;useCustomMap=false;ResetGame();enemy_cleanup(true);clear_world_voxels();clear_pickups();init_static_hash();gameState=GAME_STATE_PLAYING;enemyLab=false;creativeModeActive=false;netRequestedCreative=false;debugConsoleOpen=false;testFocused=true;activePlayers=2;
 for(int i=0;i<4;i++){playerInput[i]=INPUT_TYPE_KEYBOARD;players[i].respawn_timer=0;players[i].matter=100;players[i].isExposed=false;players[i].pos=(Vector3){i*4.f,2,0};players[i].yaw=players[i].pitch=0;players[i].meleeSwingActive=false;localSeats[i].sensitivity=1;localSeats[i].invert=false;combatantLife[i]=LIFE_ALIVE;}
 eq_reset_all();
}
static void behavior(void){fixture();NetInputCommand toggle={.pressed=NET_INPUT_AIM};for(int g=0;g<EQ_GUN_COUNT;g++){
 eq_reset(0);equipment[0].gun[0]=g;aim_command(0,&toggle);aim_tick(.075f);assert(fabsf(aimStates[0].blend-.5f)<.001f);aim_tick(.075f);assert(fabsf(aim_zoom(0)-eqGuns[g].zoom)<.001f);assert(fabsf(tanf(30*DEG2RAD)/tanf(aim_fov(0)*.5f*DEG2RAD)-eqGuns[g].zoom)<.001f);assert(aim_fov(1)==60);aim_command(0,&toggle);aim_tick(.15f);assert(aimStates[0].blend==0);}
 equipment[0].gun[0]=EQ_RAIL;players[0].matter=0;players[0].isExposed=true;aim_command(0,&toggle);aim_tick(.15f);assert(aim_zoom(0)==4);Player p={0};mouse_angles(&p,(Vector2){100,0},0,0,.016f,0);assert(fabsf(p.yaw+3.75f)<.001f);
 const uint32_t cancels[]={NET_INPUT_SWITCH,NET_INPUT_EQUIP,NET_INPUT_MELEE,NET_INPUT_BUILD,NET_INPUT_GRENADE,NET_INPUT_HOLOGRAM,NET_INPUT_AIM_CANCEL};for(int n=0;n<7;n++){aimStates[0].active=true;NetInputCommand c={.pressed=NET_INPUT_AIM|cancels[n]};aim_command(0,&c);assert(!aimStates[0].active);}
 aimStates[0].active=true;NetInputCommand c={.held=NET_INPUT_TETHER};aim_command(0,&c);assert(!aimStates[0].active);aimStates[0].active=true;c.held=NET_INPUT_REVIVE;aim_command(0,&c);assert(!aimStates[0].active);
 aimStates[0].active=true;testFocused=false;aim_focus_tick();assert(!aimStates[0].active&&aimStates[0].blend==0);testFocused=true;aim_command(0,&toggle);players[0].respawn_timer=9999;aim_tick(.016f);assert(aimStates[0].blend==0&&!aimStates[0].active);
}
static void interaction(void){fixture();NetInputCommand c={0};aim_interact(0,true,.3f,&c);assert(!aimStates[0].revive&&!c.pressed);aim_interact(0,false,.01f,&c);aim_interact(0,true,.24f,&c);aim_interact(0,false,.01f,&c);assert(c.pressed==NET_INPUT_SWITCH && !c.held);c=(NetInputCommand){0};aim_interact(0,true,.25f,&c);assert(c.held==NET_INPUT_REVIVE);c=(NetInputCommand){0};aim_interact(0,false,.01f,&c);assert(!c.pressed&&!c.held);}
static void pickup_hold(void){
 fixture();aimStates[0].tapArmed=true;pickups[0]=(Pickup){.active=true,.type=PICKUP_RAILGUN,.pos=players[0].pos};
 NetInputCommand c={0};aim_interact(0,true,.26f,&c);assert(c.pressed==NET_INPUT_EQUIP && !c.held);
 c=(NetInputCommand){0};aim_interact(0,true,.3f,&c);assert(!c.pressed && !c.held);
 aim_interact(0,false,.01f,&c);assert(!c.pressed);clear_pickups();
}
static void bindings(void){fixture();NetInputCommand c={0};playerInput[0]=INPUT_TYPE_GAMEPAD;padAvailable=true;aim_sample(0,0,.016f,&c);padDown[GAMEPAD_BUTTON_RIGHT_THUMB]=padPressed[GAMEPAD_BUTTON_RIGHT_THUMB]=true;aim_sample(0,0,.016f,&c);assert((c.pressed&NET_INPUT_AIM)&&!(c.held&NET_INPUT_REVIVE));memset(padDown,0,sizeof(padDown));memset(padPressed,0,sizeof(padPressed));c=(NetInputCommand){0};padDown[GAMEPAD_BUTTON_RIGHT_FACE_LEFT]=true;aim_sample(0,0,.26f,&c);assert(c.held&NET_INPUT_REVIVE);padDown[GAMEPAD_BUTTON_RIGHT_FACE_LEFT]=false;c=(NetInputCommand){0};aim_sample(0,0,.01f,&c);assert(!(c.pressed&NET_INPUT_EQUIP));padAvailable=false;playerInput[0]=INPUT_TYPE_KEYBOARD;c=(NetInputCommand){0};aim_sample(0,0,.016f,&c);keysPressed[KEY_LEFT_SHIFT]=keysDown[KEY_LEFT_SHIFT]=true;aim_sample(0,0,.016f,&c);assert(c.pressed&NET_INPUT_AIM);keysPressed[KEY_LEFT_SHIFT]=keysDown[KEY_LEFT_SHIFT]=false;c=(NetInputCommand){0};aim_sample(1,1,.016f,&c);keysPressed[KEY_ENTER]=keysDown[KEY_ENTER]=true;aim_sample(1,1,.016f,&c);assert(c.pressed&NET_INPUT_AIM);keysPressed[KEY_ENTER]=keysDown[KEY_ENTER]=false;}
static void receive_toggle(uint32_t sequence){uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,77,0);net_write_u8(&w,1);NetInputCommand c={.sequence=sequence,.pressed=NET_INPUT_AIM};net_write_input(&w,&c);net_on_receive(&netTransport,1,FPS_NET_CHANNEL_INPUT,packet,w.length,NULL);}
static void duplicate(void){fixture();netTransport.role=NET_ROLE_HOST;netTransport.session_id=77;netPeerPlayerMask[1]=2;netLastProcessedInput[1]=0;netHasPendingInput[1]=false;receive_toggle(10);receive_toggle(10);assert(netPendingInput[1].pressed&NET_INPUT_AIM);eq_command(1,&netPendingInput[1],.016f);assert(aimStates[1].active);netLastProcessedInput[1]=10;netHasPendingInput[1]=false;receive_toggle(10);assert(!(netPendingInput[1].pressed&NET_INPUT_AIM));receive_toggle(11);receive_toggle(12);assert(!(netPendingInput[1].pressed&NET_INPUT_AIM));netTransport.role=NET_ROLE_OFFLINE;}
static void captures(void){fixture();weaponRenderingReady=weapon_renderer_init();assert(weaponRenderingReady);character_renderer_init();const int sizes[3][2]={{960,540},{480,540},{480,270}};for(int size=0;size<3;size++)for(int gun=0;gun<EQ_GUN_COUNT;gun++){
 equipment[0].gun[0]=gun;aimStates[0].active=true;aim_tick(.15f);equipment[0].charge=.4f;int w=sizes[size][0],h=sizes[size][1];weaponViewAspect=(float)w/h;weaponFrameSampled=false;RenderTexture2D rt=LoadRenderTexture(w,h);Camera3D cam={players[0].pos,{0,2,-10},{0,1,0},aim_fov(0),CAMERA_PERSPECTIVE};BeginTextureMode(rt);ClearBackground((Color){26,36,51,255});BeginMode3D(cam);DrawPlane((Vector3){0},(Vector2){100,100},DARKGRAY);for(int n=-4;n<5;n++)DrawCube((Vector3){n*3.f,1,-20},1,2,1,n==0?ORANGE:GRAY);EndMode3D();cam.fovy=60;weapon_clear_depth();BeginMode3D(cam);float hide=gun==EQ_RAIL?1:0;weapon_view_hide(hide);character_view_hide(hide);WeaponPose pose=player_weapon_pose(0,true);weapon_draw(pose,&weaponVisuals[0],cam.position,1,false,WEAPON_COLOR);character_draw_paws(pose,BLUE,false);weapon_view_hide(0);character_view_hide(0);EndMode3D();aim_overlay(0,w,h);eq_hud(0,w,h);EndTextureMode();assert(aimShader.id && aimShader.id!=rlGetShaderIdDefault());Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/aim-%d-%d.png",gun,size));UnloadImage(im);UnloadRenderTexture(rt);}
 aim_shutdown();aimShaderTried=true;RenderTexture2D rt=LoadRenderTexture(480,270);BeginTextureMode(rt);ClearBackground(DARKGRAY);aim_overlay(0,480,270);EndTextureMode();Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,"artifacts/aim-fallback.png");UnloadImage(im);UnloadRenderTexture(rt);aim_shutdown();
}
static void unobstructed_hud(void){
 fixture();const int sizes[3][2]={{960,540},{480,540},{480,270}};
 for(int n=0;n<3;n++){
  int w=sizes[n][0],h=sizes[n][1];RenderTexture2D rt=LoadRenderTexture(w,h);Image hud[2];
  for(int mode=0;mode<2;mode++){
   aimStates[0].active=mode!=0;aim_tick(.15f);
   BeginTextureMode(rt);ClearBackground(BLANK);eq_hud(0,w,h);combat_ui_abilities(0,w,h);EndTextureMode();hud[mode]=LoadImageFromTexture(rt.texture);
  }
  Color *a=LoadImageColors(hud[0]),*b=LoadImageColors(hud[1]);assert(!memcmp(a,b,w*h*sizeof(Color)));
  UnloadImageColors(a);UnloadImageColors(b);UnloadImage(hud[0]);UnloadImage(hud[1]);
  for(int g=0;g<EQ_GUN_COUNT;g++){
   equipment[0].gun[0]=g;aimStates[0].active=true;aim_tick(.15f);
   BeginTextureMode(rt);ClearBackground(BLANK);aim_overlay(0,w,h);EndTextureMode();Image im=LoadImageFromTexture(rt.texture);Color *pixels=LoadImageColors(im);
   for(int y=0;y<h;y++)for(int x=0;x<w;x++)if(x<w/5||x>=w*4/5||y<h/5||y>=h*4/5)assert(pixels[y*w+x].a==0);
   UnloadImageColors(pixels);UnloadImage(im);
  }
  UnloadRenderTexture(rt);
 }
}
static void hud_gameplay_capture(void){
 fixture();weaponBloomEnabled=true;RenderTexture2D views[4]={0};int count=0,w=0,h=0;
 for(int seats=1;seats<=4;seats*=2){
  activePlayers=multiplayerActivePlayers=seats;
  for(int i=0;i<seats;i++)equipment[i].gun[0]=EQ_RAIL;
  for(int mode=0;mode<2;mode++){
   for(int i=0;i<seats;i++)aimStates[i].active=mode!=0;
   aim_tick(.15f);for(int frame=0;frame<3;frame++)render_gameplay_view(views,&count,&w,&h,false);
   Image frame=LoadImageFromScreen();ExportImage(frame,TextFormat("artifacts/aim-hud-%d-%d.png",seats,mode));UnloadImage(frame);
  }
 }
 for(int i=0;i<count;i++)UnloadRenderTexture(views[i]);
}
static void four_views(void){fixture();activePlayers=multiplayerActivePlayers=4;winningScore=999;for(int i=0;i<4;i++){equipment[i].gun[0]=i+2;players[i].pos=(Vector3){i*2.f,2,0};}for(int x=-4;x<=4;x++)add_static_box_at_grid(x*6,x*6+2,0,5,-25,-25,GRAY);init_static_hash();rebuild_all_voxel_surfaces();RenderTexture2D views[4]={0};int count=0,w=0,h=0;drawWeaponModels=true;
 for(int bloom=0;bloom<2;bloom++){weaponBloomEnabled=bloom!=0;
 for(int mode=0;mode<2;mode++){for(int i=0;i<4;i++)aimStates[i].active=mode!=0;aim_tick(.15f);for(int warm=0;warm<15;warm++)render_gameplay_view(views,&count,&w,&h,false);double start=GetTime();for(int n=0;n<45;n++){render_gameplay_view(views,&count,&w,&h,false);rlDrawRenderBatchActive();}assert(count==4);fprintf(stderr,"Four views aim=%d pass=%d %.3f ms/frame\n",mode,bloom,(GetTime()-start)*1000/45);if(mode){for(int i=0;i<4;i++){Image im=LoadImageFromTexture(views[i].texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/aim-four-%d.png",i));UnloadImage(im);}}}}
 for(int i=0;i<count;i++)UnloadRenderTexture(views[i]);
}
int main(int argc,char**argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Aiming checks");init_pbd_thread_pool();assert(initialize_physics_backend());behavior();interaction();pickup_hold();bindings();duplicate();captures();unobstructed_hud();hud_gameplay_capture();four_views();aim_shutdown();assert(!aimShader.id);for(int i=0;i<4;i++)weapon_bloom_unload(&weaponBlooms[i]);shutdown_world_visuals();weapon_renderer_shutdown();character_renderer_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Aiming checks passed");return 0;}
