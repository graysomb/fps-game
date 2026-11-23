#define FIREFIGHT_TESTING
#include "raylib.h"
#include <assert.h>
static bool tk[512],tp[512],tr[512],tm[7],tmp[7],pads[4],pb[4][32],pp[4][32],focused=true;
static float axes[4][8],wheel;static Vector2 motion,position;
static bool keyDown(int k){return k>=0&&k<512&&tk[k];}static bool keyPress(int k){return k>=0&&k<512&&tp[k];}static bool keyRelease(int k){return k>=0&&k<512&&tr[k];}
static bool mouseDown(int k){return k>=0&&k<7&&tm[k];}static bool mousePress(int k){return k>=0&&k<7&&tmp[k];}
static bool padAvailable(int d){return d>=0&&d<4&&pads[d];}static bool padDown(int d,int b){return padAvailable(d)&&b>=0&&b<32&&pb[d][b];}static bool padPress(int d,int b){return padAvailable(d)&&b>=0&&b<32&&pp[d][b];}static float padAxis(int d,int a){return padAvailable(d)&&a>=0&&a<8?axes[d][a]:0;}
static bool focus(void){return focused;}static Vector2 delta(void){return motion;}static Vector2 cursor(void){return position;}static float scroll(void){return wheel;}
#define IsKeyDown keyDown
#define IsKeyPressed keyPress
#define IsKeyReleased keyRelease
#define IsMouseButtonDown mouseDown
#define IsMouseButtonPressed mousePress
#define IsGamepadAvailable padAvailable
#define IsGamepadButtonDown padDown
#define IsGamepadButtonPressed padPress
#define GetGamepadAxisMovement padAxis
#define IsWindowFocused focus
#define GetMouseDelta delta
#define GetMousePosition cursor
#define GetMouseWheelMove scroll
#define main fps_application_main
#include "../fps_ray.c"
#undef main
bool ffTestReviveHeld[4];
static void neutral(void){memset(tk,0,sizeof(tk));memset(tp,0,sizeof(tp));memset(tr,0,sizeof(tr));memset(tm,0,sizeof(tm));memset(tmp,0,sizeof(tmp));memset(pb,0,sizeof(pb));memset(pp,0,sizeof(pp));memset(axes,0,sizeof(axes));for(int d=0;d<4;d++){axes[d][GAMEPAD_AXIS_LEFT_TRIGGER]=-1;axes[d][GAMEPAD_AXIS_RIGHT_TRIGGER]=-1;}motion=(Vector2){0};wheel=0;}
static void fixture(void){
    neutral();ct_defaults(&ctLibrary);ctReady=true;ctMenuOpen=false;ctWasBlocked=false;focused=true;gameState=GAME_STATE_PLAYING;gameMode=GAME_MODE_DEATHMATCH;netTransport.role=NET_ROLE_OFFLINE;
    activePlayers=multiplayerActivePlayers=4;enemyLab=enemyLabPaused=equipmentLab=creativeModeActive=debugConsoleOpen=false;cmTesting=false;worldPeriodic=false;
    memset(players,0,sizeof(players));memset(aimStates,0,sizeof(aimStates));memset(pickups,0,sizeof(pickups));eq_reset_all();
    for(int i=0;i<4;i++){playerInput[i]=multiplayerPlayerInput[i]=i<2?INPUT_TYPE_KEYBOARD:INPUT_TYPE_GAMEPAD;players[i].pos=(Vector3){i*20,2,0};players[i].matter=players[i].matterMax=100;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;}
    pads[0]=pads[2]=false;pads[1]=pads[3]=true;localSeats[2].device=3;localSeats[3].device=1;
    mouseCaptured=mouseButtonsArmed=true;mouseSeat=mouseSlot=0;ct_cancel_all();ct_frame(.016f);ct_frame(.016f);
}
static void persistence(void){
    CtLibrary a,b;ct_defaults(&a);a.profiles[0].tuning[1].x=2.3f;a.profiles[0].bindings[0][CT_FIRE][0]=(CtBinding){1,KEY_K,3};a.profiles[0].verticalFov=85;a.assigned[3]=0;
    assert(ct_save_file("roundtrip.cfg",&a));assert(ct_read_file("roundtrip.cfg",&b));assert(!memcmp(&a,&b,sizeof(a)));
    FILE *f=fopen("bad.cfg","w");fputs("CONTROLS 1 1\nSEATS 0 0 0 0\nPROFILE 0 0\nNAME Broken\nTUNE 1 x nan\nENDPROFILE\nEND\n",f);fclose(f);CtLibrary before=b;assert(!ct_read_file("bad.cfg",&b));assert(!memcmp(&before,&b,sizeof(b)));
    assert(!ct_save_file("missing/folder/settings.cfg",&a));assert(ct_read_file("roundtrip.cfg",&b));assert(!memcmp(&a,&b,sizeof(a)));
    f=fopen("couch_settings.cfg","w");fputs("1.8 0.22 1\n1 0.15 0\n1 0.15 0\n1 0.15 0\n",f);fclose(f);remove("controls.cfg");ct_load();assert(ctLibrary.profiles[0].tuning[1].x==1.8f&&ctLibrary.profiles[0].tuning[2].inner==.22f&&ctLibrary.profiles[0].tuning[1].invertY);assert(FileExists("couch_settings.cfg")&&FileExists("controls.cfg"));
    f = fopen("legacy.cfg", "w");
    assert(f);
    fprintf(f, "CONTROLS 1 3\nSEATS 0 1 2 2\n"
               "PROFILE 0 0\nNAME Legacy\nTUNE 1 x 2\nOPTIONS 2 1 0 1 0 0 0 0\nBIND 0 GRENADE 0 1 %d 0\nENDPROFILE\n"
               "PROFILE 1 1\nBIND 0 EQUIP 0 1 %d 0\nENDPROFILE\n"
               "PROFILE 2 0\nBIND 0 GRENADE 0 1 %d 0\nENDPROFILE\nEND\n", KEY_C, KEY_Y, KEY_N);
    fclose(f);
    assert(ct_read_file("legacy.cfg", &b));
    assert(b.profiles[0].verticalFov == 60 && b.profiles[0].tuning[1].x == 2);
    assert(b.profiles[0].bindings[0][CT_GRENADE][0].code == KEY_Y);
    assert(b.profiles[1].bindings[0][CT_EQUIP][0].code == KEY_C);
    assert(b.profiles[2].bindings[0][CT_GRENADE][0].code == KEY_N);
    f = fopen("bad-fov.cfg", "w");
    assert(f);
    fputs("CONTROLS 1 1\nSEATS 0 0 0 0\nPROFILE 0 0\nNAME Invalid FOV\nFOV 95\nENDPROFILE\nEND\n", f);
    fclose(f);
    assert(!ct_read_file("bad-fov.cfg", &b));
    assert(b.profiles[0].verticalFov == 60 && b.profiles[0].tuning[1].x == 2);
    puts("Profiles: round trip, migration, invalid input and atomic save failure PASS");
}
static void actions(void){
    fixture();
    tk[KEY_Y]=true;ct_frame(.016f);
    assert(ct_command(0,0,.016f).pressed&NET_INPUT_GRENADE);
    assert(!(ct_command(1,1,.016f).pressed&NET_INPUT_EQUIP));
    assert(!strcmp(ct_hint(0,CT_GRENADE),"Y"));
    neutral();ct_frame(.016f);
    tk[KEY_C]=tk[KEY_PERIOD]=true;pb[3][GAMEPAD_BUTTON_LEFT_FACE_UP]=true;ct_frame(.016f);
    assert(!(ct_command(0,0,.016f).pressed&NET_INPUT_GRENADE));
    assert((ct_command(1,1,.016f).pressed&(NET_INPUT_EQUIP|NET_INPUT_GRENADE))==(NET_INPUT_EQUIP|NET_INPUT_GRENADE));
    assert(ct_command(2,2,.016f).pressed&NET_INPUT_GRENADE);

    fixture();CtProfile *p=ct_profile(0);p->bindings[0][CT_FIRE][0]=(CtBinding){1,KEY_K,0};tk[KEY_LEFT_CONTROL]=true;ct_frame(.016f);assert(!(ct_command(0,0,.016f).pressed&NET_INPUT_FIRE));tk[KEY_LEFT_CONTROL]=false;tk[KEY_K]=true;ct_frame(.016f);assert(ct_command(0,0,.016f).pressed&NET_INPUT_FIRE);assert(ct_command(0,0,.016f).pressed&NET_INPUT_FIRE);
    neutral();ct_frame(.016f);tm[0]=true;ct_frame(.016f);assert(ct_command(0,0,.016f).held&NET_INPUT_FIRE);assert(!(ct_command(1,1,.016f).held&NET_INPUT_FIRE));
    neutral();ct_frame(.016f);wheel=1;ct_frame(.016f);assert(ct_pressed(0,CT_SWITCH));ct_frame(.016f);assert(ct_pressed(0,CT_SWITCH));
    neutral();ct_frame(.016f);axes[3][GAMEPAD_AXIS_RIGHT_TRIGGER]=.3f;ct_frame(.016f);assert(ct_command(2,2,.016f).held&NET_INPUT_FIRE);assert(!(ct_command(3,3,.016f).held&NET_INPUT_FIRE));
    focused=false;ct_frame(.016f);assert(ct_command(2,2,.016f).pressed&NET_INPUT_TETHER_CANCEL);focused=true;ct_frame(.016f);assert(!(ct_command(2,2,.016f).held&NET_INPUT_FIRE));neutral();ct_frame(.016f);axes[3][GAMEPAD_AXIS_RIGHT_TRIGGER]=1;ct_frame(.016f);assert(ct_command(2,2,.016f).held&NET_INPUT_FIRE);
    fixture();creativeModeActive=true;gameState=GAME_STATE_CREATIVE;ct_frame(.016f);tk[KEY_LEFT_CONTROL]=tk[KEY_Z]=true;ct_frame(.016f);assert(ct_pressed(0,CT_UNDO));assert(!ct_down(0,CT_PLACE));
    fixture();tm[0]=true;ct_frame(.016f);assert(ct_down(0,CT_FIRE));players[0].respawn_timer=2;ct_frame(.016f);assert(!ct_down(0,CT_FIRE));players[0].respawn_timer=0;ct_frame(.016f);assert(!ct_down(0,CT_FIRE));tm[0]=false;ct_frame(.016f);tm[0]=true;ct_frame(.016f);assert(ct_pressed(0,CT_FIRE));
    puts("Shared actions: remaps, mouse ownership, wheel, device IDs, focus and editor chords PASS");
}
static void modes(void){
    fixture();ct_profile(0)->aimToggle[0]=0;tk[KEY_LEFT_SHIFT]=true;ct_frame(.016f);NetInputCommand c=ct_command(0,0,.016f);aim_command(0,&c);assert(aimStates[0].active);tk[KEY_LEFT_SHIFT]=false;ct_frame(.016f);c=ct_command(0,0,.016f);aim_command(0,&c);assert(!aimStates[0].active);
    fixture();ct_profile(0)->tetherToggle[1]=1;tm[1]=true;ct_frame(.016f);assert(ct_command(0,0,.016f).held&NET_INPUT_TETHER);tm[1]=false;ct_frame(.016f);assert(ct_command(0,0,.016f).held&NET_INPUT_TETHER);tm[1]=true;ct_frame(.016f);assert(!(ct_command(0,0,.016f).held&NET_INPUT_TETHER));
    fixture();CtProfile *p=ct_profile(0);p->combined[0]=1;p->bindings[0][CT_INTERACT][0]=(CtBinding){1,KEY_B,0};p->bindings[0][CT_EQUIP][0]=(CtBinding){0};ct_command(0,0,.016f);tk[KEY_B]=true;ct_frame(.1f);ct_command(0,0,.1f);tk[KEY_B]=false;ct_frame(.016f);assert(ct_command(0,0,.016f).pressed&NET_INPUT_SWITCH);
    fixture();gameMode=GAME_MODE_FIREFIGHT;firefightHumans=2;firefightWaveState=FIREFIGHT_STATE_ACTIVE;combatantTeam[0]=combatantTeam[1]=0;combatantLife[1]=LIFE_DOWNED;players[1].pos=players[1].death_pos=players[0].pos;ct_profile(0)->combined[0]=1;ct_profile(0)->bindings[0][CT_INTERACT][0]=(CtBinding){1,KEY_B,0};ct_profile(0)->bindings[0][CT_EQUIP][0]=(CtBinding){0};pickups[0]=(Pickup){.active=true,.type=PICKUP_GLOCK,.pos=players[0].pos};ct_command(0,0,.016f);tk[KEY_B]=true;for(int n=0;n<4;n++){ct_frame(.1f);c=ct_command(0,0,.1f);}assert(c.held&NET_INPUT_REVIVE);assert(!(c.pressed&NET_INPUT_EQUIP));tk[KEY_B]=false;ct_frame(.016f);assert(!(ct_command(0,0,.016f).pressed&NET_INPUT_SWITCH));
    fixture();equipment[0].gun[0]=EQ_RAIL;equipment[0].charge=.4f;equipment[0].held=true;players[0].matter=40;equipment[0].buildHeld=true;ct_open();assert(!equipment[0].charge&&!equipment[0].held&&!equipment[0].buildHeld&&players[0].matter>40);ct_close(false);
    puts("Hold/toggle, tap/hold, revive priority and safe charge cancellation PASS");
}
static void tuning(void){
    float results[3];int rates[]={30,60,144};for(int r=0;r<3;r++){fixture();axes[3][GAMEPAD_AXIS_RIGHT_X]=.7f;ct_profile(2)->tuning[2].assist=0;float dt=1.f/rates[r];for(int n=0;n<rates[r];n++){ct_frame(dt);ct_move(2,dt);}results[r]=players[2].yaw;}
    assert(fabsf(results[0]-results[2])<.15f&&fabsf(results[1]-results[2])<.15f);
    fixture();ct_profile(0)->tuning[1].x=2;motion=(Vector2){10,0};ct_frame(.016f);NetInputCommand c=ct_command(0,0,.016f);assert(fabsf(c.aim_yaw+3)<.001f);
    uint8_t data[128];NetWriter w;NetReader reader;net_writer_init(&w,data,sizeof(data));net_write_input(&w,&c);net_reader_init(&reader,data,w.length);NetInputCommand decoded;assert(net_read_input(&reader,&decoded));Player host=players[0];net_apply_movement_command(&host,&decoded,.016f);assert(host.yaw==c.aim_yaw);
    puts("30/60/144 Hz turning and serialized LAN aim equivalence PASS");
}
static void lan(void){
    fixture();ct_profile(0)->bindings[0][CT_SWITCH][0]=(CtBinding){1,KEY_N,0};
    netTransport.role=NET_ROLE_CLIENT;netLocalPlayerCount=1;netLocalPlayerSlots[0]=1;
    tk[KEY_N]=true;motion=(Vector2){14,2};ct_frame(.016f);NetInputCommand command=ct_command(0,1,.016f);command.sequence=17;
    assert(command.pressed&NET_INPUT_SWITCH);
    netTransport=(NetTransport){0};assert(net_transport_global_init());assert(net_transport_host(&netTransport,27947));NetTransport client={0};assert(net_transport_connect(&client,"127.0.0.1",27947));
    double start=GetTime();while(GetTime()-start<3&&(!client.connected||!net_transport_connected_clients(&netTransport))){net_transport_pump(&netTransport,NULL,NULL,NULL);net_transport_pump(&client,NULL,NULL,NULL);WaitTime(.002);}
    assert(client.connected&&net_transport_connected_clients(&netTransport)==1);netPeerPlayerMask[1]=2;netLastProcessedInput[1]=0;netHasPendingInput[1]=false;
    uint8_t packet[128];NetWriter writer;net_writer_init(&writer,packet,sizeof(packet));net_write_header(&writer,NET_MSG_INPUT,0,netTransport.session_id,0);net_write_u8(&writer,1);net_write_input(&writer,&command);
    assert(net_transport_send(&client,0,FPS_NET_CHANNEL_INPUT,packet,writer.length,true));
    start=GetTime();while(GetTime()-start<3&&!netHasPendingInput[1]){net_transport_pump(&client,NULL,NULL,NULL);net_transport_pump(&netTransport,net_on_receive,NULL,NULL);WaitTime(.002);}
    assert(netHasPendingInput[1]);assert(netPendingInput[1].pressed&NET_INPUT_SWITCH);assert(netPendingInput[1].aim_yaw==command.aim_yaw);
    net_transport_shutdown(&client);net_transport_shutdown(&netTransport);net_transport_global_shutdown();netTransport.role=NET_ROLE_OFFLINE;
    puts("Loopback LAN: remapped client input reaches authoritative host PASS");
}
static void camera_fov(void) {
    fixture();
    ct_profile(1)->verticalFov = 90;
    assert(aim_fov(0) == 60 && aim_fov(1) == 90);
    for (int seat = 0; seat < 2; ++seat) {
        float baseFov = aim_fov(seat);
        aimLensGun[seat] = EQ_RAIL;
        aimStates[seat].blend = 1;
        float magnification = tanf(baseFov * 0.5f * DEG2RAD) / tanf(aim_fov(seat) * 0.5f * DEG2RAD);
        assert(fabsf(magnification - eqGuns[EQ_RAIL].zoom) < .001f);
        aimStates[seat].blend = 0;
        assert(aim_fov(seat) == baseFov);
    }
    netTransport.role = NET_ROLE_CLIENT;
    netLocalPlayerCount = 2;
    netLocalPlayerSlots[0] = 3;
    netLocalPlayerSlots[1] = 2;
    assert(aim_fov(3) == 60 && aim_fov(2) == 90);
    assert(aim_fov(1) == 60);
    netTransport.role = NET_ROLE_OFFLINE;
    puts("FOV: split-screen profiles, LAN seat mapping and zoom magnification PASS");
}

static void fov_menu(void) {
    fixture();
    ctTab = ctEditSeat = 0;
    ct_open();
    ctRow = CT_PROFILE_FOV;
    tp[KEY_LEFT] = true;
    for (int press = 0; press < 5; ++press) ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 45 && aim_fov(0) == 60);

    neutral();
    int width = GetScreenWidth();
    int listWidth = (int)(width * .62f);
    position = (Vector2){30 + listWidth - 24, 170 + CT_PROFILE_FOV * 32 + 8};
    tmp[MOUSE_BUTTON_LEFT] = true;
    ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 50);
    position.x -= 34;
    ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 45);

    neutral();
    couchMenuDevice = 3;
    pp[3][GAMEPAD_BUTTON_LEFT_FACE_RIGHT] = true;
    ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 50);
    neutral();
    tp[KEY_RIGHT] = true;
    for (int press = 0; press < 10; ++press) ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 90);
    neutral();
    tp[KEY_F5] = true;
    ct_menu_frame();
    assert(!ctMenuOpen && aim_fov(0) == 90 && aim_fov(1) == 60);
    CtLibrary saved;
    assert(ct_read_file("controls.cfg", &saved) && saved.profiles[0].verticalFov == 90);

    neutral();
    ct_open();
    ctRow = CT_PROFILE_RESET_FIRST;
    tp[KEY_ENTER] = true;
    ct_menu_frame();
    assert(ct_edit_profile()->verticalFov == 60);
    neutral();
    tp[KEY_ESCAPE] = true;
    ct_menu_frame();
    assert(!ctMenuOpen && aim_fov(0) == 90);
    assert(ct_read_file("controls.cfg", &saved) && saved.profiles[0].verticalFov == 90);
    neutral();
    puts("FOV menu: keyboard, mouse, controller, limits, reset, apply and cancel PASS");
}

static void menus(void){
    fixture();gameState=GAME_STATE_PAUSED;ctPauseOpened=true;couchSelection=0;tp[KEY_PAUSE]=true;couch_frame();assert(gameState==GAME_STATE_PAUSED);tp[KEY_PAUSE]=true;couch_frame();assert(gameState==GAME_STATE_PLAYING);neutral();
    ct_open();CtLibrary saved=ctLibrary;ctDraft.profiles[0].tuning[1].x=3;ct_close(false);assert(!memcmp(&saved,&ctLibrary,sizeof(saved)));
    ct_open();ct_edit_profile()->verticalFov=75;ctNameAction=1;strcpy(ctNameBuffer,"My layout");ct_name_finish();assert(ctDraft.count==5&&!strcmp(ct_edit_profile()->name,"My layout")&&ct_edit_profile()->verticalFov==75);
    ctTab=1;ctDevice=0;ctCapture=CT_FIRE;ctBindingSlot=0;ctCandidate=ct_edit_profile()->bindings[0][CT_JUMP][0];ctConflict=ct_conflict(ct_edit_profile(),CT_FIRE,ctCandidate);assert(ctConflict==CT_JUMP);CtBinding old=ct_edit_profile()->bindings[0][CT_FIRE][0];ct_assign_binding(0);assert(ct_equal_binding(ct_edit_profile()->bindings[0][CT_JUMP][0],old));
    ct_close(true);CtLibrary loaded;assert(ct_read_file("controls.cfg",&loaded));assert(!memcmp(&loaded,&ctLibrary,sizeof(loaded)));
    ct_open();couchMenuDevice=3;ctDevice=2;ctCapture=CT_JUMP;ctCaptureArmed=true;pp[3][GAMEPAD_BUTTON_RIGHT_FACE_RIGHT]=pb[3][GAMEPAD_BUTTON_RIGHT_FACE_RIGHT]=true;ct_frame(.016f);ct_menu_frame();assert(ctCapture<0);neutral();ct_close(false);
    fixture();ct_open();for(int size=0;size<2;size++){SetWindowSize(size?640:1280,size?360:720);for(int tab=0;tab<5;tab++){ctTab=tab;ctRow=ctScroll=0;ct_menu_frame();ct_menu_frame();Image image=LoadImageFromScreen();assert(ExportImage(image,TextFormat("controls-%d-%d.png",size,tab)));UnloadImage(image);}}
    ct_close(false);puts("Menu draft/apply/cancel, profile duplicate, conflict swap, B capture and screenshots PASS");
}
int main(void){SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Controls tests");MakeDirectory("artifacts");MakeDirectory("artifacts/controls-test");assert(ChangeDirectory("artifacts/controls-test"));persistence();actions();modes();tuning();lan();camera_fov();fov_menu();menus();CloseWindow();puts("Controls regression suite PASS");return 0;}
