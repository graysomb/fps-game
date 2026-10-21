/* Aiming is cosmetic: simulation rays and weapon damage remain unchanged. */
static Shader aimShader;
static bool aimShaderTried;
static int aimResolutionLoc,aimBlendLoc,aimStyleLoc,aimChargeLoc;
static int aimLensGun[MAX_PLAYERS];
static int aimDevices[MAX_PLAYERS]={-2,-2,-2,-2},aimInputs[MAX_PLAYERS]={-2,-2,-2,-2};
static bool aim_local(int i){return netTransport.role==NET_ROLE_OFFLINE || net_player_is_local(i);}
static void aim_cancel(int i){if(i>=0 && i<MAX_PLAYERS)aimStates[i].active=false;}
static float aim_zoom(int i){if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return 1;float a=aimStates[i].blend;a=a*a*(3-2*a);int gun=aimLensGun[i];return 1+(eqGuns[gun].zoom-1)*a;}
static float aim_fov(int i){return 2*atanf(tanf(30*DEG2RAD)/aim_zoom(i))*RAD2DEG;}
static float aim_scale(int i){return 1/aim_zoom(i);}
static float aim_player_scale(const Player *p){for(int i=0;i<MAX_PLAYERS;i++)if(p==&players[i])return aim_scale(i);return 1;}
static float aim_seat_scale(int seat){int i=netTransport.role==NET_ROLE_OFFLINE?seat:netLocalPlayerSlots[seat];return aim_scale(i);}
static void aim_command(int i,const NetInputCommand *c){
    if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return;
    AimState *a=&aimStates[i];
    uint32_t cancel=NET_INPUT_AIM_CANCEL|NET_INPUT_SWITCH|NET_INPUT_EQUIP|NET_INPUT_MELEE|NET_INPUT_BUILD|NET_INPUT_GRENADE|NET_INPUT_HOLOGRAM;
    bool blocked=gameState!=GAME_STATE_PLAYING || players[i].respawn_timer>0 || ff_input_blocked(i) || creativeModeActive || netRequestedCreative || (enemyLab&&enemyLabPaused);
    bool interrupted=(c->pressed&cancel)!=0 || (c->held&(NET_INPUT_TETHER|NET_INPUT_REVIVE))!=0;
    if(blocked || interrupted || players[i].meleeSwingActive || equipment[i].switchTime>0){a->active=false;}
    else if(c->pressed&NET_INPUT_AIM){a->active=!a->active;if(a->active)aimLensGun[i]=eq_gun(i);}
    if(netTransport.role==NET_ROLE_CLIENT && (interrupted || (c->pressed&NET_INPUT_AIM)))a->pending=c->sequence;
}
/* A hold never becomes a pickup tap on release, even if revival was unavailable. */
static void aim_interact(int i,bool down,float dt,NetInputCommand *c){
    AimState *a=&aimStates[i];a->revive=false;
    if(!a->tapArmed){if(!down)a->tapArmed=true;return;}
    if(down){a->tapTime+=fmaxf(0,dt);a->tapDown=true;if(a->tapTime>=.25f){a->revive=true;c->held|=NET_INPUT_REVIVE;}}
    else {if(a->tapDown && a->tapTime<.25f)c->pressed|=NET_INPUT_EQUIP;a->tapDown=false;a->tapTime=0;}
}
static void aim_sample(int seat,int i,float dt,NetInputCommand *c){
    if(seat<0 || i<0 || i>=MAX_PLAYERS || gameState!=GAME_STATE_PLAYING || creativeModeActive || (enemyLab&&!equipmentLab))return;
    bool pad=playerInput[netTransport.role==NET_ROLE_CLIENT?seat:i]==INPUT_TYPE_GAMEPAD;int d=seat_device(seat);
    if(!IsWindowFocused() || debugConsoleOpen || players[i].respawn_timer>0 || (pad&&!IsGamepadAvailable(d))){c->pressed|=NET_INPUT_AIM_CANCEL;aimStates[i].toggleArmed=false;aimStates[i].tapArmed=false;aimStates[i].tapDown=false;aimStates[i].tapTime=0;aimStates[i].revive=false;return;}
    if(pad){
        aim_interact(i,IsGamepadButtonDown(d,GAMEPAD_BUTTON_RIGHT_FACE_LEFT),dt,c);
        if(!IsGamepadButtonDown(d,GAMEPAD_BUTTON_RIGHT_THUMB))aimStates[i].toggleArmed=true;
        if(aimStates[i].toggleArmed && IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_THUMB))c->pressed|=NET_INPUT_AIM;
        if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_FACE_DOWN)||IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_TRIGGER_1)||IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_TRIGGER_1)||IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_UP)||IsGamepadButtonDown(d,GAMEPAD_BUTTON_LEFT_TRIGGER_2))c->pressed|=NET_INPUT_AIM_CANCEL;
    }else{
        bool first=seat==0;
        if(!IsKeyDown(first?KEY_LEFT_SHIFT:KEY_ENTER))aimStates[i].toggleArmed=true;
        if(aimStates[i].toggleArmed && IsKeyPressed(first?KEY_LEFT_SHIFT:KEY_ENTER))c->pressed|=NET_INPUT_AIM;
        if(IsKeyPressed(first?KEY_Z:KEY_M)||IsKeyPressed(first?KEY_E:KEY_O)||IsKeyPressed(first?KEY_C:KEY_PERIOD)||IsKeyPressed(first?KEY_V:KEY_SLASH)||IsKeyDown(first?KEY_R:KEY_P)||IsKeyDown(first?KEY_X:KEY_COMMA))c->pressed|=NET_INPUT_AIM_CANCEL;
        if(mouse_owns(i)&&mouseCaptured&&mouseButtonsArmed&&(IsMouseButtonDown(MOUSE_BUTTON_RIGHT)||IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)))c->pressed|=NET_INPUT_AIM_CANCEL;
    }
}
static void aim_tick(float dt){
    for(int i=0;i<MAX_PLAYERS;i++){
        AimState *a=&aimStates[i];if(players[i].respawn_timer>0 || is_player_bot(i)){a->active=false;a->blend=0;a->tapArmed=false;a->tapDown=false;a->revive=false;a->tapTime=0;}
        if(a->active)aimLensGun[i]=eq_gun(i);
        a->blend=clampf(a->blend+(a->active?1:-1)*dt/.15f,0,1);
    }
}
static void aim_focus_tick(void){
    for(int i=0;i<MAX_PLAYERS;i++)if(aim_local(i)){
        int seat=netTransport.role==NET_ROLE_OFFLINE?i:net_local_index_for_slot(i);if(seat<0)continue;
        bool pad=playerInput[netTransport.role==NET_ROLE_CLIENT?seat:i]==INPUT_TYPE_GAMEPAD;
        int device=pad?seat_device(seat):-1,input=pad?1:0;bool changed=aimDevices[i]!=device || aimInputs[i]!=input;aimDevices[i]=device;aimInputs[i]=input;
        if(!changed && gameState==GAME_STATE_PLAYING && IsWindowFocused() && !debugConsoleOpen && !(enemyLab&&enemyLabPaused) && (!pad||IsGamepadAvailable(seat_device(seat))))continue;
        AimState *a=&aimStates[i];bool notify=a->active||a->blend>0;uint32_t pending=a->pending;*a=(AimState){.pending=pending};
        if(notify && netTransport.role==NET_ROLE_CLIENT && netTransport.connected){
            NetInputCommand c={.sequence=++netInputSequence[i],.pressed=NET_INPUT_AIM_CANCEL};a->pending=c.sequence;
            uint8_t data[128];NetWriter w;net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,netServerTick);net_write_u8(&w,i);net_write_input(&w,&c);net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,data,w.length,true);
        }
    }
}
static void aim_shutdown(void){if(aimShader.id && aimShader.id!=rlGetShaderIdDefault())UnloadShader(aimShader);aimShader=(Shader){0};aimShaderTried=false;}
static void aim_overlay(int i,int w,int h){
    if(i<0 || i>=MAX_PLAYERS || aimStates[i].blend<=.001f || players[i].respawn_timer>0)return;
    float blend=aimStates[i].blend,charge=equipment[i].charge/eqRailCharge;int style=aimLensGun[i];
    if(!aimShaderTried){aimShaderTried=true;aimShader=LoadShader(NULL,"shaders/aim_sight.frag");aimResolutionLoc=GetShaderLocation(aimShader,"resolution");aimBlendLoc=GetShaderLocation(aimShader,"blend");aimStyleLoc=GetShaderLocation(aimShader,"style");aimChargeLoc=GetShaderLocation(aimShader,"railCharge");}
    if(aimShader.id && aimShader.id!=rlGetShaderIdDefault()){
        Vector2 size={(float)w,(float)h};SetShaderValue(aimShader,aimResolutionLoc,&size,SHADER_UNIFORM_VEC2);SetShaderValue(aimShader,aimBlendLoc,&blend,SHADER_UNIFORM_FLOAT);SetShaderValue(aimShader,aimStyleLoc,&style,SHADER_UNIFORM_INT);SetShaderValue(aimShader,aimChargeLoc,&charge,SHADER_UNIFORM_FLOAT);
        BeginShaderMode(aimShader);DrawRectangle(0,0,w,h,WHITE);EndShaderMode();
    }else{Color c=Fade(style==EQ_RAIL?VIOLET:SKYBLUE,blend);int x=w/2,y=h/2;DrawLine(x-12,y,x-3,y,c);DrawLine(x+3,y,x+12,y,c);DrawLine(x,y-12,x,y-3,c);DrawLine(x,y+3,x,y+12,c);}
    if(style==EQ_RAIL)DrawText(TextFormat("4x  CHARGE %d%%",(int)(charge*100)),w/2-55,h/2+(int)(h*.11f),h<400?12:16,Fade(VIOLET,blend));
}
