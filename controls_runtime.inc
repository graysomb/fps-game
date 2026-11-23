/* Physical inputs are sampled once per application frame; consumers read actions. */
static struct {bool keys[512],mouse[7],pad[4][32],available[4];float axes[4][8],wheel;Vector2 delta;} ctRaw;
static bool ctArmed[4][CT_COUNT],ctToggle[4],ctTetherWas[4],ctCommandValid[4];
static int ctLastDevice[4]={-99,-99,-99,-99},ctLastContext[4];
static Vector2 ctVelocity[4],ctSmooth[4];
static NetInputCommand ctCommands[4];
static int ctActionDevice[4][CT_COUNT];
static bool ctWasBlocked;
static void ct_mouse_reset(void){ctRaw.delta=(Vector2){0};memset(ctSmooth,0,sizeof(ctSmooth));if(mouseSeat>=0&&mouseSeat<4){CtProfile *p=ct_profile(mouseSeat);for(int a=0;a<CT_COUNT;a++)for(int b=0;b<2;b++){CtBinding binding=p->bindings[1][a][b];if(binding.kind==2&&ctRaw.mouse[binding.code]){ctArmed[mouseSeat][a]=false;ctState[mouseSeat][a]=(CtActionState){0};}}}}
static bool ct_pad(int seat){int slot=netTransport.role==NET_ROLE_OFFLINE?seat:netLocalPlayerSlots[seat];int input=netTransport.role==NET_ROLE_CLIENT?seat:slot;return input>=0&&input<MAX_PLAYERS&&playerInput[input]==INPUT_TYPE_GAMEPAD;}
static int ct_seat(int actor){return netTransport.role==NET_ROLE_OFFLINE?actor:net_local_index_for_slot(actor);}
static int ct_context(void){if(creativeModeActive&&gameState==GAME_STATE_CREATIVE)return CT_COMMON|CT_EDIT;int context=CT_COMMON|CT_GAME;if(enemyLab&&!cmTesting)context|=CT_LAB;if(cmTesting)context|=CT_TEST;if(gameMode==GAME_MODE_FIREFIGHT&&(firefightWaveState==FIREFIGHT_STATE_LOOT||firefightWaveState==FIREFIGHT_STATE_INTERMISSION))context|=CT_PREP;return context;}
static bool ct_binding_down(CtBinding b,int seat,int device){
    int mods=(ctRaw.keys[KEY_LEFT_CONTROL]||ctRaw.keys[KEY_RIGHT_CONTROL]?1:0)|(ctRaw.keys[KEY_LEFT_SHIFT]||ctRaw.keys[KEY_RIGHT_SHIFT]?2:0)|(ctRaw.keys[KEY_LEFT_ALT]||ctRaw.keys[KEY_RIGHT_ALT]?4:0);
    // A modifier used as the key itself does not require itself in the chord mask.
    if(b.kind==1){if(b.code==KEY_LEFT_CONTROL||b.code==KEY_RIGHT_CONTROL)mods&=~1;if(b.code==KEY_LEFT_SHIFT||b.code==KEY_RIGHT_SHIFT)mods&=~2;if(b.code==KEY_LEFT_ALT||b.code==KEY_RIGHT_ALT)mods&=~4;}
    if(b.kind==1)return b.code>0&&b.code<512&&ctRaw.keys[b.code]&&(mods&b.mods)==b.mods;
    if(b.kind==2)return seat==mouseSeat&&mouseCaptured&&b.code>=0&&b.code<7&&ctRaw.mouse[b.code];
    if(b.kind==3)return seat==mouseSeat&&mouseCaptured&&ctRaw.wheel*b.code>0;
    if(b.kind==4 && device>=0&&device<4&&ctRaw.available[device]){
        if(b.code==GAMEPAD_BUTTON_LEFT_TRIGGER_2||b.code==GAMEPAD_BUTTON_RIGHT_TRIGGER_2){int axis=b.code==GAMEPAD_BUTTON_LEFT_TRIGGER_2?GAMEPAD_AXIS_LEFT_TRIGGER:GAMEPAD_AXIS_RIGHT_TRIGGER;return (ctRaw.axes[device][axis]+1)*.5f>=ct_profile(seat)->tuning[2].trigger;}
        return b.code>=0&&b.code<32&&ctRaw.pad[device][b.code];
    }return false;
}
static bool ct_down(int seat,int action){return seat>=0&&seat<4&&action>=0&&action<CT_COUNT&&ctState[seat][action].down;}
static bool ct_pressed(int seat,int action){return seat>=0&&seat<4&&action>=0&&action<CT_COUNT&&ctState[seat][action].pressed;}
static bool ct_actor_pressed(int actor,int action){return ct_pressed(ct_seat(actor),action);}
static void ct_cancel_all(void){
    memset(ctToggle,0,sizeof(ctToggle));memset(ctTetherWas,0,sizeof(ctTetherWas));memset(ctVelocity,0,sizeof(ctVelocity));memset(ctSmooth,0,sizeof(ctSmooth));memset(ctArmed,0,sizeof(ctArmed));memset(ctState,0,sizeof(ctState));memset(ctCommandValid,0,sizeof(ctCommandValid));
    for(int seat=0;seat<4;seat++){int i=netTransport.role==NET_ROLE_OFFLINE?seat:netLocalPlayerSlots[seat];if(i<0||i>=MAX_PLAYERS)continue;
        if(netTransport.role!=NET_ROLE_CLIENT){Vector3 zero={0};release_tether_velocity(i,&zero);eq_cancel_charge(i);}if(netTransport.role==NET_ROLE_CLIENT&&netTransport.connected){NetInputCommand c={.sequence=++netInputSequence[i],.client_tick=netServerTick,.pressed=NET_INPUT_AIM_CANCEL|NET_INPUT_BUILD_CANCEL|NET_INPUT_TETHER_CANCEL,.mouse_look=1,.aim_yaw=players[i].yaw,.aim_pitch=players[i].pitch};uint8_t packet[128];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_INPUT,0,netTransport.session_id,netServerTick);net_write_u8(&w,i);net_write_input(&w,&c);net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,w.length,true);}
        equipment[i].buildHeld=false;aim_cancel(i);aimStates[i].tapDown=false;aimStates[i].tapTime=0;aimStates[i].tapArmed=false;aimStates[i].holdUsed=false;}
}
static void ct_frame(float dt){
    (void)dt;if(!ctReady){ct_defaults(&ctLibrary);ctReady=true;}memset(ctCommandValid,0,sizeof(ctCommandValid));
    for(int k=0;k<512;k++)ctRaw.keys[k]=IsKeyDown(k);
    for(int k=0;k<7;k++)ctRaw.mouse[k]=IsMouseButtonDown(k);
    ctRaw.delta=GetMouseDelta();ctRaw.wheel=GetMouseWheelMove();
    for(int d=0;d<4;d++){ctRaw.available[d]=IsGamepadAvailable(d);for(int b=0;b<32;b++)ctRaw.pad[d][b]=ctRaw.available[d]&&IsGamepadButtonDown(d,b);for(int a=0;a<8;a++)ctRaw.axes[d][a]=ctRaw.available[d]?GetGamepadAxisMovement(d,a):0;}
    int context=ct_context();
    bool inputBlocked=ctMenuOpen||!IsWindowFocused()||(gameState!=GAME_STATE_PLAYING&&gameState!=GAME_STATE_CREATIVE);
    bool blocked=inputBlocked||debugConsoleOpen;
    if(blocked&&!ctWasBlocked)ct_cancel_all();ctWasBlocked=blocked;
    for(int s=0;s<4;s++){
        bool pad=ct_pad(s);int device=pad?seat_device(s):-1;CtProfile *p=ct_profile(s);
        bool changed=ctLastDevice[s]!=device||ctLastContext[s]!=context;ctLastDevice[s]=device;ctLastContext[s]=context;
        int actor=netTransport.role==NET_ROLE_OFFLINE?s:netLocalPlayerSlots[s];
        bool dead=gameState==GAME_STATE_PLAYING&&actor>=0&&actor<MAX_PLAYERS&&players[actor].respawn_timer>0;
        bool disabled=blocked||(pad&&(device<0||device>=4||!ctRaw.available[device]));
        if(disabled||changed){ctToggle[s]=false;ctVelocity[s]=(Vector2){0};ctSmooth[s]=(Vector2){0};}
        for(int a=0;a<CT_COUNT;a++){
            // Player 1 owns the console shortcut even when playing with a controller.
            bool console=s==0&&a==CT_CONSOLE;
            bool down=false,pulse=false;for(int d=pad&&!console?2:0;d<=(pad?2:1);d++)for(int b=0;b<2;b++){
                CtBinding binding=p->bindings[d][a][b];bool active=ct_binding_down(binding,s,device);
                // Modified editor shortcuts take precedence over their plain-key actions.
                if(active&&binding.kind==1)for(int other=0;other<CT_COUNT;other++)if(ctDefinitions[other].context&context)for(int n=0;n<2;n++){CtBinding candidate=p->bindings[0][other][n];int modifier=(binding.code==KEY_LEFT_CONTROL||binding.code==KEY_RIGHT_CONTROL)?1:(binding.code==KEY_LEFT_SHIFT||binding.code==KEY_RIGHT_SHIFT)?2:(binding.code==KEY_LEFT_ALT||binding.code==KEY_RIGHT_ALT)?4:0;
                    if(((candidate.code==binding.code&&candidate.mods!=binding.mods&&(candidate.mods&binding.mods)==binding.mods)||(modifier&&(candidate.mods&modifier)))&&ct_binding_down(candidate,s,device))active=false;}
                if(active&&(context&CT_LAB)&&ctDefinitions[a].context==CT_GAME)for(int other=0;other<CT_COUNT;other++)if(ctDefinitions[other].context==CT_LAB)for(int n=0;n<2;n++){CtBinding candidate=p->bindings[d][other][n];if(candidate.kind==binding.kind&&candidate.code==binding.code&&candidate.mods==binding.mods&&ct_binding_down(candidate,s,device))active=false;}
                if(active){ctActionDevice[s][a]=d;pulse|=binding.kind==3;}down|=active;
            }
            bool capture=(console?inputBlocked:disabled)||(dead&&ctDefinitions[a].context==CT_GAME)||(ce.active&&!ce.testing&&(ce.menu||ceSeats[s].panel)&&(a==CT_PLACE||a==CT_ERASE||a==CT_SAMPLE));
            if(capture||changed)ctArmed[s][a]=false;if(!down||pulse)ctArmed[s][a]=true;
            if(!(ctDefinitions[a].context&context)||capture||!ctArmed[s][a])down=false;
            CtActionState old=ctState[s][a];ctState[s][a]=(CtActionState){down,down&&(!old.down||pulse),!down&&old.down};
        }
    }
}
static float ct_axis(float v,CtTuning *t){float a=fabsf(v);if(a<=t->inner)return 0;if(t->curve<.5f)return Clamp(v/(1-t->outer),-1,1);float x=Clamp((a-t->inner)/(1-t->inner-t->outer),0,1);if(t->curve>1.5f)x=x*x;return copysignf(x,v);}
static Vector2 ct_move_axes(int s){if(s<0||s>=4)return (Vector2){0};if(ctMenuOpen||debugConsoleOpen||!IsWindowFocused())return (Vector2){0};if(ct_pad(s)){CtTuning movement=ct_profile(s)->tuning[2];movement.inner=movement.moveInner;movement.outer=movement.moveOuter;CtTuning *t=&movement;int d=seat_device(s);if(d<0||d>=4||!ctRaw.available[d])return (Vector2){0};int x=t->swap?GAMEPAD_AXIS_RIGHT_X:GAMEPAD_AXIS_LEFT_X,y=t->swap?GAMEPAD_AXIS_RIGHT_Y:GAMEPAD_AXIS_LEFT_Y;Vector2 move={ct_axis(ctRaw.axes[d][x],t),-ct_axis(ctRaw.axes[d][y],t)};if(t->curve<.5f && Vector2Length(move)>0)move=Vector2Normalize(move);return move;}return (Vector2){ct_down(s,CT_RIGHT)-ct_down(s,CT_LEFT),ct_down(s,CT_FORWARD)-ct_down(s,CT_BACK)};}
static Vector2 ct_look_axes(int s){if(ct_pad(s)){CtTuning *t=&ct_profile(s)->tuning[2];int d=seat_device(s);if(d<0||d>=4||!ctRaw.available[d])return (Vector2){0};CtTuning look=*t;return (Vector2){ct_axis(ctRaw.axes[d][t->swap?GAMEPAD_AXIS_LEFT_X:GAMEPAD_AXIS_RIGHT_X],&look),ct_axis(ctRaw.axes[d][t->swap?GAMEPAD_AXIS_LEFT_Y:GAMEPAD_AXIS_RIGHT_Y],&look)};}return (Vector2){ct_down(s,CT_LOOK_RIGHT)-ct_down(s,CT_LOOK_LEFT),ct_down(s,CT_LOOK_DOWN)-ct_down(s,CT_LOOK_UP)};}
static float ct_turn(float input,float *velocity,float speed,float acceleration,float brake,float dt,bool direct){
    if(direct){*velocity=input*speed;return *velocity*dt;}
    // Integrate in bounded substeps for consistent ramp and stopping distances.
    float result=0;int steps=(int)ceilf(dt*240);if(steps<1)steps=1;float h=dt/steps;
    for(int n=0;n<steps;n++){float before=*velocity;if(input!=0)*velocity=Clamp(*velocity+input*acceleration*h,-speed,speed);else{float next=fmaxf(0,fabsf(*velocity)-brake*h);*velocity=copysignf(next,*velocity);}result+=(before+*velocity)*.5f*h;}return result;
}
static void ct_look(int s,int slot,Player *p,float dt,bool editor){
    CtProfile *profile=ct_profile(s);bool pad=ct_pad(s);CtTuning *t=&profile->tuning[pad?2:0];Vector2 look=ct_look_axes(s);float zoom=editor?1:aim_scale(slot);zoom*=1+(t->zoom-1)*(1-zoom);
    if(t->invertX)look.x=-look.x;if(t->invertY)look.y=-look.y;
    bool mouse=!pad&&s==mouseSeat&&mouseCaptured;float speed=editor?120:TURN_SPEED;
    p->yaw-=ct_turn(look.x,&ctVelocity[s].x,speed*t->x*zoom,TURN_ACCELERATION*t->x*t->accel*zoom,TURN_FRICTION*t->brake,dt,t->direct||t->accel<=0||mouse||editor);
    p->pitch-=ct_turn(look.y,&ctVelocity[s].y,speed*t->y*zoom,TURN_ACCELERATION*t->y*t->accel*zoom,TURN_FRICTION*t->brake,dt,t->direct||t->accel<=0||mouse||editor);
    if(mouse){t=&profile->tuning[1];Vector2 delta=ctRaw.delta;float smooth=t->smooth;if(smooth>0){float blend=1-expf(-dt/smooth);Vector2 velocity=dt>0?Vector2Scale(delta,1/dt):(Vector2){0};ctSmooth[s]=Vector2Lerp(ctSmooth[s],velocity,blend);delta=Vector2Scale(ctSmooth[s],dt);}float rate=dt>0?Vector2Length(ctRaw.delta)/dt:0;float gain=1+t->accel*fminf(rate/1000,2);float scale=.15f*gain*(editor?1:aim_scale(slot))*(editor?1:1+(t->zoom-1)*(1-aim_scale(slot)));p->yaw-=delta.x*scale*t->x*(t->invertX?-1:1);p->pitch-=delta.y*scale*t->y*(t->invertY?-1:1);}
    p->yaw=remainderf(p->yaw,360);p->pitch=Clamp(p->pitch,-89,89);p->yaw_vel=p->pitch_vel=0;
}
static float ct_assist(int slot){int s=ct_seat(slot);if(!ctReady||s<0||s>=4)return 1;return ct_pad(s)?ct_profile(s)->tuning[2].assist:1;}
static bool ct_safe(int s,int slot){return !ctMenuOpen&&!debugConsoleOpen&&IsWindowFocused()&&gameState==GAME_STATE_PLAYING&&players[slot].respawn_timer<=0&&!ff_input_blocked(slot)&&!(enemyLab&&enemyLabPaused)&&(!ct_pad(s)||IsGamepadAvailable(seat_device(s)));}
static NetInputCommand ct_command(int s,int slot,float dt){
    NetInputCommand c={0};if(s<0||s>=4||slot<0||slot>=MAX_PLAYERS)return c;if(ctCommandValid[s])return ctCommands[s];ctCommandValid[s]=true;
    c.mouse_look=1;c.aim_yaw=players[slot].yaw;c.aim_pitch=players[slot].pitch;
    if(!ct_safe(s,slot)){ctToggle[s]=false;c.pressed=NET_INPUT_AIM_CANCEL|NET_INPUT_BUILD_CANCEL|NET_INPUT_TETHER_CANCEL;}
    else {
        Vector2 move=ct_move_axes(s);c.move_x=net_quantize_axis(move.x);c.move_y=net_quantize_axis(move.y);
        Player aim=players[slot];ct_look(s,slot,&aim,dt,false);
        if(ct_pad(s)||!mouse_owns(slot)){InputType savedInput=playerInput[slot];playerInput[slot]=ct_pad(s)?INPUT_TYPE_GAMEPAD:INPUT_TYPE_KEYBOARD;float yaw=players[slot].yaw,pitch=players[slot].pitch;players[slot].yaw=aim.yaw;players[slot].pitch=aim.pitch;Vector2 look=ct_look_axes(s);apply_aim_assist(slot,dt,fmaxf(fabsf(look.x),fabsf(look.y)));aim.yaw=players[slot].yaw;aim.pitch=players[slot].pitch;players[slot].yaw=yaw;players[slot].pitch=pitch;playerInput[slot]=savedInput;}
        c.aim_yaw=aim.yaw;c.aim_pitch=aim.pitch;
        static const struct {int action;uint32_t flag;} map[]={{CT_FIRE,NET_INPUT_FIRE},{CT_JUMP,NET_INPUT_JUMP},{CT_MELEE,NET_INPUT_MELEE},{CT_BUILD,NET_INPUT_BUILD},{CT_SWITCH,NET_INPUT_SWITCH},{CT_EQUIP,NET_INPUT_EQUIP},{CT_REVIVE,NET_INPUT_REVIVE},{CT_GRENADE,NET_INPUT_GRENADE},{CT_CYCLE,NET_INPUT_GRENADE_CYCLE},{CT_HOLOGRAM,NET_INPUT_HOLOGRAM}};
        for(unsigned n=0;n<sizeof(map)/sizeof(map[0]);n++){if(ct_pressed(s,map[n].action))c.pressed|=map[n].flag;if(ct_down(s,map[n].action))c.held|=map[n].flag;}
        CtProfile *p=ct_profile(s);int d=ct_pad(s)?2:ctActionDevice[s][CT_AIM];
        if(p->aimToggle[d]){if(ct_pressed(s,CT_AIM))c.pressed|=NET_INPUT_AIM;}else{if(ct_down(s,CT_AIM)&&!aimStates[slot].active)c.pressed|=NET_INPUT_AIM;if(!ct_down(s,CT_AIM))c.pressed|=NET_INPUT_AIM_CANCEL;}
        int td=ct_pad(s)?2:ctActionDevice[s][CT_TETHER];if(p->tetherToggle[td]){if(ct_pressed(s,CT_TETHER))ctToggle[s]=!ctToggle[s];}else ctToggle[s]=ct_down(s,CT_TETHER);if(ctToggle[s])c.held|=NET_INPUT_TETHER;
        aimStates[slot].revive=false;
        if(p->combined[ct_pad(s)?2:ctActionDevice[s][CT_INTERACT]])aim_interact(slot,ct_down(s,CT_INTERACT),dt,&c);
        if(ct_down(s,CT_REVIVE))aimStates[slot].revive=true;
    }
    if(ct_context()&CT_PREP){if(ct_pressed(s,CT_READY))c.pressed|=NET_INPUT_PREP_READY;if(ct_pressed(s,CT_PREP_LEFT))c.pressed|=NET_INPUT_PREP_LEFT;if(ct_pressed(s,CT_PREP_RIGHT))c.pressed|=NET_INPUT_PREP_RIGHT;if(ct_pressed(s,CT_PREP_MENU))c.pressed|=NET_INPUT_PREP_MENU;}
    ctCommands[s]=c;return c;
}
static void ct_move(int slot,float dt){int s=ct_seat(slot);if(s<0)return;NetInputCommand c=ct_command(s,slot,dt);c.pressed&=~NET_INPUT_JUMP;net_apply_movement_command(&players[slot],&c,dt);}
static void ct_local_actions(float dt){
    if(netTransport.role==NET_ROLE_CLIENT)return;
    for(int i=0;i<MAX_PLAYERS&&i<activePlayers;i++)if(!is_player_bot(i)&&(netTransport.role==NET_ROLE_OFFLINE||net_player_is_local(i))){int s=ct_seat(i);if(s<0)continue;NetInputCommand c=ct_command(s,i,dt);bool tether=(c.held&NET_INPUT_TETHER)!=0;
        if(c.pressed&NET_INPUT_TETHER_CANCEL){Vector3 zero={0};release_tether_velocity(i,&zero);eq_cancel_charge(i);}else if(tether&&!ctTetherWas[s])start_tether(i);else if(!tether&&ctTetherWas[s])release_tether(i);ctTetherWas[s]=tether;
        if(c.pressed&NET_INPUT_MELEE)perform_melee(i);if(c.pressed&NET_INPUT_GRENADE)clarity_throw(i);if(c.pressed&NET_INPUT_HOLOGRAM)hologram_deploy(i);if((c.pressed&NET_INPUT_JUMP)&&players[i].onGround){players[i].vel.y=JUMP_SPEED;players[i].onGround=false;}}
}
static void ct_edit_look(int seat,float dt){ct_look(seat,seat,&players[seat],dt,true);}
static const char *ct_binding_name(CtBinding b){
    static char names[16][80];static int next;char *out=names[next++&15];const char *key="Unbound";
    if(b.kind==1){switch(b.code){case KEY_SPACE:key="Space";break;case KEY_ENTER:key="Enter";break;case KEY_TAB:key="Tab";break;case KEY_ESCAPE:key="Esc";break;case KEY_LEFT_SHIFT:key="LShift";break;case KEY_RIGHT_SHIFT:key="RShift";break;case KEY_LEFT_CONTROL:key="LCtrl";break;case KEY_RIGHT_CONTROL:key="RCtrl";break;case KEY_UP:key="Up";break;case KEY_DOWN:key="Down";break;case KEY_LEFT:key="Left";break;case KEY_RIGHT:key="Right";break;default:key=b.code>=32&&b.code<=96?TextFormat("%c",b.code):TextFormat("Key %d",b.code);}}
    if(b.kind==2)key=b.code==0?"Mouse Left":b.code==1?"Mouse Right":b.code==2?"Mouse Middle":TextFormat("Mouse %d",b.code+1);
    if(b.kind==3)key=b.code>0?"Wheel Up":"Wheel Down";
    if(b.kind==4){static const char *buttons[]={"?","DUp","DRight","DDown","DLeft","Y","B","A","X","LB","LT","RB","RT","Back","Guide","Start","L3","R3"};key=b.code<18?buttons[b.code]:"Pad";}
    snprintf(out,80,"%s%s%s%s",b.mods&1?"Ctrl+":"",b.mods&2?"Shift+":"",b.mods&4?"Alt+":"",key);return out;
}
static const char *ct_hint(int actor,int action){int s=ct_seat(actor);if(s<0||s>=4)s=0;CtProfile *p=ct_profile(s);int d=ct_pad(s)?2:0;CtBinding b=p->bindings[d][action][0];if(!ct_pad(s)&&mouse_owns(actor)&&p->bindings[1][action][0].kind)b=p->bindings[1][action][0];if(!b.kind&&(action==CT_EQUIP||action==CT_REVIVE)&&p->combined[d])return TextFormat("Hold %s",ct_binding_name(p->bindings[d][CT_INTERACT][0]));return ct_binding_name(b);}
