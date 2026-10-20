/* One physical mouse belongs to the first local keyboard seat. */
static void mouse_angles(Player *p,Vector2 delta,float keysX,float keysY,float dt,int seat){
    float scale=.15f*localSeats[seat].sensitivity;
    p->yaw=remainderf(p->yaw-delta.x*scale-keysX*TURN_SPEED*dt,360);
    p->pitch=clampf(p->pitch-delta.y*scale*(localSeats[seat].invert?-1:1)-keysY*TURN_SPEED*dt,-89,89);
    p->yaw_vel=p->pitch_vel=0;
}
static void mouse_sync(void){
    int seat=-1,slot=-1;
    int count=netTransport.role==NET_ROLE_OFFLINE?multiplayerActivePlayers:netLocalPlayerCount;
    for(int s=0;s<count && s<MAX_PLAYERS;s++){
        int actor=netTransport.role==NET_ROLE_OFFLINE?s:netLocalPlayerSlots[s];
        int input=netTransport.role==NET_ROLE_CLIENT?s:actor;
        if(actor>=0 && input>=0 && playerInput[input]==INPUT_TYPE_KEYBOARD){seat=s;slot=actor;break;}
    }
    bool capture=slot>=0 && IsWindowFocused() && !debugConsoleOpen &&
        (gameState==GAME_STATE_PLAYING || gameState==GAME_STATE_CREATIVE) && !(enemyLab && enemyLabPaused);
    if(mouseCaptured && (!capture || slot!=mouseSlot)){
        if(netTransport.role==NET_ROLE_CLIENT && mouseSlot>=0){
            NetInputCommand cancel={.sequence=++netInputSequence[mouseSlot],.client_tick=netServerTick,
                .pressed=NET_INPUT_TETHER_CANCEL,.mouse_look=1,.aim_yaw=remainderf(players[mouseSlot].yaw,360),.aim_pitch=clampf(players[mouseSlot].pitch,-89,89)};
            uint8_t packet[128];NetWriter writer;net_writer_init(&writer,packet,sizeof(packet));
            net_write_header(&writer,NET_MSG_INPUT,0,netTransport.session_id,netServerTick);
            net_write_u8(&writer,(uint8_t)mouseSlot);net_write_input(&writer,&cancel);
            if(netTransport.connected)net_transport_send(&netTransport,0,FPS_NET_CHANNEL_INPUT,packet,writer.length,true);
        }
        else if(mouseTetherHeld && mouseSlot>=0){Vector3 zero={0};release_tether_velocity(mouseSlot,&zero);}
        EnableCursor();mouseCaptured=false;mouseTetherHeld=false;mouseButtonsArmed=false;
    }
    mouseSeat=seat;mouseSlot=slot;
    if(capture && !mouseCaptured){DisableCursor();mouseCaptured=true;mouseButtonsArmed=false;mouseFrameDelta=(Vector2){0};}
}
static void mouse_command(NetInputCommand *c,int seat,int slot,float dt){
    if(slot!=mouseSlot)return;
    Player aim=players[slot];
    bool playable=mouseCaptured && players[slot].respawn_timer<=0 && !ff_input_blocked(slot);
    if(playable)mouse_angles(&aim,mouseFrameDelta,net_axis_value(c->look_x),net_axis_value(c->look_y),dt,seat);
    c->mouse_look=1;c->aim_yaw=remainderf(aim.yaw,360);c->aim_pitch=clampf(aim.pitch,-89,89);c->look_x=c->look_y=0;
    if(!playable){c->pressed=c->held=0;c->move_x=c->move_y=0;}
    if(playable && mouseButtonsArmed && !netRequestedCreative && gameState!=GAME_STATE_CREATIVE){
        if(IsMouseButtonPressed(MOUSE_BUTTON_LEFT))c->pressed|=NET_INPUT_FIRE;
        if(IsMouseButtonDown(MOUSE_BUTTON_LEFT))c->held|=NET_INPUT_FIRE;
        if(GetMouseWheelMove()!=0)c->pressed|=NET_INPUT_SWITCH;
        if(IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE))c->pressed|=NET_INPUT_MELEE;
        if(IsMouseButtonDown(MOUSE_BUTTON_RIGHT))c->held|=NET_INPUT_TETHER;
    }
}
static void mouse_frame(float dt){
    bool wasCaptured=mouseCaptured;int previous=mouseSlot;mouse_sync();
    mouseFrameDelta=mouseCaptured && wasCaptured && previous==mouseSlot?GetMouseDelta():(Vector2){0};
    if(!mouseButtonsArmed && wasCaptured && previous==mouseSlot && mouseCaptured &&
       !IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !IsMouseButtonDown(MOUSE_BUTTON_RIGHT) && !IsMouseButtonDown(MOUSE_BUTTON_MIDDLE) &&
       !IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && !IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE))mouseButtonsArmed=true;
    if(netTransport.role==NET_ROLE_CLIENT || mouseSlot<0)return;
    NetInputCommand c={0};int i=mouseSlot;
    bool first=mouseSeat==0;
    c.look_x=net_quantize_axis((IsKeyDown(first?KEY_H:KEY_RIGHT)?1.f:0)-(IsKeyDown(first?KEY_F:KEY_LEFT)?1.f:0));
    c.look_y=net_quantize_axis((IsKeyDown(first?KEY_G:KEY_DOWN)?1.f:0)-(IsKeyDown(first?KEY_T:KEY_UP)?1.f:0));
    if(IsKeyDown(first?KEY_R:KEY_P))c.held|=NET_INPUT_TETHER;
    if(IsKeyPressed(first?KEY_LEFT_CONTROL:KEY_RIGHT_CONTROL))c.pressed|=NET_INPUT_FIRE;
    if(IsKeyPressed(first?KEY_Z:KEY_M))c.pressed|=NET_INPUT_MELEE;
    mouse_command(&c,mouseSeat,i,dt);
    players[i].yaw=c.aim_yaw;players[i].pitch=c.aim_pitch;players[i].yaw_vel=players[i].pitch_vel=0;
    if(gameState!=GAME_STATE_PLAYING || !mouseCaptured)return;

    if(c.pressed&NET_INPUT_MELEE)perform_melee(i);
    bool held=(c.held&NET_INPUT_TETHER)!=0;
    if(held && !mouseTetherHeld)start_tether(i);
    if(!held && mouseTetherHeld)release_tether(i);
    mouseTetherHeld=held;
}
