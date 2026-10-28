static void combat_ui_reset(void) {
    memset(deathNotices,0,sizeof(deathNotices));deathContext=DEATH_UNKNOWN;
    memset(damageCues,0,sizeof(damageCues));memset(reviveComplete,0,sizeof(reviveComplete));
    combatUiReceiveTick=0;memset(chargeReceiveTick,0,sizeof(chargeReceiveTick));damageOriginKnown=false;
}
static const char *combat_ui_death_text(int victim) {
    DeathNotice *n=&deathNotices[victim];
    const char *causes[]={"Unknown cause","Fell out of the arena","Voxel shot","Crushed by debris","Explosion","Melee","Shotgun","Cutter beam","Railgun","Rocket","Sword","Glob gun","Glock","DMR"};
    if(n->cause==DEATH_FALL || n->attacker<0)return causes[n->cause];
    if(n->attacker==victim)return TextFormat("%s - self-inflicted",causes[n->cause]);
    const char *names[]={"Rifleman","Swarmer","Goliath","Scrapper","Shieldbearer","Sapper","Lobber","Builder","Anchor","Foundry","Golem"};
    return n->enemy>=0?TextFormat("%s - %s",causes[n->cause],names[n->enemy]):TextFormat("%s - Player %d",causes[n->cause],n->attacker+1);
}
static void combat_ui_death(int victim,int attacker,bool debris,bool downed) {
    if(victim<0 || victim>=MAX_PLAYERS || is_player_bot(victim) || netTransport.role==NET_ROLE_CLIENT)return;
    DeathCause cause=deathContext;
    if(cause==DEATH_UNKNOWN)cause=players[victim].pos.y< -10?DEATH_FALL:debris?DEATH_DEBRIS:explosionDamageContext?DEATH_EXPLOSION:attacker>=0?DEATH_MELEE:DEATH_UNKNOWN;
    int enemy=attacker>=0 && attacker<activePlayers && is_player_bot(attacker)?players[attacker].enemyType:-1;
    if(go_combat_actor(attacker))enemy=ENEMY_TYPE_COUNT;
    deathNotices[victim]=(DeathNotice){cause,attacker,enemy,downed,(float)GetTime()+5};
    if(netTransport.role==NET_ROLE_HOST){
        uint8_t packet[48];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_COMBAT_UI,0,netTransport.session_id,netServerTick);
        net_write_u8(&w,2);net_write_u8(&w,victim);net_write_u8(&w,cause);net_write_u8(&w,attacker+1);net_write_u8(&w,enemy+1);net_write_u8(&w,downed);
        net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_CONTROL,packet,w.length,true);
    }
}
static void combat_ui_results(void) {
    int w=GetScreenWidth();bool victory=firefightWaveState==FIREFIGHT_STATE_VICTORY;
    ClearBackground((Color){15,18,26,255});
    Color accent=victory?GOLD:(Color){235,75,65,255};
    DrawRectangle(0,0,w,8,accent);
    const char *title=victory?"FIREFIGHT VICTORY":"FIREFIGHT DEFEAT";
    int font=w<800?32:52;DrawText(title,(w-MeasureText(title,font))/2,30,font,accent);
    const char *subtitle=victory?"The arena is clear.":"TEAM WIPED - No teammates or reserves remaining.";
    int small=w<800?16:20;DrawText(subtitle,(w-MeasureText(subtitle,small))/2,95,small,RAYWHITE);
    const char *stats=TextFormat("Wave %d  |  Score %d  |  Enemies defeated %d",firefightWave,firefightScore,firefightKillsTotal);
    DrawText(stats,(w-MeasureText(stats,small))/2,130,small,GOLD);
    for(int i=0;i<firefightHumans;i++){
        int y=180+i*60;DrawRectangle(30,y-6,w-60,54,(Color){27,32,44,255});
        DrawText(TextFormat("P%d   %d kills   %d physics kills   %d revives",i+1,players[i].kills,players[i].debrisKills,reviveCount[i]),45,y,small,player_palette_color(i));
        if(deathNotices[i].until>0)DrawText(TextFormat("Last hit: %s",combat_ui_death_text(i)),45,y+25,small-2,LIGHTGRAY);
    }
}
static void combat_ui_add_damage(int victim,Vector3 origin,bool known) {
    if(victim<0 || victim>=MAX_PLAYERS)return;
    int slot=0;float oldest=1e30f;
    for(int n=0;n<8;n++){
        DamageCue *c=&damageCues[victim][n];
        if(c->known==known && (!known || v_length(world_delta(c->origin,origin))<1)){slot=n;break;}
        if(c->until<oldest){oldest=c->until;slot=n;}
    }
    damageCues[victim][slot]=(DamageCue){origin,(float)GetTime()+1,known};
}
static void combat_ui_damage(int victim,int attacker) {
    if(victim<0 || victim>=MAX_PLAYERS || netTransport.role==NET_ROLE_CLIENT)return;
    bool known=damageOriginKnown || aiProjectileDamage || (attacker>=0 && attacker<activePlayers && attacker!=victim);
    Vector3 origin=damageOriginKnown?damageOrigin:aiProjectileDamage?aiProjectileOrigin:known?players[attacker].pos:(Vector3){0};
    combat_ui_add_damage(victim,origin,known);
    if(netTransport.role==NET_ROLE_HOST){
        uint8_t packet[48];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_COMBAT_UI,0,netTransport.session_id,netServerTick);
        net_write_u8(&w,1);net_write_u8(&w,victim);net_write_u8(&w,known);
        net_write_f32(&w,origin.x);net_write_f32(&w,origin.y);net_write_f32(&w,origin.z);
        net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_CONTROL,packet,w.length,true);
    }
}
static void combat_ui_shutdown(void) {
    if(reviveShader.id && reviveShader.id!=rlGetShaderIdDefault())UnloadShader(reviveShader);
    reviveShader=(Shader){0};reviveShaderAttempted=false;
}
static void combat_ui_world(void) {
    if(gameMode!=GAME_MODE_FIREFIGHT)return;
    for(int i=0;i<firefightHumans;i++){
        bool active=combatantLife[i]==LIFE_DOWNED && reviveOwner[i]>=0 && reviveProgress[i]>0;
        if(!active && reviveComplete[i]<=0)continue;
        Vector3 base=active?players[i].death_pos:players[i].pos;base.y-=BASE_EYE_HEIGHT-.03f;
        Color color=player_palette_color(active?reviveOwner[i]:i);
        float progress=active?clampf(reviveProgress[i]/FF_REVIVE_SECONDS,0,1):1;
        if(!reviveShaderAttempted){
            reviveShaderAttempted=true;
            reviveShader=LoadShaderFromMemory(
                "#version 330\nin vec3 vertexPosition;in vec4 vertexColor;uniform mat4 mvp;out vec3 p;out vec4 c;void main(){p=vertexPosition;c=vertexColor;gl_Position=mvp*vec4(vertexPosition,1);}",
                "#version 330\nin vec3 p;in vec4 c;uniform float clock;out vec4 finalColor;void main(){float band=pow(.5+.5*sin(p.y*25.-clock*12.),8.);finalColor=vec4(c.rgb*(.7+band*.8),c.a*(.15+.85*band));}");
        }
        bool shader=reviveShader.id && reviveShader.id!=rlGetShaderIdDefault();
        float clock=(float)GetTime();
        rlDrawRenderBatchActive();rlDisableDepthMask();
        if(shader){SetShaderValue(reviveShader,GetShaderLocation(reviveShader,"clock"),&clock,SHADER_UNIFORM_FLOAT);BeginShaderMode(reviveShader);}
        DrawCylinder(base,.65f,.65f,1.2f,24,Fade(color,active?.35f:reviveComplete[i]));
        if(shader)EndShaderMode();
        for(int n=0;n<48;n++){
            float a=n*2*PI/48,b=(n+1)*2*PI/48;
            DrawLine3D(v_add(base,(Vector3){cosf(a)*.8f,0,sinf(a)*.8f}),v_add(base,(Vector3){cosf(b)*.8f,0,sinf(b)*.8f}),n<progress*48?color:Fade(color,.15f));
        }
        rlDrawRenderBatchActive();rlEnableDepthMask();
    }
}
static Vector2 combat_ui_direction(Vector3 delta,Camera3D camera) {
    Vector3 f=v_sub(camera.target,camera.position);f.y=0;f=v_norm(f);
    Vector3 right=v_cross(f,(Vector3){0,1,0});
    Vector2 d={v_dot(delta,right),-v_dot(delta,f)};
    return Vector2Length(d)>.001f?Vector2Normalize(d):(Vector2){0,-1};
}
static void combat_ui_hud(int viewer,Camera3D camera,int w,int h) {
    if(viewer<0 || viewer>=MAX_PLAYERS)return;
    DeathNotice *notice=&deathNotices[viewer];
    if(notice->until>0 && players[viewer].respawn_timer>0){
        int font=w<700?16:24,y=(int)(h*.32f);
        DrawRectangle(12,y-12,w-24,font*2+32,Fade(BLACK,.8f));
        const char *title=notice->downed && combatantLife[viewer]==LIFE_DOWNED?"DOWNED":"YOU DIED";
        DrawText(title,(w-MeasureText(title,font))/2,y,font,RED);
        const char *cause=combat_ui_death_text(viewer);int detail=font-2;
        DrawText(cause,(w-MeasureText(cause,detail))/2,y+font+8,detail,RAYWHITE);
    }
    for(int n=0;n<8;n++){
        DamageCue *cue=&damageCues[viewer][n];float alpha=clampf(cue->until-(float)GetTime(),0,1);if(alpha<=0)continue;
        Color color=Fade(RED,alpha*.85f);
        if(!cue->known){DrawRectangleLinesEx((Rectangle){8,8,w-16,h-16},5,color);continue;}
        Vector2 d=combat_ui_direction(world_delta(cue->origin,camera.position),camera);
        float sx=(w*.5f-22)/fmaxf(.001f,fabsf(d.x)),sy=(h*.5f-22)/fmaxf(.001f,fabsf(d.y));
        Vector2 at={w*.5f+d.x*fminf(sx,sy),h*.5f+d.y*fminf(sx,sy)},side={-d.y*20,d.x*20};
        Vector2 tip=Vector2Subtract(at,Vector2Scale(d,12));
        DrawLineEx(Vector2Add(at,side),tip,7,color);DrawLineEx(tip,Vector2Subtract(at,side),7,color);
    }
    for(int i=0;i<firefightHumans && gameMode==GAME_MODE_FIREFIGHT;i++)if(reviveProgress[i]>0 && (i==viewer || reviveOwner[i]==viewer)){
        Vector2 at={w*.5f,h*.72f};float fraction=clampf(reviveProgress[i]/FF_REVIVE_SECONDS,0,1);
        DrawCircleV(at,26,Fade(BLACK,.75f));DrawRing(at,21,25,-90,-90+fraction*360,48,player_palette_color(reviveOwner[i]));
        const char *label=i==viewer?"BEING REVIVED":TextFormat("REVIVING P%d",i+1);int font=w<700?14:18;
        DrawText(label,(int)at.x-MeasureText(label,font)/2,(int)at.y+31,font,WHITE);
    }
}
static void combat_ui_abilities(int viewer,int w,int h) {
    if(viewer<0 || viewer>=MAX_PLAYERS || creativeModeActive)return;
    int local=netTransport.role==NET_ROLE_OFFLINE?viewer:net_local_index_for_slot(viewer);if(local<0)return;
    int size=w<700?40:56;float y=h-size-80;
    bool usable=players[viewer].respawn_timer<=0 && combatantLife[viewer]==LIFE_ALIVE && !ff_input_blocked(viewer);
    for(int n=0;n<2;n++){
        float x=w-18-(2-n)*(size+12);Vector2 c={x+size*.5f,y+size*.5f};
        bool ready=usable && (n?ex_grenade_count(viewer,exSelectedGrenade[viewer])>0 && grenadeCooldown[viewer]<=0:holograms[viewer].cooldown<=0);
        Color ink=ready?(n?(exSelectedGrenade[viewer]==GRENADE_GRAVITY?PURPLE:exSelectedGrenade[viewer]==GRENADE_IGLOO?SKYBLUE:ORANGE):SKYBLUE):(Color){100,106,114,255};
        DrawRectangleRounded((Rectangle){x,y,size,size},.2f,6,Fade(BLACK,.8f));
        if(n){
            if(exSelectedGrenade[viewer]==GRENADE_GRAVITY){DrawRing(c,size*.16f,size*.24f,0,360,24,ink);DrawEllipseLines((int)c.x,(int)c.y,size*.32f,size*.12f,ink);}
            else if(exSelectedGrenade[viewer]==GRENADE_IGLOO){DrawCircleSector((Vector2){c.x-3,c.y+8},size*.28f,180,360,8,ink);DrawRectangle(c.x-6,c.y,8,12,BLACK);DrawRectangle(c.x+4,c.y-size*.28f,6,12,ink);}
            else {DrawEllipse((int)c.x-5,(int)c.y+2,size*.19f,size*.25f,ink);DrawRectangle((int)c.x-10,(int)y+6,12,6,ink);DrawLineEx((Vector2){c.x,y+8},(Vector2){c.x+7,y+18},3,ink);}
            DrawText(TextFormat("%d",ex_grenade_count(viewer,exSelectedGrenade[viewer])),(int)x+size-17,(int)y+size-25,size<50?22:28,ink);
        }else{
            DrawCircleV((Vector2){c.x-9,c.y-8},size*.13f,ink);DrawCircleV((Vector2){c.x+8,c.y-8},size*.13f,ink);
            DrawEllipse((int)c.x,(int)c.y+3,size*.18f,size*.23f,ink);DrawCircleV((Vector2){c.x,c.y+10},3,BLACK);
            if(holograms[viewer].cooldown>0)DrawRing(c,size*.43f,size*.47f,-90,-90+360*(1-clampf(holograms[viewer].cooldown/15,0,1)),32,SKYBLUE);
            if(holograms[viewer].life>0)DrawCircle((int)x+size-5,(int)y+5,5,SKYBLUE);
        }
        bool pad=playerInput[local]==INPUT_TYPE_GAMEPAD;
        const char *key=n?(pad?"D-UP":local==0?"C":"."):(pad?"LB":local==0?"V":"/");
        DrawText(key,(int)c.x-MeasureText(key,12)/2,(int)y+size+3,12,LIGHTGRAY);
    }
}
static void combat_ui_snapshot(int peer) {
    if(netTransport.role!=NET_ROLE_HOST)return;
    uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
    net_write_header(&w,NET_MSG_COMBAT_UI,0,netTransport.session_id,netServerTick);net_write_u8(&w,0);
    for(int i=0;i<MAX_PLAYERS;i++){
        net_write_u8(&w,combatantLife[i]);net_write_u8(&w,reviveOwner[i]+1);net_write_f32(&w,reviveProgress[i]);net_write_f32(&w,reviveComplete[i]);
        net_write_f32(&w,players[i].death_pos.x);net_write_f32(&w,players[i].death_pos.y);net_write_f32(&w,players[i].death_pos.z);
    }
    if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    for(int part=0;part<4;part++){
        net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_CHARGES,0,netTransport.session_id,netServerTick);net_write_u8(&w,part);
        for(int n=part*16;n<(part+1)*16;n++){
            EnemyCharge *c=&enemyCharges[n];net_write_u8(&w,c->active);net_write_u8(&w,c->owner);net_write_u64(&w,c->id);
            net_write_f32(&w,c->pos.x);net_write_f32(&w,c->pos.y);net_write_f32(&w,c->pos.z);net_write_f32(&w,c->fuse);
        }
        if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    }
}
static bool combat_ui_receive(NetMessageType type,NetReader *r,uint32_t tick) {
    if(type!=NET_MSG_COMBAT_UI && type!=NET_MSG_CHARGES)return false;
    if(netTransport.role!=NET_ROLE_CLIENT)return true;
    if(type==NET_MSG_CHARGES){
        int part=net_read_u8(r);if(part>=4 || tick<chargeReceiveTick[part])return true;
        EnemyCharge incoming[16]={0};
        for(int n=0;n<16;n++){
            EnemyCharge *c=&incoming[n];c->active=net_read_u8(r)!=0;c->owner=net_read_u8(r);c->id=net_read_u64(r);
            c->pos=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};c->fuse=net_read_f32(r);
            if(r->failed || c->owner>=MAX_COMBATANTS || !v_isfinite(c->pos) || !isfinite(c->fuse) || c->fuse>4 || c->fuse< -1)return true;
        }
        for(int n=0;n<16;n++){EnemyCharge old=enemyCharges[part*16+n];if(old.id==incoming[n].id)incoming[n].beep=old.beep;}
        memcpy(enemyCharges+part*16,incoming,sizeof(incoming));chargeReceiveTick[part]=tick;return true;
    }
    int kind=net_read_u8(r);
    if(kind==2){
        int victim=net_read_u8(r),cause=net_read_u8(r),attacker=(int)net_read_u8(r)-1,enemy=(int)net_read_u8(r)-1,downed=net_read_u8(r);
        if(!r->failed && victim<MAX_PLAYERS && cause<DEATH_CAUSE_COUNT && attacker<MAX_COMBATANTS && enemy<=ENEMY_TYPE_COUNT && downed<=1)
            deathNotices[victim]=(DeathNotice){cause,attacker,enemy,downed,(float)GetTime()+5};
        return true;
    }
    if(kind==1){int victim=net_read_u8(r),known=net_read_u8(r);Vector3 origin={net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        if(!r->failed && victim<MAX_PLAYERS && known<=1 && v_isfinite(origin))combat_ui_add_damage(victim,origin,known);return true;}
    if(kind!=0 || tick<combatUiReceiveTick)return true;
    int life[4],owner[4];float progress[4],complete[4];Vector3 pos[4];
    for(int i=0;i<MAX_PLAYERS;i++){
        life[i]=net_read_u8(r);owner[i]=(int)net_read_u8(r)-1;progress[i]=net_read_f32(r);complete[i]=net_read_f32(r);
        pos[i]=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        if(r->failed || life[i]>LIFE_UNUSED || owner[i]>=MAX_PLAYERS || !isfinite(progress[i]) || progress[i]<0 || progress[i]>FF_REVIVE_SECONDS || (progress[i]>0 && owner[i]<0) || !isfinite(complete[i]) || complete[i]<0 || complete[i]>.6f || !v_isfinite(pos[i]))return true;
    }
    for(int i=0;i<MAX_PLAYERS;i++){combatantLife[i]=life[i];reviveOwner[i]=owner[i];reviveProgress[i]=progress[i];reviveComplete[i]=complete[i];players[i].death_pos=pos[i];}
    combatUiReceiveTick=tick;return true;
}
