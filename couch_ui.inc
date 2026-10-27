/* Offline couch shell. LAN keeps its existing UI and four-seat protocol. */
static bool couchLobby=false, couchAdvanced=false, couchSettings=false;
static bool couchReady[MAX_PLAYERS];
static int couchSelection=0, couchSettingSeat=0, couchSettingRow=0;
static float couchNavDelay=0;

static void couch_load_settings(void) {
    FILE *f=fopen("couch_settings.cfg","r");
    if(!f) return;
    for(int i=0;i<MAX_PLAYERS;i++) {
        float sensitivity,deadzone; int invert;
        if(fscanf(f,"%f %f %d",&sensitivity,&deadzone,&invert)!=3) break;
        if(isfinite(sensitivity)&&isfinite(deadzone)&&sensitivity>=.5f&&sensitivity<=2&&deadzone>=.05f&&deadzone<=.3f&&(invert==0||invert==1)) {
            localSeats[i].sensitivity=sensitivity; localSeats[i].deadzone=deadzone; localSeats[i].invert=invert!=0;
        }
    }
    fclose(f);
}
static void couch_save_settings(void) {
    FILE *f=fopen("couch_settings.cfg","w"); if(!f) return;
    for(int i=0;i<MAX_PLAYERS;i++) fprintf(f,"%.3f %.3f %d\n",localSeats[i].sensitivity,localSeats[i].deadzone,localSeats[i].invert);
    fclose(f);
}
static bool couch_button(int seat, int button) {
    return seat>=0 && seat<MAX_PLAYERS && multiplayerPlayerInput[seat]==INPUT_TYPE_GAMEPAD &&
        IsGamepadAvailable(seat_device(seat)) && IsGamepadButtonPressed(seat_device(seat),button);
}
/* Menu ownership is independent of gameplay seats and survives hot-plugging. */
static int couchMenuDevice=-1;
static int couch_menu_device(void) {
    if(couchMenuDevice>=0 && IsGamepadAvailable(couchMenuDevice)) return couchMenuDevice;
    couchMenuDevice=-1;couchNavDelay=0;
    for(int d=0;d<4;d++) if(IsGamepadAvailable(d)){couchMenuDevice=d;break;}
    return couchMenuDevice;
}
static bool couch_menu_button(int button) {int d=couch_menu_device();return d>=0 && IsGamepadButtonPressed(d,button);}
static bool couch_confirm(void) { return couch_menu_device()<0?IsKeyPressed(KEY_ENTER):couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_DOWN); }
static bool couch_back(void) { return couch_menu_device()<0?IsKeyPressed(KEY_ESCAPE):couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_RIGHT); }
static int couch_axis(bool horizontal) {
    int d=couch_menu_device();
    if(d<0) {
        if(IsKeyPressed(horizontal?KEY_LEFT:KEY_UP)) return -1;
        if(IsKeyPressed(horizontal?KEY_RIGHT:KEY_DOWN)) return 1;
        return 0;
    }
    if(IsGamepadButtonPressed(d,horizontal?GAMEPAD_BUTTON_LEFT_FACE_LEFT:GAMEPAD_BUTTON_LEFT_FACE_UP)) return -1;
    if(IsGamepadButtonPressed(d,horizontal?GAMEPAD_BUTTON_LEFT_FACE_RIGHT:GAMEPAD_BUTTON_LEFT_FACE_DOWN)) return 1;
    float a=GetGamepadAxisMovement(d,horizontal?GAMEPAD_AXIS_LEFT_X:GAMEPAD_AXIS_LEFT_Y);
    if(fabsf(a)>.6f && couchNavDelay<=0) {couchNavDelay=.2f;return a<0?-1:1;}
    return 0;
}
static void couch_title(const char *title) {
    ClearBackground((Color){16,21,30,255});
    DrawText(title,48,40,36,RAYWHITE);
}
static void couch_row(int row,const char *label,bool selected) {
    int y=116+row*48;
    if(selected) DrawRectangle(40,y-6,(gameState==GAME_STATE_MENU && !couchLobby && !couchSettings && !couchAdvanced?GetScreenWidth()/2-64:GetScreenWidth()-80),42,(Color){40,77,94,255});
    DrawText(label,56,y,24,selected?GOLD:LIGHTGRAY);
}
static void couch_begin_lobby(bool practice) {
    enemyLab=false;
    if(gameMode==GAME_MODE_DEATHMATCH) {currentWorldType=WORLD_TYPE_CROSSCUT;useCustomMap=false;}
    firefightPractice=practice;
    couchLobby=true; couchSelection=0;
    multiplayerActivePlayers=1;
    int device=couch_menu_device();
    multiplayerPlayerInput[0]=device>=0?INPUT_TYPE_GAMEPAD:INPUT_TYPE_KEYBOARD;
    if(device>=0)localSeats[0].device=device;
    memset(couchReady,0,sizeof(couchReady));
}
static void couch_launch(void) {
    if(firefightPractice) multiplayerActivePlayers=1;
    if(gameMode==GAME_MODE_FIREFIGHT && (firefightFoundry||firefightPractice)) useCustomMap=false;
    ResetGame(); gameState=GAME_STATE_PLAYING; couchLobby=false;
}
static void ff_draw_hud(int i,int w,int h) {
    if(enemyLab) {enemy_lab_hud(w,h);return;}
    if(gameMode!=GAME_MODE_FIREFIGHT || i>=firefightHumans) return;
    if(firefightWaveState==FIREFIGHT_STATE_LOOT && combatantLife[i]==LIFE_ALIVE){
        int ready=0,total=0,fs=w<700?12:16;
        for(int n=0;n<firefightHumans;n++)if(combatantLife[n]!=LIFE_UNUSED){total++;ready+=prepReady[n];}
        int seat=netTransport.role==NET_ROLE_CLIENT?net_local_index_for_slot(i):i;
        if(seat<0)seat=i;
        const char *key=playerInput[seat]==INPUT_TYPE_GAMEPAD?"A":seat==0?"Space":"Right Shift";
        const char *line=TextFormat("%s: %s | Skip loot (%d/%d ready)",key,prepReady[i]?"Cancel":"Ready",ready,total);
        int width=MeasureText(line,fs),x=(w-width)/2;
        DrawRectangle(x-5,62,width+10,fs+4,Fade(BLACK,.65f));DrawText(line,x,64,fs,prepReady[i]?GOLD:RAYWHITE);
    }
    if(firefightWaveState==FIREFIGHT_STATE_INTERMISSION && !firefightPractice && combatantLife[i]==LIFE_ALIVE){
        int fs=w<700?12:16,step=fs+4,y=64;
        bool pad=playerInput[i]==INPUT_TYPE_GAMEPAD;
        const WaveDefinition *next=&ff_waves[firefightWave%FF_WAVE_COUNT];int total=0;
        for(int t=0;t<FF_TYPE_COUNT;t++)total+=ff_count(firefightWave+1,fp_team_size(),t);
        const char *lines[3];
        lines[0]=TextFormat("%s %s",ff_upgrades_maxed()?"All upgrades maxed":ff_upgrade_label(prepVotes[i]<0?0:prepVotes[i]),prepReady[i]?"[READY]":"");
        lines[1]=TextFormat("Next: %s (%d enemies)",next->name,total);
        lines[2]=prepMenuOpen[i]?(pad?"D-pad vote | A ready | Y build":"Look vote | Jump ready | V/N build"):(pad?"Y: vote menu":"V/N: vote menu");
        int rows=prepMenuOpen[i]?3:1;
        for(int n=0;n<rows;n++){const char *line=rows==1?lines[2]:lines[n];int width=MeasureText(line,fs),x=(w-width)/2;
            DrawRectangle(x-5,y+n*step-2,width+10,step,Fade(BLACK,.65f));DrawText(line,x,y+n*step,fs,n==0?GOLD:RAYWHITE);}
        return;
    }
    int font=w<700?16:20;
    for(int boss=MAX_PLAYERS;boss<activePlayers;boss++)if(players[boss].respawn_timer<=0 && players[boss].enemyType==ENEMY_TYPE_FOUNDRY) {
        EnemySpecialState *s=&enemySpecial[boss];int lost=0;for(int bit=0;bit<4;bit++)lost+=(s->lostMask>>bit)&1;
        int barWidth=w<700?230:340,x=(w-barWidth)/2,y=70;
        DrawRectangle(x-6,y-5,barWidth+12,44,Fade(BLACK,.85f));
        DrawText(enemy_core_open(boss)?TextFormat("CORE EXPOSED %.0f / %.0f",s->coreHealth,s->coreMax):TextFormat("FOUNDRY: Break 2 supports (%d/2)",lost<2?lost:2),x,y,w<700?12:16,enemy_core_open(boss)?GOLD:ORANGE);
        if(s->blockedCue>0)DrawText("CORE PROTECTED - SHOOT ORANGE SUPPORTS",x,y+38,w<700?10:14,ORANGE);
        if(s->coreCue>0)DrawText("SUPPORTS BROKEN - ATTACK THE GOLD CORE",x,y+38,w<700?10:14,GOLD);
        DrawRectangle(x,y+24,barWidth,7,DARKGRAY);
        DrawRectangle(x,y+24,(int)(barWidth*clampf(s->coreHealth/fmaxf(1,s->coreMax),0,1)),7,GOLD);break;
    }
    const char *text=NULL;
    if(combatantLife[i]==LIFE_DOWNED) text=TextFormat("DOWNED | Awaiting revive %.1f / 3s",reviveProgress[i]);
    else if(combatantLife[i]==LIFE_SPECTATING) text="Spectating until wave clear";
    else if(combatantLife[i]==LIFE_WAITING) text=TextFormat("Respawn in %.0fs",players[i].respawn_timer);
    else if(firefightPractice) {
        const char *steps[]={"Melee a block to harvest Matter", "Shoot the enemy until its shield breaks", "Finish the exposed enemy with melee", "Expose the next enemy; throw debris to finish", "Training complete!"};
        text=steps[practiceStep<5?practiceStep:4];
    } else {
        for(int j=0;j<firefightHumans;j++) if(combatantLife[j]==LIFE_DOWNED) {
            text=TextFormat("P%d: Hold %s to revive",j+1,playerInput[i]==INPUT_TYPE_GAMEPAD?"X":i==0?"X":","); break;
        }
        if(!text && players[i].matter<1) text="No Matter: melee blocks to harvest";
        if(!text) {
            Vector3 dir; float angle;
            if(find_aim_assist_target(i,&dir,&angle) && angle<12) {
                for(int j=MAX_PLAYERS;j<MAX_COMBATANTS;j++) if(players[j].respawn_timer<=0 && players[j].isExposed) {
                    Vector3 delta=v_norm(v_sub(players[j].pos,players[i].pos));
                    if(v_dot(delta,player_forward(&players[i]))>.97f) {
                        text=playerInput[i]==INPUT_TYPE_GAMEPAD?"EXPOSED: melee (A/bottom face) or throw (LT)":i==0?"EXPOSED: melee Z or throw R":"EXPOSED: melee M or throw P"; break;
                    }
                }
            }
        }
        if(!text && firefightEnemiesSpawned<firefightWaveEnemiesTotal) text=TextFormat("Incoming: %s entrance",ffEntrance?"EAST":"WEST");
    }
    if(text) {
        /* Keep the hint between the radar and ability icons, above the grenade row. */
        int panelX=HUD_PADDING_X+2*hud_radar_radius(h)+16;
        int panelW=w-18-2*((w<700?40:56)+12)-8-panelX-6;
        if(panelW<1)return;
        int fs=font;
        while(fs>12 && MeasureText(text,fs)>panelW) fs--;
        char lines[4][128]={{0}};int count=0,width=0;
        const char *next=text;
        while(*next && count<4){
            int length=0,space=0;
            while(next[length] && length<(int)sizeof(lines[0])-1){
                lines[count][length]=next[length];lines[count][length+1]=0;
                if(MeasureText(lines[count],fs)>panelW && length>0)break;
                if(next[length]==' ')space=length;
                length++;
            }
            if(next[length] && space>0)length=space;
            lines[count][length]=0;
            int lineWidth=MeasureText(lines[count],fs);if(lineWidth>width)width=lineWidth;
            next+=length;while(*next==' ')next++;
            count++;
        }
        int y=h-155;
        DrawRectangle(panelX-6,y-4,width+12,count*(fs+2)+6,Fade(BLACK,.8f));
        for(int n=0;n<count;n++)DrawText(lines[n],panelX,y+n*(fs+2),fs,GOLD);
    }
    /* Equipped gun costs are shown beside the two-slot equipment HUD. */
}

static bool couch_frame(void) {
    if(netTransport.role!=NET_ROLE_OFFLINE){
        if(gameMode!=GAME_MODE_FIREFIGHT||gameState!=GAME_STATE_GAMEOVER)return false;
        BeginDrawing();combat_ui_results();DrawText(netTransport.role==NET_ROLE_HOST?"Enter / A: new run   Esc / B: leave":"Waiting for host   Esc / B: leave",48,GetScreenHeight()-45,22,GOLD);EndDrawing();
        if(couch_back()){net_transport_shutdown(&netTransport);netRequestedRole=NET_ROLE_OFFLINE;netLobbyStarted=false;gameState=GAME_STATE_MENU;couchLobby=false;}
        else if(netTransport.role==NET_ROLE_HOST&&couch_confirm()){ResetGame();gameState=GAME_STATE_PLAYING;fp_network_reset();}
        return true;
    }
    bool confirm=couch_confirm(), back=couch_back(), pausePressed=IsKeyPressed(KEY_PAUSE);
    couchNavDelay-=GetFrameTime();
    if(couchAdvanced) {
        if(back) couchAdvanced=false;
        return false;
    }
    if(gameState==GAME_STATE_PLAYING) {
        for(int i=0;i<multiplayerActivePlayers;i++) if(multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD &&
            (!IsGamepadAvailable(seat_device(i))||couch_button(i,GAMEPAD_BUTTON_MIDDLE_RIGHT))) {
            gameState=GAME_STATE_PAUSED; couchSelection=0; break;
        }
        return false;
    }
    if(gameState!=GAME_STATE_MENU && gameState!=GAME_STATE_PAUSED && gameState!=GAME_STATE_GAMEOVER) return false;
    if(couchSettings) {
        couchSettingRow=(couchSettingRow+couch_axis(false)+4)%4;
        int change=couch_axis(true); LocalPlayerSeat *seat=&localSeats[couchSettingSeat];
        if(couchSettingRow==0) couchSettingSeat=(couchSettingSeat+change+4)%4;
        if(couchSettingRow==1) seat->sensitivity=clampf(seat->sensitivity+.1f*change,.5f,2);
        if(couchSettingRow==2) seat->deadzone=clampf(seat->deadzone+.01f*change,.05f,.3f);
        if(couchSettingRow==3 && (change||confirm)) seat->invert=!seat->invert;
        BeginDrawing();couch_title("CONTROLLER SETTINGS");
        couch_row(0,TextFormat("Player %d",couchSettingSeat+1),couchSettingRow==0);
        seat=&localSeats[couchSettingSeat];
        couch_row(1,TextFormat("Look sensitivity: %.1fx",seat->sensitivity),couchSettingRow==1);
        couch_row(2,TextFormat("Deadzone: %.2f",seat->deadzone),couchSettingRow==2);
        couch_row(3,TextFormat("Invert Y: %s",seat->invert?"Yes":"No"),couchSettingRow==3);
        DrawText("Left/right: change | Back: save",56,350,22,LIGHTGRAY);EndDrawing();
        if(back) {couch_save_settings();couchSettings=false;}
        return true;
    }
    if(gameState==GAME_STATE_PAUSED) {
        int missing=-1;
        for(int i=0;i<multiplayerActivePlayers;i++) if(multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD&&!IsGamepadAvailable(seat_device(i))) {missing=i;break;}
        if(missing>=0) for(int d=0;d<4;d++) if(IsGamepadAvailable(d)&&IsGamepadButtonPressed(d,GAMEPAD_BUTTON_MIDDLE_RIGHT)) {
            bool used=false;for(int j=0;j<multiplayerActivePlayers;j++) if(j!=missing&&multiplayerPlayerInput[j]==INPUT_TYPE_GAMEPAD&&seat_device(j)==d) used=true;
            if(!used) localSeats[missing].device=d;
        }
        couchSelection=(couchSelection+couch_axis(false)+3)%3;
        BeginDrawing();couch_title("PAUSED");
        couch_row(0,"Resume",couchSelection==0);couch_row(1,"Controller settings",couchSelection==1);couch_row(2,"Return to lobby",couchSelection==2);
        if(missing>=0) DrawText(TextFormat("P%d disconnected. Reconnect, or press Start on a free controller to reassign.",missing+1),48,310,20,GOLD);
        EndDrawing();
        if(confirm || pausePressed) {
            if(couchSelection==0&&missing<0) gameState=GAME_STATE_PLAYING;
            if(couchSelection==1) couchSettings=true;
            if(couchSelection==2) {if(enemyLab) enemy_lab_exit();else {gameState=GAME_STATE_MENU;couchLobby=true;memset(couchReady,0,sizeof(couchReady));}}
        }
        return true;
    }
    if(gameState==GAME_STATE_GAMEOVER) {
        couchSelection=(couchSelection+couch_axis(false)+2)%2;
        if(gameMode==GAME_MODE_FIREFIGHT){
            BeginDrawing();combat_ui_results();
            int y=GetScreenHeight()-65;
            DrawText(couchSelection==0?"> Retry":"  Retry",56,y,26,couchSelection==0?GOLD:LIGHTGRAY);
            DrawText(couchSelection==1?"> Lobby":"  Lobby",240,y,26,couchSelection==1?GOLD:LIGHTGRAY);
            DrawText(couch_menu_device()<0?"Arrows: select   Enter: confirm":"D-pad: select   A: confirm",56,y+34,16,LIGHTGRAY);
            EndDrawing();
            if(confirm){if(couchSelection==0)couch_launch();else {gameState=GAME_STATE_MENU;couchLobby=true;memset(couchReady,0,sizeof(couchReady));}}
            return true;
        }
        BeginDrawing();couch_title(gameMode==GAME_MODE_FIREFIGHT?(firefightWaveState==FIREFIGHT_STATE_VICTORY?"VICTORY":"FIREFIGHT DEFEAT"):TextFormat("PLAYER %d WINS",winnerId));
        couch_row(0,"Rematch",couchSelection==0);couch_row(1,"Lobby",couchSelection==1);
        DrawText(gameMode==GAME_MODE_FIREFIGHT?TextFormat("Wave %d | Team score %d",firefightWave,firefightScore):TextFormat("First to %d points",winningScore),56,230,24,GOLD);
        for(int i=0;i<multiplayerActivePlayers;i++) {
            const char *stats=gameMode==GAME_MODE_FIREFIGHT?
                TextFormat("P%d: %d kills | %d physics kills | %d revives",i+1,players[i].kills,players[i].debrisKills,reviveCount[i]):
                TextFormat("P%d: %d points | %d kills | %d deaths",i+1,player_points(&players[i]),players[i].kills,players[i].deaths);
            DrawText(stats,56,278+36*i,22,RAYWHITE);
        }
        EndDrawing();
        if(confirm) {
            if(couchSelection==0) couch_launch();
            else {gameState=GAME_STATE_MENU;couchLobby=true;memset(couchReady,0,sizeof(couchReady));}
        }
        return true;
    }
    if(couchLobby) {
        for(int d=0;d<4;d++) if(IsGamepadAvailable(d)&&IsGamepadButtonPressed(d,GAMEPAD_BUTTON_MIDDLE_RIGHT)) {
            bool used=false;for(int i=0;i<multiplayerActivePlayers;i++) if(multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD&&seat_device(i)==d) used=true;
            if(!used && multiplayerActivePlayers<4 && !firefightPractice) {int i=multiplayerActivePlayers++;multiplayerPlayerInput[i]=INPUT_TYPE_GAMEPAD;localSeats[i].device=d;couchReady[i]=false;}
        }
        if(IsKeyPressed(KEY_K) && multiplayerActivePlayers<2 && !firefightPractice) {int i=multiplayerActivePlayers++;multiplayerPlayerInput[i]=INPUT_TYPE_KEYBOARD;couchReady[i]=false;}
        int count=firefightPractice?1:multiplayerActivePlayers;
        bool all=true;
        for(int i=0;i<count;i++) {
            if(couch_button(i,GAMEPAD_BUTTON_RIGHT_FACE_DOWN)||(multiplayerPlayerInput[i]==INPUT_TYPE_KEYBOARD&&IsKeyPressed(i==0?KEY_ENTER:KEY_RIGHT_SHIFT))) couchReady[i]=!couchReady[i];
            if(i>0&&(couch_button(i,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT)||(multiplayerPlayerInput[i]==INPUT_TYPE_KEYBOARD&&IsKeyPressed(KEY_BACKSPACE)))) {
                for(int j=i;j<count-1;j++) {multiplayerPlayerInput[j]=multiplayerPlayerInput[j+1];localSeats[j]=localSeats[j+1];couchReady[j]=false;}
                multiplayerActivePlayers--;all=false;break;
            }
            all &= couchReady[i];
        }
        int change=couch_axis(true);
        if(change && gameMode==GAME_MODE_DEATHMATCH) {
            currentWorldType=(WorldType)((currentWorldType+change+WORLD_TYPE_COUNT)%WORLD_TYPE_COUNT);
            useCustomMap=false;memset(couchReady,0,sizeof(couchReady));all=false;
        }
        if(gameMode==GAME_MODE_DEATHMATCH && (IsKeyPressed(KEY_TAB)||couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_UP))) {
            useCustomMap=!useCustomMap;memset(couchReady,0,sizeof(couchReady));all=false;
        }
        if(change && gameMode==GAME_MODE_FIREFIGHT && !firefightPractice) {
            firefightEndless=!firefightEndless;firefightFoundry=true;memset(couchReady,0,sizeof(couchReady));all=false;
        }
        if(gameMode==GAME_MODE_FIREFIGHT && firefightEndless && (IsKeyPressed(KEY_TAB)||couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_UP))) {firefightFoundry=!firefightFoundry;memset(couchReady,0,sizeof(couchReady));all=false;}
        if(gameMode==GAME_MODE_FIREFIGHT && !firefightPractice && (IsKeyPressed(KEY_C)||couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_LEFT))){firefightCheckpoints=!firefightCheckpoints;memset(couchReady,0,sizeof(couchReady));all=false;}
        if(gameMode==GAME_MODE_FIREFIGHT && !firefightPractice && IsKeyPressed(KEY_H)){
            netRequestedRole=NET_ROLE_HOST;netRequestedCreative=false;netRequestedLocalPlayers=count;
            if(net_transport_host(&netTransport,netRequestedPort)){net_set_host_local_players(count);netLobbyStarted=false;couchLobby=false;gameState=GAME_STATE_LOBBY;return true;}
        }
        BeginDrawing();couch_title(firefightPractice?"LEARN COMBAT":gameMode==GAME_MODE_FIREFIGHT?"FIREFIGHT LOBBY":"LOCAL VERSUS LOBBY");
        couch_row(0,firefightPractice?"Harvest > expose > melee > throw":gameMode==GAME_MODE_FIREFIGHT?(firefightEndless?"Endless Survival (left/right changes mode)":"Twelve-wave Firefight (left/right changes mode)"):"Classic deathmatch",false);
        if(gameMode==GAME_MODE_FIREFIGHT)couch_row(1,firefightFoundry?"Arena: Foundry Holdout":"Arena: current world / custom map",false);
        else {
            const char *names[]={"Greek Temple","Megalith","Sun-Henge","Forerunner","Citadel","Test Arena","Blood Arena","Procedural","Crosscut"};
            couch_row(1,useCustomMap?TextFormat("Custom map slot %d",creativeMapSlot+1):TextFormat("Arena: %s | Left/right: change",names[currentWorldType]),false);
        }
        if(gameMode==GAME_MODE_FIREFIGHT && !firefightPractice)couch_row(2,TextFormat("Checkpoints: %s | X / C: toggle",firefightCheckpoints?"ON (retry current wave)":"OFF"),false);
        for(int i=0;i<count;i++) couch_row(i+2+(gameMode==GAME_MODE_FIREFIGHT&&!firefightPractice),TextFormat("P%d: %s %d | %s",i+1,multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD?"Controller":"Keyboard",multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD?seat_device(i)+1:i+1,couchReady[i]?"READY":"Not ready"),couchReady[i]);
        DrawText("Start: join | Bottom face / Enter / RShift: ready | K: add keyboard P2",48,465,18,LIGHTGRAY);
        DrawText(gameMode==GAME_MODE_DEATHMATCH?"Back: leave | Top face / Tab: use custom map":"Back: leave | Y / Tab: endless arena | H: host LAN",48,492,18,LIGHTGRAY);EndDrawing();
        if(back) {couchLobby=false;memset(couchReady,0,sizeof(couchReady));}
        else if(all) couch_launch();
        return true;
    }
    couchSelection=(couchSelection+couch_axis(false)+6)%6;
    BeginDrawing();couch_title("CONSTRUCTOR BRAWLER");clarity_menu();
    const char *items[]={"Firefight","Local versus","Learn combat","Controller settings","World builder / LAN","Enemy Lab"};
    for(int i=0;i<6;i++) couch_row(i,items[i],couchSelection==i);
    DrawText(couch_menu_device()>=0?"D-pad / stick: navigate":"Arrow keys: navigate",48,425,18,LIGHTGRAY);DrawText(couch_menu_device()>=0?"Bottom face: confirm":"Enter: confirm",48,450,18,LIGHTGRAY);EndDrawing();
    if(confirm) {
        if(couchSelection==0||couchSelection==1||couchSelection==2) {gameMode=couchSelection==1?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;couch_begin_lobby(couchSelection==2);}
        else if(couchSelection==3)couchSettings=true;
        else if(couchSelection==5)enemy_lab_start(ENEMY_TYPE_SCRAPPER);
        else {couchAdvanced=true;firefightPractice=false;gameMode=GAME_MODE_DEATHMATCH;activePlayers=clamp_active_players(multiplayerActivePlayers);}
    }
    return true;
}
