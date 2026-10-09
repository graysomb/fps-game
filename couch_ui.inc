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
static bool couch_confirm(void) { return IsKeyPressed(KEY_ENTER)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_DOWN); }
static bool couch_back(void) { return IsKeyPressed(KEY_ESCAPE)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT); }
static int couch_axis(bool horizontal) {
    int minus=horizontal?KEY_LEFT:KEY_UP, plus=horizontal?KEY_RIGHT:KEY_DOWN;
    if(IsKeyPressed(minus)) return -1;
    if(IsKeyPressed(plus)) return 1;
    if(multiplayerPlayerInput[0]!=INPUT_TYPE_GAMEPAD || !IsGamepadAvailable(seat_device(0))) return 0;
    int d=seat_device(0);
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
    if(selected) DrawRectangle(40,y-6,GetScreenWidth()-80,42,(Color){40,77,94,255});
    DrawText(label,56,y,24,selected?GOLD:LIGHTGRAY);
}
static void couch_begin_lobby(bool practice) {
    if(gameMode==GAME_MODE_DEATHMATCH) {currentWorldType=WORLD_TYPE_CROSSCUT;useCustomMap=false;}
    firefightPractice=practice;
    couchLobby=true; couchSelection=0;
    multiplayerActivePlayers=1;
    memset(couchReady,0,sizeof(couchReady));
}
static void couch_launch(void) {
    if(firefightPractice) multiplayerActivePlayers=1;
    if(gameMode==GAME_MODE_FIREFIGHT && (firefightFoundry||firefightPractice)) useCustomMap=false;
    ResetGame(); gameState=GAME_STATE_PLAYING; couchLobby=false;
}
static void ff_draw_hud(int i,int w,int h) {
    if(gameMode!=GAME_MODE_FIREFIGHT || i>=firefightHumans) return;
    int font=w<700?16:20;
    /* Reserve the lower-left corner for the existing radar. */
    int panelX=160, panelW=w-panelX-12;
    DrawRectangle(panelX-6,h-154,panelW+6,100,Fade(BLACK,.8f));
    const char *text=NULL;
    if(combatantLife[i]==LIFE_DOWNED) text=TextFormat("DOWNED %.0fs | Revive %.1f / 3s",bleedTimer[i],reviveProgress[i]);
    else if(combatantLife[i]==LIFE_SPECTATING) text="Spectating until wave clear";
    else if(combatantLife[i]==LIFE_WAITING) text=TextFormat("Respawn in %.0fs",players[i].respawn_timer);
    else if(firefightPractice) {
        const char *steps[]={"Melee a block to harvest Matter", "Shoot the enemy until its shield breaks", "Finish the exposed enemy with melee", "Expose the next enemy; throw debris to finish", "Training complete!"};
        text=steps[practiceStep<5?practiceStep:4];
    } else if(firefightWaveState==FIREFIGHT_STATE_INTERMISSION) {
        int v=prepVotes[i]<0?0:prepVotes[i];
        text=TextFormat("%s %s",ff_upgrade_names[v],prepReady[i]?"[READY]":"");
        DrawText(playerInput[i]==INPUT_TYPE_GAMEPAD?(prepMenuOpen[i]?"D-pad: vote | Bottom: ready | Top: build":"Top face: vote menu"):(prepMenuOpen[i]?"Look: vote | Jump: ready | V/N: build":"V/N: vote menu"),panelX,h-112,font,RAYWHITE);
        const WaveDefinition *next=&ff_waves[firefightWave%5];
        DrawText(TextFormat("Next: %s (%d standard, %d swarm, %d giants)",next->name,ff_count(firefightWave+1,firefightHumans,0),ff_count(firefightWave+1,firefightHumans,1),ff_count(firefightWave+1,firefightHumans,2)),panelX,h-142,w<700?14:font,RAYWHITE);
    } else {
        for(int j=0;j<firefightHumans;j++) if(combatantLife[j]==LIFE_DOWNED) {
            text=TextFormat("P%d DOWN | Hold %s nearby: revive (25 Matter)",j+1,playerInput[i]==INPUT_TYPE_GAMEPAD?"LB":i==0?"X":","); break;
        }
        if(!text && players[i].matter<1) text="No Matter: melee blocks to harvest";
        if(!text) {
            Vector3 dir; float angle;
            if(find_aim_assist_target(i,&dir,&angle) && angle<12) {
                for(int j=MAX_PLAYERS;j<MAX_COMBATANTS;j++) if(players[j].respawn_timer<=0 && players[j].isExposed) {
                    Vector3 delta=v_norm(v_sub(players[j].pos,players[i].pos));
                    if(v_dot(delta,player_forward(&players[i]))>.97f) {
                        text=playerInput[i]==INPUT_TYPE_GAMEPAD?"EXPOSED: melee (bottom face) or throw (LT)":i==0?"EXPOSED: melee Z or throw R":"EXPOSED: melee M or throw P"; break;
                    }
                }
            }
        }
        if(!text && firefightEnemiesSpawned<firefightWaveEnemiesTotal) text=TextFormat("Incoming: %s entrance",ffEntrance?"EAST":"WEST");
    }
    if(text) {
        int fs=font;
        while(fs>12 && MeasureText(text,fs)>panelW) fs--;
        DrawText(text,panelX,h-82,fs,GOLD);
    }
    if(combatantLife[i]==LIFE_ALIVE && firefightWaveState!=FIREFIGHT_STATE_INTERMISSION)
        DrawText(TextFormat("Shot 1 | Build %d | Harvest +%d",ff_bonus(i,UPGRADE_BUILD)?5:10,ff_bonus(i,UPGRADE_HARVEST)?15:10),panelX,h-112,font,RAYWHITE);
}

static bool couch_frame(void) {
    if(netTransport.role!=NET_ROLE_OFFLINE) return false;
    bool confirm=couch_confirm(), back=couch_back(), pausePressed=IsKeyPressed(KEY_PAUSE);
    couchNavDelay-=GetFrameTime();
    /* Controller-only startup: Start assigns the first seat. */
    if(gameState==GAME_STATE_MENU && !couchLobby) for(int d=0;d<4;d++)
        if(IsGamepadAvailable(d)&&IsGamepadButtonPressed(d,GAMEPAD_BUTTON_MIDDLE_RIGHT)) {
            multiplayerPlayerInput[0]=INPUT_TYPE_GAMEPAD;localSeats[0].device=d;
        }
    if(couchAdvanced) {
        if(IsKeyPressed(KEY_ESCAPE)) couchAdvanced=false;
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
            if(couchSelection==2) {gameState=GAME_STATE_MENU;couchLobby=true;memset(couchReady,0,sizeof(couchReady));}
        }
        return true;
    }
    if(gameState==GAME_STATE_GAMEOVER) {
        couchSelection=(couchSelection+couch_axis(false)+2)%2;
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
        if(gameMode==GAME_MODE_DEATHMATCH && (IsKeyPressed(KEY_TAB)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_UP))) {
            useCustomMap=!useCustomMap;memset(couchReady,0,sizeof(couchReady));all=false;
        }
        if(change && gameMode==GAME_MODE_FIREFIGHT && !firefightPractice) {
            firefightEndless=!firefightEndless;firefightFoundry=true;memset(couchReady,0,sizeof(couchReady));all=false;
        }
        if(gameMode==GAME_MODE_FIREFIGHT && firefightEndless && (IsKeyPressed(KEY_TAB)||couch_button(0,GAMEPAD_BUTTON_RIGHT_FACE_UP))) {firefightFoundry=!firefightFoundry;memset(couchReady,0,sizeof(couchReady));all=false;}
        BeginDrawing();couch_title(firefightPractice?"LEARN COMBAT":gameMode==GAME_MODE_FIREFIGHT?"FIREFIGHT LOBBY":"LOCAL VERSUS LOBBY");
        couch_row(0,firefightPractice?"Harvest > expose > melee > throw":gameMode==GAME_MODE_FIREFIGHT?(firefightEndless?"Endless Survival (left/right changes mode)":"Five-wave Firefight (left/right changes mode)"):"Classic deathmatch",false);
        if(gameMode==GAME_MODE_FIREFIGHT)couch_row(1,firefightFoundry?"Arena: Foundry Holdout":"Arena: current world / custom map",false);
        else {
            const char *names[]={"Greek Temple","Megalith","Sun-Henge","Forerunner","Citadel","Test Arena","Blood Arena","Procedural","Crosscut"};
            couch_row(1,useCustomMap?TextFormat("Custom map slot %d",creativeMapSlot+1):TextFormat("Arena: %s | Left/right: change",names[currentWorldType]),false);
        }
        for(int i=0;i<count;i++) couch_row(i+2,TextFormat("P%d: %s %d | %s",i+1,multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD?"Controller":"Keyboard",multiplayerPlayerInput[i]==INPUT_TYPE_GAMEPAD?seat_device(i)+1:i+1,couchReady[i]?"READY":"Not ready"),couchReady[i]);
        DrawText("Start: join | Bottom face / Enter / RShift: ready | K: add keyboard P2",48,440,20,LIGHTGRAY);
        DrawText(gameMode==GAME_MODE_DEATHMATCH?"Back: leave | Top face / Tab: use custom map":"Back: leave | Endless: top face / Tab changes arena",48,470,20,LIGHTGRAY);EndDrawing();
        if(back) {couchLobby=false;memset(couchReady,0,sizeof(couchReady));}
        else if(all) couch_launch();
        return true;
    }
    couchSelection=(couchSelection+couch_axis(false)+5)%5;
    BeginDrawing();couch_title("CONSTRUCTOR BRAWLER");
    const char *items[]={"Firefight","Local versus","Learn combat","Controller settings","World builder / LAN / advanced setup"};
    for(int i=0;i<5;i++) couch_row(i,items[i],couchSelection==i);
    DrawText("Start on a controller to use it | D-pad: navigate | Bottom face: confirm",48,405,20,LIGHTGRAY);EndDrawing();
    if(confirm) {
        if(couchSelection==0||couchSelection==1||couchSelection==2) {gameMode=couchSelection==1?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;couch_begin_lobby(couchSelection==2);}
        else if(couchSelection==3)couchSettings=true;
        else {couchAdvanced=true;firefightPractice=false;gameMode=GAME_MODE_DEATHMATCH;activePlayers=clamp_active_players(multiplayerActivePlayers);}
    }
    return true;
}
