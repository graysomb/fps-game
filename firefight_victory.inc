/* Local presentation flow: GAMEOVER remains the authoritative match state.
 * No epilogue state, input or presentation geometry is sent over the network. */
enum FfEndPhase { FFEND_NONE,FFEND_BANNER,FFEND_RESULTS,FFEND_EPILOGUE };
static struct {
    enum FfEndPhase phase;
    uint32_t generation;
    float runSeconds,banner,time,skip,fade;
    int selection,humans,score,kills,wave;
    int playerKills[4],physicsKills[4],revives[4];
    Color colors[4];
    bool checkpoints,leaving;
    Music music;
    Texture2D backdrop;
    RenderTexture2D scenes[2];
    WeaponBloom bloom;
    const char *error;
} ffend;
static void ffend_release_movie(void){
    if(IsMusicValid(ffend.music)){StopMusicStream(ffend.music);UnloadMusicStream(ffend.music);}
    ffend.music=(Music){0};
    for(int i=0;i<2;i++){if(ffend.scenes[i].id)UnloadRenderTexture(ffend.scenes[i]);ffend.scenes[i]=(RenderTexture2D){0};}
    weapon_bloom_unload(&ffend.bloom);victory_town_shutdown();
    ffend.time=ffend.skip=ffend.fade=0;ffend.leaving=false;
}
static void ffend_reset(void){
    ffend_release_movie();
    if(ffend.backdrop.id)UnloadTexture(ffend.backdrop);
    memset(&ffend,0,sizeof(ffend));ffend.generation=fpGeneration;
}
static bool ffend_eligible(void){
    return gameMode==GAME_MODE_FIREFIGHT && gameState==GAME_STATE_GAMEOVER &&
        firefightWaveState==FIREFIGHT_STATE_VICTORY && firefightWave>=FF_WAVE_COUNT &&
        !firefightPractice && !firefightEndless && !enemyLab;
}
static void ffend_capture_results(void){
    ffend.humans=clampi(firefightHumans,1,4);ffend.wave=firefightWave;
    ffend.score=firefightScore;ffend.kills=firefightKillsTotal;ffend.checkpoints=firefightCheckpoints;
    for(int i=0;i<ffend.humans;i++){
        ffend.playerKills[i]=players[i].kills;ffend.physicsKills[i]=players[i].debrisKills;
        ffend.revives[i]=reviveCount[i];ffend.colors[i]=player_palette_color(i);
    }
}
static void ffend_snapshot(void){
    ffend_capture_results();
    Image image=LoadImageFromScreen();if(image.data){ffend.backdrop=LoadTextureFromImage(image);UnloadImage(image);}
    ffend.phase=FFEND_BANNER;ffend.banner=0;ffend.selection=0;
}
static bool ffend_targets(void){
    int w=GetScreenWidth(),h=GetScreenHeight();if(w<=0||h<=0)return false;
    if(ffend.scenes[0].texture.width==w&&ffend.scenes[0].texture.height==h)return true;
    for(int i=0;i<2;i++){
        if(ffend.scenes[i].id)UnloadRenderTexture(ffend.scenes[i]);
        ffend.scenes[i]=LoadRenderTexture(w,h);
    }
    return ffend.scenes[0].id&&ffend.scenes[1].id;
}
static bool ffend_start_movie(void){
    ffend.error=NULL;
    if(!victory_town_init()||!ffend_targets()){
        ffend_release_movie();ffend.error="Epilogue assets unavailable. You can still play again or leave.";return false;
    }
    for(int i=0;i<ffend.humans;i++)townPlayerColors[i]=ffend.colors[i];
    duck_renderer_init(NULL);
    if(IsAudioDeviceReady()){
        ffend.music=LoadMusicStream(TOWN_MUSIC_PATH);
        if(IsMusicValid(ffend.music)){ffend.music.looping=false;SetMusicVolume(ffend.music,.75f);PlayMusicStream(ffend.music);}
    }
    ffend.time=ffend.skip=ffend.fade=0;ffend.leaving=false;ffend.phase=FFEND_EPILOGUE;return true;
}
static void ffend_render_scene(int index,TownSample sample){
    RenderTexture2D target=ffend.scenes[index];
    BeginTextureMode(target);ClearBackground((Color){19,31,43,255});town_sky(sample);
    BeginMode3D(sample.camera);victory_town_draw(sample,false);EndMode3D();EndTextureMode();
    if(weaponBloomEnabled&&weapon_bloom_resize(&ffend.bloom,target.texture.width,target.texture.height)){
        weapon_bloom_begin(&ffend.bloom,target);BeginMode3D(sample.camera);victory_town_draw(sample,true);EndMode3D();
        weapon_bloom_composite(&ffend.bloom,target);EndTextureMode();
    }
}
static void ffend_draw_movie(float dt){
    ffend.time+=dt;
    if(IsMusicValid(ffend.music)){
        UpdateMusicStream(ffend.music);
        if(IsMusicStreamPlaying(ffend.music))ffend.time=GetMusicTimePlayed(ffend.music);
    }
    int device=couch_menu_device();
    bool held=IsKeyDown(KEY_ESCAPE)||(device>=0&&IsGamepadButtonDown(device,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT));
    ffend.skip=held?ffend.skip+dt:0;
    if(ffend.skip>=.65f||ffend.time>=TOWN_MUSIC_SECONDS-.4f)ffend.leaving=true;
    if(ffend.leaving)ffend.fade=fminf(1,ffend.fade+dt/.4f);
    if(IsMusicValid(ffend.music))SetMusicVolume(ffend.music,.75f*(1-ffend.fade));
    if(!ffend_targets()){ffend_release_movie();ffend.phase=FFEND_RESULTS;return;}
    float time=fminf(ffend.time,TOWN_MUSIC_SECONDS-.001f)/town_playback_scale(),local;
    int scene=town_timeline(time,&local);
    float blend=scene==SC_COUNT-1?0:fmaxf(0,local-(townScenes[scene].duration-1));
    ffend_render_scene(0,victory_town_sample(scene,local));
    if(blend>0)ffend_render_scene(1,victory_town_sample((scene+1)%SC_COUNT,0));
    int w=GetScreenWidth(),h=GetScreenHeight();
    BeginDrawing();ClearBackground(BLACK);
    for(int i=0;i<(blend>0?2:1);i++)DrawTextureRec(ffend.scenes[i].texture,(Rectangle){0,0,w,-h},(Vector2){0},Fade(WHITE,i?blend:1));
    town_credits_draw(time,w,h);
    const char *hint=device>=0?"Hold B to skip":"Hold Escape to skip";
    int fs=h<600?14:18,x=w-MeasureText(hint,fs)-24;
    DrawRectangle(x-8,18,MeasureText(hint,fs)+16,fs+12,Fade(BLACK,.55f));DrawText(hint,x,24,fs,RAYWHITE);
    if(ffend.skip>0)DrawRectangle(x,fs+32,(int)(MeasureText(hint,fs)*fminf(ffend.skip/.65f,1)),3,GOLD);
    if(ffend.fade>0)DrawRectangle(0,0,w,h,Fade(BLACK,ffend.fade));
    EndDrawing();
    if(ffend.fade>=1){ffend_release_movie();ffend.phase=FFEND_RESULTS;ffend.selection=0;}
}
static void ffend_leave(void){
    ffend_reset();
    if(netTransport.role!=NET_ROLE_OFFLINE){net_transport_shutdown(&netTransport);netRequestedRole=NET_ROLE_OFFLINE;netLobbyStarted=false;}
    gameState=GAME_STATE_MENU;couchLobby=false;couchAdvanced=false;
}
static void ffend_draw_results(void){
    float s=fminf(GetScreenWidth()/1280.f,GetScreenHeight()/720.f);int x=(int)(48*s),fs=(int)(26*s);if(fs<14)fs=14;
    BeginDrawing();ClearBackground((Color){15,18,26,255});DrawRectangle(0,0,GetScreenWidth(),8,GOLD);
    DrawText("FIREFIGHT COMPLETE",x,(int)(35*s),(int)(44*s),GOLD);
    DrawText("The arena is yours. A new home begins.",x,(int)(100*s),fs,RAYWHITE);
    int seconds=(int)ffend.runSeconds;
    DrawText(TextFormat("Wave %d  |  Score %d  |  %d enemies defeated",ffend.wave,ffend.score,ffend.kills),x,(int)(150*s),fs,RAYWHITE);
    DrawText(TextFormat("%s %02d:%02d  |  Wave-scaled difficulty  |  Checkpoints %s",netTransport.role==NET_ROLE_CLIENT?"Time in run":"Completion time",seconds/60,seconds%60,ffend.checkpoints?"ON":"OFF"),x,(int)(190*s),fs-2,LIGHTGRAY);
    for(int i=0;i<ffend.humans;i++){
        const char *row=netTransport.role==NET_ROLE_CLIENT?
            TextFormat("P%d   %d kills   %d physics kills",i+1,ffend.playerKills[i],ffend.physicsKills[i]):
            TextFormat("P%d   %d kills   %d physics kills   %d revives",i+1,ffend.playerKills[i],ffend.physicsKills[i],ffend.revives[i]);
        DrawText(row,x,(int)((245+i*43)*s),fs,ffend.colors[i]);
    }
    const char *labels[]={"Watch epilogue",netTransport.role==NET_ROLE_CLIENT?"Play again (host starts next run)":"Play again","Main menu"};
    for(int i=0;i<3;i++){
        int y=GetScreenHeight()-(int)((175-i*42)*s);
        if(i==ffend.selection)DrawRectangle(x-8,y-4,GetScreenWidth()-2*x+16,fs+8,(Color){40,65,72,255});
        DrawText(labels[i],x,y,fs,i==ffend.selection?GOLD:RAYWHITE);
    }
    if(ffend.error)DrawText(ffend.error,x,(int)(440*s),fs-2,ORANGE);
    DrawText(couch_menu_device()<0?"Up/down: select   Enter: confirm   Escape: main menu":"D-pad: select   A: confirm   B: main menu",x,GetScreenHeight()-(int)(30*s),fs-4,LIGHTGRAY);
    EndDrawing();
}
static bool ffend_frame(float dt){
    dt=clampf(dt,0,.25f);
    if(ffend.generation!=fpGeneration)ffend_reset();
    if(!ffend_eligible()){
        if(ffend.phase!=FFEND_NONE)ffend_reset();
        if(gameMode==GAME_MODE_FIREFIGHT&&gameState==GAME_STATE_PLAYING)ffend.runSeconds+=dt;
        return false;
    }
    if(netTransport.role==NET_ROLE_CLIENT&&!netTransport.connected){ffend_leave();return false;}
    couchNavDelay-=dt;
    if(ffend.phase==FFEND_NONE)ffend_snapshot();
    if(ffend.phase==FFEND_BANNER){
        ffend.banner+=dt;int w=GetScreenWidth(),h=GetScreenHeight();
        BeginDrawing();ClearBackground((Color){15,18,26,255});
        if(ffend.backdrop.id)DrawTexturePro(ffend.backdrop,(Rectangle){0,0,ffend.backdrop.width,ffend.backdrop.height},(Rectangle){0,0,w,h},(Vector2){0},0,WHITE);
        DrawRectangle(0,0,w,h,Fade(BLACK,.45f));int fs=w<800?30:52;
        DrawText("FIREFIGHT COMPLETE",(w-MeasureText("FIREFIGHT COMPLETE",fs))/2,h/2-30,fs,GOLD);EndDrawing();
        if(ffend.banner>=3){ffend_capture_results();if(ffend.backdrop.id)UnloadTexture(ffend.backdrop);ffend.backdrop=(Texture2D){0};ffend.phase=FFEND_RESULTS;}
        return true;
    }
    if(ffend.phase==FFEND_EPILOGUE){ffend_draw_movie(dt);return true;}
    ffend.selection=(ffend.selection+couch_axis(false)+3)%3;
    bool confirm=couch_confirm(),back=couch_back();ffend_draw_results();
    if(back){ffend_leave();return true;}
    if(confirm){
        if(ffend.selection==0)ffend_start_movie();
        else if(ffend.selection==2)ffend_leave();
        else if(netTransport.role==NET_ROLE_CLIENT)ffend.error="Waiting for the host to start the next run.";
        else if(netTransport.role==NET_ROLE_HOST){ResetGame();gameState=GAME_STATE_PLAYING;fp_network_reset();}
        else couch_launch();
    }
    return true;
}
