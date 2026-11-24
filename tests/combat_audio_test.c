#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[MAX_PLAYERS];

/* Capture actual stereo device frames, then silence the test output. */
static pthread_mutex_t captureLock=PTHREAD_MUTEX_INITIALIZER;
static double energy[2];
static void capture_mix(void *samples,unsigned frames) {
    float *out=samples;
    pthread_mutex_lock(&captureLock);
    for(unsigned i=0;i<frames*2;i++){energy[i%2]+=out[i]*out[i];out[i]=0;}
    pthread_mutex_unlock(&captureLock);
}
static double pan_ratio(float x,float yaw){
    players[0].yaw=yaw;
    WaitTime(.25);
    pthread_mutex_lock(&captureLock);energy[0]=energy[1]=0;pthread_mutex_unlock(&captureLock);
    clarityClock+=1;
    play_sfx_at(SFX_FIRE,(Vector3){x,2,x==0?-3:0},4,1,1);
    WaitTime(.2);
    pthread_mutex_lock(&captureLock);
    assert(energy[0]+energy[1]>0);
    double db=10*log10(energy[1]/energy[0]);
    pthread_mutex_unlock(&captureLock);
    return db;
}
int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(320,240,"Combat audio regression");
    init_sfx();assert(IsAudioDeviceReady());
    AttachAudioMixedProcessor(capture_mix);
    gameState=GAME_STATE_PLAYING;gameMode=GAME_MODE_FIREFIGHT;
    activePlayers=firefightHumans=1;playerInput[0]=INPUT_TYPE_KEYBOARD;
    players[0].pos=(Vector3){0,2,0};combatantLife[0]=LIFE_ALIVE;
    double right=pan_ratio(4,0),left=pan_ratio(-4,0),center=pan_ratio(0,0),turned=pan_ratio(4,180);
    printf("Right-minus-left energy, dB: right %.2f, left %.2f, front %.2f, turned %.2f\n",right,left,center,turned);
    assert(right>15 && left< -15 && fabs(center)<.1 && turned< -15);

    players[0].matter=80;
    apply_matter_damage(0,-1,5);
    float first=sfxLastPlay[SFX_PLAYER_HURT];assert(first>=0);
    apply_matter_damage(0,-1,5);assert(first==sfxLastPlay[SFX_PLAYER_HURT]);
    WaitTime(.16);players[0].invuln_timer=1;apply_matter_damage(0,-1,5);
    players[0].invuln_timer=0;eq_spend(0,1);assert(first==sfxLastPlay[SFX_PLAYER_HURT]);
    apply_matter_damage(0,-1,100);update_exposed_audio();
    assert(warningActive && first==sfxLastPlay[SFX_PLAYER_HURT]);
    double pulse=nextWarning;update_exposed_audio();assert(nextWarning==pulse);
    WaitTime(1.27);update_exposed_audio();assert(nextWarning>pulse);
    activePlayers=firefightHumans=2;playerInput[1]=INPUT_TYPE_KEYBOARD;
    combatantLife[1]=LIFE_ALIVE;players[1].isExposed=true;players[0].isExposed=false;
    pulse=nextWarning;update_exposed_audio();assert(nextWarning==pulse);
    combatantLife[1]=LIFE_DOWNED;update_exposed_audio();assert(!warningActive);
    netTransport.role=NET_ROLE_HOST;net_set_host_local_players(1);netPlayerPresent[1]=true;
    ClarityEffect hurt={.kind=CLARITY_PLAYER_HURT,.attacker=1};
    clarity_effect(hurt,false);assert(first==sfxLastPlay[SFX_PLAYER_HURT]);
    netTransport.role=NET_ROLE_CLIENT;netLocalPlayerSlots[0]=1;activePlayers=2;
    clarity_effect(hurt,false);assert(sfxLastPlay[SFX_PLAYER_HURT]>first);
    combatantLife[1]=LIFE_ALIVE;update_exposed_audio();assert(warningActive);
    gameState=GAME_STATE_MENU;update_exposed_audio();assert(!warningActive);
    netTransport.role=NET_ROLE_OFFLINE;
    uint32_t initial=clarityNextId;
    ai_sighting_audio(4,true,.01f);ai_sighting_audio(4,true,.01f);
    ai_sighting_audio(4,false,2.9f);ai_sighting_audio(4,true,.01f);
    assert(clarityNextId==initial+1);
    ai_sighting_audio(4,false,3);ai_sighting_audio(4,true,.01f);
    assert(clarityNextId==initial+2);
    gameState=GAME_STATE_PLAYING;activePlayers=6;firefightHumans=1;
    for(int i=1;i<4;i++)players[i].respawn_timer=9999;
    players[0].pos=(Vector3){0,2,-5};players[0].isExposed=false;
    players[4].pos=(Vector3){0,2,0};players[4].yaw=180;
    playerInput[4]=INPUT_TYPE_BOT_MEDIUM;combatantTeam[4]=1;
    players[5].respawn_timer=9999;ai_reset(4);
    initial=clarityNextId;ai_alert(4,0,players[0].pos,AI_ALERT_HEARD);
    ai_tick(4,.01f);assert(clarityNextId==initial);
    players[4].yaw=0;ai_tick(4,.01f);assert(clarityNextId==initial+1);
    WaitTime(.36);play_sfx_at(SFX_SIGHT_GROWL,players[4].pos,4,2,1);
    double firstGrowl=lastGrowl;assert(GetTime()-firstGrowl<.05);
    play_sfx_at(SFX_SIGHT_GROWL,players[4].pos,5,2,1);assert(lastGrowl==firstGrowl);
    activePlayers=firefightHumans=2;
    for(int i=0;i<2;i++){
        eq_reset(i);equipment[i].gun[0]=EQ_RAIL;equipment[i].charge=.35f;
        players[i].respawn_timer=0;players[i].matter=50;players[i].isExposed=false;
    }
    lab_update_audio();
    assert(IsAudioStreamPlaying(labChargeStream[0]) && IsAudioStreamPlaying(labChargeStream[1]));
    netTransport.role=NET_ROLE_CLIENT;netPlayerPresent[0]=netPlayerPresent[1]=true;
    equipment[1].gun[0]=EQ_CUTTER;equipment[1].held=equipment[1].beamContact=true;
    equipment[1].beamEnd=players[1].pos;
    lab_update_audio();
    assert(atomic_load(&labChargeLevel[0])>0 && !atomic_load(&labChargeLevel[1]));
    assert(IsSoundPlaying(labBeam[1]) && IsSoundPlaying(labContact[1]));
    equipment[1].held=false;gameState=GAME_STATE_MENU;WaitTime(.06);lab_update_audio();
    assert(!labBeam[1].frameCount && !labContact[1].frameCount && !IsAudioStreamPlaying(labChargeStream[0]));
    DetachAudioMixedProcessor(capture_mix);shutdown_sfx();CloseWindow();
    puts("Stereo, damage, warning, local ownership and growl rearm passed");
}
