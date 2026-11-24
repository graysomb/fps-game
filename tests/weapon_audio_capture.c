// Run with "original" or "prototype" and an output path for stereo f32le at 48 kHz.
#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[MAX_PLAYERS];
enum { RATE=48000, SECONDS=20 };
static float reel[SECONDS*RATE*2];
static atomic_uint recorded;

static void capture(void *samples,unsigned frames)
{
    unsigned start=atomic_load(&recorded);
    unsigned remaining=SECONDS*RATE-start;
    unsigned count=frames<remaining?frames:remaining;
    memcpy(reel+2*start,samples,count*2*sizeof(float));
    memset(samples,0,frames*2*sizeof(float));
    atomic_store(&recorded,start+count);
}

int main(int argc,char **argv)
{
    assert(argc==3);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(320,240,"Weapon Lab capture");
    SetRandomSeed(1835);
    srand(1835);
    init_sfx();
    assert(IsAudioDeviceReady());
    gameState=GAME_STATE_PLAYING;
    gameMode=GAME_MODE_FIREFIGHT;
    activePlayers=firefightHumans=1;
    playerInput[0]=INPUT_TYPE_KEYBOARD;
    players[0].pos=(Vector3){0,2,0};
    combatantLife[0]=LIFE_ALIVE;
    eq_reset(0);
    enemyLab=equipmentLab=true;
    enemyLabPaused=false;
    bool prototype=!strcmp(argv[1],"prototype");
    debug_console_exec(prototype?"weaponaudio prototype":"weaponaudio original");
    assert(labPrototype==prototype);
    double started=GetTime();
    AttachAudioMixedProcessor(capture);
    for(int tick=0;tick<SECONDS*60;tick++){
        while(atomic_load(&recorded)<tick*(RATE/60)) WaitTime(.001);
        clarityClock=1+tick/60.f;
        Equipment *e=&equipment[0];
        bool shotgun=tick==24 || (tick>=384 && tick<=480 && (tick-384)%48==0) || tick==990;
        bool glob=tick==302 || (tick>=780 && tick<870 && (tick-780)%10==0) ||
                  (tick>=1032 && tick<1082 && (tick-1032)%10==0);
        if(shotgun||glob){
            int gun=shotgun?EQ_SHOTGUN:EQ_GLOB;
            e->gun[e->selected]=gun;
            e->cooldown[gun]=0;players[0].matter=100;
            players[0].pitch=-14;
            players[0].yaw=tick%20<10?-55:55;
            assert(eq_fire(0));
            players[0].yaw=0;
        }
        // Two incomplete pulls (released early, then underfunded), then three full shots.
        const int pulls[]={84,122,182,540,690};
        for(int n=0;n<5;n++){
            int age=tick-pulls[n],length=n<2?18:40;
            if(age<0 || age>length)continue;
            if(age==0){
                e->gun[e->selected]=EQ_RAIL;e->cooldown[EQ_RAIL]=0;
                players[0].matter=n==1?2:100;players[0].isExposed=false;
            }
            unsigned before=labTake[LAB_RAIL];
            eq_command(0,&(NetInputCommand){.held=age<length?NET_INPUT_FIRE:0},1.f/60);
            lab_update_charge();
            if(prototype && age==length)assert(labTake[LAB_RAIL]==before+(n>=2));
            if(n==1 && age==10)assert(!atomic_load(&labChargeLevel));
        }
        ex_tick(1.f/60);
        if(tick>=924 && tick<1160 && (tick-924)%22==0)
            play_sfx_at(SFX_STEP_1+(tick/22)%3,(Vector3){tick%44<22?-3:3,2,-1},2,1,1);
        players[0].isExposed=tick>=1075 && tick<1160;
        update_exposed_audio();
    }
    while(atomic_load(&recorded)<SECONDS*RATE)WaitTime(.001);
    DetachAudioMixedProcessor(capture);
    assert(fabs(GetTime()-started-SECONDS)<.35);
    if(prototype)assert(labTake[LAB_SPLAT]==labTake[LAB_GLOB] && labTake[LAB_GLOB]==15);
    FILE *file=fopen(argv[2],"wb");assert(file);
    assert(fwrite(reel,sizeof(float),SECONDS*RATE*2,file)==SECONDS*RATE*2);fclose(file);
    shutdown_sfx();CloseWindow();
}
