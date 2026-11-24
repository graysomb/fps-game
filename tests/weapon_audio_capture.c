// Run with "live", "original", or "prototype" and an output path for stereo f32le at 48 kHz.
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

// Sweep through open air and a voxel wall, fire diagonally across the listener's
// forward view, then alternate the two hitscan weapons at a fixed half-second cadence.
static void second_weapon_sequence(int tick,bool prototype)
{
    Equipment *e=&equipment[0];
    Player *p=&players[0];
    if(tick==0){
        init_static_hash();
        add_static_box_at_grid(-2,2,3,5,-40,-10,GRAY);
        init_static_hash(); rebuild_voxel_hash();
        p->matter=p->matterMax=100;
        e->gun[e->selected]=EQ_CUTTER;
    }
    if(tick<=140){
        p->yaw=tick>=54 && tick<110 ? 0 : 60;
        eq_command(0,&(NetInputCommand){.held=tick>=15 && tick<140 ? NET_INPUT_FIRE : 0},1.f/60);
    }
    if(tick==182){
        e->gun[e->selected]=EQ_ROCKET;
        p->yaw=-55;p->pitch=-8;
        assert(eq_fire(0));
    }
    if(tick>=182){p->yaw=0;p->pitch=0;}
    if(tick>=305 && tick<=515 && (tick-305)%30==0){
        e->gun[e->selected]=((tick-305)/30)%2 ? EQ_DMR : EQ_GLOCK;
        assert(eq_fire(0));
    }
    eq_tick(1.f/60);
    lab_update_audio();
    if(prototype){
        if(tick==42)assert(labBeam[0].frameCount && !labContact[0].frameCount);
        if(tick==66)assert(labContact[0].frameCount);
        if(tick==125)assert(!labContact[0].frameCount);
        if(tick==142)assert(!labBeam[0].frameCount);
        if(tick==194)assert(labMotors[0].frameCount);
        if(tick==285)assert(!labMotors[0].frameCount && labTake[LAB_BOOM]==1);
    }
}

// Alternating launch sounds share the same cadence. The sword contacts are
// resolved by its normal sweep, including the miss in the empty first scene.
static void third_weapon_sequence(int tick,bool prototype)
{
    Player *p=&players[0];
    Equipment *e=&equipment[0];
    if (tick==0) {
        init_static_hash();
        p->last_shot_time=p->last_melee_time=-10;
        p->matter=p->matterMax=100;
    }
    bool burst=tick>=18 && tick<=58 && (tick-18)%10==0;
    bool alternating=tick>=96 && tick<=246 && (tick-96)%30==0;
    if (burst || alternating) {
        bool voxel=alternating && ((tick-96)/30)%2;
        e->gun[e->selected]=voxel ? EQ_VOXEL : EQ_DEFAULT;
        p->dynamicShotActive=voxel;
        p->pitch=voxel ? -15 : 12;
        p->yaw=voxel ? -45 : 0;
        assert(FireVoxel(0));
        p->yaw=p->pitch=0;
    }
    if (tick<282 && voxel_count) {
        update_projectiles(1.f/60);
        simulate_voxel_pbd_steps(1.f/60,1);
        handle_pbd_projectile_hits();
    }
    if (tick==282 || tick==386 || tick==490) {
        clear_world_voxels();init_static_hash();
        if (tick==386) add_static_box_at_grid(-2,2,2,6,-5,-5,GRAY);
        if (tick==490) {
            activePlayers=5;
            for (int i=1;i<4;i++)players[i].respawn_timer=9999;
            enemy_lab_spawn(4,ENEMY_TYPE_STANDARD,(Vector3){0,2,-2.5f});
        }
        rebuild_voxel_hash();
        e->gun[e->selected]=EQ_BLADE;
        assert(eq_fire(0));
    }
    eq_sword_tick(1.f/60);
    if (prototype && tick==570) {
        assert(labTake[LAB_DEFAULT]==8 && labTake[LAB_VOXEL]==3 && labTake[LAB_FRACTURE]==3);
        assert(labTake[LAB_SWISH]==3 && labTake[LAB_SHEAR]==1 && labTake[LAB_STRIKE]==1);
    }
}

int main(int argc,char **argv)
{
    assert(argc==3 || (argc==4 && (!strcmp(argv[3],"second") || !strcmp(argv[3],"third"))));
    bool second=argc==4 && !strcmp(argv[3],"second");
    bool third=argc==4 && !second;
    int seconds=argc==4 ? 10 : SECONDS;
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(320,240,"Weapon Lab capture");
    SetRandomSeed(1835);
    srand(1835);
    if (third) {init_pbd_thread_pool();assert(initialize_physics_backend());}
    init_sfx();
    assert(IsAudioDeviceReady());
    gameState=GAME_STATE_PLAYING;
    gameMode=GAME_MODE_FIREFIGHT;
    activePlayers=firefightHumans=1;
    playerInput[0]=INPUT_TYPE_KEYBOARD;
    players[0].pos=(Vector3){0,2,0};
    combatantLife[0]=LIFE_ALIVE;
    eq_reset(0);
    bool live=!strcmp(argv[1],"live");
    enemyLab=equipmentLab=!live;
    enemyLabPaused=false;
    bool prototype=live || !strcmp(argv[1],"prototype");
    if(!live) debug_console_exec(prototype?"weaponaudio prototype":"weaponaudio original");
    assert(weapon_audio_enabled()==prototype);
    double started=GetTime();
    AttachAudioMixedProcessor(capture);
    for(int tick=0;tick<seconds*60;tick++){
        while(atomic_load(&recorded)<tick*(RATE/60)) WaitTime(.001);
        clarityClock=1+tick/60.f;
        if(third){third_weapon_sequence(tick,prototype);continue;}
        if(second){second_weapon_sequence(tick,prototype);continue;}
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
            if(n==1 && age==10)assert(!atomic_load(&labChargeLevel[0]));
        }
        ex_tick(1.f/60);
        if(tick>=924 && tick<1160 && (tick-924)%22==0)
            play_sfx_at(SFX_STEP_1+(tick/22)%3,(Vector3){tick%44<22?-3:3,2,-1},2,1,1);
        players[0].isExposed=tick>=1075 && tick<1160;
        update_exposed_audio();
    }
    while(atomic_load(&recorded)<seconds*RATE)WaitTime(.001);
    DetachAudioMixedProcessor(capture);
    assert(fabs(GetTime()-started-seconds)<.35);
    if(prototype && !second && !third)assert(labTake[LAB_SPLAT]==labTake[LAB_GLOB] && labTake[LAB_GLOB]==15);
    FILE *file=fopen(argv[2],"wb");assert(file);
    size_t samples=(size_t)seconds*RATE*2;
    assert(fwrite(reel,sizeof(float),samples,file)==samples);fclose(file);
    if (third) {clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();}
    shutdown_sfx();CloseWindow();
}
