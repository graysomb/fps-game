#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

static void fresh(void) {
    clear_world_voxels();activePlayers=0;netTransport.role=NET_ROLE_OFFLINE;
    memset(players,0,sizeof(players));worldPeriodic=false;init_static_hash();
}
static void advance(float seconds,float dt) {
    int n=(int)lroundf(seconds/dt);
    for(int i=0;i<n;i++)recycle_dead_voxels(dt);
}
static void check_rate(int fps) {
    float dt=1.f/fps;
    fresh();
    int v=addVoxel(4.25f,2.25f,4.25f,false,true,GRAY,0);assert(v>=0);
    voxels[v].owner=0;voxels[v].sleeping=true;
    advance(9,dt);assert(voxel_count==1);
    advance(2,dt);assert(voxel_count==0);
    assert(abs(recycleFrameCounter-660)<=1);

    fresh();v=add_static_voxel_at_grid(8,4,8,GRAY,0);assert(v>=0);voxels[v].owner=0;
    advance(29,dt);assert(voxel_count==1);
    advance(2,dt);assert(voxel_count==0);

    fresh();v=add_static_voxel_at_grid(8,4,8,GRAY,0);assert(v>=0);
    assert(recycle_queue_push(&voxels[v]));remove_voxel_index(v);
    advance(14,dt);assert(voxel_count==0 && !regeneration.active);
    advance(3,dt);assert(voxel_count==1 && recycleQueueCount==0);

    fresh();v=addVoxel(4.25f,2.25f,4.25f,false,true,GRAY,0);assert(v>=0);
    voxels[v].owner=0;voxels[v].restorationQueued=true;
    voxels[v].restorationDeadlineFrame=RECYCLE_STATIC_RESTORE_DELAY;
    activePlayers=1;players[0].tetherHolding=true;players[0].tetherVoxel=v;
    players[0].tetherVoxelIdentity=voxels[v].identity;tetherTag[v]=1;
    advance(20,dt);assert(voxel_count==1 && voxels[0].lifeFrames==0);
    players[0].tetherHolding=false;recycle_dead_voxels(dt);assert(voxel_count==0);
    printf("elapsed cleanup %d FPS: PASS\n",fps);
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Elapsed cleanup regression");
    check_rate(4);check_rate(30);check_rate(60);check_rate(144);
    fresh();
    for(int i=0;i<100;i++){recycle_dead_voxels(.004f);recycle_dead_voxels(.246f);}
    assert(abs(recycleFrameCounter-1500)<=1);
    int ticks=recycleFrameCounter,updates=recycleUpdateCounter;
    recycle_dead_voxels(0);assert(recycleFrameCounter==ticks && recycleUpdateCounter==updates+1);
    fresh();assert(recycleFrameCounter==0 && recycleTickRemainder==0);
    CloseWindow();puts("Variable frame time, pause, world reset: PASS");return 0;
}

