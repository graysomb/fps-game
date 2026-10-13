#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

static int find_identity(uint64_t identity) {
    for (int i=0;i<voxel_count;i++) if (voxels[i].identity==identity) return i;
    return -1;
}
static void fire_and_check(bool dynamic) {
    Player *p=&players[0];
    p->dynamicShotActive=dynamic;
    p->last_shot_time=-100;
    float matter=p->matter;
    int before=voxel_count;
    FireVoxel(0);
    assert(voxel_count==before+1);
    uint64_t identity=voxels[before].identity;
    fprintf(stderr,"shot dynamic=%d matter=%.0f->%.0f queued=%d deadline=%d frame=%d\n",
            dynamic,matter,p->matter,voxels[before].restorationQueued,
            voxels[before].restorationDeadlineFrame,recycleFrameCounter);
    assert(p->matter==matter-MATTER_SHOT_COST);
    Vector3 start=voxels[before].pos;
    update_projectiles(1.0f/60);
    simulate_voxel_pbd_steps(1.0f/120,2);
    recycle_dead_voxels(1.0f/60);
    int shot=find_identity(identity);
    if(shot<0){fprintf(stderr,"FAIL: launcher fired, but fresh cube was recycled before rendering\n");exit(1);}
    assert(!voxels[shot].restorationQueued && voxels[shot].restorationDeadlineFrame==0);
    assert(voxels[shot].pos.z<start.z-.1f);
    voxel_render_frame++;
    prepare_dynamic_voxel_transforms();
    bool visible=false;
    if(dynamic) {
        for(int n=0;n<instanceTransformsCount;n++) {
            float *t=instanceTransforms[n].v;
            visible |= v_length(v_sub((Vector3){t[12],t[13],t[14]},voxels[shot].pos))<.01f;
        }
    } else {
        for(int n=0;n<bullet_voxel_count;n++)visible |= bullet_voxel_indices[n]==shot;
    }
    assert(visible);
    remove_voxel_index(shot);
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Launcher recycling regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    gameMode=GAME_MODE_FIREFIGHT;firefightHumans=1;activePlayers=1;
    clear_world_voxels();
    players[0]=(Player){.pos={0,5,0},.matter=100,.matterMax=100};
    // A real terrain activation creates a queued-original deadline. Expiring
    // that debris leaves its former dense-array slot available for the gun.
    int source=add_static_voxel_at_grid(20,5,20,GRAY,0);assert(source>=0);
    assert(activate_static_voxel_for_tether(source,0,ACTIVATION_TETHER_BELIEF));
    assert(voxel_count==1 && voxels[0].restorationQueued);
    recycleFrameCounter=voxels[0].restorationDeadlineFrame;
    recycle_dead_voxels(1.0f/60);
    assert(voxel_count==0);
    // Before the fix every shot reused the expired deadline and vanished.
    for(int shot=0;shot<8;shot++)fire_and_check(true);
    fire_and_check(false);
    // Terrain descendants must still retain their intentional deadline.
    int other=add_static_voxel_at_grid(22,5,20,GRAY,0);assert(other>=0);
    assert(!voxels[other].restorationQueued && voxels[other].restorationDeadlineFrame==0);
    assert(activate_static_voxel_for_tether(other,0,ACTIVATION_TETHER_BELIEF));
    int debris=voxel_count-1;uint64_t id=voxels[debris].identity;
    assert(voxels[debris].restorationQueued && voxels[debris].restorationDeadlineFrame>recycleFrameCounter);
    recycleFrameCounter=voxels[debris].restorationDeadlineFrame;
    recycle_dead_voxels(1.0f/60);assert(find_identity(id)<0);
    // Repeat on the actual Firefight map, with its resident supply physics.
    firefightFoundry=true;firefightPractice=false;enemyLab=false;
    multiplayerActivePlayers=1;multiplayerPlayerInput[0]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;activePlayers=1;
    players[0].pos=(Vector3){0,15,0};players[0].yaw=players[0].pitch=0;
    players[0].matter=100;players[0].isExposed=false;
    source=table_get_static_only(41,5,-8);assert(source>=0);
    assert(activate_static_voxel_for_tether(source,0,ACTIVATION_TETHER_BELIEF));
    assert(voxels[voxel_count-1].restorationQueued);
    recycleFrameCounter=voxels[voxel_count-1].restorationDeadlineFrame;
    recycle_dead_voxels(1.0f/60);
    for(int shot=0;shot<8;shot++)fire_and_check(true);
    fire_and_check(false);
    shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    puts("Launcher recycling regression passed");return 0;
}
