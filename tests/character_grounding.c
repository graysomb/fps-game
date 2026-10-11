#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void expect_near(float actual,float expected) {if(fabsf(actual-expected)>.005f){fprintf(stderr,"got %f expected %f\n",actual,expected);abort();}}
static void reset_human(Vector3 pos) {
    players[0]=(Player){.pos=pos,.onGround=true,.matter=100,.matterMax=100,.tetherVoxel=-1};
    playerInput[0]=INPUT_TYPE_KEYBOARD;activePlayers=1;netTransport.role=NET_ROLE_OFFLINE;
}
static void cell(int x,int y,int z) {assert(add_static_voxel_at_grid(x,y,z,GRAY,0)>=0);}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,480,"Grounding checks");assert(IsWindowReady());
    init_static_hash();init_pbd_thread_pool();assert(initialize_physics_backend());
    clear_world_voxels();reset_human((Vector3){0,1,0});
    BoundingBox b=combatant_bounds(0,.25f);expect_near(b.min.y,0);expect_near(b.max.y,1.12f);expect_near(b.max.x,.39f);
    assert(segment_intersects_aabb((Vector3){0,.5f,-2},(Vector3){0,.5f,2},b.min,b.max));
    assert(!segment_intersects_aabb((Vector3){0,1.25f,-2},(Vector3){0,1.25f,2},b.min,b.max));
    CharacterPose pose=character_pose(players[0].pos,0,0,0,true,false,SKYBLUE,0);
    expect_near(Vector3Transform((Vector3){0,.287f,0},character_body_transform(pose)).y,1);
    CharacterAnimation animation={0};
    for(int frame=0;frame<120;frame++) {
        pose=character_pose(players[0].pos,0,0,0,true,false,SKYBLUE,frame/60.0f);
        character_animate(&animation,&pose,1.0f/60,character_ground_query,0);
        expect_near(pose.soles[0].y,0);expect_near(pose.soles[1].y,0);
    }
    for(int frame=0;frame<120;frame++) {
        Vector3 before[2]={animation.soles[0],animation.soles[1]};bool stance[2]={!animation.swinging[0],!animation.swinging[1]};
        players[0].pos.x+=.015f;
        pose=character_pose(players[0].pos,0,0,.9f,true,false,SKYBLUE,frame/60.0f);
        character_animate(&animation,&pose,1.0f/60,character_ground_query,0);
        for(int f=0;f<2;f++)if(stance[f]&&!animation.swinging[f])assert(Vector3Distance(before[f],pose.soles[f])<.02f);
    }
    // Stop, turn, reverse and teleport: feet remain finite and within leg reach.
    for(int frame=0;frame<240;frame++) {
        if(frame>=120)players[0].pos.x-=.02f;
        float yaw=fminf(180,frame*2);
        pose=character_pose(players[0].pos,yaw,0,frame>=120?1.2f:0,true,false,SKYBLUE,frame/60.0f);
        character_animate(&animation,&pose,1.0f/60,character_ground_query,0);
        for(int f=0;f<2;f++) {
            Vector3 nominal=Vector3Transform((Vector3){(f?1:-1)*.10f*CHARACTER_MODEL_SCALE,0,-.045f*CHARACTER_MODEL_SCALE},MatrixRotateY(yaw*DEG2RAD));
            nominal=Vector3Add(nominal,(Vector3){players[0].pos.x,0,players[0].pos.z});
            assert(isfinite(pose.soles[f].x));assert(Vector3Distance(nominal,pose.soles[f])<.42f);
        }
    }
    players[0].pos.x+=10;
    pose=character_pose(players[0].pos,0,0,0,true,false,SKYBLUE,0);
    character_animate(&animation,&pose,1.0f/60,character_ground_query,0);
    assert(fabsf(pose.soles[0].x-players[0].pos.x)<.3f);
    reset_human((Vector3){0,2,0});players[0].onGround=false;players[0].vel.y=-30;
    net_integrate_player(&players[0],.1f,true);expect_near(players[0].pos.y,1);assert(players[0].onGround);
    players[0].vel.y=JUMP_SPEED;players[0].onGround=false;
    for(int f=0;f<180;f++)net_integrate_player(&players[0],1.0f/60,true);
    expect_near(players[0].pos.y,1);assert(players[0].onGround);
    for(int x=-3;x<=3;x++)for(int z=-3;z<=3;z++)for(int y=0;y<2;y++)cell(x,y,z);
    rebuild_static_hash_if_dirty();reset_human((Vector3){0,3,0});players[0].vel.y=-40;players[0].onGround=false;
    net_integrate_player(&players[0],.1f,true);expect_near(players[0].pos.y,2);assert(players[0].onGround);
    CharacterSupport support=character_ground_query((Vector3){0,1,0},0);assert(support.valid);expect_near(support.height,1);
    clear_world_voxels();rebuild_static_hash_if_dirty();net_integrate_player(&players[0],1.0f/60,true);assert(!players[0].onGround);
    for(int f=0;f<120;f++)net_integrate_player(&players[0],1.0f/60,true);expect_near(players[0].pos.y,1);
    // Swept wall and ceiling: do not tunnel even under knockback.
    for(int z=-2;z<=2;z++)for(int y=0;y<8;y++)cell(2,y,z);
    rebuild_static_hash_if_dirty();reset_human((Vector3){0,1,0});players[0].vel.x=50;players[0].meleeKnockbackActive=true;
    net_integrate_player(&players[0],.1f,true);expect_near(players[0].pos.x,.61f);
    clear_world_voxels();for(int x=-2;x<=2;x++)for(int z=-2;z<=2;z++)cell(x,4,z);
    rebuild_static_hash_if_dirty();reset_human((Vector3){0,1,0});players[0].vel.y=20;players[0].onGround=false;
    net_integrate_player(&players[0],.1f,true);expect_near(players[0].pos.y,1.88f);
    // Enlarged body rejects a half-unit corridor; bots retain their original box.
    clear_world_voxels();for(int y=0;y<4;y++){cell(-1,y,0);cell(1,y,0);}rebuild_static_hash_if_dirty();
    assert(!spawn_position_clear_ex((Vector3){.25f,1,.25f},.25f,0,0));
    playerInput[0]=INPUT_TYPE_BOT_EASY;b=combatant_bounds(0,.25f);expect_near(b.max.y-b.min.y,.5f);
    clear_world_voxels();reset_human((Vector3){0,1,0});
    // Exercise the current backend's body contact uniforms with a resting cube.
    int id=addVoxel(.65f,.35f,0,false,true,BLUE,0);assert(id>=0);
    rebuild_voxel_hash();rebuild_static_hash_if_dirty();
    assert(sim_particle_count>0);
    float beforeX=0;for(int i=0;i<sim_particle_count;i++)beforeX+=sim_particles[i]->pos.x;
    prepare_tether_forces();
    simulate_voxel_pbd_steps(1.0f/120,1);
    for(int i=0;i<voxel_count;i++)assert(isfinite(voxels[i].pos.x)&&isfinite(voxels[i].pos.y));
    float afterX=0;for(int i=0;i<sim_particle_count;i++)afterX+=sim_particles[i]->pos.x;
    assert(afterX>beforeX+.01f); // lower torso contact, below the old eye-centered box
    if(physics_backend_is_gpu(physicsBackend.requested))assert(physicsBackend.active==physicsBackend.requested);
    shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Grounding checks passed");return 0;
}
