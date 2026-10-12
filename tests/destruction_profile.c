/* Destruction stress profiler. Run from repo root with --physics=gpu.
 * FPS_STRESS_GRENADES (0..24), FPS_STRESS_FRAMES (>=120) configure the workload.
 * Gameplay AI is disabled to isolate terrain activation, physics and rendering. */
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static int profile_compare(const void *a,const void *b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static void profile_summary(const char *label,double *values,int count) {
    if(count<=0)return;
    double sum=0;for(int i=0;i<count;i++)sum+=values[i];qsort(values,count,sizeof(*values),profile_compare);
    fprintf(stderr,"timing %s count=%d mean=%.3f p50=%.3f p95=%.3f p99=%.3f max=%.3f ms\n",label,count,sum/count,values[count/2],values[(int)((count-1)*.95)],values[(int)((count-1)*.99)],values[count-1]);
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Grenade stress profile");SetTargetFPS(0);
    init_pbd_thread_pool();assert(initialize_physics_backend());weaponRenderingReady=weapon_renderer_init();character_renderer_init();
    bool versus=getenv("FPS_STRESS_MAP") && !strcmp(getenv("FPS_STRESS_MAP"),"versus");
    gameMode=versus?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=!versus;
    if(versus)currentWorldType=WORLD_TYPE_CROSSCUT;
multiplayerActivePlayers=1;multiplayerPlayerInput[0]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;
    for(int i=0;i<activePlayers;i++)players[i].respawn_timer=9999;
    int limit=getenv("FPS_STRESS_GRENADES")?atoi(getenv("FPS_STRESS_GRENADES")):24;
    int frames=getenv("FPS_STRESS_FRAMES")?atoi(getenv("FPS_STRESS_FRAMES")):3600;
    if(limit<0 || limit>24 || frames<120)return 2;
    double *frameTimes=calloc(frames,sizeof(double));assert(frameTimes);
    double commitTimes[2048];int commitCount=0;
    Vector3 targets[24];int count=0;
    for(int i=0;i<voxel_count && count<limit;i+=37)if(!voxels[i].simulate && voxels[i].pos.y>.4f && voxels[i].pos.y<5){targets[count++]=v_add(voxels[i].pos,(Vector3){0,2,0});}
    Camera3D cam={{24,18,24},{0,2,0},{0,1,0},65,CAMERA_PERSPECTIVE};
    double phases[9]={0};int ticks=0;
    fprintf(stderr,"frame,voxels,awake,sleep,sim_particles,pool,collision,chunks,queue,invisible,fixed_particles,blast_ms,activate_ms,physics_ms,recycle_ms,sleep_ms,hash_ms,mesh_ms,render_ms,total_ms\n");
    for(int f=0;f<frames;f++) {
        double t=GetTime(),start=t;uint64_t commitsBefore=regenerationCommits;
        if(getenv("FPS_STRESS_ACTIVATE") && f==180) {
            double begin=GetTime();bool activated=activate_static_voxels_near_region(39,46,0,9,-10,10,0);
            double elapsed=(GetTime()-begin)*1000;int dynamic=0;
            for(int i=0;i<voxel_count;i++)dynamic+=voxels[i].simulate;
            fprintf(stderr,"activation changed=%d dynamic_voxels=%d synchronous_ms=%.3f\n",activated,dynamic,elapsed);
        }
        if(f>=180 && (f-180)%30==0 && (f-180)/30<count){int n=(f-180)/30;grenades[n]=(GrenadeState){.id=100+n,.owner=-1,.pos=targets[n],.vel={0,-2,0},.fuse=-1};}
        clarity_tick(1.0f/60);phases[0]+=GetTime()-t;t=GetTime();
        activate_static_voxels_near_dynamic();prepare_tether_forces();wake_sleeping_voxels_near_awake();rebuild_static_hash_if_dirty();phases[1]+=GetTime()-t;t=GetTime();
        simulate_voxel_pbd_steps(1.0f/120,2);phases[2]+=GetTime()-t;t=GetTime();
        recycle_dead_voxels(1.0f/60);phases[3]+=GetTime()-t;t=GetTime();
        deactivate_sleeping_voxels();phases[4]+=GetTime()-t;t=GetTime();
        rebuild_voxel_hash();rebuild_static_hash_if_dirty();phases[5]+=GetTime()-t;t=GetTime();
        if(meshDirty || (voxel_chunk_count==0 && voxel_count>0)){mark_all_chunks_dirty();meshDirty=false;}
        for(int c=0;c<voxel_chunk_count;c++)if(voxel_chunks[c].in_use && voxel_chunks[c].dirty)rebuild_chunk(&voxel_chunks[c]);
        phases[6]+=GetTime()-t;t=GetTime();
        voxel_render_frame++;
        BeginDrawing();ClearBackground(SKYBLUE);BeginMode3D(cam);draw_world_surfaces();DrawVoxels(cam);draw_pickups(cam);EndMode3D();EndDrawing();phases[7]+=GetTime()-t;
        frameTimes[f]=(GetTime()-start)*1000;
        if(regenerationCommits!=commitsBefore) {
            if(commitCount<2048)commitTimes[commitCount++]=regenerationLastCommitMs;
            fprintf(stderr,"commit frame=%d blocks_remaining=%d transaction_ms=%.3f frame_ms=%.3f\n",f+1,recycleQueueCount,regenerationLastCommitMs,frameTimes[f]);
        }
        phases[8]+=GetTime()-start;ticks++;
        if(f%120==119){int awake=0,sleep=0,invisible=0,fixed=0;for(int v=0;v<voxel_count;v++)invisible+=voxels[v].simulate && !voxels[v].vgs_active;for(int p=0;p<sim_particle_count;p++)fixed+=sim_particles[p]->inv_mass==0;for(int i=0;i<voxel_count;i++)if(voxels[i].simulate){if(voxels[i].sleeping)sleep++;else awake++;}
            fprintf(stderr,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",f+1,voxel_count,awake,sleep,sim_particle_count,particle_pool_count,collision_particle_count,voxel_chunk_count,recycleQueueCount,invisible,fixed);
            for(int k=0;k<9;k++){fprintf(stderr,",%.3f",phases[k]*1000/ticks);phases[k]=0;}fprintf(stderr,"\n");fflush(stderr);ticks=0;
        }
    }
    if(frames>2400)profile_summary("quiet-tail",frameTimes+2400,frames-2400);
    profile_summary("all-frames",frameTimes,frames);profile_summary("transactions",commitTimes,commitCount);free(frameTimes);
    if(getenv("FPS_STRESS_RESTORE_ONLY")) {
        /* Controlled isolation: queue 300 real terrain voxels with no explosion/AI. */
        int queued=0;
        for(int n=voxel_count-1;n>=0 && queued<300;n--)if(!voxels[n].simulate && voxels[n].owner==-1){
            assert(recycle_queue_push(&voxels[n]));remove_voxel_index(n);queued++;
        }
        recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;
        rebuild_voxel_hash();rebuild_static_hash_if_dirty();meshDirty=true;
        BeginDrawing();BeginMode3D(cam);DrawVoxels(cam);EndMode3D();EndDrawing();
    }
    /* Hold the final world still and isolate the restoration queue's render cost. */
    for(int mode=0;mode<2;mode++) {
        double start=GetTime();int queued=recycleQueueCount;
        for(int f=0;f<120;f++) {
            if(!mode)recycle_dead_voxels(1.0f/60);
            voxel_render_frame++;
            BeginDrawing();BeginMode3D(cam);DrawVoxels(cam);EndMode3D();EndDrawing();
        }
        fprintf(stderr,"render-control restore=%d queue=%d->%d mean_ms=%.3f\n",!mode,queued,recycleQueueCount,(GetTime()-start)*1000/120);
    }
    regeneration_shutdown();weapon_renderer_shutdown();character_renderer_shutdown();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();
    return 0;
}
