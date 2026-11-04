/* Reproduction only: production gameplay functions, GPU timings and a fixed replay.
 * Run from repo root: .build/wave7_tether_profile.exe --physics=gpu
 * FPS_TETHER_CONTROL=1 omits grabs; FPS_TETHER_CYCLES defaults to 15.
 * FPS_TETHER_ROUND=3 selects round three; FPS_TETHER_LARGE=1 re-grabs large debris.
 * FPS_TETHER_REALTIME=1 schedules actions by wall time and uses actual frame dt.
 * Otherwise timers use fixed simulated time; measurements always use wall time.
 * FPS_TETHER_TARGET=shield|tower|scrap targets enemy-owned geometry.
 * FPS_TETHER_QUIET sets idle seconds (default 100).
 * FPS_TETHER_DIAGNOSTICS=1 prints optional counters in diagnostic builds.
 * The human has unlimited Matter/invulnerability to keep the workload repeatable. */
#include "raylib.h"
static double wall_time(void){return GetTime();}
static double replayTime=1000;
static double replay_time(void){return replayTime;}
#define GetTime replay_time
#define FIREFIGHT_TESTING
#define main fps_application_main
#ifdef TETHER_PROFILE_BASELINE
#include "../.build/fps_tether_baseline.c"
#elif defined(TETHER_PROFILE_INSTRUMENTED)
#include "../.build/fps_tether_instrumented.c"
#else
#include "../fps_ray.c"
#endif
#undef main
#include <assert.h>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif
bool ffTestReviveHeld[4];
static double phases[15];
static unsigned long long debrisWork[12];
static int compare_ms(const void *a,const void *b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static void report(int frame,const char *stage,int ticks){
    if(getenv("FPS_TETHER_DIAGNOSTICS")){
        printf("DEBRIS time=%.2f frames=%d searches=%llu cached=%llu rejected=%llu visited=%llu groups=%llu blocked=%llu wake_queries=%llu transitions=%llu protected=%llu mixed=%llu not_calm=%llu unsupported=%llu\n",frame/60.,ticks,debrisWork[0],debrisWork[1],debrisWork[2],debrisWork[3],debrisWork[4],debrisWork[5],debrisWork[6],debrisWork[7],debrisWork[8],debrisWork[9],debrisWork[10],debrisWork[11]);
        memset(debrisWork,0,sizeof(debrisWork));
    }
#ifdef TETHER_PROFILE_INSTRUMENTED
    printf("WORK time=%.2f frames=%d collect=%llu visited=%llu rejects=%llu rejected_cells=%llu wakes=%llu wake_cells=%llu freezes=%llu freeze_cells=%llu\n",frame/60.,ticks,
        atomic_exchange(&profileCollect,0),atomic_exchange(&profileVisited,0),atomic_exchange(&profileReject,0),atomic_exchange(&profileRejectedCells,0),
        profileWakeCalls,profileWakeCells,profileFreezeCalls,profileFreezeCells);
    profileWakeCalls=profileWakeCells=profileFreezeCalls=profileFreezeCells=0;
#endif
    int awake=0,sleep=0,tagged=0,tower=0,owned=0,live=0,held=0;
    for(int v=0;v<voxel_count;v++)if(voxels[v].simulate){awake+=!voxels[v].sleeping;sleep+=voxels[v].sleeping;tagged+=voxels[v].wasTethered;tower+=voxels[v].towerTag!=0;owned+=voxels[v].owner>=0;held+=recycle_voxel_held(v);}
    for(int i=MAX_PLAYERS;i<activePlayers;i++)live+=players[i].respawn_timer<=0;
    printf("%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",frame,stage,firefightWave,firefightRound,live,voxel_count,awake,sleep,sim_particle_count,particle_pool_count,free_particle_count,collision_particle_count,tagged,tower,owned,held);
    for(int k=0;k<15;k++){printf(",%.4f",phases[k]*1000/ticks);phases[k]=0;}puts("");
}
static void diagnose(int schedule){
    build_recycle_physical_island_graph();
    int groups=0,mixed=0,old=0,heldGroups=0,persistent=0,scrap=0;
    for(int r=0;r<voxel_count;r++)if(recycleLifetimeParent[r]==r){
        int n=0,a=0,sl=0,lo=2147483647,hi=0,h=0,pc=0,sc=0,tw=0,due=0;
        for(int v=recycleLifetimeHead[r];v>=0;v=recycleLifetimeNext[v]){
            Voxel *x=&voxels[v];n++;a+=!x->sleeping;sl+=x->sleeping;
            if(x->lifeFrames<lo)lo=x->lifeFrames;if(x->lifeFrames>hi)hi=x->lifeFrames;
            h+=recycle_voxel_held(v);pc+=x->persistentConstruction;sc+=x->owner==ENEMY_SCRAP_OWNER;tw+=x->towerTag!=0;
            due+=x->restorationQueued && x->restorationDeadlineFrame>0 && recycleFrameCounter>=x->restorationDeadlineFrame;
        }
        groups++;mixed+=a>0&&sl>0;old+=hi>600;heldGroups+=h>0;persistent+=pc;scrap+=sc;
        if(n>=8)
            printf("ISLAND time=%.1f root=%d n=%d awake=%d sleep=%d age_min=%d age_max=%d held=%d persistent=%d scrap=%d tower=%d due=%d mixed_blocks_cleanup=%d\n",schedule/60.,r,n,a,sl,lo,hi,h,pc,sc,tw,due,a>0&&sl>0&&!due&&!h);
    }
    printf("DIAG time=%.1f groups=%d mixed=%d old=%d held=%d persistent=%d scrap=%d queue=%d\n",schedule/60.,groups,mixed,old,heldGroups,persistent,scrap,recycleQueueCount);
}
static void grab(int cycle,bool control){
    Player *p=&players[0];if(p->tetherHolding)release_tether(0);
    float angle=cycle*2.399963f;
    p->pos=nudge_spawn_clear((Vector3){27*cosf(angle),BASE_EYE_HEIGHT,27*sinf(angle)},.3f,4,0);p->vel=(Vector3){0};
    if(control)return;
    int selected=-1;float score=1e30f;
    for(int v=0;v<voxel_count;v++){
        Voxel *x=&voxels[v];if(x->isBullet || x->type || x->pos.y<1 || x->pos.y>3.5f || x->persistentConstruction)continue;
        float distance=v_length(world_delta(x->pos,p->pos));if(distance<2 || distance>GOLD_TETHER_RANGE-1)continue;
        float cost=distance+(x->simulate?2:0);if(cost>=score)continue;
        Vector3 delta=world_delta(x->pos,p->pos);VoxelHit hit;
        if(!first_voxel_hit_detailed((Ray){p->pos,v_norm(delta)},distance+.4f,-1,false,&hit))continue;
        selected=hit.id;score=cost;
    }
    if(getenv("FPS_TETHER_LARGE") && cycle%3!=0){
        static unsigned char visited[MAX_VOXELS];memset(visited,0,(size_t)voxel_count);
        int best=0,seed=-1;
        for(int v=0;v<voxel_count;v++){
            if(visited[v] || !voxels[v].simulate || voxels[v].isBullet || voxels[v].type || voxels[v].persistentConstruction || enemy_equipment_owner(v)>=0)continue;
            int n=build_glue_cluster_indices(v,glueClusterIndices);
            for(int k=0;k<n;k++)visited[glueClusterIndices[k]]=1;
            if(n>best && voxels[v].pos.y<8){best=n;seed=v;}
        }
        if(seed>=0 && best>8){
            Vector3 at=voxels[seed].pos;p->pos=nudge_spawn_clear(world_wrap(v_add(at,(Vector3){6,BASE_EYE_HEIGHT,0})),.3f,4,0);
            p->pos.y=fmaxf(BASE_EYE_HEIGHT,p->pos.y);Vector3 delta=world_delta(at,p->pos);VoxelHit hit;
            if(first_voxel_hit_detailed((Ray){p->pos,v_norm(delta)},v_length(delta)+.5f,-1,false,&hit))selected=hit.id;
            printf("large candidate cells=%d\n",best);
        }
    }
    const char *target=getenv("FPS_TETHER_TARGET");
    if(target){
        int seed=-1;
        for(int v=0;v<voxel_count;v++){
            bool wanted=!strcmp(target,"tower")?voxels[v].towerTag!=0:
                !strcmp(target,"shield")?enemy_equipment_owner(v)>=0:
                voxels[v].owner==ENEMY_SCRAP_OWNER;
            if(wanted && !voxels[v].isBullet){seed=v;break;}
        }
        if(seed>=0){
            Vector3 at=voxels[seed].pos;
            p->pos=nudge_spawn_clear(world_wrap(v_add(at,(Vector3){4,BASE_EYE_HEIGHT,0})),.3f,4,0);
            Vector3 delta=world_delta(at,p->pos);VoxelHit hit;
            if(first_voxel_hit_detailed((Ray){p->pos,v_norm(delta)},v_length(delta)+.5f,-1,false,&hit)){
                selected=hit.id;printf("target %s requested=%d actual=%d tower=%u equipment=%d\n",target,seed,selected,voxels[selected].towerTag,enemy_equipment_owner(selected));
            }
        }
    }
    if(selected<0){printf("grab cycle=%d no target\n",cycle);return;}
    Vector3 chosen=voxels[selected].pos;Vector3 delta=world_delta(chosen,p->pos);dir_to_yaw_pitch(delta,&p->yaw,&p->pitch);
    acquire_tether_hit(0,selected,(Vector3){0,0,1});
    int count=p->tetherHolding?build_glue_cluster_indices(p->tetherVoxel,glueClusterIndices):0;
    printf("grab cycle=%d held=%d gold=%d cells=%d target=(%.1f %.1f %.1f)\n",cycle,p->tetherHolding,p->goldTetherHolding,count,chosen.x,chosen.y,chosen.z);
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Wave seven tether reproduction");SetTargetFPS(0);setvbuf(stdout,NULL,_IONBF,0);
    init_pbd_thread_pool();assert(initialize_physics_backend());weaponRenderingReady=weapon_renderer_init();character_renderer_init();drawWeaponModels=true;weaponBloomEnabled=true;
    SetRandomSeed(1234);gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;multiplayerActivePlayers=1;multiplayerPlayerInput[0]=INPUT_TYPE_KEYBOARD;
    ResetGame();gameState=GAME_STATE_PLAYING;start_firefight_wave(7);int round=getenv("FPS_TETHER_ROUND")?atoi(getenv("FPS_TETHER_ROUND")):1;while(firefightRound<round && firefightRound<3)ff_wave_complete();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,8};
    int cycles=getenv("FPS_TETHER_CYCLES")?atoi(getenv("FPS_TETHER_CYCLES")):15;
    bool control=getenv("FPS_TETHER_CONTROL")!=NULL,realtime=getenv("FPS_TETHER_REALTIME")!=NULL;
    if(cycles<1 || cycles>100)return 2;
    int quiet=getenv("FPS_TETHER_QUIET")?atoi(getenv("FPS_TETHER_QUIET")):100;
    int stressEnd=600+cycles*120,quietEnd=stressEnd+quiet*60,frames=quietEnd+1200;
    int capacity=realtime?60000:frames;double *totals=calloc(capacity,sizeof(double)),*sorted=calloc(capacity,sizeof(double));int *sampleStage=malloc(capacity*sizeof(int));assert(totals && sorted && sampleStage);
    RenderTexture2D views[4]={0};int count=0,w=0,h=0;
    printf("GPU %s | %s | control=%d cycles=%d round=%d large=%d\n",glGetString(GL_VENDOR),glGetString(GL_RENDERER),control,cycles,round,getenv("FPS_TETHER_LARGE")!=NULL);
    puts("frame,stage,wave,round,enemies,voxels,awake,sleep,sim_particles,pool,free_particles,collision_particles,tethered_history,tower_dynamic,owned_dynamic,held,events_ms,ai_ms,activation_ms,tether_ms,wake_ms,physics_ms,recycle_ms,hits_ms,sleep_ms,hash_ms,render_ms,total_ms,projectiles_ms,pickups_ms,drain_ms");
    for(int warm=0;warm<5;warm++){voxel_render_frame++;render_gameplay_view(views,&count,&w,&h,false);glFinish();}
    double replayStart=wall_time(),lastFrame=replayStart;int samples=0,sampleTicks=0,lastReport=0,lastCycle=-1;bool released=false,reset=false,cycleReleased=false;int lastDiag=-1;
    for(int f=0;f<capacity;f++){
        double start=wall_time(),t=start;
        int schedule=realtime?(int)((start-replayStart)*60):f;if(schedule>=frames)break;
        float dt=realtime&&f?(float)(start-lastFrame):1.f/60;lastFrame=start;replayTime+=dt;
        sampleStage[f]=schedule<600?0:schedule<stressEnd?1:schedule<quietEnd?2:3;
        if(schedule<120 || (sampleStage[f]==2 && schedule<stressEnd+600) || (sampleStage[f]==3 && schedule<quietEnd+120))sampleStage[f]=-1;
        players[0].invuln_timer=1000;players[0].matter=players[0].matterMax=100;players[0].isExposed=false;players[0].goldTetherCharged=true;
        if(schedule>=stressEnd && !released){if(players[0].tetherHolding)release_tether(0);released=true;}
        if(schedule>=quietEnd && !reset){reset=true;puts("ACTION console wave 8");start_firefight_wave(8);}
        if(schedule>=600 && schedule<stressEnd){int tick=(schedule-600)%120,cycle=(schedule-600)/120;
            if(cycle!=lastCycle){grab(cycle,control);lastCycle=cycle;cycleReleased=false;}
            if(tick>0 && tick<75 && players[0].tetherHolding){players[0].yaw+=60*dt;players[0].pitch=10;}
            if(tick>=75 && !cycleReleased){if(players[0].tetherHolding)release_tether(0);cycleReleased=true;}
        }
        clarity_tick(dt);eq_tick(dt);hologram_tick(dt);ducks_tick(dt);update_firefight_logic(dt);phases[0]+=wall_time()-t;t=wall_time();
        update_melee_swings();
        for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0){if(is_player_bot(i))UpdateBot(i,dt);net_integrate_player(&players[i],dt,true);}
        for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0){if(players[i].enemyType==ENEMY_TYPE_GOLIATH)update_goliath_armor_positions(i);if(enemy_is_special(i))enemy_update_equipment(i);}
        update_enemy_contact_damage(dt);phases[1]+=wall_time()-t;t=wall_time();
        activate_static_voxels_near_dynamic();phases[2]+=wall_time()-t;t=wall_time();
        update_projectiles(dt);phases[12]+=wall_time()-t;t=wall_time();
        update_pickups(dt);phases[13]+=wall_time()-t;t=wall_time();
        update_tether_drains(dt);phases[14]+=wall_time()-t;t=wall_time();
        prepare_tether_forces();phases[3]+=wall_time()-t;t=wall_time();
        wake_sleeping_voxels_near_awake();rebuild_static_hash_if_dirty();phases[4]+=wall_time()-t;t=wall_time();
        pbdTimeAccumulator=fminf(pbdTimeAccumulator+dt,PBD_MAX_STEP_DT*PBD_MAX_ACCUM_STEPS);
        int due=(int)(pbdTimeAccumulator/PBD_MAX_STEP_DT);if(due){simulate_voxel_pbd_steps(PBD_MAX_STEP_DT,due);pbdTimeAccumulator-=due*PBD_MAX_STEP_DT;}phases[5]+=wall_time()-t;t=wall_time();
        recycle_dead_voxels(dt);phases[6]+=wall_time()-t;t=wall_time();
        handle_pbd_projectile_hits();phases[7]+=wall_time()-t;t=wall_time();
        deactivate_sleeping_voxels();phases[8]+=wall_time()-t;t=wall_time();
#if !defined(TETHER_PROFILE_BASELINE) && !defined(TETHER_PROFILE_INSTRUMENTED)
        if(getenv("FPS_TETHER_DIAGNOSTICS")){
            for(int worker=0;worker<=PBD_MAX_THREADS;worker++){
                ActivationCandidates *c=&activationCandidates[worker];
                debrisWork[0]+=c->searches;debrisWork[1]+=c->hits;debrisWork[2]+=c->rejects;debrisWork[3]+=c->visited;
            }
            debrisWork[4]+=sleepGroupsChecked;debrisWork[5]+=sleepGroupsBlocked;debrisWork[6]+=wakeQueries;debrisWork[7]+=wakeTransitions;
            debrisWork[8]+=sleepProtected;debrisWork[9]+=sleepMixed;debrisWork[10]+=sleepNotCalm;debrisWork[11]+=sleepUnsupported;
        }
#endif
        if(++voxelHashFramesSinceRebuild>=VOXEL_HASH_REBUILD_INTERVAL)rebuild_voxel_hash();rebuild_static_hash_if_dirty();phases[9]+=wall_time()-t;t=wall_time();
        voxel_render_frame++;render_gameplay_view(views,&count,&w,&h,false);glFinish();phases[10]+=wall_time()-t;
        totals[f]=(wall_time()-start)*1000;phases[11]+=totals[f]/1000;
        samples=f+1;sampleTicks++;
        if(schedule/600!=lastDiag){diagnose(schedule);lastDiag=schedule/600;}
        if(realtime){double remaining=1.0/60-(wall_time()-start);if(remaining>0)WaitTime(remaining);}
        if(schedule/120>lastReport){report(schedule,schedule<600?"baseline":schedule<stressEnd?"tether":schedule<quietEnd?"released":"wave-reset",sampleTicks);sampleTicks=0;lastReport=schedule/120;}
    }
    const char *names[]={"baseline","tether","quiet-tail","after-wave-reset"};
    for(int k=0;k<4;k++){int n=0;for(int f=0;f<samples;f++)if(sampleStage[f]==k)sorted[n++]=totals[f];if(!n)continue;qsort(sorted,n,sizeof(double),compare_ms);printf("SUMMARY %s median_ms=%.3f p95_ms=%.3f fps=%.1f samples=%d realtime=%d\n",names[k],sorted[n/2],sorted[(n-1)*95/100],1000/sorted[n/2],n,realtime);}
    free(totals);free(sorted);free(sampleStage);for(int i=0;i<count;i++)UnloadRenderTexture(views[i]);clear_world_voxels();shutdown_world_visuals();weapon_renderer_shutdown();character_renderer_shutdown();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
