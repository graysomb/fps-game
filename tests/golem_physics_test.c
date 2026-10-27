#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static int attached_count(void){int n=0;for(int c=0;c<goCount;c++)n+=goCells[c].index>=0;return n;}
static void healthy(void){
    assert(attached_count()==goCount);
    for(int c=0;c<goCount;c++)assert(voxels[goCells[c].index].vgs_active);
}
static int compare_time(const void *a,const void *b){float x=*(const float*)a,y=*(const float*)b;return (x>y)-(x<y);}
static void lab_frame(float dt){
    activate_static_voxels_near_dynamic();update_projectiles(dt);prepare_tether_forces();
    wake_sleeping_voxels_near_awake();rebuild_static_hash_if_dirty();
    pbdTimeAccumulator+=dt;
    int steps=(int)(pbdTimeAccumulator/PBD_MAX_STEP_DT);
    if(steps>PBD_MAX_ACCUM_STEPS)steps=PBD_MAX_ACCUM_STEPS;
    if(steps){simulate_voxel_pbd_steps(PBD_MAX_STEP_DT,steps);pbdTimeAccumulator-=steps*PBD_MAX_STEP_DT;}
    recycle_dead_voxels(dt);handle_pbd_projectile_hits();
    if(++voxelHashFramesSinceRebuild>=VOXEL_HASH_REBUILD_INTERVAL)rebuild_voxel_hash();
    deactivate_sleeping_voxels();rebuild_static_hash_if_dirty();
}
static void lab_check(void){
    const int rates[]={30,60,144};
    for(int r=0;r<3;r++){
        enemy_lab_start(ENEMY_TYPE_COUNT+1);go_combat_reset();pbdTimeAccumulator=0;
        Vector3 inactive=players[MAX_PLAYERS].pos;
        int probe=addVoxel(inactive.x,inactive.y,inactive.z,false,true,GRAY,0);assert(probe>=0);
        uint64_t probeId=voxels[probe].identity;handle_pbd_projectile_hits();
        probe=go_find(probeId,probe);assert(probe>=0);remove_voxel_index(probe);
        for(int f=0;f<rates[r]*60;f++){
            lab_frame(1.f/rates[r]);
            int live=0;for(int c=0;c<goCount;c++)live+=go_find(goCells[c].id,goCells[c].index)>=0;
            if(live!=goCount){printf("LAB LOSS fps=%d frame=%d live=%d/%d\n",rates[r],f,live,goCount);fflush(stdout);}
            assert(live==goCount);healthy();
        }
        for(int f=0;f<20;f++){lab_frame(.10f);healthy();}
        printf("Lab startup: 60 seconds at %d FPS passed\n",rates[r]);
    }
}
static void fixture(void){
    clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();
    memset(players,0,sizeof(players));netTransport.role=NET_ROLE_OFFLINE;activePlayers=1;worldPeriodic=false;worldFloorHalfSize=32;
    players[0].pos=(Vector3){-6,BASE_EYE_HEIGHT,0};players[0].matter=100;players[0].respawn_timer=9999;
    assert(golem_physics_start());goRegen=false;
    assert(goCount==291);
    for(int i=0;i<goCount;i++){
        Voxel *v=&voxels[goCells[i].index];assert(v->simulate&&!v->fixed&&v->golemOwned);
        assert(v->rest_volume>0 && v->rest_dimensions.x>0 && v->rest_dimensions.y>0 && v->rest_dimensions.z>0);
        for(int k=0;k<8;k++)assert(Vector3Distance(v->particles[k]->pos,Vector3Transform(go_corner(&goCells[i],k),go_matrix(goCells[i].part)))<.0001f);
    }
    printf("Body: %d cells, %d connections\n",goCount,goWeldCount);
}
static void capture(const char *name){
    voxel_render_frame++;
    RenderTexture2D rt=LoadRenderTexture(800,600);Camera3D cam={{-9,6,0},{0,2.8f,10},{0,1,0},43,CAMERA_PERSPECTIVE};
    BeginTextureMode(rt);ClearBackground((Color){19,28,39,255});BeginMode3D(cam);DrawPlane((Vector3){0,-.01f,10},(Vector2){50,50},(Color){40,51,62,255});DrawVoxels(cam);EndMode3D();
    DrawText(TextFormat("PHYSICAL GOLEM | %.2fs | error %.3fm",goClock,goError),15,15,20,SKYBLUE);EndTextureMode();
    Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);ExportImage(image,TextFormat("artifacts/golem/physics/%s.png",name));UnloadImage(image);UnloadRenderTexture(rt);
}
static void accuracy_check(void){
    const float ranges[]={6,12,20,30},angles[]={-45,0,45};
    puts("ACCURACY backend,range,angle,closest_center_m,closest_box_m,ground_shortfall_m");
    for(int r=0;r<4;r++)for(int a=0;a<3;a++){
        fixture();goAction=GOLEM_THROW;playerInput[0]=INPUT_TYPE_KEYBOARD;
        combatantLife[0]=LIFE_ALIVE;players[0].respawn_timer=0;
        float angle=angles[a]*DEG2RAD;
        Vector3 target={sinf(angle)*ranges[r],1.2f,10-cosf(angle)*ranges[r]};
        players[0].pos=v_add(target,(Vector3){0,.45f,0});
        uint64_t ids[64];int count=0;
        for(int c=0;c<goCount;c++)if(goCells[c].part==GOLEM_HAND_R)ids[count++]=goCells[c].id;
        float closest=FLT_MAX,boxMiss=FLT_MAX,shortfall=0;int launched=-1;
        for(int f=0;f<360;f++){
            simulate_voxel_pbd_steps(1.f/60,1);
            if(goCells[goMount[GOLEM_HAND_R]].index>=0)continue;
            if(launched<0){launched=f;players[0].respawn_timer=9999;}
            Vector3 center={0},lo={FLT_MAX,FLT_MAX,FLT_MAX},hi={-FLT_MAX,-FLT_MAX,-FLT_MAX};
            for(int c=0;c<count;c++){
                int index=go_find(ids[c],-1);assert(index>=0);
                center=v_add(center,voxels[index].pos);
                for(int k=0;k<8;k++){
                    Vector3 p=voxels[index].particles[k]->pos;assert(v_isfinite(p));
                    lo=Vector3Min(lo,p);hi=Vector3Max(hi,p);
                }
            }
            center=v_mul(center,1.f/count);
            closest=fminf(closest,Vector3Distance(center,target));
            boxMiss=fminf(boxMiss,Vector3Distance(Vector3Clamp(target,lo,hi),target));
            shortfall=(target.x-center.x)*sinf(angle)-(target.z-center.z)*cosf(angle);
            if(f>launched+2&&lo.y<.06f)break;
        }
        assert(launched>=0);
        printf("ACCURACY %s,%.0f,%.0f,%.3f,%.3f,%.3f\n",physics_backend_name(physicsBackend.active),ranges[r],angles[a],closest,boxMiss,shortfall);
    }
}
static void repair_check(void){
    const int losses[]={GOLEM_SHIN_L,GOLEM_THIGH_L,GOLEM_ARM_R,GOLEM_HAND_R,GM_CHEST,GOLEM_HEAD};
    for(int t=0;t<6;t++){
        fixture();goRegen=true;go_release_part(losses[t]);
        if(t==1)go_release_part(GOLEM_THIGH_R);
        int before=attached_count();bool restored=false;
        for(int f=0;f<480;f++){
            simulate_voxel_pbd_steps(1.f/60,1);recycle_dead_voxels(1.f/60);
            int now=attached_count();assert(now==before||now==goCount);
            if(now==goCount)restored=true;
            if(restored)healthy();
            if(getenv("GOLEM_CAPTURE")&&t==1&&f<360&&f%3==0)capture(TextFormat("leg-regrowth/%03d",f/3));
        }
        healthy();assert(goPose.kneel<.01f);
        for(int part=1;part<GM_COUNT;part++){
            int cluster=build_glue_cluster_indices(goCells[goMount[part]].index,glueClusterIndices);
            for(int c=0;c<goCount;c++)if(goCells[c].part==part){
                bool found=false;for(int n=0;n<cluster;n++)found|=glueClusterIndices[n]==goCells[c].index;
                assert(found);
            }
        }
        printf("Atomic repair with debris: part %d passed\n",losses[t]);
    }
    /* An allocation failure after some cells were created must not publish
     * any part of the replacement, consume debris, or advance recovery. */
    fixture();goRegen=true;go_release_part(GOLEM_ARM_R);goClock=2;
    int before=attached_count(),total=voxel_count;
    goRepairFailAfter=3;
    for(int f=0;f<30;f++)simulate_voxel_pbd_steps(1.f/60,1);
    assert(attached_count()==before&&voxel_count==total);
    goRepairFailAfter=-1;
    for(int f=0;f<300;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();
    /* Remove just one cell and verify the complete authored graph is restored. */
    fixture();goRegen=true;remove_voxel_index(goCells[goMount[GOLEM_HAND_L]].index);
    for(int f=0;f<360;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();
    /* A held piece cannot be reclaimed even when it occupies the exact target. */
    fixture();goRegen=true;
    int cell=goMount[GOLEM_HAND_R],held=goCells[cell].index;uint64_t heldId=voxels[held].identity;
    go_release_part(GOLEM_HAND_R);goClock=2;
    players[0].tetherHolding=true;players[0].tetherVoxel=held;players[0].tetherVoxelIdentity=heldId;
    for(int i=0;i<20;i++)golem_physics_prepare();
    assert(goCells[cell].index<0&&go_find(heldId,held)>=0&&voxels[held].golemReclaimStart==0);
    players[0].tetherHolding=false;
    for(int f=0;f<300;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();
    /* An actor blocks the repair; clearing it doesn't restart the delay. */
    fixture();goRegen=true;cell=goMount[GOLEM_HAND_R];
    Vector3 target=Vector3Transform(goCells[cell].center,go_matrix(GOLEM_HAND_R));
    go_release_part(GOLEM_HAND_R);goClock=2;
    players[0].pos=target;players[0].respawn_timer=0;
    golem_physics_prepare();assert(goCells[cell].index<0);
    assert(goRepair[GOLEM_ARM_R].phase==GO_REPAIR_BLOCKED);
    players[0].respawn_timer=9999;
    for(int f=0;f<300;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();puts("Allocation rollback, isolated cell, held debris and actor obstruction passed");
}
static void repair_edges(void){
    fixture();goRegen=true;go_release_part(GOLEM_SHIN_L);
    for(int f=0;f<300&&attached_count()<goCount;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();assert(goRepair[GOLEM_THIGH_L].phase==GO_REPAIR_RECOVER);
    remove_voxel_index(goCells[goMount[GOLEM_SHIN_R]].index);
    for(int f=0;f<360;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();assert(goPose.kneel<.01f);
    fixture();goRegen=true;go_release_part(GOLEM_HAND_R);goClock=2;
    golem_physics_prepare();assert(goRepair[GOLEM_ARM_R].phase==GO_REPAIR_BUILD);
    goRegen=false;golem_physics_prepare();
    for(int n=0;n<voxel_count;n++)assert(voxels[n].golemReclaimStart==0);
    goRegen=true;
    for(int f=0;f<300;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();puts("Damage during recovery and canceled dissolve passed");
}
int main(int argc,char **argv){
    printf("Velocity damping: %.6f per 1/120 second\n",(double)VELOCITY_DAMPING);
    setvbuf(stdout,NULL,_IONBF,0);assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,600,"Physical golem checks");init_pbd_thread_pool();if(!initialize_physics_backend()){fprintf(stderr,"Backend init: %s\n",physicsBackend.fallback_reason);return 2;}
    if(getenv("GOLEM_REPAIR_EDGES")){repair_edges();goto cleanup;}
    if(getenv("GOLEM_LAB_CHECK")){lab_check();goto cleanup;}
    if(getenv("GOLEM_REPAIR_CHECK")){repair_check();repair_edges();goto cleanup;}
    if(getenv("GOLEM_ACCURACY")){accuracy_check();goto cleanup;}
    fixture();capture("initial");float timings[360];
    for(int frame=0;frame<360;frame++){
        simulate_voxel_pbd_steps(1.f/60,1);
        timings[frame]=(float)physicsBackend.last_step_ms;
        if(getenv("GOLEM_CAPTURE")&&frame%3==0)capture(TextFormat("idle/%03d",frame/3));
        if(frame<6||frame%60==59){int dead=0;for(int q=0;q<voxel_count;q++)dead+=!voxels[q].vgs_active;printf("inactive %d ",dead);int live=0;for(int i=0;i<goCount;i++)live+=goCells[i].index>=0;printf("Frame %d attached %d constraints %d error %.4f backend %s\n",frame,live,skeletonConstraintCount,goError,physics_backend_name(physicsBackend.active));}
        if(frame==2)for(int c=0;c<goCount;c++)if(goCells[c].index>=0&&!voxels[goCells[c].index].vgs_active)printf("broken cell %d part %d\n",c,goCells[c].part);
        for(int i=0;i<sim_particle_count;i++)assert(v_isfinite(sim_particles[i]->pos));
    }
    capture("settled");
    qsort(timings,360,sizeof(float),compare_time);
    printf("Physics step median %.3fms p95 %.3fms (%s)\n",timings[180],timings[342],physics_backend_name(physicsBackend.active));
    assert(fabsf(goClock-6.f)<.005f);
    for(int c=0;c<goCount;c++)assert(goCells[c].index>=0 && voxels[goCells[c].index].vgs_active);
    if(getenv("GOLEM_IDLE_ONLY"))goto cleanup;
    for(int action=1;action<=4;action++){
        fixture();goAction=action;
        uint64_t fist=goCells[goMount[GOLEM_HAND_R]].id;
        Vector3 launchOrigin={0};int launchFrame=-1;
        for(int f=0;f<360;f++){
            simulate_voxel_pbd_steps(1.f/60,1);
            if(action==4){
                int index=go_find(fist,-1);assert(index>=0);
                if(launchFrame<0&&goCells[goMount[GOLEM_HAND_R]].index<0){
                    launchFrame=f;launchOrigin=voxels[index].pos;
                    assert(voxels[index].vel.z < -15);
                    assert(goCells[goMount[GOLEM_FOREARM_R]].index>=0);
                }
                if(launchFrame>=0&&f==launchFrame+20){
                    float advance=launchOrigin.z-voxels[index].pos.z;
                    printf("Launched fist travelled %.2fm forward in 0.33s\n",advance);
                    assert(advance>4);
                }
            }
            if(getenv("GOLEM_CAPTURE")&&f%3==0)capture(TextFormat("motion-%d/%03d",action,f/3));
        }
        int attached=0;for(int c=0;c<goCount;c++)attached+=goCells[c].index>=0;
        printf("Action %d: attached %d / %d error %.3f\n",action,attached,goCount,goError);
        if(action<4)healthy();else {assert(attached==goCount-12);assert(launchFrame>=0);}
        for(int p=1;p<GM_COUNT;p++){int missing=0,total=0;for(int c=0;c<goCount;c++)if(goCells[c].part==p){total++;missing+=goCells[c].index<0;}if(missing)printf("Part %d missing %d/%d\n",p,missing,total);}
        capture(TextFormat("action-%d",action));
    }
    fixture();
    go_release_part(GOLEM_ARM_R);for(int i=0;i<60;i++)simulate_voxel_pbd_steps(1.f/60,1);capture("arm-off");
    assert(goCells[goMount[GOLEM_ARM_R]].index<0);
    assert(goCells[goMount[GOLEM_FOREARM_R]].index<0);
    assert(goCells[goMount[GOLEM_HAND_R]].index<0);
    /* Removal compacts the voxel pool; identities must still resolve correctly. */
    uint64_t head=goCells[goMount[GOLEM_HEAD]].id;
    remove_voxel_index(0);golem_physics_prepare();
    assert(go_find(head,goCells[goMount[GOLEM_HEAD]].index)>=0);
    fixture();go_release_part(GOLEM_SHIN_L);
    for(int f=0;f<360;f++){
        simulate_voxel_pbd_steps(1.f/60,1);
        if(getenv("GOLEM_CAPTURE")&&f%3==0)capture(TextFormat("leg-loss/%03d",f/3));
    }
    assert(goInput.damage[GOLEM_LEFT_LEG].removed);capture("kneeling");
    /* Regrow while the fallen arm remains in the world. */
    fixture();go_release_part(GOLEM_ARM_R);
    goRegen=true;
    for(int f=0;f<360;f++){
        simulate_voxel_pbd_steps(1.f/60,1);
        if(getenv("GOLEM_CAPTURE")&&f%3==0)capture(TextFormat("regeneration/%03d",f/3));
    }
    printf("Regeneration: attached %d / %d\n",attached_count(),goCount);capture("regenerated");
    healthy();
    if(getenv("GOLEM_STRESS_SECONDS")){
        int seconds=atoi(getenv("GOLEM_STRESS_SECONDS"));assert(seconds>0&&seconds<=3600);
        fixture();goRegen=true;int peak=voxel_count;
        for(int f=0;f<seconds*60;f++){
            if(f%1200==600)go_release_part((f/1200)%2?GOLEM_SHIN_L:GOLEM_ARM_R);
            simulate_voxel_pbd_steps(1.f/60,1);recycle_dead_voxels(1.f/60);
            if(voxel_count>peak)peak=voxel_count;
            assert(recycleQueueCount==0);
            for(int i=0;i<sim_particle_count;i++)assert(v_isfinite(sim_particles[i]->pos));
            if(f%3600==3599)printf("Stress %ds: %d attached, %d cells, %d peak\n",(f+1)/60,attached_count(),voxel_count,peak);
        }
        printf("Stress completed: %ds; %d peak cells\n",seconds,peak);
    }
    /* A contacted face cut is the engine's normal structural-damage path. */
    fixture();int cut=goMount[GOLEM_ARM_R];
    for(int f=0;f<6;f++)break_face_link(&voxels[goCells[cut].index],f);
    golem_physics_prepare();assert(attached_count()<goCount);
    /* Real tether acquisition and solver forces must overcome a mount. */
    fixture();players[0].respawn_timer=0;players[0].goldTetherCharged=true;
    playerInput[0]=INPUT_TYPE_KEYBOARD;combatantLife[0]=LIFE_ALIVE;
    acquire_tether_hit(0,goCells[goMount[GOLEM_ARM_R]].index,(Vector3){1,0,0});
    assert(players[0].tetherHolding&&players[0].goldTetherHolding);
    for(int f=0;f<120;f++){prepare_tether_forces();simulate_voxel_pbd_steps(1.f/60,1);}
    assert(goCells[goMount[GOLEM_ARM_R]].index<0);release_tether(0);
    fixture();go_release_part(GOLEM_ARM_R);
    for(int i=voxel_count-1;i>=0;i--)if(voxels[i].golemOwned&&!voxels[i].persistentConstruction)remove_voxel_index(i);
    golem_physics_prepare();
    Vector3 blocked=Vector3Transform(goCells[goMount[GOLEM_ARM_R]].center,go_matrix(GOLEM_ARM_R));
    int blocker=addVoxel(blocked.x,blocked.y,blocked.z,true,false,GRAY,0);assert(blocker>=0);
    uint64_t blockerId=voxels[blocker].identity;goRegen=true;
    for(int f=0;f<240;f++)simulate_voxel_pbd_steps(1.f/60,1);
    assert(goCells[goMount[GOLEM_ARM_R]].index<0);
    remove_voxel_index(go_find(blockerId,blocker));
    for(int f=0;f<240;f++)simulate_voxel_pbd_steps(1.f/60,1);
    healthy();
cleanup:
    clear_world_voxels();assert(!goActive&&skeletonConstraintCount==0&&voxel_count==0&&skeletonConstraints==NULL&&goRollback.cells==NULL);
    shutdown_pbd_thread_pool();gpu_physics_shutdown();
    shutdown_shard_rendering();clear_all_chunks();
    if(voxelInstanceVboId)rlUnloadVertexBuffer(voxelInstanceVboId);
    if(greedyMaterialInit){greedyMaterial.maps[MATERIAL_MAP_DIFFUSE].texture=(Texture2D){0};UnloadMaterial(greedyMaterial);}
    shutdown_world_visuals();CloseWindow();puts("Physical golem test completed");return 0;
}
