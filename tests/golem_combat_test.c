#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void frame(void){
    activate_static_voxels_near_dynamic();update_projectiles(1.f/60);eq_tick(1.f/60);
    prepare_tether_forces();wake_sleeping_voxels_near_awake();rebuild_static_hash_if_dirty();
    simulate_voxel_pbd_steps(1.f/60,1);recycle_dead_voxels(1.f/60);handle_pbd_projectile_hits();
    rebuild_voxel_hash();deactivate_sleeping_voxels();
}
static void fixture(void){enemy_lab_start(ENEMY_TYPE_COUNT+1);enemyLabGod=true;assert(goCombatEnabled&&!goCombatDead);}
static int attached(void){int n=0;for(int c=0;c<goCount;c++)n+=go_find(goCells[c].id,goCells[c].index)>=0;return n;}
static void capture(const char *name){
    char file[256];snprintf(file,sizeof(file),"artifacts/golem/combat/%s.png",name);
    voxel_render_frame++;RenderTexture2D rt=LoadRenderTexture(960,720);
    Camera3D camera={{-12,8,-6},{0,3.3f,10},{0,1,0},45,CAMERA_PERSPECTIVE};
    BeginTextureMode(rt);ClearBackground((Color){19,28,39,255});BeginMode3D(camera);
    DrawPlane((Vector3){0,-.01f,10},(Vector2){64,64},(Color){40,51,62,255});DrawVoxels(camera);eq_world();clarity_draw();EndMode3D();enemy_lab_hud(960,720);EndTextureMode();
    Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);ExportImage(image,file);UnloadImage(image);UnloadRenderTexture(rt);
}

static void arsenal_tests(void){
    fixture();
    int cell=goCells[goMount[GOLEM_HAND_R]].index,other=goCells[goMount[GOLEM_HAND_L]].index;
    assert(go_shield_hit(cell,35));assert(goShields[GOLEM_HAND_R].value==15);
    assert(goShields[GOLEM_HAND_L].value==50&&goFight.shield==150);
    assert(go_shield_hit(cell,35)&&goShields[GOLEM_HAND_R].value==0);
    assert(go_find(goCells[goMount[GOLEM_HAND_R]].id,cell)>=0);
    recycleFrameCounter++;assert(!go_shield_hit(cell,1));
    goFight.time=7.9f;go_shields_tick(.1f);assert(goShields[GOLEM_HAND_R].value==0);
    goFight.time=8;go_shields_tick(1);assert(goShields[GOLEM_HAND_R].value==5);
    goShields[GOLEM_HAND_R].suppressed=5;go_shields_tick(1);assert(goShields[GOLEM_HAND_R].value==5);
    assert(go_drain_start(0,other));assert(!players[0].tetherHolding);release_tether(0);assert(goDrainPart[0]<0);
    go_release_part(GOLEM_ARM_R);go_shields_tick(.1f);assert(!go_part_alive(GM_TURRET_R));assert(goShields[GOLEM_HAND_R].value==0);
    fixture();goFight.target=0;goFight.aim=go_target_center(0);goArsenal.gunTimer=0;
    for(int n=0;n<7;n++){goFight.time+=1.f/3;go_arsenal_tick(1.f/3,10,false);}
    assert(goArsenal.bullets==3&&goArsenal.side==1);
    goArsenal.rocketTimer=0;
    for(int n=0;n<61;n++){goFight.time+=1.f/60;go_arsenal_tick(1.f/60,10,false);}
    assert(goArsenal.rockets==1);int rockets=0;for(int n=0;n<32;n++)rockets+=eqRockets[n].active;assert(rockets==1);
    go_release_part(GOLEM_HEAD);goArsenal.rocketTimer=0;go_arsenal_tick(2,10,false);assert(goArsenal.rockets==1);
    fixture();goFight.target=0;goFight.aim=go_target_center(0);
    assert(go_arsenal_grab());golem_physics_prepare();assert(attached()==goCount);assert(players[MAX_PLAYERS].goldTetherHolding);
    int payload=players[MAX_PLAYERS].tetherVoxel;uint64_t id=voxels[payload].identity;
    prepare_tether_forces();assert(build_glue_cluster_indices(payload,glueClusterIndices)>1);
    float matter=goFight.shield;go_arsenal_throw();payload=enemy_voxel(id);
    assert(payload>=0&&!players[MAX_PLAYERS].tetherHolding&&go_part_alive(GOLEM_HAND_R));
    assert(v_isfinite(voxels[payload].vel)&&v_length(voxels[payload].vel)>=TETHER_THROW_IMPULSE-.01f);
    assert(goFight.shield==matter&&goArsenal.throws==1);
    puts("Independent shields, regen delay, drain lock, turret bursts, head rocket and shared gold tether passed");
}

static void shield_weapon_tests(void){
    fixture();int part=GOLEM_ARM_R,v=goCells[goMount[part]].index;uint64_t id=voxels[v].identity;
    ShotgunImpact impacts[EQ_PELLETS];for(int n=0;n<EQ_PELLETS;n++)impacts[n]=(ShotgunImpact){.id=id,.direction={0,0,1},.distance=2};
    eq_shotgun_impacts(0,impacts,EQ_PELLETS);assert(goShields[part].value==25);
    SwordWave wave={0};VoxelHit hit={.id=v};
    assert(!sword_cut_hit(0,&wave,hit));assert(goShields[part].value==0);
    recycleFrameCounter++;assert(!sword_cut_hit(0,&wave,hit));assert(enemy_voxel(id)>=0);
    fixture();v=goCells[goMount[GOLEM_HAND_R]].index;id=voxels[v].identity;
    eq_rocket_blast((EquipmentRocket){.owner=0,.pos=voxels[v].pos});
    assert(enemy_voxel(id)>=0);assert(goShields[GOLEM_HAND_R].value==0);
    fixture();goFight.target=0;goFight.aim=go_target_center(0);goArsenal.rocketTimer=0;
    for(int n=0;n<32;n++)eqRockets[n].active=true;
    go_arsenal_tick(1.1f,10,false);assert(goArsenal.rockets==0&&goArsenal.rocketTimer==1);
    memset(eqRockets,0,sizeof(eqRockets));
    puts("Shotgun aggregation, sword dedupe, shield-breaking rocket protection and pool exhaustion passed");
}
static void shield_tests(void){
    fixture();assert(goFight.shield==150&&goFight.health==200);
    apply_matter_damage(MAX_PLAYERS,0,175);assert(goFight.shield==0&&goFight.health==175);
    go_combat_prepare(7.9f);assert(goFight.shield==0);
    go_combat_prepare(.1f);assert(fabsf(goFight.shield-.5f)<.01f);
    eq_damage(0,MAX_PLAYERS,20,false,DEATH_GLOB,players[0].pos);assert(goFight.shield==0&&goFight.health==175);
    go_combat_prepare(7.9f);assert(goFight.shield==0);go_combat_prepare(.2f);assert(goFight.shield>0);
    goFight.shield=5;tetherDrainDamage=true;apply_matter_damage(MAX_PLAYERS,0,20);tetherDrainDamage=false;assert(goFight.shield==0&&goFight.health==175);
    eq_damage(0,MAX_PLAYERS,50,true,DEATH_RAIL,players[0].pos);assert(goFight.health==125&&!goCombatDead);
    eq_damage(0,MAX_PLAYERS,100,true,DEATH_ROCKET,players[0].pos);assert(goFight.health==25&&!goCombatDead);
    go_combat_damage(30,false);assert(goCombatDead&&!goRegen&&skeletonConstraintCount==0&&attached()==0);
    for(int n=0;n<120;n++)frame();assert(goCombatDead&&attached()==0);capture("defeated");
    fixture();assert(goFight.shield==150&&goFight.health==200);puts("Shield, overflow, suppression, drain, defeat/reset passed");
}
static void ray_tests(void){
    fixture();for(int n=0;n<30;n++)frame();
    Vector3 core=go_core_position();Ray ray={v_add(core,(Vector3){0,0,-8}),{0,0,1}};VoxelHit hit;
    capture("intact");assert(first_voxel_hit_detailed(ray,8,-1,false,&hit));
    printf("Chest blocker part/voxel: %d\n",hit.id);
    for(int c=goCount-1;c>=0;c--)if(goCells[c].part==GM_CHEST){int n=go_find(goCells[c].id,goCells[c].index);if(n>=0)remove_voxel_index(n);}
    golem_physics_prepare();
    assert(!first_voxel_hit_detailed(ray,7.7f,-1,false,&hit));float distance=8;
    assert(eq_actor_ray(0,ray,&distance)==MAX_PLAYERS);puts("Physical chest occlusion and exposed reactor targeting passed");
}
static void combat_tests(void){
    fixture();float initial=players[0].matter;
    for(int f=0;f<1800;f++){
        if(f==450)players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,7};
        if(f==900)players[0].pos=(Vector3){-2,BASE_EYE_HEIGHT,6};
        if(f==1350)players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,-10};
        frame();
        for(int i=0;i<sim_particle_count;i++)assert(v_isfinite(sim_particles[i]->pos));
        if(f%120==0)printf("COMBAT %.1fs state %s actions %d/%d/%d cells %d/%d position %.1f %.1f\n",f/60.f,go_combat_status(),goFight.actions[0],goFight.actions[1],goFight.actions[2],attached(),goCount,goFight.position.x,goFight.position.z);
        if(getenv("GOLEM_COMBAT_CAPTURE")&&f<600&&f%6==0)capture(TextFormat("fight-%03d",f/6));
    }
    assert(goFight.actions[0]>0&&goFight.actions[1]>0&&goFight.actions[2]>0);assert(players[0].matter==initial);
    enemyLabGod=false;players[0].matter=100;players[0].isExposed=false;players[0].invuln_timer=0;
    bool hit[MAX_COMBATANTS+MAX_PLAYERS]={0};Vector3 at=players[0].pos;
    go_sweep_hits(at,at,.65f,40,hit);go_sweep_hits(at,at,.65f,40,hit);assert(players[0].matter==60);
    players[0].invuln_timer=2;memset(hit,0,sizeof(hit));go_sweep_hits(at,at,.65f,40,hit);assert(players[0].matter==60);
    puts("Autonomous attacks, god mode, hit deduplication and invulnerability passed");
    fixture();go_release_part(GOLEM_SHIN_L);goRegen=false;for(int f=0;f<90;f++)frame();
    assert(goFight.state==GC_KNEEL);capture("kneeling");goRegen=true;
    for(int f=0;f<600;f++)frame();assert(go_part_alive(GOLEM_SHIN_L));puts("Kneeling and repair passed");
}
static void attack_hits(void){
    for(int scenario=0;scenario<2;scenario++){
        int attack=scenario<2?scenario:2;
        fixture();enemyLabGod=false;goRegen=false;
        const float distances[]={2.5f,3,6,12,20,25};float distance=distances[scenario];
        players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,10-distance};players[0].matter=100;players[0].invuln_timer=0;
        goFight.target=0;goFight.aim=go_target_center(0);go_begin_attack(attack);
        for(int a=0;a<3;a++)goFight.cooldown[a]=100;goArsenal.gunTimer=goArsenal.rocketTimer=100;
        float closest=100;Vector3 best={0};
        for(int f=0;f<240;f++){
            frame();Vector3 at=go_limb_center(attack==1?0:1);
            if(goFight.state==GC_STRIKE&&Vector3Distance(at,go_target_center(0))<closest){closest=Vector3Distance(at,go_target_center(0));best=at;}
            if(getenv("GOLEM_COMBAT_CAPTURE")&&f%4==0)capture(TextFormat("attack-%d-%03d",scenario,f/4));
        }
        printf("Attack %d actual player Matter %.1f, nearest fist %.2f\n",attack,players[0].matter,closest);printf("Best fist %.2f %.2f %.2f player %.2f %.2f %.2f\n",best.x,best.y,best.z,players[0].pos.x,players[0].pos.y,players[0].pos.z);
        assert(players[0].matter==100-goAttacks[attack].damage);
        if(attack<2)assert(attached()==goCount);
    }
}
static int compare_ms(const void *a,const void *b){float x=*(const float*)a,y=*(const float*)b;return (x>y)-(x<y);}
static void stress_test(int seconds){
    fixture();float *times=malloc((size_t)seconds*60*sizeof(float));assert(times);int peak=voxel_count;
    for(int f=0;f<seconds*60;f++){
        if(f%600==0){float distance=(f/600)%2?3:18,angle=(f/600)*1.1f;players[0].pos=(Vector3){Clamp(goFight.position.x+sinf(angle)*distance,-24,24),BASE_EYE_HEIGHT,Clamp(goFight.position.z-cosf(angle)*distance,-24,24)};}
        if(f%1200==600)go_release_part((f/1200)%2?GOLEM_SHIN_L:GOLEM_ARM_L);
        double start=GetTime();frame();times[f]=(GetTime()-start)*1000;
        peak=peak>voxel_count?peak:voxel_count;
        for(int p=0;p<sim_particle_count;p++)assert(v_isfinite(sim_particles[p]->pos));
        if(f%3600==3599)printf("STRESS %ds cells=%d peak=%d actions=%d/%d/%d\n",(f+1)/60,voxel_count,peak,goFight.actions[0],goFight.actions[1],goFight.actions[2]);
    }
    int count=seconds*60,window=count/5;float early=0,late=0;
    for(int i=0;i<window;i++){early+=times[i]/window;late+=times[count-window+i]/window;}
    qsort(times,count,sizeof(float),compare_ms);printf("STRESS %s %ds frame median=%.3fms p95=%.3fms early=%.3fms late=%.3fms peak=%d\n",physics_backend_name(physicsBackend.active),seconds,times[count/2],times[(int)(count*.95f)],early,late,peak);
    printf("ARSENAL bullets=%d rockets=%d throws=%d\n",goArsenal.bullets,goArsenal.rockets,goArsenal.throws);
    free(times);assert(peak<2000); /* Regression guard, not a production voxel cap. */
}
static void integration_tests(void){
    fixture();for(int f=0;f<20;f++)frame();
    eq_rocket_blast((EquipmentRocket){.owner=0,.pos=v_add(go_core_position(),(Vector3){0,0,-1.4f})});assert(goFight.shield==150);
    puts("Rocket protection evaluated before cover destruction passed");
    fixture();players[MAX_PLAYERS].invuln_timer=1;apply_matter_damage(MAX_PLAYERS,0,50);assert(goFight.shield==150);
    players[MAX_PLAYERS].invuln_timer=0;apply_matter_damage(MAX_PLAYERS,MAX_PLAYERS,50);assert(goFight.shield==150);
    goFight.target=0;goFight.aim=go_target_center(0);go_begin_attack(0);go_release_part(GOLEM_ARM_R);goRegen=false;frame();assert(goFight.state==GC_RECOVER);
    fixture();go_release_part(GOLEM_SHIN_L);goRegen=false;for(int f=0;f<90;f++)frame();
    Vector3 core=go_core_position(),eye={core.x,BASE_EYE_HEIGHT,core.z-8};Vector3 d=v_sub(core,eye);VoxelHit hit;
    assert(!first_voxel_hit_detailed((Ray){eye,v_norm(d)},v_length(d)-.35f,-1,false,&hit));
    puts("Core invulnerability, friendly fire, limb cancellation and kneeling opening passed");
    fixture();goFight.target=-1;players[0].respawn_timer=9999;go_combat_prepare(.02f);assert(goFight.state==GC_IDLE);
    holograms[0].life=10;holograms[0].body.pos=(Vector3){0,BASE_EYE_HEIGHT,5};holograms[0].body.respawn_timer=0;go_combat_prepare(.02f);assert(target_is_hologram(goFight.target));
    holograms[0].life=0;go_combat_prepare(.02f);assert(goFight.state==GC_IDLE);
    enemy_lab_exit();assert(!goCombatEnabled&&!goActive);puts("Target lifecycle and lab exit passed");
}
static void fallback_attack(void){
    fixture();goFight.target=0;goFight.aim=go_target_center(0);assert(go_arsenal_grab());go_begin_attack(2);
    for(int a=0;a<3;a++)goFight.cooldown[a]=100;goArsenal.gunTimer=goArsenal.rocketTimer=100;
    bool injected=false;
    for(int f=0;f<240;f++){
        if(!injected&&goFight.state==GC_STRIKE){gpu_forced_validation_failure_consumed=false;injected=true;}
        float before=goFight.time;frame();assert(fabsf(goFight.time-before-1.f/60)<.0001f);
    }
    assert(injected&&gpu_forced_validation_failure_consumed&&goArsenal.throws==1&&go_part_alive(GOLEM_HAND_R));
    puts("Forced GPU rollback: one shared tether throw, retained hand and one clock advance");
}
static void previews(void){
    const char *names[]={"part-shields","shoulder-guns","head-rocket","structure-throw"};
    for(int mode=0;mode<4;mode++){
        fixture();goRegen=false;players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,-2};
        goFight.target=0;goFight.aim=go_target_center(0);
        for(int a=0;a<3;a++)goFight.cooldown[a]=100;
        goArsenal.gunTimer=goArsenal.rocketTimer=100;
        if(mode==0)go_shield_hit(goCells[goMount[GOLEM_ARM_R]].index,60);
        if(mode==1)goArsenal.gunTimer=0;
        if(mode==2)goArsenal.rocketTimer=.5f;
        if(mode==3){assert(go_arsenal_grab());go_begin_attack(2);}
        int frames=mode==0?900:300;
        for(int f=0;f<frames;f++){
            goFight.position=(Vector3){0,.04f,10};frame();clarity_tick(1.f/60);
            if(f%6==0)capture(TextFormat("%s-%03d",names[mode],f/6));
        }
    }
}
int main(int argc,char **argv){
    setvbuf(stdout,NULL,_IONBF,0);assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,720,"Golem combat tests");init_pbd_thread_pool();assert(initialize_physics_backend());weaponRenderingReady=weapon_renderer_init();
    if(getenv("GOLEM_COMBAT_PREVIEW"))previews();
    else if(getenv("GOLEM_COMBAT_STRESS"))stress_test(atoi(getenv("GOLEM_COMBAT_STRESS")));
    else {if(getenv("FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE"))fallback_attack();arsenal_tests();shield_weapon_tests();shield_tests();ray_tests();attack_hits();combat_tests();integration_tests();}
    clear_world_voxels();assert(!goActive&&!goCombatEnabled&&skeletonConstraintCount==0);
    shutdown_pbd_thread_pool();gpu_physics_shutdown();shutdown_shard_rendering();clear_all_chunks();shutdown_world_visuals();weapon_renderer_shutdown();CloseWindow();puts("Golem combat checks passed");return 0;
}
