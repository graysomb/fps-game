/* Small, deterministic fixtures. No lifecycle transitions can conceal motion. */
static struct {
    FILE *trace, *positions;
    int expected_particles,expected_independent;
    float expected_mass,initial_y,reference_y,reference_v;
    float mass,penetration,strain,speed,angle,drift,dependency,freefall_error;
    Vector3 center;
} sizedStats;

static bool debug_setup_sized_drop(DebugRunState *state) {
    debug_reset_world();memset(&sizedStats,0,sizeof(sizedStats));
    if(physics_backend_is_gpu(physicsBackend.requested)||physics_backend_is_gpu(physicsBackend.active)) {
        debug_add_failure(state,"sized fixtures require --physics=cpu-st or cpu-mt");return false;
    }
    sizedFixture=true;activePlayers=0;
    bool mixed=strstr(debugRunConfig.scenario,"mixed")!=NULL;
    bool group=mixed||strstr(debugRunConfig.scenario,"double")!=NULL;
    int child[8];
    if(group) {
        int base_y=mixed?5:4;
        for(int k=0;k<8;k++) {
            child[k]=debug_add_cell(-1+(k&1),base_y+((k>>1)&1),-1+((k>>2)&1),false,true,
                (Color){245,160,60,255},DEBUG_TEST_TAG);
            if(child[k]<0)return false;
        }
        if(mixed)for(int z=-1;z<1;z++)for(int x=-1;x<1;x++)
            if(debug_add_cell(x,4,z,false,true,(Color){60,180,240,255},DEBUG_TEST_TAG)<0)return false;
    } else {
        if(add_voxel_sized(0,4.5f*VOXEL_SIZE,0,false,true,(Color){245,160,60,255},0,VOXEL_SIZE)<0)return false;
        voxels[0].debugClusterTag=DEBUG_TEST_TAG;
    }
    debug_rebuild_world_state();
    if(group&&!coarse_share_child_corners())return false;
    rebuild_particle_collision_metadata();
    if(group&&!coarse_group_bind(child)) {debug_add_failure(state,"child group binding failed");return false;}
    sizedStats.expected_particles=group?(mixed?36:27):8;
    sizedStats.expected_independent=group?(mixed?17:8):8;
    sizedStats.expected_mass=8*(mixed?12:group?8:1);
    state->expected_voxels=mixed?12:group?8:1;
    state->activation_succeeded=true;
    debug_set_initial_centroid(state);
    state->camera.position=(Vector3){4,3.5f,6};state->camera.target=(Vector3){0,1.7f,0};
    char path[DEBUG_PATH_CAPACITY];
    if(!debug_join_path(path,sizeof(path),state->output_dir,"sized-trace.csv"))return false;
    sizedStats.trace=fopen(path,"w");if(!sizedStats.trace)return false;
    if(!debug_join_path(path,sizeof(path),state->output_dir,"sized-positions.csv"))return false;
    sizedStats.positions=fopen(path,"w");if(!sizedStats.positions){fclose(sizedStats.trace);sizedStats.trace=NULL;return false;}
    fputs("step,particle,x,y,z,vx,vy,vz\n",sizedStats.positions);
    fputs("step,cx,cy,cz,mass,particles,independent,penetration,strain,speed,angle,drift,dependency,freefallError\n",sizedStats.trace);
    return true;
}
static void debug_validate_sized(DebugRunState *state,int step) {
    debug_validate_all_sim_particles(state,step);
    float mass=0,speed=0,min_y=FLT_MAX,strain=0,angle=0,dependency=0;
    Vector3 center={0};
    for(int i=0;i<sim_particle_count;i++) {
        Particle *p=sim_particles[i];float m=1/p->inv_mass;
        mass+=m;center=v_add(center,v_mul(p->pos,m));speed=fmaxf(speed,v_length(p->vel));
    }
    center=v_mul(center,1/mass);
    for(int i=0;i<active_particle_count;i++) {
        Particle *p=active_particles[i];min_y=fminf(min_y,p->pos.y);
        if(sizedStats.positions)fprintf(sizedStats.positions,"%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
            step,i,p->pos.x,p->pos.y,p->pos.z,p->vel.x,p->vel.y,p->vel.z);
        if(!v_isfinite(p->pos)||!v_isfinite(p->vel))debug_add_failure(state,"nonfinite sized particle");
    }
    for(int v=0;v<voxel_count;v++) {
        float st,sh;voxel_measure_strain(&voxels[v],&st,&sh);strain=fmaxf(strain,st);
        for(int axis=0;axis<3;axis++) {
            Vector3 a={0};for(int k=0;k<8;k++)if(!(k&(1<<axis)))
                a=v_add(a,v_sub(voxels[v].particles[k|(1<<axis)]->pos,voxels[v].particles[k]->pos));
            float component=axis==0?a.x:axis==1?a.y:a.z;
            angle=fmaxf(angle,acosf(clampf(component/fmaxf(v_length(a),1e-8f),-1,1))*RAD2DEG);
        }
    }
    for(int i=0;i<coarseGroup.sample_count;i++)if(coarseGroup.source_count[i]>1) {
        if(coarseGroup.sample[i]->sim_index!=-1 || coarseGroup.sample[i]->base_inv_mass!=0)
            debug_add_failure(state,"dependent particle was independently integrated");
        Vector3 goal={0};for(int k=0;k<8;k++)goal=v_add(goal,v_mul(coarseGroup.source[i][k]->pos,coarseGroup.weight[i][k]));
        dependency=fmaxf(dependency,v_length(v_sub(goal,coarseGroup.sample[i]->pos)));
    }
    if(step==0) {
        sizedStats.initial_y=sizedStats.reference_y=center.y;
        if(fabsf(min_y-4*VOXEL_SIZE)>1e-6f)debug_add_failure(state,"incorrect initial bottom height");
    } else for(int sub=0;sub<PBD_SUBSTEPS;sub++) {
        float dt=PBD_MAX_STEP_DT/PBD_SUBSTEPS;
        sizedStats.reference_v=sizedStats.reference_v*VELOCITY_DAMPING-GRAVITY*dt;
        sizedStats.reference_y+=sizedStats.reference_v*dt;
    }
    float drift=sqrtf(center.x*center.x+center.z*center.z);
    if(step<=20)sizedStats.freefall_error=fmaxf(sizedStats.freefall_error,fabsf(center.y-sizedStats.reference_y));
    sizedStats.mass=mass;sizedStats.center=center;sizedStats.speed=speed;sizedStats.strain=strain;
    sizedStats.angle=angle;sizedStats.drift=drift;
    sizedStats.penetration=fmaxf(sizedStats.penetration,fmaxf(0,-min_y));
    sizedStats.dependency=fmaxf(sizedStats.dependency,dependency);
    if(active_particle_count!=sizedStats.expected_particles||sim_particle_count!=sizedStats.expected_independent)
        debug_add_failure(state,"sized particle counts changed (%d actual, %d independent)",active_particle_count,sim_particle_count);
    if(fabsf(mass-sizedStats.expected_mass)>1e-4f)debug_add_failure(state,"incorrect sized mass %.6g",mass);
    if(dependency>1e-5f*VOXEL_SIZE)debug_add_failure(state,"dependent particle left parent shape");
    if(sizedStats.trace)fprintf(sizedStats.trace,"%d,%.9g,%.9g,%.9g,%.9g,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
        step,center.x,center.y,center.z,mass,active_particle_count,sim_particle_count,fmaxf(0,-min_y),strain,speed,angle,drift,dependency,sizedStats.freefall_error);
}
static void debug_finalize_sized(DebugRunState *state) {
    if(sizedStats.trace){fclose(sizedStats.trace);sizedStats.trace=NULL;}
    if(sizedStats.positions){fclose(sizedStats.positions);sizedStats.positions=NULL;}
    if(physicsBackend.active==PHYSICS_BACKEND_CPU_MT && atomic_load(&sizedWorkerCalls)==0)
        debug_add_failure(state,"CPU MT fixture did not execute worker tasks");
    if(debugRunConfig.steps<600)debug_add_failure(state,"sized drops require 600 steps");
    if(sizedStats.freefall_error>1e-4f*VOXEL_SIZE)debug_add_failure(state,"freefall error %.6g",sizedStats.freefall_error);
    if(sizedStats.penetration>.02f*VOXEL_SIZE)debug_add_failure(state,"floor penetration %.6g",sizedStats.penetration);
    if(sizedStats.strain>.02f)debug_add_failure(state,"final strain %.6g",sizedStats.strain);
    if(sizedStats.speed>.01f)debug_add_failure(state,"final speed %.6g",sizedStats.speed);
    if(sizedStats.angle>1)debug_add_failure(state,"spurious rotation %.6g degrees",sizedStats.angle);
    if(sizedStats.drift>.02f*VOXEL_SIZE)debug_add_failure(state,"lateral drift %.6g",sizedStats.drift);
}
static void debug_report_sized(FILE *file) {
    if(!sizedFixture)return;
    fprintf(file,"  \"sizedWorkerCalls\": %u,\n",atomic_load(&sizedWorkerCalls));
    fprintf(file,"  \"sized\": {\"particles\": %d, \"independent\": %d, \"dependent\": %d, \"mass\": %.9g, \"centroid\": [%.9g,%.9g,%.9g], \"maxPenetration\": %.9g, \"finalStrain\": %.9g, \"finalSpeed\": %.9g, \"finalAngle\": %.9g, \"finalDrift\": %.9g, \"dependencyError\": %.9g, \"freefallError\": %.9g},\n",
        active_particle_count,sim_particle_count,active_particle_count-sim_particle_count,sizedStats.mass,
        sizedStats.center.x,sizedStats.center.y,sizedStats.center.z,sizedStats.penetration,sizedStats.strain,
        sizedStats.speed,sizedStats.angle,sizedStats.drift,sizedStats.dependency,sizedStats.freefall_error);
}
