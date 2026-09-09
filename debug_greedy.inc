static struct {
    UnitVoxelBuffer buffer;
    int ticks,activation_step,original_particles;
    float max_dependency,max_penetration,max_strain,final_strain;
    int histogram[32];
    FILE *positions;
    bool fine;
    double physics_seconds;
    int physics_steps;
} greedyDebug;

static bool debug_setup_greedy_irregular(DebugRunState *state){
    debug_reset_world();memset(&greedyDebug,0,sizeof(greedyDebug));memset(&greedyStageProfile,0,sizeof(greedyStageProfile));greedyStageProfile.enabled=true;greedyFixture=true;greedyDebug.fine=strstr(debugRunConfig.scenario,"-fine")!=NULL;greedyFineFixture=greedyDebug.fine;greedyActivationCubes=!greedyDebug.fine;activePlayers=0;greedyDebug.activation_step=20;
    if(physics_backend_is_gpu(physicsBackend.requested)||physics_backend_is_gpu(physicsBackend.active)){debug_add_failure(state,"greedy activation requires CPU ST or MT");return false;}
    for(int z=-6;z<6;z++)for(int y=8;y<20;y++)for(int x=-6;x<6;x++){
        if(x>=-2&&x<2&&y>=12&&y<16)continue;if(debug_add_cell(x,y,z,true,false,(Color){115,125,138,255},DEBUG_TEST_TAG)<0)return false;}
    for(int z=-2;z<=2;z++)for(int y=8;y<13;y++)for(int x=6;x<11;x++)if(debug_add_cell(x,y,z,true,false,(Color){115,125,138,255},DEBUG_TEST_TAG)<0)return false;
    for(int z=-1;z<=1;z++)for(int y=8;y<11;y++)for(int x=11;x<14;x++)if(debug_add_cell(x,y,z,true,false,(Color){115,125,138,255},DEBUG_TEST_TAG)<0)return false;
    for(int x=14;x<18;x++)if(debug_add_cell(x,8,0,true,false,(Color){115,125,138,255},DEBUG_TEST_TAG)<0)return false;
    for(int i=0;i<voxel_count;i++){Voxel *v=&voxels[i];v->pendingActivation=true;if(!unit_voxel_buffer_push(&greedyDebug.buffer,v->gx,v->gy,v->gz,v->color,v->type,v->fixed,i,v->debugClusterTag,-1))return false;}
    greedyDebug.original_particles=greedyDebug.buffer.count*8;greedyFixtureCellCount=greedyDebug.buffer.count;state->expected_voxels=greedyDebug.buffer.count;state->activation_succeeded=true;
    char path[DEBUG_PATH_CAPACITY];if(!debug_join_path(path,sizeof(path),state->output_dir,"greedy-positions.csv"))return false;greedyDebug.positions=fopen(path,"w");if(!greedyDebug.positions)return false;
    fputs("step,particle,x,y,z,vx,vy,vz\n",greedyDebug.positions);debug_set_initial_centroid(state);state->camera.position=(Vector3){13,10,17};state->camera.target=(Vector3){2,6,0};return true;
}
static void debug_greedy_step(DebugRunState *state){(void)state;greedyDebug.ticks++;if(greedyDebug.ticks==greedyDebug.activation_step){
    if(greedyDebug.fine){remove_buffered_static_voxels(&greedyDebug.buffer);emit_unit_voxels_from_units(&greedyDebug.buffer,false,true,-1);rebuild_voxel_hash();rebuild_all_voxel_surfaces();rebuild_glue_constraints();}
    else if(!emit_greedy_dynamic_from_buffer(&greedyDebug.buffer))return;memset(greedyDebug.histogram,0,sizeof(greedyDebug.histogram));
    for(int c=0;c<greedyCoarse.group_count;c++){int size=greedyCoarse.sizes[c];if(size<32)greedyDebug.histogram[size]++;}
    for(int v=greedyCoarse.child_begin;v<greedyCoarse.child_begin+greedyCoarse.child_count;v++){int c=greedyCoarse.voxel_group[v];if(c>=0)voxels[v].color=greedy_size_color(greedyCoarse.sizes[c]);}
    rebuild_all_voxel_surfaces();}if((greedyCoarse.active||(greedyDebug.fine&&greedyDebug.ticks>=greedyDebug.activation_step))&&greedyDebug.ticks>greedyDebug.activation_step){physics_force_single_cpu=physicsBackend.active==PHYSICS_BACKEND_CPU_ST;double started=GetTime();simulate_voxel_pbd_cpu_steps(PBD_MAX_STEP_DT/PBD_SUBSTEPS,PBD_SUBSTEPS);greedyDebug.physics_seconds+=GetTime()-started;greedyDebug.physics_steps++;physics_force_single_cpu=false;}}
static void debug_validate_greedy(DebugRunState *state,int step){if(greedyDebug.fine){if(greedyDebug.ticks>=greedyDebug.activation_step)debug_validate_all_sim_particles(state,step);return;}if(!greedyCoarse.active)return;debug_validate_all_sim_particles(state,step);float dependency=0,miny=FLT_MAX,strain=0;
    for(int i=0;i<greedyCoarse.sample_count;i++){Particle *p=greedyCoarse.samples[i];miny=fminf(miny,p->pos.y);if(greedyDebug.positions)fprintf(greedyDebug.positions,"%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",step,i,p->pos.x,p->pos.y,p->pos.z,p->vel.x,p->vel.y,p->vel.z);if(greedyCoarse.source_count[i]>1){Vector3 goal={0};for(int k=0;k<greedyCoarse.source_count[i];k++)goal=v_add(goal,v_mul(greedyCoarse.sources[i][k]->pos,greedyCoarse.weights[i][k]));dependency=fmaxf(dependency,v_length(v_sub(goal,p->pos)));}}
    for(int c=0;c<greedyCoarse.group_count;c++){float a,b;voxel_measure_strain(&greedyCoarse.shapes[c],&a,&b);strain=fmaxf(strain,a);}
    greedyDebug.final_strain=strain;greedyDebug.max_dependency=fmaxf(greedyDebug.max_dependency,dependency);greedyDebug.max_penetration=fmaxf(greedyDebug.max_penetration,fmaxf(0,-miny));greedyDebug.max_strain=fmaxf(greedyDebug.max_strain,strain);
    if(dependency>VOXEL_SIZE*1e-5f)debug_add_failure(state,"greedy dependency error %.9g",dependency);
    if(step==debugRunConfig.steps){int distinct=0;for(int i=1;i<32;i++)distinct+=greedyDebug.histogram[i]>0;if(distinct<3)debug_add_failure(state,"greedy cover produced only %d sizes",distinct);
        if(greedyCoarse.independent_count>=greedyDebug.original_particles)debug_add_failure(state,"greedy cover did not reduce independent particles");if(greedyDebug.max_penetration>VOXEL_SIZE*.02f)debug_add_failure(state,"greedy floor penetration %.9g",greedyDebug.max_penetration);if(greedyDebug.max_strain>.5f)debug_add_failure(state,"greedy transient strain %.9g",greedyDebug.max_strain);if(greedyDebug.final_strain>.02f)debug_add_failure(state,"greedy final strain %.9g",greedyDebug.final_strain);if(physicsBackend.active==PHYSICS_BACKEND_CPU_MT&&atomic_load(&sizedWorkerCalls)==0)debug_add_failure(state,"greedy CPU MT fixture did not execute a worker job");if(greedyDebug.positions){fclose(greedyDebug.positions);greedyDebug.positions=NULL;}}}
static void debug_report_greedy(FILE *f){if(!greedyFixture)return;fprintf(f,"  \"greedy\": {\"mode\": \"%s\", \"activated\": %s, \"cells\": %d, \"groups\": %d, \"samples\": %d, \"independent\": %d, \"dependent\": %d, \"coverMs\": %.9g, \"physicsSteps\": %d, \"physicsTotalMs\": %.9g, \"physicsMeanMs\": %.9g, \"integrationMs\": %.9g, \"vgsMs\": %.9g, \"dynamicContactMs\": %.9g, \"staticContactMs\": %.9g, \"interfaceRefreshMs\": %.9g, \"passiveRefreshMs\": %.9g, \"maxDependency\": %.9g, \"maxPenetration\": %.9g, \"maxStrain\": %.9g, \"finalStrain\": %.9g, \"workerCalls\": %u, \"histogram\": {",
    greedyDebug.fine?"fine":"greedy",greedyCoarse.active?"true":"false",greedyDebug.buffer.count,greedyCoarse.group_count,greedyDebug.fine?active_particle_count:greedyCoarse.sample_count,greedyDebug.fine?sim_particle_count:greedyCoarse.independent_count,greedyDebug.fine?0:greedyCoarse.sample_count-greedyCoarse.independent_count,greedyCoarse.cover_ms,greedyDebug.physics_steps,greedyDebug.physics_seconds*1000,greedyDebug.physics_steps?greedyDebug.physics_seconds*1000/greedyDebug.physics_steps:0,greedyStageProfile.integration_ms,greedyStageProfile.vgs_ms,greedyStageProfile.dynamic_ms,greedyStageProfile.static_ms,greedyStageProfile.interface_ms,greedyStageProfile.child_ms,greedyDebug.max_dependency,greedyDebug.max_penetration,greedyDebug.max_strain,greedyDebug.final_strain,atomic_load(&sizedWorkerCalls));
    bool first=true;for(int i=1;i<32;i++)if(greedyDebug.histogram[i]){fprintf(f,"%s\"%d\": %d",first?"":", ",i,greedyDebug.histogram[i]);first=false;}fputs("}},\n",f);}
