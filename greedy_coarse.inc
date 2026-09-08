/* Scalable CPU coarse groups produced by GreedyCubeCover. Unit child voxels
 * remain the source of rendering and collision samples. */
typedef struct {
    bool active;
    Voxel *shapes; int *sizes; int group_count;
    int *voxel_group; int voxel_capacity;
    Particle **samples; int sample_count;
    Particle *(*sources)[64]; float (*weights)[64]; unsigned char *source_count;
    int *owner;
    int *dependent_order,dependent_count;
    Particle **collision_samples;int collision_count;
    Vector3 *delta; float *load;
    int child_begin,child_count,independent_count;
    double cover_ms;
} GreedyCoarseSystem;
static GreedyCoarseSystem greedyCoarse;
static bool greedyFineFixture;
static int greedyFixtureCellCount;
static Color greedy_size_color(int size){static const Color p[]={{80,185,255,255},{255,190,65,255},{105,220,125,255},{220,95,235,255},{255,105,95,255},{90,225,215,255},{245,235,95,255},{150,125,255,255}};return p[(size-1)&7];}
static int greedy_sample_index(Particle *p);
static void greedy_refresh(void);

static void greedy_coarse_reset(void){free(greedyCoarse.shapes);free(greedyCoarse.sizes);free(greedyCoarse.voxel_group);
    free(greedyCoarse.samples);free(greedyCoarse.sources);free(greedyCoarse.weights);free(greedyCoarse.source_count);free(greedyCoarse.owner);free(greedyCoarse.dependent_order);free(greedyCoarse.collision_samples);
    free(greedyCoarse.delta);free(greedyCoarse.load);memset(&greedyCoarse,0,sizeof(greedyCoarse));}
static bool greedy_solve3(float a[3][3],Vector3 b,Vector3 *out){float m[3][4]={{a[0][0],a[0][1],a[0][2],b.x},{a[1][0],a[1][1],a[1][2],b.y},{a[2][0],a[2][1],a[2][2],b.z}};
    for(int c=0;c<3;c++){int p=c;for(int r=c+1;r<3;r++)if(fabsf(m[r][c])>fabsf(m[p][c]))p=r;if(fabsf(m[p][c])<1e-9f)return false;for(int q=c;q<4;q++){float t=m[c][q];m[c][q]=m[p][q];m[p][q]=t;}float d=m[c][c];for(int q=c;q<4;q++)m[c][q]/=d;for(int r=0;r<3;r++)if(r!=c){float f=m[r][c];for(int q=c;q<4;q++)m[r][q]-=f*m[c][q];}}
    *out=(Vector3){m[0][3],m[1][3],m[2][3]};return true;}
static void greedy_coarse_materialize(void){GreedyCoarseSystem *g=&greedyCoarse;if(!g->active)return;greedy_refresh();Vector3 target_p={0},target_l={0};
    for(int i=0;i<sim_particle_count;i++){Particle *p=sim_particles[i];float m=1/p->inv_mass;Vector3 q=v_mul(p->vel,m);target_p=v_add(target_p,q);target_l=v_add(target_l,v_cross(p->pos,q));}
    float *mass=calloc((size_t)g->sample_count,sizeof(*mass));if(!mass){greedy_coarse_reset();return;}
    for(int v=0;v<voxel_count;v++)if(voxels[v].simulate)for(int k=0;k<8;k++){int i=greedy_sample_index(voxels[v].particles[k]);if(i>=0)mass[i]+=1;}
    float total=0;Vector3 center={0},now_p={0},now_l={0};for(int i=0;i<g->sample_count;i++)if(mass[i]>0){Particle *p=g->samples[i];p->base_inv_mass=p->inv_mass=1/mass[i];if(p->sim_index<0)sim_particles_add(p);total+=mass[i];center=v_add(center,v_mul(p->pos,mass[i]));Vector3 q=v_mul(p->vel,mass[i]);now_p=v_add(now_p,q);now_l=v_add(now_l,v_cross(p->pos,q));}center=v_mul(center,1/fmaxf(total,1e-9f));
    Vector3 dp=v_sub(target_p,now_p),u=v_mul(dp,1/fmaxf(total,1e-9f));Vector3 dl=v_sub(target_l,v_add(now_l,v_cross(center,dp)));float I[3][3]={{0}};
    for(int i=0;i<g->sample_count;i++)if(mass[i]>0){Vector3 r=v_sub(g->samples[i]->pos,center);float m=mass[i],r2=v_dot(r,r);I[0][0]+=m*(r2-r.x*r.x);I[1][1]+=m*(r2-r.y*r.y);I[2][2]+=m*(r2-r.z*r.z);I[0][1]-=m*r.x*r.y;I[0][2]-=m*r.x*r.z;I[1][2]-=m*r.y*r.z;}I[1][0]=I[0][1];I[2][0]=I[0][2];I[2][1]=I[1][2];Vector3 omega={0};greedy_solve3(I,dl,&omega);
    for(int i=0;i<g->sample_count;i++)if(mass[i]>0){Particle *p=g->samples[i];p->vel=v_add(p->vel,v_add(u,v_cross(omega,v_sub(p->pos,center))));}free(mass);greedy_coarse_reset();collisionTopologyDirty=true;}
static int greedy_sample_index(Particle *p){for(int i=0;i<greedyCoarse.sample_count;i++)if(greedyCoarse.samples[i]==p)return i;return -1;}
static void greedy_refresh(void){GreedyCoarseSystem *g=&greedyCoarse;for(int q=0;q<g->dependent_count;q++){int i=g->dependent_order[q];
    Particle *p=g->samples[i];p->pos=p->prev_pos=p->predicted_pos=p->vel=(Vector3){0};float inv=0;
    for(int k=0;k<g->source_count[i];k++){Particle *s=g->sources[i][k];float a=g->weights[i][k];p->pos=v_add(p->pos,v_mul(s->pos,a));
        p->prev_pos=v_add(p->prev_pos,v_mul(s->prev_pos,a));p->predicted_pos=v_add(p->predicted_pos,v_mul(s->predicted_pos,a));p->vel=v_add(p->vel,v_mul(s->vel,a));inv+=a*a*s->inv_mass;}p->inv_mass=inv;}}
static bool greedy_share_range(int begin,int count){int cap=count*8,n=0;Particle **list=malloc((size_t)cap*sizeof(*list));int *types=malloc((size_t)cap*sizeof(*types));if(cap&&(!list||!types)){free(list);free(types);return false;}
    for(int v=begin;v<begin+count;v++)for(int k=0;k<8;k++){Particle *old=voxels[v].particles[k],*same=NULL;for(int i=0;i<n;i++)
        if(types[i]==voxels[v].type&&list[i]->rest_cell_x==old->rest_cell_x&&list[i]->rest_cell_y==old->rest_cell_y&&list[i]->rest_cell_z==old->rest_cell_z){same=list[i];break;}
        if(same&&same!=old){particle_retain(same);voxels[v].particles[k]=same;particle_release(old);}else if(!same){list[n]=old;types[n++]=voxels[v].type;}}
    free(list);free(types);collisionTopologyDirty=true;return true;}
static int greedy_voxel_at(int begin,int count,int x,int y,int z){for(int v=begin;v<begin+count;v++)if(voxels[v].orig_min_gx==x&&voxels[v].orig_min_gy==y&&voxels[v].orig_min_gz==z)return v;return -1;}
static bool greedy_group_precedes(int a,int b){if(b<0)return true;GreedyCoarseSystem *g=&greedyCoarse;if(g->sizes[a]!=g->sizes[b])return g->sizes[a]>g->sizes[b];
    Voxel *x=&g->shapes[a],*y=&g->shapes[b];if(x->orig_min_gx!=y->orig_min_gx)return x->orig_min_gx<y->orig_min_gx;
    if(x->orig_min_gy!=y->orig_min_gy)return x->orig_min_gy<y->orig_min_gy;return x->orig_min_gz<y->orig_min_gz;}
static void greedy_add_mass(float *mass,int sample,float amount,int depth){GreedyCoarseSystem *g=&greedyCoarse;if(sample<0||depth>g->group_count+1)return;
    if(g->source_count[sample]<=1){mass[sample]+=amount;return;}for(int k=0;k<g->source_count[sample];k++)greedy_add_mass(mass,greedy_sample_index(g->sources[sample][k]),amount*g->weights[sample][k],depth+1);}
static bool greedy_coarse_bind(const GreedyCubeCover *cover,const UnitVoxelBuffer *buffer,int begin){
    greedy_coarse_reset();GreedyCoarseSystem *g=&greedyCoarse;if(!cover||!buffer||cover->cube_count<=0)return false;
    g->group_count=cover->cube_count;g->child_begin=begin;g->child_count=buffer->count;g->voxel_capacity=voxel_count;
    g->shapes=calloc((size_t)g->group_count,sizeof(*g->shapes));g->sizes=calloc((size_t)g->group_count,sizeof(*g->sizes));g->voxel_group=malloc((size_t)voxel_count*sizeof(*g->voxel_group));
    if(!g->shapes||!g->sizes||!g->voxel_group)goto fail;for(int i=0;i<voxel_count;i++)g->voxel_group[i]=-1;
    for(int c=0;c<g->group_count;c++){GreedyCube q=cover->cubes[c];g->sizes[c]=q.size;int first=greedy_voxel_at(begin,buffer->count,q.x,q.y,q.z);if(first<0)goto fail;
        Voxel shape=voxels[first];shape.rest_edge=q.size*VOXEL_SIZE;shape.rest_volume=shape.rest_edge*shape.rest_edge*shape.rest_edge;shape.simulate_dofs=true;
        shape.orig_min_gx=shape.rest_min_gx=q.x;shape.orig_min_gy=shape.rest_min_gy=q.y;shape.orig_min_gz=shape.rest_min_gz=q.z;
        shape.orig_max_gx=shape.rest_max_gx=q.x+q.size-1;shape.orig_max_gy=shape.rest_max_gy=q.y+q.size-1;shape.orig_max_gz=shape.rest_max_gz=q.z+q.size-1;
        for(int k=0;k<8;k++){int v=greedy_voxel_at(begin,buffer->count,q.x+((k&1)?q.size-1:0),q.y+((k&2)?q.size-1:0),q.z+((k&4)?q.size-1:0));if(v<0)goto fail;shape.particles[k]=voxels[v].particles[k];}g->shapes[c]=shape;}
    for(int i=0;i<buffer->count;i++){int c=cover->cell_cube[i];if(c<0||c>=g->group_count)goto fail;g->voxel_group[begin+i]=c;}
    g->sample_count=active_particle_count;size_t n=(size_t)g->sample_count;g->samples=malloc(n*sizeof(*g->samples));g->sources=calloc(n,sizeof(*g->sources));g->weights=calloc(n,sizeof(*g->weights));
    g->source_count=calloc(n,sizeof(*g->source_count));g->delta=calloc(n,sizeof(*g->delta));g->load=calloc(n,sizeof(*g->load));g->owner=malloc(n*sizeof(*g->owner));g->dependent_order=malloc(n*sizeof(*g->dependent_order));if(!g->samples||!g->sources||!g->weights||!g->source_count||!g->delta||!g->load||!g->owner||!g->dependent_order)goto fail;
    memcpy(g->samples,active_particles,n*sizeof(*g->samples));int *owner=g->owner;for(int i=0;i<g->sample_count;i++)owner[i]=-1;
    for(int v=begin;v<begin+buffer->count;v++){int c=g->voxel_group[v];for(int k=0;k<8;k++){int i=greedy_sample_index(voxels[v].particles[k]);if(i>=0&&greedy_group_precedes(c,owner[i]))owner[i]=c;}}
    for(int i=0;i<g->sample_count;i++){Particle *p=g->samples[i];int c=owner[i];bool control=c<0;if(c>=0)for(int k=0;k<8;k++)if(g->shapes[c].particles[k]==p)control=true;
        if(control){g->source_count[i]=1;g->sources[i][0]=p;g->weights[i][0]=1;}else{g->source_count[i]=8;Voxel *s=&g->shapes[c];Vector3 uv=v_mul(v_sub(p->pos,s->particles[0]->pos),1/s->rest_edge);
            for(int k=0;k<8;k++){g->sources[i][k]=s->particles[k];g->weights[i][k]=((k&1)?uv.x:1-uv.x)*((k&2)?uv.y:1-uv.y)*((k&4)?uv.z:1-uv.z);}}}
    /* Flatten the acyclic ownership graph. Corrections and mass then reach
     * independent controls in one sparse pass instead of recursively branching
     * for every contact. */
    for(int c=0;c<g->group_count;c++)for(int i=0;i<g->sample_count;i++)if(g->owner[i]==c&&g->source_count[i]>1){
        Particle *flat[64]={0};float fw[64]={0};int fn=0;
        for(int k=0;k<8;k++){int si=greedy_sample_index(g->sources[i][k]);float outer=g->weights[i][k];if(si<0)continue;
            for(int q=0;q<g->source_count[si];q++){Particle *p=g->sources[si][q];float w=outer*g->weights[si][q];int j=0;while(j<fn&&flat[j]!=p)j++;if(j==fn){if(fn>=64)goto fail;flat[fn]=p;fw[fn++]=0;}fw[j]+=w;}}
        g->source_count[i]=(unsigned char)fn;for(int k=0;k<fn;k++){g->sources[i][k]=flat[k];g->weights[i][k]=fw[k];}g->dependent_order[g->dependent_count++]=i;}
    float *mass=calloc(n,sizeof(*mass));if(!mass)goto fail;for(int v=0;v<voxel_count;v++)if(voxels[v].simulate)for(int k=0;k<8;k++){int i=greedy_sample_index(voxels[v].particles[k]);if(i>=0)greedy_add_mass(mass,i,1,0);}
    for(int i=0;i<g->sample_count;i++)if(g->source_count[i]>1){sim_particles_remove(g->samples[i]);g->samples[i]->base_inv_mass=0;}else if(mass[i]>0){g->samples[i]->base_inv_mass=g->samples[i]->inv_mass=1/mass[i];g->independent_count++;}free(mass);
    g->active=true;greedy_refresh();collisionTopologyDirty=true;return true;
fail:greedy_coarse_reset();return false;}
static void greedy_begin(void){memset(greedyCoarse.delta,0,(size_t)greedyCoarse.sample_count*sizeof(*greedyCoarse.delta));memset(greedyCoarse.load,0,(size_t)greedyCoarse.sample_count*sizeof(*greedyCoarse.load));}
static void greedy_correction_recursive(int sample,Vector3 d,float coefficient,int depth){GreedyCoarseSystem *g=&greedyCoarse;if(sample<0||depth>g->group_count+1)return;
    if(g->source_count[sample]<=1){g->delta[sample]=v_add(g->delta[sample],v_mul(d,coefficient));g->load[sample]+=fabsf(coefficient);return;}
    for(int k=0;k<g->source_count[sample];k++)greedy_correction_recursive(greedy_sample_index(g->sources[sample][k]),d,coefficient*g->weights[sample][k],depth+1);}
static void greedy_correction(int sample,Vector3 d){if(sample>=0&&v_length(d)>0)greedy_correction_recursive(sample,d,1,0);}
static void greedy_point_impulse(int sample,Vector3 impulse){GreedyCoarseSystem *g=&greedyCoarse;if(sample<0)return;for(int k=0;k<g->source_count[sample];k++){Particle *p=g->sources[sample][k];p->vel=v_add(p->vel,v_mul(impulse,g->weights[sample][k]*p->inv_mass));}greedy_refresh();}
static void greedy_apply(void){GreedyCoarseSystem *g=&greedyCoarse;for(int i=0;i<g->sample_count;i++)if(g->source_count[i]==1&&g->load[i]>0)g->samples[i]->predicted_pos=v_add(g->samples[i]->predicted_pos,v_mul(g->delta[i],1/g->load[i]));greedy_refresh();}
static void greedy_shape_delta(Voxel *shape,Vector3 d[8]){Voxel local=*shape;Particle p[8];for(int k=0;k<8;k++){p[k]=*local.particles[k];local.particles[k]=&p[k];atomic_store(&p[k].corr_sum_x,0);atomic_store(&p[k].corr_sum_y,0);atomic_store(&p[k].corr_sum_z,0);atomic_store(&p[k].corr_weight_bits,0);}gather_voxel_shape_constraints(&local);for(int k=0;k<8;k++){float w=bits_to_float(atomic_load(&p[k].corr_weight_bits));d[k]=w>0?v_mul((Vector3){bits_to_float(atomic_load(&p[k].corr_sum_x)),bits_to_float(atomic_load(&p[k].corr_sum_y)),bits_to_float(atomic_load(&p[k].corr_sum_z))},PBD_SOR_FACTOR/w):(Vector3){0};}}
static void greedy_shapes(void){GreedyCoarseSystem *g=&greedyCoarse;greedy_begin();for(int c=0;c<g->group_count;c++){Vector3 d[8];greedy_shape_delta(&g->shapes[c],d);for(int k=0;k<8;k++)greedy_correction(greedy_sample_index(g->shapes[c].particles[k]),d[k]);}
    for(int v=0;v<voxel_count;v++)if(voxels[v].simulate&&(v>=g->voxel_capacity||g->voxel_group[v]<0)){Vector3 d[8];greedy_shape_delta(&voxels[v],d);for(int k=0;k<8;k++)greedy_correction(greedy_sample_index(voxels[v].particles[k]),d[k]);}greedy_apply();}
static void greedy_build_collision_samples(void){GreedyCoarseSystem *g=&greedyCoarse;free(g->collision_samples);g->collision_samples=malloc((size_t)g->sample_count*sizeof(*g->collision_samples));g->collision_count=0;if(!g->collision_samples)return;
    unsigned char *seen=calloc((size_t)g->sample_count,1);if(!seen)return;for(int v=g->child_begin;v<g->child_begin+g->child_count;v++){bool shell=false;for(int f=0;f<6;f++)shell|=voxels[v].surface[f];if(!shell)continue;for(int k=0;k<8;k++){int i=greedy_sample_index(voxels[v].particles[k]);if(i>=0&&!seen[i]){seen[i]=1;g->collision_samples[g->collision_count++]=g->samples[i];}}}free(seen);}
static void greedy_collisions(void){GreedyCoarseSystem *g=&greedyCoarse;int count=g->collision_count;Vector3 *before=malloc((size_t)count*sizeof(*before));if(!before)return;for(int i=0;i<count;i++){particle_snapshot[i]=g->collision_samples[i];before[i]=particle_snapshot[i]->predicted_pos;}
    build_particle_hash(particle_snapshot,count);gather_particle_collisions(PBD_MAX_STEP_DT,particle_snapshot,count);greedy_begin();for(int i=0;i<count;i++){Vector3 d=v_sub(particle_snapshot[i]->predicted_pos,before[i]);particle_snapshot[i]->predicted_pos=before[i];greedy_correction(greedy_sample_index(particle_snapshot[i]),d);}free(before);greedy_apply();}
static void greedy_scene(float dt){GreedyCoarseSystem *g=&greedyCoarse;int count=g->collision_count;Vector3 *before=malloc((size_t)count*sizeof(*before));if(!before)return;for(int i=0;i<count;i++){collision_particles[i]=g->collision_samples[i];before[i]=collision_particles[i]->predicted_pos;}collision_particle_count=count;solve_static_collisions(dt);greedy_begin();for(int i=0;i<count;i++){Particle *p=collision_particles[i];Vector3 d=v_sub(p->predicted_pos,before[i]);p->predicted_pos=before[i];greedy_correction(greedy_sample_index(p),d);}free(before);greedy_apply();}
static void greedy_coarse_cpu_step(float dt){integrate_particles(dt);greedy_refresh();greedy_collisions();for(int i=0;i<PBD_CONSTRAINT_ITERS;i++)greedy_shapes();for(int i=0;i<3;i++)greedy_scene(dt);VelocityUpdateJob job={1/dt};pbd_parallel_for(0,sim_particle_count,update_particle_velocity_range,&job);greedy_refresh();pbd_parallel_for(0,voxel_count,update_voxel_velocity_range,&job);}

static bool emit_greedy_dynamic_from_buffer(UnitVoxelBuffer *buffer){if(!buffer||buffer->count<=0)return false;if(physics_backend_is_gpu(physicsBackend.requested)||physics_backend_is_gpu(physicsBackend.active)){fprintf(stderr,"--activation-cubes=greedy requires CPU ST or MT physics\n");return false;}
    GreedyCell *cells=malloc((size_t)buffer->count*sizeof(*cells));if(!cells)return false;for(int i=0;i<buffer->count;i++)cells[i]=(GreedyCell){buffer->voxels[i].gx,buffer->voxels[i].gy,buffer->voxels[i].gz,buffer->voxels[i].type,i};
    GreedyCubeCover cover={0};double started=GetTime();bool ok=greedy_cube_cover_build(cells,buffer->count,&cover);free(cells);if(!ok)return false;cover.milliseconds=(GetTime()-started)*1000;
    remove_buffered_static_voxels(buffer);int begin=voxel_count;if(emit_unit_voxels_from_units(buffer,false,true,-1)!=buffer->count||!greedy_share_range(begin,buffer->count)||!greedy_coarse_bind(&cover,buffer,begin)){while(voxel_count>begin)remove_voxel_index(voxel_count-1);emit_static_voxels_from_units(buffer);greedy_cube_cover_free(&cover);return false;}
    greedyCoarse.cover_ms=cover.milliseconds;greedy_cube_cover_free(&cover);rebuild_voxel_hash();rebuild_all_voxel_surfaces();greedy_build_collision_samples();rebuild_glue_constraints();refresh_static_voxel_beliefs();meshDirty=true;return true;}
