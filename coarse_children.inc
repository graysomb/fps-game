/* One-level, CPU-only child group for the sized-voxel experiments. The children
 * are ordinary voxels. Only the parent's eight existing outer particles move
 * independently; the other child particles are trilinearly dependent. */
#define COARSE_SAMPLE_CAP 64
#define COARSE_FINE_CAP 16

typedef struct {
    bool active;
    Voxel shape;
    int child[8];
    Particle *sample[COARSE_SAMPLE_CAP];
    int sample_count;
    /* A row has one entry for an independent particle, or eight for a child. */
    Particle *source[COARSE_SAMPLE_CAP][8];
    float weight[COARSE_SAMPLE_CAP][8];
    int source_count[COARSE_SAMPLE_CAP];
    Vector3 delta[COARSE_SAMPLE_CAP];
    float load[COARSE_SAMPLE_CAP];
} CoarseVoxelGroup;
static CoarseVoxelGroup coarseGroup;

static void coarse_group_reset(void) { memset(&coarseGroup,0,sizeof(coarseGroup)); sizedFixture=false; atomic_store(&sizedWorkerCalls,0); }
static int coarse_sample_index(const Particle *p) {
    for(int i=0;i<coarseGroup.sample_count;i++)if(coarseGroup.sample[i]==p)return i;
    return -1;
}
static bool coarse_is_child(int id) {
    for(int i=0;i<8;i++)if(coarseGroup.child[i]==id)return true;
    return false;
}
static void coarse_refresh(void) {
    CoarseVoxelGroup *g=&coarseGroup;
    for(int i=0;i<g->sample_count;i++)if(g->source_count[i]>1) {
        Particle *p=g->sample[i];
        p->pos=p->prev_pos=p->predicted_pos=p->vel=(Vector3){0};
        float inverse=0;
        for(int k=0;k<g->source_count[i];k++) {
            Particle *s=g->source[i][k];float a=g->weight[i][k];
            p->pos=v_add(p->pos,v_mul(s->pos,a));
            p->prev_pos=v_add(p->prev_pos,v_mul(s->prev_pos,a));
            p->predicted_pos=v_add(p->predicted_pos,v_mul(s->predicted_pos,a));
            p->vel=v_add(p->vel,v_mul(s->vel,a));
            inverse+=a*a*s->inv_mass;
        }
        p->inv_mass=inverse; /* constraint effective mass, never integrated */
    }
}
/* Face-by-face legacy glue can leave multiple pointers at an edge where
 * three or more freshly created cells meet. Complete the pointer equivalence
 * for this bounded construction, before assigning any dependency or mass. */
static bool coarse_share_child_corners(void) {
    if(voxel_count>COARSE_FINE_CAP)return false;
    Particle *canonical[COARSE_SAMPLE_CAP];int count=0;
    for(int v=0;v<voxel_count;v++)for(int k=0;k<8;k++) {
        Particle *old=voxels[v].particles[k],*shared=NULL;
        for(int i=0;i<count;i++)if(canonical[i]->rest_cell_x==old->rest_cell_x &&
            canonical[i]->rest_cell_y==old->rest_cell_y && canonical[i]->rest_cell_z==old->rest_cell_z) {
            shared=canonical[i];break;
        }
        if(shared && shared!=old) {
            particle_retain(shared);voxels[v].particles[k]=shared;particle_release(old);
        } else if(!shared) {
            if(count==COARSE_SAMPLE_CAP)return false;
            canonical[count++]=old;
        }
    }
    collisionTopologyDirty=true;return true;
}
static bool coarse_group_bind(const int child[8]) {
    CoarseVoxelGroup *g=&coarseGroup;
    memset(g,0,sizeof(*g));memcpy(g->child,child,sizeof(g->child));
    if(voxel_count>COARSE_FINE_CAP || active_particle_count>COARSE_SAMPLE_CAP)return false;
    g->shape=voxels[child[0]];g->shape.rest_edge=2*VOXEL_SIZE;
    g->shape.rest_volume=8*VOXEL_SIZE*VOXEL_SIZE*VOXEL_SIZE;
    g->shape.simulate_dofs=true;
    for(int k=0;k<8;k++)g->shape.particles[k]=voxels[child[k]].particles[k];
    Vector3 origin=g->shape.particles[0]->pos;
    for(int i=0;i<active_particle_count;i++)g->sample[g->sample_count++]=active_particles[i];
    for(int i=0;i<g->sample_count;i++) {
        Particle *p=g->sample[i];bool inside=false,control=false;
        for(int j=0;j<8;j++) {
            if(p==g->shape.particles[j])control=true;
            for(int k=0;k<8;k++)if(p==voxels[child[j]].particles[k])inside=true;
        }
        if(!inside||control) {g->source_count[i]=1;g->source[i][0]=p;g->weight[i][0]=1;}
        else {
            g->source_count[i]=8;
            Vector3 uv=v_mul(v_sub(p->pos,origin),1/g->shape.rest_edge);
            for(int k=0;k<8;k++) {
                g->source[i][k]=g->shape.particles[k];
                g->weight[i][k]=((k&1)?uv.x:1-uv.x)*((k&2)?uv.y:1-uv.y)*((k&4)?uv.z:1-uv.z);
            }
        }
    }
    /* Assemble material corner masses once, after ordinary pointer glue. */
    float mass[COARSE_SAMPLE_CAP]={0};
    for(int v=0;v<voxel_count;v++)for(int k=0;k<8;k++) {
        int i=coarse_sample_index(voxels[v].particles[k]);if(i<0)return false;
        float r=voxels[v].rest_edge/VOXEL_SIZE;
        for(int j=0;j<g->source_count[i];j++) {
            int dst=coarse_sample_index(g->source[i][j]);
            mass[dst]+=r*r*r*g->weight[i][j];
        }
    }
    for(int i=0;i<g->sample_count;i++) {
        Particle *p=g->sample[i];
        if(g->source_count[i]>1) {sim_particles_remove(p);p->base_inv_mass=0;}
        else {if(mass[i]<=0)return false;p->base_inv_mass=p->inv_mass=1/mass[i];}
    }
    g->active=true;coarse_refresh();return true;
}
static void coarse_point_impulse(int sample,Vector3 impulse) {
    CoarseVoxelGroup *g=&coarseGroup;
    for(int k=0;k<g->source_count[sample];k++) {
        Particle *p=g->source[sample][k];
        p->vel=v_add(p->vel,v_mul(impulse,g->weight[sample][k]*p->inv_mass));
    }
    coarse_refresh();
}
static void coarse_begin_corrections(void) {
    memset(coarseGroup.delta,0,sizeof(coarseGroup.delta));
    memset(coarseGroup.load,0,sizeof(coarseGroup.load));
}
/* Expand both endpoints before computing the denominator. Shared sources cancel
 * here rather than receiving two inconsistent endpoint corrections. */
static void coarse_contact(int a,int b,Vector3 normal,float depth) {
    CoarseVoxelGroup *g=&coarseGroup;
    float gradient[COARSE_SAMPLE_CAP]={0};
    int row[2]={a,b};
    for(int side=0;side<2;side++)if(row[side]>=0) {
        int i=row[side];
        for(int k=0;k<g->source_count[i];k++)
            gradient[coarse_sample_index(g->source[i][k])]+=(side?-1:1)*g->weight[i][k];
    }
    float denom=0;
    for(int i=0;i<g->sample_count;i++)denom+=gradient[i]*gradient[i]*g->sample[i]->inv_mass;
    if(denom<=0)return;
    for(int i=0;i<g->sample_count;i++)if(gradient[i]!=0) {
        float coefficient=g->sample[i]->inv_mass*gradient[i]/denom;
        g->delta[i]=v_add(g->delta[i],v_mul(normal,depth*coefficient));
        g->load[i]+=fabsf(coefficient);
    }
}
static void coarse_apply_corrections(void) {
    CoarseVoxelGroup *g=&coarseGroup;
    /* One common relaxation retains the equal/opposite impulse of each contact.
       Independently normalizing the endpoints would introduce net momentum. */
    float load=1;
    for(int i=0;i<g->sample_count;i++)load=fmaxf(load,g->load[i]);
    for(int i=0;i<g->sample_count;i++)if(g->source_count[i]==1)
        g->sample[i]->predicted_pos=v_add(g->sample[i]->predicted_pos,v_mul(g->delta[i],1/load));
    coarse_refresh();
}
/* Evaluate unchanged legacy VGS on private particles, allowing fine voxels to
 * gather in parallel without racing on either shared corners or parent sources. */
typedef struct { Vector3 delta[COARSE_FINE_CAP][8]; } CoarseShapeJob;
static void coarse_shape_range(int begin,int end,int worker,void *user) {
    (void)worker;CoarseShapeJob *job=user;
    for(int v=begin;v<end;v++)if(!coarse_is_child(v)) {
        Voxel local=voxels[v];Particle p[8];
        for(int k=0;k<8;k++) {
            p[k]=*local.particles[k];local.particles[k]=&p[k];
            atomic_store(&p[k].corr_sum_x,0);atomic_store(&p[k].corr_sum_y,0);
            atomic_store(&p[k].corr_sum_z,0);atomic_store(&p[k].corr_weight_bits,0);
        }
        gather_voxel_shape_constraints(&local);
        for(int k=0;k<8;k++) {
            float w=bits_to_float(atomic_load(&p[k].corr_weight_bits));
            if(w>0)job->delta[v][k]=v_mul((Vector3){bits_to_float(atomic_load(&p[k].corr_sum_x)),
                bits_to_float(atomic_load(&p[k].corr_sum_y)),bits_to_float(atomic_load(&p[k].corr_sum_z))},PBD_SOR_FACTOR/w);
        }
    }
}
static void coarse_shapes(void) {
    reset_particle_accumulators();
    gather_voxel_shape_constraints(&coarseGroup.shape);
    apply_particle_accumulators();coarse_refresh();
    CoarseShapeJob job={0};
    pbd_parallel_for(0,voxel_count,coarse_shape_range,&job);
    coarse_begin_corrections();
    for(int v=0;v<voxel_count;v++)if(!coarse_is_child(v))for(int k=0;k<8;k++) {
        Vector3 d=job.delta[v][k];float length=v_length(d);
        if(length>0)coarse_contact(coarse_sample_index(voxels[v].particles[k]),-1,v_mul(d,1/length),length);
    }
    coarse_apply_corrections();
}
static bool coarse_internal_pair(int a,int b) {
    CoarseVoxelGroup *g=&coarseGroup;
    bool in_a=false,in_b=false;
    for(int j=0;j<8;j++)for(int k=0;k<8;k++) {
        Particle *p=voxels[g->child[j]].particles[k];
        in_a|=p==g->sample[a];in_b|=p==g->sample[b];
    }
    return (in_a&&in_b)||particles_are_local_structural_neighbors(g->sample[a],g->sample[b]);
}
static void coarse_collisions(void) {
    CoarseVoxelGroup *g=&coarseGroup;coarse_begin_corrections();
    /* The fixture contains at most 64 samples. Use the same sphere contact law
       as the legacy hash solver; no new collision geometry is introduced. */
    for(int a=0;a<g->sample_count;a++)for(int b=a+1;b<g->sample_count;b++) {
        if(coarse_internal_pair(a,b))continue;
        Particle *pa=g->sample[a],*pb=g->sample[b];
        Vector3 d=v_sub(pa->predicted_pos,pb->predicted_pos);float length=v_length(d);
        float depth=pa->radius+pb->radius-length;
        if(depth>0)coarse_contact(a,b,length>1e-6f?v_mul(d,1/length):(Vector3){1,0,0},.5f*COLLISION_RELAXATION*depth);
    }
    coarse_apply_corrections();
}
static void coarse_scene_contacts(float dt) {
    CoarseVoxelGroup *g=&coarseGroup;
    Vector3 before[COARSE_SAMPLE_CAP],previous[COARSE_SAMPLE_CAP];
    for(int i=0;i<g->sample_count;i++) {
        before[i]=g->sample[i]->predicted_pos;previous[i]=g->sample[i]->prev_pos;
        collision_particles[i]=g->sample[i];
    }
    collision_particle_count=g->sample_count;
    solve_static_collisions(dt); /* Existing floor, static patches and players. */
    Vector3 correction[COARSE_SAMPLE_CAP];
    for(int i=0;i<g->sample_count;i++) {
        correction[i]=v_sub(g->sample[i]->predicted_pos,before[i]);
        g->sample[i]->predicted_pos=before[i];g->sample[i]->prev_pos=previous[i];
    }
    coarse_begin_corrections();
    for(int i=0;i<g->sample_count;i++) {
        float length=v_length(correction[i]);
        if(length>0)coarse_contact(i,-1,v_mul(correction[i],1/length),length);
    }
    coarse_apply_corrections();
    /* Apply the existing floor friction to independent points only. */
    for(int i=0;i<g->sample_count;i++)if(g->source_count[i]==1) {
        Particle *p=g->sample[i];
        if(p->predicted_pos.y<=0) {
            p->prev_pos.x=p->predicted_pos.x-(p->predicted_pos.x-p->prev_pos.x)*FLOOR_TANGENTIAL_VELOCITY_RETENTION;
            p->prev_pos.z=p->predicted_pos.z-(p->predicted_pos.z-p->prev_pos.z)*FLOOR_TANGENTIAL_VELOCITY_RETENTION;
        }
    }
    coarse_refresh();
}
static void coarse_cpu_step(float dt) {
    integrate_particles(dt);coarse_refresh();coarse_collisions();
    for(int it=0;it<PBD_CONSTRAINT_ITERS;it++)coarse_shapes();
    coarse_scene_contacts(dt);
    VelocityUpdateJob job={1/dt};
    pbd_parallel_for(0,sim_particle_count,update_particle_velocity_range,&job);
    coarse_refresh();
    pbd_parallel_for(0,voxel_count,update_voxel_velocity_range,&job);
}
