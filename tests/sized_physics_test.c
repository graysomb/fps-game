/* Exercise the production implementation without starting a window. */
#define main fps_game_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

static void near_vector(Vector3 a,Vector3 b,float eps) {assert(v_length(v_sub(a,b))<eps);}
static void reset_test(void) {debug_reset_world();activePlayers=0;physics_force_single_cpu=true;}
static void momenta(Vector3 *p,Vector3 *l) {
    *p=*l=(Vector3){0};
    for(int i=0;i<sim_particle_count;i++) {
        Particle *s=sim_particles[i];Vector3 impulse=v_mul(s->vel,1/s->inv_mass);
        *p=v_add(*p,impulse);*l=v_add(*l,v_cross(s->pos,impulse));
    }
}
static void make_group(bool mixed) {
    reset_test();int child[8];
    for(int k=0;k<8;k++)child[k]=debug_add_cell(-1+(k&1),5+((k>>1)&1),-1+((k>>2)&1),false,true,WHITE,63);
    if(mixed)for(int z=-1;z<1;z++)for(int x=-1;x<1;x++)debug_add_cell(x,4,z,false,true,BLUE,63);
    debug_rebuild_world_state();assert(coarse_share_child_corners());rebuild_particle_collision_metadata();
    assert(coarse_group_bind(child));
    assert(active_particle_count==(mixed?36:27));assert(sim_particle_count==(mixed?17:8));
}
static void construction_test(void) {
    reset_test();assert(add_voxel_sized(0,3,0,false,true,WHITE,0,2*VOXEL_SIZE)==0);
    Voxel *v=&voxels[0];assert(sim_particle_count==8);
    assert(v->rest_edge==1 && v->rest_volume==1 && v->particle_radius==.5f);
    for(int k=0;k<8;k++) {
        near_vector(v->particles[k]->pos,(Vector3){corner_signs[k][0]*.5f,3+corner_signs[k][1]*.5f,corner_signs[k][2]*.5f},1e-6f);
        assert(v->particles[k]->inv_mass==.125f);
    }
    assert(v->rest_max_gx-v->rest_min_gx+1==2);
    assert(add_voxel_sized(0,0,0,false,true,WHITE,0,NAN)==-1);assert(voxel_count==1);
}
static void transfer_test(void) {
    make_group(true);
    int point=coarse_sample_index(voxels[8].particles[3]);
    assert(coarseGroup.source_count[point]>1);
    Vector3 position=coarseGroup.sample[point]->pos,J={.02f,.01f,-.03f},p,l;
    coarse_point_impulse(point,J);momenta(&p,&l);
    near_vector(p,J,1e-6f);near_vector(l,v_cross(position,J),1e-6f);
    assert(v_length(coarseGroup.shape.particles[0]->vel)>0);

    /* A contact between dependent points can have overlapping sources. Its
       internal impulse must still have zero linear and angular resultant. */
    make_group(true);
    int a=-1,b=-1;
    for(int i=0;i<coarseGroup.sample_count;i++)if(coarseGroup.source_count[i]>1) {
        if(a<0)a=i;else {b=i;break;}
    }
    Vector3 n=v_norm(v_sub(coarseGroup.sample[a]->pos,coarseGroup.sample[b]->pos));
    coarse_begin_corrections();coarse_contact(a,b,n,.001f);coarse_apply_corrections();
    for(int i=0;i<sim_particle_count;i++)sim_particles[i]->vel=v_sub(sim_particles[i]->predicted_pos,sim_particles[i]->pos);
    momenta(&p,&l);near_vector(p,(Vector3){0},2e-5f);near_vector(l,(Vector3){0},2e-5f);

    /* Impulse at a free bottom corner must travel through fine VGS to parent. */
    make_group(true);voxels[8].particles[0]->vel.x=.01f;
    for(int i=0;i<10;i++)coarse_cpu_step(PBD_MAX_STEP_DT/PBD_SUBSTEPS);
    float movement=0;for(int k=0;k<8;k++)movement+=fabsf(coarseGroup.shape.particles[k]->vel.x);
    assert(movement>1e-5f);
    reset_test();assert(!coarseGroup.active && coarseGroup.sample_count==0 && active_particle_count==0);
}
static void rigid_pose_test(void) {
    make_group(true);
    const float angle=.4f;Vector3 before[COARSE_SAMPLE_CAP];
    for(int i=0;i<sim_particle_count;i++) {
        Particle *p=sim_particles[i];Vector3 q=p->pos;
        p->pos=p->predicted_pos=p->prev_pos=(Vector3){cosf(angle)*q.x-sinf(angle)*q.y+.3f,
            sinf(angle)*q.x+cosf(angle)*q.y+.7f,q.z-.2f};
    }
    coarse_refresh();
    for(int i=0;i<coarseGroup.sample_count;i++)before[i]=coarseGroup.sample[i]->predicted_pos;
    for(int pass=0;pass<3;pass++)coarse_shapes();
    for(int i=0;i<coarseGroup.sample_count;i++)near_vector(before[i],coarseGroup.sample[i]->predicted_pos,1e-5f);
}
static void broadphase_test(void) {
    reset_test();add_voxel_sized(0,3,0,false,true,WHITE,0,1);addVoxel(2,3,0,false,true,BLUE,0);
    Particle *a=voxels[0].particles[0],*b=voxels[1].particles[0];
    a->predicted_pos=(Vector3){.49f,3,0};b->predicted_pos=(Vector3){1.01f,3,0};
    a->collision_group=b->collision_group=-1;
    Particle *pair[2]={a,b};Vector3 old_a=a->predicted_pos,old_b=b->predicted_pos;
    build_particle_hash(pair,2);assert(abs(a->cell_x-b->cell_x)==2);
    gather_particle_collisions(PBD_MAX_STEP_DT,pair,2);
    assert(v_length(v_sub(a->predicted_pos,b->predicted_pos))>.52f);
    near_vector(v_add(v_mul(v_sub(a->predicted_pos,old_a),1/a->inv_mass),
                      v_mul(v_sub(b->predicted_pos,old_b),1/b->inv_mass)),(Vector3){0},1e-5f);
}
static void static_contact_test(void) {
    reset_test();int box=addVoxel(1.25f,3.25f,.25f,true,false,WHITE,0);
    int ids[16];int n=gather_static_voxels_near_point((Vector3){.01f,3.25f,.25f},1.1f,ids,16);
    assert(n==1 && ids[0]==box); /* old fixed one-cell search missed it */
    int large=add_voxel_sized(0,3,0,false,true,BLUE,0,1);
    Particle *p=voxels[large].particles[0];p->predicted_pos=(Vector3){.8f,3.25f,.25f};
    assert(push_particle_out_of_static(&voxels[box],p,.5f*p->radius));
    assert(p->predicted_pos.x<.8f);
}
int main(void) {
    construction_test();make_group(false);transfer_test();rigid_pose_test();broadphase_test();static_contact_test();reset_test();
    puts("sized construction, sharing, transfer, teardown and contact tests passed");return 0;
}
