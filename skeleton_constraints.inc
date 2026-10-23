/* Dynamic mounting/weld constraints. Included after the CPU PBD helpers.
 * A constrained corner follows either a world target or an affine combination
 * of another cell's corners. Every participating particle retains its mass. */
typedef struct SkeletonConstraint {
    Particle *corner, *support[8];
    float weight[8];
    Vector3 target,previous_target;
    float compliance, limit;
} SkeletonConstraint;
typedef struct SkeletonGpuConstraint {
    uint32_t corner[12];
    float weight[8], target[4], parameters[4], previous_target[4];
} SkeletonGpuConstraint;
_Static_assert(sizeof(SkeletonGpuConstraint)==128,"Skeleton constraint GPU layout");
static SkeletonConstraint *skeletonConstraints;
static int skeletonConstraintCount,skeletonConstraintCapacity;
static void golem_physics_tick(float dt);
static void golem_physics_prepare(void);
static void golem_physics_reset(void);
static bool skeleton_reserve(int count){
    if(count<=skeletonConstraintCapacity)return true;
    SkeletonConstraint *next=realloc(skeletonConstraints,(size_t)count*sizeof(*next));
    if(!next)return false;
    skeletonConstraints=next;skeletonConstraintCapacity=count;return true;
}
static void skeleton_solve_cpu(float dt){
    if(!skeletonConstraintCount)return;
    reset_particle_accumulators();
    for(int n=0;n<skeletonConstraintCount;n++){
        SkeletonConstraint *c=&skeletonConstraints[n];Particle *p=c->corner;
        if(!p||!p->active)continue;
        Vector3 target=c->target,previous=c->previous_target;float denom=p->inv_mass+c->compliance/(dt*dt);
        for(int k=0;k<8;k++)if(c->support[k]){
            target=v_add(target,v_mul(c->support[k]->predicted_pos,c->weight[k]));
            previous=v_add(previous,v_mul(c->support[k]->pos,c->weight[k]));
            denom+=c->support[k]->inv_mass*c->weight[k]*c->weight[k];
        }
        if(denom<=0)continue;
        Vector3 relative=v_sub(v_sub(p->predicted_pos,p->pos),v_sub(target,previous));
        Vector3 delta=v_sub(v_sub(target,p->predicted_pos),v_mul(relative,20*dt/12));float length=v_length(delta);
        if(length>c->limit)delta=v_mul(delta,c->limit/length);
        delta=v_mul(delta,1/denom);
        accumulate_particle_correction(p,v_mul(delta,p->inv_mass),1);
        for(int k=0;k<8;k++)if(c->support[k])
            accumulate_particle_correction(c->support[k],v_mul(delta,-c->weight[k]*c->support[k]->inv_mass),1);
    }
    apply_particle_accumulators();
}
