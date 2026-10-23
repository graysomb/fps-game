/* Local-only articulated body fixture, backed by the normal particle solver. */
#include "golem_renderer.c"
typedef struct GoCell {Vector3 center,size;Color color;int part,index;uint64_t id;float missing;} GoCell;
typedef struct GoLink {int a,b,face;float weight[4][8];int corner[4];} GoLink;
typedef struct GoWeld {GoLink definition;uint64_t a,b;int ia,ib;} GoWeld;
static GoCell *goCells;
static GoLink *goLinks;
static GoWeld *goWelds;
static int goCount,goPartStart,goLinkCount,goWeldCount,goWeldCapacity;
static bool goActive,goRegen=true,goOverlay,goAllocationFailed,goConnectNeeded;
static GolemPose goPose,goPreviousPose;
static float goLastDt;
static GolemInput goInput;
static float goClock,goError,goOverload[GM_COUNT],goLost[GM_COUNT],goLoad;
static int goMount[GM_COUNT],goAnchor[GM_COUNT][4],goAnchorCorner[GM_COUNT][4];
static int goAction;
static float goActionStart;
static float goLastHandLaunch=-1;
static const float goScale=1.5f;
typedef enum GoRepairPhase { GO_REPAIR_IDLE, GO_REPAIR_WAIT, GO_REPAIR_BLOCKED, GO_REPAIR_BUILD, GO_REPAIR_RECOVER } GoRepairPhase;
typedef struct GoRepair {GoRepairPhase phase;float since;} GoRepair;
static GoRepair goRepair[GM_COUNT];
#ifdef FIREFIGHT_TESTING
static int goRepairFailAfter=-1;
#endif
static struct {
    GoCell *cells;int capacity;
    GolemInput input;GolemPose pose,previous;
    float clock,lastDt,overload[GM_COUNT],lost[GM_COUNT];
    GoRepair repair[GM_COUNT];
} goRollback;
static bool golem_checkpoint_save(void){
    if(!goActive)return true;
    if(goRollback.capacity<goCount){
        GoCell *cells=realloc(goRollback.cells,(size_t)goCount*sizeof(*cells));
        if(!cells)return false;goRollback.cells=cells;goRollback.capacity=goCount;
    }
    memcpy(goRollback.cells,goCells,(size_t)goCount*sizeof(*goCells));
    goRollback.input=goInput;goRollback.pose=goPose;goRollback.previous=goPreviousPose;
    memcpy(goRollback.repair,goRepair,sizeof(goRepair));
    goRollback.clock=goClock;goRollback.lastDt=goLastDt;
    memcpy(goRollback.overload,goOverload,sizeof(goOverload));memcpy(goRollback.lost,goLost,sizeof(goLost));
    return true;
}
static void golem_checkpoint_restore(void){
    if(!goActive)return;
    memcpy(goCells,goRollback.cells,(size_t)goCount*sizeof(*goCells));
    goInput=goRollback.input;goPose=goRollback.pose;goPreviousPose=goRollback.previous;
    memcpy(goRepair,goRollback.repair,sizeof(goRepair));
    goClock=goRollback.clock;goLastDt=goRollback.lastDt;
    memcpy(goOverload,goRollback.overload,sizeof(goOverload));memcpy(goLost,goRollback.lost,sizeof(goLost));
    for(int i=0;i<goCount;i++)if(goCells[i].index>=0)voxels[goCells[i].index].persistentConstruction=true;
    collisionTopologyDirty=true;
}
static void go_collect(Vector3 c,Vector3 size,Color color){
    GoCell *next=realloc(goCells,(size_t)(goCount+1)*sizeof(*next));
    if(!next){goAllocationFailed=true;return;}goCells=next;
    goCells[goCount++]=(GoCell){.center=c,.size=size,.color=color,.index=-1};
}
static void go_finish(int part){for(int i=goPartStart;i<goCount;i++)goCells[i].part=part;goPartStart=goCount;}
static int go_find(uint64_t id,int cached){
    if(!id)return -1;
    if(cached>=0&&cached<voxel_count&&voxels[cached].identity==id)return cached;
    for(int i=0;i<voxel_count;i++)if(voxels[i].identity==id)return i;
    return -1;
}
static Matrix go_matrix(int part){return part<GOLEM_JOINT_COUNT?goPose.joint[part]:part==GM_CHEST?goPose.chest:part==GM_CORE?goPose.core:goPose.eye;}
static Matrix go_smooth_matrix(Matrix a,Matrix b,float amount){
    Vector3 pa,pb,sa,sb;Quaternion qa,qb;
    MatrixDecompose(a,&pa,&qa,&sa);MatrixDecompose(b,&pb,&qb,&sb);
    Vector3 position=Vector3Lerp(pa,pb,amount);
    /* Compose scale, rotation and translation with the raylib 5.5 API. */
    Matrix result=MatrixMultiply(MatrixScale(sb.x,sb.y,sb.z),QuaternionToMatrix(QuaternionSlerp(qa,qb,amount)));
    result.m12=position.x;result.m13=position.y;result.m14=position.z;
    return result;
}
static Vector3 go_corner(GoCell *c,int corner){return v_add(c->center,(Vector3){corner_signs[corner][0]*c->size.x*.5f,corner_signs[corner][1]*c->size.y*.5f,corner_signs[corner][2]*c->size.z*.5f});}
static void golem_physics_reset(void){
    goActive=false;skeletonConstraintCount=0;
    free(goCells);free(goLinks);free(goWelds);free(skeletonConstraints);
    free(goRollback.cells);memset(&goRollback,0,sizeof(goRollback));
    goCells=NULL;goLinks=NULL;goWelds=NULL;skeletonConstraints=NULL;
    goCount=goPartStart=goLinkCount=goWeldCount=goWeldCapacity=skeletonConstraintCapacity=0;
    goAllocationFailed=false;goConnectNeeded=false;goLastDt=0;goClock=goError=goLoad=0;goAction=0;goActionStart=0;
    goLastHandLaunch=-1;
    memset(goRepair,0,sizeof(goRepair));
    memset(goOverload,0,sizeof(goOverload));memset(goLost,0,sizeof(goLost));
}
static bool go_spawn_cell(int i,bool initial){
    GoCell *c=&goCells[i];Matrix transform=go_matrix(c->part);Vector3 position=Vector3Transform(c->center,transform);
    int id=addVoxel(position.x,position.y,position.z,false,true,c->color,0);if(id<0)return false;
    Voxel *v=&voxels[id];v->golemOwned=true;v->persistentConstruction=true;v->glueEligible=false;
    v->owner=STATIC_DEBRIS_OWNER; /* A manufactured body has no terrain to restore. */
    v->rest_dimensions=v_mul(c->size,goScale);v->rest_edge=fmaxf(v->rest_dimensions.x,fmaxf(v->rest_dimensions.y,v->rest_dimensions.z));
    v->rest_volume=v->rest_dimensions.x*v->rest_dimensions.y*v->rest_dimensions.z;
    v->particle_radius=.008f; /* Contact skin, not an extra half-cell around each corner. */
    for(int k=0;k<8;k++){
        Particle *p=v->particles[k];p->pos=Vector3Transform(go_corner(c,k),transform);p->prev_pos=p->predicted_pos=p->pos;
        p->radius=v->particle_radius;p->authored_shape=true;p->vel=(Vector3){0};
        if(!initial&&goLastDt>0){
            Matrix previous=c->part<GOLEM_JOINT_COUNT?goPreviousPose.joint[c->part]:c->part==GM_CHEST?goPreviousPose.chest:c->part==GM_CORE?goPreviousPose.core:goPreviousPose.eye;
            p->vel=v_mul(v_sub(p->pos,Vector3Transform(go_corner(c,k),previous)),1/goLastDt);
            p->prev_pos=v_sub(p->pos,v_mul(p->vel,goLastDt));
        }
    }
    c->id=v->identity;c->index=id;c->missing=0;collisionTopologyDirty=true;goConnectNeeded=true;return true;
}
static bool go_link_definition(int a,int b,int face){
    GoLink *next=realloc(goLinks,(size_t)(goLinkCount+1)*sizeof(*next));if(!next)return false;goLinks=next;
    GoLink *link=&goLinks[goLinkCount++];*link=(GoLink){.a=a,.b=b,.face=face};
    for(int k=0;k<4;k++){
        int corner=face_corner_indices[opposite_face[face]][k];link->corner[k]=corner;
        Vector3 relative=v_sub(go_corner(&goCells[b],corner),goCells[a].center);
        for(int j=0;j<8;j++)link->weight[k][j]=
            (.5f+corner_signs[j][0]*relative.x/goCells[a].size.x)*
            (.5f+corner_signs[j][1]*relative.y/goCells[a].size.y)*
            (.5f+corner_signs[j][2]*relative.z/goCells[a].size.z);
    }return true;
}
static void go_connect(void){
    if(!goConnectNeeded)return;goConnectNeeded=false;
    for(int n=0;n<goLinkCount;n++){
        GoLink *d=&goLinks[n];GoCell *a=&goCells[d->a],*b=&goCells[d->b];
        if(a->index<0||b->index<0)continue;
        Voxel *va=&voxels[a->index],*vb=&voxels[b->index];int face=d->face;
        bool known=false;for(int w=0;w<goWeldCount;w++)if(goWelds[w].a==a->id&&goWelds[w].b==b->id){known=true;break;}
        if(known)continue;
        if(va->glued_faces[face]||vb->glued_faces[opposite_face[face]])continue;
        if(goWeldCount==goWeldCapacity){
            int capacity=goWeldCapacity+goLinkCount;
            GoWeld *next=realloc(goWelds,(size_t)capacity*sizeof(*next));if(!next){goConnectNeeded=true;return;}goWelds=next;goWeldCapacity=capacity;
        }
        goWelds[goWeldCount++]=(GoWeld){*d,a->id,b->id,a->index,b->index};
        va->glued_faces[face]=vb->glued_faces[opposite_face[face]]=true;
        va->glued_neighbors[face]=b->index;vb->glued_neighbors[opposite_face[face]]=a->index;collisionTopologyDirty=true;
    }
}
static void go_release_part(int part){
    if(part<=0||part>=GM_COUNT)return;
    goLost[part]=goClock;goOverload[part]=0;
    for(int i=0;i<goCount;i++)if(goCells[i].part==part&&goCells[i].index>=0){
        int index=go_find(goCells[i].id,goCells[i].index);
        if(index>=0)voxels[index].persistentConstruction=false;goCells[i].id=0;goCells[i].index=-1;goCells[i].missing=goClock;
    }
    collisionTopologyDirty=true;
    if(part==GOLEM_ARM_L||part==GOLEM_FOREARM_L||part==GOLEM_ARM_R||part==GOLEM_FOREARM_R||part==GOLEM_THIGH_L||part==GOLEM_SHIN_L||part==GOLEM_THIGH_R||part==GOLEM_SHIN_R)go_release_part(part+1);
}
#include "golem_repair.inc"
static void go_launch_hand(void){
    Vector3 origin={0};int count=0;
    for(int i=0;i<goCount;i++)if(goCells[i].part==GOLEM_HAND_R&&goCells[i].index>=0){
        origin=v_add(origin,voxels[goCells[i].index].pos);count++;
    }
    if(!count)return;
    origin=v_mul(origin,1.f/count);
    Vector3 forward={-sinf(goInput.yaw*DEG2RAD),0,-cosf(goInput.yaw*DEG2RAD)};
    Vector3 target=v_add(origin,v_mul(forward,16));float nearest=FLT_MAX;
    for(int p=0;p<MAX_PLAYERS&&p<activePlayers;p++)if(players[p].respawn_timer<=0&&!is_player_bot(p)){
        float distance=Vector3DistanceSqr(origin,players[p].pos);
        if(distance<nearest){nearest=distance;target=v_add(players[p].pos,(Vector3){0,-.45f,0});}
    }
    Vector3 delta=v_sub(target,origin);
    float flight=fmaxf(.15f,sqrtf(delta.x*delta.x+delta.z*delta.z)/22.f);
    Vector3 velocity=v_mul(delta,1/flight);velocity.y+=.5f*GRAVITY*flight;
    /* Apply before GPU packing, not during a GPU substep: both backends must
     * receive the impulse once. Equal corner velocities preserve the fist. */
    for(int i=0;i<goCount;i++)if(goCells[i].part==GOLEM_HAND_R&&goCells[i].index>=0)
        set_voxel_velocity(&voxels[goCells[i].index],velocity);
    go_release_part(GOLEM_HAND_R);
}
static void golem_physics_prepare(void){
    if(!goActive)return;
    for(int i=0;i<goCount;i++){
        GoCell *c=&goCells[i];int old=c->index;c->index=go_find(c->id,c->index);
        if(c->index>=0&&!voxels[c->index].vgs_active){voxels[c->index].persistentConstruction=false;c->index=-1;}
        if(c->index<0&&old>=0){if(getenv("GOLEM_DEBUG"))printf("Golem detach: removed/inactive cell=%d part=%d time=%.3f\n",i,c->part,goClock);c->id=0;c->missing=goClock;if(goMount[c->part]==i)go_release_part(c->part);}
        if(c->index>=0){Voxel *v=&voxels[c->index];if(v->sleeping)wake_sleeping_cluster_in_place(c->index);v->sleepFrames=0;v->lifeFrames=0;v->wake_source=true;}
    }
    /* Connectivity is checked against live face links, including gun/sword cuts. */
    bool *seen=calloc((size_t)goCount,sizeof(bool));if(!seen)return;
    for(int p=1;p<GM_COUNT;p++)if(goMount[p]>=0&&goCells[goMount[p]].index>=0)seen[goMount[p]]=true;
    bool changed=true;while(changed){changed=false;for(int n=0;n<goLinkCount;n++){
        GoLink *d=&goLinks[n];GoCell *a=&goCells[d->a],*b=&goCells[d->b];if(a->index<0||b->index<0)continue;
        Voxel *v=&voxels[a->index];if(!v->glued_faces[d->face]||v->glued_neighbors[d->face]!=b->index)continue;
        if(seen[d->a]!=seen[d->b]){seen[d->a]=seen[d->b]=true;changed=true;}
    }}
    for(int i=0;i<goCount;i++)if(goCells[i].index>=0&&!seen[i]){
        if(getenv("GOLEM_DEBUG"))printf("Golem detach: broken connection cell=%d part=%d time=%.3f\n",i,goCells[i].part,goClock);
        voxels[goCells[i].index].persistentConstruction=false;goCells[i].index=-1;goCells[i].id=0;goCells[i].missing=goClock;
    }free(seen);
    go_mount_loads();
    go_repair_prepare();
    go_connect();
    float cycleStart=goActionStart+floorf((goClock-goActionStart)/6)*6;
    if(goAction==GOLEM_THROW&&goClock-cycleStart>=.7f+.6f*golem_action_seconds(GOLEM_THROW)&&goLastHandLaunch<cycleStart){
        goLastHandLaunch=cycleStart;go_launch_hand();
    }
}
static void golem_physics_tick(float dt){
    skeletonConstraintCount=0;if(!goActive)return;
    goClock+=dt;goInput.position.y=.04f;goInput.time=goClock;goInput.action=(GolemAction)goAction;goInput.actionProgress=Clamp((fmodf(goClock-goActionStart,6)-.7f)/golem_action_seconds(goInput.action),0,1);
    if(goAction==GOLEM_WALK){goInput.travel=goClock*.45f;goInput.position.x=sinf(goClock*.25f)*1.0f;goInput.yaw=sinf(goClock*.25f)*20;}
    for(int side=0;side<2;side++){
        int arm=side?GOLEM_ARM_R:GOLEM_ARM_L,shin=side?GOLEM_SHIN_R:GOLEM_SHIN_L;
        int parts[2]={arm,shin};for(int j=0;j<2;j++){
            int part=parts[j],section=j?GOLEM_LEFT_LEG+side:side;bool absent=false;
            int mount=goMount[part];absent=mount<0||goCells[mount].index<0;
            GolemDamage *d=&goInput.damage[section];if(absent){if(!d->removed){d->removed=true;d->age=0;}d->age+=dt;d->recovery=0;}else if(d->removed){d->recovery=fminf(1,d->recovery+dt/.8f);if(d->recovery>=1)*d=(GolemDamage){0};}
        }
    }
    GolemPose previousPose=goPose;goPreviousPose=goPose;goLastDt=dt;
    goPose=golem_sample(&goInput);
    float blend=1-expf(-dt/0.15f);
    for(int p=0;p<GOLEM_JOINT_COUNT;p++)goPose.joint[p]=go_smooth_matrix(previousPose.joint[p],goPose.joint[p],blend);
    goPose.chest=go_smooth_matrix(previousPose.chest,goPose.chest,blend);
    goPose.core=go_smooth_matrix(previousPose.core,goPose.core,blend);
    goPose.eye=go_smooth_matrix(previousPose.eye,goPose.eye,blend);
    /* The authored poses assume a point foot. Keep the physical cuboid soles
     * above the lab floor, including crouches and the planted slam fist. */
    go_floor_pose();
    if(!skeleton_reserve(goWeldCount*4+GM_COUNT*4+goCount*8))return;
    for(int n=0;n<goWeldCount;){
        GoWeld *w=&goWelds[n];w->ia=go_find(w->a,w->ia);w->ib=go_find(w->b,w->ib);
        if(w->ia<0||w->ib<0){goWelds[n]=goWelds[--goWeldCount];continue;}n++;
        Voxel *a=&voxels[w->ia],*b=&voxels[w->ib];int face=w->definition.face;
        if(!a->glued_faces[face]||a->glued_neighbors[face]!=w->ib)continue;
        for(int k=0;k<4;k++){
            SkeletonConstraint *c=&skeletonConstraints[skeletonConstraintCount++];*c=(SkeletonConstraint){.corner=b->particles[w->definition.corner[k]],.compliance=1e-9f,.limit=a->rest_edge*.25f};
            for(int j=0;j<8;j++){c->support[j]=a->particles[j];c->weight[j]=w->definition.weight[k][j];}
        }
    }
    /* A compliant distributed motor supplies pose acceleration to each
     * connected cell. Four segment mounts carry the attachment; the motor
     * avoids asking a single weld chain to transmit an instantaneous turn.
     * Disconnected cells immediately lose this drive and retain momentum. */
    for(int i=0;i<goCount;i++)if(goCells[i].index>=0){
        GoCell *cell=&goCells[i];int part=cell->part;Voxel *v=&voxels[cell->index];
        Matrix previous=part<GOLEM_JOINT_COUNT?previousPose.joint[part]:part==GM_CHEST?previousPose.chest:part==GM_CORE?previousPose.core:previousPose.eye;
        for(int k=0;k<8;k++)skeletonConstraints[skeletonConstraintCount++]=(SkeletonConstraint){
            .corner=v->particles[k],.target=Vector3Transform(go_corner(cell,k),go_matrix(part)),
            .previous_target=Vector3Transform(go_corner(cell,k),previous),.compliance=1e-3f,.limit=v->rest_edge*.25f};
    }
    for(int part=1;part<GM_COUNT;part++){
        int mount=goMount[part];if(mount<0||goCells[mount].index<0)continue;GoCell *cell=&goCells[mount];Voxel *v=&voxels[cell->index];
        Matrix previous=part<GOLEM_JOINT_COUNT?previousPose.joint[part]:part==GM_CHEST?previousPose.chest:part==GM_CORE?previousPose.core:previousPose.eye;
        for(int k=0;k<4;k++){
            int anchor=goAnchor[part][k];if(anchor<0||goCells[anchor].index<0)continue;
            GoCell *ac=&goCells[anchor];int corner=goAnchorCorner[part][k];
            Vector3 target=Vector3Transform(go_corner(ac,corner),go_matrix(part));
            Particle *p=voxels[ac->index].particles[corner];
            skeletonConstraints[skeletonConstraintCount++]=(SkeletonConstraint){.corner=p,.target=target,.previous_target=Vector3Transform(go_corner(ac,corner),previous),.compliance=1e-7f,.limit=v->rest_edge*.25f};
        }
    }
}
static bool golem_physics_start(void){
    golem_physics_reset();goRegen=true;goInput=golem_input();goInput.scale=goScale;goInput.position=(Vector3){0,.04f,10};goInput.side=1;goPose=golem_sample(&goInput);
    golem_body_generate(go_collect,go_finish);if(goAllocationFailed){golem_physics_reset();return false;}
    for(int p=0;p<GM_COUNT;p++)goMount[p]=-1;
    for(int i=0;i<goCount;i++){int p=goCells[i].part;int old=goMount[p];if(old<0||v_length(goCells[i].center)<v_length(goCells[old].center))goMount[p]=i;}
    /* Four widely separated material points mount each connected segment.
     * Nothing is made kinematic: these are compliant solver constraints. */
    for(int part=1;part<GM_COUNT;part++){
        Vector3 low={FLT_MAX,FLT_MAX,FLT_MAX},high={-FLT_MAX,-FLT_MAX,-FLT_MAX};
        for(int i=0;i<goCount;i++)if(goCells[i].part==part)for(int k=0;k<8;k++){
            Vector3 q=go_corner(&goCells[i],k);low=Vector3Min(low,q);high=Vector3Max(high,q);
        }
        const int tetra[4]={0,3,5,6};
        for(int k=0;k<4;k++){
            Vector3 target={tetra[k]&1?high.x:low.x,tetra[k]&2?high.y:low.y,tetra[k]&4?high.z:low.z};float best=FLT_MAX;goAnchor[part][k]=-1;
            for(int i=0;i<goCount;i++)if(goCells[i].part==part)for(int c=0;c<8;c++){
                float distance=Vector3DistanceSqr(target,go_corner(&goCells[i],c));
                if(distance<best){best=distance;goAnchor[part][k]=i;goAnchorCorner[part][k]=c;}
            }
        }
    }
    /* Each face selects its nearest overlapping neighbor within the segment. */
    for(int a=0;a<goCount;a++)for(int face=0;face<6;face++){
        int axis=face/2;float sign=face%2?-1:1;int best=-1;float score=100;
        for(int b=0;b<goCount;b++)if(a!=b&&goCells[a].part==goCells[b].part){
            Vector3 d=v_sub(goCells[b].center,goCells[a].center);float *delta=&d.x,*sa=&goCells[a].size.x,*sb=&goCells[b].size.x;
            if(delta[axis]*sign<=0)continue;
            float gap=fabsf(delta[axis])-(sa[axis]+sb[axis])*.5f;
            bool overlap=true;for(int j=0;j<3;j++)if(j!=axis&&fabsf(delta[j])>=(sa[j]+sb[j])*.5f-.01f)overlap=false;
            if(overlap&&gap>=-.12f&&gap<.04f&&fabsf(gap)<score){best=b;score=fabsf(gap);}
        }
        if(best>a&&!go_link_definition(a,best,face)){golem_physics_reset();return false;}
    }
    go_floor_pose();goPreviousPose=goPose;
    goActive=true;
    for(int i=0;i<goCount;i++)if(!go_spawn_cell(i,true)){goAllocationFailed=true;break;}
    if(goAllocationFailed){
        for(int i=voxel_count-1;i>=0;i--)if(voxels[i].golemOwned)remove_voxel_index(i);
        golem_physics_reset();return false;
    }
    go_connect();rebuild_voxel_hash();goClock=0;return true;
}
static void golem_lab_draw(void){
    if(!goActive)return;
    if(goOverlay)golem_draw_skeleton(&goPose);
    if(goRegen)for(int i=0;i<goCount;i++)if(goCells[i].index<0&&goClock-goCells[i].missing>1.25f){
        Vector3 corners[8];go_cell_points(i,corners);
        Color color=goRepair[go_repair_group(goCells[i].part)].phase==GO_REPAIR_BLOCKED?ORANGE:SKYBLUE;
        for(int k=0;k<8;k++)for(int axis=1;axis<=4;axis*=2)if(!(k&axis))DrawLine3D(corners[k],corners[k|axis],Fade(color,.45f));
    }
}

/* Only the reactor and visor emit. Use their simulated corners, never the
 * animated bone transform, so glow follows impacts and detached cells. */
static void golem_emission(void){
    if(!goActive)return;
    for(int i=0;i<voxel_count;i++){
        Voxel *v=&voxels[i];
        if(!v->golemOwned||!v->vgs_active||v->color.r!=255)continue;
        float remaining=v->golemReclaimStart>0?Clamp(1-(goClock-v->golemReclaimStart)/.25f,0,1):1;
        Vector3 center={0};for(int k=0;k<8;k++)center=v_add(center,v_mul(v->particles[k]->pos,.125f));
        for(int face=0;face<6;face++){
            Vector3 a=v->particles[face_corner_indices[face][0]]->pos;
            Vector3 b=v->particles[face_corner_indices[face][1]]->pos;
            Vector3 c=v->particles[face_corner_indices[face][3]]->pos;
            Vector3 d=v->particles[face_corner_indices[face][2]]->pos;
            a=Vector3Lerp(center,a,remaining);b=Vector3Lerp(center,b,remaining);
            c=Vector3Lerp(center,c,remaining);d=Vector3Lerp(center,d,remaining);
            DrawTriangle3D(a,b,c,v->color);DrawTriangle3D(c,b,a,v->color);
            DrawTriangle3D(a,c,d,v->color);DrawTriangle3D(d,c,a,v->color);
        }
    }
}

static void golem_collision_groups(void){
 if(!goActive)return;
 for(int i=0;i<goCount;i++)if(goCells[i].index>=0)for(int k=0;k<8;k++)voxels[goCells[i].index].particles[k]->collision_group=MAX_VOXELS+goCells[i].part;
}
