/* Repairs are committed between physics steps. Each articulated segment keeps
 * its own weld graph; a transaction may contain several dependent segments. */
static const int goParent[GM_COUNT]={-1,0,1,2,2,4,5,2,7,8,1,10,11,1,13,14,2,2,3};
static int go_repair_group(int part){
    if(part>=GOLEM_ARM_L&&part<=GOLEM_HAND_L)return GOLEM_ARM_L;
    if(part>=GOLEM_ARM_R&&part<=GOLEM_HAND_R)return GOLEM_ARM_R;
    if(part>=GOLEM_THIGH_L&&part<=GOLEM_FOOT_R)return GOLEM_THIGH_L;
    return part;
}
static void go_mount_loads(void){
    goError=goLoad=0;if(goLastDt<=0)return;
    for(int part=1;part<GM_COUNT;part++){
        if(goMount[part]<0||goCells[goMount[part]].index<0)continue;
        float mass=0,load=0,error=0;
        for(int i=0;i<goCount;i++)if(goCells[i].part==part&&goCells[i].index>=0)
            for(int k=0;k<8;k++){float im=voxels[goCells[i].index].particles[k]->inv_mass;if(im>0)mass+=1/im;}
        for(int k=0;k<4;k++){
            int cell=goAnchor[part][k];if(cell<0||goCells[cell].index<0)continue;
            Particle *p=voxels[goCells[cell].index].particles[goAnchorCorner[part][k]];
            float distance=v_length(v_sub(Vector3Transform(go_corner(&goCells[cell],goAnchorCorner[part][k]),go_matrix(part)),p->pos));
            error=fmaxf(error,distance);if(p->inv_mass>0)load+=distance/(goLastDt*goLastDt*p->inv_mass);
        }
        goError=fmaxf(goError,error);goLoad=fmaxf(goLoad,load/fmaxf(mass*GRAVITY,1));
        bool recovering=goRepair[go_repair_group(part)].phase==GO_REPAIR_RECOVER;
        goOverload[part]=!recovering&&load>8*mass*GRAVITY?goOverload[part]+2*goLastDt:0;
        if(goOverload[part]>.12f){
            if(getenv("GOLEM_DEBUG"))printf("Golem detach: overload part=%d time=%.3f error=%.4f load=%.2fx\n",part,goClock,error,load/fmaxf(mass*GRAVITY,1));
            go_release_part(part);
        }
    }
}
static void go_floor_pose(void){
    float bottom=0;
    /* Include missing geometry: regrowth must not change the floor offset. */
    for(int i=0;i<goCount;i++)for(int k=0;k<8;k++)
        bottom=fminf(bottom,Vector3Transform(go_corner(&goCells[i],k),go_matrix(goCells[i].part)).y-.012f);
    if(bottom<0){
        for(int p=0;p<GOLEM_JOINT_COUNT;p++)goPose.joint[p].m13-=bottom;
        goPose.chest.m13-=bottom;goPose.core.m13-=bottom;goPose.eye.m13-=bottom;
    }
}
static void go_cell_points(int cell,Vector3 points[8]){
    for(int k=0;k<8;k++)points[k]=Vector3Transform(go_corner(&goCells[cell],k),go_matrix(goCells[cell].part));
}
static bool go_separates(Vector3 axis,const Vector3 a[8],const Vector3 b[8]){
    float length=v_length(axis);if(length<1e-6f)return false;axis=v_mul(axis,1/length);
    float alo=FLT_MAX,ahi=-FLT_MAX,blo=FLT_MAX,bhi=-FLT_MAX;
    for(int k=0;k<8;k++){float x=v_dot(a[k],axis),y=v_dot(b[k],axis);alo=fminf(alo,x);ahi=fmaxf(ahi,x);blo=fminf(blo,y);bhi=fmaxf(bhi,y);}
    return ahi<=blo+.002f || bhi<=alo+.002f;
}
static bool go_boxes_overlap(const Vector3 a[8],const Vector3 b[8]){
    BoundingBox ba={a[0],a[0]},bb={b[0],b[0]};
    for(int k=1;k<8;k++){ba.min=Vector3Min(ba.min,a[k]);ba.max=Vector3Max(ba.max,a[k]);bb.min=Vector3Min(bb.min,b[k]);bb.max=Vector3Max(bb.max,b[k]);}
    if(!CheckCollisionBoxes(ba,bb))return false;
    Vector3 ea[3]={v_sub(a[1],a[0]),v_sub(a[2],a[0]),v_sub(a[4],a[0])};
    Vector3 eb[3]={v_sub(b[1],b[0]),v_sub(b[2],b[0]),v_sub(b[4],b[0])};
    for(int j=0;j<3;j++){
        if(go_separates(Vector3CrossProduct(ea[j],ea[(j+1)%3]),a,b) ||
           go_separates(Vector3CrossProduct(eb[j],eb[(j+1)%3]),a,b))return false;
        for(int k=0;k<3;k++)if(go_separates(Vector3CrossProduct(ea[j],eb[k]),a,b))return false;
    }
    return true;
}
static void go_box_points(BoundingBox box,Vector3 points[8]){
    for(int k=0;k<8;k++)points[k]=(Vector3){k&1?box.max.x:box.min.x,k&2?box.max.y:box.min.y,k&4?box.max.z:box.min.z};
}
static bool go_held(int index){
    if(tetherTag[index]>0||tetherHoldClusterTag[index]>0)return true;
    for(int p=0;p<activePlayers;p++)if(players[p].tetherHolding&&players[p].tetherVoxel==index)return true;
    return false;
}
static bool go_reclaimable(int index){
    Voxel *v=&voxels[index];return v->golemOwned&&!v->persistentConstruction&&!go_held(index);
}
static void go_repair_prepare(void){
    /* Index storage is scratch only; persistent references use identities. */
    int *pending=malloc((size_t)goCount*sizeof(int));
    GoCell *backup=malloc((size_t)goCount*sizeof(GoCell));
    uint64_t *reclaim=malloc((size_t)(voxel_count+1)*sizeof(uint64_t));
    if(!pending||!backup||!reclaim){free(pending);free(backup);free(reclaim);return;}
    for(int n=0;n<voxel_count;n++)if(voxels[n].golemReclaimStart>0&&!go_reclaimable(n))voxels[n].golemReclaimStart=0;
    for(int group=1;group<GM_COUNT;group++){
        int count=0;float earliest=goClock;
        for(int i=0;i<goCount;i++)if(goCells[i].index<0&&go_repair_group(goCells[i].part)==group){pending[count++]=i;earliest=fminf(earliest,goCells[i].missing);}
        GoRepair *repair=&goRepair[group];
        if(!count){if(repair->phase!=GO_REPAIR_RECOVER||goClock-repair->since>=.8f)repair->phase=GO_REPAIR_IDLE;continue;}
        if(repair->phase==GO_REPAIR_IDLE||repair->phase==GO_REPAIR_RECOVER)*repair=(GoRepair){GO_REPAIR_WAIT,earliest};
        if(!goRegen||goClock-repair->since<1.25f)continue;
        repair->phase=GO_REPAIR_BLOCKED;
        bool blocked=false;int reclaimCount=0;
        for(int k=0;k<count&&!blocked;k++){
            int cell=pending[k],parent=goParent[goCells[cell].part];
            if(parent>0&&goCells[goMount[parent]].index<0&&go_repair_group(parent)!=group){blocked=true;break;}
            Vector3 a[8];go_cell_points(cell,a);
            for(int p=0;p<activePlayers&&!blocked;p++)if(players[p].respawn_timer<=0){
                Vector3 b[8];go_box_points(combatant_bounds_at(p,players[p].pos,PLAYER_SIZE*.5f),b);
                if(go_boxes_overlap(a,b))blocked=true;
            }
            for(int n=0;n<voxel_count&&!blocked;n++){
                Voxel *v=&voxels[n];if(v->golemOwned&&v->persistentConstruction)continue;
                Vector3 b[8];
                for(int c=0;c<8;c++)b[c]=v->particles[c]?v->particles[c]->pos:v_add(v->pos,(Vector3){corner_signs[c][0]*v->rest_edge*.5f,corner_signs[c][1]*v->rest_edge*.5f,corner_signs[c][2]*v->rest_edge*.5f});
                if(!go_boxes_overlap(a,b))continue;
                if(!go_reclaimable(n)){blocked=true;break;}
                bool known=false;for(int j=0;j<reclaimCount;j++)if(reclaim[j]==v->identity)known=true;
                if(!known)reclaim[reclaimCount++]=v->identity;
            }
        }
        if(blocked)continue;
        /* Reserve all fallible storage before changing the world or reclaiming
         * debris. No joint can be published with only part of its cells. */
        if(voxel_count+count>MAX_VOXELS || count*8>free_particle_count+MAX_PARTICLES-particle_pool_count)continue;
        GoWeld *reserved=realloc(goWelds,(size_t)(goWeldCount+goLinkCount)*sizeof(*reserved));
        if(!reserved)continue;goWelds=reserved;goWeldCapacity=goWeldCount+goLinkCount;
        if(!skeleton_reserve((goWeldCount+goLinkCount)*4+GM_COUNT*4+goCount*8))continue;
        bool dissolving=false;
        for(int j=0;j<reclaimCount;j++){
            int n=go_find(reclaim[j],-1);if(n<0)continue;
            if(voxels[n].golemReclaimStart<=0)voxels[n].golemReclaimStart=fmaxf(goClock,.0001f);
            if(goClock-voxels[n].golemReclaimStart<.25f)dissolving=true;
        }
        repair->phase=GO_REPAIR_BUILD;if(dissolving)continue;
        int made=0;
        for(int k=0;k<count;k++){
            int i=pending[k];backup[k]=goCells[i];
#ifdef FIREFIGHT_TESTING
            if(goRepairFailAfter>=0&&made>=goRepairFailAfter)break;
#endif
            /* Obstruction was checked for the complete assembly above. */
            if(!go_spawn_cell(i,false))break;
            made++;
        }
        if(made!=count){
            while(made>0){--made;int i=pending[made];remove_voxel_index(go_find(goCells[i].id,goCells[i].index));goCells[i]=backup[made];}
            repair->phase=GO_REPAIR_BLOCKED;continue;
        }
        /* Remove only debris still overlapping this transaction. Compaction
         * is resolved before publishing connections. */
        for(int j=0;j<reclaimCount;j++){int n=go_find(reclaim[j],-1);if(n>=0)remove_voxel_index(n);}
        for(int i=0;i<goCount;i++)goCells[i].index=go_find(goCells[i].id,goCells[i].index);
        go_connect();collisionTopologyDirty=true;
        for(int k=0;k<count;k++)goOverload[goCells[pending[k]].part]=0;
        *repair=(GoRepair){GO_REPAIR_RECOVER,goClock};
    }
    /* A moving/held obstruction or a disabled repair cancels its cosmetic
     * dissolve. Reclamation never leaves a surviving voxel invisible. */
    for(int n=0;n<voxel_count;n++)if(voxels[n].golemReclaimStart>0){
        bool needed=false;
        if(goRegen&&go_reclaimable(n)){
            Vector3 b[8];for(int k=0;k<8;k++)b[k]=voxels[n].particles[k]->pos;
            for(int i=0;i<goCount&&!needed;i++)if(goCells[i].index<0&&goRepair[go_repair_group(goCells[i].part)].phase==GO_REPAIR_BUILD){
                Vector3 a[8];go_cell_points(i,a);needed=go_boxes_overlap(a,b);
            }
        }
        if(!needed)voxels[n].golemReclaimStart=0;
    }
    free(pending);free(backup);free(reclaim);
}
static const char *go_repair_status(void){
    static const char *names[]={"ready","waiting","blocked","rebuilding","recovering"};
    int phase=0;for(int p=1;p<GM_COUNT;p++)if(goRepair[p].phase>phase)phase=goRepair[p].phase;
    return names[phase];
}
