static void foundry_traps_reset(void) {
    for(int n=0;n<foundryTrapCount;n++)free(foundryTraps[n].original);
    memset(foundryTraps,0,sizeof(foundryTraps));foundryTrapCount=0;
    foundryFragileReferenceCells=foundryFragileCells=foundryFragileReferenceCount=0;
}
static bool foundry_trap_hit(int id,int actor,bool payload) {
    if(id<0 || id>=voxel_count || voxels[id].simulate)return false;
    Voxel *v=&voxels[id];
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];if(t->active || t->pending)continue;
        if(foundry_trap_trigger(t,v->gx,v->gy,v->gz) || (payload && foundry_trap_contains(t,v->gx,v->gy,v->gz))) {
            t->pending=true;t->owner=actor;return true;
        }
    }
    return false;
}
static bool foundry_trap_complete(const FoundryTrap *t) {
    for(int k=0;k<t->originalCount;k++){const Voxel *v=&t->original[k];int id=table_get_static_only(v->gx,v->gy,v->gz);if(id<0 || voxels[id].owner!=-1)return false;}
    return true;
}
static void foundry_trap_release(int n) {
    FoundryTrap *t=&foundryTraps[n];
    int total=t->kind==FOUNDRY_BAY?t->originalCount:(t->x1-t->x0+1)*(t->y1-t->y0+1)*(t->z1-t->z0+1);
    /* Capacity preflight keeps exhausted pools from starting a release. */
    if(total>FOUNDRY_TRAP_PAYLOAD || voxel_count+total>=MAX_VOXELS || particle_pool_count-free_particle_count+8*total>=MAX_PARTICLES || recycleQueueCount+t->originalCount>MAX_VOXELS)return;
    t->active=true;t->pending=false;t->victims=0;t->activation++;
    static UnitVoxelBuffer buffer; buffer.count=0;
    for(int z=t->z0;z<=t->z1;z++)for(int y=t->fourPosts?t->ry+1:t->y0;y<=t->y1;y++)for(int x=t->x0;x<=t->x1;x++) {
        if(!foundry_trap_contains(t,x,y,z))continue;
        int id=table_get_static_only(x,y,z);
        if(id>=0)unit_voxel_buffer_push(&buffer,x,y,z,voxels[id].color,0,false,id,0,t->owner);
    }
    for(int k=0;k<t->originalCount;k++){const Voxel *v=&t->original[k];if(!foundry_trap_trigger(t,v->gx,v->gy,v->gz) && !(t->fourPosts && v->gy>=t->ry && v->gy<t->y0))continue;
        int release=table_get_static_only(v->gx,v->gy,v->gz);
        if(release>=0){recycle_queue_push(&voxels[release]);remove_voxel_index(release);}
    }
    // A blast may have activated orange supports before the linked release.
    if(t->fourPosts || t->kind==FOUNDRY_BAY)for(int id=voxel_count-1;id>=0;id--){
        Voxel *v=&voxels[id];
        if(v->simulate && !v->isBullet && v->type==0 &&
           ((t->kind==FOUNDRY_BAY && foundry_trap_trigger(t,v->orig_min_gx,v->orig_min_gy,v->orig_min_gz)) ||
           (t->fourPosts && v->orig_min_gy>=t->ry && v->orig_min_gy<t->y0 &&
           (v->orig_min_gx==t->x0 || v->orig_min_gx==t->x1) &&
           (v->orig_min_gz==t->z0 || v->orig_min_gz==t->z1))))remove_voxel_index(id);
    }
    // Removing the switch may compact the voxel array; refresh seed indices.
    for(int k=0;k<buffer.count;k++)buffer.voxels[k].voxelIndex=table_get_static_only(buffer.voxels[k].gx,buffer.voxels[k].gy,buffer.voxels[k].gz);
    remove_buffered_static_voxels(&buffer);
    for(int k=0;k<buffer.count;k++) {
        UnitVoxelSeed *s=&buffer.voxels[k];
        int id=addVoxel((s->gx+.5f)*VOXEL_SIZE,(s->gy+.5f)*VOXEL_SIZE,(s->gz+.5f)*VOXEL_SIZE,false,true,s->color,0);
        if(id<0)continue;
        Voxel *v=&voxels[id];v->foundryTrapId=n+1;v->owner=v->activator=t->owner;
        v->restorationQueued=s->restorationQueued;v->restorationDeadlineFrame=s->restorationDeadlineFrame;
        glue_neighbor_faces_for_voxel(id);
    }
    // An explosion can already have activated part of the payload. Adopt those
    // survivors by their original cells, without resurrecting destroyed cells.
    int members[FOUNDRY_TRAP_PAYLOAD],count=0;
    for(int id=0;id<voxel_count;id++) {
        Voxel *v=&voxels[id];
        if(!v->simulate || v->isBullet || v->type || !foundry_trap_contains(t,v->orig_min_gx,v->orig_min_gy,v->orig_min_gz))continue;
        if(count>=FOUNDRY_TRAP_PAYLOAD)break;
        v->foundryTrapId=n+1;v->owner=v->activator=t->owner;members[count++]=id;
        for(int f=0;f<6;f++) {
            int other=v->glued_neighbors[f];
            if(v->glued_faces[f] && other>=0 && other<voxel_count && !foundry_trap_contains(t,voxels[other].orig_min_gx,voxels[other].orig_min_gy,voxels[other].orig_min_gz))break_face_link(v,f);
        }
    }
    // Destructive mining/blasts may remove cells before this release. Preserve
    // the original assembly in the normal occupied-space regeneration queue.
    for(int k=0;k<t->originalCount;k++) {
        Voxel *original=&t->original[k];bool queued=false;
        for(int q=0;q<recycleQueueCount;q++) {
            Voxel *snap=&recycleQueue[(recycleQueueHead+q)%MAX_VOXELS];
            if(snap->orig_min_gx==original->gx && snap->orig_min_gy==original->gy && snap->orig_min_gz==original->gz){queued=true;break;}
        }
        if(!queued)recycle_queue_push(original);
    }
    for(int k=0;k<count;k++) {
        Voxel *v=&voxels[members[k]];v->restorationQueued=true;
        v->restorationDeadlineFrame=recycleFrameCounter+RECYCLE_STATIC_RESTORE_DELAY;
    }
    detach_gold_tether_static_boundary(members,count);
    /* Removed crack cells separate the actual irregular physical fragments. */
    int stamp=++particle_sync_stamp;
    for(int k=0;k<count;k++) {
        Voxel *v=&voxels[members[k]];v->sleeping=false;v->sleepFrames=0;v->wake_source=true;v->wake_timer=COARSENING_WAKE_FRAMES;
        for(int c=0;c<8;c++) {
            Particle *p=v->particles[c];if(!p || p->sync_stamp==stamp)continue;p->sync_stamp=stamp;
            if(t->kind==FOUNDRY_BAY && v->orig_min_gy<t->y1){
                float push=.8f+1.2f*clampf((p->pos.y-2)/3,0,1);
                p->vel.z+=v->orig_min_gz<(t->z0+t->z1)/2?push:-push;
            }
            if(t->kind==FOUNDRY_TIP)p->vel.z+=1.2f+1.8f*clampf(p->pos.y-t->y0*VOXEL_SIZE,0,3);
        }
    }
    refresh_activation_surfaces(&buffer);buffer.count=0;
    collisionTopologyDirty=true;dynamic_voxels_dirty=true;rebuild_voxel_hash();rebuild_glue_constraints();
    play_sfx_at(SFX_GLUE_BREAK,(Vector3){(t->rx+.5f)*VOXEL_SIZE,(t->ry+.5f)*VOXEL_SIZE,(t->rz+.5f)*VOXEL_SIZE},t->owner,2,.65f);
}
static void foundry_traps_blast(Vector3 pos,int owner) {
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];if(t->active)continue;
        Vector3 release={(t->rx+.5f)*VOXEL_SIZE,(t->ry+.5f)*VOXEL_SIZE,(t->rz+.5f)*VOXEL_SIZE};
        Vector3 at=world_image(pos,release);
        Vector3 nearest={clampf(at.x,t->x0*VOXEL_SIZE,(t->x1+1)*VOXEL_SIZE),clampf(at.y,t->y0*VOXEL_SIZE,(t->y1+1)*VOXEL_SIZE),clampf(at.z,t->z0*VOXEL_SIZE,(t->z1+1)*VOXEL_SIZE)};
        if(t->fourPosts)for(int x=t->x0;x<=t->x1;x+=t->x1-t->x0)for(int z=t->z0;z<=t->z1;z+=t->z1-t->z0)
            if(v_length(world_delta((Vector3){(x+.5f)*VOXEL_SIZE,(t->ry+.5f)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE},pos))<4){t->pending=true;t->owner=owner;}
        if(v_length(world_delta(release,pos))<4 || v_length(v_sub(at,nearest))<4){t->pending=true;t->owner=owner;}
    }
}
static void foundry_traps_prepare(void) {
    if(!foundryTrapCount)return;
    bool live[FOUNDRY_TRAPS]={0};
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId && voxels[i].foundryTrapId<=foundryTrapCount)live[voxels[i].foundryTrapId-1]=true;
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];bool complete=foundry_trap_complete(t);
        if(t->active) {
            if(!live[n] && complete){t->active=false;t->owner=-1;t->victims=0;}
        } else if(t->pending || !complete)foundry_trap_release(n);
    }
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId) {
        voxels[i].foundryTrapPrevious=voxels[i].pos;
        voxels[i].foundryTrapPreviousSpeed=v_length(voxels[i].vel);
    }
}
static void foundry_traps_hits(void) {
    if(!foundryTrapCount)return;
    bool any=false;for(int n=0;n<foundryTrapCount;n++)if(foundryTraps[n].active){any=true;break;}
    if(!any)return;
    static int sizes[MAX_VOXELS];memset(sizes,0,(size_t)voxel_count*sizeof(int));
    rebuild_glue_topology_neighbors();build_glue_cluster_ids(collisionVoxelClusterIds);
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId && voxels[i].vgs_active) {
        int group=collisionVoxelClusterIds[i];if(group>=0)sizes[group]++;
    }
    for(int i=0;i<voxel_count;i++) {
        Voxel *v=&voxels[i];if(!v->foundryTrapId || v->foundryTrapId>foundryTrapCount)continue;
        int group=collisionVoxelClusterIds[i];
        v->foundryTrapHeavy=v->vgs_active && group>=0 && sizes[group]>=8;
        if(!v->foundryTrapHeavy || v->sleeping || fmaxf(v_length(v->vel),v->foundryTrapPreviousSpeed)<2.5f)continue;
        hologram_projectile(i,world_image(v->foundryTrapPrevious,v->pos),v->pos,100);
        FoundryTrap *t=&foundryTraps[v->foundryTrapId-1];
        Vector3 start=world_image(v->foundryTrapPrevious,v->pos),path=v_sub(v->pos,start);float length=v_length(path);
        for(int j=0;j<activePlayers;j++) {
            if(tetherHoldClusterTag[i]==j+1 || tetherTag[i]==j+1)continue;
            if(t->victims&(1u<<j) || players[j].respawn_timer>0 || players[j].invuln_timer>0 || !enemy_core_open(j))continue;
            BoundingBox b=expand_bounds(combatant_bounds_at(j,world_image(players[j].pos,start),PLAYER_SIZE*.5f),v->rest_edge*.5f+voxel_particle_radius(v)+.03f);
            bool hit=(start.x>=b.min.x && start.x<=b.max.x && start.y>=b.min.y && start.y<=b.max.y && start.z>=b.min.z && start.z<=b.max.z);
            if(v->tetherReleaseOverlap==j+1){if(hit)continue;v->tetherReleaseOverlap=0;}
            if(!hit && length>.00001f){RayCollision c=GetRayCollisionBox((Ray){start,v_mul(path,1/length)},b);hit=c.hit && c.distance<=length;}
            if(!hit)continue;
            t->victims|=1u<<j;
            bool previous=explosionDamageContext;explosionDamageContext=true;
            if(players[j].enemyType==ENEMY_TYPE_FOUNDRY)apply_matter_damage(j,t->owner,100);
            else kill_player(j,t->owner,false,t->owner>=0);
            explosionDamageContext=previous;smushBannerTimer=1;play_sfx_at(SFX_SMUSH,players[j].pos,j,1,.8f);
        }
    }
}
static void foundry_traps_draw(void) {
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];Vector3 center={(t->x0+t->x1+1)*.25f,.015f,(t->z0+t->z1+1)*.25f};
        float width=(t->x1-t->x0+1)*.5f,depth=(t->z1-t->z0+1)*.5f;
        if(t->kind==FOUNDRY_TIP){depth+=3;center.z+=1.5f;}
        if(!periodic_visible(center,fmaxf(width,depth)))continue;
        Color c=t->active?Fade(ORANGE,.3f):ORANGE;
        DrawCubeWires(center,width,.02f,depth,c);
        for(float x=-width*.5f;x<width*.5f;x+=.75f)DrawLine3D(v_add(center,(Vector3){x,0,-depth*.5f}),v_add(center,(Vector3){x+.35f,0,-depth*.5f+.35f}),c);
        if(t->active)continue;
        Vector3 release={(t->rx+.5f)*.5f,(t->ry+.5f)*.5f,(t->rz+.5f)*.5f};
        if(t->kind==FOUNDRY_BAY){
            for(int k=0;k<t->originalCount;k++){const Voxel *v=&t->original[k];
                if(v->gx==t->rx && v->gy==t->ry && v->gz==t->rz)DrawCubeWires(v->pos,.6f,.6f,.6f,YELLOW);
            }
        }else if(t->fourPosts){
            for(int x=t->x0;x<=t->x1;x+=t->x1-t->x0)for(int z=t->z0;z<=t->z1;z+=t->z1-t->z0)
                DrawCubeWires((Vector3){(x+.5f)*VOXEL_SIZE,(t->ry+.5f)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE},.6f,.6f,.6f,YELLOW);
        }else DrawCubeWires(release,.6f,.6f,.6f,YELLOW);
        if(t->kind==FOUNDRY_SWITCH){
            DrawCubeWires(release,.75f,.75f,.75f,ORANGE);
            for(int x=t->x0;x<=t->x1;x+=t->x1-t->x0)for(int z=t->z0;z<=t->z1;z+=t->z1-t->z0)
                DrawLine3D((Vector3){(x+.5f)*VOXEL_SIZE,(t->y1+1)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE},(Vector3){(x+.5f)*VOXEL_SIZE,11*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE},ORANGE);
            DrawLine3D(release,(Vector3){center.x,(t->y0+.5f)*VOXEL_SIZE,center.z},ORANGE);
        }
    }
}
