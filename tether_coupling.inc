static bool weak_tether_enabled(int i){
    if(i<0||i>=MAX_PLAYERS||i>=activePlayers||is_player_bot(i))return false;
    Player *p=&players[i];
    return netTransport.role!=NET_ROLE_CLIENT&&gameState==GAME_STATE_PLAYING&&!creativeModeActive&&
        !ff_input_blocked(i)&&combatantLife[i]==LIFE_ALIVE&&p->respawn_timer<=0&&
        !(enemyLab&&enemyLabPaused)&&p->tetherHolding&&!p->goldTetherHolding&&!p->drainTarget;
}
static bool weak_tether_any(void){for(int i=0;i<MAX_PLAYERS;i++)if(weak_tether_enabled(i))return true;return false;}
static bool weak_tether_sample(int i,WeakTetherCoupling *c){
    c->active=false;c->force=c->impulse=(Vector3){0};c->count=0;
    if(!weak_tether_enabled(i)||!resolve_tether_voxel(&players[i]))return false;
    Player *p=&players[i];Voxel *v=&voxels[p->tetherVoxel];
    if(v->isBullet||v->type!=0||voxel_is_fluid(v)||(v->simulate&&!v->vgs_active))return false;
    if(v->sleeping){wake_sleeping_cluster_in_place(p->tetherVoxel);c->identity=0;}
    if(c->identity!=v->identity||c->generation!=collisionMetadataGeneration||c->voxelCount!=voxel_count||collisionTopologyDirty){
        int count=build_glue_cluster_indices(p->tetherVoxel,glueClusterIndices);
        if(count<=0){glueClusterIndices[0]=p->tetherVoxel;count=1;}
        c->mass=0;c->anchored=!v->simulate||v->fixed;int stamp=++particle_sync_stamp;
        for(int n=0;n<count;n++){
            Voxel *member=&voxels[glueClusterIndices[n]];c->anchored|=member->fixed||!member->simulate;
            for(int k=0;k<8;k++){Particle *q=member->particles[k];if(!q||q->sync_stamp==stamp)continue;q->sync_stamp=stamp;
                // A welded corner represents material from every owning cell.
                // Count it once, with that ownership weight; welding must not
                // make a sixteen-cube group lighter than sixteen loose cubes.
                if(q->inv_mass>0)c->mass+=fmaxf(1,q->glue_count)/q->inv_mass;else c->anchored=true;
            }
        }
        c->identity=v->identity;c->generation=collisionMetadataGeneration;c->voxelCount=voxel_count;
    }
    c->face=p->tetherFace;c->target=c->velocity=(Vector3){0};
    if(c->face<0||c->face>=6){
        Vector3 toward=world_delta(player_tether_hold_position(p),v->pos);float best=-FLT_MAX;
        for(int f=0;f<6;f++){float dot=toward.x*face_offsets[f][0]+toward.y*face_offsets[f][1]+toward.z*face_offsets[f][2];if(dot>best){best=dot;c->face=f;}}
    }
    for(int k=0;k<4;k++){
        Particle *q=v->particles[face_corner_indices[c->face][k]];if(!q)continue;
        c->target=v_add(c->target,v_mul(q->pos,.25f));c->velocity=v_add(c->velocity,v_mul(q->vel,.25f));
    }
    c->target=world_image(c->target,p->pos);
    for(int k=0;k<4;k++){
        Particle *q=v->particles[face_corner_indices[c->face][k]];if(!q||q->inv_mass<=0)continue;
        bool duplicate=false;for(int j=0;j<c->count;j++)duplicate|=c->corner[j]==q;
        if(!duplicate)c->corner[c->count++]=q;
    }
    c->active=v_isfinite(c->target)&&v_isfinite(c->velocity);
    return c->active;
}
/* Impulse on the player; the object receives its negative. Damping is relative,
 * so translating both endpoints together never changes the spring response. */
static Vector3 weak_tether_impulse(Vector3 error,Vector3 relative,float inverseMass,float dt){
    if(dt<=0||!v_isfinite(error)||!v_isfinite(relative)||!isfinite(inverseMass)||inverseMass<0)return (Vector3){0};
    float mobility=1/TETHER_GRAPPLE_PLAYER_MASS+inverseMass,mu=1/mobility;
    Vector3 force=v_mul(v_sub(v_mul(error,100),v_mul(relative,20)),mu);
    float cap=TETHER_GRAPPLE_PLAYER_MASS*TETHER_GRAPPLE_ACCELERATION;
    if(inverseMass>0)cap=fminf(cap,TETHER_GRAPPLE_VOXEL_ACCELERATION/inverseMass);
    float length=v_length(force);if(length>cap)force=v_mul(force,cap/length);
    Vector3 impulse=v_mul(force,dt);
    float distance=v_length(error);
    if(distance>1e-5f){
        Vector3 direction=v_mul(error,1/distance);float closing=v_dot(relative,direction),added=v_dot(impulse,direction)*mobility;
        if(added>0){float allowed=fmaxf(0,TETHER_GRAPPLE_MAX_APPROACH-closing);if(added>allowed)impulse=v_mul(impulse,allowed/added);}
    }
    return impulse;
}
static void weak_tether_apply_player(int i,Vector3 impulse){
    Player *p=&players[i];p->vel=v_add(p->vel,v_mul(impulse,1/TETHER_GRAPPLE_PLAYER_MASS));
    if(impulse.y>0&&p->vel.y>0)p->onGround=false;
    if(!p->onGround&&v_length(impulse)>1e-6f)p->grappleAirborne=true;
}
static void weak_tether_prepare(float dt){
    if(dt<=0||!isfinite(dt)){for(int i=0;i<MAX_PLAYERS;i++){weakTethers[i].active=false;weakTethers[i].force=weakTethers[i].impulse=(Vector3){0};}return;}
    for(int i=0;i<MAX_PLAYERS;i++){
        WeakTetherCoupling *c=&weakTethers[i];if(!weak_tether_sample(i,c))continue;
        float inverseMass=c->anchored?0:1/fmaxf(c->mass,TETHER_REFERENCE_MASS);
        Vector3 error=world_delta(c->target,player_tether_hold_position(&players[i]));
        c->impulse=weak_tether_impulse(error,v_sub(players[i].vel,c->velocity),inverseMass,dt);
        if(!c->anchored&&!c->count){c->impulse=(Vector3){0};continue;}
        c->force=v_mul(c->impulse,-1/dt);
        weak_tether_apply_player(i,c->impulse);
    }
}
static void weak_tether_integrate(float dt){
    for(int i=0;i<MAX_PLAYERS;i++){
        WeakTetherCoupling *c=&weakTethers[i];if(!c->active||!weak_tether_enabled(i)||c->count<=0)continue;
        for(int k=0;k<c->count;k++){Particle *q=c->corner[k];if(!q||q->inv_mass<=0)continue;
            Vector3 dv=v_mul(c->force,dt*q->inv_mass/(c->count*fmaxf(1,q->glue_count)));
            q->vel=v_add(q->vel,dv);q->predicted_pos=v_add(q->predicted_pos,v_mul(dv,dt));
        }
    }
}
