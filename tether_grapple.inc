/* Player-side response to a default voxel tether. Never run in voxel substeps. */
#define TETHER_GRAPPLE_STOP_DISTANCE 1.5f
#define TETHER_GRAPPLE_SPRING 20.0f
#define TETHER_GRAPPLE_DAMPING 4.0f
#define TETHER_GRAPPLE_ACCELERATION 30.0f
#define TETHER_GRAPPLE_MAX_APPROACH 12.0f

static bool tether_grapple_anchor(int i,Vector3 *target,Vector3 *velocity){
    if(i<0 || i>=MAX_PLAYERS || i>=activePlayers || is_player_bot(i) ||
       gameState!=GAME_STATE_PLAYING || creativeModeActive || ff_input_blocked(i) ||
       combatantLife[i]!=LIFE_ALIVE || players[i].respawn_timer>0 || (enemyLab&&enemyLabPaused))return false;
    Player *p=&players[i];
    if(netTransport.role==NET_ROLE_CLIENT){
        if(!p->grappleInputHeld || !p->netGrappleAnchorTick || !p->netTetherVisualActive ||
           p->netGoldTetherVisualActive || p->netDrainVisualActive)return false;
        *target=world_image(p->netTetherVisualTarget,p->pos);
        *velocity=p->netGrappleAnchorVelocity;
    }else{
        if(!p->tetherHolding || p->goldTetherHolding || p->drainTarget || !resolve_tether_voxel(p))return false;
        Voxel *v=&voxels[p->tetherVoxel];
        if(v->isBullet || v->type!=0 || voxel_is_fluid(v) || (v->simulate&&!v->vgs_active))return false;
        if(!player_tether_visual_target(i,target))return false;
        *velocity=(Vector3){0};
        if(v->simulate && p->tetherFace>=0 && p->tetherFace<6){
            for(int k=0;k<4;k++){
                Particle *corner=v->particles[face_corner_indices[p->tetherFace][k]];
                if(corner)*velocity=v_add(*velocity,v_mul(corner->vel,.25f));
            }
        }
    }
    return v_isfinite(*target) && v_isfinite(*velocity);
}

static bool tether_grapple_active(int i){
    Vector3 target,velocity;
    return tether_grapple_anchor(i,&target,&velocity) &&
        v_length(world_delta(target,players[i].pos))>TETHER_GRAPPLE_STOP_DISTANCE;
}

static void tether_grapple_pull(int i,float dt){
    Vector3 target,anchorVelocity;
    if(dt<=0 || !tether_grapple_anchor(i,&target,&anchorVelocity))return;
    Player *p=&players[i];Vector3 delta=world_delta(target,p->pos);float distance=v_length(delta);
    if(distance<=TETHER_GRAPPLE_STOP_DISTANCE)return;
    Vector3 direction=v_mul(delta,1/distance);
    float approach=v_dot(v_sub(p->vel,anchorVelocity),direction);
    float acceleration=clampf(TETHER_GRAPPLE_SPRING*(distance-TETHER_GRAPPLE_STOP_DISTANCE)-
        TETHER_GRAPPLE_DAMPING*approach,0,TETHER_GRAPPLE_ACCELERATION);
    float impulse=fminf(acceleration*dt,fmaxf(0,TETHER_GRAPPLE_MAX_APPROACH-approach));
    p->vel=v_add(p->vel,v_mul(direction,impulse));
    if(direction.y*impulse>0 && p->vel.y>0)p->onGround=false;
}

/* Existing target snapshots supply anchor velocity without a new wire field. */
static void tether_grapple_receive(int i,unsigned flags,Vector3 target,uint32_t tick){
    Player *p=&players[i];
    bool active=(flags&NET_PLAYER_VISUAL_TETHER) && !(flags&(NET_PLAYER_VISUAL_GOLD_TETHER|NET_PLAYER_VISUAL_DRAIN));
    Vector3 velocity={0};
    if(active && p->netTetherVisualActive && !p->netGoldTetherVisualActive && !p->netDrainVisualActive &&
       p->netGrappleAnchorTick && tick>p->netGrappleAnchorTick){
        float dt=(tick-p->netGrappleAnchorTick)/(float)FPS_NET_TICK_RATE;
        if(dt<=.25f)velocity=v_mul(world_delta(target,p->netTetherVisualTarget),1/dt);
        // A replacement anchor or correction must not manufacture launch speed.
        if(!v_isfinite(velocity) || v_length(velocity)>30)velocity=(Vector3){0};
    }
    p->netGrappleAnchorVelocity=velocity;p->netGrappleAnchorTick=active?tick:0;
    p->netTetherVisualActive=(flags&NET_PLAYER_VISUAL_TETHER)!=0;
    p->netGoldTetherVisualActive=(flags&NET_PLAYER_VISUAL_GOLD_TETHER)!=0;
    p->netDrainVisualActive=(flags&NET_PLAYER_VISUAL_DRAIN)!=0;
    p->netTetherVisualTarget=target;
}
