/* Presentation/prediction adapter. Authoritative recoil runs in fixed physics. */
static bool tether_grapple_anchor(int i,Vector3 *target,Vector3 *velocity){
    if(i<0||i>=MAX_PLAYERS||i>=activePlayers||is_player_bot(i)||gameState!=GAME_STATE_PLAYING||creativeModeActive||
       ff_input_blocked(i)||combatantLife[i]!=LIFE_ALIVE||players[i].respawn_timer>0||(enemyLab&&enemyLabPaused))return false;
    Player *p=&players[i];
    if(netTransport.role==NET_ROLE_CLIENT){
        if(!p->grappleInputHeld||!p->netGrappleAnchorTick||!p->netTetherVisualActive||p->netGoldTetherVisualActive||p->netDrainVisualActive)return false;
        *target=world_image(p->netTetherVisualTarget,p->pos);*velocity=p->netGrappleAnchorVelocity;
    }else{
        WeakTetherCoupling c=weakTethers[i];if(!weak_tether_sample(i,&c))return false;
        *target=c.target;*velocity=c.velocity;
    }
    return v_isfinite(*target)&&v_isfinite(*velocity);
}
static bool tether_grapple_active(int i){Vector3 target,velocity;return tether_grapple_anchor(i,&target,&velocity);}
static bool tether_grapple_airborne(int i){
    return i>=0&&i<MAX_PLAYERS&&!players[i].onGround&&(players[i].grappleAirborne||tether_grapple_active(i));
}
static void player_walk_accelerate(Player *p,Vector3 direction,float dt){
    float previousSpeed=hypotf(p->vel.x,p->vel.z);
    p->vel=v_add(p->vel,v_mul(direction,ACCELERATION*dt));
    int i=(int)(p-players);
    if(!tether_grapple_active(i)&&!tether_grapple_airborne(i))return;
    /* The higher grapple ceiling belongs to physical impulses, not the walk
     * motor. Steering may preserve existing momentum, but cannot manufacture
     * extra speed above walking speed. Braking input can still reduce it. */
    float limit=fmaxf(MOVE_SPEED,previousSpeed),speed=hypotf(p->vel.x,p->vel.z);
    if(speed>limit){p->vel.x*=limit/speed;p->vel.z*=limit/speed;}
}
static float tether_grapple_ground_friction(int i){
    if(i<0||i>=MAX_PLAYERS||!players[i].onGround)return FRICTION;
    float inverseMass=players[i].netGrappleInverseMass;
    Vector3 target,velocity;if(!tether_grapple_anchor(i,&target,&velocity))return FRICTION;
    if(netTransport.role!=NET_ROLE_CLIENT){WeakTetherCoupling c=weakTethers[i];if(!weak_tether_sample(i,&c))return FRICTION;inverseMass=c.anchored?0:1/fmaxf(c.mass,TETHER_REFERENCE_MASS);}
    // Walking brakes at 400 m/s2. Keep them for light pickups; heavy anchors
    // need finite physical drag below the grapple's 30 m/s2 acceleration.
    return inverseMass<=1/TETHER_GRAPPLE_PLAYER_MASS?fminf(FRICTION,12):FRICTION;
}
static void tether_grapple_pull(int i,float dt){
    if(netTransport.role!=NET_ROLE_CLIENT)return;
    Vector3 target,velocity;if(dt<=0||!tether_grapple_anchor(i,&target,&velocity))return;
    Vector3 error=world_delta(target,player_tether_hold_position(&players[i]));
    weak_tether_apply_player(i,weak_tether_impulse(error,v_sub(players[i].vel,velocity),players[i].netGrappleInverseMass,dt));
}
static void tether_grapple_receive(int i,unsigned flags,Vector3 target,Vector3 velocity,float inverseMass,uint32_t tick){
    Player *p=&players[i];
    bool active=(flags&NET_PLAYER_VISUAL_TETHER)&&!(flags&(NET_PLAYER_VISUAL_GOLD_TETHER|NET_PLAYER_VISUAL_DRAIN));
    p->netGrappleAnchorVelocity=active?velocity:(Vector3){0};p->netGrappleInverseMass=active?inverseMass:0;
    p->netGrappleAnchorTick=active?tick:0;
    p->netTetherVisualActive=(flags&NET_PLAYER_VISUAL_TETHER)!=0;
    p->netGoldTetherVisualActive=(flags&NET_PLAYER_VISUAL_GOLD_TETHER)!=0;
    p->netDrainVisualActive=(flags&NET_PLAYER_VISUAL_DRAIN)!=0;p->netTetherVisualTarget=target;
}
