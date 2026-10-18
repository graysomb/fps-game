/* Actor tether locks are separate from physical voxel constraints. */
static bool drain_target_valid(int from,int to) {
    if(to<0 || to>=activePlayers || to==from)return false;
    Player *p=&players[from],*t=&players[to];
    return p->respawn_timer<=0 && t->respawn_timer<=0 &&
        combatantLife[from]==LIFE_ALIVE && combatantLife[to]==LIFE_ALIVE &&
        (gameMode!=GAME_MODE_FIREFIGHT || combatantTeam[from]!=combatantTeam[to]) &&
        enemy_core_open(to);
}
static bool drain_clear(int from,int to) {
    Vector3 delta=world_delta(players[to].pos,players[from].pos);float distance=v_length(delta);
    if(distance>TETHER_RANGE || distance<.001f)return false;
    VoxelHit hit;
    if(!first_voxel_hit_detailed((Ray){players[from].pos,v_mul(delta,1/distance)},distance,-1,false,&hit))return true;
    // An enemy's own armor is part of the target, not intervening world cover.
    return enemy_equipment_owner(hit.id)==to || get_goliath_owner_of_voxel(hit.id)==to;
}
static bool start_tether_drain(int from) {
    if(netTransport.role==NET_ROLE_CLIENT)return false;
    Player *p=&players[from];Ray aim={p->pos,player_forward(p)};
    int best=-1;float nearest=TETHER_RANGE;
    for(int to=0;to<activePlayers;to++)if(drain_target_valid(from,to)){
        BoundingBox bounds=combatant_bounds_at(to,world_image(players[to].pos,p->pos),PLAYER_SIZE*.5f);
        RayCollision hit=GetRayCollisionBox(aim,expand_bounds(bounds,.15f));
        if(hit.hit && hit.distance<nearest && drain_clear(from,to)){best=to;nearest=hit.distance;}
    }
    if(best<0)return false;
    p->drainTarget=best+1;play_sfx_at(SFX_TETHER,p->pos,from,1,1);return true;
}
static void update_tether_drains(float dt) {
    if(netTransport.role==NET_ROLE_CLIENT || dt<=0)return;
    for(int from=0;from<activePlayers;from++){
        Player *p=&players[from];int to=p->drainTarget-1;
        if(to<0)continue;
        Vector3 delta=to<activePlayers?world_delta(players[to].pos,p->pos):(Vector3){0};
        if(!drain_target_valid(from,to) || p->tetherHolding || !drain_clear(from,to) ||
           v_dot(player_forward(p),v_norm(delta))<cosf(65*DEG2RAD)) {p->drainTarget=0;continue;}
        Player *target=&players[to];
        if(target->invuln_timer>0 || (enemyLab && enemyLabGod && to<firefightHumans))continue;
        // Foundry Matter is its lethal core health: leave one point rather than kill.
        float minimum=target->enemyType==ENEMY_TYPE_FOUNDRY?1:0;
        float amount=fminf(2.f*player_bullet_damage(from)*dt,fmaxf(0,target->matter-minimum));
        if(amount<=0)continue;
        tetherDrainDamage=true;apply_matter_damage(to,from,amount);tetherDrainDamage=false;
        if(target->matter<=minimum+.0001f){target->matter=minimum;}
    }
}
