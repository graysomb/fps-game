static int ai_level(int i) {return playerInput[i]==INPUT_TYPE_BOT_EASY?0:playerInput[i]==INPUT_TYPE_BOT_HARD?2:1;}
static float ai_reaction(int i) {const float value[]={.6f,.4f,.3f};return value[ai_level(i)];}
static float ai_pause(int i) {const float value[]={1.2f,1,.8f};return value[ai_level(i)];}
static float ai_random(int i) {
    uint32_t *s=&aiCombat[i].random;*s=*s*1664525u+1013904223u;return (*s>>8)/16777216.f;
}
static void ai_reset(int i) {
    if(i<0 || i>=MAX_COMBATANTS)return;
    players[i].drainTarget=0;players[i].netDrainVisualActive=false;
    for(int j=0;j<MAX_COMBATANTS;j++)if(players[j].drainTarget==i+1)players[j].drainTarget=0;
    aiCombat[i]=(AiCombat){.target=-1,.pendingTarget=-1,.intent=-1,.lastFire=-1000,.initialized=true,.random=0x51a700u+(uint32_t)i*977u,.lastPosition=players[i].pos};
    botStates[i]=(BotState){.targetIndex=-1};
    // Invalidate memory when a slot is reused, including another bot's target.
    for(int j=0;j<MAX_COMBATANTS;j++)if(j!=i) {
        if(aiCombat[j].target==i){aiCombat[j].target=-1;aiCombat[j].memory=0;}
        if(aiCombat[j].pendingTarget==i)aiCombat[j].pendingTarget=-1;
    }
}
static bool ai_hostile(int i,int j) {
    return i>=0 && i<activePlayers && j>=0 && j<activePlayers && i!=j && players[j].respawn_timer<=0 &&
        (gameMode!=GAME_MODE_FIREFIGHT || combatantTeam[i]!=combatantTeam[j]);
}
static void ai_cancel_burst(int i) {
    AiCombat *a=&aiCombat[i];if(a->burstShots)a->recovery=fmaxf(a->recovery,ai_pause(i));a->burstShots=0;
}
static void ai_select(int i,int target,Vector3 position,AiAlert source) {
    AiCombat *a=&aiCombat[i];
    if(a->target!=target){ai_cancel_burst(i);a->reaction=ai_reaction(i);a->tracking=0;a->visible=false;a->wasVisible=false;a->commit=1.5f;a->intentTime=0;}
    a->target=target;a->known=position;a->memory=4;a->source=source;
}
static void ai_alert(int victim,int attacker,Vector3 position,AiAlert source) {
    if(netTransport.role==NET_ROLE_CLIENT || !ai_hostile(victim,attacker) || !is_player_bot(victim) || players[victim].respawn_timer>0)return;
    AiCombat *a=&aiCombat[victim];if(!a->initialized)ai_reset(victim);
    if(source==AI_ALERT_HIT)a->tracking=0;
    if(enemy_is_special(victim) && enemySpecial[victim].windup>0) {
        if(a->pendingTarget>=0 && a->pendingSource==AI_ALERT_HIT && source!=AI_ALERT_HIT)return;
        a->pendingTarget=attacker;a->pendingPosition=position;a->pendingSource=source;return;
    }
    bool override=source==AI_ALERT_HIT && a->overrideDelay<=0;
    if(a->target==attacker || !ai_hostile(victim,a->target) || a->memory<=0 || override) {
        if(override && a->target!=attacker)a->overrideDelay=1;
        ai_select(victim,attacker,position,source);
    }
}
static void ai_gunfire(int shooter,Vector3 origin) {
    if(netTransport.role==NET_ROLE_CLIENT)return;
    for(int i=0;i<activePlayers;i++)if(is_player_bot(i) && ai_hostile(i,shooter)) {
        float d=v_length(world_delta(origin,players[i].pos));
        if(d<=8 || (d<=24 && !is_view_occluded_by_voxels(players[i].pos,origin)))ai_alert(i,shooter,origin,AI_ALERT_HEARD);
    }
}
static void ai_near_miss(int projectile,Vector3 start,Vector3 end) {
    Voxel *v=&voxels[projectile];if(!v->isBullet || netTransport.role==NET_ROLE_CLIENT)return;
    Vector3 path=v_sub(end,start);float distance=v_length(path);if(distance<.0001f)return;
    VoxelHit wall;if(first_voxel_hit_detailed((Ray){start,v_mul(path,1/distance)},distance,projectile,false,&wall))end=v_add(start,v_mul(path,wall.t/distance));
    // Stop the audible path at the first actor, too.
    path=v_sub(end,start);distance=v_length(path);
    for(int j=0;j<activePlayers;j++)if(j!=v->owner && players[j].respawn_timer<=0 && distance>.0001f) {
        BoundingBox b=combatant_bounds_at(j,world_image(players[j].pos,start),PLAYER_SIZE*.5f);
        RayCollision c=GetRayCollisionBox((Ray){start,v_mul(path,1/distance)},b);
        if(c.hit && c.distance<distance){end=v_add(start,v_mul(path,c.distance/distance));path=v_sub(end,start);distance=c.distance;}
    }
    for(int i=0;i<activePlayers;i++)if(is_player_bot(i) && ai_hostile(i,v->owner)) {
        Vector3 point=world_image(players[i].pos,start),segment=v_sub(end,start);
        float t=clampf(v_dot(v_sub(point,start),segment)/fmaxf(v_dot(segment,segment),.000001f),0,1);
        Vector3 closest=v_add(start,v_mul(segment,t));
        if(v_length(v_sub(point,closest))<=1.5f && !is_view_occluded_by_voxels(point,closest))
            ai_alert(i,v->owner,v->shotOrigin,AI_ALERT_NEAR);
    }
}
static void ai_tick(int i,float dt) {
    AiCombat *a=&aiCombat[i];if(!a->initialized)ai_reset(i);
    dt=fmaxf(0,dt);a->clock+=dt;a->memory=fmaxf(0,a->memory-dt);a->commit=fmaxf(0,a->commit-dt);
    a->overrideDelay=fmaxf(0,a->overrideDelay-dt);a->reaction=fmaxf(0,a->reaction-dt);a->recovery=fmaxf(0,a->recovery-dt);
    a->intentTime=fmaxf(0,a->intentTime-dt);a->jumpDelay=fmaxf(0,a->jumpDelay-dt);
    bool locked=enemy_is_special(i) && enemySpecial[i].windup>0;
    if(!locked && a->pendingTarget>=0){int target=a->pendingTarget;a->pendingTarget=-1;ai_alert(i,target,a->pendingPosition,a->pendingSource);}
    int nearest=find_nearest_aware_enemy(i,NULL);
    bool valid=ai_hostile(i,a->target),visible=valid && check_ai_awareness_of_player(i,a->target,NULL);
    if(!locked && nearest>=0 && (!valid || a->memory<=0 || (a->commit<=0 &&
        v_length(world_delta(players[nearest].pos,players[i].pos))<.7f*v_length(world_delta(a->known,players[i].pos))))) {
        ai_select(i,nearest,players[nearest].pos,AI_ALERT_SEEN);valid=true;visible=true;
    }
    if(!valid || (!visible && a->memory<=0 && !locked)){a->target=-1;a->memory=0;visible=false;}
    if(visible) {
        a->known=players[a->target].pos;a->memory=4;
        if(!a->wasVisible){a->reaction=fmaxf(a->reaction,ai_reaction(i));a->tracking=0;}
        Vector3 velocity=players[a->target].vel;velocity.y=0;
        if(v_length(velocity)>1 && v_length(a->lastTargetVelocity)>1 && v_dot(v_norm(velocity),v_norm(a->lastTargetVelocity))<0)a->tracking=0;
        a->lastTargetVelocity=velocity;a->tracking=fminf(1,a->tracking+dt);
    } else {ai_cancel_burst(i);a->tracking=0;a->reaction=fmaxf(a->reaction,ai_reaction(i));}
    a->visible=visible;a->wasVisible=visible;
    BotState *b=&botStates[i];b->targetIndex=a->target;b->hasTarget=a->target>=0;b->lastKnownTargetPos=a->known;b->lostTrackTimer=a->memory;
    Vector3 travel=world_delta(players[i].pos,a->lastPosition);travel.y=0;
    if(v_length(travel)<dt*.2f && v_length(players[i].vel)>.5f)a->stuck+=dt;else a->stuck=0;
    a->lastPosition=players[i].pos;
}
static bool ai_shot_allowed(int i) {
    if(!aiCombat[i].initialized)ai_reset(i);
    return aiCombat[i].clock-aiCombat[i].lastFire>=1.0/3.0-1e-7;
}
static void ai_shot_fired(int i) {
    AiCombat *a=&aiCombat[i];a->lastFire=a->clock;
    if(++a->burstShots>=(ai_level(i)==0?2:3)){a->burstShots=0;a->recovery=ai_pause(i);}
}
static bool ai_ready_to_fire(int i) {
    AiCombat *a=&aiCombat[i];if(!a->visible || a->reaction>0 || a->recovery>0 || !ai_shot_allowed(i))return false;
    Vector3 to=v_norm(world_delta(a->known,players[i].pos));return v_dot(player_forward(&players[i]),to)>=cosf(8*DEG2RAD);
}
static void ai_aim(int i,Vector3 target,float dt) {
    AiCombat *a=&aiCombat[i];const float initial[]={4,3,2.5f},settled[]={2,1.5f,1};int level=ai_level(i);
    Vector3 velocity=players[i].vel;velocity.y=0;
    a->error=initial[level]+(settled[level]-initial[level])*a->tracking+(v_length(velocity)>.5f?1.5f:0);
    a->aimTimer-=dt;
    if(a->aimTimer<=0){float angle=ai_random(i)*2*PI,radius=sqrtf(ai_random(i));a->goalYaw=cosf(angle)*radius;a->goalPitch=sinf(angle)*radius;a->aimTimer+=.25f;if(a->aimTimer<=0)a->aimTimer=.25f;}
    float blend=1-expf(-dt*8);a->aimYaw+=(a->goalYaw-a->aimYaw)*blend;a->aimPitch+=(a->goalPitch-a->aimPitch)*blend;
    float yaw,pitch;dir_to_yaw_pitch(world_delta(target,players[i].pos),&yaw,&pitch);
    players[i].yaw+=clampf(angle_diff_deg(players[i].yaw,yaw+a->aimYaw*a->error),-180*dt,180*dt);
    players[i].pitch+=clampf(pitch+a->aimPitch*a->error-players[i].pitch,-120*dt,120*dt);
}
