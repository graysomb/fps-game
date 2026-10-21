/* Authoritative equipment behavior. All simulation runs once, outside viewports. */
static void eq_reset(int i){if(i<0 || i>=MAX_PLAYERS)return;aimStates[i]=(AimState){0};netLastAimInput[i]=0;equipment[i]=(Equipment){.gun={EQ_DEFAULT,EQ_EMPTY},.initialized=true};players[i].dynamicShotActive=false;}
static void eq_reset_all(void){for(int i=0;i<MAX_PLAYERS;i++)eq_reset(i);memset(eqRockets,0,sizeof(eqRockets));memset(eqTraces,0,sizeof(eqTraces));eqReceiveTick=0;memset(eqPartTick,0,sizeof(eqPartTick));}
static int eq_gun(int i){if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return EQ_DEFAULT;if(!equipment[i].initialized)eq_reset(i);return equipment[i].gun[equipment[i].selected];}
static bool eq_upgrade(int i,int flag){return i>=0 && i<MAX_PLAYERS && !is_player_bot(i) && (equipment[i].upgrades&flag)!=0;}
static bool eq_usable(int i){return i>=0 && i<MAX_PLAYERS && !is_player_bot(i) && netTransport.role!=NET_ROLE_CLIENT && gameState==GAME_STATE_PLAYING && players[i].respawn_timer<=0 && !ff_input_blocked(i) && !creativeModeActive && !(enemyLab&&enemyLabPaused);}
static void eq_spend(int i,float cost){players[i].matter=fmaxf(0,players[i].matter-cost);players[i].isExposed=players[i].matter<=0;}
static int eq_pickup_gun(int type){switch(type){case PICKUP_DYNAMIC_SHOT:return EQ_VOXEL;case PICKUP_SHOTGUN:return EQ_SHOTGUN;case PICKUP_CUTTER:return EQ_CUTTER;case PICKUP_RAILGUN:return EQ_RAIL;case PICKUP_ROCKET:return EQ_ROCKET;case PICKUP_SWORD:return EQ_BLADE;default:return -1;}}
static int eq_pickup_kind(int type){int gun=eq_pickup_gun(type);if(gun>=0)return eqGuns[gun].visual;return type==PICKUP_SWORD?WEAPON_SWORD:type==PICKUP_ADHESIVE?WEAPON_ADHESIVE_GRENADE:type==PICKUP_REINFORCED?WEAPON_REINFORCED_BUILDER:-1;}
static bool eq_pickup(int i,int type,bool replace){
    if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return false;
    eq_gun(i);Equipment *e=&equipment[i];int gun=eq_pickup_gun(type);
    if(gun>=0){
        if(e->gun[0]==gun || e->gun[1]==gun)return false;
        int slot=e->gun[0]==EQ_EMPTY?0:e->gun[1]==EQ_EMPTY?1:replace?e->selected:-1;
        if(slot<0)return false;aim_cancel(i);e->gun[slot]=gun;e->selected=slot;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;
        players[i].dynamicShotActive=gun==EQ_VOXEL;return true;
    }
    int flag=type==PICKUP_ADHESIVE?EQ_ADHESIVE:type==PICKUP_REINFORCED?EQ_BUILD:0;
    if(!flag || (e->upgrades&flag))return false;e->upgrades|=flag;return true;
}
static int eq_near_pickup(int i){int best=-1;float distance=1.4f;for(int p=0;p<MAX_PICKUPS;p++)if(pickups[p].active && eq_pickup_kind(pickups[p].type)>=0){float d=v_length(world_delta(pickups[p].pos,players[i].pos));if(d<distance){distance=d;best=p;}}return best;}
static bool eq_hostile(int i,int j){return i!=j && players[j].respawn_timer<=0 && (gameMode!=GAME_MODE_FIREFIGHT || combatantTeam[i]!=combatantTeam[j]);}
static void eq_damage(int i,int target,float amount,bool lethal,DeathCause cause,Vector3 origin){
    if(amount<=0 || !eq_hostile(i,target) || players[target].invuln_timer>0 || !enemy_core_open(target))return;
    bool exposed=players[target].isExposed;float matter=players[target].matter;
    DeathCause prev=deathContext;deathContext=cause;bool known=damageOriginKnown;Vector3 old=damageOrigin;damageOriginKnown=true;damageOrigin=origin;
    if(exposed)kill_player(target,i,true,false);
    else {apply_matter_damage(target,i,amount);if(lethal && matter<=amount && players[target].respawn_timer<=0)kill_player(target,i,true,false);}
    damageOriginKnown=known;damageOrigin=old;deathContext=prev;
}
static int eq_actor_ray(int owner,Ray ray,float *distance){
    int found=-1;for(int j=0;j<activePlayers;j++)if(j!=owner && players[j].respawn_timer<=0){
        BoundingBox box=combatant_bounds_at(j,world_image(players[j].pos,ray.position),PLAYER_SIZE*.5f);
        RayCollision hit=GetRayCollisionBox(ray,box);if(hit.hit && hit.distance>=0 && hit.distance<*distance){*distance=hit.distance;found=j;}
    }return found;
}
static void eq_trace(Vector3 a,Vector3 b,int gun,int owner){int slot=0;for(int n=0;n<32;n++)if(eqTraces[n].life<=0){slot=n;break;}eqTraces[slot]=(EquipmentTrace){a,b,.22f,gun,owner};}
static bool eq_cut(int owner,Ray ray,float reach){
    VoxelHit hit;if(!first_voxel_hit_detailed(ray,reach,-1,false,&hit))return false;
    int v=first_voxel_hit(ray,reach,-1);if(v<0)return false;
    if(clarity_shield_voxel(v)){clarity_shield_melee(v,owner);return true;}
    if(!voxels[v].simulate){
        Vector3 at=world_wrap(v_add(ray.position,v_mul(ray.direction,hit.t+.01f)));
        int x=(int)floorf(at.x/VOXEL_SIZE),y=(int)floorf(at.y/VOXEL_SIZE),z=(int)floorf(at.z/VOXEL_SIZE);
        foundry_trap_hit(v,owner,false);
        remove_static_voxels_in_region_recycle(x,x,y,y,z,z);mark_surface_neighbors(at);
        activate_static_voxels_near_region(x,x,y,y,z,z,owner);
    }else {Vector3 at=voxels[v].pos;enemy_detach(v);detach_goliath_armor_voxel(v);remove_voxel_index(v);mark_surface_neighbors(at);}
    return true;
}
// Bore a bounded, one-voxel-wide path. Actor hits are resolved first, so
// compaction and equipment removal cannot change which actors the shot hits.
static void eq_rail_terrain(int owner,Ray ray){
    for(int n=0;n<512;n++){
        VoxelHit hit;if(!first_voxel_hit_detailed(ray,eqGuns[EQ_RAIL].range,-1,false,&hit))break;
        int v=hit.id;if(v<0 || v>=voxel_count)break;
        Vector3 at=world_wrap(v_add(ray.position,v_mul(ray.direction,hit.t+.01f)));
        uint64_t id=voxels[v].identity;
        if(!voxels[v].simulate){
            if(foundry_trap_hit(v,owner,false))continue;
            remove_static_voxels_in_region_recycle(hit.x,hit.x,hit.y,hit.y,hit.z,hit.z);
        }else {enemy_detach(v);detach_goliath_armor_voxel(v);remove_voxel_index(v);}
        mark_surface_neighbors(at);
        // Fluids are activated, not deleted; do not repeatedly revisit them.
        if(enemy_voxel(id)>=0)break;
    }
    rebuild_voxel_hash();
}
typedef struct { uint64_t id; VoxelHit hit; Vector3 direction; float distance; bool wasStatic; } ShotgunImpact;
static void eq_shotgun_impacts(int owner,ShotgunImpact *hits,int count){
    if(!count)return;
    static UnitVoxelBuffer buffer;unit_voxel_buffer_clear(&buffer);
    refresh_static_voxel_beliefs();
    for(int v=0;v<voxel_count;v++)atomic_store_explicit(&activationClaims[v],0,memory_order_relaxed);
    // Collect once for the entire volley; use the same resistance and bounded
    // connected-voxel activation as an impacting dynamic voxel.
    for(int n=0;n<count;n++){
        ShotgunImpact *h=&hits[n];int v=enemy_voxel(h->id);if(v<0)continue;
        int gear=enemy_equipment_owner(v),armor=get_goliath_owner_of_voxel(v);
        if(gear>=0)ai_alert(gear,owner,players[owner].pos,AI_ALERT_HIT);
        if(armor>=0)ai_alert(armor,owner,players[owner].pos,AI_ALERT_HIT);
        if(gear>=0 || armor>=0 || voxel_is_fluid(&voxels[v]) || voxels[v].isBullet){h->id=0;continue;}
        if(!h->wasStatic)continue;
        if(foundry_trap_hit(v,owner,false)){
            uint64_t id=h->id;for(int k=n;k<count;k++)if(hits[k].id==id)hits[k].id=0;
            continue;
        }
        if(voxels[v].pendingActivation || voxels[v].activationCooldownFrames>0)continue;
        int start=buffer.count;
        if(collect_static_activation_cluster(v,owner,h->hit.x,h->hit.y,h->hit.z,
            (float)(VOXEL_ACTIVATION_RADIUS*VOXEL_ACTIVATION_RADIUS),&buffer)<=0)continue;
        if(!dynamic_belief_overcomes_static(eqShotgunActivationBelief,compute_cluster_freeze_belief(&buffer,start))){rollback_activation_buffer(&buffer,start);continue;}
        expand_activation_cluster_unbounded(&buffer,start,eqShotgunActivationBelief,owner);
    }
    if(buffer.count && (voxel_count+buffer.count>=MAX_VOXELS || particle_pool_count-free_particle_count+8*buffer.count>=MAX_PARTICLES || recycleQueueCount+buffer.count>=MAX_VOXELS))rollback_activation_buffer(&buffer,0);
    if(buffer.count){
        qsort(buffer.voxels,(size_t)buffer.count,sizeof(UnitVoxelSeed),compare_unit_voxel_seed);
        remove_buffered_static_voxels(&buffer);
        emit_unit_voxels_from_units(&buffer,false,true,-1);
        wake_fluid_near_unit_buffer(&buffer);rebuild_voxel_hash();refresh_activation_surfaces(&buffer);
        glue_dynamic_voxel_to_static_neighbors();refresh_static_voxel_beliefs();
    }
    int targets[EQ_PELLETS];
    for(int n=0;n<count;n++){ShotgunImpact *h=&hits[n];targets[n]=!h->id?-1:h->wasStatic?table_get(h->hit.x,h->hit.y,h->hit.z):enemy_voxel(h->id);}
    for(int n=0;n<count;n++){
        int v=targets[n];
        if(v<0 || !voxels[v].simulate)continue;
        if(voxels[v].sleeping)wake_sleeping_cluster_in_place(v);
        int members=build_glue_cluster_indices(v,glueClusterIndices);Vector3 impulse={0};
        for(int k=0;k<members;k++){
            int member=glueClusterIndices[k];voxels[member].activator=owner;
            for(int hit=n;hit<count;hit++)if(targets[hit]==member){
                float falloff=1-.75f*clampf(hits[hit].distance/eqGuns[EQ_SHOTGUN].range,0,1);
                impulse=v_add(impulse,v_mul(hits[hit].direction,eqShotgunImpulse/EQ_PELLETS*falloff));targets[hit]=-1;
            }
        }
        // Shared impact code conserves total momentum across connected mass and
        // updates shared particles once per volley, not once per pellet/corner.
        apply_tether_impact_to_dynamic_island(v,0,impulse);
    }
    play_sfx_at(SFX_IMPACT,world_wrap(v_add(players[owner].pos,v_mul(hits[0].direction,hits[0].distance))),owner,1,.8f);
}
static void eq_rocket_blast(EquipmentRocket rocket){
    float damages[MAX_COMBATANTS]={0};
    for(int j=0;j<activePlayers;j++)if(players[j].respawn_timer<=0){Vector3 d=world_delta(players[j].pos,rocket.pos);float dist=v_length(d);VoxelHit hit;
        if(dist<eqRocketOuter && (dist<.15f || !first_voxel_hit_detailed((Ray){rocket.pos,v_norm(d)},fmaxf(0,dist-.3f),-1,false,&hit)))damages[j]=eqGuns[EQ_ROCKET].damage*(1-clampf((dist-eqRocketInner)/(eqRocketOuter-eqRocketInner),0,1));
    }
    // Protection is evaluated before the terrain becomes debris.
    enemy_blast_terrain((EnemyCharge){.pos=rocket.pos,.owner=rocket.owner});
    for(int j=0;j<activePlayers;j++)if(damages[j]>0){
        if(j==rocket.owner && players[j].invuln_timer<=0){bool previous=explosionDamageContext;explosionDamageContext=true;DeathCause old=deathContext;deathContext=DEATH_ROCKET;
            if(players[j].matter<=damages[j] || players[j].isExposed)kill_player(j,j,false,false);else apply_matter_damage(j,j,damages[j]);deathContext=old;explosionDamageContext=previous;
        }else eq_damage(rocket.owner,j,damages[j],true,DEATH_ROCKET,rocket.pos);
        if(eq_hostile(rocket.owner,j))players[j].vel=v_add(players[j].vel,v_mul(v_norm(world_delta(players[j].pos,rocket.pos)),damages[j]*.1f));
    }
    hologram_blast(rocket.pos,4,100);clarity_effect((ClarityEffect){.id=clarityNextId++,.kind=2,.pos=rocket.pos,.attacker=rocket.owner},true);
}
static bool eq_fire(int i){
    if(!eq_usable(i))return false;int gun=eq_gun(i);if(gun<EQ_SHOTGUN || gun==EQ_CUTTER || gun==EQ_RAIL)return false;
    Equipment *e=&equipment[i];Player *p=&players[i];const GunDefinition *def=&eqGuns[gun];
    if(gun==EQ_BLADE){if(e->switchTime>0)return false;perform_melee(i);return p->meleeSwingActive;}
    if(e->switchTime>0 || e->cooldown[gun]>0 || p->isExposed || p->matter<def->cost)return false;
    Vector3 f=player_forward(p);Ray ray={p->pos,f};
    if(gun==EQ_SHOTGUN){
        ShotgunImpact hits[EQ_PELLETS];int hitCount=0;
        float damage[MAX_COMBATANTS]={0};Vector3 r=v_norm(v_cross(f,(Vector3){0,1,0})),up=v_cross(r,f);
        for(int n=0;n<EQ_PELLETS;n++){float angle=n*2*PI/EQ_PELLETS;Vector3 d=v_norm(v_add(f,v_add(v_mul(r,cosf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)),v_mul(up,sinf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)))));Ray pellet={p->pos,d};float length=def->range;VoxelHit hit;
            bool blocked=first_voxel_hit_detailed(pellet,length,-1,false,&hit);if(blocked)length=hit.t;
            int target=eq_actor_ray(i,pellet,&length);if(target>=0)damage[target]+=def->damage;
            else if(blocked)hits[hitCount++]=(ShotgunImpact){voxels[hit.id].identity,hit,d,length,!voxels[hit.id].simulate};
            eq_trace(p->pos,v_add(p->pos,v_mul(d,length)),gun,i);
        }
        for(int j=0;j<activePlayers;j++)if(damage[j]>0)eq_damage(i,j,damage[j],false,DEATH_SHOTGUN,p->pos);
        eq_shotgun_impacts(i,hits,hitCount);
    }else {
        int free=-1;for(int n=0;n<32;n++)if(!eqRockets[n].active){free=n;break;}if(free<0)return false;
        // Start at the eye and sweep from it: no spawning through nearby walls.
        eqRockets[free]=(EquipmentRocket){true,i,p->pos,v_mul(f,eqRocketSpeed),eqRocketLife};
    }
    eq_spend(i,def->cost);e->cooldown[gun]=def->cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,gun==EQ_ROCKET?.55f:.8f);ai_gunfire(i,p->pos);return true;
}
static void eq_command(int i,const NetInputCommand *c,float dt){
    if(i<0 || i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];eq_gun(i);
    aim_command(i,c);
    if(!eq_usable(i)){e->held=false;e->charge=0;return;}
    if((c->pressed&NET_INPUT_SWITCH) && !players[i].meleeSwingActive){int other=1-e->selected;if(e->gun[other]!=EQ_EMPTY){e->selected=other;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;}}
    if((c->pressed&NET_INPUT_EQUIP) && !players[i].meleeSwingActive){int n=eq_near_pickup(i);if(n>=0 && eq_pickup(i,pickups[n].type,true)){pickups[n].active=false;pickups[n].respawnTimer=pickups[n].type==PICKUP_RAILGUN||pickups[n].type==PICKUP_ROCKET?60:30;}}
    bool held=(c->held&NET_INPUT_FIRE)!=0;int gun=eq_gun(i);Player *p=&players[i];
    if(e->switchTime>0 || (p->isExposed && gun!=EQ_BLADE) || p->meleeSwingActive || e->cooldown[gun]>0){e->held=false;e->charge=0;return;}
    if((c->pressed&NET_INPUT_FIRE) && gun!=EQ_CUTTER && gun!=EQ_RAIL)FireVoxel(i);
    Ray ray={p->pos,player_forward(p)};
    if(gun==EQ_RAIL){
        if(held && p->matter>=eqGuns[gun].cost)e->charge=fminf(eqRailCharge,e->charge+dt);
        if(!held && e->held && e->charge>=eqRailCharge && p->matter>=eqGuns[gun].cost){
            for(int j=0;j<activePlayers;j++)if(eq_hostile(i,j)){
                bool hitActor=false;float period=worldFloorHalfSize*2;int images=worldPeriodic?(int)ceilf(eqGuns[gun].range/period)+1:0;
                Vector3 nearest=world_image(players[j].pos,p->pos);
                for(int x=-images;x<=images && !hitActor;x++)for(int z=-images;z<=images && !hitActor;z++){
                    Vector3 at=v_add(nearest,(Vector3){x*period,0,z*period});RayCollision hit=GetRayCollisionBox(ray,combatant_bounds_at(j,at,PLAYER_SIZE*.5f));
                    hitActor=hit.hit && hit.distance>=0 && hit.distance<=eqGuns[gun].range;
                }
                if(hitActor)eq_damage(i,j,eqGuns[gun].damage,true,DEATH_RAIL,p->pos);
            }
            eq_rail_terrain(i,ray);
            eq_trace(player_weapon_pose(i,false).muzzle,v_add(p->pos,v_mul(ray.direction,100)),EQ_RAIL,i);eq_spend(i,eqGuns[gun].cost);e->cooldown[gun]=eqGuns[gun].cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,.65f);ai_gunfire(i,p->pos);
        }
        if(!held)e->charge=0;
    }else if(gun==EQ_CUTTER && held && p->matter>0){
        float step=fminf(dt,p->matter/eqGuns[gun].cost),length=eqGuns[gun].range;VoxelHit hit;bool wall=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(wall)length=hit.t;
        int target=eq_actor_ray(i,ray,&length);e->beamEnd=v_add(p->pos,v_mul(ray.direction,length));
        if(target>=0)eq_damage(i,target,eqGuns[gun].damage*step,false,DEATH_CUTTER,p->pos);
        e->cutTime+=step;bool cut=false;
        for(int cuts=0;cuts<12 && e->cutTime>=eqCutterUnitSeconds;cuts++){
            e->cutTime-=eqCutterUnitSeconds;float stop=eqGuns[gun].range;VoxelHit surface;
            bool blocked=first_voxel_hit_detailed(ray,stop,-1,false,&surface);if(blocked)stop=surface.t;
            if(eq_actor_ray(i,ray,&stop)<0 && blocked)cut|=eq_cut(i,ray,eqGuns[gun].range);
        }
        if(cut){play_sfx_at(SFX_FIRE,p->pos,i,1,1.6f);ai_gunfire(i,p->pos);}
        eq_spend(i,eqGuns[gun].cost*step);
    }else e->cutTime=0;
    e->held=held;
}
static void eq_local(float dt){
    if(netTransport.role==NET_ROLE_CLIENT)return;
    for(int i=0;i<MAX_PLAYERS && i<activePlayers;i++)if(!is_player_bot(i) && (netTransport.role==NET_ROLE_OFFLINE || net_player_is_local(i))){
        NetInputCommand c={0};bool pad=playerInput[i]==INPUT_TYPE_GAMEPAD;int d=seat_device(i);
        if(pad){if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_TRIGGER_2))c.pressed|=NET_INPUT_FIRE;if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_LEFT))c.pressed|=NET_INPUT_SWITCH;if(IsGamepadButtonDown(d,GAMEPAD_BUTTON_RIGHT_TRIGGER_2))c.held|=NET_INPUT_FIRE;}
        else {if(!mouse_owns(i) && IsKeyPressed(i==0?KEY_LEFT_CONTROL:KEY_RIGHT_CONTROL))c.pressed|=NET_INPUT_FIRE;if(IsKeyPressed(i==0?KEY_Q:KEY_U))c.pressed|=NET_INPUT_SWITCH;if(IsKeyPressed(i==0?KEY_B:KEY_Y))c.pressed|=NET_INPUT_EQUIP;if(IsKeyDown(i==0?KEY_LEFT_CONTROL:KEY_RIGHT_CONTROL))c.held|=NET_INPUT_FIRE;}
        if(mouse_owns(i)){if(mouseCaptured&&mouseButtonsArmed&&IsMouseButtonPressed(MOUSE_BUTTON_LEFT))c.pressed|=NET_INPUT_FIRE;if(mouseCaptured&&mouseButtonsArmed&&IsMouseButtonDown(MOUSE_BUTTON_LEFT))c.held|=NET_INPUT_FIRE;if(mouseCaptured&&GetMouseWheelMove()!=0)c.pressed|=NET_INPUT_SWITCH;}
        aim_sample(netTransport.role==NET_ROLE_OFFLINE?i:net_local_index_for_slot(i),i,dt,&c);
        if(debugConsoleOpen || !IsWindowFocused()){c.held=c.pressed=0;equipment[i].charge=0;equipment[i].held=false;}
        eq_command(i,&c,dt);
    }
}
static void eq_tick(float dt){
    if(dt<=0)return;
    aim_tick(dt);
    for(int i=0;i<MAX_PLAYERS;i++){Equipment *e=&equipment[i];for(int g=0;g<EQ_GUN_COUNT;g++)e->cooldown[g]=fmaxf(0,e->cooldown[g]-dt);e->switchTime=fmaxf(0,e->switchTime-dt);e->buildTime=fmaxf(0,e->buildTime-dt);
        if(e->initialized && !is_player_bot(i))players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;
        weaponVisuals[i].charge=e->charge/.65f;weaponVisuals[i].activity=e->held||e->buildTime>0;
        if(players[i].respawn_timer>0){e->held=false;e->charge=0;}
    }
    for(int n=0;n<32;n++)eqTraces[n].life=fmaxf(0,eqTraces[n].life-dt);
    if(netTransport.role==NET_ROLE_CLIENT){for(int n=0;n<32;n++)if(eqRockets[n].active)eqRockets[n].pos=world_wrap(v_add(eqRockets[n].pos,v_mul(eqRockets[n].vel,dt)));return;}
    for(int n=0;n<32;n++){EquipmentRocket *r=&eqRockets[n];if(!r->active)continue;float remaining=fminf(dt,.25f);
        while(remaining>0 && r->active){float step=fminf(remaining,1.f/120);remaining-=step;float length=v_length(r->vel)*step;Ray ray={r->pos,v_norm(r->vel)};VoxelHit hit;
            bool collision=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(collision)length=hit.t;int actor=eq_actor_ray(r->owner,ray,&length);
            r->pos=world_wrap(v_add(r->pos,v_mul(ray.direction,length)));r->life-=step;
            if(collision || actor>=0 || r->pos.y<.1f){EquipmentRocket copy=*r;r->active=false;eq_rocket_blast(copy);}else if(r->life<=0)r->active=false;
        }
    }
}
static void eq_sword(int i){
    Ray center={players[i].pos,player_forward(&players[i])};
    for(int j=0;j<activePlayers;j++)if(eq_hostile(i,j)){Vector3 d=world_delta(players[j].pos,center.position);float length=v_length(d);if(length>eqSwordReach || v_dot(v_norm(d),center.direction)<cosf(eqSwordSweep*.5f*DEG2RAD))continue;
        VoxelHit wall;if(first_voxel_hit_detailed((Ray){center.position,v_norm(d)},fmaxf(0,length-.3f),-1,false,&wall))continue;
        eq_damage(i,j,eqSwordDamage,true,DEATH_SWORD,center.position);
    }
    // Sample a diagonal cutting plane densely, then cut through its depth.
    Vector3 right=v_norm(v_cross(center.direction,(Vector3){0,1,0})),up=v_cross(right,center.direction);
    Vector3 diagonal=v_norm(v_sub(right,up));int cuts=0;
    for(int depth=0;depth<3;depth++)for(int n=0;n<16 && cuts<EQ_SWORD_CUTS;n++){
        float angle=(-35+70*n/15.f)*DEG2RAD;
        Vector3 dir=v_norm(v_add(v_mul(center.direction,cosf(angle)),v_mul(diagonal,sinf(angle))));
        if(eq_cut(i,(Ray){center.position,dir},eqSwordReach))cuts++;
    }
    if(cuts)add_player_matter(&players[i],ff_harvest_amount(i));play_sfx_at(SFX_MELEE,center.position,i,1,.7f);
}
static bool eq_panel_positions(int i,Vector3 out[EQ_PANEL_CUBES]){
    Ray ray={players[i].pos,player_forward(&players[i])};VoxelHit hit;Vector3 at;
    if(first_voxel_hit_detailed(ray,TETHER_RANGE,-1,false,&hit))at=v_add(ray.position,v_add(v_mul(ray.direction,hit.t),v_mul(hit.normal,VOXEL_SIZE*.51f)));
    else {if(ray.direction.y>=-.001f)return false;float t=-ray.position.y/ray.direction.y;if(t>TETHER_RANGE)return false;at=v_add(ray.position,v_mul(ray.direction,t));at.y=VOXEL_SIZE*.5f;}
    at=world_wrap(at);at=(Vector3){(floorf(at.x/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.y/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.z/VOXEL_SIZE)+.5f)*VOXEL_SIZE};
    Vector3 side=fabsf(ray.direction.z)>=fabsf(ray.direction.x)?(Vector3){VOXEL_SIZE,0,0}:(Vector3){0,0,VOXEL_SIZE};bool valid=true;
    for(int n=0;n<EQ_PANEL_CUBES;n++){out[n]=world_wrap(v_add(at,v_add(v_mul(side,n%EQ_PANEL_SIDE-EQ_PANEL_SIDE/2),(Vector3){0,(n/EQ_PANEL_SIDE)*VOXEL_SIZE,0})));int x=floorf(out[n].x/VOXEL_SIZE),y=floorf(out[n].y/VOXEL_SIZE),z=floorf(out[n].z/VOXEL_SIZE);
        if(y<0 || occupied(x,y,z) || (!worldPeriodic && (fabsf(out[n].x)>FLOOR_SIZE || fabsf(out[n].z)>FLOOR_SIZE)))valid=false;
        BoundingBox cube={v_sub(out[n],(Vector3){.25f,.25f,.25f}),v_add(out[n],(Vector3){.25f,.25f,.25f})};
        for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && !voxels[v].isBullet){float radius=voxels[v].rest_edge*.5f;Vector3 center=world_image(voxels[v].pos,out[n]);if(CheckCollisionBoxes(cube,(BoundingBox){v_sub(center,(Vector3){radius,radius,radius}),v_add(center,(Vector3){radius,radius,radius})}))valid=false;}
        for(int j=0;j<activePlayers;j++)if((players[j].respawn_timer<=0 || combatantLife[j]==LIFE_DOWNED) && CheckCollisionBoxes(cube,combatant_bounds_at(j,world_image(players[j].pos,out[n]),PLAYER_SIZE*.5f)))valid=false;
    }return valid;
}
static bool eq_build(int i){
    if(!eq_usable(i))return false;Player *p=&players[i];float cost=3*ff_build_cost(i);Vector3 at[EQ_PANEL_CUBES];
    if(p->isExposed || p->matter<cost || GetTime()-p->last_build_time<BUILD_COOLDOWN_SECONDS || !eq_panel_positions(i,at) || voxel_count+EQ_PANEL_CUBES>=MAX_VOXELS || MAX_PARTICLES-particle_pool_count+free_particle_count<8*EQ_PANEL_CUBES)return false;
    uint64_t ids[EQ_PANEL_CUBES]={0};int made=0;
    for(int n=0;n<EQ_PANEL_CUBES;n++){int v=add_static_voxel_at_grid(floorf(at[n].x/VOXEL_SIZE),floorf(at[n].y/VOXEL_SIZE),floorf(at[n].z/VOXEL_SIZE),player_palette_color(i),0);if(v<0)break;ids[made++]=voxels[v].identity;voxels[v].reinforced=true;tag_owned_static_voxel(v,i);mark_surface(v);mark_surface_neighbors(at[n]);}
    if(made!=EQ_PANEL_CUBES){for(int n=0;n<made;n++){int v=enemy_voxel(ids[n]);if(v>=0)remove_voxel_index(v);}return false;}
    eq_spend(i,cost);p->last_build_time=GetTime();equipment[i].buildTime=.6f;play_sfx(SFX_SHIELD);return true;
}

static Shader eqBeamShader;
static bool eqBeamTried;
static void eq_beam_shutdown(void){if(eqBeamShader.id && eqBeamShader.id!=rlGetShaderIdDefault())UnloadShader(eqBeamShader);eqBeamShader=(Shader){0};eqBeamTried=false;}
static void eq_world(void){
    if(!eqBeamTried){eqBeamTried=true;eqBeamShader=LoadShader(NULL,"shaders/equipment_beam.frag");}
    bool shader=eqBeamShader.id && eqBeamShader.id!=rlGetShaderIdDefault();
    if(shader){float time=(float)GetTime();SetShaderValue(eqBeamShader,GetShaderLocation(eqBeamShader,"effectTime"),&time,SHADER_UNIFORM_FLOAT);BeginShaderMode(eqBeamShader);}

    for(int n=0;n<32;n++){EquipmentTrace *t=&eqTraces[n];if(t->life>0){
            Vector3 start=t->start;if(t->owner==drawViewPlayerIndex)start=player_weapon_pose(t->owner,true).muzzle;
            float width=t->gun==EQ_RAIL?.065f:.008f;
            DrawCylinderEx(start,t->end,width,width,8,Fade(t->gun==EQ_RAIL?PURPLE:ORANGE,t->life/.22f));
            if(t->gun==EQ_RAIL){DrawSphere(start,.10f*t->life/.22f,WHITE);DrawCylinderEx(start,t->end,.018f,.018f,6,Fade(WHITE,t->life/.22f));}
        }
        EquipmentRocket *r=&eqRockets[n];if(r->active){DrawCube(r->pos,.22f,.22f,.3f,ORANGE);DrawCylinderEx(r->pos,v_sub(r->pos,v_mul(v_norm(r->vel),.8f)),.1f,.02f,5,RED);}}
    for(int i=0;i<MAX_PLAYERS;i++){Equipment *e=&equipment[i];if(e->initialized && e->held && eq_gun(i)==EQ_CUTTER && players[i].respawn_timer<=0){WeaponPose pose=player_weapon_pose(i,i==drawViewPlayerIndex);DrawCylinderEx(pose.muzzle,e->beamEnd,.055f,.055f,8,LIME);DrawSphere(e->beamEnd,.11f,LIME);}
        if(e->charge>0 && eq_gun(i)==EQ_RAIL){Vector3 end=v_add(players[i].pos,v_mul(player_forward(&players[i]),100));DrawCylinderEx(player_weapon_pose(i,i==drawViewPlayerIndex).muzzle,end,.006f,.006f,4,Fade(PURPLE,.5f));}
        if(!equipmentEmission && i==drawViewPlayerIndex && eq_upgrade(i,EQ_BUILD) && players[i].respawn_timer<=0){Vector3 pos[EQ_PANEL_CUBES]={{0}};bool valid=eq_panel_positions(i,pos)&&players[i].matter>=3*ff_build_cost(i);if(v_length(pos[1])>.01f)for(int n=0;n<EQ_PANEL_CUBES;n++)DrawCubeWires(pos[n],VOXEL_SIZE,VOXEL_SIZE,VOXEL_SIZE,valid?GOLD:RED);}
    }
    if(shader)EndShaderMode();
}
static void eq_icon(int gun,int x,int y,Color c){
    if(gun==EQ_BLADE){DrawLineEx((Vector2){x-5,y+22},(Vector2){x+22,y-2},4,c);DrawLine(x-3,y+9,x+10,y+22,c);return;}
    DrawRectangle(x,y+5,25,9,c);DrawRectangle(x+17,y+13,5,10,c);
    if(gun==EQ_SHOTGUN){for(int n=0;n<3;n++)DrawRectangle(x-4,y+1+n*6,5,4,c);}
    if(gun==EQ_RAIL||gun==EQ_CUTTER){DrawRectangle(x-12,y+3,18,3,c);DrawRectangle(x-12,y+12,18,3,c);}
    if(gun==EQ_ROCKET)DrawRectangleLines(x-6,y,34,19,c);
    if(gun==EQ_VOXEL)DrawRectangleLines(x+22,y-2,12,17,c);
}
static void eq_hud(int i,int w,int h){
    if(i<0||i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];if(!e->initialized)return;int sz=w<700?12:16;
    int radarCenter=h-HUD_PADDING_Y-RADAR_PIXEL_RADIUS-HUD_BAR_THICKNESS-6;
    radarCenter=radarCenter-RADAR_PIXEL_RADIUS<HUD_PADDING_Y+HUD_BAR_HEIGHT?HUD_PADDING_Y+HUD_BAR_HEIGHT+RADAR_PIXEL_RADIUS+4:radarCenter;
    int y=radarCenter-RADAR_PIXEL_RADIUS-8-46;
    for(int s=0;s<2;s++){int x=14+s*105,g=e->gun[s];Color c=s==e->selected?GOLD:GRAY;DrawRectangle(x,y,98,46,Fade(BLACK,.7f));DrawRectangleLines(x,y,98,46,c);if(g!=EQ_EMPTY){eq_icon(g,x+20,y+2,c);DrawText(eqGuns[g].name,x+5,y+27,sz,c);}else DrawText("Empty",x+5,y+27,sz,GRAY);}
    int g=eq_gun(i);DrawText(TextFormat("%.0f Matter%s  %.1fs",eqGuns[g].cost,g==EQ_CUTTER?"/s":"",e->cooldown[g]),14,y-19,sz,RAYWHITE);
    if(g==EQ_RAIL && e->charge>0 && aimStates[i].blend<.02f)DrawRectangle(w/2-40,h/2+35,(int)(80*e->charge/.65f),5,PURPLE);
    DrawText(TextFormat("%s %s %s","",e->upgrades&EQ_ADHESIVE?"STICKY":"",e->upgrades&EQ_BUILD?"PANEL":""),230,y+27,sz,GOLD);
    int n=eq_near_pickup(i);if(n>=0 && eq_pickup_gun(pickups[n].type)>=0 && e->gun[0]!=eq_pickup_gun(pickups[n].type) && e->gun[1]!=eq_pickup_gun(pickups[n].type)){const char *key=playerInput[i]==INPUT_TYPE_GAMEPAD?"Hold X":i==0?"B":"Y";const char *label=TextFormat("%s: equip %s",key,pickup_type_label(pickups[n].type));DrawText(label,(w-MeasureText(label,sz))/2,h/2+55,sz,WHITE);}
}
static void eq_layout(void){
    if(useCustomMap||firefightPractice)return;bool foundry=gameMode==GAME_MODE_FIREFIGHT&&firefightFoundry;if(!foundry&&!clarity_crosscut())return;
    const int types[]={PICKUP_SHOTGUN,PICKUP_CUTTER,PICKUP_RAILGUN,PICKUP_ROCKET,PICKUP_SWORD,PICKUP_ADHESIVE,PICKUP_REINFORCED};
    const Vector3 anchors[]={{20,0,14},{28,0,14},{6,0,3},{3,0,6},{20,0,25},{28,0,25},{32,0,7}};
    for(int t=0;t<7;t++){Vector3 at;if(!clarity_supply_position(anchors[t],0,&at))continue;for(int q=0;q<(foundry?2:4);q++){int free=-1;for(int n=0;n<MAX_PICKUPS;n++)if(!pickups[n].active && pickups[n].respawnTimer<=0){free=n;break;}if(free<0)return;Vector3 pos=crosscut_rotate(at,foundry?q*2:q);Vector3 clear;if(!clarity_supply_position(pos,0,&clear))continue;pickups[free]=(Pickup){.pos=clear,.active=true,.type=types[t]};}}
}
#include "equipment_net.inc"

#include "aiming.inc"
