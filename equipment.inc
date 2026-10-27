/* Authoritative equipment behavior. All simulation runs once, outside viewports. */
static void eq_reset(int i){if(i<0 || i>=MAX_COMBATANTS)return;swordNextDirection[i]=swordDirection[i]=0;if(!swordWaves[i].emitted)swordWaves[i].active=false;aimStates[i]=(AimState){0};netLastAimInput[i]=0;equipment[i]=(Equipment){.gun={EQ_DEFAULT,EQ_EMPTY},.initialized=true};players[i].dynamicShotActive=false;players[i].grappleInputHeld=false;players[i].grappleAirborne=false;players[i].netGrappleInverseMass=0;if(i<MAX_PLAYERS)weakTethers[i]=(WeakTetherCoupling){0};players[i].netGrappleAnchorTick=0;players[i].netGrappleAnchorVelocity=(Vector3){0};weaponVisuals[i].deployment=weaponVisuals[i].build_pulse=weaponVisuals[i].cooldown=0;}
static void eq_reset_all(void){for(int i=0;i<MAX_COMBATANTS;i++)eq_reset(i);memset(swordWaves,0,sizeof(swordWaves));memset(swordSequences,0,sizeof(swordSequences));memset(eqRockets,0,sizeof(eqRockets));memset(eqTraces,0,sizeof(eqTraces));eqReceiveTick=0;memset(eqPartTick,0,sizeof(eqPartTick));}
static int eq_gun(int i){if(i<0 || i>=MAX_COMBATANTS || (is_player_bot(i)&&!fp_carrier(i)))return EQ_DEFAULT;if(!equipment[i].initialized)eq_reset(i);return equipment[i].gun[equipment[i].selected];}
static bool eq_upgrade(int i,int flag){return i>=0 && i<MAX_COMBATANTS && (!is_player_bot(i)||fp_carrier(i)) && (equipment[i].upgrades&flag)!=0;}
static bool eq_usable(int i){return i>=0 && i<MAX_COMBATANTS && (!is_player_bot(i)||fp_carrier(i)) && netTransport.role!=NET_ROLE_CLIENT && gameState==GAME_STATE_PLAYING && players[i].respawn_timer<=0 && !ff_input_blocked(i) && !creativeModeActive && !(enemyLab&&enemyLabPaused);}
static void eq_spend(int i,float cost){players[i].matter=fmaxf(0,players[i].matter-cost);players[i].isExposed=players[i].matter<=0;}
// Charge time represents Matter already transferred into the railgun capacitor.
// Refunds return owned Matter, so Glob's restoration suppression does not apply.
static void eq_cancel_charge(int i){
    Equipment *e=&equipment[i];Player *p=&players[i];
    if(netTransport.role!=NET_ROLE_CLIENT && e->charge>0 && p->respawn_timer<=0 && combatantLife[i]==LIFE_ALIVE){
        p->matter=fminf(player_max_matter(p),p->matter+eqGuns[EQ_RAIL].cost*e->charge/eqRailCharge);
        p->isExposed=p->matter<=0;
    }
    e->charge=0;e->held=false;
}
static int eq_pickup_gun(int type){switch(type){case PICKUP_DYNAMIC_SHOT:return EQ_VOXEL;case PICKUP_SHOTGUN:return EQ_SHOTGUN;case PICKUP_CUTTER:return EQ_CUTTER;case PICKUP_RAILGUN:return EQ_RAIL;case PICKUP_ROCKET:return EQ_ROCKET;case PICKUP_GLOB:return EQ_GLOB;case PICKUP_GLOCK:return EQ_GLOCK;case PICKUP_DMR:return EQ_DMR;case PICKUP_SWORD:return EQ_BLADE;default:return -1;}}
static int eq_pickup_kind(int type){int gun=eq_pickup_gun(type);if(gun>=0)return eqGuns[gun].visual;return type==PICKUP_BLACK_HOLE?WEAPON_BLACK_HOLE_BOMB:type==PICKUP_IGLOO?WEAPON_IGLOO_GRENADE:type==PICKUP_SWORD?WEAPON_SWORD:type==PICKUP_ADHESIVE?WEAPON_ADHESIVE_GRENADE:type==PICKUP_REINFORCED?WEAPON_REINFORCED_BUILDER:-1;}
static bool eq_pickup(int i,int type,bool replace){
    if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return false;
    eq_gun(i);Equipment *e=&equipment[i];int gun=eq_pickup_gun(type);
    if(gun>=0){
        if(e->gun[0]==gun || e->gun[1]==gun)return false;
        int slot=e->gun[0]==EQ_EMPTY?0:e->gun[1]==EQ_EMPTY?1:replace?e->selected:-1;
        if(slot<0)return false;aim_cancel(i);eq_cancel_charge(i);e->gun[slot]=gun;e->selected=slot;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;
        players[i].dynamicShotActive=gun==EQ_VOXEL;return true;
    }
    int flag=type==PICKUP_ADHESIVE?EQ_ADHESIVE:type==PICKUP_REINFORCED?EQ_BUILD:0;
    if(!flag || (e->upgrades&flag))return false;e->upgrades|=flag;return true;
}
static int eq_near_pickup(int i){int best=-1;float distance=1.4f;for(int p=0;p<MAX_PICKUPS;p++)if(pickups[p].active && eq_pickup_kind(pickups[p].type)>=0 && pickups[p].type!=PICKUP_BLACK_HOLE && pickups[p].type!=PICKUP_IGLOO){float d=v_length(world_delta(pickups[p].pos,players[i].pos));if(d<distance){distance=d;best=p;}}return best;}
static bool eq_hostile(int i,int j){return i!=j && players[j].respawn_timer<=0 && (gameMode!=GAME_MODE_FIREFIGHT || combatantTeam[i]!=combatantTeam[j]);}
static void eq_damage(int i,int target,float amount,bool lethal,DeathCause cause,Vector3 origin){
    if(amount<=0 || !eq_hostile(i,target) || players[target].invuln_timer>0 || !enemy_core_open(target))return;
    if(go_combat_actor(target)){go_combat_damage(amount,cause==DEATH_GLOB);if(cause==DEATH_GLOB)goFight.suppressed=5;return;}
    bool exposed=players[target].isExposed;float matter=players[target].matter;
    DeathCause prev=deathContext;deathContext=cause;bool known=damageOriginKnown;Vector3 old=damageOrigin;damageOriginKnown=true;damageOrigin=origin;
    if(cause==DEATH_GLOB){
        // Suppression never executes exposed actors, including the Foundry's
        // special core-health path which normally kills at zero Matter.
        if(!exposed){if(players[target].enemyType==ENEMY_TYPE_FOUNDRY)amount=fminf(amount,fmaxf(0,enemySpecial[target].coreHealth-.01f));apply_matter_damage(target,i,amount);}
    }
    else if(exposed)kill_player(target,i,true,false);
    else {apply_matter_damage(target,i,amount);if(lethal && matter<=amount && players[target].respawn_timer<=0)kill_player(target,i,true,false);}
    damageOriginKnown=known;damageOrigin=old;deathContext=prev;
}
static int eq_actor_ray(int owner,Ray ray,float *distance){
    int found=-1;float period=worldFloorHalfSize*2;int images=worldPeriodic?(int)ceilf(*distance/period):0;
    for(int j=0;j<activePlayers;j++)if(j!=owner && players[j].respawn_timer<=0){
        Vector3 nearest=world_image(players[j].pos,ray.position);
        for(int x=-images;x<=images;x++)for(int z=-images;z<=images;z++){
            Vector3 at=v_add(nearest,(Vector3){x*period,0,z*period});
            RayCollision hit=GetRayCollisionBox(ray,combatant_bounds_at(j,at,PLAYER_SIZE*.5f));
            if(hit.hit && hit.distance>=0 && hit.distance<*distance){*distance=hit.distance;found=j;}
        }
    }return found;
}
static void eq_trace(Vector3 a,Vector3 b,int gun,int owner){int slot=0;for(int n=0;n<32;n++)if(eqTraces[n].life<=0){slot=n;break;}eqTraces[slot]=(EquipmentTrace){a,b,.22f,gun,owner};}
static bool eq_cut(int owner,Ray ray,float reach){
    VoxelHit hit;if(!first_voxel_hit_detailed(ray,reach,-1,false,&hit))return false;
    int v=first_voxel_hit(ray,reach,-1);if(v<0)return false;
    if(go_shield_hit(v,0))return true;
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
        if(go_shield_hit(v,eqGuns[EQ_RAIL].damage))break;
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
    float shieldDamage[GM_COUNT]={0};int shieldIndex[GM_COUNT];for(int p=0;p<GM_COUNT;p++)shieldIndex[p]=-1;
    for(int n=0;n<count;n++){int v=enemy_voxel(hits[n].id),p=go_shield_part(v);
        if(p>=0 && go_shield_hit(v,0)){shieldDamage[p]+=eqGuns[EQ_SHOTGUN].damage;shieldIndex[p]=v;hits[n].id=0;}}
    for(int p=0;p<GM_COUNT;p++)if(shieldIndex[p]>=0)go_shield_hit(shieldIndex[p],shieldDamage[p]);
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
                float falloff=1;
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
    ducks_blast(rocket.pos,4);
    float damages[MAX_COMBATANTS]={0};
    for(int j=0;j<activePlayers;j++)if(players[j].respawn_timer<=0){Vector3 d=world_delta(players[j].pos,rocket.pos);float dist=v_length(d);VoxelHit hit;
        if(dist<eqRocketOuter && (dist<.15f || !first_voxel_hit_detailed((Ray){rocket.pos,v_norm(d)},fmaxf(0,dist-.3f),-1,false,&hit)))damages[j]=eqGuns[EQ_ROCKET].damage*(1-clampf((dist-eqRocketInner)/(eqRocketOuter-eqRocketInner),0,1));
    }
    go_shield_blast(rocket.pos,eqRocketOuter,eqGuns[EQ_ROCKET].damage);
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
#include "equipment_expansion.inc"
static bool eq_fire(int i){
    if(!eq_usable(i))return false;int gun=eq_gun(i);if(gun<EQ_SHOTGUN || gun==EQ_CUTTER || gun==EQ_RAIL)return false;
    Equipment *e=&equipment[i];Player *p=&players[i];const GunDefinition *def=&eqGuns[gun];
    if(gun==EQ_BLADE){if(e->switchTime>0 || p->meleeSwingActive)return false;return perform_sword(i);}
    if(e->switchTime>0 || e->cooldown[gun]>0 || p->isExposed || p->matter<def->cost)return false;
    Vector3 f=player_forward(p);Ray ray={p->pos,f};
    if(gun==EQ_SHOTGUN){
        ShotgunImpact hits[EQ_PELLETS];int hitCount=0;
        float damage[MAX_COMBATANTS]={0};Vector3 r=v_norm(v_cross(f,(Vector3){0,1,0})),up=v_cross(r,f);
        for(int n=0;n<EQ_PELLETS;n++){float angle=n*2*PI/EQ_PELLETS;Vector3 d=v_norm(v_add(f,v_add(v_mul(r,cosf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)),v_mul(up,sinf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)))));Ray pellet={p->pos,d};float length=def->range;VoxelHit hit;
            bool blocked=first_voxel_hit_detailed(pellet,length,-1,false,&hit);if(blocked)length=hit.t;
            int target=eq_actor_ray(i,pellet,&length);if(target>=0)damage[target]+=def->damage;
            else if(blocked)hits[hitCount++]=(ShotgunImpact){voxels[hit.id].identity,hit,d,length,!voxels[hit.id].simulate};
            ducks_ray(pellet,length,0,-1,false);
            eq_trace(p->pos,v_add(p->pos,v_mul(d,length)),gun,i);
        }
        for(int j=0;j<activePlayers;j++)if(damage[j]>0)eq_damage(i,j,damage[j],true,DEATH_SHOTGUN,p->pos);
        eq_shotgun_impacts(i,hits,hitCount);
    }else if(gun==EQ_GLOB||gun==EQ_GLOCK||gun==EQ_DMR){if(!ex_fire(i,gun,ray))return false;
    }else {
        int free=-1;for(int n=0;n<32;n++)if(!eqRockets[n].active){free=n;break;}if(free<0)return false;
        // Start at the eye and sweep from it: no spawning through nearby walls.
        eqRockets[free]=(EquipmentRocket){true,i,p->pos,v_mul(f,eqRocketSpeed),eqRocketLife};
    }
    eq_spend(i,def->cost);e->cooldown[gun]=def->cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,gun==EQ_ROCKET?.55f:gun==EQ_GLOCK?1.25f:gun==EQ_GLOB?1.6f:gun==EQ_DMR?.95f:.8f);ai_gunfire(i,p->pos);return true;
}
static void eq_command(int i,const NetInputCommand *c,float dt){
    if(i<0 || i>=MAX_COMBATANTS)return;Equipment *e=&equipment[i];eq_gun(i);
    aim_command(i,c);
    eq_build_command(i,c);
    if(!eq_usable(i)){eq_cancel_charge(i);return;}
    if(c->pressed&NET_INPUT_GRENADE_CYCLE)ex_cycle(i);
    if((c->pressed&NET_INPUT_SWITCH) && !players[i].meleeSwingActive){int other=1-e->selected;if(e->gun[other]!=EQ_EMPTY){eq_cancel_charge(i);e->selected=other;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;}}
    if((c->pressed&NET_INPUT_EQUIP) && !players[i].meleeSwingActive && !fp_interact(i)){int n=eq_near_pickup(i);if(n>=0 && eq_pickup(i,pickups[n].type,true)){pickups[n].active=false;pickups[n].respawnTimer=pickups[n].type==PICKUP_RAILGUN||pickups[n].type==PICKUP_ROCKET||pickups[n].type==PICKUP_DMR?60:30;}}
    bool held=(c->held&NET_INPUT_FIRE)!=0;int gun=eq_gun(i);Player *p=&players[i];
    if(e->switchTime>0 || (p->isExposed && gun!=EQ_BLADE && !(gun==EQ_RAIL && e->charge>0)) || p->meleeSwingActive || e->cooldown[gun]>0){eq_cancel_charge(i);return;}
    if((c->pressed&NET_INPUT_FIRE) && gun!=EQ_CUTTER && gun!=EQ_RAIL)FireVoxel(i);
    if(gun==EQ_GLOB && held)eq_fire(i);
    Ray ray={p->pos,player_forward(p)};
    if(gun==EQ_RAIL){
        if(held && e->charge<eqRailCharge){
            float next=fminf(eqRailCharge,e->charge+fmaxf(0,dt));
            float needed=(next-e->charge)*eqGuns[gun].cost/eqRailCharge;
            float loaded=fminf(p->matter,needed);
            e->charge=loaded>=needed?next:e->charge+loaded*eqRailCharge/eqGuns[gun].cost;
            eq_spend(i,loaded);
            // Repeated float transfers must not strand a fully paid shot just
            // below the charge threshold (or leave a microscopic Matter shield).
            if(e->charge>=eqRailCharge-1e-5f){
                e->charge=eqRailCharge;
                if(p->matter<1e-4f)eq_spend(i,p->matter);
            }
        }
        if(!held && e->held && e->charge>=eqRailCharge){
            ducks_ray(ray,eqGuns[gun].range,0,-1,true);
            float shieldStop=eqGuns[gun].range;
            uint64_t passed[goCount>0?goCount:1];int passCount=0;VoxelHit barrier;
            while(goCombatEnabled&&passCount<goCount&&go_body_ray(ray,shieldStop,-1,passed,passCount,&barrier)){
                if(go_shielded(barrier.id)){shieldStop=barrier.t;break;}
                passed[passCount++]=voxels[barrier.id].identity;
            }
            for(int j=0;j<activePlayers;j++)if(eq_hostile(i,j)){
                bool hitActor=false;float period=worldFloorHalfSize*2;int images=worldPeriodic?(int)ceilf(eqGuns[gun].range/period)+1:0;
                Vector3 nearest=world_image(players[j].pos,p->pos);
                for(int x=-images;x<=images && !hitActor;x++)for(int z=-images;z<=images && !hitActor;z++){
                    Vector3 at=v_add(nearest,(Vector3){x*period,0,z*period});RayCollision hit=GetRayCollisionBox(ray,combatant_bounds_at(j,at,PLAYER_SIZE*.5f));
                    hitActor=hit.hit && hit.distance>=0 && hit.distance<=shieldStop;
                }
                if(hitActor)eq_damage(i,j,eqGuns[gun].damage,true,DEATH_RAIL,p->pos);
            }
            eq_rail_terrain(i,ray);
            eq_trace(player_weapon_pose(i,false).muzzle,v_add(p->pos,v_mul(ray.direction,shieldStop)),EQ_RAIL,i);e->charge=0;e->cooldown[gun]=eqGuns[gun].cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,.65f);ai_gunfire(i,p->pos);
        }
        if(!held)eq_cancel_charge(i);
    }else if(gun==EQ_CUTTER && held && p->matter>0){
        float step=fminf(dt,p->matter/eqGuns[gun].cost),length=eqGuns[gun].range;VoxelHit hit;bool wall=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(wall)length=hit.t;
        int target=eq_actor_ray(i,ray,&length);e->beamEnd=v_add(p->pos,v_mul(ray.direction,length));
        ducks_ray(ray,length,.08f,-1,false);
        if(target>=0)eq_damage(i,target,eqGuns[gun].damage*step,true,DEATH_CUTTER,p->pos);
        bool shieldContact=target<0&&wall&&go_shield_hit(hit.id,eqGuns[gun].damage*step);
        e->cutTime+=step;bool cut=false;
        for(int cuts=0;!shieldContact && cuts<24 && e->cutTime>=eqCutterUnitSeconds;cuts++){
            e->cutTime-=eqCutterUnitSeconds;float stop=eqGuns[gun].range;VoxelHit surface;
            bool blocked=first_voxel_hit_detailed(ray,stop,-1,false,&surface);if(blocked)stop=surface.t;
            if(eq_actor_ray(i,ray,&stop)<0 && blocked)cut|=eq_cut(i,ray,eqGuns[gun].range);
        }
        if(shieldContact)e->cutTime=0;
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
        ex_input(i,&c);
        eq_build_sample(i,i,&c);
        aim_sample(netTransport.role==NET_ROLE_OFFLINE?i:net_local_index_for_slot(i),i,dt,&c);
        if(debugConsoleOpen || !IsWindowFocused()){c.held=c.pressed=0;eq_cancel_charge(i);}
        eq_command(i,&c,dt);
    }
}
static void eq_tick(float dt){
    if(dt<=0)return;
    fp_tick(dt);
    ex_tick(dt);
    aim_tick(dt);
    eq_sword_tick(dt);
    for(int i=0;i<MAX_COMBATANTS;i++){Equipment *e=&equipment[i];for(int g=0;g<EQ_GUN_COUNT;g++)e->cooldown[g]=fmaxf(0,e->cooldown[g]-dt);e->switchTime=fmaxf(0,e->switchTime-dt);e->buildTime=fmaxf(0,e->buildTime-dt);
        if(netTransport.role!=NET_ROLE_CLIENT && !eq_usable(i))eq_cancel_charge(i);
        if(e->initialized && (!is_player_bot(i)||fp_carrier(i)))players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;
        if(gameState!=GAME_STATE_PLAYING || ff_input_blocked(i) || players[i].respawn_timer>0)e->buildHeld=false;
        float empty=players[i].matter<=0 && players[i].respawn_timer<=0?1:0;
        if(empty && !e->wasEmpty)e->depletionPulse=.45f;
        e->wasEmpty=empty!=0;e->depletionPulse=empty?fmaxf(0,e->depletionPulse-dt):0;
        e->depleted+=(empty-e->depleted)*fminf(1,dt*10);
        weaponVisuals[i].depleted=eq_gun(i)==EQ_RAIL && e->charge>0?0:e->depleted;
        weaponVisuals[i].charge=e->charge/.65f;weaponVisuals[i].activity=e->held||e->buildTime>0;
        float duration=eqGuns[eq_gun(i)].cooldown;
        weaponVisuals[i].cooldown=duration>0?Clamp(e->cooldown[eq_gun(i)]/duration,0,1):0;
        float deploy=e->buildHeld || e->buildTime>.35f?1:0;
        float step=dt/.15f;weaponVisuals[i].deployment+=Clamp(deploy-weaponVisuals[i].deployment,-step,step);
        weaponVisuals[i].build_pulse=e->buildTime/.6f;
        if(players[i].respawn_timer>0){e->held=false;e->charge=0;}
    }
    for(int n=0;n<32;n++)eqTraces[n].life=fmaxf(0,eqTraces[n].life-dt);
    if(netTransport.role==NET_ROLE_CLIENT){for(int n=0;n<32;n++)if(eqRockets[n].active)eqRockets[n].pos=world_wrap(v_add(eqRockets[n].pos,v_mul(eqRockets[n].vel,dt)));return;}
    for(int n=0;n<32;n++){EquipmentRocket *r=&eqRockets[n];if(!r->active)continue;float remaining=fminf(dt,.25f);
        while(remaining>0 && r->active){float step=fminf(remaining,1.f/120);remaining-=step;float length=v_length(r->vel)*step;Ray ray={r->pos,v_norm(r->vel)};VoxelHit hit;
            bool collision=go_combat_actor(r->owner)?go_world_ray(ray,length,&hit):first_voxel_hit_detailed(ray,length,-1,false,&hit);if(collision)length=hit.t;int actor=eq_actor_ray(r->owner,ray,&length);
            bool duck=ducks_ray(ray,length,.12f,-1,false);
            r->pos=world_wrap(v_add(r->pos,v_mul(ray.direction,length)));r->life-=step;
            if(collision || actor>=0 || duck || r->pos.y<.1f){EquipmentRocket copy=*r;r->active=false;eq_rocket_blast(copy);}else if(r->life<=0)r->active=false;
        }
    }
}
static Shader eqSlashShader;
static bool eqSlashTried;
#include "sword_wave.inc"
static int eq_build_side(int i){return eq_upgrade(i,EQ_BUILD)?5:3;}
static int eq_build_count(int i){int side=eq_build_side(i);return side*side*side;}
static bool eq_panel_positions(int i,Vector3 out[EQ_PANEL_CUBES]){
    Ray ray={players[i].pos,player_forward(&players[i])};VoxelHit hit;Vector3 at,normal={0,1,0};
    int side=eq_build_side(i),half=side/2;
    if(first_voxel_hit_detailed(ray,TETHER_RANGE,-1,false,&hit)){
        normal=hit.normal;at=(Vector3){(hit.x+.5f)*VOXEL_SIZE,(hit.y+.5f)*VOXEL_SIZE,(hit.z+.5f)*VOXEL_SIZE};
        at=v_add(at,v_mul(normal,(half+1)*VOXEL_SIZE));
    }else {if(ray.direction.y>=-.001f)return false;float t=-ray.position.y/ray.direction.y;if(t<0 || t>TETHER_RANGE)return false;
        at=v_add(ray.position,v_mul(ray.direction,t));at.y=(half+.5f)*VOXEL_SIZE;}
    at=world_wrap(at);at=(Vector3){(floorf(at.x/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.y/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.z/VOXEL_SIZE)+.5f)*VOXEL_SIZE};
    bool valid=true;int count=eq_build_count(i);
    for(int n=0;n<count;n++){
        out[n]=world_wrap(v_add(at,(Vector3){(n%side-half)*VOXEL_SIZE,((n/side)%side-half)*VOXEL_SIZE,(n/(side*side)-half)*VOXEL_SIZE}));
        int x=floorf(out[n].x/VOXEL_SIZE),y=floorf(out[n].y/VOXEL_SIZE),z=floorf(out[n].z/VOXEL_SIZE);
        int existing=table_get(x,y,z);
        bool blocked=existing>=0 && !voxel_is_fluid(&voxels[existing]) && enemy_equipment_owner(existing)<0 && get_goliath_owner_of_voxel(existing)<0;
        if(y<0 || out[n].y>FLOOR_SIZE*2 || blocked || (!worldPeriodic && (fabsf(out[n].x)>FLOOR_SIZE || fabsf(out[n].z)>FLOOR_SIZE)))valid=false;
        Vector3 extent={VOXEL_SIZE*.499f,VOXEL_SIZE*.499f,VOXEL_SIZE*.499f};BoundingBox cube={v_sub(out[n],extent),v_add(out[n],extent)};
        for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && !voxels[v].isBullet){VoxelWorldBounds b;voxel_particle_world_bounds(&voxels[v],&b);Vector3 center={(b.minx+b.maxx)*.5f,(b.miny+b.maxy)*.5f,(b.minz+b.maxz)*.5f};Vector3 shift=v_sub(world_image(center,out[n]),center);if(CheckCollisionBoxes(cube,(BoundingBox){v_add((Vector3){b.minx,b.miny,b.minz},shift),v_add((Vector3){b.maxx,b.maxy,b.maxz},shift)}) && get_goliath_owner_of_voxel(v)<0)valid=false;}
    }return valid;
}
static bool eq_build_ready(int i){int count=eq_build_count(i);return !players[i].isExposed && players[i].matter>=ff_build_cost(i) && GetTime()-players[i].last_build_time>=BUILD_COOLDOWN_SECONDS && voxel_count+count<MAX_VOXELS && MAX_PARTICLES-particle_pool_count+free_particle_count>=8*count;}
static bool eq_build_lifts(int count,const Vector3 *at,float lifts[MAX_COMBATANTS]){
    float roof=0;for(int n=0;n<count;n++)roof=fmaxf(roof,at[n].y+VOXEL_SIZE*.5f);
    memset(lifts,0,sizeof(float)*MAX_COMBATANTS);
    for(int j=0;j<activePlayers;j++)if(players[j].respawn_timer<=0 || combatantLife[j]==LIFE_DOWNED){
        BoundingBox body=combatant_bounds_at(j,players[j].pos,player_collision_half(&players[j]));bool overlap=false;
        for(int n=0;n<count&&!overlap;n++){Vector3 center=world_image(at[n],players[j].pos),extent={VOXEL_SIZE*.499f,VOXEL_SIZE*.499f,VOXEL_SIZE*.499f};overlap=CheckCollisionBoxes(body,(BoundingBox){v_sub(center,extent),v_add(center,extent)});}
        if(overlap){
            lifts[j]=fmaxf(0,roof-body.min.y+.02f);
            Vector3 raised=v_add(players[j].pos,(Vector3){0,lifts[j],0});
            if(player_aabb_hits_world_ex(raised,player_collision_half(&players[j]),j))return false;
        }
    }return true;
}
static bool eq_build(int i){
    if(!eq_usable(i) || !eq_build_ready(i))return false;Player *p=&players[i];Vector3 at[EQ_PANEL_CUBES];
    if(!eq_panel_positions(i,at))return false;
    float lifts[MAX_COMBATANTS];if(!eq_build_lifts(eq_build_count(i),at,lifts))return false;
    int count=eq_build_count(i);uint64_t ids[EQ_PANEL_CUBES]={0};int made=0;
    for(int n=0;n<count;n++){int v=add_static_voxel_at_grid(floorf(at[n].x/VOXEL_SIZE),floorf(at[n].y/VOXEL_SIZE),floorf(at[n].z/VOXEL_SIZE),player_palette_color(i),0);if(v<0)break;ids[made++]=voxels[v].identity;voxels[v].reinforced=eq_upgrade(i,EQ_BUILD);tag_owned_static_voxel(v,i);mark_surface(v);mark_surface_neighbors(at[n]);}
    if(made!=count){for(int n=0;n<made;n++){int v=enemy_voxel(ids[n]);if(v>=0)remove_voxel_index(v);}return false;}
    for(int j=0;j<activePlayers;j++)if(lifts[j]>0){players[j].pos.y+=lifts[j];players[j].vel.y=fmaxf(0,players[j].vel.y);players[j].onGround=false;}
    eq_spend(i,ff_build_cost(i));p->last_build_time=GetTime();equipment[i].buildTime=.6f;play_sfx(SFX_SHIELD);return true;
}
static void eq_build_command(int i,const NetInputCommand *c){
    if(i<0 || i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];
    if((c->pressed&NET_INPUT_BUILD_CANCEL) || gameState!=GAME_STATE_PLAYING || creativeModeActive || ff_input_blocked(i) || players[i].respawn_timer>0){e->buildHeld=false;return;}
    if(c->pressed&NET_INPUT_BUILD){e->buildHeld=true;aim_cancel(i);}
    if(e->buildHeld && !(c->held&NET_INPUT_BUILD)){e->buildHeld=false;if(netTransport.role!=NET_ROLE_CLIENT)eq_build(i);}
}
static void eq_build_sample(int seat,int slot,NetInputCommand *c){
    bool pad=playerInput[seat]==INPUT_TYPE_GAMEPAD;int device=seat_device(seat);
    if(debugConsoleOpen || !IsWindowFocused() || gameState!=GAME_STATE_PLAYING || (pad && !IsGamepadAvailable(device))){c->pressed|=NET_INPUT_BUILD_CANCEL;equipment[slot].buildHeld=false;return;}
    bool down=pad?IsGamepadButtonDown(device,GAMEPAD_BUTTON_RIGHT_TRIGGER_1):IsKeyDown(seat==0?KEY_E:KEY_O);
    bool press=pad?IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_TRIGGER_1):IsKeyPressed(seat==0?KEY_E:KEY_O);
    if(down)c->held|=NET_INPUT_BUILD;if(press)c->pressed|=NET_INPUT_BUILD;
    if(netTransport.role==NET_ROLE_CLIENT)eq_build_command(slot,c);
}

static Shader eqBeamShader;
static bool eqBeamTried;
static void eq_beam_shutdown(void){if(eqSlashShader.id && eqSlashShader.id!=rlGetShaderIdDefault())UnloadShader(eqSlashShader);eqSlashShader=(Shader){0};eqSlashTried=false;if(eqBeamShader.id && eqBeamShader.id!=rlGetShaderIdDefault())UnloadShader(eqBeamShader);eqBeamShader=(Shader){0};eqBeamTried=false;}
static void eq_world(void){
    if(!eqBeamTried){eqBeamTried=true;eqBeamShader=LoadShader(NULL,"shaders/equipment_beam.frag");}
    bool shader=eqBeamShader.id && eqBeamShader.id!=rlGetShaderIdDefault();
    float time=weaponFrameSampled?weaponFrameTime:(float)GetTime();
    if(shader){SetShaderValue(eqBeamShader,GetShaderLocation(eqBeamShader,"effectTime"),&time,SHADER_UNIFORM_FLOAT);BeginShaderMode(eqBeamShader);}

    for(int n=0;n<32;n++){EquipmentTrace *t=&eqTraces[n];if(t->life>0){
            Vector3 start=t->start;if(t->owner==drawViewPlayerIndex)start=player_weapon_pose(t->owner,true).muzzle;
            float width=t->gun==EQ_RAIL?.065f:.008f;
            DrawCylinderEx(start,t->end,width,width,8,Fade(t->gun==EQ_RAIL?PURPLE:t->gun==EQ_DMR?SKYBLUE:t->gun==EQ_GLOCK?LIGHTGRAY:ORANGE,t->life/.22f));
            if(t->gun==EQ_RAIL){DrawSphere(start,.10f*t->life/.22f,WHITE);DrawCylinderEx(start,t->end,.018f,.018f,6,Fade(WHITE,t->life/.22f));}
        }
        EquipmentRocket *r=&eqRockets[n];if(r->active){
            Vector3 axis=v_norm(r->vel),tail=v_sub(r->pos,v_mul(axis,.35f));
            DrawCylinderEx(tail,v_sub(tail,v_mul(axis,.7f)),.08f,.012f,5,ORANGE);
            Vector3 right=v_norm(v_cross(axis,fabsf(axis.y)>.9f?(Vector3){1,0,0}:(Vector3){0,1,0})),up=v_cross(axis,right);
            for(int k=0;k<4;k++){float age=fmodf(time*2+k*.25f+n*.13f,1),angle=k*2.4f+n;
                Vector3 point=v_add(v_sub(tail,v_mul(axis,age*1.7f)),v_add(v_mul(right,cosf(angle)*age*.18f),v_mul(up,sinf(angle)*age*.18f)));
                float size=.055f*(1-age);DrawCube(point,size,size,size,Fade(ORANGE,1-age));}
        }}
    for(int i=0;i<MAX_COMBATANTS;i++){Equipment *e=&equipment[i];if(e->initialized && e->held && eq_gun(i)==EQ_CUTTER && players[i].respawn_timer<=0){WeaponPose pose=player_weapon_pose(i,i==drawViewPlayerIndex);DrawCylinderEx(pose.muzzle,e->beamEnd,.075f,.06f,8,Fade(LIME,.40f));DrawCylinderEx(pose.muzzle,e->beamEnd,.026f,.025f,6,(Color){220,255,185,255});
            /* Full-range misses have no contact effect. A shorter actor/voxel ray does. */
            if(v_length(v_sub(e->beamEnd,players[i].pos))<eqGuns[EQ_CUTTER].range-.01f){
                DrawSphere(e->beamEnd,.085f,(Color){210,255,160,255});
                for(int k=0;k<6;k++){float age=fmodf(time*3+k/6.f,1),angle=k*2.4f;
                    Vector3 offset={cosf(angle)*age*.34f,age*.30f-age*age*.4f,sinf(angle)*age*.34f};
                    Vector3 point=v_add(e->beamEnd,offset);DrawLine3D(point,v_add(point,v_mul(offset,.25f)),Fade(LIME,1-age));}
            }}
        if(e->charge>0 && eq_gun(i)==EQ_RAIL){Vector3 end=v_add(players[i].pos,v_mul(player_forward(&players[i]),100));DrawCylinderEx(player_weapon_pose(i,i==drawViewPlayerIndex).muzzle,end,.006f,.006f,4,Fade(PURPLE,.5f));}
        if(!equipmentEmission && i==drawViewPlayerIndex && e->buildHeld && players[i].respawn_timer<=0){Vector3 pos[EQ_PANEL_CUBES]={{0}};float lifts[MAX_COMBATANTS];bool valid=eq_panel_positions(i,pos)&&eq_build_ready(i)&&eq_build_lifts(eq_build_count(i),pos,lifts);if(v_length(pos[1])>.01f)for(int n=0;n<eq_build_count(i);n++)DrawCubeWires(pos[n],VOXEL_SIZE,VOXEL_SIZE,VOXEL_SIZE,valid?player_palette_color(i):RED);}
    }
    if(shader)EndShaderMode();
    Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};
    for(int n=0;n<32;n++)if(eqRockets[n].active && periodic_visible(eqRockets[n].pos,2))weapon_draw_rocket(eqRockets[n].pos,eqRockets[n].vel,eye,time,equipmentEmission?WEAPON_EMISSION:WEAPON_COLOR);
    fp_draw(equipmentEmission?WEAPON_EMISSION:WEAPON_COLOR);
    eq_sword_draw();ex_draw(equipmentEmission?WEAPON_EMISSION:WEAPON_COLOR);
}
static void eq_icon(int gun,int x,int y,Color c){
    if(gun==EQ_BLADE){DrawLineEx((Vector2){x-5,y+22},(Vector2){x+22,y-2},4,c);DrawLine(x-3,y+9,x+10,y+22,c);return;}
    if(gun==EQ_GLOB){DrawCircle(x+6,y+10,7,c);DrawRectangleLines(x+17,y+2,10,14,c);}
    if(gun==EQ_DMR){DrawRectangleLines(x+7,y-2,13,5,c);DrawLine(x-12,y+9,x+30,y+9,c);}
    if(gun==EQ_GLOCK){DrawRectangle(x,y+5,24,6,c);DrawRectangle(x+16,y+11,6,12,c);return;}
    DrawRectangle(x,y+5,25,9,c);DrawRectangle(x+17,y+13,5,10,c);
    if(gun==EQ_SHOTGUN){for(int n=0;n<3;n++)DrawRectangle(x-4,y+1+n*6,5,4,c);}
    if(gun==EQ_RAIL||gun==EQ_CUTTER){DrawRectangle(x-12,y+3,18,3,c);DrawRectangle(x-12,y+12,18,3,c);}
    if(gun==EQ_ROCKET)DrawRectangleLines(x-6,y,34,19,c);
    if(gun==EQ_VOXEL)DrawRectangleLines(x+22,y-2,12,17,c);
}
static void eq_hud(int i,int w,int h){
    fp_hud(i,w,h);
    if(i<0||i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];if(!e->initialized)return;int sz=w<700?12:16;
    ex_hud(i,w,h);
    int radius=hud_radar_radius(h);
    int radarCenter=h-HUD_PADDING_Y-radius-HUD_BAR_THICKNESS-6;
    radarCenter=radarCenter-radius<HUD_PADDING_Y+HUD_BAR_HEIGHT?HUD_PADDING_Y+HUD_BAR_HEIGHT+radius+4:radarCenter;
    int y=radarCenter-radius-8-46;
    for(int s=0;s<2;s++){int x=14+s*105,g=e->gun[s];Color c=s==e->selected?GOLD:GRAY;DrawRectangle(x,y,98,46,Fade(BLACK,.7f));DrawRectangleLines(x,y,98,46,c);if(g!=EQ_EMPTY){eq_icon(g,x+20,y+2,c);DrawText(eqGuns[g].name,x+5,y+27,sz,c);}else DrawText("Empty",x+5,y+27,sz,GRAY);}
    int g=eq_gun(i);DrawText(TextFormat("%.0f Matter%s  %.1fs",eqGuns[g].cost,g==EQ_CUTTER?"/s":"",e->cooldown[g]),14,y-19,sz,RAYWHITE);
    if(g==EQ_RAIL && e->charge>0 && aimStates[i].blend<.02f)DrawRectangle(w/2-40,h/2+35,(int)(80*e->charge/.65f),5,PURPLE);
    DrawText(TextFormat("%s %s %s","",e->upgrades&EQ_ADHESIVE?"STICKY":"",e->upgrades&EQ_BUILD?"5x5x5":""),230,y+27,sz,GOLD);
    int n=eq_near_pickup(i);if(n>=0 && eq_pickup_gun(pickups[n].type)>=0 && e->gun[0]!=eq_pickup_gun(pickups[n].type) && e->gun[1]!=eq_pickup_gun(pickups[n].type)){const char *key=playerInput[i]==INPUT_TYPE_GAMEPAD?"Hold X":i==0?"B":"Y";const char *label=TextFormat("%s: equip %s",key,pickup_type_label(pickups[n].type));DrawText(label,(w-MeasureText(label,sz))/2,h/2+55,sz,WHITE);}
}
static void eq_layout(void){
    if(fp_active())return;
    if(useCustomMap||firefightPractice)return;bool foundry=gameMode==GAME_MODE_FIREFIGHT&&firefightFoundry;if(!foundry&&!clarity_crosscut())return;
    const int types[]={PICKUP_SHOTGUN,PICKUP_CUTTER,PICKUP_RAILGUN,PICKUP_ROCKET,PICKUP_SWORD,PICKUP_ADHESIVE,PICKUP_REINFORCED};
    const Vector3 anchors[]={{20,0,14},{28,0,14},{6,0,3},{3,0,6},{20,0,25},{28,0,25},{32,0,7}};
    for(int t=0;t<7;t++){Vector3 at;if(!clarity_supply_position(anchors[t],0,&at))continue;for(int q=0;q<(foundry?2:4);q++){int free=-1;for(int n=0;n<MAX_PICKUPS;n++)if(!pickups[n].active && pickups[n].respawnTimer<=0){free=n;break;}if(free<0)return;Vector3 pos=crosscut_rotate(at,foundry?q*2:q);Vector3 clear;if(!clarity_supply_position(pos,0,&clear))continue;pickups[free]=(Pickup){.pos=clear,.active=true,.type=types[t]};}}
    ex_layout();
}
#include "equipment_net.inc"
#include "equipment_expansion_net.inc"

#include "aiming.inc"
