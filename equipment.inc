/* Authoritative equipment behavior. All simulation runs once, outside viewports. */
static void eq_reset(int i){if(i<0 || i>=MAX_PLAYERS)return;equipment[i]=(Equipment){.gun={EQ_DEFAULT,EQ_EMPTY},.initialized=true};players[i].dynamicShotActive=false;}
static void eq_reset_all(void){for(int i=0;i<MAX_PLAYERS;i++)eq_reset(i);memset(eqRockets,0,sizeof(eqRockets));memset(eqTraces,0,sizeof(eqTraces));eqReceiveTick=0;memset(eqPartTick,0,sizeof(eqPartTick));}
static int eq_gun(int i){if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return EQ_DEFAULT;if(!equipment[i].initialized)eq_reset(i);return equipment[i].gun[equipment[i].selected];}
static bool eq_upgrade(int i,int flag){return i>=0 && i<MAX_PLAYERS && !is_player_bot(i) && (equipment[i].upgrades&flag)!=0;}
static bool eq_usable(int i){return i>=0 && i<MAX_PLAYERS && !is_player_bot(i) && netTransport.role!=NET_ROLE_CLIENT && gameState==GAME_STATE_PLAYING && players[i].respawn_timer<=0 && !ff_input_blocked(i) && !creativeModeActive && !(enemyLab&&enemyLabPaused);}
static void eq_spend(int i,float cost){players[i].matter=fmaxf(0,players[i].matter-cost);players[i].isExposed=players[i].matter<=0;}
static int eq_pickup_gun(int type){switch(type){case PICKUP_DYNAMIC_SHOT:return EQ_VOXEL;case PICKUP_SHOTGUN:return EQ_SHOTGUN;case PICKUP_CUTTER:return EQ_CUTTER;case PICKUP_RAILGUN:return EQ_RAIL;case PICKUP_ROCKET:return EQ_ROCKET;default:return -1;}}
static int eq_pickup_kind(int type){int gun=eq_pickup_gun(type);if(gun>=0)return eqGuns[gun].visual;return type==PICKUP_SWORD?WEAPON_SWORD:type==PICKUP_ADHESIVE?WEAPON_ADHESIVE_GRENADE:type==PICKUP_REINFORCED?WEAPON_REINFORCED_BUILDER:-1;}
static bool eq_pickup(int i,int type,bool replace){
    if(i<0 || i>=MAX_PLAYERS || is_player_bot(i))return false;
    eq_gun(i);Equipment *e=&equipment[i];int gun=eq_pickup_gun(type);
    if(gun>=0){
        if(e->gun[0]==gun || e->gun[1]==gun)return false;
        int slot=e->gun[0]==EQ_EMPTY?0:e->gun[1]==EQ_EMPTY?1:replace?e->selected:-1;
        if(slot<0)return false;e->gun[slot]=gun;e->selected=slot;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;
        players[i].dynamicShotActive=gun==EQ_VOXEL;return true;
    }
    int flag=type==PICKUP_SWORD?EQ_SWORD:type==PICKUP_ADHESIVE?EQ_ADHESIVE:type==PICKUP_REINFORCED?EQ_BUILD:0;
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
static void eq_trace(Vector3 a,Vector3 b,int gun){int slot=0;for(int n=0;n<32;n++)if(eqTraces[n].life<=0){slot=n;break;}eqTraces[slot]=(EquipmentTrace){a,b,.22f,gun};}
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
static void eq_impulse(Vector3 at,Vector3 forward,float radius,float strength){
    for(int n=0;n<voxel_count;n++){Voxel *v=&voxels[n];if(!v->simulate || v->isBullet || v_length(world_delta(v->pos,at))>radius)continue;
        Vector3 d=world_delta(v->pos,at);if(v_dot(v_norm(d),forward)<.75f)continue;
        if(v->sleeping)wake_sleeping_cluster_in_place(n);
        Vector3 impulse=v_mul(forward,strength*(1-v_length(d)/radius));v->vel=v_add(v->vel,impulse);
        for(int p=0;p<8;p++)v->particles[p]->vel=v_add(v->particles[p]->vel,impulse);
    }
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
    if(e->switchTime>0 || e->cooldown[gun]>0 || p->isExposed || p->matter<def->cost)return false;
    Vector3 f=player_forward(p);Ray ray={p->pos,f};
    if(gun==EQ_SHOTGUN){
        float damage[MAX_COMBATANTS]={0};Vector3 r=v_norm(v_cross(f,(Vector3){0,1,0})),up=v_cross(r,f);
        for(int n=0;n<EQ_PELLETS;n++){float angle=n*2*PI/EQ_PELLETS;Vector3 d=v_norm(v_add(f,v_add(v_mul(r,cosf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)),v_mul(up,sinf(angle)*tanf(eqShotgunCone*.5f*DEG2RAD)))));Ray pellet={p->pos,d};float length=def->range;VoxelHit hit;
            if(first_voxel_hit_detailed(pellet,length,-1,false,&hit))length=hit.t;
            int target=eq_actor_ray(i,pellet,&length);if(target>=0)damage[target]+=def->damage;
            eq_trace(p->pos,v_add(p->pos,v_mul(d,length)),gun);
        }
        for(int j=0;j<activePlayers;j++)if(damage[j]>0)eq_damage(i,j,damage[j],false,DEATH_SHOTGUN,p->pos);
        eq_impulse(p->pos,f,8,18);
    }else {
        int free=-1;for(int n=0;n<32;n++)if(!eqRockets[n].active){free=n;break;}if(free<0)return false;
        // Start at the eye and sweep from it: no spawning through nearby walls.
        eqRockets[free]=(EquipmentRocket){true,i,p->pos,v_mul(f,eqRocketSpeed),eqRocketLife};
    }
    eq_spend(i,def->cost);e->cooldown[gun]=def->cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,gun==EQ_ROCKET?.55f:.8f);ai_gunfire(i,p->pos);return true;
}
static void eq_command(int i,const NetInputCommand *c,float dt){
    if(i<0 || i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];eq_gun(i);
    if(!eq_usable(i)){e->held=false;e->charge=0;return;}
    if(c->pressed&NET_INPUT_SWITCH){int other=1-e->selected;if(e->gun[other]!=EQ_EMPTY){e->selected=other;e->switchTime=eqSwitchSeconds;e->charge=0;e->held=false;players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;}}
    if(c->pressed&NET_INPUT_EQUIP){int n=eq_near_pickup(i);if(n>=0 && eq_pickup(i,pickups[n].type,true)){pickups[n].active=false;pickups[n].respawnTimer=pickups[n].type==PICKUP_RAILGUN||pickups[n].type==PICKUP_ROCKET?60:30;}}
    bool held=(c->held&NET_INPUT_FIRE)!=0;int gun=eq_gun(i);Player *p=&players[i];
    if(e->switchTime>0 || p->isExposed || p->meleeSwingActive || e->cooldown[gun]>0){e->held=false;e->charge=0;return;}
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
            eq_trace(p->pos,v_add(p->pos,v_mul(ray.direction,100)),EQ_RAIL);eq_spend(i,eqGuns[gun].cost);e->cooldown[gun]=eqGuns[gun].cooldown;weapon_visual_shot(&weaponVisuals[i],(float)GetTime());play_sfx_at(SFX_FIRE,p->pos,i,2,.65f);ai_gunfire(i,p->pos);
        }
        if(!held)e->charge=0;
    }else if(gun==EQ_CUTTER && held && p->matter>0){
        float step=fminf(dt,p->matter/eqGuns[gun].cost),length=eqGuns[gun].range;VoxelHit hit;bool wall=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(wall)length=hit.t;
        int target=eq_actor_ray(i,ray,&length);e->beamEnd=v_add(p->pos,v_mul(ray.direction,length));
        if(target>=0)eq_damage(i,target,eqGuns[gun].damage*step,false,DEATH_CUTTER,p->pos);
        e->cutTime+=step;if(e->cutTime>=eqCutterUnitSeconds){e->cutTime=fmodf(e->cutTime,eqCutterUnitSeconds);if(target<0 && wall)eq_cut(i,ray,12);play_sfx_at(SFX_FIRE,p->pos,i,1,1.6f);ai_gunfire(i,p->pos);}
        eq_spend(i,eqGuns[gun].cost*step);
    }else e->cutTime=0;
    e->held=held;
}
static void eq_local(float dt){
    if(netTransport.role==NET_ROLE_CLIENT)return;
    for(int i=0;i<MAX_PLAYERS && i<activePlayers;i++)if(!is_player_bot(i) && (netTransport.role==NET_ROLE_OFFLINE || net_player_is_local(i))){
        NetInputCommand c={0};bool pad=playerInput[i]==INPUT_TYPE_GAMEPAD;int d=seat_device(i);
        if(pad){if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_TRIGGER_2))c.pressed|=NET_INPUT_FIRE;if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_LEFT))c.pressed|=NET_INPUT_SWITCH;if(IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_RIGHT))c.pressed|=NET_INPUT_EQUIP;if(IsGamepadButtonDown(d,GAMEPAD_BUTTON_RIGHT_TRIGGER_2))c.held|=NET_INPUT_FIRE;}
        else {if(!mouse_owns(i) && IsKeyPressed(i==0?KEY_LEFT_CONTROL:KEY_RIGHT_CONTROL))c.pressed|=NET_INPUT_FIRE;if(IsKeyPressed(i==0?KEY_Q:KEY_U))c.pressed|=NET_INPUT_SWITCH;if(IsKeyPressed(i==0?KEY_B:KEY_Y))c.pressed|=NET_INPUT_EQUIP;if(IsKeyDown(i==0?KEY_LEFT_CONTROL:KEY_RIGHT_CONTROL))c.held|=NET_INPUT_FIRE;}
        if(mouse_owns(i)){if(mouseCaptured&&mouseButtonsArmed&&IsMouseButtonPressed(MOUSE_BUTTON_LEFT))c.pressed|=NET_INPUT_FIRE;if(mouseCaptured&&mouseButtonsArmed&&IsMouseButtonDown(MOUSE_BUTTON_LEFT))c.held|=NET_INPUT_FIRE;if(mouseCaptured&&GetMouseWheelMove()!=0)c.pressed|=NET_INPUT_SWITCH;}
        if(debugConsoleOpen || !IsWindowFocused()){c.held=c.pressed=0;equipment[i].charge=0;equipment[i].held=false;}
        eq_command(i,&c,dt);
    }
}
static void eq_tick(float dt){
    if(dt<=0)return;
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
    int cuts=0;for(int n=0;n<EQ_SWORD_CUTS;n++){float a=(-30+n*20)*DEG2RAD;Vector3 dir=Vector3RotateByAxisAngle(center.direction,(Vector3){0,1,0},a);if(eq_cut(i,(Ray){center.position,dir},eqSwordReach))cuts++;}
    if(cuts)add_player_matter(&players[i],ff_harvest_amount(i));play_sfx_at(SFX_MELEE,center.position,i,1,.7f);
}
static bool eq_panel_positions(int i,Vector3 out[3]){
    Ray ray={players[i].pos,player_forward(&players[i])};VoxelHit hit;Vector3 at;
    if(first_voxel_hit_detailed(ray,TETHER_RANGE,-1,false,&hit))at=v_add(ray.position,v_add(v_mul(ray.direction,hit.t),v_mul(hit.normal,VOXEL_SIZE*.51f)));
    else {if(ray.direction.y>=-.001f)return false;float t=-ray.position.y/ray.direction.y;if(t>TETHER_RANGE)return false;at=v_add(ray.position,v_mul(ray.direction,t));at.y=VOXEL_SIZE*.5f;}
    at=world_wrap(at);at=(Vector3){(floorf(at.x/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.y/VOXEL_SIZE)+.5f)*VOXEL_SIZE,(floorf(at.z/VOXEL_SIZE)+.5f)*VOXEL_SIZE};
    Vector3 side=fabsf(ray.direction.z)>=fabsf(ray.direction.x)?(Vector3){VOXEL_SIZE,0,0}:(Vector3){0,0,VOXEL_SIZE};bool valid=true;
    for(int n=0;n<3;n++){out[n]=world_wrap(v_add(at,v_mul(side,n-1)));int x=floorf(out[n].x/VOXEL_SIZE),y=floorf(out[n].y/VOXEL_SIZE),z=floorf(out[n].z/VOXEL_SIZE);
        if(y<0 || occupied(x,y,z) || (!worldPeriodic && (fabsf(out[n].x)>FLOOR_SIZE || fabsf(out[n].z)>FLOOR_SIZE)))valid=false;
        BoundingBox cube={v_sub(out[n],(Vector3){.25f,.25f,.25f}),v_add(out[n],(Vector3){.25f,.25f,.25f})};
        for(int v=0;v<voxel_count;v++)if(voxels[v].simulate && !voxels[v].isBullet){float radius=voxels[v].rest_edge*.5f;Vector3 center=world_image(voxels[v].pos,out[n]);if(CheckCollisionBoxes(cube,(BoundingBox){v_sub(center,(Vector3){radius,radius,radius}),v_add(center,(Vector3){radius,radius,radius})}))valid=false;}
        for(int j=0;j<activePlayers;j++)if((players[j].respawn_timer<=0 || combatantLife[j]==LIFE_DOWNED) && CheckCollisionBoxes(cube,combatant_bounds_at(j,world_image(players[j].pos,out[n]),PLAYER_SIZE*.5f)))valid=false;
    }return valid;
}
static bool eq_build(int i){
    if(!eq_usable(i))return false;Player *p=&players[i];float cost=3*ff_build_cost(i);Vector3 at[3];
    if(p->isExposed || p->matter<cost || GetTime()-p->last_build_time<BUILD_COOLDOWN_SECONDS || !eq_panel_positions(i,at) || voxel_count+3>=MAX_VOXELS || MAX_PARTICLES-particle_pool_count+free_particle_count<24)return false;
    uint64_t ids[3]={0};int made=0;
    for(int n=0;n<3;n++){int v=add_static_voxel_at_grid(floorf(at[n].x/VOXEL_SIZE),floorf(at[n].y/VOXEL_SIZE),floorf(at[n].z/VOXEL_SIZE),player_palette_color(i),0);if(v<0)break;ids[made++]=voxels[v].identity;voxels[v].reinforced=true;tag_owned_static_voxel(v,i);mark_surface(v);mark_surface_neighbors(at[n]);}
    if(made!=3){for(int n=0;n<made;n++){int v=enemy_voxel(ids[n]);if(v>=0)remove_voxel_index(v);}return false;}
    eq_spend(i,cost);p->last_build_time=GetTime();equipment[i].buildTime=.6f;play_sfx(SFX_SHIELD);return true;
}

static void eq_world(void){
    for(int n=0;n<32;n++){EquipmentTrace *t=&eqTraces[n];if(t->life>0)DrawCylinderEx(t->start,t->end,t->gun==EQ_RAIL?.025f:.008f,t->gun==EQ_RAIL?.025f:.008f,5,Fade(t->gun==EQ_RAIL?PURPLE:ORANGE,t->life/.22f));
        EquipmentRocket *r=&eqRockets[n];if(r->active){DrawCube(r->pos,.22f,.22f,.3f,ORANGE);DrawCylinderEx(r->pos,v_sub(r->pos,v_mul(v_norm(r->vel),.8f)),.1f,.02f,5,RED);}}
    for(int i=0;i<MAX_PLAYERS;i++){Equipment *e=&equipment[i];if(e->initialized && e->held && eq_gun(i)==EQ_CUTTER && players[i].respawn_timer<=0){WeaponPose pose=player_weapon_pose(i,i==drawViewPlayerIndex);DrawCylinderEx(pose.muzzle,e->beamEnd,.018f,.018f,5,LIME);DrawSphere(e->beamEnd,.06f,LIME);}
        if(e->charge>0 && eq_gun(i)==EQ_RAIL){Vector3 end=v_add(players[i].pos,v_mul(player_forward(&players[i]),100));DrawCylinderEx(player_weapon_pose(i,i==drawViewPlayerIndex).muzzle,end,.006f,.006f,4,Fade(PURPLE,.5f));}
        if(!equipmentEmission && i==drawViewPlayerIndex && eq_upgrade(i,EQ_BUILD) && players[i].respawn_timer<=0){Vector3 pos[3]={{0}};bool valid=eq_panel_positions(i,pos)&&players[i].matter>=3*ff_build_cost(i);if(v_length(pos[1])>.01f)for(int n=0;n<3;n++)DrawCubeWires(pos[n],VOXEL_SIZE,VOXEL_SIZE,VOXEL_SIZE,valid?GOLD:RED);}
    }
}
static void eq_icon(int gun,int x,int y,Color c){
    DrawRectangle(x,y+5,25,9,c);DrawRectangle(x+17,y+13,5,10,c);
    if(gun==EQ_SHOTGUN){for(int n=0;n<3;n++)DrawRectangle(x-4,y+1+n*6,5,4,c);}
    if(gun==EQ_RAIL||gun==EQ_CUTTER){DrawRectangle(x-12,y+3,18,3,c);DrawRectangle(x-12,y+12,18,3,c);}
    if(gun==EQ_ROCKET)DrawRectangleLines(x-6,y,34,19,c);
    if(gun==EQ_VOXEL)DrawRectangleLines(x+22,y-2,12,17,c);
}
static void eq_hud(int i,int w,int h){
    if(i<0||i>=MAX_PLAYERS)return;Equipment *e=&equipment[i];if(!e->initialized)return;int y=h-105,sz=w<700?12:16;
    for(int s=0;s<2;s++){int x=14+s*105,g=e->gun[s];Color c=s==e->selected?GOLD:GRAY;DrawRectangle(x,y,98,46,Fade(BLACK,.7f));DrawRectangleLines(x,y,98,46,c);if(g!=EQ_EMPTY){eq_icon(g,x+20,y+2,c);DrawText(eqGuns[g].name,x+5,y+27,sz,c);}else DrawText("Empty",x+5,y+27,sz,GRAY);}
    int g=eq_gun(i);DrawText(TextFormat("%.0f Matter%s  %.1fs",eqGuns[g].cost,g==EQ_CUTTER?"/s":"",e->cooldown[g]),14,y-19,sz,RAYWHITE);
    if(g==EQ_RAIL && e->charge>0)DrawRectangle(w/2-40,h/2+35,(int)(80*e->charge/.65f),5,PURPLE);
    DrawText(TextFormat("%s %s %s",e->upgrades&EQ_SWORD?"SWORD":"",e->upgrades&EQ_ADHESIVE?"STICKY":"",e->upgrades&EQ_BUILD?"PANEL":""),14,y-36,sz,GOLD);
    int n=eq_near_pickup(i);if(n>=0 && eq_pickup_gun(pickups[n].type)>=0 && e->gun[0]!=eq_pickup_gun(pickups[n].type) && e->gun[1]!=eq_pickup_gun(pickups[n].type)){const char *key=playerInput[i]==INPUT_TYPE_GAMEPAD?"D-pad Right":i==0?"B":"Y";const char *label=TextFormat("%s: equip %s",key,pickup_type_label(pickups[n].type));DrawText(label,(w-MeasureText(label,sz))/2,h/2+55,sz,WHITE);}
}
static void eq_layout(void){
    if(useCustomMap||firefightPractice)return;bool foundry=gameMode==GAME_MODE_FIREFIGHT&&firefightFoundry;if(!foundry&&!clarity_crosscut())return;
    const int types[]={PICKUP_SHOTGUN,PICKUP_CUTTER,PICKUP_RAILGUN,PICKUP_ROCKET,PICKUP_SWORD,PICKUP_ADHESIVE,PICKUP_REINFORCED};
    const Vector3 anchors[]={{20,0,14},{28,0,14},{6,0,3},{3,0,6},{20,0,25},{28,0,25},{32,0,7}};
    for(int t=0;t<7;t++){Vector3 at;if(!clarity_supply_position(anchors[t],0,&at))continue;for(int q=0;q<(foundry?2:4);q++){int free=-1;for(int n=0;n<MAX_PICKUPS;n++)if(!pickups[n].active && pickups[n].respawnTimer<=0){free=n;break;}if(free<0)return;Vector3 pos=crosscut_rotate(at,foundry?q*2:q);Vector3 clear;if(!clarity_supply_position(pos,0,&clear))continue;pickups[free]=(Pickup){.pos=clear,.active=true,.type=types[t]};}}
}
#include "equipment_net.inc"
