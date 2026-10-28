/* Shared builder/rifleman watchtowers. Geometry and routes use the same grid blueprint. */
typedef struct {int x,y,z;} TowerCell;
static TowerCell towerCells[ENEMY_STRUCTURE_CELLS];
static int towerCellCount;
static void tower_blueprint(void){
    if(towerCellCount)return;
    for(int z=0;z<18;z++)for(int y=0;y<10;y++)for(int x=0;x<6;x++){
        bool stair=z<12 && x>=1 && x<=3 && y<=z*8/12;
        bool deck=z>=12 && y==7;
        bool corner=z>=12 && (z==12 || z==17) && (x==0 || x==5);
        bool shelter=z>=12 && (y==4 || (y<4 && (x==0 || x==5 || z==12 || (z==17 && (x<2 || x>3)))));
        if(stair || deck || corner || shelter){if(towerCellCount>=ENEMY_STRUCTURE_CELLS)return;towerCells[towerCellCount++]=(TowerCell){x,y,z};}
    }
}
static Vector3 tower_point(const EnemyStructure *t,Vector3 local){return world_wrap(v_add(t->origin,crosscut_rotate(local,t->quarter)));}
static Vector3 tower_cell(const EnemyStructure *t,int k){TowerCell c=towerCells[k];return tower_point(t,(Vector3){(c.x+.5f)*.5f,(c.y+.5f)*.5f,(c.z+.5f)*.5f});}
static int tower_at(Vector3 p){return table_get_static_only((int)floorf(p.x/.5f),(int)floorf(p.y/.5f),(int)floorf(p.z/.5f));}
static EnemyStructure *tower_find(uint32_t tag){for(int b=0;b<ENEMY_STRUCTURE_MAX;b++)if(enemyStructures[b].tower && enemyStructures[b].tag==tag)return &enemyStructures[b];return NULL;}
static bool tower_owner(const EnemyStructure *t){return t->owner>=MAX_PLAYERS && t->owner<activePlayers && players[t->owner].respawn_timer<=0 && players[t->owner].enemyType==ENEMY_TYPE_BUILDER && enemySpecial[t->owner].generation==t->generation;}
static bool tower_eligible(int i){
    if(i<MAX_PLAYERS || i>=activePlayers || players[i].respawn_timer>0 || players[i].enemyType!=ENEMY_TYPE_STANDARD)return false;
    int gun=fpCarrierType[i]>0?eq_gun(i):EQ_DEFAULT;
    return gun==EQ_DEFAULT || gun==EQ_GLOB || gun==EQ_GLOCK || gun==EQ_DMR || gun==EQ_RAIL;
}
static void tower_release(int i){EnemyStructure *t=tower_find(towerRifles[i].tag);if(t && t->occupant==i)t->occupant=-1;memset(&towerRifles[i],0,sizeof(towerRifles[i]));towerRifles[i].retry=2;}
static bool tower_stationed(int i){return i>=0 && i<MAX_COMBATANTS && towerRifles[i].tag && towerRifles[i].state==2;}
static void tower_reset(void){
    /* Tags follow activated pieces even when their identity or array index changes. */
    for(int v=voxel_count-1;v>=0;v--)if(voxels[v].towerTag)remove_voxel_index(v);
    memset(towerRifles,0,sizeof(towerRifles));towerNetTick=0;
}
static bool tower_actor_blocked(Vector3 p){
    BoundingBox cell={v_sub(p,(Vector3){.249f,.249f,.249f}),v_add(p,(Vector3){.249f,.249f,.249f})};
    for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0 && combatantLife[i]!=LIFE_UNUSED){
        BoundingBox actor=combatant_bounds_at(i,world_image(players[i].pos,p),player_collision_half(&players[i]));
        if(CheckCollisionBoxes(cell,actor))return true;
    }return false;
}
static bool tower_site(const EnemyStructure *t){
    /* Actor clearance is checked once in local space, not once per empty cell. */
    for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0 && combatantLife[i]!=LIFE_UNUSED){
        Vector3 p=crosscut_rotate(world_delta(players[i].pos,t->origin),(4-t->quarter)%4);
        if(p.x>-1.5f && p.x<4.5f && p.z>-1.5f && p.z<10.5f && p.y<7)return false;
    }
    /* Entire walkable envelope, including stair headroom, must be clear. */
    for(int z=-2;z<20;z++)for(int x=-2;x<8;x++)for(int y=0;y<13;y++){
        Vector3 p=tower_point(t,(Vector3){(x+.5f)*.5f,(y+.5f)*.5f,(z+.5f)*.5f});
        if(fabsf(p.x)>FLOOR_SIZE-1 || fabsf(p.z)>FLOOR_SIZE-1 || tower_at(p)>=0)return false;

    }
    Vector3 center=tower_point(t,(Vector3){1.5f,0,4.5f});
    for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active && v_length(world_delta(center,pickups[n].pos))<6)return false;
    for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner==-1 && v_length(world_delta(center,fpItems[n].pos))<6)return false;
    for(int n=0;n<foundryTrapCount;n++){
        FoundryTrap *f=&foundryTraps[n];Vector3 c={(f->x0+f->x1+1)*.25f,0,(f->z0+f->z1+1)*.25f};
        Vector3 d=world_delta(center,c);if(fabsf(d.x)<(f->x1-f->x0+1)*.25f+5 && fabsf(d.z)<(f->z1-f->z0+1)*.25f+5)return false;
    }
    for(int n=0;n<ENEMY_STRUCTURE_MAX;n++){EnemyStructure *other=&enemyStructures[n];if(other==t || !other->tower)continue;if(v_length(world_delta(center,tower_point(other,(Vector3){1.5f,0,4.5f})))<11)return false;}
    if(gameMode==GAME_MODE_FIREFIGHT && firefightFoundry){
        if(fabsf(center.x)<13 && fabsf(center.z)<12)return false;
        if(fabsf(center.x)>FLOOR_SIZE-8 && fabsf(center.z)<9)return false;
    }
    return true;
}
static bool tower_place(EnemyStructure *t){
    tower_blueprint();int missing=0;
    for(int k=0;k<towerCellCount;k++){
        Vector3 p=tower_cell(t,k);int id=tower_at(p);
        if(id>=0 && voxels[id].towerTag==t->tag){t->ids[k]=voxels[id].identity;continue;}
        if(id>=0 || tower_actor_blocked(p))return false;
        int occupied=table_get((int)floorf(p.x/.5f),(int)floorf(p.y/.5f),(int)floorf(p.z/.5f));
        if(occupied>=0 && voxels[occupied].towerTag!=t->tag)return false;
        missing++;
    }
    if(!missing){t->built=true;return true;}
    if(voxel_count+missing>=MAX_VOXELS || particle_pool_count-free_particle_count+8*missing>=MAX_PARTICLES)return false;
    for(int i=0;i<activePlayers;i++)if(players[i].tetherHolding){int v=resolve_tether_voxel(&players[i])?players[i].tetherVoxel:-1;if(v>=0 && voxels[v].towerTag==t->tag){players[i].tetherHolding=players[i].goldTetherHolding=false;players[i].tetherVoxel=-1;players[i].tetherVoxelIdentity=0;}}
    for(int v=voxel_count-1;v>=0;v--)if(voxels[v].towerTag==t->tag && voxels[v].simulate)remove_voxel_index(v);
    uint64_t created[ENEMY_STRUCTURE_CELLS];int count=0;
    for(int k=0;k<towerCellCount;k++){
        Vector3 p=tower_cell(t,k);int id=tower_at(p);
        if(id<0){
            Color c=towerCells[k].z<12?(Color){100,151,165,255}:(Color){88,110,130,255};
            id=add_static_voxel_at_grid((int)floorf(p.x/.5f),(int)floorf(p.y/.5f),(int)floorf(p.z/.5f),c,0);
            if(id<0){for(int j=0;j<count;j++){int v=enemy_voxel(created[j]);if(v>=0)remove_voxel_index(v);}return false;}
            voxels[id].owner=ENEMY_CONSTRUCTION_OWNER;voxels[id].towerTag=t->tag;created[count++]=voxels[id].identity;
        }
        t->ids[k]=voxels[id].identity;
    }
    t->count=towerCellCount;t->built=true;t->pulse=.35f;collisionTopologyDirty=true;
    return true;
}
static bool tower_search(int i){
    int slot=-1;for(int b=0;b<ENEMY_STRUCTURE_MAX;b++)if(!enemyStructures[b].tower && !enemyStructures[b].count){slot=b;break;}
    if(slot<0)return false;
    Vector3 desired=players[i].pos;float nearest=16;
    for(int j=MAX_PLAYERS;j<activePlayers;j++)if(tower_eligible(j)){float d=v_length(world_delta(players[j].pos,players[i].pos));if(d<nearest){nearest=d;desired=players[j].pos;}}
    Vector3 dir=world_delta(desired,players[i].pos);float base=atan2f(dir.z,dir.x)+enemySpecial[i].buildQuarter*PI/4;
    int start=(int)(enemySpecial[i].selected%352);enemySpecial[i].selected=(start+16)%352;
    for(int attempt=0;attempt<16;attempt++){
        int candidateIndex=(start+attempt)%352,n=candidateIndex/4,q=candidateIndex%4;
        float a=base+(n%8)*PI/4,r=6+(n/8)*2;
        Vector3 center=v_add(players[i].pos,(Vector3){cosf(a)*r,0,sinf(a)*r});center.y=0;
        Vector3 origin=v_sub(center,crosscut_rotate((Vector3){1.5f,0,4.5f},q));origin.x=floorf(origin.x*2)*.5f;origin.z=floorf(origin.z*2)*.5f;
        EnemyStructure candidate={.tower=true,.owner=i,.generation=enemySpecial[i].generation,.origin=world_wrap(origin),.quarter=q,.occupant=-1};
        if(tower_site(&candidate)){
            candidate.bestDistance=FLT_MAX;candidate.tag=towerNextTag++;if(!candidate.tag)candidate.tag=towerNextTag++;
            enemyStructures[slot]=candidate;return true;
        }
    }return false;
}
static void tower_builder(int i,float dt){
    if(netTransport.role==NET_ROLE_CLIENT || players[i].respawn_timer>0)return;
    tower_blueprint();ai_tick(i,dt);EnemySpecialState *s=&enemySpecial[i];s->searchTimer-=dt;
    EnemyStructure *t=NULL;for(int b=0;b<ENEMY_STRUCTURE_MAX;b++)if(enemyStructures[b].tower && enemyStructures[b].owner==i && enemyStructures[b].generation==s->generation){t=&enemyStructures[b];break;}
    if(!t){
        s->reason="seeking tower site";players[i].vel.x=players[i].vel.z=0;
        if(s->searchTimer<=0){s->searchTimer=.25f;tower_search(i);}return;
    }
    Vector3 service=tower_point(t,(Vector3){-1,BASE_EYE_HEIGHT,6.5f});
    if(t->built){
        /* Approach the rear door around the outside, then enter the sheltered room. */
        const Vector3 route[3]={{-.65f,BASE_EYE_HEIGHT,9.5f},{1.5f,BASE_EYE_HEIGHT,9.5f},{1.5f,BASE_EYE_HEIGHT,7.5f}};
        service=tower_point(t,route[t->shelterStep]);
        Vector3 delta=world_delta(service,players[i].pos);delta.y=0;
        if(v_length(delta)<.2f && t->shelterStep<2)service=tower_point(t,route[++t->shelterStep]);
    }
    Vector3 to=world_delta(service,players[i].pos);to.y=0;
    float distance=v_length(to);s->buildOrigin=t->origin;s->buildKind=2;
    if(!t->built && distance>4){
        if(distance<t->bestDistance-.1f){t->bestDistance=distance;t->approach=0;}else t->approach+=dt;
        if(t->approach>5){memset(t,0,sizeof(*t));s->searchTimer=2;s->buildQuarter=(s->buildQuarter+1)%8;return;}
    }
    if(distance>(t->built?.15f:1)){bot_apply_move(&players[i],&botStates[i],to,fminf(MOVE_SPEED*.65f,distance/fmaxf(dt,.001f)),dt);s->reason=t->built?"entering tower shelter":"approaching watchtower";}
    else {players[i].vel.x=players[i].vel.z=0;s->reason=t->built?"maintaining watchtower":"constructing watchtower";}
    ai_aim(i,tower_point(t,(Vector3){1.5f,3,7.5f}),dt);
    Vector3 local=crosscut_rotate(world_delta(players[i].pos,t->origin),(4-t->quarter)%4);
    float dx=fmaxf(0,fmaxf(-local.x,local.x-3)),dz=fmaxf(0,fmaxf(-local.z,local.z-9));
    if(t->built?sqrtf(dx*dx+dz*dz)>6:distance>4)return;
    if(!t->built){
        t->construction=fminf(2,t->construction+dt);s->windup=2-t->construction;
        if(t->construction>=2){if(!tower_place(t)){s->reason="tower space blocked";t->blocked+=dt;if(t->blocked>5){memset(t,0,sizeof(*t));s->searchTimer=2;s->buildQuarter=(s->buildQuarter+1)%8;}}else {s->actions++;s->windup=0;}}
    }else if(t->repairClock<=0){t->repairClock=.25f;if(tower_place(t)&&t->pulse>0){s->repairs++;s->actions++;}}
}
static void tower_tick(float dt){
    for(int b=0;b<ENEMY_STRUCTURE_MAX;b++){
        EnemyStructure *t=&enemyStructures[b];if(!t->tower)continue;
        t->repairClock-=dt;t->pulse=fmaxf(0,t->pulse-dt);
        if(!tower_owner(t)&&!t->built){memset(t,0,sizeof(*t));continue;}
        if(t->occupant>=0 && (!tower_eligible(t->occupant)||towerRifles[t->occupant].tag!=t->tag)){if(t->occupant<MAX_COMBATANTS)tower_release(t->occupant);t->occupant=-1;}
    }
}
static Vector3 tower_waypoint(const EnemyStructure *t,int step){
    if(step==0)return tower_point(t,(Vector3){1.25f,BASE_EYE_HEIGHT,-.7f});
    if(step<=8){int s=step-1,start=(s*12+7)/8,end=((s+1)*12+7)/8;return tower_point(t,(Vector3){1.25f,step*.5f+.27f,(start+end)*.25f});}
    return tower_point(t,(Vector3){1.5f,4.27f,7.5f});
}
static bool tower_floor(Vector3 waypoint,uint32_t tag){Vector3 at=waypoint;at.y-=.3f;int id=tower_at(at);return id>=0 && voxels[id].towerTag==tag;}
static bool tower_move(int i,Vector3 goal,float dt,bool stairs){
    Player *p=&players[i];TowerRifle *r=&towerRifles[i];Vector3 to=world_delta(goal,p->pos);float height=to.y;to.y=0;float d=v_length(to);
    if(d<(stairs?.22f:.06f) && height<.2f)return true;
    float speed=fminf(MOVE_SPEED,d/fmaxf(dt,.001f));
    r->jump=fmaxf(0,r->jump-dt);
    if(stairs){
        if(height>.12f && p->onGround && r->jump<=0){p->vel.y=sqrtf(2*GRAVITY*.85f);p->onGround=false;r->jump=.3f;}
        Vector3 move=d>.01f?v_mul(to,speed/d):(Vector3){0};p->vel.x=move.x;p->vel.z=move.z;
    }else bot_apply_move(p,&botStates[i],to,speed,dt);
    float yaw,pitch;dir_to_yaw_pitch(to,&yaw,&pitch);p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-180*dt,180*dt);
    return false;
}
static bool tower_rifle(int i,float dt){
    TowerRifle *r=&towerRifles[i];if(!tower_eligible(i)){if(r->tag)tower_release(i);return false;}
    r->retry=fmaxf(0,r->retry-dt);EnemyStructure *t=tower_find(r->tag);
    if(!t){
        if(r->tag)tower_release(i);if(r->retry>0)return false;
        float nearest=24;
        for(int b=0;b<ENEMY_STRUCTURE_MAX;b++){EnemyStructure *candidate=&enemyStructures[b];if(!candidate->tower||!candidate->built||candidate->occupant>=0)continue;
            Vector3 at=tower_waypoint(candidate,0);float d=v_length(world_delta(at,players[i].pos));
            if(d<nearest && !is_view_occluded_by_voxels(players[i].pos,at)){nearest=d;t=candidate;}
        }
        if(!t)return false;t->occupant=i;*r=(TowerRifle){.tag=t->tag,.state=1,.target=-1,.last=players[i].pos};
    }
    Player *p=&players[i];int gun=fpCarrierType[i]>0?eq_gun(i):EQ_DEFAULT;
    if(p->matter<eqGuns[gun].cost && equipment[i].charge<=0){tower_release(i);return false;}
    ai_tick(i,dt);AiCombat *a=&aiCombat[i];
    if(r->state==3){
        r->jump=fmaxf(0,r->jump-dt);
        bool landed=p->onGround && p->vel.y<=0 && r->jump<=0;
        bool valid=ai_hostile(i,r->target);Vector3 goal=valid?ai_target_body(r->target)->pos:p->pos;
        if(!valid || !a->visible)r->lost+=dt;else r->lost=0;
        bool exposed=valid && ai_target_body(r->target)->matter<=0;
        Vector3 from=world_delta(goal,tower_point(t,(Vector3){1.5f,0,7.5f}));from.y=0;
        if(landed && (!exposed||r->lost>2||v_length(from)>10)){r->state=1;r->waypoint=0;return true;}
        Vector3 delta=world_delta(goal,p->pos);delta.y=0;
        if(landed)bot_apply_move(p,&botStates[i],delta,MOVE_SPEED*.85f,dt);
        ai_aim(i,goal,dt);if(valid && a->visible && v_length(world_delta(goal,p->pos))<MELEE_RANGE)perform_melee(i);
        return true;
    }
    if(r->state==1){
        Vector3 goal=tower_waypoint(t,r->waypoint);
        if(r->waypoint>0 && !tower_floor(goal,t->tag)){
            r->blocked+=dt;p->vel.x=p->vel.z=0;if(r->blocked>2)tower_release(i);return true;
        }r->blocked=0;
        if(v_length(world_delta(p->pos,r->last))<.01f)r->stuck+=dt;else r->stuck=0;r->last=p->pos;
        if(r->stuck>5){tower_release(i);return true;}
        if(tower_move(i,goal,dt,r->waypoint>0)){r->waypoint++;r->stuck=0;if(r->waypoint>9){r->state=2;p->vel.x=p->vel.z=0;}}
        return true;
    }
    if(p->pos.y<3.9f || !tower_floor(tower_waypoint(t,9),t->tag)){p->vel.x=p->vel.z=0;if(p->onGround){r->state=1;r->waypoint=0;}return true;}
    if(a->visible && ai_hostile(i,a->target) && !target_is_hologram(a->target) && ai_target_body(a->target)->matter<=0){
        Vector3 local=world_delta(ai_target_body(a->target)->pos,tower_point(t,(Vector3){1.5f,0,7.5f}));local.y=0;
        Vector3 relative=crosscut_rotate(world_delta(ai_target_body(a->target)->pos,t->origin),(4-t->quarter)%4);
        float dx=fmaxf(0,fmaxf(-relative.x,relative.x-3)),dz=fmaxf(0,fmaxf(-relative.z,relative.z-9));
        if(sqrtf(dx*dx+dz*dz)<=6){
            Vector3 direction=v_norm(local),goal=world_wrap(v_add(tower_point(t,(Vector3){1.5f,4.27f,7.5f}),v_mul(direction,2)));
            if(tower_move(i,goal,dt,false) || v_length(world_delta(p->pos,tower_waypoint(t,9)))>1){
                Vector3 probe=v_add(p->pos,v_add(v_mul(direction,.7f),(Vector3){0,.7f,0}));
                if(!player_aabb_hits_world_ex(probe,player_collision_half(p),i)){
                    p->vel=v_add(v_mul(direction,MOVE_SPEED*.8f),(Vector3){0,JUMP_SPEED*.55f,0});p->onGround=false;r->state=3;r->target=a->target;r->lost=0;r->jump=.4f;eq_cancel_charge(i);
                }
            }return true;
        }
    }
    p->vel.x=p->vel.z=0;
    if(!a->visible && a->target>=0){
        const Vector3 posts[4]={{.05f,4.27f,7.5f},{2.95f,4.27f,7.5f},{1.5f,4.27f,6.05f},{1.5f,4.27f,8.95f}};
        for(int n=0;n<4;n++){Vector3 post=tower_point(t,posts[n]);if(!is_view_occluded_by_voxels(post,a->known) && tower_floor(post,t->tag)){tower_move(i,post,dt,false);ai_aim(i,a->known,dt);eq_cancel_charge(i);return true;}}
    }
    if(!a->visible && a->target<0){
        r->lost+=dt;
        if(r->lost>3){
            const Vector3 posts[4]={{.05f,4.27f,7.5f},{1.5f,4.27f,6.05f},{2.95f,4.27f,7.5f},{1.5f,4.27f,8.95f}};
            Vector3 goal=tower_point(t,posts[r->post%4]);
            if(tower_floor(goal,t->tag) && tower_move(i,goal,dt,false)){r->post++;r->lost=0;}
            eq_cancel_charge(i);return true;
        }
    }else if(a->visible)r->lost=0;
    if(!a->visible){p->yaw+=45*dt;p->pitch=-15;eq_cancel_charge(i);equipment[i].held=false;return true;}
    ai_aim(i,a->known,dt);
    if(is_view_occluded_by_voxels(p->pos,a->known))return true;
    float distance=v_length(world_delta(a->known,p->pos));
    if(distance>eqGuns[gun].range || !ai_ready_to_fire(i))return true;
    if(gun==EQ_RAIL){NetInputCommand c={0};if(equipment[i].charge<eqRailCharge)c.held=NET_INPUT_FIRE;eq_command(i,&c,dt);}
    else if(gun==EQ_GLOB){NetInputCommand c={.held=NET_INPUT_FIRE};eq_command(i,&c,dt);}
    else if(gun==EQ_DEFAULT)FireVoxel(i);
    else if(eq_fire(i))ai_shot_fired(i);
    return true;
}
static void tower_draw(void){
    tower_blueprint();
    for(int b=0;b<ENEMY_STRUCTURE_MAX;b++){
        EnemyStructure *t=&enemyStructures[b];if(!t->tower || (!t->built && t->construction<=0))continue;
        if(!t->built)for(int k=0;k<towerCellCount;k++){Vector3 p=tower_cell(t,k);DrawCubeWires(p,.49f,.49f,.49f,Fade(SKYBLUE,.25f+.5f*t->construction/2));}
        if(t->pulse>0){Vector3 top=tower_point(t,(Vector3){1.5f,4.05f,7.5f});DrawCubeWires(top,3,.1f,3,Fade(SKYBLUE,t->pulse/.35f));if(netTransport.role==NET_ROLE_CLIENT || tower_owner(t))DrawCylinderEx(players[t->owner].pos,top,.025f,.025f,6,SKYBLUE);}
    }
}
