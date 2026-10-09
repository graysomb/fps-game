/* Specialist behavior uses the normal combat/voxel paths; the lab calls these too. */
static int enemy_voxel(uint64_t id) {
    if(!id) return -1;
    for(int v=0;v<voxel_count;v++) if(voxels[v].identity==id) return v;
    return -1;
}
static Color enemy_color(EnemyType t) {
    const Color c[7]={{170,120,70,255},{65,110,205,255},{245,140,30,255},{185,90,215,255},{225,200,60,255},{70,220,200,255},{160,70,40,255}};
    return t>=ENEMY_TYPE_SCRAPPER && t<ENEMY_TYPE_COUNT?c[t-ENEMY_TYPE_SCRAPPER]:GRAY;
}
static void enemy_cue(int i) {
    const int sounds[7]={SFX_TETHER,SFX_IMPACT,SFX_SHIELD,SFX_FIRE,SFX_GLUE_BREAK,SFX_SHIELD_BREAK,SFX_EXPLOSION};
    play_sfx(sounds[players[i].enemyType-ENEMY_TYPE_SCRAPPER]);
}
static float enemy_matter(EnemyType t) {
    const float m[FF_TYPE_COUNT]={100,25,100,50,75,50,60,50,60,200};
    return t==ENEMY_TYPE_FOUNDRY?200*(1+.5f*(firefightHumans-1)):m[t<ENEMY_TYPE_COUNT?t:0];
}
static void enemy_pin(int v,Vector3 pos,Vector3 vel) {
    Voxel *x=&voxels[v];
    x->pos=pos;x->vel=vel;x->sleeping=false;x->sleepFrames=x->lifeFrames=0;
    x->glueEligible=false;
    /* Moving an already pinned part does not change collision topology. */
    if(x->particles[0] && x->particles[0]->inv_mass!=0) goliath_set_armor_attached(x,true);
    x->gx=(int)floorf(pos.x/VOXEL_SIZE);x->gy=(int)floorf(pos.y/VOXEL_SIZE);x->gz=(int)floorf(pos.z/VOXEL_SIZE);
    for(int c=0;c<8;c++) if(x->particles[c]) {
        Particle *p=x->particles[c];
        p->pos=v_add(pos,(Vector3){corner_signs[c][0]*.25f,corner_signs[c][1]*.25f,corner_signs[c][2]*.25f});
        p->prev_pos=p->predicted_pos=p->pos;p->vel=vel;
    }
}
static Vector3 enemy_local(int i,Vector3 offset) {
    float a=-players[i].yaw*DEG2RAD;
    return v_add(players[i].pos,(Vector3){offset.x*cosf(a)+offset.z*sinf(a),offset.y,-offset.x*sinf(a)+offset.z*cosf(a)});
}
static bool enemy_attach(int i,int v,Vector3 offset,int group) {
    EnemySpecialState *s=&enemySpecial[i];
    if(v<0 || s->equipmentCount>=ENEMY_EQUIPMENT_MAX) return false;
    for(int face=0;face<6;face++) break_face_link(&voxels[v],face);
    s->equipment[s->equipmentCount++]=(EnemyEquipment){voxels[v].identity,offset,group};
    voxels[v].owner=voxels[v].activator=i;
    enemy_pin(v,enemy_local(i,offset),players[i].vel);
    return true;
}
static void enemy_plate(int i,Vector3 offset,int group) {
    Vector3 p=enemy_local(i,offset);
    int v=addVoxel(p.x,p.y,p.z,false,true,enemy_color(players[i].enemyType),0);
    enemy_attach(i,v,offset,group);
}
static int enemy_equipment_owner(int v) {
    if(v<0 || v>=voxel_count) return -1;
    uint64_t id=voxels[v].identity;
    /* A planted charge outlives its actor and must remain shootable after slot reuse. */
    for(int c=MAX_PLAYERS;c<MAX_COMBATANTS;c++)
        if(enemyCharges[c].active && enemyCharges[c].id==id) return c;
    int i=voxels[v].owner;
    if(i<MAX_PLAYERS || i>=MAX_COMBATANTS) return -1;
    for(int n=0;n<enemySpecial[i].equipmentCount;n++) if(enemySpecial[i].equipment[n].id==id) return i;
    return -1;
}
static void enemy_detach(int v) {
    if(v<0 || v>=voxel_count) return;
    uint64_t id=voxels[v].identity;
    goliath_set_armor_attached(&voxels[v],false);voxels[v].glueEligible=true;
    voxels[v].owner=voxels[v].activator=-1;
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        if(enemyCharges[i].active && enemyCharges[i].id==id) enemyCharges[i].fuse=-1;
        EnemySpecialState *s=&enemySpecial[i];
        for(int n=s->equipmentCount-1;n>=0;n--) if(s->equipment[n].id==id) {
            s->equipment[n]=s->equipment[--s->equipmentCount];
            if(players[i].enemyType==ENEMY_TYPE_ANCHOR) s->disabled=8;
        }
    }
}
static void enemy_release(int i) {
    EnemySpecialState *s=&enemySpecial[i];
    while(s->equipmentCount) {
        int v=enemy_voxel(s->equipment[s->equipmentCount-1].id);
        if(v>=0) enemy_detach(v);else s->equipmentCount--;
    }
    s->windup=0;s->state=0;s->disabled=8;
}
static void enemy_on_death(int i) {
    enemy_release(i);
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY && !enemyLab) {
        int pending=ffRemaining[ENEMY_TYPE_SWARMER];
        ffRemaining[ENEMY_TYPE_SWARMER]=0;
        firefightWaveEnemiesTotal-=pending;firefightWaveEnemiesRemaining-=pending;
    }
}
static void enemy_cleanup(bool removeEquipment) {
    enemyCleanupPending=false;
    if(removeEquipment) for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        for(int n=0;n<enemySpecial[i].equipmentCount;n++) {
            int v=enemy_voxel(enemySpecial[i].equipment[n].id);if(v>=0) remove_voxel_index(v);
        }
        memset(&enemySpecial[i],0,sizeof(enemySpecial[i]));
    }
    for(int i=0;i<MAX_COMBATANTS;i++) {
        int v=enemy_voxel(enemyCharges[i].id);if(v>=0) remove_voxel_index(v);
    }
    for(int i=0;i<MAX_COMBATANTS*2;i++) {int v=enemy_voxel(enemyLobs[i].id);if(v>=0) remove_voxel_index(v);}
    for(int i=0;i<ENEMY_STRUCTURE_MAX;i++) for(int n=0;n<enemyStructures[i].count;n++) {
        int v=enemy_voxel(enemyStructures[i].ids[n]);if(v>=0) remove_voxel_index(v);
    }
    memset(enemyCharges,0,sizeof(enemyCharges));memset(enemyLobs,0,sizeof(enemyLobs));memset(enemyStructures,0,sizeof(enemyStructures));
}
static void enemy_init(int i) {
    EnemySpecialState *s=&enemySpecial[i];memset(s,0,sizeof(*s));
    s->generation=++enemyGeneration;s->cooldown=2;s->spawnTimer=15;s->target=-1;
    players[i].matter=players[i].matterMax=enemy_matter(players[i].enemyType);
    if(players[i].enemyType==ENEMY_TYPE_SHIELDBEARER)
        for(int x=0;x<4;x++) for(int y=0;y<5;y++) enemy_plate(i,(Vector3){(x-1.5f)*.5f,(y-1.5f)*.5f,-1},0);
    if(players[i].enemyType==ENEMY_TYPE_ANCHOR) enemy_plate(i,(Vector3){0,.9f,.65f},0);
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) {
        for(int g=0;g<4;g++) for(int y=0;y<4;y++)
            enemy_plate(i,(Vector3){g%2?1.25f:-1.25f,-BASE_EYE_HEIGHT+.25f+y*.5f,g/2?1.25f:-1.25f},g);
        for(int x=0;x<4;x++) for(int z=0;z<2;z++) enemy_plate(i,(Vector3){(x-1.5f)*.5f,1.6f,(z-.5f)*.5f},4);
    }
}
static void enemy_update_equipment(int i) {
    EnemySpecialState *s=&enemySpecial[i];int before=s->equipmentCount,write=0,groups=0;
    for(int n=0;n<s->equipmentCount;n++) {
        EnemyEquipment e=s->equipment[n];int v=enemy_voxel(e.id);
        if(v<0 || !voxels[v].simulate || tetherTag[v]>0 || voxels[v].owner!=i) continue;
        enemy_pin(v,enemy_local(i,e.offset),players[i].vel);
        s->equipment[write++]=e;if(e.group<4) groups|=1<<e.group;
    }
    s->equipmentCount=write;
    if(players[i].enemyType==ENEMY_TYPE_ANCHOR && before && !write) s->disabled=8;
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) {
        int lost=(~groups)&15, count=0;for(int g=0;g<4;g++) count+=(lost>>g)&1;
        if(count>=2 && lost!=s->lostMask) s->exposure=10;
        s->lostMask=lost;
    }
}
static bool enemy_core_open(int i) {
    return players[i].enemyType!=ENEMY_TYPE_FOUNDRY || enemySpecial[i].lostMask==15 || enemySpecial[i].exposure>0;
}
static bool enemy_clear_line(Vector3 a,Vector3 b) {
    Vector3 d=v_sub(b,a);float len=v_length(d);
    return len<.1f || first_voxel_hit((Ray){a,v_mul(d,1/len)},len,-1)<0;
}
static bool enemy_braces(int anchor,int target) {
    if(anchor==target || players[anchor].respawn_timer>0 || players[target].respawn_timer>0 ||
       players[anchor].enemyType!=ENEMY_TYPE_ANCHOR || enemySpecial[anchor].disabled>0 ||
       !enemySpecial[anchor].equipmentCount || combatantTeam[anchor]!=combatantTeam[target]) return false;
    if(v_length(v_sub(players[anchor].pos,players[target].pos))>8 ||
       !enemy_clear_line(players[anchor].pos,players[target].pos)) return false;
    int rank=0;
    for(int j=MAX_PLAYERS;j<target;j++) if(j!=anchor && players[j].respawn_timer<=0 &&
        v_length(v_sub(players[anchor].pos,players[j].pos))<=8 && enemy_clear_line(players[anchor].pos,players[j].pos)) rank++;
    return rank<3;
}
static float enemy_knockback_scale(int i) {
    for(int a=MAX_PLAYERS;a<activePlayers;a++) if(enemy_braces(a,i)) return .25f;
    return 1;
}
static bool enemy_scrap_eligible(int v,int i) {
    Voxel *x=&voxels[v];
    if(!x->simulate || x->isBullet || x->goliathProjectile || x->type || x->owner>=0 ||
       x->activator>=0 || tetherTag[v]>0 || tetherHoldClusterTag[v]>0 || enemy_equipment_owner(v)>=0 ||
       get_goliath_owner_of_voxel(v)>=0 || v_length(x->vel)>.5f || x->rest_edge>VOXEL_SIZE*1.1f) return false;
    return v_length(v_sub(x->pos,players[i].pos))<=6;
}
static uint64_t enemy_find_scrap(int i) {
    for(int v=0;v<voxel_count;v++) if(enemy_scrap_eligible(v,i)) return voxels[v].identity;
    return 0;
}
static bool enemy_collect(int i) {
    EnemySpecialState *s=&enemySpecial[i];int v=enemy_voxel(s->selected);
    if(v<0 || !enemy_scrap_eligible(v,i)) return false;
    int n=s->equipmentCount;
    Vector3 offset={(n%4-1.5f)*.5f,(n/4-1)*.5f,.85f};
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) offset=(Vector3){(s->repairs%4-1.5f)*.5f,1.1f,(s->repairs/4-1)*.5f};
    return enemy_attach(i,v,offset,4);
}
static uint64_t enemy_sapper_target(int i,int target) {
    float best=-1;uint64_t found=0;
    Vector3 threat=v_norm(v_sub(players[target].pos,players[i].pos));
    for(int v=0;v<voxel_count;v++) {
        Voxel *x=&voxels[v];if(x->simulate || x->pos.y>players[i].pos.y+1.5f) continue;
        Vector3 d=v_sub(x->pos,players[i].pos);float dist=v_length(d);
        if(dist>10 || dist<.2f || v_dot(v_norm(d),threat)<.2f) continue;
        float score=12-dist+(table_get_static_only(x->gx,x->gy+1,x->gz)>=0?10:0);
        if(score>best) {best=score;found=x->identity;}
    }
    return found;
}
static bool enemy_plant(int i) {
    EnemySpecialState *s=&enemySpecial[i];int support=enemy_voxel(s->selected);
    if(support<0 || voxels[support].simulate || enemyCharges[i].active) return false;
    Vector3 d=v_sub(players[i].pos,voxels[support].pos);
    if(v_length(d)>2.5f) return false;
    Vector3 p=v_add(voxels[support].pos,v_mul(v_norm(d),.55f));
    int v=addVoxel(p.x,p.y,p.z,false,true,ORANGE,0);if(v<0) return false;
    enemyCharges[i]=(EnemyCharge){true,i,voxels[v].identity,s->selected,p,4,0};
    voxels[v].owner=i;enemy_pin(v,p,(Vector3){0});return true;
}
static void enemy_blast(EnemyCharge c) {
    int gx=(int)floorf(c.pos.x/.5f),gy=(int)floorf(c.pos.y/.5f),gz=(int)floorf(c.pos.z/.5f);
    /* Activate first so surviving overhead pieces enter ordinary debris physics. */
    activate_static_voxels_near_region(gx-6,gx+6,gy-6,gy+6,gz-6,gz+6,c.owner);
    for(int v=voxel_count-1;v>=0;v--) if(v_length(v_sub(voxels[v].pos,c.pos))<1.25f) remove_voxel_index(v);
    for(int j=0;j<firefightHumans;j++) if(players[j].respawn_timer<=0) {
        float d=v_length(v_sub(players[j].pos,c.pos));if(d>3) continue;
        if(players[j].isExposed) kill_player(j,c.owner,true,false);
        else apply_matter_damage(j,c.owner,30*(1-d/3));
    }
    play_sfx(SFX_EXPLOSION);meshDirty=true;
}
static bool enemy_lob_velocity(Vector3 from,Vector3 target,Vector3 *velocity) {
    float d=v_length(v_sub(target,from));float t=clampf(d/12,1.1f,2.4f);
    *velocity=v_mul(v_sub(target,from),1/t);velocity->y+=.5f*GRAVITY*t;
    Vector3 last=from;
    for(int n=1;n<=24;n++) {
        float at=t*n/24;Vector3 next=v_add(from,v_mul(*velocity,at));next.y-=.5f*GRAVITY*at*at;
        Vector3 segment=v_sub(next,last);float length=v_length(segment);VoxelHit hit;
        if(fabsf(next.x)>FLOOR_SIZE-.5f || fabsf(next.z)>FLOOR_SIZE-.5f ||
            first_voxel_hit_detailed((Ray){last,v_norm(segment)},length,-1,true,&hit)) return false;
        last=next;
    }
    return true;
}
static bool enemy_launch(int i) {
    EnemySpecialState *s=&enemySpecial[i];int count=0,free=-1;
    for(int n=0;n<MAX_COMBATANTS*2;n++) {if(enemyLobs[n].id && enemyLobs[n].owner==i) count++;if(!enemyLobs[n].id) free=n;}
    if(count>=2 || free<0) return false;
    Vector3 from=v_add(players[i].pos,(Vector3){0,2.3f,0}),vel;
    if(!enemy_lob_velocity(from,s->aim,&vel)) return false;
    int v=addVoxel(from.x,from.y,from.z,false,true,enemy_color(players[i].enemyType),0);if(v<0) return false;
    voxels[v].owner=voxels[v].activator=i;voxels[v].goliathProjectile=true;voxels[v].glueEligible=false;
    voxels[v].tetherThrowCcdFrames=TETHER_THROW_CCD_FRAMES;set_voxel_velocity(&voxels[v],vel);
    enemyLobs[free]=(EnemyLob){voxels[v].identity,i,0,s->aim};play_sfx(SFX_TETHER);return true;
}
static bool enemy_build_valid(int i,Vector3 origin,int kind) {
    if(fabsf(origin.x)>FLOOR_SIZE-5 || fabsf(origin.z)>FLOOR_SIZE-5) return false;
    /* A clear apron guarantees routes around both ends; never seal a doorway. */
    for(float x=-3;x<=3;x+=.5f) for(float z=-3;z<=3;z+=.5f)
        for(int y=0;y<4;y++) if(table_get_static_only((int)floorf((origin.x+x)/.5f),y,(int)floorf((origin.z+z)/.5f))>=0) return false;
    for(int j=0;j<activePlayers;j++) if(players[j].respawn_timer<=0 &&
        fabsf(players[j].pos.x-origin.x)<2.5f && fabsf(players[j].pos.z-origin.z)<2.5f) return false;
    for(int p=0;p<MAX_PICKUPS;p++) if(pickups[p].active && v_length(v_sub(pickups[p].pos,origin))<4) return false;
    for(int p=0;p<8;p++) if(v_length(v_sub(crosscut_spawn_point(p),origin))<4) return false;
    if(fabsf(origin.x)>24 && fabsf(origin.z)<9) return false; /* Foundry entrance lanes. */
    int count=0,live=0;
    for(int b=0;b<ENEMY_STRUCTURE_MAX;b++) if(enemyStructures[b].owner==i && enemyStructures[b].generation==enemySpecial[i].generation) {
        bool any=false;for(int n=0;n<enemyStructures[b].count;n++) if(enemy_voxel(enemyStructures[b].ids[n])>=0) {live++;any=true;}
        count+=any;
    }
    return count<2 && live+(kind?40:24)<=96;
}
static bool enemy_build(int i) {
    EnemySpecialState *s=&enemySpecial[i];
    if(!enemy_build_valid(i,s->buildOrigin,s->buildKind)) return false;
    int slot=-1;
    for(int b=0;b<ENEMY_STRUCTURE_MAX;b++) {
        bool alive=false;for(int n=0;n<enemyStructures[b].count;n++) if(enemy_voxel(enemyStructures[b].ids[n])>=0) alive=true;
        if(!alive) {slot=b;break;}
    }
    if(slot<0) return false;
    EnemyStructure *b=&enemyStructures[slot];memset(b,0,sizeof(*b));b->owner=i;b->generation=s->generation;
    int ox=(int)floorf(s->buildOrigin.x/.5f),oz=(int)floorf(s->buildOrigin.z/.5f);
    int depth=s->buildKind?4:1,width=s->buildKind?4:6;
    for(int z=0;z<depth;z++) for(int x=0;x<width;x++) for(int y=0;y<(s->buildKind?z+1:4);y++) {
        Vector3 cell=crosscut_rotate((Vector3){x-width/2,y,z},s->buildQuarter);
        int v=add_static_voxel_at_grid(ox+(int)cell.x,y,oz+(int)cell.z,enemy_color(ENEMY_TYPE_BUILDER),0);
        if(v>=0) {b->ids[b->count++]=voxels[v].identity;voxels[v].owner=ENEMY_CONSTRUCTION_OWNER;}
    }
    rebuild_all_voxel_surfaces();rebuild_glue_constraints();meshDirty=true;return b->count>0;
}
static int enemy_target(int i) {
    int best=-1;float distance=1000;
    for(int j=0;j<firefightHumans;j++) if(players[j].respawn_timer<=0) {
        float d=v_length(v_sub(players[j].pos,players[i].pos));
        if(d<distance) {distance=d;best=j;}
    }
    return best;
}
static void enemy_tick(float dt) {
    if(enemyCleanupPending) enemy_cleanup(true);
    if(enemyLab && enemyLabPaused && !enemyLabStep) return;
    if(enemyLab && players[0].respawn_timer>0) {
        players[0].respawn_timer-=dt;
        if(players[0].respawn_timer<=0) {ff_restore_human(0,100);players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};}
    }
    for(int i=MAX_PLAYERS;i<activePlayers;i++) if(enemy_is_special(i) && players[i].respawn_timer<=0) {
        enemySpecial[i].exposure=fmaxf(0,enemySpecial[i].exposure-dt);
        enemySpecial[i].disabled=fmaxf(0,enemySpecial[i].disabled-dt);
        enemy_update_equipment(i);
    }
    for(int i=0;i<MAX_COMBATANTS;i++) {
        EnemyCharge *c=&enemyCharges[i];if(!c->active) continue;
        int v=enemy_voxel(c->id),support=enemy_voxel(c->support);
        bool disarmed=v<0 || support<0 || (support>=0 && voxels[support].simulate) || c->fuse<0;
        c->fuse-=dt;c->beep-=dt;
        if(disarmed || c->fuse<=0) {
            EnemyCharge copy=*c;c->active=false;c->id=0;
            if(v>=0) remove_voxel_index(v);
            if(!disarmed) enemy_blast(copy);
        } else {
            enemy_pin(v,c->pos,(Vector3){0});
            if(c->beep<=0) {play_sfx(SFX_SHIELD);c->beep=fmaxf(.12f,c->fuse*.2f);}
        }
    }
    for(int n=0;n<MAX_COMBATANTS*2;n++) if(enemyLobs[n].id) {
        EnemyLob *l=&enemyLobs[n];int v=enemy_voxel(l->id);l->age+=dt;
        if(v<0) l->id=0;
        else if(l->age>10 && tetherTag[v]<=0) {remove_voxel_index(v);l->id=0;}
        else if(tetherTag[v]>0) { /* Catching transfers projectile attribution through the normal tether path. */
            voxels[v].owner=voxels[v].activator=tetherTag[v]-1;
        }
    }
}
static void enemy_update(int i,float dt) {
    if(enemyLab && enemyLabPaused && !enemyLabStep) {players[i].vel=(Vector3){0};return;}
    EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];EnemyType t=p->enemyType;
    s->cooldown-=dt;s->melee-=dt;s->repair-=dt;s->target=enemy_target(i);
    if(s->target<0) {p->vel.x=p->vel.z=0;return;}
    Vector3 to=v_sub(players[s->target].pos,p->pos);float distance=v_length(to);
    bool visible=enemy_clear_line(p->pos,players[s->target].pos);
    float yaw=atan2f(-to.x,-to.z)*RAD2DEG;
    float turn=t==ENEMY_TYPE_SHIELDBEARER?60:180;
    p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-turn*dt,turn*dt);
    if(p->isExposed) {s->windup=0;s->state=0;}
    if(s->windup>0) {
        s->windup-=dt;p->vel.x=p->vel.z=0;
        if(s->windup<=0) {
            if(s->state==1) enemy_collect(i);
            if(s->state==2) enemy_plant(i);
            if(s->state==3) enemy_launch(i);
            if(s->state==4) enemy_build(i);
            s->state=0;
        }
        return;
    }
    Vector3 move=to;float speed=1;
    if(t==ENEMY_TYPE_SHIELDBEARER) {
        speed=.7f;
        /* Interpose toward the threat, ahead of the nearest unshielded ally. */
        for(int j=MAX_PLAYERS;j<activePlayers;j++) if(j!=i && players[j].respawn_timer<=0 && v_length(v_sub(p->pos,players[j].pos))<8) {
            Vector3 destination=v_add(players[j].pos,v_mul(v_norm(v_sub(players[s->target].pos,players[j].pos)),2));
            move=v_sub(destination,p->pos);break;
        }
    }
    if(t==ENEMY_TYPE_LOBBER || t==ENEMY_TYPE_FOUNDRY) {
        speed=t==ENEMY_TYPE_FOUNDRY?.35f:.75f;
        if(distance<10) move=v_mul(to,-1);else if(distance<24) move=(Vector3){0};
        if(!p->isExposed && s->cooldown<=0 && visible && distance<=24) {
            s->aim=players[s->target].pos;Vector3 vel;
            if(enemy_lob_velocity(v_add(p->pos,(Vector3){0,2.3f,0}),s->aim,&vel)) {
                s->windup=1.5f;s->state=3;enemy_cue(i);
            }
            s->cooldown=t==ENEMY_TYPE_FOUNDRY?6:5;
        }
    }
    if(t==ENEMY_TYPE_SCRAPPER && !p->isExposed && s->cooldown<=0 && s->equipmentCount<12) {
        s->selected=enemy_find_scrap(i);s->cooldown=3;
        if(s->selected) {s->state=1;s->windup=1;enemy_cue(i);}
    }
    if(t==ENEMY_TYPE_SAPPER) {
        if(enemyCharges[i].active) move=v_sub(p->pos,enemyCharges[i].pos);
        else if(!p->isExposed && s->cooldown<=0) {
            int v=enemy_voxel(s->selected);
            if(v<0 || voxels[v].simulate) {s->selected=enemy_sapper_target(i,s->target);v=enemy_voxel(s->selected);}
            if(v>=0) {
                move=v_sub(voxels[v].pos,p->pos);
                if(v_length(move)<=2.5f) {s->state=2;s->windup=1.5f;s->cooldown=8;enemy_cue(i);}
            } else s->cooldown=1;
        }
    }
    if(t==ENEMY_TYPE_BUILDER) {
        speed=.85f;
        if(!p->isExposed && s->cooldown<=0) {
            s->buildKind=0;
            for(int j=MAX_PLAYERS;j<activePlayers;j++) if(j!=i && players[j].respawn_timer<=0 &&
                players[s->target].pos.y-players[j].pos.y>.5f && players[s->target].pos.y-players[j].pos.y<=2 &&
                v_length(v_sub(p->pos,players[j].pos))<8) s->buildKind=1;
            Vector3 direction=to;direction.y=0;direction=v_norm(direction);
            s->buildQuarter=fabsf(direction.x)>fabsf(direction.z)?(direction.x>0?3:1):(direction.z>0?0:2);
            s->buildOrigin=v_add(p->pos,v_mul(direction,4));s->buildOrigin.y=0;
            if(enemy_build_valid(i,s->buildOrigin,s->buildKind)) {s->windup=2;s->state=4;enemy_cue(i);}
            s->cooldown=8;
        }
    }
    if(t==ENEMY_TYPE_ANCHOR) {
        speed=.65f;
        if(!s->equipmentCount && s->disabled<=0 && !p->isExposed) enemy_plate(i,(Vector3){0,.9f,.65f},0);
        if(distance<7) move=v_mul(to,-1);else if(distance<12) move=(Vector3){0};
    }
    if(t==ENEMY_TYPE_FOUNDRY) {
        if(s->repair<=0 && s->repairs<12 && s->equipmentCount<ENEMY_EQUIPMENT_MAX && !p->isExposed) {
            s->selected=enemy_find_scrap(i);if(enemy_collect(i)) s->repairs++;s->repair=3;
        }
        s->spawnTimer-=dt;
        if(s->spawnTimer<=0 && s->spawned<6) {
            int added=0;
            for(int j=MAX_PLAYERS;j<activePlayers && added<2;j++) if(players[j].respawn_timer>0) {
                if(!enemyLab && ffRemaining[ENEMY_TYPE_SWARMER]<=0) break;
                Vector3 at=v_add(p->pos,(Vector3){added?3:-3,0,3});
                at=nudge_spawn_clear(at,PLAYER_SIZE*.25f,1,j);
                if(!spawn_position_clear_ex(at,PLAYER_SIZE*.25f,1,j)) continue;
                init_firefight_bot_entity(j,ENEMY_TYPE_SWARMER,INPUT_TYPE_BOT_EASY);players[j].pos=at;
                if(!enemyLab) {ffRemaining[ENEMY_TYPE_SWARMER]--;firefightEnemiesSpawned++;}
                s->spawned++;added++;
            }
            s->spawnTimer=added?15:1;
        }
    }
    if(distance<1.4f && visible && s->melee<=0) {
        if(players[s->target].isExposed) kill_player(s->target,i,true,false);
        else apply_matter_damage(s->target,i,10);
        s->melee=1.5f;
    }
    move.y=0;
    if(s->windup>0 || v_length(move)<.25f) {p->vel.x=p->vel.z=0;return;}
    bot_apply_move(p,&botStates[i],move,MOVE_SPEED*speed,dt);
}
static void enemy_draw(void) {
    for(int i=MAX_PLAYERS;i<activePlayers;i++) if(enemy_is_special(i) && players[i].respawn_timer<=0) {
        EnemySpecialState *s=&enemySpecial[i];Vector3 p=players[i].pos;Color c=enemy_color(players[i].enemyType);
        /* Distinct silhouettes in addition to color: pack, mast, mortar, hardhat. */
        if(players[i].enemyType==ENEMY_TYPE_SCRAPPER) DrawCube(v_add(p,(Vector3){0,.1f,.5f}),.8f,.8f,.4f,c);
        if(players[i].enemyType==ENEMY_TYPE_SAPPER) {DrawCube(v_add(p,(Vector3){0,.4f,0}),.7f,.3f,.7f,c);DrawSphere(v_add(p,(Vector3){0,.7f,0}),.15f,ORANGE);}
        if(players[i].enemyType==ENEMY_TYPE_BUILDER) DrawCube(v_add(p,(Vector3){0,.4f,0}),1,.2f,1,c);
        if(players[i].enemyType==ENEMY_TYPE_LOBBER) DrawCylinder(v_add(p,(Vector3){0,.25f,0}),.4f,.3f,1,8,c);
        if(players[i].enemyType==ENEMY_TYPE_ANCHOR) {
            DrawLine3D(p,v_add(p,(Vector3){0,1.1f,0}),c);
            for(int j=MAX_PLAYERS;j<activePlayers;j++) if(enemy_braces(i,j)) DrawLine3D(v_add(p,(Vector3){0,.9f,0}),players[j].pos,c);
        }
        if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) DrawCubeWires(p,2.5f,2.5f,2.5f,enemy_core_open(i)?GOLD:c);
        if(s->windup>0) {
            DrawCircle3D(p,.9f,(Vector3){1,0,0},90,c);
            if(s->state==3) {DrawCircle3D(s->aim,1.5f,(Vector3){1,0,0},90,ORANGE);DrawLine3D(p,s->aim,Fade(c,.5f));}
            if(s->state==4) DrawCubeWires(v_add(s->buildOrigin,(Vector3){0,1,s->buildKind?1:0}),3,2,s->buildKind?2:.5f,YELLOW);
        }
    }
    for(int i=0;i<MAX_COMBATANTS;i++) if(enemyCharges[i].active) {
        EnemyCharge *c=&enemyCharges[i];Color col=((int)(c->fuse*8)%2)?RED:YELLOW;
        DrawSphere(c->pos,.35f,col);DrawCircle3D(c->pos,3,(Vector3){1,0,0},90,Fade(col,.7f));
    }
    for(int i=0;i<MAX_COMBATANTS*2;i++) if(enemyLobs[i].id) DrawCircle3D(enemyLobs[i].target,1.5f,(Vector3){1,0,0},90,ORANGE);
}
