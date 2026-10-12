/* Specialist behavior uses the normal combat/voxel paths; the lab calls these too. */
static int enemy_voxel(uint64_t id) {
    if(!id)return -1;
    /* Validate cached indices after compaction; identities are never recycled. */
    static uint64_t keys[4096];static int indices[4096];unsigned slot=(unsigned)(id^(id>>32))&4095;
    int cached=indices[slot];
    if(keys[slot]==id && cached>=0 && cached<voxel_count && voxels[cached].identity==id)return cached;
    for(int v=0;v<voxel_count;v++)if(voxels[v].identity==id){keys[slot]=id;indices[slot]=v;return v;}
    return -1;
}
static Color enemy_color(EnemyType t) {
    const Color c[7]={{170,120,70,255},{65,110,205,255},{245,140,30,255},{185,90,215,255},{225,200,60,255},{70,220,200,255},{160,70,40,255}};
    return t>=ENEMY_TYPE_SCRAPPER && t<ENEMY_TYPE_COUNT?c[t-ENEMY_TYPE_SCRAPPER]:GRAY;
}
static void enemy_cue(int i) {
    const int sounds[7]={SFX_TETHER,SFX_IMPACT,SFX_SHIELD,SFX_FIRE,SFX_GLUE_BREAK,SFX_SHIELD_BREAK,SFX_EXPLOSION};
    play_sfx_at(sounds[players[i].enemyType-ENEMY_TYPE_SCRAPPER],players[i].pos,i,2,1);
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
    float a=players[i].yaw*DEG2RAD;
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
    enemy_drop_payload(i);s->windup=0;s->state=0;s->disabled=8;
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
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++)enemy_drop_payload(i);
    if(removeEquipment) for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        for(int n=0;n<enemySpecial[i].equipmentCount;n++) {
            int v=enemy_voxel(enemySpecial[i].equipment[n].id);if(v>=0) remove_voxel_index(v);
        }
        memset(&enemySpecial[i],0,sizeof(enemySpecial[i]));
    }
    for(int i=0;i<MAX_COMBATANTS;i++) {
        int v=enemy_voxel(enemyCharges[i].id);if(v>=0) remove_voxel_index(v);
    }
    for(int i=0;i<MAX_COMBATANTS*2;i++)for(int k=0;k<enemyLobs[i].count;k++) {
        int v=enemy_voxel(enemyLobs[i].ids[k]);if(v>=0 && voxels[v].owner==enemyLobs[i].owner) {
            voxels[v].goliathProjectile=false;voxels[v].glueEligible=true;voxels[v].owner=voxels[v].activator=-1;
        }
    }
    for(int i=0;i<ENEMY_STRUCTURE_MAX;i++) for(int n=0;n<enemyStructures[i].count;n++) {
        int v=enemy_voxel(enemyStructures[i].ids[n]);if(v>=0) remove_voxel_index(v);
    }
    memset(enemyCharges,0,sizeof(enemyCharges));memset(enemyLobs,0,sizeof(enemyLobs));memset(enemyStructures,0,sizeof(enemyStructures));
}
static void enemy_init(int i) {
    EnemySpecialState *s=&enemySpecial[i];memset(s,0,sizeof(*s));
    s->generation=++enemyGeneration;s->cooldown=2;s->spawnTimer=15;s->target=-1;s->reason="seeking";
    s->coreMax=s->coreHealth=enemy_matter(players[i].enemyType);
    players[i].matter=players[i].matterMax=enemy_matter(players[i].enemyType);
    if(players[i].enemyType==ENEMY_TYPE_SHIELDBEARER)
        for(int x=0;x<4;x++) for(int y=0;y<5;y++) enemy_plate(i,(Vector3){(x-1.5f)*.5f,(y-1.5f)*.5f,-1},x);
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
        if(count>=2)s->coreOpen=true;
        s->lostMask=lost;
    }
}
static bool enemy_core_open(int i) {
    return players[i].enemyType!=ENEMY_TYPE_FOUNDRY || enemySpecial[i].coreOpen || enemySpecial[i].lostMask==15;
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
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY)return .25f;
    for(int a=MAX_PLAYERS;a<activePlayers;a++) if(enemy_braces(a,i)) return .25f;
    return 1;
}
static bool enemy_scrap_eligible(int v,int i) {
    Voxel *x=&voxels[v];
    if(!x->simulate || x->isBullet || x->goliathProjectile || x->type || tetherTag[v]>0 || tetherHoldClusterTag[v]>0 || enemy_equipment_owner(v)>=0 ||
       get_goliath_owner_of_voxel(v)>=0 || v_length(x->vel)>1.0f || x->rest_edge>VOXEL_SIZE*1.1f) return false;
    return v_length(v_sub(x->pos,players[i].pos))<=16;
}
static uint64_t enemy_find_scrap(int i) {
    uint64_t result=0;float best=17;
    for(int v=0;v<voxel_count;v++)if(enemy_scrap_eligible(v,i)) {
        float d=v_length(v_sub(voxels[v].pos,players[i].pos));
        if(d<best){best=d;result=voxels[v].identity;}
    }
    return result;
}
static bool enemy_collect(int i) {
    EnemySpecialState *s=&enemySpecial[i];int v=enemy_voxel(s->selected);
    if(v<0 || !enemy_scrap_eligible(v,i) || v_length(v_sub(voxels[v].pos,players[i].pos))>2.5f) return false;
    int n=s->equipmentCount;
    Vector3 offset={(n%3-1)*.5f,(n/6)*.5f-.25f,n%6<3?-.7f:.7f};
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
static void enemy_blast_terrain(EnemyCharge c) {
    /* Excavation activation skips its interior; a blast must activate the radial shell explicitly. */
    static UnitVoxelBuffer shell;unit_voxel_buffer_clear(&shell);
    for(int v=0;v<voxel_count && shell.count<VOXEL_ACTIVATION_UNIT_BUDGET;v++) {
        Voxel *x=&voxels[v];float distance=v_length(v_sub(x->pos,c.pos));
        if(!x->simulate && !x->type && distance>=2 && distance<4)
            unit_voxel_buffer_push(&shell,x->gx,x->gy,x->gz,x->color,x->type,false,v,x->debugClusterTag,c.owner);
    }
    remove_buffered_static_voxels(&shell);emit_unit_voxels_from_units(&shell,false,true,-1);
    init_static_hash();rebuild_voxel_hash();
    for(int v=voxel_count-1;v>=0;v--)if(v_length(v_sub(voxels[v].pos,c.pos))<2 && enemy_equipment_owner(v)<0)remove_voxel_index(v);
    int stamp=++particle_sync_stamp;
    for(int v=0;v<voxel_count;v++) {
        Voxel *x=&voxels[v];float distance=v_length(v_sub(x->pos,c.pos));
        if(!x->simulate || distance>=4 || enemy_equipment_owner(v)>=0)continue;
        if(x->sleeping)wake_sleeping_cluster_in_place(v);
        for(int f=0;f<6;f++)break_face_link(x,f);
        Vector3 direction=v_sub(x->pos,c.pos);direction.y+=1;
        Vector3 impulse=v_mul(v_norm(direction),12*(1-distance/4));
        x->owner=x->activator=c.owner;x->explosionUntil=clarityClock+10;x->sleeping=false;x->vel=v_add(x->vel,impulse);
        for(int n=0;n<8;n++)if(x->particles[n] && x->particles[n]->sync_stamp!=stamp) {
            Particle *part=x->particles[n];part->sync_stamp=stamp;
            if(part->inv_mass>0)part->vel=v_add(part->vel,impulse);
        }
    }
    rebuild_voxel_hash();rebuild_all_voxel_surfaces();
    meshDirty=true;
}
/* Supplies are real glued voxel bundles, shared by players and enemies. */
static void enemy_scrap_bundle(Vector3 origin) {
    int first=voxel_count;
    for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int z=0;z<2;z++) {
        int v=addVoxel(origin.x+x*.5f+.25f,origin.y+y*.5f+.25f,origin.z+z*.5f+.25f,false,true,BROWN,0);
        if(v>=0){voxels[v].owner=ENEMY_SCRAP_OWNER;voxels[v].activator=-1;}
    }
    int members[8],count=0;
    for(int v=first;v<voxel_count;v++){glue_neighbor_faces_for_voxel(v);members[count++]=v;}
    if(count)freeze_dynamic_cluster_in_place(members,count);
}
static void enemy_stock_scrap(bool initial) {
    for(int n=0;n<8;n++) {
        const float lanes[4]={-28,-4,4,28};
        Vector3 at={n<4?-29:29,0,lanes[n%4]};bool clear=true;
        for(int v=0;v<voxel_count && clear;v++)if(fabsf(voxels[v].pos.x-at.x-.5f)<1.25f && fabsf(voxels[v].pos.z-at.z-.5f)<1.25f && voxels[v].pos.y<1.5f)clear=false;
        for(int i=0;!initial && i<activePlayers && clear;i++)if(players[i].respawn_timer<=0 && v_length(v_sub(players[i].pos,at))<2)clear=false;
        if(clear)enemy_scrap_bundle(at);
    }
}
static bool enemy_payload(int i,int v,int *count,Vector3 *center,Vector3 *extent) {
    if(v<0 || v>=voxel_count || !voxels[v].simulate || voxels[v].isBullet || voxels[v].type || enemy_equipment_owner(v)>=0)return false;
    int n=build_glue_cluster_indices(v,glueClusterIndices);
    if(n<8 || n>ENEMY_LOB_MAX_VOXELS)return false;
    Vector3 lo={1e6f,1e6f,1e6f},hi={-1e6f,-1e6f,-1e6f};*center=(Vector3){0};
    for(int k=0;k<n;k++) {
        int id=glueClusterIndices[k];Voxel *x=&voxels[id];
        if(!x->simulate || x->type || x->isBullet || enemy_equipment_owner(id)>=0 ||
            (tetherTag[id]>0 && tetherTag[id]!=i+1) || (tetherHoldClusterTag[id]>0 && tetherHoldClusterTag[id]!=i+1))return false;
        if(!players[i].tetherHolding && (x->goliathProjectile || v_length(x->vel)>1))return false;
        *center=v_add(*center,x->pos);
        lo.x=fminf(lo.x,x->pos.x-.25f);lo.y=fminf(lo.y,x->pos.y-.25f);lo.z=fminf(lo.z,x->pos.z-.25f);
        hi.x=fmaxf(hi.x,x->pos.x+.25f);hi.y=fmaxf(hi.y,x->pos.y+.25f);hi.z=fmaxf(hi.z,x->pos.z+.25f);
    }
    *center=v_mul(*center,1.0f/n);*extent=v_mul(v_sub(hi,lo),.5f);*count=n;
    return extent->x<=1.5f && extent->y<=1.5f && extent->z<=1.5f;
}
static uint64_t enemy_find_payload(int i) {
    float best=20;uint64_t found=0;
    for(int v=0;v<voxel_count;v++) {
        if(!voxels[v].simulate || voxels[v].type || voxels[v].isBullet)continue;
        float d=v_length(v_sub(voxels[v].pos,players[i].pos));if(d>=best)continue;
        int count;Vector3 center,extent;
        if(enemy_payload(i,v,&count,&center,&extent)){found=voxels[v].identity;best=d;}
    }
    return found;
}
static void enemy_drop_payload(int i) {
    Player *p=&players[i];
    if(p->tetherHolding){Vector3 zero={0};release_tether_velocity(i,&zero);}
    p->goldTetherCharged=false;
    for(int v=0;v<voxel_count;v++)if(tetherTag[v]==i+1 || tetherHoldClusterTag[v]==i+1){tetherTag[v]=tetherHoldClusterTag[v]=0;}
}
static bool enemy_lob_arc(Vector3 from,Vector3 target,Vector3 extent,Vector3 *velocity) {
    /* Match the damped fixed-step PBD integrator, not vacuum ballistics. */
    const float dt=PBD_MAX_STEP_DT/PBD_SUBSTEPS,q=VELOCITY_DAMPING;
    int steps=0;Vector3 d=v_sub(target,from);
    for(int n=42*PBD_SUBSTEPS;n<=240*PBD_SUBSTEPS;n+=3*PBD_SUBSTEPS) {
        float response=0,weight=1,gravityVelocity=0,gravityOffset=0;
        for(int k=0;k<n;k++){weight*=q;response+=weight*dt;gravityVelocity=q*gravityVelocity-GRAVITY*dt;gravityOffset+=gravityVelocity*dt;}
        Vector3 candidate={d.x/response,(d.y-gravityOffset)/response,d.z/response};
        if(v_length(candidate)>40)continue;
        float y=from.y,vy=candidate.y,apex=y;for(int k=0;k<n;k++){vy=q*vy-GRAVITY*dt;y+=vy*dt;apex=fmaxf(apex,y);}
        if(apex>=fmaxf(from.y,target.y)+4){*velocity=candidate;steps=n;break;}
    }
    if(!steps)return false;
    Vector3 position=from,vel=*velocity,last=from;
    for(int n=1;n<=steps;n++) {
        vel=v_mul(vel,q);vel.y-=GRAVITY*dt;position=v_add(position,v_mul(vel,dt));
        if(n%4 && n!=steps)continue;
        for(int corner=0;corner<9;corner++) {
            Vector3 offset=corner==8?(Vector3){0}:(Vector3){corner_signs[corner][0]*extent.x,corner_signs[corner][1]*extent.y,corner_signs[corner][2]*extent.z};
            Vector3 segment=v_sub(position,last);VoxelHit hit;float length=v_length(segment);
            if(fabsf(position.x)>FLOOR_SIZE-.5f || fabsf(position.z)>FLOOR_SIZE-.5f || first_voxel_hit_detailed((Ray){v_add(last,offset),v_norm(segment)},length,-1,true,&hit))return false;
        }
        last=position;
    }
    return true;
}
static bool enemy_lob_velocity(Vector3 from,Vector3 target,Vector3 *velocity) {
    return enemy_lob_arc(from,target,(Vector3){.5f,.5f,.5f},velocity);
}
static bool enemy_launch(int i) {
    EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];int live=0,free=-1;
    for(int n=0;n<MAX_COMBATANTS*2;n++){if(enemyLobs[n].id && enemyLobs[n].owner==i)live++;if(!enemyLobs[n].id)free=n;}
    if(live>=2 || free<0 || !p->tetherHolding || !resolve_tether_voxel(p))return false;
    int count;Vector3 from,extent,velocity;
    if(!enemy_payload(i,p->tetherVoxel,&count,&from,&extent) || !enemy_lob_arc(from,s->aim,extent,&velocity))return false;
    EnemyLob lob={.id=voxels[p->tetherVoxel].identity,.owner=i,.count=count,.target=s->aim};
    for(int k=0;k<count;k++)lob.ids[k]=voxels[glueClusterIndices[k]].identity;
    release_tether_velocity(i,&velocity);
    for(int k=0;k<count;k++) {
        int v=enemy_voxel(lob.ids[k]);if(v<0)continue;
        voxels[v].goliathProjectile=true;voxels[v].glueEligible=false;
        voxels[v].tetherThrowCcdFrames=TETHER_THROW_CCD_FRAMES;track_tether_throw_ccd(&voxels[v]);
    }
    enemyLobs[free]=lob;return true;
}
/* A whole payload causes one bounded hit, not damage multiplied by its voxel count. */
static bool enemy_lob_hit(int voxel,int victim) {
    uint64_t id=voxels[voxel].identity;
    for(int n=0;n<MAX_COMBATANTS*2;n++)if(enemyLobs[n].id) {
        EnemyLob *l=&enemyLobs[n];bool member=false;
        for(int k=0;k<l->count;k++)if(l->ids[k]==id)member=true;
        if(!member)continue;
        if(l->hitCooldown[victim]<=0 && combatantTeam[victim]!=combatantTeam[l->owner]) {
            l->hitCooldown[victim]=.35f;
            if(players[victim].isExposed)kill_player(victim,l->owner,true,true);else apply_matter_damage(victim,l->owner,30);
        }
        return true;
    }
    return false;
}
static bool enemy_build_valid(int i,Vector3 origin,int kind) {
    if(fabsf(origin.x)>FLOOR_SIZE-5 || fabsf(origin.z)>FLOOR_SIZE-5) return false;
    /* A clear apron guarantees routes around both ends; never seal a doorway. */
    for(int x=-4;x<=4;x++)for(int z=-1;z<=(kind?4:1);z++)for(int y=0;y<4;y++) {
        Vector3 cell=crosscut_rotate((Vector3){x,y,z},enemySpecial[i].buildQuarter);
        if(table_get((int)floorf(origin.x/.5f)+(int)cell.x,y,(int)floorf(origin.z/.5f)+(int)cell.z)>=0)return false;
    }
    for(int j=0;j<activePlayers;j++) if(players[j].respawn_timer<=0 &&
        fabsf(players[j].pos.x-origin.x)<2.5f && fabsf(players[j].pos.z-origin.z)<2.5f) return false;
    for(int p=0;p<MAX_PICKUPS;p++) if(pickups[p].active && v_length(v_sub(pickups[p].pos,origin))<4) return false;
    if(fabsf(origin.x)>FLOOR_SIZE-8 && fabsf(origin.z)<9) return false; /* Foundry entrance lanes. */
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
static bool enemy_find_build_site(int i,Vector3 threat) {
    EnemySpecialState *s=&enemySpecial[i];float base=atan2f(threat.z,threat.x);
    for(int n=0;n<16;n++) {
        float a=base+(n%8)*PI/4,range=n<8?4:n<12?6:8;
        s->buildOrigin=v_add(players[i].pos,(Vector3){cosf(a)*range,0,sinf(a)*range});s->buildOrigin.y=0;
        s->buildQuarter=fabsf(cosf(a))>fabsf(sinf(a))?(cosf(a)>0?3:1):(sinf(a)>0?0:2);
        if(enemy_build_valid(i,s->buildOrigin,s->buildKind))return true;
    }
    return false;
}
static void enemy_stomp(int i) {
    for(int j=0;j<firefightHumans;j++)if(players[j].respawn_timer<=0 && v_length(v_sub(players[j].pos,players[i].pos))<=3) {
        apply_matter_damage(j,i,20);
        Vector3 away=v_sub(players[j].pos,players[i].pos);away.y=.6f;
        players[j].vel=v_add(players[j].vel,v_mul(v_norm(away),6));players[j].onGround=false;players[j].meleeKnockbackActive=true;
    }
    play_sfx_at(SFX_EXPLOSION,players[i].pos,i,3,.7f);
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
            if(c->beep<=0) {play_sfx_at(SFX_SHIELD,c->pos,c->owner,3,1);c->beep=fmaxf(.12f,c->fuse*.2f);}
        }
    }
    for(int n=0;n<MAX_COMBATANTS*2;n++)if(enemyLobs[n].id) {
        EnemyLob *l=&enemyLobs[n];l->age+=dt;int present=0,holder=-1;
        for(int j=0;j<MAX_COMBATANTS;j++)l->hitCooldown[j]=fmaxf(0,l->hitCooldown[j]-dt);
        for(int k=0;k<l->count;k++) {int v=enemy_voxel(l->ids[k]);if(v<0)continue;present++;if(tetherTag[v]>0)holder=tetherTag[v]-1;}
        if(holder>=0 || l->age>10 || !present) {
            for(int k=0;k<l->count;k++) {int v=enemy_voxel(l->ids[k]);if(v<0)continue;
                voxels[v].owner=voxels[v].activator=holder;voxels[v].goliathProjectile=false;voxels[v].glueEligible=true;
            }
            memset(l,0,sizeof(*l));
        }
    }
}
static void enemy_update(int i,float dt) {
    if(enemyLab && enemyLabPaused && !enemyLabStep){players[i].vel=(Vector3){0};return;}
    EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];EnemyType t=p->enemyType;
    if(t==ENEMY_TYPE_FOUNDRY && !s->supplied) {
        int stocked=0;
        for(int n=0;n<16 && stocked<2;n++) {
            float angle=(n%8)*PI/4,range=n<8?4:6;
            Vector3 at={floorf(p->pos.x+cosf(angle)*range),0,floorf(p->pos.z+sinf(angle)*range)};
            if(fabsf(at.x)>FLOOR_SIZE-2 || fabsf(at.z)>FLOOR_SIZE-2)continue;
            bool clear=true;
            for(int x=0;x<2;x++)for(int z=0;z<2;z++)for(int y=0;y<2;y++)
                if(table_get((int)(at.x*2)+x,y,(int)(at.z*2)+z)>=0)clear=false;
            for(int j=0;j<firefightHumans;j++)if(players[j].respawn_timer<=0 && v_length(v_sub(players[j].pos,at))<2)clear=false;
            if(clear){enemy_scrap_bundle(at);stocked++;}
        }
        s->supplied=true;
    }
    s->cooldown-=dt;s->melee-=dt;s->repair-=dt;s->searchTimer-=dt;s->stompCooldown-=dt;s->repositionTimer=fmaxf(0,s->repositionTimer-dt);s->target=enemy_target(i);
    if(s->target<0){enemy_drop_payload(i);p->vel.x=p->vel.z=0;s->reason="no living target";return;}
    Vector3 to=v_sub(players[s->target].pos,p->pos);float distance=v_length(to);
    float yaw=atan2f(-to.x,-to.z)*RAD2DEG,turn=t==ENEMY_TYPE_SHIELDBEARER?120:180;
    p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-turn*dt,turn*dt);
    if(t==ENEMY_TYPE_FOUNDRY)p->isExposed=false;
    if(p->isExposed){enemy_drop_payload(i);s->windup=0;s->state=0;s->selected=0;s->reason="exposed";}
    if(s->windup>0) {
        s->windup-=dt;p->vel.x=p->vel.z=0;s->work=2;s->reason="wind-up";
        if(s->windup<=0) {
            bool done=false;
            if(s->state==1)done=enemy_collect(i);
            if(s->state==2)done=enemy_plant(i);
            if(s->state==3){done=enemy_launch(i);if(!done){enemy_drop_payload(i);s->repositionTimer=2;}}
            if(s->state==4)done=enemy_build(i);
            if(s->state==5){enemy_stomp(i);done=true;}
            s->actions+=done;s->reason=done?"recovering":"target invalid or blocked";
            s->cooldown=done?(t==ENEMY_TYPE_BUILDER?8:t==ENEMY_TYPE_SAPPER?8:t==ENEMY_TYPE_SCRAPPER?1:5):1;
            s->searchTimer=done?.5f:1;s->selected=0;s->state=0;s->work=3;
        }
        return;
    }
    Vector3 move=to;float speed=1;s->work=0;
    if(t==ENEMY_TYPE_SHIELDBEARER) {
        speed=.7f;
        for(int j=MAX_PLAYERS;j<activePlayers;j++)if(j!=i && players[j].respawn_timer<=0 && v_length(v_sub(p->pos,players[j].pos))<8) {
            move=v_sub(v_add(players[j].pos,v_mul(v_norm(v_sub(players[s->target].pos,players[j].pos)),2)),p->pos);break;
        }
    }
    if(t==ENEMY_TYPE_SCRAPPER && !p->isExposed && s->equipmentCount<12) {
        if(s->searchTimer<=0 && s->cooldown<=0){s->selected=enemy_find_scrap(i);s->searchTimer=s->selected?.5f:1;s->reason=s->selected?"seeking scrap":"no settled scrap";}
        int v=enemy_voxel(s->selected);
        if(v>=0 && enemy_scrap_eligible(v,i)) {
            move=v_sub(voxels[v].pos,p->pos);s->work=1;
            if(v_length(move)<=2.5f && s->cooldown<=0){s->state=1;s->windup=1;enemy_cue(i);}
        }
    }
    if(t==ENEMY_TYPE_SAPPER) {
        if(enemyCharges[i].active)move=v_sub(p->pos,enemyCharges[i].pos);
        else if(!p->isExposed && s->cooldown<=0) {
            if(s->searchTimer<=0){s->selected=enemy_sapper_target(i,s->target);s->searchTimer=s->selected?.5f:1;}
            int v=enemy_voxel(s->selected);
            if(v>=0 && !voxels[v].simulate) {move=v_sub(voxels[v].pos,p->pos);s->work=1;if(v_length(move)<=2.5f){s->state=2;s->windup=1.5f;enemy_cue(i);}}
        }
    }
    if(t==ENEMY_TYPE_BUILDER) {
        speed=.85f;if(distance<6)move=v_mul(to,-1);else if(distance<10)move=(Vector3){0};
        if(!p->isExposed && s->cooldown<=0) {
            if(s->searchTimer<=0) {
                s->buildKind=players[s->target].pos.y-p->pos.y>.5f?1:0;
                s->work=enemy_find_build_site(i,to)?1:0;s->searchTimer=s->work?.5f:1;
                s->reason=s->work?"seeking build site":"no clear build site";
                if(s->work)s->state=4;else s->state=0;
            }
            if(s->state==4) {
                Vector3 d=v_sub(s->buildOrigin,p->pos);d.y=0;
                if(v_length(d)>4.1f){move=d;s->work=1;}
                else if(enemy_build_valid(i,s->buildOrigin,s->buildKind)){s->windup=2;enemy_cue(i);}
                else {s->state=0;s->searchTimer=1;}
            }
        }
    }
    if(t==ENEMY_TYPE_LOBBER || t==ENEMY_TYPE_FOUNDRY) {
        speed=t==ENEMY_TYPE_FOUNDRY?.35f:.75f;
        if(distance<10)move=v_mul(to,-1);else if(distance<24)move=(Vector3){0};
        if(s->repositionTimer>0){move=distance>24?to:(Vector3){to.z,0,-to.x};s->reason="repositioning for clear arc";}
        if(!p->isExposed && s->cooldown<=0 && s->repositionTimer<=0) {
            if(s->searchTimer<=0 && !p->tetherHolding){s->selected=enemy_find_payload(i);s->searchTimer=s->selected?.5f:1;s->reason=s->selected?"seeking payload":"no throwable group";}
            int v=enemy_voxel(s->selected);
            if(v>=0 && !p->tetherHolding) {
                move=v_sub(voxels[v].pos,p->pos);s->work=1;
                if(v_length(move)<=5 && enemy_clear_line(v_add(p->pos,(Vector3){0,1,0}),v_add(voxels[v].pos,(Vector3){0,1,0}))) {
                    int count;Vector3 center,extent;
                    if(enemy_payload(i,v,&count,&center,&extent)) {p->goldTetherCharged=true;acquire_tether_hit(i,v,(Vector3){0,1,0});}
                }
            }
            if(p->tetherHolding) {
                s->aim=players[s->target].pos;Vector3 velocity;
                if(enemy_lob_velocity(player_tether_hold_position(p),s->aim,&velocity)) {s->state=3;s->windup=1.5f;enemy_cue(i);}
                else {enemy_drop_payload(i);s->selected=0;s->searchTimer=s->cooldown=1;s->reason="blocked arc";s->repositionTimer=2;move=distance>24?to:(Vector3){to.z,0,-to.x};}
            }
        }
        if(t==ENEMY_TYPE_FOUNDRY && distance<3 && s->stompCooldown<=0 && s->windup<=0) {s->state=5;s->windup=1.25f;s->stompCooldown=6;enemy_cue(i);}
    }
    if(t==ENEMY_TYPE_ANCHOR) {
        speed=.65f;if(!s->equipmentCount && s->disabled<=0 && !p->isExposed)enemy_plate(i,(Vector3){0,.9f,.65f},0);
        if(distance<7)move=v_mul(to,-1);else if(distance<12)move=(Vector3){0};
    }
    if(t==ENEMY_TYPE_FOUNDRY) {
        if(s->repair<=0 && s->repairs<12 && s->equipmentCount<ENEMY_EQUIPMENT_MAX && !p->isExposed) {
            uint64_t work=s->selected;s->selected=enemy_find_scrap(i);if(enemy_collect(i))s->repairs++;s->selected=work;s->repair=3;
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
    if(distance<1.4f && enemy_clear_line(p->pos,players[s->target].pos) && s->melee<=0 && s->windup<=0 && t!=ENEMY_TYPE_FOUNDRY) {
        if(players[s->target].isExposed)kill_player(s->target,i,true,false);else apply_matter_damage(s->target,i,10);s->melee=1.5f;
    }
    move.y=0;
    if(s->windup>0 || v_length(move)<.25f){p->vel.x=p->vel.z=0;return;}
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
        if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) {
            DrawCube(p,1.5f,1.5f,1.5f,enemy_core_open(i)?GOLD:DARKGRAY);
            DrawCubeWires(p,1.53f,1.53f,1.53f,enemy_core_open(i)?YELLOW:c);
            for(int n=0;n<s->equipmentCount;n++)if(s->equipment[n].group<4)DrawCubeWires(enemy_local(i,s->equipment[n].offset),.52f,.52f,.52f,ORANGE);
        }
        if(s->windup>0) {
            DrawCircle3D(p,.9f,(Vector3){1,0,0},90,c);
            if(s->state==3) {DrawCircle3D(s->aim,1.5f,(Vector3){1,0,0},90,ORANGE);DrawLine3D(p,s->aim,Fade(c,.5f));}
            if(s->state==5)DrawCircle3D(p,3,(Vector3){1,0,0},90,RED);
            if(s->state==4) DrawCubeWires(v_add(s->buildOrigin,(Vector3){0,1,s->buildKind?1:0}),3,2,s->buildKind?2:.5f,YELLOW);
        }
    }
    for(int i=0;i<MAX_COMBATANTS;i++) if(enemyCharges[i].active) {
        EnemyCharge *c=&enemyCharges[i];Color col=((int)(c->fuse*8)%2)?RED:YELLOW;
        DrawSphere(c->pos,.35f,col);DrawCircle3D(c->pos,4,(Vector3){1,0,0},90,Fade(col,.7f));
    }
    for(int i=0;i<MAX_COMBATANTS*2;i++) if(enemyLobs[i].id) DrawCircle3D(enemyLobs[i].target,1.5f,(Vector3){1,0,0},90,ORANGE);
}
