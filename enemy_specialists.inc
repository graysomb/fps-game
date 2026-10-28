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
    for(int c=0;c<ENEMY_CHARGE_MAX;c++)
        if(enemyCharges[c].active && enemyCharges[c].id==id) return enemyCharges[c].owner;
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
    for(int c=0;c<ENEMY_CHARGE_MAX;c++)if(enemyCharges[c].active && enemyCharges[c].id==id)enemyCharges[c].fuse=-1;
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
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
static int enemy_live_charge(int owner) {
    for(int c=0;c<ENEMY_CHARGE_MAX;c++)if(enemyCharges[c].active && enemyCharges[c].owner==owner)return c;
    return -1;
}
static bool enemy_create_charge(int owner,Vector3 pos,uint64_t support) {
    int slot=!enemyCharges[owner].active?owner:-1;
    if(slot<0)for(int c=0;c<ENEMY_CHARGE_MAX;c++)if(!enemyCharges[c].active){slot=c;break;}
    if(slot<0)return false;
    int v=addVoxel(pos.x,pos.y,pos.z,false,true,ORANGE,0);if(v<0)return false;
    enemyCharges[slot]=(EnemyCharge){true,owner,voxels[v].identity,support,pos,4,0};
    voxels[v].owner=owner;voxels[v].glueEligible=false;
    if(support)enemy_pin(v,pos,(Vector3){0});
    return true;
}
static void enemy_on_death(int i) {
    if(players[i].enemyType==ENEMY_TYPE_SAPPER)
        enemy_create_charge(i,v_add(players[i].pos,(Vector3){0,.3f,0}),0);
    enemy_release(i);
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY && !enemyLab) {
        int pending=ffRemaining[ENEMY_TYPE_SWARMER];
        ffRemaining[ENEMY_TYPE_SWARMER]=0;
        firefightWaveEnemiesTotal-=pending;firefightWaveEnemiesRemaining-=pending;
    }
}
static void enemy_cleanup(bool removeEquipment) {
    goRemoteActive=false;goRemoteTick=0;
    if(goActive){for(int v=voxel_count-1;v>=0;v--)if(voxels[v].golemOwned)remove_voxel_index(v);golem_physics_reset();}
    tower_reset();
    enemyCleanupPending=false;
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++)enemy_drop_payload(i);
    if(removeEquipment) for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        for(int n=0;n<enemySpecial[i].equipmentCount;n++) {
            int v=enemy_voxel(enemySpecial[i].equipment[n].id);if(v>=0) remove_voxel_index(v);
        }
        memset(&enemySpecial[i],0,sizeof(enemySpecial[i]));
    }
    for(int i=0;i<ENEMY_CHARGE_MAX;i++) {
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
/* Authored equipment layout shared by simulation and isolated menu previews. */
static int enemy_equipment_layout(EnemyType type,EnemyEquipment *out) {
    int count=0;
    if(type==ENEMY_TYPE_SHIELDBEARER)
        for(int x=0;x<4;x++)for(int y=0;y<5;y++)out[count++]=(EnemyEquipment){.offset={(x-1.5f)*.5f,(y-1.5f)*.5f,-1},.group=x};
    if(type==ENEMY_TYPE_ANCHOR)out[count++]=(EnemyEquipment){.offset={0,.9f,.65f},.group=0};
    if(type==ENEMY_TYPE_FOUNDRY){
        for(int g=0;g<4;g++)for(int y=0;y<4;y++)out[count++]=(EnemyEquipment){.offset={g%2?1.25f:-1.25f,-BASE_EYE_HEIGHT+.25f+y*.5f,g/2?1.25f:-1.25f},.group=g};
        for(int x=0;x<4;x++)for(int z=0;z<2;z++)out[count++]=(EnemyEquipment){.offset={(x-1.5f)*.5f,1.6f,(z-.5f)*.5f},.group=4};
    }
    return count;
}
static void enemy_init(int i) {
    EnemySpecialState *s=&enemySpecial[i];memset(s,0,sizeof(*s));
    for(int g=0;g<4;g++)s->shieldHits[g]=6;
    s->generation=++enemyGeneration;s->cooldown=2;s->spawnTimer=15;s->target=-1;s->reason="seeking";
    s->coreMax=s->coreHealth=enemy_matter(players[i].enemyType);
    players[i].matter=players[i].matterMax=enemy_matter(players[i].enemyType);
    EnemyEquipment layout[ENEMY_EQUIPMENT_MAX];int count=enemy_equipment_layout(players[i].enemyType,layout);
    for(int n=0;n<count;n++)enemy_plate(i,layout[n].offset,layout[n].group);
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
        if(count>=2 && !s->coreOpen){s->coreOpen=true;s->coreCue=3;play_sfx_at(SFX_SHIELD_BREAK,players[i].pos,i,2,.7f);}
        s->lostMask=lost;
    }
}
static bool enemy_core_open(int i) {
    return players[i].enemyType!=ENEMY_TYPE_FOUNDRY || enemySpecial[i].coreOpen || enemySpecial[i].lostMask==15;
}
static bool enemy_clear_line(Vector3 a,Vector3 b) {
    Vector3 d=world_delta(b,a);float len=v_length(d);
    return len<.1f || first_voxel_hit((Ray){a,v_mul(d,1/len)},len,-1)<0;
}
static bool enemy_braces(int anchor,int target) {
    if(anchor==target || players[anchor].respawn_timer>0 || players[target].respawn_timer>0 ||
       players[anchor].enemyType!=ENEMY_TYPE_ANCHOR || enemySpecial[anchor].disabled>0 ||
       !enemySpecial[anchor].equipmentCount || combatantTeam[anchor]!=combatantTeam[target]) return false;
    if(v_length(world_delta(players[anchor].pos,players[target].pos))>8 ||
       !enemy_clear_line(players[anchor].pos,players[target].pos)) return false;
    int rank=0;
    for(int j=MAX_PLAYERS;j<target;j++) if(j!=anchor && players[j].respawn_timer<=0 && combatantTeam[j]==combatantTeam[anchor] &&
        v_length(world_delta(players[anchor].pos,players[j].pos))<=8 && enemy_clear_line(players[anchor].pos,players[j].pos)) rank++;
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
    return v_length(world_delta(x->pos,players[i].pos))<=16;
}
static uint64_t enemy_find_scrap(int i) {
    uint64_t result=0;float best=17;
    for(int v=0;v<voxel_count;v++)if(enemy_scrap_eligible(v,i)) {
        float d=v_length(world_delta(voxels[v].pos,players[i].pos));
        if(d<best){best=d;result=voxels[v].identity;}
    }
    return result;
}
static bool enemy_collect(int i) {
    EnemySpecialState *s=&enemySpecial[i];int v=enemy_voxel(s->selected);
    if(v<0 || !enemy_scrap_eligible(v,i) || v_length(world_delta(voxels[v].pos,players[i].pos))>2.5f) return false;
    int n=s->equipmentCount;
    Vector3 offset={(n%3-1)*.5f,(n/6)*.5f-.25f,n%6<3?-.7f:.7f};
    if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) offset=(Vector3){(s->repairs%4-1.5f)*.5f,1.1f,(s->repairs/4-1)*.5f};
    return enemy_attach(i,v,offset,4);
}
static uint64_t enemy_sapper_target(int i,int target) {
    float best=-1;uint64_t found=0;
    Vector3 threat=v_norm(world_delta(ai_target_body(target)->pos,players[i].pos));
    for(int v=0;v<voxel_count;v++) {
        Voxel *x=&voxels[v];if(x->simulate || x->pos.y>players[i].pos.y+1.5f) continue;
        Vector3 d=world_delta(x->pos,players[i].pos);float dist=v_length(d);
        if(dist>10 || dist<.2f || v_dot(v_norm(d),threat)<.2f) continue;
        float score=12-dist+(table_get_static_only(x->gx,x->gy+1,x->gz)>=0?10:0);
        if(score>best) {best=score;found=x->identity;}
    }
    return found;
}
static bool enemy_plant(int i) {
    EnemySpecialState *s=&enemySpecial[i];int support=enemy_voxel(s->selected);
    if(support<0 || voxels[support].simulate || enemy_live_charge(i)>=0) return false;
    Vector3 d=world_delta(players[i].pos,voxels[support].pos);
    if(v_length(d)>2.5f) return false;
    Vector3 p=v_add(voxels[support].pos,v_mul(v_norm(d),.55f));
    return enemy_create_charge(i,p,s->selected);
}
static void enemy_blast_terrain(EnemyCharge c) {
    foundry_traps_blast(c.pos,c.owner);
    /* Excavation activation skips its interior; a blast must activate the radial shell explicitly. */
    static UnitVoxelBuffer shell;unit_voxel_buffer_clear(&shell);
    for(int v=0;v<voxel_count && shell.count<VOXEL_ACTIVATION_UNIT_BUDGET;v++) {
        Voxel *x=&voxels[v];float distance=v_length(world_delta(x->pos,c.pos));
        if(!x->simulate && !x->type && distance>=2 && distance<4)
            unit_voxel_buffer_push(&shell,x->gx,x->gy,x->gz,x->color,x->type,false,v,x->debugClusterTag,c.owner);
    }
    remove_buffered_static_voxels(&shell);emit_unit_voxels_from_units(&shell,false,true,-1);
    init_static_hash();rebuild_voxel_hash();
    for(int v=voxel_count-1;v>=0;v--)if(v_length(world_delta(voxels[v].pos,c.pos))<2 && enemy_equipment_owner(v)<0 && !go_shield_hit(v,0))remove_voxel_index(v);
    int stamp=++particle_sync_stamp;
    for(int v=0;v<voxel_count;v++) {
        Voxel *x=&voxels[v];float distance=v_length(world_delta(x->pos,c.pos));
        if(!x->simulate || distance>=4 || enemy_equipment_owner(v)>=0 || go_shield_hit(v,0))continue;
        if(x->sleeping)wake_sleeping_cluster_in_place(v);
        for(int f=0;f<6;f++)break_face_link(x,f);
        Vector3 direction=world_delta(x->pos,c.pos);direction.y+=1;
        Vector3 impulse=v_mul(v_norm(direction),12*(1-distance/4));
        x->owner=x->activator=c.owner;x->explosionUntil=clarityClock+10;x->sleeping=false;x->vel=v_add(x->vel,impulse);
        for(int n=0;n<8;n++)if(x->particles[n] && x->particles[n]->sync_stamp!=stamp) {
            Particle *part=x->particles[n];part->sync_stamp=stamp;
            if(part->inv_mass>0)part->vel=v_add(part->vel,impulse);
        }
    }
    rebuild_voxel_hash();rebuild_all_voxel_surfaces();
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
        for(int i=0;!initial && i<activePlayers && clear;i++)if(players[i].respawn_timer<=0 && v_length(world_delta(players[i].pos,at))<2)clear=false;
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
    // A two-meter cube can have a wider AABB after rotating in the tether.
    return extent->x<=2.f && extent->y<=2.f && extent->z<=2.f;
}
/* Extract one connected 4x4x4 patch; never activate the surrounding building. */
static int enemy_extract_payload(int owner,int seed) {
    if(seed<0 || voxels[seed].simulate)return seed;
    Voxel origin=voxels[seed];int queue[64][3],head=0,tail=1;bool seen[64]={0};
    int seedCell[3]={origin.gx,origin.gy,origin.gz},low[3];
    // Fit the box inward from a surface rather than wasting half on empty air.
    for(int axis=0;axis<3;axis++){
        int ends[2]={0,0};
        for(int side=0;side<2;side++)for(int step=1;step<=3;step++){
            int cell[3]={origin.gx,origin.gy,origin.gz};cell[axis]+=(side?1:-1)*step;
            int part=table_get_static_only(cell[0],cell[1],cell[2]);
            if(part<0 || voxels[part].type || voxels[part].pendingActivation)break;
            ends[side]=step;
        }
        int min=seedCell[axis]-ends[0],max=seedCell[axis]+ends[1];
        low[axis]=max-min>=3?(int)clampf(seedCell[axis]-1,min,max-3):min;
    }
    queue[0][0]=origin.gx;queue[0][1]=origin.gy;queue[0][2]=origin.gz;
    seen[origin.gx-low[0]+4*(origin.gy-low[1])+16*(origin.gz-low[2])]=true;
    static UnitVoxelBuffer units;unit_voxel_buffer_clear(&units);
    const int delta[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    while(head<tail){
        int x=queue[head][0],y=queue[head][1],z=queue[head++][2];
        int v=table_get_static_only(x,y,z);
        if(v<0 || voxels[v].type || voxels[v].pendingActivation || voxels[v].rest_edge>VOXEL_SIZE+.001f)continue;
        Voxel *part=&voxels[v];
        unit_voxel_buffer_push(&units,x,y,z,part->color,0,false,v,part->debugClusterTag,owner);
        for(int d=0;d<6;d++){
            int nx=x+delta[d][0],ny=y+delta[d][1],nz=z+delta[d][2];
            int ix=nx-low[0],iy=ny-low[1],iz=nz-low[2];
            if(ix<0 || ix>=4 || iy<0 || iy>=4 || iz<0 || iz>=4)continue;
            int key=ix+4*iy+16*iz;if(seen[key])continue;seen[key]=true;
            queue[tail][0]=nx;queue[tail][1]=ny;queue[tail++][2]=nz;
        }
    }
    if(units.count<8 || voxel_count+units.count>=MAX_VOXELS || MAX_PARTICLES-particle_pool_count+free_particle_count<units.count*16)return -1;
    remove_buffered_static_voxels(&units);
    int first=voxel_count,spawned=emit_unit_voxels_from_units(&units,false,true,-1);
    if(!spawned)return -1;
    // Isolate the patch from any neighboring dynamic material glued during creation.
    for(int v=first;v<voxel_count;v++)for(int face=0;face<6;face++)break_face_link(&voxels[v],face);
    rebuild_voxel_hash();
    for(int v=first;v<voxel_count;v++)glue_neighbor_faces_for_voxel(v);
    for(int v=first;v<voxel_count;v++)for(int face=0;face<6;face++)
        if(voxels[v].glued_neighbors[face]<first)break_face_link(&voxels[v],face);
    refresh_activation_surfaces(&units);
    return first;
}
static bool enemy_payload_reachable(int owner,int v) {
    Vector3 delta=world_delta(voxels[v].pos,players[owner].pos);float distance=v_length(delta);
    if(distance>5)return false;
    if(distance<.01f)return true;
    VoxelHit hit;
    return !first_voxel_hit_detailed((Ray){players[owner].pos,v_mul(delta,1/distance)},distance,-1,true,&hit) || hit.id==v;
}
static uint64_t enemy_find_payload(int i) {
    float best=20;uint64_t found=0;
    if(players[i].enemyType==ENEMY_TYPE_LOBBER && enemySpecial[i].payloadFallback<=0)for(int v=0;v<voxel_count;v++){
        Voxel *x=&voxels[v];if(x->simulate || x->type || x->pendingActivation || x->rest_edge>VOXEL_SIZE+.001f)continue;
        float d=v_length(world_delta(x->pos,players[i].pos));
        if(d<best){found=x->identity;best=d;}
    }
    if(found)return found;
    static unsigned visited[MAX_VOXELS],stamp; if(!++stamp){memset(visited,0,sizeof(visited));stamp=1;}
    for(int v=0;v<voxel_count;v++) {
        if(visited[v]==stamp || !voxels[v].simulate || voxels[v].type || voxels[v].isBullet)continue;
        float d=v_length(world_delta(voxels[v].pos,players[i].pos));if(d>=best)continue;
        int members=build_glue_cluster_indices(v,glueClusterIndices);
        for(int k=0;k<members;k++)visited[glueClusterIndices[k]]=stamp;
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
    const float dt=PBD_MAX_STEP_DT/PBD_SUBSTEPS,q=VELOCITY_DAMPING;
    Vector3 d=world_delta(target,from);
    for(int steps=60*PBD_SUBSTEPS;steps<=480*PBD_SUBSTEPS;steps+=15*PBD_SUBSTEPS){
        float response=0,weight=1,gv=0,gy=0;
        for(int k=0;k<steps;k++){weight*=q;response+=weight*dt;gv=q*gv-GRAVITY*dt;gy+=gv*dt;}
        Vector3 candidate={d.x/response,(d.y-gy)/response,d.z/response};
        if(v_length(candidate)>fmaxf(80,4*FLOOR_SIZE))continue;
        Vector3 position=from,vel=candidate,last=from;float apex=from.y;bool blocked=false;
        for(int n=1;n<=steps && !blocked;n++){
            vel=v_mul(vel,q);vel.y-=GRAVITY*dt;position=v_add(position,v_mul(vel,dt));apex=fmaxf(apex,position.y);
            if(n%4 && n!=steps)continue;
            Vector3 segment=v_sub(position,last);float length=v_length(segment);
            if(!worldPeriodic && (fabsf(position.x)>FLOOR_SIZE-.5f || fabsf(position.z)>FLOOR_SIZE-.5f)){blocked=true;break;}
            for(int corner=0;corner<9 && !blocked;corner++){
                Vector3 offset=corner==8?(Vector3){0}:(Vector3){corner_signs[corner][0]*extent.x,corner_signs[corner][1]*extent.y,corner_signs[corner][2]*extent.z};
                VoxelHit hit;if(length>.0001f && first_voxel_hit_detailed((Ray){v_add(last,offset),v_norm(segment)},length,-1,true,&hit))blocked=true;
            }
            last=position;
        }
        if(!blocked && apex>=fmaxf(from.y,target.y)+4){*velocity=candidate;return true;}
    }
    return false;
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
            bool prior=damageOriginKnown;Vector3 origin=damageOrigin;damageOriginKnown=true;damageOrigin=v_sub(voxels[voxel].pos,v_mul(v_norm(voxels[voxel].vel),2));
            if(players[victim].isExposed)kill_player(victim,l->owner,true,true);else apply_matter_damage(victim,l->owner,30);
            damageOriginKnown=prior;damageOrigin=origin;
        }
        return true;
    }
    return false;
}
static void enemy_stomp(int i) {
    hologram_blast(players[i].pos,3,20);
    for(int j=0;j<firefightHumans;j++)if(players[j].respawn_timer<=0 && v_length(world_delta(players[j].pos,players[i].pos))<=3) {
        apply_matter_damage(j,i,20);
        Vector3 away=world_delta(players[j].pos,players[i].pos);away.y=.6f;
        players[j].vel=v_add(players[j].vel,v_mul(v_norm(away),6));players[j].onGround=false;players[j].meleeKnockbackActive=true;
    }
    play_sfx_at(SFX_EXPLOSION,players[i].pos,i,3,.7f);
}
static int enemy_target(int i) {return aiCombat[i].target;}

static void enemy_tick(float dt) {
    if(netTransport.role==NET_ROLE_CLIENT)return;
    if(enemyCleanupPending) enemy_cleanup(true);
    if(enemyLab && enemyLabPaused && !enemyLabStep) return;
    tower_tick(dt);
    if(enemyLab && players[0].respawn_timer>0) {
        players[0].respawn_timer-=dt;
        if(players[0].respawn_timer<=0) {ff_restore_human(0,100);players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};}
    }
    for(int i=MAX_PLAYERS;i<activePlayers;i++) if(enemy_is_special(i) && players[i].respawn_timer<=0) {
        enemySpecial[i].exposure=fmaxf(0,enemySpecial[i].exposure-dt);
        enemySpecial[i].disabled=fmaxf(0,enemySpecial[i].disabled-dt);
        enemySpecial[i].coreCue=fmaxf(0,enemySpecial[i].coreCue-dt);
        enemySpecial[i].blockedCue=fmaxf(0,enemySpecial[i].blockedCue-dt);
        enemy_update_equipment(i);
    }
    for(int i=0;i<ENEMY_CHARGE_MAX;i++) {
        EnemyCharge *c=&enemyCharges[i];if(!c->active) continue;
        int v=enemy_voxel(c->id),support=enemy_voxel(c->support);
        bool disarmed=v<0 || (c->support && (support<0 || voxels[support].simulate)) || c->fuse<0;
        if(!c->support && v>=0)c->pos=voxels[v].pos;
        c->fuse-=dt;c->beep-=dt;
        if(disarmed || c->fuse<=0) {
            EnemyCharge copy=*c;c->active=false;c->id=0;
            if(v>=0) remove_voxel_index(v);
            if(!disarmed) enemy_blast(copy);
        } else {
            if(c->support)enemy_pin(v,c->pos,(Vector3){0});
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
    if(players[i].enemyType==ENEMY_TYPE_BUILDER){tower_builder(i,dt);return;}
    if(netTransport.role==NET_ROLE_CLIENT)return;
    if(enemyLab && enemyLabPaused && !enemyLabStep){players[i].vel=(Vector3){0};return;}
    ai_tick(i,dt);
    EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];EnemyType t=p->enemyType;
    if(t==ENEMY_TYPE_FOUNDRY)p->isExposed=false;
    if(p->isExposed){enemy_drop_payload(i);s->windup=0;s->state=0;s->selected=0;s->reason="exposed";}
    if(clarity_evade(i,dt))return;
    if(t==ENEMY_TYPE_FOUNDRY && !s->supplied) {
        int stocked=0;
        for(int n=0;n<16 && stocked<2;n++) {
            float angle=(n%8)*PI/4,range=n<8?4:6;
            Vector3 at={floorf(p->pos.x+cosf(angle)*range),0,floorf(p->pos.z+sinf(angle)*range)};
            if(!worldPeriodic && (fabsf(at.x)>FLOOR_SIZE-2 || fabsf(at.z)>FLOOR_SIZE-2))continue;
            bool clear=true;
            for(int x=0;x<2;x++)for(int z=0;z<2;z++)for(int y=0;y<2;y++)
                if(table_get((int)(at.x*2)+x,y,(int)(at.z*2)+z)>=0)clear=false;
            for(int j=0;j<firefightHumans;j++)if(players[j].respawn_timer<=0 && v_length(world_delta(players[j].pos,at))<2)clear=false;
            if(clear){enemy_scrap_bundle(at);stocked++;}
        }
        s->supplied=true;
    }
    s->payloadFallback=fmaxf(0,s->payloadFallback-dt);s->cooldown-=dt;s->melee-=dt;s->repair-=dt;s->searchTimer-=dt;s->stompCooldown-=dt;s->repositionTimer=fmaxf(0,s->repositionTimer-dt);s->target=enemy_target(i);
    // Gathering nearby ammunition does not require sight of a player. Keep this
    // preparation separate from the visible-target requirement for a throw.
    if(!aiCombat[i].visible && !p->isExposed && !p->tetherHolding && s->windup<=0 &&
       s->cooldown<=0 && (t==ENEMY_TYPE_LOBBER || t==ENEMY_TYPE_FOUNDRY)) {
        if(s->searchTimer<=0){s->selected=enemy_find_payload(i);s->searchTimer=1;}
        int v=enemy_voxel(s->selected);
        if(v>=0) {
            Vector3 gather=world_delta(voxels[v].pos,p->pos);
            if(v_length(gather)<=5 && enemy_clear_line(v_add(p->pos,(Vector3){0,1,0}),v_add(voxels[v].pos,(Vector3){0,1,0}))) {
                int count;Vector3 center,extent;
                if(enemy_payload(i,v,&count,&center,&extent)){p->goldTetherCharged=true;acquire_tether_hit(i,v,(Vector3){0,1,0});}
            }
            if(!p->tetherHolding){
                s->reason="seeking payload";
                float yaw=atan2f(-gather.x,-gather.z)*RAD2DEG;
                p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-180*dt,180*dt);
                bot_apply_move(p,&botStates[i],gather,MOVE_SPEED*(t==ENEMY_TYPE_FOUNDRY?.35f:.75f),dt);return;
            }
        }
    }
    if(s->target<0){
        s->windup=0;s->state=0;s->reason="patrolling";
        BotState *b=&botStates[i];b->pathTimer-=dt;
        if(b->pathTimer<=0 || aiCombat[i].stuck>.6f){
            // Search the arena, never an unseen player's current coordinates.
            float radius=fminf(FLOOR_SIZE*.6f,24)*ai_random(i),angle=ai_random(i)*2*PI;
            b->moveTarget=(Vector3){cosf(angle)*radius,p->pos.y,sinf(angle)*radius};b->pathTimer=4;
        }
        Vector3 search=world_delta(b->moveTarget,p->pos);search.y=0;
        float yaw=atan2f(-search.x,-search.z)*RAD2DEG;
        p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-120*dt,120*dt);
        if(t!=ENEMY_TYPE_FOUNDRY)bot_apply_move(p,b,search,MOVE_SPEED*.7f,dt);
        else p->vel.x=p->vel.z=0;
        return;
    }
    Vector3 to=world_delta(aiCombat[i].known,p->pos);float distance=v_length(to);
    float yaw=atan2f(-to.x,-to.z)*RAD2DEG,turn=t==ENEMY_TYPE_SHIELDBEARER?120:180;
    p->yaw+=clampf(angle_diff_deg(p->yaw,yaw),-turn*dt,turn*dt);
    if(s->windup>0) {
        s->windup-=dt;p->vel.x=p->vel.z=0;s->work=2;s->reason="wind-up";
        if(s->windup<=0) {
            bool done=false;
            if(s->state==1)done=enemy_collect(i);
            if(s->state==2)done=enemy_plant(i);
            if(s->state==3){done=enemy_launch(i);if(!done){enemy_drop_payload(i);s->repositionTimer=2;}}

            if(s->state==5){enemy_stomp(i);done=true;}
            s->actions+=done;s->reason=done?"recovering":"target invalid or blocked";
            s->cooldown=done?(t==ENEMY_TYPE_BUILDER?8:t==ENEMY_TYPE_SAPPER?8:t==ENEMY_TYPE_SCRAPPER?1:5):1;
            s->searchTimer=done?.5f:1;s->selected=0;s->state=0;s->work=3;
        }
        return;
    }
    if(!aiCombat[i].visible) {
        Vector3 investigate=to;investigate.y=0;
        if(v_length(investigate)>.5f)bot_apply_move(p,&botStates[i],investigate,MOVE_SPEED*.7f,dt);
        else p->vel.x=p->vel.z=0;
        s->reason="investigating last known position";return;
    }
    Vector3 move=to;float speed=1;s->work=0;
    if(t==ENEMY_TYPE_SHIELDBEARER) {
        speed=.7f;
        for(int j=MAX_PLAYERS;j<activePlayers;j++)if(j!=i && players[j].respawn_timer<=0 && v_length(world_delta(p->pos,players[j].pos))<8) {
            move=world_delta(v_add(players[j].pos,v_mul(v_norm(world_delta(ai_target_body(s->target)->pos,players[j].pos)),2)),p->pos);break;
        }
    }
    if(t==ENEMY_TYPE_SCRAPPER && !p->isExposed && s->equipmentCount<12) {
        if(s->searchTimer<=0 && s->cooldown<=0){s->selected=enemy_find_scrap(i);s->searchTimer=s->selected?.5f:1;s->reason=s->selected?"seeking scrap":"no settled scrap";}
        int v=enemy_voxel(s->selected);
        if(v>=0 && enemy_scrap_eligible(v,i)) {
            move=world_delta(voxels[v].pos,p->pos);s->work=1;
            if(v_length(move)<=2.5f && s->cooldown<=0){s->state=1;s->windup=1;enemy_cue(i);}
        }
    }
    if(t==ENEMY_TYPE_SAPPER) {
        int charge=enemy_live_charge(i);
        if(charge>=0)move=world_delta(p->pos,enemyCharges[charge].pos);
        else if(!p->isExposed && s->cooldown<=0) {
            if(s->searchTimer<=0){s->selected=enemy_sapper_target(i,s->target);s->searchTimer=s->selected?.5f:1;}
            int v=enemy_voxel(s->selected);
            if(v>=0 && !voxels[v].simulate) {move=world_delta(voxels[v].pos,p->pos);s->work=1;if(v_length(move)<=2.5f){s->state=2;s->windup=1.5f;enemy_cue(i);}}
        }
    }
    if(t==ENEMY_TYPE_LOBBER || t==ENEMY_TYPE_FOUNDRY) {
        speed=t==ENEMY_TYPE_FOUNDRY?.35f:.75f;
        if(distance<10)move=v_mul(to,-1);else if(distance<24)move=(Vector3){0};
        if(s->repositionTimer>0){move=distance>24?to:(Vector3){to.z,0,-to.x};s->reason="repositioning for clear arc";}
        if(!p->isExposed && s->cooldown<=0 && s->repositionTimer<=0) {
            if(s->searchTimer<=0 && !p->tetherHolding){s->selected=enemy_find_payload(i);s->searchTimer=1;s->reason=s->selected?"seeking payload":"no throwable group";}
            int v=enemy_voxel(s->selected);
            if(v>=0 && !p->tetherHolding) {
                move=world_delta(voxels[v].pos,p->pos);s->work=1;
                if(enemy_payload_reachable(i,v)) {
                    int count;Vector3 center,extent;
                    if(!voxels[v].simulate){v=enemy_extract_payload(i,v);if(v<0){s->payloadFallback=3;s->selected=0;s->searchTimer=s->cooldown=1;}}
                    if(enemy_payload(i,v,&count,&center,&extent)) {p->goldTetherCharged=true;acquire_tether_hit(i,v,(Vector3){0,1,0});}
                }
            }
            if(p->tetherHolding) {
                s->aim=ai_target_body(s->target)->pos;Vector3 velocity,center,extent;int count;
                if(resolve_tether_voxel(p) && enemy_payload(i,p->tetherVoxel,&count,&center,&extent) && enemy_lob_arc(player_tether_hold_position(p),s->aim,extent,&velocity)) {s->state=3;s->windup=1.5f;enemy_cue(i);}
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
                clarity_spawn(j);
                if(!enemyLab) {ffRemaining[ENEMY_TYPE_SWARMER]--;firefightEnemiesSpawned++;}
                s->spawned++;added++;
            }
            s->spawnTimer=added?15:1;
        }
    }
    if(distance<1.4f && enemy_clear_line(p->pos,ai_target_body(s->target)->pos) && s->melee<=0 && s->windup<=0 && t!=ENEMY_TYPE_FOUNDRY) {
        if(target_is_hologram(s->target))hologram_damage(s->target,10);else if(players[s->target].isExposed)kill_player(s->target,i,true,false);else apply_matter_damage(s->target,i,10);s->melee=1.5f;
    }
    move.y=0;
    if(s->windup>0 || v_length(move)<.25f){p->vel.x=p->vel.z=0;return;}
    bot_apply_move(p,&botStates[i],move,MOVE_SPEED*speed,dt);
}
/* Shared silhouette details, independent of simulation state and usable in menus. */
static void enemy_draw_details(EnemyType type,Vector3 p,Color c) {
        if(type==ENEMY_TYPE_SCRAPPER) DrawCube(v_add(p,(Vector3){0,.1f,.5f}),.8f,.8f,.4f,c);
        if(type==ENEMY_TYPE_SAPPER) {DrawCube(v_add(p,(Vector3){0,.4f,0}),.7f,.3f,.7f,c);DrawSphere(v_add(p,(Vector3){0,.7f,0}),.15f,ORANGE);}
        if(type==ENEMY_TYPE_BUILDER) DrawCube(v_add(p,(Vector3){0,.4f,0}),1,.2f,1,c);
        if(type==ENEMY_TYPE_LOBBER) DrawCylinder(v_add(p,(Vector3){0,.25f,0}),.4f,.3f,1,8,c);
}
static void enemy_draw(void) {
    tower_draw();
    for(int i=MAX_PLAYERS;i<activePlayers;i++) if(enemy_is_special(i) && players[i].respawn_timer<=0) {
        if(!periodic_visible(players[i].pos,5))continue;
        EnemySpecialState *s=&enemySpecial[i];Vector3 p=players[i].pos;Color c=enemy_color(players[i].enemyType);
        /* Distinct silhouettes in addition to color: pack, mast, mortar, hardhat. */
        enemy_draw_details(players[i].enemyType,p,c);
        if(players[i].enemyType==ENEMY_TYPE_ANCHOR) {
            DrawLine3D(p,v_add(p,(Vector3){0,1.1f,0}),c);
            for(int j=MAX_PLAYERS;j<activePlayers;j++) if(enemy_braces(i,j)) {
                Vector3 start=v_add(p,(Vector3){0,.9f,0}),end=world_image(players[j].pos,p);DrawLine3D(start,end,c);
                if(enemySpecial[j].healingPulse>0)DrawCube(Vector3Lerp(start,end,fmodf(clarityClock*1.5f,1)),.12f,.12f,.12f,WHITE);
            }
        }
        if(players[i].enemyType==ENEMY_TYPE_FOUNDRY) {
            DrawCube(p,1.5f,1.5f,1.5f,enemy_core_open(i)?GOLD:DARKGRAY);
            DrawCubeWires(p,1.53f,1.53f,1.53f,enemy_core_open(i)?YELLOW:c);
            for(int n=0;n<s->equipmentCount;n++)if(s->equipment[n].group<4)DrawCubeWires(enemy_local(i,s->equipment[n].offset),.52f,.52f,.52f,ORANGE);
        }
        if(s->windup>0) {
            DrawCircle3D(p,.9f,(Vector3){1,0,0},90,c);
            if(s->state==3) {Vector3 aim=world_image(s->aim,p);DrawCircle3D(aim,1.5f,(Vector3){1,0,0},90,ORANGE);DrawLine3D(p,aim,Fade(c,.5f));}
            if(s->state==5)DrawCircle3D(p,3,(Vector3){1,0,0},90,RED);
            if(s->state==4) DrawCubeWires(v_add(s->buildOrigin,(Vector3){0,1,s->buildKind?1:0}),3,2,s->buildKind?2:.5f,YELLOW);
        }
    }
    for(int i=0;i<ENEMY_CHARGE_MAX;i++) if(enemyCharges[i].active) {
        EnemyCharge *c=&enemyCharges[i];Color col=((int)(c->fuse*8)%2)?RED:YELLOW;
        DrawSphere(c->pos,.35f,col);DrawCircle3D(c->pos,4,(Vector3){1,0,0},90,Fade(col,.7f));
    }
    for(int i=0;i<MAX_COMBATANTS*2;i++) if(enemyLobs[i].id) DrawCircle3D(enemyLobs[i].target,1.5f,(Vector3){1,0,0},90,ORANGE);
}
