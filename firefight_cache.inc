static uint32_t cacheVisited[MAX_PLAYERS];

static void fp_build_spawn_glass(void) {
    Vector3 low={0},high={0};
    for(int n=0;n<FP_RULE_COUNT;n++) {
        Vector3 d=world_delta(fpCachePositions[n],fpCacheOrigin);
        low=Vector3Min(low,d); high=Vector3Max(high,d);
    }
    for(int i=1;i<firefightHumans;i++) {
        if(combatantLife[i]!=LIFE_ALIVE)continue;
        Vector3 d=world_delta(players[i].pos,fpCacheOrigin);
        // Separate spawn sites belong to their own part of the map.
        if(v_length(d)<18){low=Vector3Min(low,d);high=Vector3Max(high,d);}
    }
    int xmin=floorf((fpCacheOrigin.x+low.x-2.5f)/VOXEL_SIZE);
    int xmax=ceilf((fpCacheOrigin.x+high.x+2.5f)/VOXEL_SIZE);
    int zmin=floorf((fpCacheOrigin.z+low.z-2.5f)/VOXEL_SIZE);
    int zmax=ceilf((fpCacheOrigin.z+high.z+2.5f)/VOXEL_SIZE);
    for(int z=zmin;z<=zmax;z++)for(int x=xmin;x<=xmax;x++) {
        bool edgeX=x==xmin || x==xmax, edgeZ=z==zmin || z==zmax;
        if(!edgeX && !edgeZ)continue;
        if((edgeX && abs(2*z-zmin-zmax)<=8) || (edgeZ && abs(2*x-xmin-xmax)<=8))continue;
        Vector3 at=world_wrap((Vector3){(x+.5f)*VOXEL_SIZE,fpCacheOrigin.y+2,(z+.5f)*VOXEL_SIZE});
        if(!worldPeriodic && (fabsf(at.x)>=FLOOR_SIZE || fabsf(at.z)>=FLOOR_SIZE))continue;
        VoxelHit floorHit;
        float ground=first_voxel_hit_detailed((Ray){at,{0,-1,0}},at.y,-1,true,&floorHit)?at.y-floorHit.t:0;
        if(ground>fpCacheOrigin.y+.25f)continue;
        int y0=lroundf(ground/VOXEL_SIZE);
        // Leave a clearance buffer beside existing structures and beneath lintels.
        bool clear=true;
        for(int y=y0;y<y0+6 && clear;y++)for(int dz=-1;dz<=1 && clear;dz++)for(int dx=-1;dx<=1;dx++) {
            int id=table_get_static_only(x+dx,y,z+dz);
            if(id>=0 && !spawn_glass(id)){clear=false;break;}
        }
        if(!clear)continue;
        for(int y=y0;y<y0+4;y++) {
            int id=add_static_voxel_at_grid(x,y,z,spawnGlassColor,0);
            if(id>=0){voxels[id].reinforced=true;voxels[id].persistentConstruction=true;}
        }
    }
}
static void fp_cache_place(void){
 fpCacheReady=false;if(!fp_active() || netTransport.role==NET_ROLE_CLIENT)return;
 fpCacheOrigin=players[0].pos;
 for(int r=0;r<FP_RULE_COUNT;r++){
  FfItem stand={.pos=world_wrap(v_add(fpCacheOrigin,(Vector3){((r%4)-1.5f)*2.8f,0,4+(r/4)*2.8f}))};
  fp_ground(&stand);fpCachePositions[r]=stand.pos;
 }
 memset(cacheVisited,0,sizeof(cacheVisited));
 fp_build_spawn_glass();
 fpCacheReady=true;
}
static bool fp_cache_eligible(int i,int type){
 Equipment *e=&equipment[i];int gun=eq_pickup_gun(type);
 if(gun>=0)return e->gun[0]!=gun && e->gun[1]!=gun;
 if(type==PICKUP_ADHESIVE)return !(e->upgrades&EQ_ADHESIVE);
 if(type==PICKUP_REINFORCED)return !(e->upgrades&EQ_BUILD);
 if(type==PICKUP_SHIELD_REGEN)return !(e->upgrades&EQ_SHIELD);
 if(type==PICKUP_GOLD_TETHER)return !players[i].goldTetherCharged;
 int kind=type==PICKUP_BLACK_HOLE?GRENADE_GRAVITY:type==PICKUP_IGLOO?GRENADE_IGLOO:GRENADE_REGULAR;
 return ex_grenade_count(i,kind)<2;
}
static bool fp_cache_reachable(int i,int r,float distance) {
    if(distance>=1.8f)return false;
    Vector3 delta=world_delta(fpCachePositions[r],players[i].pos);
    return distance<.15f || !first_voxel_hit_detailed((Ray){players[i].pos,v_norm(delta)},distance-.15f,-1,false,NULL);
}
static int fp_cache_near(int i){
    if(!fp_active() || !fpCacheReady || i<0 || i>=firefightHumans || players[i].respawn_timer>0 || combatantLife[i]!=LIFE_ALIVE)return -1;
    int best=-1;float nearest=1.8f;
    for(int r=0;r<FP_RULE_COUNT;r++){
        int type=fpRules[r].type;
        if(!(fpUnlocked&(1u<<r)) || eq_pickup_gun(type)<0 || !fp_cache_eligible(i,type))continue;
        float distance=v_length(world_delta(fpCachePositions[r],players[i].pos));
        if(distance<nearest && fp_cache_reachable(i,r,distance)){best=r;nearest=distance;}
    }
    return best;
}
static bool fp_cache_grant(int i,int r){
    if(r<0 || !fp_reserve(fpCount+1))return false;
    size_t n=fpCount++;
    fpItems[n]=(FfItem){.id=fpNextId++,.type=fpRules[r].type,.owner=-1,.slot=-1,.pos=fpCachePositions[r]};
    if(!fp_collect(i,n,true)){fpCount--;return false;}
    return true;
}
static bool fp_cache_collect(int i){
    return netTransport.role!=NET_ROLE_CLIENT && fp_cache_grant(i,fp_cache_near(i));
}
static void fp_cache_walkup(int i){
    if(!fpCacheReady || combatantLife[i]!=LIFE_ALIVE || players[i].respawn_timer>0){
        cacheVisited[i]=0;return;
    }
    for(int r=0;r<FP_RULE_COUNT;r++) {
        unsigned bit=1u<<r;int type=fpRules[r].type;
        if(!(fpUnlocked&bit) || eq_pickup_gun(type)>=0)continue;
        float distance=v_length(world_delta(fpCachePositions[r],players[i].pos));
        if(distance>=1.8f){cacheVisited[i]&=~bit;continue;}
        if((cacheVisited[i]&bit) || !fp_cache_reachable(i,r,distance))continue;
        if(!fp_cache_eligible(i,type) || fp_cache_grant(i,r))cacheVisited[i]|=bit;
    }
}

/* Console skips grant only equipment introduced before the requested wave. */
static void fp_unlock_previous_waves(int wave){
 if(!fp_active() || netTransport.role==NET_ROLE_CLIENT)return;
 for(int r=0;r<FP_RULE_COUNT;r++)if(fpRules[r].wave>0 && fpRules[r].wave<wave)fp_queue_unlock(fpRules[r].type);
}
