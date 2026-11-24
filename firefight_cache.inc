static uint32_t cacheVisited[MAX_PLAYERS];

enum { FP_GATE_CELLS=8, FP_WALL_CELLS=4, FP_LAYOUT_MARGIN=8 };
typedef struct { int x0,x1,z0,z1; } FpRect;
typedef struct { int ground; bool walkable,keepOpen; } FpSite;
typedef struct { FpRect area; int width; FpSite *sites; int *seen,*queue,stamp; } FpSurvey;

static int fp_site_index(const FpSurvey *s,int x,int z){
    return (z-s->area.z0)*s->width+x-s->area.x0;
}
static FpSite fp_site(const FpSurvey *s,int x,int z){return s->sites[fp_site_index(s,x,z)];}
static bool fp_gate(int cell,int min,int max){
    int first=min+(max-min+1-FP_GATE_CELLS)/2;
    return cell>=first && cell<first+FP_GATE_CELLS;
}

// Survey before inserting any glass. A wall or tower is an obstruction, not a hill
// that should lift the whole perimeter. Low ledges count as walkable ground.
static FpSite fp_survey_column(int x,int z,int foot){
    int y=foot;
    while(y>0 && table_get_static_only(x,y-1,z)<0)y--;
    while(y<=foot+1 && table_get_static_only(x,y,z)>=0)y++;
    FpSite site={.ground=y,.walkable=y<=foot+1};
    if(!site.walkable){site.ground=0;return site;}
    for(int k=0;k<FP_WALL_CELLS;k++)if(table_get_static_only(x,y+k,z)>=0)site.walkable=false;
    // Preserve headroom through existing lintels when a perimeter crosses a doorway.
    for(int k=FP_WALL_CELLS;k<=FP_WALL_CELLS+2;k++)
        if(site.walkable && table_get_static_only(x,y+k,z)>=0)site.keepOpen=true;
    return site;
}

static int fp_lane_obstructions(const FpSurvey *s,FpRect r){
    int blocked=0;
    for(int axis=0;axis<2;axis++){
        int lo=axis?r.x0:r.z0,hi=axis?r.x1:r.z1;
        for(int t=lo;t<=hi;t++)for(int side=0;side<2;side++){
            int wall=axis?(side?r.z1:r.z0):(side?r.x1:r.x0);
            FpSite site=fp_site(s,axis?t:wall,axis?wall:t);
            if(!fp_gate(t,lo,hi)){blocked+=site.keepOpen;continue;}
            int previous=site.ground;
            for(int d=-2;d<=2;d++){
                FpSite lane=fp_site(s,axis?t:wall+d,axis?wall+d:t);
                if(!lane.walkable || (d>-2 && abs(lane.ground-previous)>1)){blocked++;break;}
                previous=lane.ground;
            }
        }
    }
    return blocked;
}

// Flood only the interior. Every stand must have a reachable approach within
// its pickup radius, and every entrance must connect back to the spawn.
static int fp_route_obstructions(FpSurvey *s,FpRect r){
    int x=(int)floorf(fpCacheOrigin.x/VOXEL_SIZE),z=(int)floorf(fpCacheOrigin.z/VOXEL_SIZE);
    int start=fp_site_index(s,x,z),head=0,tail=0;
    s->stamp++;
    if(x>r.x0 && x<r.x1 && z>r.z0 && z<r.z1 && s->sites[start].walkable){
        s->queue[tail++]=start;s->seen[start]=s->stamp;
    }
    const int steps[]={-1,1,-s->width,s->width};
    while(head<tail){
        int at=s->queue[head++];
        for(int n=0;n<4;n++){
            int next=at+steps[n];
            x=s->area.x0+next%s->width;z=s->area.z0+next/s->width;
            if(x<=r.x0||x>=r.x1||z<=r.z0||z>=r.z1||s->seen[next]==s->stamp)continue;
            if(!s->sites[next].walkable||abs(s->sites[next].ground-s->sites[at].ground)>1)continue;
            s->seen[next]=s->stamp;s->queue[tail++]=next;
        }
    }
    int blocked=0;
    for(int n=0;n<FP_RULE_COUNT;n++){
        Vector3 stand=v_add(fpCacheOrigin,world_delta(fpCachePositions[n],fpCacheOrigin));
        x=(int)floorf(stand.x/VOXEL_SIZE);z=(int)floorf(stand.z/VOXEL_SIZE);
        bool reached=false;
        for(int dz=-1;dz<=1;dz++)for(int dx=-1;dx<=1;dx++){
            int px=x+dx,pz=z+dz;
            if(px<=r.x0||px>=r.x1||pz<=r.z0||pz>=r.z1)continue;
            int i=fp_site_index(s,px,pz);
            if(s->seen[i]==s->stamp && fabsf(stand.y-(s->sites[i].ground*VOXEL_SIZE+BASE_EYE_HEIGHT))<1.5f)reached=true;
        }
        blocked+=!reached;
    }
    int mx=r.x0+(r.x1-r.x0+1)/2,mz=r.z0+(r.z1-r.z0+1)/2;
    int exits[]={fp_site_index(s,r.x0+1,mz),fp_site_index(s,r.x1-1,mz),
                 fp_site_index(s,mx,r.z0+1),fp_site_index(s,mx,r.z1-1)};
    for(int n=0;n<4;n++)blocked+=s->seen[exits[n]]!=s->stamp;
    return blocked;
}

static FpRect fp_foundry_deck(void){
    return (FpRect){floorf(-FOUNDRY_DECK_X*1.875f),ceilf((FOUNDRY_DECK_X+1)*1.875f)-1,
                    floorf(-FOUNDRY_DECK_Z*1.875f),ceilf((FOUNDRY_DECK_Z+1)*1.875f)-1};
}
static Vector3 fp_base_center(void){
    Vector3 center=fpCacheOrigin;
    if(firefightFoundry && !useCustomMap){
        FpRect deck=fp_foundry_deck();
        center.x=(deck.x0+deck.x1+1)*VOXEL_SIZE*.5f;
        center.z=(deck.z0+deck.z1+1)*VOXEL_SIZE*.5f;
    }
    return center;
}

static void fp_build_spawn_glass(void){
    Vector3 low=fpCacheOrigin,high=low;
    for(int n=0;n<FP_RULE_COUNT;n++){
        Vector3 at=v_add(fpCacheOrigin,world_delta(fpCachePositions[n],fpCacheOrigin));
        low=Vector3Min(low,at);high=Vector3Max(high,at);
    }
    for(int i=1;i<firefightHumans;i++)if(combatantLife[i]==LIFE_ALIVE){
        Vector3 d=world_delta(players[i].pos,fpCacheOrigin);
        if(v_length(d)<18){Vector3 at=v_add(fpCacheOrigin,d);low=Vector3Min(low,at);high=Vector3Max(high,at);}
    }
    int limit=(int)floorf(FLOOR_SIZE/VOXEL_SIZE);
    Vector3 center=fp_base_center();
    if(firefightFoundry && !useCustomMap){
        FpRect deck=fp_foundry_deck();
        low.x=fminf(low.x,deck.x0*VOXEL_SIZE);high.x=fmaxf(high.x,(deck.x1+1)*VOXEL_SIZE);
        low.z=fminf(low.z,deck.z0*VOXEL_SIZE);high.z=fmaxf(high.z,(deck.z1+1)*VOXEL_SIZE);
    }
    int cx=lroundf(center.x/VOXEL_SIZE),cz=lroundf(center.z/VOXEL_SIZE);
    int xmin=floorf((low.x-2.5f)/VOXEL_SIZE)-1,xend=ceilf((high.x+2.5f)/VOXEL_SIZE)+1;
    int zmin=floorf((low.z-2.5f)/VOXEL_SIZE)-1,zend=ceilf((high.z+2.5f)/VOXEL_SIZE)+1;
    if(!worldPeriodic){
        xmin=clampi(xmin,-limit+3,limit-3);xend=clampi(xend,-limit+3,limit-3);
        zmin=clampi(zmin,-limit+3,limit-3);zend=clampi(zend,-limit+3,limit-3);
        // Move the center only as far as needed to fit the footprint and a gate approach.
        cx=clampi(cx,ceilf((xend-limit+3)*.5f),floorf((xmin+limit-3)*.5f));
        cz=clampi(cz,ceilf((zend-limit+3)*.5f),floorf((zmin+limit-3)*.5f));
    }
    int rx=(int)fmaxf(cx-xmin,xend-cx),rz=(int)fmaxf(cz-zmin,zend-cz);
    FpRect base={cx-rx,cx+rx-1,cz-rz,cz+rz-1};
    FpSurvey s={.area={base.x0-FP_LAYOUT_MARGIN-2,base.x1+FP_LAYOUT_MARGIN+2,
                       base.z0-FP_LAYOUT_MARGIN-2,base.z1+FP_LAYOUT_MARGIN+2}};
    s.width=s.area.x1-s.area.x0+1;
    int count=s.width*(s.area.z1-s.area.z0+1);
    s.sites=calloc(count,sizeof(*s.sites));s.seen=calloc(count,sizeof(int));s.queue=malloc(count*sizeof(int));
    if(!s.sites||!s.seen||!s.queue){TraceLog(LOG_ERROR,"Spawn shield survey allocation failed");goto done;}
    int foot=(int)floorf((fpCacheOrigin.y-BASE_EYE_HEIGHT)/VOXEL_SIZE);
    for(int z=s.area.z0;z<=s.area.z1;z++)for(int x=s.area.x0;x<=s.area.x1;x++)
        s.sites[fp_site_index(&s,x,z)]=fp_survey_column(x,z,clampi(foot,0,511));
    if(useCustomMap && cmRuntime)for(int n=0;n<cmRuntime->spawnCount;n++){
        Vector3 at=v_add(fpCacheOrigin,world_delta(cmRuntime->spawns[n].pos,fpCacheOrigin));
        int sx=floorf(at.x/VOXEL_SIZE),sz=floorf(at.z/VOXEL_SIZE);
        for(int z=sz-1;z<=sz+1;z++)for(int x=sx-1;x<=sx+1;x++)
            if(x>=s.area.x0&&x<=s.area.x1&&z>=s.area.z0&&z<=s.area.z1)s.sites[fp_site_index(&s,x,z)].keepOpen=true;
    }
    FpRect best=base;int bestScore=INT_MAX,bestExtra=INT_MAX;
    // Minimize blocked access first, then total outward displacement. Even side
    // lengths put an eight-cell gate exactly on each side's midpoint.
    for(int dx=0;dx<=FP_LAYOUT_MARGIN;dx++)for(int dz=0;dz<=FP_LAYOUT_MARGIN;dz++){
        int extra=dx+dz;if(bestScore==0 && extra>=bestExtra)continue;
        FpRect r={base.x0-dx,base.x1+dx,base.z0-dz,base.z1+dz};
        if(!worldPeriodic && (r.x0< -limit+3||r.x1>limit-4||r.z0< -limit+3||r.z1>limit-4))continue;
        if(worldPeriodic && (r.x1-r.x0+1>=2*limit||r.z1-r.z0+1>=2*limit))continue;
        int score=fp_lane_obstructions(&s,r);
        if(score>bestScore || (score==bestScore && extra>=bestExtra))continue;
        score+=fp_route_obstructions(&s,r);
        if(score<bestScore || (score==bestScore && extra<bestExtra)){best=r;bestScore=score;bestExtra=extra;}
    }
    if(bestScore)TraceLog(LOG_WARNING,"Firefight spawn shield: obstructed layout on %s at (%.1f, %.1f), score %d",
        useCustomMap&&cmRuntime?cmRuntime->name:"arena",fpCacheOrigin.x,fpCacheOrigin.z,bestScore);
    int top=0;
    for(int z=best.z0;z<=best.z1;z++)for(int x=best.x0;x<=best.x1;x++){
        if(x!=best.x0&&x!=best.x1&&z!=best.z0&&z!=best.z1)continue;
        for(int dz=-1;dz<=1;dz++)for(int dx=-1;dx<=1;dx++){
            FpSite site=fp_site(&s,x+dx,z+dz);
            if(site.walkable&&site.ground>top)top=site.ground;
        }
    }
    top+=FP_WALL_CELLS;
    for(int z=best.z0;z<=best.z1;z++)for(int x=best.x0;x<=best.x1;x++){
        bool edgeX=x==best.x0||x==best.x1,edgeZ=z==best.z0||z==best.z1;
        if(!edgeX&&!edgeZ)continue;
        if((edgeX&&fp_gate(z,best.z0,best.z1))||(edgeZ&&fp_gate(x,best.x0,best.x1)))continue;
        FpSite site=fp_site(&s,x,z);
        if(site.keepOpen)continue;
        for(int y=site.ground;y<top;y++){
            if(table_get_static_only(x,y,z)>=0)continue;
            int id=add_static_voxel_at_grid(x,y,z,spawnGlassColor,0);
            if(id>=0){voxels[id].reinforced=true;voxels[id].persistentConstruction=true;}
        }
    }
done:
    free(s.sites);free(s.seen);free(s.queue);
}
static void fp_cache_place(void){
 fpCacheReady=false;if(!fp_active() || netTransport.role==NET_ROLE_CLIENT)return;
 fpCacheOrigin=players[0].pos;
 Vector3 center=fp_base_center();
 // Fit the Foundry rows between the deck's raised blocks and harvest rail.
 bool onDeck=firefightFoundry && !useCustomMap;
 float columnSpacing=onDeck?2.2f:2.8f,rowSpacing=onDeck?4.2f:2.8f;
 float halfDepth=((FP_RULE_COUNT-1)/4)*rowSpacing*.5f;
 if(!worldPeriodic){
  // Shift the formation together so stands and the outside gate approaches fit.
  float edge=FLOOR_SIZE-4.5f;
  center.x=clampf(center.x,-edge+1.5f*columnSpacing,edge-1.5f*columnSpacing);
  center.z=clampf(center.z,-edge+halfDepth,edge-halfDepth);
 }
 for(int r=0;r<FP_RULE_COUNT;r++){
  Vector3 anchor=world_wrap(v_add(center,(Vector3){((r%4)-1.5f)*columnSpacing,0,(r/4)*rowSpacing-halfDepth}));
  FfItem stand={.pos=anchor};fp_ground(&stand);
  // A downward supply ray can land on top of a tower. Try the surrounding
  // floor before accepting a stand higher than the spawn's walkable surface.
  if(stand.pos.y>fpCacheOrigin.y+.5f)for(int n=0;n<32;n++){
   float angle=n*PI/4,radius=(1+n/8)*.75f;
   Vector3 at=world_wrap(v_add(anchor,(Vector3){radius*cosf(angle),0,radius*sinf(angle)})),safe;
   if(!clarity_supply_position(at,0,&safe)||safe.y>fpCacheOrigin.y+.5f)continue;
   bool separate=true;
   for(int k=0;k<r;k++)if(v_length(world_delta(safe,fpCachePositions[k]))<1.4f)separate=false;
   if(separate){stand.pos=safe;break;}
  }
  fpCachePositions[r]=stand.pos;
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
