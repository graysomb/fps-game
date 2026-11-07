/* Versioned, pointer-free wave saves. Never serialize Player/Equipment runtime state. */
#define FP_SAVE_VERSION 1u
#define FP_SAVE_MAX_ITEMS 100000u
#define FP_SAVE_MAX_WAVE 100000
/* Fixed-width fields keep CPU/GPU builds compatible. Bump the version on schema changes. */
typedef struct {
    int32_t kills,deaths,debrisKills,grenades,special[2],selected,revives,gold;
    int32_t gun[2],gunSelected,upgrades,input,device;
} FpSavedPlayer;
typedef struct {
    uint32_t version,size,seed,nextId,unlocked,count;
    int32_t wave,lives,score,kills,repeatCursor,remaining[FF_TYPE_COUNT],pending[8],pendingCount,total;
    uint32_t upgrades[3];
    int32_t humans,world,foundry,endless,checkpoints,customMap,mapSlot;
    int32_t templeStage,megalithStage,megalithType,hyperStage,forerunnerStage,forerunnerType,sanctumSteps,sanctumTopology;
    FpSavedPlayer player[MAX_PLAYERS];
} FpSavedRun;
typedef struct {uint32_t id;int32_t type,owner,slot;float x,y,z;} FpSavedItem;
static int fpSaveWave=-1,fpSaveHumans;

static bool fp_save_path(char *out,size_t size,bool create) {
    const char *override=getenv("FPS_FIREFIGHT_SAVE");
    if(override) return *override && snprintf(out,size,"%s",override)>0 && strlen(override)<size;
#ifdef FIREFIGHT_TESTING
    return false; /* Tests must opt into an isolated path, never the player's save. */
#endif
    char directory[4096];const char *base;
#ifdef _WIN32
    base=getenv("APPDATA");
    if(!base || snprintf(directory,sizeof(directory),"%s/FPS Game",base)>=(int)sizeof(directory))return false;
#elif defined(__APPLE__)
    base=getenv("HOME");
    if(!base || snprintf(directory,sizeof(directory),"%s/Library/Application Support/FPS Game",base)>=(int)sizeof(directory))return false;
#else
    base=getenv("XDG_DATA_HOME");
    if(base && *base) {if(snprintf(directory,sizeof(directory),"%s/fps-game",base)>=(int)sizeof(directory))return false;}
    else {base=getenv("HOME");if(!base || snprintf(directory,sizeof(directory),"%s/.local/share/fps-game",base)>=(int)sizeof(directory))return false;}
#endif
    if(create && !DirectoryExists(directory) && MakeDirectory(directory)!=0)return false;
    return snprintf(out,size,"%s/firefight.sav",directory)<(int)size;
}
static uint32_t fp_save_hash(uint32_t hash,const void *data,size_t size) {
    const unsigned char *bytes=data;
    for(size_t n=0;n<size;n++)hash=(hash^bytes[n])*16777619u;
    return hash;
}
static bool fp_save_valid(const FpSavedRun *s) {
    if(s->version!=FP_SAVE_VERSION || s->size!=sizeof(*s) || s->count>FP_SAVE_MAX_ITEMS ||
       s->wave<1 || s->wave>FP_SAVE_MAX_WAVE || s->humans<1 || s->humans>MAX_PLAYERS ||
       s->world<0 || s->world>=WORLD_TYPE_COUNT || s->pendingCount<0 || s->pendingCount>8 ||
       s->lives<0 || s->lives>FP_SAVE_MAX_WAVE || s->score<0 || s->kills<0 ||
       s->repeatCursor<0 || s->repeatCursor>INT_MAX-100 || !s->nextId ||
       s->unlocked>>FP_RULE_COUNT || s->foundry<0 || s->foundry>1 || s->endless<0 || s->endless>1 ||
       s->checkpoints<0 || s->checkpoints>1 || s->customMap<0 || s->customMap>1 ||
       s->mapSlot<0 || s->mapSlot>2 || s->templeStage<0 || s->templeStage>4 ||
       s->megalithStage<0 || s->megalithStage>5 || s->megalithType<0 || s->megalithType>1 ||
       s->hyperStage<0 || s->hyperStage>4 || s->forerunnerStage<0 || s->forerunnerStage>4 ||
       s->forerunnerType<0 || s->forerunnerType>=FORERUNNER_ARCHETYPE_COUNT ||
       s->sanctumSteps<1 || s->sanctumSteps>100 || s->sanctumTopology< -1 || s->sanctumTopology>3)return false;
    int total=0;
    for(int t=0;t<FF_TYPE_COUNT;t++){if(s->remaining[t]<0 || s->remaining[t]>100000)return false;total+=s->remaining[t];}
    if(total!=s->total || total<=0)return false;
    for(int n=0;n<s->pendingCount;n++)if(s->pending[n]<0 || s->pending[n]>=FP_RULE_COUNT)return false;
    for(int k=0;k<3;k++)if(s->upgrades[k]>teamUpgradeCaps[k])return false;
    for(int i=0;i<s->humans;i++) {
        const FpSavedPlayer *p=&s->player[i];
        if(p->kills<0 || p->deaths<0 || p->debrisKills<0 || p->revives<0 ||
           p->grenades<0 || p->grenades>2 || p->special[0]<0 || p->special[0]>2 || p->special[1]<0 || p->special[1]>2 ||
           p->selected<0 || p->selected>2 || p->gold<0 || p->gold>1 || p->gunSelected<0 || p->gunSelected>1 ||
           p->upgrades<0 || (p->upgrades&~EQ_UPGRADE_MASK) ||
           (p->input!=INPUT_TYPE_KEYBOARD && p->input!=INPUT_TYPE_GAMEPAD) || p->device<0 || p->device>=4)return false;
        for(int g=0;g<2;g++)if(p->gun[g]!=EQ_EMPTY && (p->gun[g]<0 || p->gun[g]>=EQ_GUN_COUNT))return false;
    }
    return true;
}
/* Read into temporary storage; invalid or incomplete files cannot mutate a running game. */
static bool fp_save_read(FpSavedRun *s,FfItem **items) {
    char path[4096];*items=NULL;
    if(!fp_save_path(path,sizeof(path),false))return false;
    FILE *f=fopen(path,"rb");if(!f)return false;
    char magic[8];uint32_t expected;
    bool ok=fread(magic,1,8,f)==8 && !memcmp(magic,"FPSFFSV1",8) &&
        fread(&expected,sizeof(expected),1,f)==1 && fread(s,sizeof(*s),1,f)==1 && fp_save_valid(s);
    uint32_t hash=ok?fp_save_hash(2166136261u,s,sizeof(*s)):0;
    FfItem *copy=ok && s->count?calloc(s->count,sizeof(*copy)):NULL;
    if(ok && s->count && !copy)ok=false;
    for(uint32_t n=0;ok && n<s->count;n++) {
        FpSavedItem item;
        ok=fread(&item,sizeof(item),1,f)==1;
        if(!ok)break;
        hash=fp_save_hash(hash,&item,sizeof(item));
        ok=item.id && item.id<s->nextId && (n==0 || item.id>copy[n-1].id) && fp_rule(item.type)>=0 &&
            item.owner>=-2 && item.owner<s->humans && item.slot>=-1 && item.slot<=1 &&
            isfinite(item.x) && isfinite(item.y) && isfinite(item.z) &&
            fabsf(item.x)<10000 && fabsf(item.y)<10000 && fabsf(item.z)<10000;
        copy[n]=(FfItem){item.id,item.type,item.owner,item.slot,{item.x,item.y,item.z}};
    }
    ok=ok && hash==expected && fgetc(f)==EOF && !ferror(f);
    fclose(f);
    if(!ok){free(copy);return false;}
    *items=copy;return true;
}
static bool fp_save_available(void) {
    if(fpSaveWave<0){FpSavedRun s;FfItem *items;if(fp_save_read(&s,&items)){fpSaveWave=s.wave;fpSaveHumans=s.humans;free(items);}else fpSaveWave=0;}
    return fpSaveWave>0;
}
static void fp_save_write(void) {
    if(!fpSaveEnabled || !fp_active() || netTransport.role!=NET_ROLE_OFFLINE || !fpCheckpoint.valid)return;
    FfCheckpoint *c=&fpCheckpoint;
    FpSavedRun s={.version=FP_SAVE_VERSION,.size=sizeof(s),.seed=get_current_world_seed(),.nextId=c->nextId,
        .unlocked=c->unlocked,.count=(uint32_t)c->count,.wave=c->wave,.lives=c->lives,.score=c->score,.kills=c->kills,
        .repeatCursor=c->repeatCursor,.pendingCount=c->pendingCount,.total=c->total,.humans=firefightHumans,
        .world=currentWorldType,.foundry=firefightFoundry,.endless=firefightEndless,.checkpoints=firefightCheckpoints,
        .customMap=useCustomMap,.mapSlot=creativeMapSlot,.templeStage=templeTargetStage,.megalithStage=megalithTargetStage,
        .megalithType=megalithArchetype,.hyperStage=hyperTargetStage,.forerunnerStage=forerunnerTargetStage,
        .forerunnerType=forerunnerArchetype,.sanctumSteps=sanctumGrowthSteps,.sanctumTopology=sanctumTopologyMode};
    memcpy(s.remaining,c->remaining,sizeof(s.remaining));memcpy(s.pending,c->pending,sizeof(s.pending));memcpy(s.upgrades,c->upgrades,sizeof(s.upgrades));
    for(int i=0;i<s.humans;i++) {
        Player *p=&c->humans[i];Equipment *e=&c->equipment[i];
        s.player[i]=(FpSavedPlayer){.kills=p->kills,.deaths=p->deaths,.debrisKills=p->debrisKills,
            .grenades=c->grenades[i],.special={c->special[i][0],c->special[i][1]},.selected=c->selected[i],.revives=c->revives[i],
            .gold=p->goldTetherCharged,.gun={e->gun[0],e->gun[1]},.gunSelected=e->selected,.upgrades=e->upgrades,
            .input=playerInput[i],.device=seat_device(i)};
    }
    fpSaveError="Could not save Firefight. Check free space and save-folder access.";
    if(c->count>FP_SAVE_MAX_ITEMS || !fp_save_valid(&s))return;
    char path[4096],temp[4100];
    if(!fp_save_path(path,sizeof(path),true))return;
    snprintf(temp,sizeof(temp),"%s.tmp",path);
    FILE *f=fopen(temp,"wb");if(!f)return;
    uint32_t hash=fp_save_hash(2166136261u,&s,sizeof(s)),placeholder=0;
    bool ok=fwrite("FPSFFSV1",1,8,f)==8 && fwrite(&placeholder,4,1,f)==1 && fwrite(&s,sizeof(s),1,f)==1;
    for(size_t n=0;ok && n<c->count;n++) {
        FfItem *v=&c->items[n];FpSavedItem item={v->id,v->type,v->owner,v->slot,v->pos.x,v->pos.y,v->pos.z};
        hash=fp_save_hash(hash,&item,sizeof(item));ok=fwrite(&item,sizeof(item),1,f)==1;
    }
    ok=ok && fseek(f,8,SEEK_SET)==0 && fwrite(&hash,4,1,f)==1 && fflush(f)==0;
    if(fclose(f)!=0)ok=false;
#ifdef _WIN32
    if(ok)ok=MoveFileExA(temp,path,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0;
#else
    if(ok)ok=rename(temp,path)==0;
#endif
    if(!ok){remove(temp);TraceLog(LOG_WARNING,"%s",fpSaveError);return;}
    fpSaveWave=s.wave;fpSaveHumans=s.humans;fpSaveError=NULL;
}
static void fp_save_clear(void) {
    if(!fpSaveEnabled || netTransport.role!=NET_ROLE_OFFLINE)return;
    char path[4096];
    if(fp_save_path(path,sizeof(path),false) && FileExists(path) && remove(path)!=0) {
        fpSaveError="Could not remove the finished Firefight save.";return;
    }
    fpSaveEnabled=false;fpSaveWave=0;
}
static void fp_start_new_game(void) {
    ResetGame();
    fpSaveEnabled=fp_active() && netTransport.role==NET_ROLE_OFFLINE;
    fp_save_write();
}
static bool fp_save_continue(void) {
    if(netTransport.role!=NET_ROLE_OFFLINE)return false;
    FpSavedRun s;FfItem *items;
    if(!fp_save_read(&s,&items)) {
        fpSaveWave=0;fpSaveError="Save unavailable or incompatible. Choose New Game to start again.";return false;
    }
    if(!fp_reserve(s.count)){free(items);fpSaveError="Not enough memory to continue Firefight.";return false;}
    free(fpCheckpoint.items);
    fpCheckpoint=(FfCheckpoint){.valid=true,.wave=s.wave,.lives=s.lives,.score=s.score,.kills=s.kills,
        .repeatCursor=s.repeatCursor,.nextId=s.nextId,.unlocked=s.unlocked,.pendingCount=s.pendingCount,.total=s.total,.items=items,.count=s.count};
    memcpy(fpCheckpoint.remaining,s.remaining,sizeof(s.remaining));memcpy(fpCheckpoint.pending,s.pending,sizeof(s.pending));memcpy(fpCheckpoint.upgrades,s.upgrades,sizeof(s.upgrades));
    gameMode=GAME_MODE_FIREFIGHT;firefightPractice=enemyLab=equipmentLab=false;
    multiplayerActivePlayers=firefightHumans=s.humans;
    currentWorldType=s.world;set_current_world_seed(s.seed);
    firefightFoundry=s.foundry;firefightEndless=s.endless;firefightCheckpoints=s.checkpoints;
    useCustomMap=s.customMap;creativeMapSlot=s.mapSlot;
    templeTargetStage=s.templeStage;megalithTargetStage=s.megalithStage;megalithArchetype=s.megalithType;
    hyperTargetStage=s.hyperStage;forerunnerTargetStage=s.forerunnerStage;forerunnerArchetype=s.forerunnerType;
    sanctumGrowthSteps=s.sanctumSteps;sanctumTopologyMode=s.sanctumTopology;
    for(int i=0;i<s.humans;i++) {
        FpSavedPlayer *p=&s.player[i];Player *v=&fpCheckpoint.humans[i];
        v->kills=p->kills;v->deaths=p->deaths;v->debrisKills=p->debrisKills;
        v->kd_ratio=((float)p->kills+1)/((float)p->deaths+1);v->goldTetherCharged=p->gold;
        fpCheckpoint.equipment[i]=(Equipment){.initialized=true,.gun={p->gun[0],p->gun[1]},.selected=p->gunSelected,.upgrades=p->upgrades};
        fpCheckpoint.grenades[i]=p->grenades;fpCheckpoint.special[i][0]=p->special[0];fpCheckpoint.special[i][1]=p->special[1];
        fpCheckpoint.selected[i]=p->selected;fpCheckpoint.revives[i]=p->revives;
        multiplayerPlayerInput[i]=p->input;localSeats[i].device=p->device;
    }
    fpSaveEnabled=true;fp_restore();fpSaveError=NULL;fpSaveWave=s.wave;fpSaveHumans=s.humans;
    return true;
}
