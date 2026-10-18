#include "combat_feedback.inc"
/* Cosmetic effects and authoritative ordnance are sampled once, outside views. */
typedef struct {Sound sound; int source,id,priority; float until;} ClarityVoice;
static ClarityVoice clarityVoices[32];
static float claritySoundLast[MAX_COMBATANTS+1][SFX_COUNT];
static void clarity_audio_shutdown(void) {
    for(int n=0;n<32;n++)if(clarityVoices[n].sound.frameCount){UnloadSoundAlias(clarityVoices[n].sound);clarityVoices[n]=(ClarityVoice){0};}
}
static void play_sfx_at(SfxId id,Vector3 pos,int source,int priority,float pitch) {
    if(!sfxReady || !sfxEnabled)return;
    int listener=-1;float distance=30;
    for(int i=0;i<activePlayers;i++)if(!is_player_bot(i) && players[i].respawn_timer<=0 && (netTransport.role==NET_ROLE_OFFLINE || net_player_is_local(i))) {
        float d=v_length(world_delta(players[i].pos,pos));if(d<distance){distance=d;listener=i;}
    }
    if(listener<0)return;
    int src=source>=0 && source<MAX_COMBATANTS?source:MAX_COMBATANTS;
    if(clarityClock-claritySoundLast[src][id]<.12f)return;
    int voice=-1;for(int n=0;n<32;n++)if(!clarityVoices[n].sound.frameCount || !IsSoundPlaying(clarityVoices[n].sound)){voice=n;break;}
    if(voice<0)for(int n=0;n<32;n++)if(clarityVoices[n].priority<priority){voice=n;break;}
    if(voice<0)return;
    ClarityVoice *v=&clarityVoices[voice];if(v->sound.frameCount)UnloadSoundAlias(v->sound);
    v->sound=LoadSoundAlias(sfxSounds[id]);v->priority=priority;v->source=src;v->id=id;
    Vector3 delta=world_delta(pos,players[listener].pos),right=v_cross(player_forward(&players[listener]),(Vector3){0,1,0});
    float gain=clampf((30-distance)/27,0,1);VoxelHit hit;
    if(distance>.1f && first_voxel_hit_detailed((Ray){players[listener].pos,v_norm(delta)},distance-.2f,-1,true,&hit))gain*=.4f;
    SetSoundVolume(v->sound,gain*(id>=SFX_STEP_1?.25f:1));SetSoundPan(v->sound,.5f+.45f*v_dot(v_norm(delta),v_norm(right)));
    SetSoundPitch(v->sound,pitch);PlaySound(v->sound);claritySoundLast[src][id]=clarityClock;
}
static void clarity_reset(void) {
    for(int n=0;n<32;n++)if(clarityVoices[n].sound.frameCount)StopSound(clarityVoices[n].sound);
    pendingRemainsCount=0;memset(enemyRemains,0,sizeof(enemyRemains));memset(footTravel,0,sizeof(footTravel));memset(footVariant,0,sizeof(footVariant));
    memset(grenades,0,sizeof(grenades));memset(clarityEffects,0,sizeof(clarityEffects));
    memset(grenadeCooldown,0,sizeof(grenadeCooldown));memset(killConfirm,0,sizeof(killConfirm));
    memset(swarmHintActive,0,sizeof(swarmHintActive));
    swarmHintTimer=swarmPulseTimer=0;explosionDamageContext=false;
    for(int i=0;i<MAX_COMBATANTS;i++)grenadeInventory[i]=1;
}
static void clarity_effect(ClarityEffect e,bool broadcast) {
    int slot=0;for(int n=0;n<CLARITY_EFFECTS;n++){if(!clarityEffects[n].id){slot=n;break;}if(clarityEffects[n].age>clarityEffects[slot].age)slot=n;}
    clarityEffects[slot]=e;
    if(e.kind==1){play_sfx_at(SFX_DEATH,e.pos,-1,3,e.type==ENEMY_TYPE_SWARMER?1.6f:.8f);if(e.attacker>=0 && e.attacker<MAX_COMBATANTS){killConfirm[e.attacker]=.45f;if(net_player_is_local(e.attacker) && !is_player_bot(e.attacker))play_sfx(SFX_KILL);}}
    else play_sfx_at(SFX_EXPLOSION,e.pos,-1,3,1);
    if(broadcast)clarity_net_event(e);
}
static void clarity_death(int slot,int attacker) {
    swarmHintActive[slot]=false;
    if(netTransport.role==NET_ROLE_CLIENT)return;
    feedback_body(slot);
    clarity_effect((ClarityEffect){.id=clarityNextId++,.kind=1,.type=players[slot].enemyType,.attacker=attacker,.pos=players[slot].pos,.yaw=players[slot].yaw},true);
}
static bool clarity_shield_voxel(int v) {
    if(v<0 || v>=voxel_count)return false;
    int owner=voxels[v].owner;
    if(owner<MAX_PLAYERS || owner>=activePlayers || players[owner].enemyType!=ENEMY_TYPE_SHIELDBEARER)return false;
    EnemySpecialState *s=&enemySpecial[owner];
    for(int n=0;n<s->equipmentCount;n++)if(s->equipment[n].id==voxels[v].identity)return true;
    return false;
}
static bool clarity_shield_melee(int v,int attacker) {
    if(!clarity_shield_voxel(v))return false;
    int owner=enemy_equipment_owner(v),group=-1;EnemySpecialState *s=&enemySpecial[owner];
    for(int n=0;n<s->equipmentCount;n++)if(s->equipment[n].id==voxels[v].identity)group=s->equipment[n].group;
    Vector3 at=voxels[v].pos;
    for(int n=s->equipmentCount-1;n>=0;n--)if(s->equipment[n].group==group){int part=enemy_voxel(s->equipment[n].id);if(part>=0){enemy_detach(part);remove_voxel_index(part);}}
    play_sfx_at(SFX_SHIELD_BREAK,at,owner,3,1);players[owner].matter_flash_timer=.2f;
    if(attacker>=0 && attacker<activePlayers)ai_alert(owner,attacker,players[attacker].pos,AI_ALERT_HIT);return true;
}
static void clarity_actor_damage(int target,int owner,float damage,bool debris) {
    if(target<0 || target>=activePlayers || players[target].respawn_timer>0 || players[target].invuln_timer>0)return;
    bool prior=explosionDamageContext;explosionDamageContext=true;
    if(players[target].isExposed)kill_player(target,owner,true,debris);else apply_matter_damage(target,owner,damage);
    explosionDamageContext=prior;
}
static bool clarity_debris_hit(int v,int target) {
    Voxel *x=&voxels[v];
    if(!x->wasTethered && clarity_is_remains(x->identity))return true;
    if(x->explosionUntil<=clarityClock || x->sleeping)return false;
    if(v_length(x->vel)>2 && players[target].contactDamageTimer<=0) {
        clarity_actor_damage(target,x->owner,20,true);players[target].contactDamageTimer=.35f;
    }
    return true;
}
static void clarity_explode(Vector3 pos,int owner) {
    hologram_blast(pos,4,60);
    if(netTransport.role==NET_ROLE_CLIENT)return;
    float damage[MAX_COMBATANTS]={0};
    for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0) {
        Vector3 d=world_delta(players[i].pos,pos);float distance=v_length(d);VoxelHit hit;
        if(distance<4 && (distance<.1f || !first_voxel_hit_detailed((Ray){pos,v_norm(d)},fmaxf(0,distance-.25f),-1,true,&hit)))damage[i]=60*(1-distance/4);
    }
    enemy_blast_terrain((EnemyCharge){.pos=pos,.owner=owner});
    for(int i=0;i<activePlayers;i++)if(damage[i]>0) {
        clarity_actor_damage(i,owner,damage[i],false);
        Vector3 d=world_delta(players[i].pos,pos);d.y+=1;players[i].vel=v_add(players[i].vel,v_mul(v_norm(d),damage[i]/10));players[i].onGround=false;players[i].meleeKnockbackActive=true;
    }
    clarity_effect((ClarityEffect){.id=clarityNextId++,.kind=2,.pos=pos,.attacker=owner},true);
}
static void enemy_blast(EnemyCharge c){clarity_explode(c.pos,c.owner);}
static bool clarity_throw(int slot) {
    if(netTransport.role==NET_ROLE_CLIENT || slot<0 || slot>=activePlayers || creativeModeActive || (netTransport.role!=NET_ROLE_OFFLINE && netRequestedCreative) || gameState!=GAME_STATE_PLAYING || players[slot].respawn_timer>0 || ff_input_blocked(slot) || (enemyLab && enemyLabPaused) || grenadeInventory[slot]<=0 || grenadeCooldown[slot]>0)return false;
    int free=-1,owned=0;for(int n=0;n<CLARITY_GRENADES;n++){if(!grenades[n].id)free=n;else if(grenades[n].owner==slot)owned++;}
    if(free<0 || owned>=2)return false;
    Vector3 forward=player_forward(&players[slot]),pos=players[slot].pos;VoxelHit hit;
    float reach=first_voxel_hit_detailed((Ray){pos,forward},.6f,-1,false,&hit)?fmaxf(0,hit.t-.2f):.6f;
    pos=v_add(pos,v_mul(forward,reach));
    grenades[free]=(GrenadeState){.id=clarityNextId++,.owner=slot,.pos=pos,.vel=v_add(v_mul(forward,14),(Vector3){0,2,0}),.fuse=-1};
    grenadeInventory[slot]--;grenadeCooldown[slot]=.4f;play_sfx_at(SFX_MELEE,pos,slot,2,1.4f);return true;
}
static void clarity_grenade_step(GrenadeState *g,float dt) {
    const float radius=.15f;g->vel.y-=GRAVITY*dt;
    Vector3 movement=v_mul(g->vel,dt),normal={0};float length=v_length(movement),fraction=1;bool collision=false;
    if(length>1e-5f) {
        for(int c=0;c<9;c++) {
            Vector3 offset=c==8?(Vector3){0}:(Vector3){corner_signs[c][0]*radius,corner_signs[c][1]*radius,corner_signs[c][2]*radius};VoxelHit hit;
            if(first_voxel_hit_detailed((Ray){v_add(g->pos,offset),v_norm(movement)},length,-1,false,&hit) && hit.t/length<=fraction){fraction=fmaxf(0,hit.t/length-.001f);normal=hit.normal;collision=true;}
        }
        if(g->pos.y+movement.y<radius && movement.y<0){float f=clampf((radius-g->pos.y)/movement.y,0,1);if(f<=fraction){fraction=f;normal=(Vector3){0,1,0};collision=true;}}
        g->pos=v_add(g->pos,v_mul(movement,fraction));
        if(collision){float inward=v_dot(g->vel,normal);if(inward<0)g->vel=v_sub(g->vel,v_mul(normal,1.4f*inward));g->vel.x*=.8f;g->vel.z*=.8f;if(v_length(g->vel)<.6f)g->vel=(Vector3){0};g->pos=v_add(g->pos,v_mul(normal,.002f));}
    }
    g->pos=world_wrap(g->pos);
    /* Only the first contact starts the timer. Account for contact within this step. */
    if(collision && g->fuse==-1)
        g->fuse=fmaxf(1.5f-g->age,dt*fraction+.5f);
    g->age=fminf(1.5f,g->age+dt);
    if(g->fuse>=0)g->fuse=fmaxf(0,g->fuse-dt);
}
static bool clarity_evade(int slot,float dt) {
    if(slot<0 || !is_player_bot(slot) || players[slot].respawn_timer>0 || players[slot].enemyType==ENEMY_TYPE_FOUNDRY)return false;
    Vector3 threats[CLARITY_GRENADES+MAX_COMBATANTS];int count=0;
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id && v_length(world_delta(players[slot].pos,grenades[n].pos))<6)threats[count++]=grenades[n].pos;
    for(int n=0;n<MAX_COMBATANTS;n++)if(enemyCharges[n].active && v_length(world_delta(players[slot].pos,enemyCharges[n].pos))<6)threats[count++]=enemyCharges[n].pos;
    if(!count)return false;
    Vector3 best={0};float score=-1;
    for(int n=0;n<16;n++) {
        Vector3 direction={cosf(n*PI/8),0,sinf(n*PI/8)},at=v_add(players[slot].pos,v_mul(direction,3));VoxelHit hit;
        if((!worldPeriodic && (fabsf(at.x)>FLOOR_SIZE-1 || fabsf(at.z)>FLOOR_SIZE-1)) || first_voxel_hit_detailed((Ray){players[slot].pos,direction},2,-1,true,&hit))continue;
        float clearance=1000;for(int j=0;j<count;j++)clearance=fminf(clearance,v_length(world_delta(at,threats[j])));
        if(clearance>score){score=clearance;best=direction;}
    }
    if(enemy_is_special(slot)){enemy_drop_payload(slot);enemySpecial[slot].windup=0;enemySpecial[slot].state=0;enemySpecial[slot].reason="evading bomb";}
    bot_apply_move(&players[slot],&botStates[slot],best,MOVE_SPEED,dt);return true;
}
static void clarity_tick(float dt) {
    if(dt<=0)return;clarityClock+=dt;feedback_tick(dt);
    for(int n=0;n<CLARITY_EFFECTS;n++)if(clarityEffects[n].id){clarityEffects[n].age+=dt;if(clarityEffects[n].age>(clarityEffects[n].kind==2?.9f:.8f))clarityEffects[n].id=0;}
    for(int i=0;i<activePlayers;i++){grenadeCooldown[i]=fmaxf(0,grenadeCooldown[i]-dt);killConfirm[i]=fmaxf(0,killConfirm[i]-dt);}
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id) {
        int oldBeat=(int)(grenades[n].fuse*3);
        if(netTransport.role==NET_ROLE_CLIENT){if(grenades[n].fuse>=0)grenades[n].fuse=fmaxf(0,grenades[n].fuse-dt);grenades[n].pos=Vector3Lerp(grenades[n].pos,grenadeVisualTarget[n],fminf(1,dt*15));if((int)(grenades[n].fuse*3)!=oldBeat)play_sfx_at(SFX_SHIELD,grenades[n].pos,grenades[n].owner,2,1.4f);continue;}
        float remaining=fminf(dt,.25f);while(remaining>0 && grenades[n].id){float step=fminf(remaining,1.0f/120);clarity_grenade_step(&grenades[n],step);remaining-=step;if(grenades[n].fuse==0){GrenadeState g=grenades[n];grenades[n].id=0;clarity_explode(g.pos,g.owner);}}
        if(grenades[n].id && (int)(grenades[n].fuse*3)!=oldBeat)play_sfx_at(SFX_SHIELD,grenades[n].pos,grenades[n].owner,2,1.4f);
    }
    int alive=0,queued=0;
    for(int i=MAX_PLAYERS;i<activePlayers;i++)if(players[i].respawn_timer<=0) {
        alive++;
    }
    for(int t=0;t<FF_TYPE_COUNT;t++)queued+=ffRemaining[t];
    if(gameMode==GAME_MODE_FIREFIGHT && !enemyLab && firefightWaveState==FIREFIGHT_STATE_ACTIVE && queued==0 && alive>0 && alive<=3){swarmHintTimer+=dt;swarmPulseTimer+=dt;if(swarmHintTimer>=8 && (swarmPulseTimer>=4)){swarmPulseTimer=0;for(int i=MAX_PLAYERS;i<activePlayers;i++){swarmHintActive[i]=players[i].respawn_timer<=0;swarmHints[i]=players[i].pos;}}}
    else {swarmHintTimer=swarmPulseTimer=0;memset(swarmHintActive,0,sizeof(swarmHintActive));}
}
static void clarity_draw(void) {
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id && periodic_visible(grenades[n].pos,4)){GrenadeState *g=&grenades[n];DrawCube(g->pos,.22f,.3f,.22f,ORANGE);DrawCubeWires(g->pos,.26f,.34f,.26f,((int)(g->fuse*8)%2)?YELLOW:RED);DrawCircle3D(g->pos,4,(Vector3){1,0,0},90,Fade(ORANGE,.35f));}
    for(int n=0;n<CLARITY_EFFECTS;n++)if(clarityEffects[n].id){ClarityEffect *e=&clarityEffects[n];float t=e->age;Color color=e->kind==2?ORANGE:e->type==ENEMY_TYPE_SWARMER?RED:LIGHTGRAY;
        if(!periodic_visible(e->pos,4))continue;
        if(e->kind==2){feedback_blast(e);continue;}
        if(t<.12f)DrawCube(e->pos,1-t*3,1-t*3,1-t*3,YELLOW);
        for(int c=0;c<8;c++){Vector3 p=v_add(e->pos,(Vector3){corner_signs[c][0]*t*2,corner_signs[c][1]*t-t*t*3,corner_signs[c][2]*t*2});DrawCube(p,.3f*(1-t),.3f*(1-t),.3f*(1-t),Fade(color,1-t/.8f));}
    }
    /* Sort shield actors back to front and preserve the opaque depth buffer. */
    int order[MAX_COMBATANTS],count=0;float distances[MAX_COMBATANTS];
    Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};
    for(int i=MAX_PLAYERS;i<activePlayers;i++)if(players[i].enemyType==ENEMY_TYPE_SHIELDBEARER && players[i].respawn_timer<=0 && enemySpecial[i].equipmentCount) {
        if(!periodic_visible(players[i].pos,3))continue;
        float distance=v_length(world_delta(players[i].pos,eye));int at=count++;
        while(at>0 && distances[at-1]<distance){order[at]=order[at-1];distances[at]=distances[at-1];at--;}
        order[at]=i;distances[at]=distance;
    }
    if(!count)return;
    rlDrawRenderBatchActive();rlDisableDepthMask();
    for(int o=0;o<count;o++) {
        int i=order[o];EnemySpecialState *s=&enemySpecial[i];
        for(int n=0;n<s->equipmentCount;n++) {
            Vector3 at=enemy_local(i,s->equipment[n].offset);
            rlPushMatrix();rlTranslatef(at.x,at.y,at.z);rlRotatef(players[i].yaw,0,1,0);
            DrawCube((Vector3){0},.49f,.49f,.1f,Fade(SKYBLUE,.3f));
            DrawCubeWires((Vector3){0},.5f,.5f,.12f,s->shieldFlash[s->equipment[n].group]>0?YELLOW:BLUE);
            rlPopMatrix();
        }
    }
    rlDrawRenderBatchActive();rlEnableDepthMask();
}
static void clarity_hud(int slot,int w,int h) {
    if(slot<0 || slot>=MAX_PLAYERS || creativeModeActive)return;
    int local=netTransport.role==NET_ROLE_OFFLINE?slot:net_local_index_for_slot(slot);if(local<0)return;
    DrawText(TextFormat("Grenades %d/2 [%s]",grenadeInventory[slot],playerInput[local]==INPUT_TYPE_GAMEPAD?"D-pad Up":local==0?"C":"."),w-225,h-60,16,ORANGE);
    if(killConfirm[slot]>0)DrawText("ELIMINATED",w/2-55,h/2+30,16,YELLOW);
}
static bool clarity_supply_position(Vector3 anchor,int quarter,Vector3 *out) {
    for(int n=0;n<81;n++) {
        int ring=n?1+(n-1)/16:0;float a=n*PI/8;Vector3 offset=crosscut_rotate((Vector3){cosf(a)*ring*.4f,0,sinf(a)*ring*.4f},quarter);Vector3 at=v_add(anchor,offset);
        if(ring*.4f>2 || (!worldPeriodic && (fabsf(at.x)>FLOOR_SIZE-1 || fabsf(at.z)>FLOOR_SIZE-1)))continue;
        Ray ray={{at.x,12,at.z},{0,-1,0}};VoxelHit hit;float floor=0;
        if(first_voxel_hit_detailed(ray,12,-1,true,&hit))floor=12-hit.t;
        if(clarity_crosscut() && fabsf(anchor.x)+fabsf(anchor.z)>.1f && floor>.55f)continue;
        at.y=floor+.8f;
        bool clear=true;
        if(clarity_crosscut()) {
            clear=clarity_supply_clear(at);
            if(fabsf(anchor.x)+fabsf(anchor.z)>.1f)for(int q=1;q<4 && clear;q++) {
                Vector3 other=crosscut_rotate(at,q);clear=clarity_supply_clear(other);
                VoxelHit support;float ground=0;
                if(first_voxel_hit_detailed((Ray){{other.x,12,other.z},{0,-1,0}},12,-1,true,&support))ground=12-support.t;
                if(fabsf(ground+.8f-at.y)>.01f)clear=false;
            }
        }
        else for(int x=-1;x<=1;x++)for(int z=-1;z<=1;z++)for(int y=0;y<3;y++)if(table_get_static_only((int)floorf(at.x*2)+x,(int)floorf(at.y*2)+y,(int)floorf(at.z*2)+z)>=0)clear=false;
        if(clear){*out=at;return true;}
    }
    return false;
}
static void clarity_supply_layout(void) {
    if(useCustomMap || firefightPractice)return;
    bool foundry=gameMode==GAME_MODE_FIREFIGHT && firefightFoundry,crosscut=clarity_crosscut();
    if(!foundry && !crosscut)return;
    clear_pickups();int count=0;
    if(foundry) {
        Vector3 anchors[10]={{0,0,-10},{0,0,10},{-12,0,0},{12,0,0},{-14,0,-20},{-14,0,20},{14,0,-20},{14,0,20},{0,0,-28},{0,0,28}};
        PickupType types[10]={PICKUP_HEALTH,PICKUP_HEALTH,PICKUP_AMMO,PICKUP_AMMO,PICKUP_DYNAMIC_SHOT,PICKUP_DYNAMIC_SHOT,PICKUP_GOLD_TETHER,PICKUP_GOLD_TETHER,PICKUP_GRENADE,PICKUP_GRENADE};
        for(int n=0;n<10;n++){Vector3 at;if(clarity_supply_position(anchors[n],0,&at))pickups[count++]=(Pickup){.pos=at,.active=true,.type=types[n]};}
    }else {
        Vector3 positions[3];bool valid[3];
        for(int k=0;k<3;k++){Vector3 anchor=k==0?crosscut_supply_point(0):k==1?(Vector3){24,0,9}:(Vector3){9,0,6};valid[k]=clarity_supply_position(anchor,0,&positions[k]);if(!valid[k])TraceLog(LOG_ERROR,"Crosscut supply group %d has no accessible placement",k);}
        for(int q=0;q<4;q++)for(int k=0;k<3;k++)if(valid[k])pickups[count++]=(Pickup){.pos=crosscut_rotate(positions[k],q),.active=true,.type=k==0?PICKUP_AMMO:k==1?PICKUP_GRENADE:PICKUP_DYNAMIC_SHOT};
        Vector3 at;if(clarity_supply_position((Vector3){0},0,&at))pickups[count++]=(Pickup){.pos=at,.active=true,.type=PICKUP_GOLD_TETHER};
    }
}
static void clarity_supply_update(void) {
    if(useCustomMap || firefightPractice || (!firefightFoundry && currentWorldType!=WORLD_TYPE_CROSSCUT))return;
    for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active || pickups[n].respawnTimer>0) {
        Vector3 at=pickups[n].pos;VoxelHit hit;float floor=0;
        if(first_voxel_hit_detailed((Ray){at,{0,-1,0}},at.y,-1,true,&hit))floor=at.y-hit.t;
        if(at.y>floor+.85f)pickups[n].pos.y=floor+.8f;
    }
}
