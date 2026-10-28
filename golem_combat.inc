/* The combat controller runs once before/after a committed fixed tick. Pose
 * sampling is replayable; attacks and health never run inside a GPU attempt. */
typedef enum { GC_IDLE,GC_APPROACH,GC_WINDUP,GC_STRIKE,GC_RECOVER,GC_KNEEL,GC_DEAD } GoCombatState;
static const struct {float range,windup,strike,recovery,cooldown,damage;} goAttacks[3]={
    {4,1,.4f,1.2f,4,40},{5,.8f,.5f,1,3,25},{GOLD_TETHER_RANGE,1.5f,.35f,1.5f,10,50}};
static const float goShieldMax=150,goCoreMax=200,goShieldDelay=8,goShieldRate=5;
static struct {
    GoCombatState state;
    int target,attack,sequence,actions[3];
    float time,elapsed,shield,health,lastHit,suppressed,flash,broken,cooldown[3],routeTimer,travel;
    Vector3 position,waypoint,aim,lastLimb[2];
    float yaw,progress;
    bool hit[MAX_COMBATANTS+MAX_PLAYERS],launched,terrainHit;
    struct {uint64_t id[32];int count;float age;Vector3 previous[32];bool hit[MAX_COMBATANTS+MAX_PLAYERS];} fist;
} goFight;
#include "golem_arsenal.inc"
static void go_combat_reset(void){
    goRemoteActive=false;goRemoteTick=0;
    go_arsenal_reset();
    if(goCombatEnabled){players[MAX_PLAYERS].respawn_timer=9999;combatantLife[MAX_PLAYERS]=LIFE_UNUSED;}
    goInput.combatReach=false;goCombatEnabled=goCombatDead=false;memset(&goFight,0,sizeof(goFight));goFight.target=-1;
}
/* Load inferred from motor position error also includes commanded acceleration.
 * Do not fracture an unheld mount merely because its attack meets the floor.
 * A held segment still evaluates external overload; face cuts always apply. */
static bool go_combat_motor_grace(int part){
    if(!goCombatEnabled||goCombatDead)return false;
    if(goFight.state!=GC_WINDUP&&goFight.state!=GC_STRIKE&&goFight.state!=GC_RECOVER&&!(part==GM_CHEST&&goFight.state==GC_KNEEL))return false;
    for(int c=0;c<goCount;c++)if(goCells[c].part==part){int v=go_find(goCells[c].id,goCells[c].index);if(v>=0&&go_held(v))return false;}
    return true;
}
static bool go_part_alive(int part){return goMount[part]>=0&&go_find(goCells[goMount[part]].id,goCells[goMount[part]].index)>=0;}
static Vector3 go_core_position(void){return Vector3Transform((Vector3){0,.52f,-.06f},goPose.core);}
static void go_combat_sync(void){
    if(!goCombatEnabled)return;
    Player *p=&players[MAX_PLAYERS];p->pos=go_core_position();p->yaw=goFight.yaw;p->vel=(Vector3){0};
    p->matter=goFight.shield;p->matterMax=goShieldMax;p->isExposed=false;
    p->respawn_timer=goCombatDead?9999:0;combatantLife[MAX_PLAYERS]=goCombatDead?LIFE_UNUSED:LIFE_ALIVE;
}
static void go_combat_start(void){
    go_combat_reset();goCombatEnabled=true;
    goFight.position=goInput.position;goFight.yaw=goInput.yaw;goFight.target=-1;
    goFight.shield=goShieldMax;goFight.health=goCoreMax;goFight.lastHit=-goShieldDelay;
    players[MAX_PLAYERS].invuln_timer=0;players[MAX_PLAYERS].enemyType=ENEMY_TYPE_STANDARD;playerInput[MAX_PLAYERS]=INPUT_TYPE_BOT_MEDIUM;
    combatantTeam[MAX_PLAYERS]=1;goFight.lastLimb[0]=goPose.hand[0];goFight.lastLimb[1]=goPose.hand[1];
    go_arsenal_start();go_combat_sync();
}
static void go_combat_damage_from(float amount,bool shieldOnly,int attacker){
    if(!goCombatEnabled||goCombatDead||amount<=0||!isfinite(amount))return;
    if(shieldOnly&&goFight.shield<=0)return;
    float shield=goFight.shield;
    goFight.shield=fmaxf(0,shield-amount);
    if(!shieldOnly)goFight.health=fmaxf(0,goFight.health-fmaxf(0,amount-shield));
    goFight.lastHit=goFight.time;goFight.flash=.25f;
    if(shield>0&&goFight.shield<=0){goFight.broken=.8f;play_sfx_at(SFX_SHIELD_BREAK,go_core_position(),MAX_PLAYERS,2,.7f);}
    if(goFight.health<=0){
        go_arsenal_drop();
        goCombatDead=true;goFight.broken=1.2f;goFight.state=GC_DEAD;goRegen=false;goFight.fist.count=0;
        for(int p=1;p<GM_COUNT;p++)go_release_part(p);
        skeletonConstraintCount=0;play_sfx_at(SFX_EXPLOSION,go_core_position(),MAX_PLAYERS,2,.65f);
        for(int n=0;n<voxel_count;n++)if(voxels[n].golemOwned)voxels[n].golemReclaimStart=0;
        if(gameMode==GAME_MODE_FIREFIGHT && !enemyLab){
            /* Use the ordinary authoritative death path exactly once. */
            goCombatEnabled=false;
            kill_player(MAX_PLAYERS,attacker,true,false);
            goCombatEnabled=true;
        }
    }
    go_combat_sync();
}
static void go_combat_damage(float amount,bool shieldOnly){go_combat_damage_from(amount,shieldOnly,-1);}
/* Physical cuboids are not unit grid cells. Test their oriented faces, and
 * exclude the reactor's visual cells: the actor proxy owns reactor damage. */
static bool go_target_valid(int target){
    if(target_is_hologram(target))return holograms[target-MAX_COMBATANTS].life>0;
    return target>=0&&target<MAX_PLAYERS&&target<activePlayers&&!is_player_bot(target)&&players[target].respawn_timer<=0&&combatantLife[target]==LIFE_ALIVE;
}
static Vector3 go_target_center(int target){
    Player *p=ai_target_body(target);return v_add(p->pos,(Vector3){0,-.45f,0});
}
static bool go_line_clear(Vector3 a,Vector3 b){
    Vector3 d=v_sub(b,a);float distance=v_length(d);if(distance<.01f)return true;
    /* Ignore our own body for AI visibility, but never intact world cover. */
    Ray ray={a,v_mul(d,1/distance)};
    for(int n=0;n<voxel_count;n++){
        Voxel *v=&voxels[n];if(v->golemOwned||v->isBullet||(!v->vgs_active&&v->simulate))continue;
        float edge=v->rest_edge;Vector3 h={edge*.5f,edge*.5f,edge*.5f};
        RayCollision hit=GetRayCollisionBox(ray,(BoundingBox){v_sub(v->pos,h),v_add(v->pos,h)});
        if(hit.hit&&hit.distance<distance-.3f)return false;
    }return true;
}
/* An arena navigation grid, expanded by the body's shoulder clearance.
 * Regenerate on a timer because the arena's structures are destructible. */
static void go_route(Vector3 target){
    enum {N=96,CELLS=N*N};unsigned char blocked[CELLS]={0};int parent[CELLS],queue[CELLS];
    for(int n=0;n<CELLS;n++)parent[n]=-2;
    for(int n=0;n<voxel_count;n++){
        Voxel *v=&voxels[n];if(v->simulate||v->golemOwned||v->pos.y-v->rest_edge*.5f>5.3f||v->pos.y+v->rest_edge*.5f<.1f)continue;
        float r=v->rest_edge*.5f+2.1f;
        int x0=Clamp(floorf(v->pos.x-r+48),0,95),x1=Clamp(ceilf(v->pos.x+r+48),0,95);
        int z0=Clamp(floorf(v->pos.z-r+48),0,95),z1=Clamp(ceilf(v->pos.z+r+48),0,95);
        for(int z=z0;z<=z1;z++)for(int x=x0;x<=x1;x++)blocked[z*N+x]=1;
    }
    for(int z=0;z<N;z++)for(int x=0;x<N;x++)if(fabsf(x-48)>FLOOR_SIZE-3||fabsf(z-48)>FLOOR_SIZE-3)blocked[z*N+x]=1;
    int start=(int)Clamp(roundf(goFight.position.z+48),0,95)*N+(int)Clamp(roundf(goFight.position.x+48),0,95);
    int head=0,tail=0,best=start;float bestDist=FLT_MAX;parent[start]=-1;queue[tail++]=start;
    while(head<tail){int at=queue[head++],x=at%N,z=at/N;float dx=x-48-target.x,dz=z-48-target.z,d=dx*dx+dz*dz;
        if(d<bestDist){best=at;bestDist=d;}if(d<1)break;
        const int sx[4]={-1,1,0,0},sz[4]={0,0,-1,1};
        for(int k=0;k<4;k++){int nx=x+sx[k],nz=z+sz[k];if(nx<0||nx>=N||nz<0||nz>=N)continue;int next=nz*N+nx;
            if(!blocked[next]&&parent[next]==-2){parent[next]=at;queue[tail++]=next;}}
    }
    while(parent[best]>=0&&parent[best]!=start)best=parent[best];
    goFight.waypoint=(Vector3){best%N-48,.04f,best/N-48};goFight.routeTimer=.35f;
}
static bool go_move_clear(Vector3 position){
    if(fabsf(position.x)>FLOOR_SIZE-3||fabsf(position.z)>FLOOR_SIZE-3)return false;
    BoundingBox body={{position.x-1.9f,.08f,position.z-1.9f},{position.x+1.9f,5.3f,position.z+1.9f}};
    for(int n=0;n<voxel_count;n++){Voxel *v=&voxels[n];if(v->simulate||v->golemOwned)continue;
        Vector3 half={v->rest_edge*.5f,v->rest_edge*.5f,v->rest_edge*.5f};
        if(CheckCollisionBoxes(body,(BoundingBox){v_sub(v->pos,half),v_add(v->pos,half)}))return false;
    }return true;
}
static void go_begin_attack(int attack){
    goFight.attack=attack;goFight.state=GC_WINDUP;goFight.elapsed=0;goFight.progress=0;
    goFight.launched=goFight.terrainHit=false;memset(goFight.hit,0,sizeof(goFight.hit));
    goFight.sequence++;goFight.actions[attack]++;goFight.cooldown[attack]=goAttacks[attack].cooldown;
    play_sfx_at(SFX_TETHER,go_core_position(),MAX_PLAYERS,1,.6f+attack*.15f);
}
static void go_combat_prepare(float dt){
    if(!goCombatEnabled)return;
    for(int n=voxel_count-1;n>=0;n--){
        Voxel *v=&voxels[n];if(!v->golemOwned||v->persistentConstruction)continue;
        if(go_held(n))v->golemDebrisAge=0;else v->golemDebrisAge+=dt;
        if(v->golemDebrisAge>=45)remove_voxel_index(n);
    }
    if(goCombatDead){goFight.broken=fmaxf(0,goFight.broken-dt);return;}
    go_shields_tick(dt);
    goFight.time+=dt;goFight.elapsed+=dt;goFight.flash=fmaxf(0,goFight.flash-dt);goFight.broken=fmaxf(0,goFight.broken-dt);
    goFight.suppressed=fmaxf(0,goFight.suppressed-dt);goFight.routeTimer-=dt;
    for(int n=0;n<3;n++)goFight.cooldown[n]=fmaxf(0,goFight.cooldown[n]-dt);
    if(goFight.time-goFight.lastHit>=goShieldDelay&&goFight.suppressed<=0&&!ex_suppressed(MAX_PLAYERS))goFight.shield=fminf(goShieldMax,goFight.shield+goShieldRate*dt);
    bool kneel=!go_part_alive(GOLEM_SHIN_L)||!go_part_alive(GOLEM_SHIN_R)||goPose.kneel>.01f;
    if(kneel){go_arsenal_drop();goArsenal.rocketWindup=0;goFight.state=GC_KNEEL;goFight.elapsed=0;goFight.progress=0;go_combat_sync();return;}
    if(goFight.state==GC_KNEEL){go_arsenal_drop();goFight.state=GC_RECOVER;goFight.elapsed=0;}
    if(!go_target_valid(goFight.target)){
        goFight.target=-1;float nearest=FLT_MAX;
        for(int k=0;k<MAX_PLAYERS*2;k++){int t=k<MAX_PLAYERS?k:MAX_COMBATANTS+k-MAX_PLAYERS;if(!go_target_valid(t))continue;
            float distance=Vector3DistanceSqr(goFight.position,ai_target_body(t)->pos);if(distance<nearest){goFight.target=t;nearest=distance;}}
    }
    if(goFight.target<0){go_arsenal_drop();goArsenal.rocketWindup=0;goFight.state=GC_IDLE;goFight.progress=0;go_combat_sync();return;}
    Vector3 target=go_target_center(goFight.target),delta=v_sub(target,goFight.position);delta.y=0;float distance=v_length(delta);
    bool attackState=goFight.state==GC_WINDUP||goFight.state==GC_STRIKE;
    bool locked=goFight.state==GC_STRIKE||(goFight.state==GC_WINDUP&&goFight.elapsed>=goAttacks[goFight.attack].windup-.25f);
    if(!locked){goFight.aim=target;float yaw=atan2f(-delta.x,-delta.z)*RAD2DEG;
        float change=fmodf(yaw-goFight.yaw+540,360)-180;goFight.yaw+=Clamp(change,-90*dt,90*dt);}
    bool ranged=go_arsenal_tick(dt,distance,attackState);
    if(ranged){go_combat_sync();return;}
    if(attackState){
        int arm=goFight.attack==1?GOLEM_ARM_L:GOLEM_ARM_R;
        if(!go_part_alive(arm)||(!goFight.launched&&!go_part_alive(arm+2))||(goFight.attack==2&&!goFight.launched&&!players[MAX_PLAYERS].tetherHolding)){goFight.state=GC_RECOVER;goFight.elapsed=0;}
        else if(goFight.state==GC_WINDUP){
            if(goFight.attack<2&&!locked&&distance>2.2f){
                float step=fminf(2*dt,distance-2.2f);Vector3 next=v_add(goFight.position,v_mul(v_norm(delta),step));
                if(go_move_clear(next)){goFight.position=next;goFight.travel+=step;}
            }
            goFight.progress=.4f*Clamp(goFight.elapsed/goAttacks[goFight.attack].windup,0,1);
            if(goFight.elapsed>=goAttacks[goFight.attack].windup){goFight.state=GC_STRIKE;goFight.elapsed=0;play_sfx_at(SFX_IMPACT,go_core_position(),MAX_PLAYERS,2,.6f);}
        }else{
            goFight.progress=Lerp(.4f,.7f,Clamp(goFight.elapsed/goAttacks[goFight.attack].strike,0,1));
            if(goFight.attack==2&&!goFight.launched&&goFight.progress>=.6f){goFight.launched=true;go_arsenal_throw();}
            if(goFight.elapsed>=goAttacks[goFight.attack].strike){goFight.state=GC_RECOVER;goFight.elapsed=0;}
        }
    }else if(goFight.state==GC_RECOVER){
        goFight.progress=Lerp(.7f,1,Clamp(goFight.elapsed/goAttacks[goFight.attack].recovery,0,1));
        if(goFight.elapsed>=goAttacks[goFight.attack].recovery){goFight.state=GC_IDLE;goFight.progress=0;}
    }else{
        bool clear=go_line_clear(go_core_position(),target);int choice=-1;
        for(int a=0;a<3;a++){int arm=a==1?GOLEM_ARM_L:GOLEM_ARM_R;
            if(goFight.cooldown[a]<=0&&go_part_alive(arm)&&go_part_alive(arm+2)&&clear&&distance<=goAttacks[a].range&&(a!=2||distance>=5)) {choice=a;break;}}
        float angle=fabsf(fmodf(atan2f(-delta.x,-delta.z)*RAD2DEG-goFight.yaw+540,360)-180);
        if(choice==2 && angle<20 && !go_arsenal_grab()){goFight.cooldown[2]=2;choice=-1;}
        if(choice>=0&&angle<20)go_begin_attack(choice);
        else if(choice==2)go_arsenal_drop();
        else if(distance>3.2f||!clear){
            goFight.state=GC_APPROACH;
            if(goFight.routeTimer<=0||Vector3Distance(goFight.position,goFight.waypoint)<.15f)go_route(target);
            Vector3 move=v_sub(goFight.waypoint,goFight.position);move.y=0;float step=fminf(2*dt,v_length(move));Vector3 next=v_add(goFight.position,v_mul(v_norm(move),step));
            if(go_move_clear(next)){goFight.position=next;goFight.travel+=step;}else goFight.routeTimer=0;
        }else goFight.state=GC_IDLE;
    }
    go_combat_sync();
}
static void go_combat_pose(void){
    goInput.combatReach=true;goInput.attackTarget=goFight.aim;
    goInput.position=goFight.position;goInput.yaw=goFight.yaw;goInput.travel=goFight.travel;goInput.side=goFight.attack==1?0:1;
    goInput.action=goFight.state==GC_APPROACH?GOLEM_WALK:(goFight.state==GC_WINDUP||goFight.state==GC_STRIKE||goFight.state==GC_RECOVER)?(GolemAction)(GOLEM_SLAM+goFight.attack):GOLEM_IDLE;
    goInput.actionProgress=goFight.progress;
}
static Vector3 go_ballistic_velocity(Vector3 delta,float flight){
    /* Match semi-implicit integration at the solver's actual substep damping. */
    const float dt=PBD_MAX_STEP_DT/PBD_SUBSTEPS;int steps=(int)ceilf(flight/dt);
    float displacement=0,factor=1,vy=0,gravityDisplacement=0;
    for(int s=0;s<steps;s++){factor*=VELOCITY_DAMPING;displacement+=factor*dt;vy=vy*VELOCITY_DAMPING-GRAVITY*dt;gravityDisplacement+=vy*dt;}
    Vector3 velocity=v_mul(delta,1/fmaxf(.001f,displacement));velocity.y-=gravityDisplacement/fmaxf(.001f,displacement);return velocity;
}
static void go_combat_launch(Vector3 origin,Vector3 *target,Vector3 *velocity){
    if(goCombatEnabled){
        *target=goFight.aim;memset(&goFight.fist,0,sizeof(goFight.fist));
        for(int c=0;c<goCount;c++)if(goCells[c].part==GOLEM_HAND_R&&goCells[c].index>=0&&goFight.fist.count<32){
            int n=goFight.fist.count++;goFight.fist.id[n]=goCells[c].id;goFight.fist.previous[n]=voxels[goCells[c].index].pos;}
    }
    Vector3 delta=v_sub(*target,origin);float flight=fmaxf(.2f,sqrtf(delta.x*delta.x+delta.z*delta.z)/22);
    *velocity=go_ballistic_velocity(delta,flight);
}
static void go_hit_target(int target,float damage,Vector3 origin){
    if(target_is_hologram(target)){hologram_damage(target,damage);return;}
    Player *p=&players[target];if(p->invuln_timer>0||(enemyLab&&enemyLabGod))return;
    bool known=damageOriginKnown;Vector3 old=damageOrigin;damageOriginKnown=true;damageOrigin=origin;
    if(p->isExposed)kill_player(target,MAX_PLAYERS,true,true);else apply_matter_damage(target,MAX_PLAYERS,damage);
    damageOriginKnown=known;damageOrigin=old;
    Vector3 push=v_norm(v_sub(p->pos,origin));push.y=fmaxf(.25f,push.y);p->vel=v_add(p->vel,v_mul(v_norm(push),6));p->meleeKnockbackActive=true;
}
static void go_sweep_hits(Vector3 a,Vector3 b,float radius,float damage,bool *hit){
    for(int k=0;k<MAX_PLAYERS*2;k++){
        int t=k<MAX_PLAYERS?k:MAX_COMBATANTS+k-MAX_PLAYERS;if(hit[t]||!go_target_valid(t))continue;
        BoundingBox body=target_is_hologram(t)?combatant_bounds_at(t-MAX_COMBATANTS,ai_target_body(t)->pos,PLAYER_SIZE*.5f):combatant_bounds(t,PLAYER_SIZE*.5f);
        body=expand_bounds(body,radius);
        if(segment_intersects_aabb(a,b,body.min,body.max)){hit[t]=true;go_hit_target(t,damage,a);}
    }
}
static Vector3 go_limb_center(int side){
    int part=side?GOLEM_HAND_R:GOLEM_HAND_L;Vector3 sum={0};int count=0;
    for(int c=0;c<goCount;c++)if(goCells[c].part==part){int n=go_find(goCells[c].id,goCells[c].index);if(n>=0){sum=v_add(sum,voxels[n].pos);count++;}}
    return count?v_mul(sum,1.f/count):goPose.hand[side];
}
static void go_combat_commit(float dt){
    if(!goCombatEnabled||goCombatDead)return;
    Vector3 limb[2]={go_limb_center(0),go_limb_center(1)};
    if(goFight.state==GC_STRIKE&&goFight.attack<2){int side=goFight.attack==1?0:1;
        go_sweep_hits(goFight.lastLimb[side],limb[side],.65f,goAttacks[goFight.attack].damage,goFight.hit);
        if(!goFight.terrainHit&&goFight.progress>=.6f){
            goFight.terrainHit=true;Vector3 at=limb[side];int x=floorf(at.x/VOXEL_SIZE),y=floorf(at.y/VOXEL_SIZE),z=floorf(at.z/VOXEL_SIZE);
            activate_static_voxels_near_region(x-1,x+1,y-1,y+1,z-1,z+1,MAX_PLAYERS);
        }
    }
    for(int n=0;n<2;n++)goFight.lastLimb[n]=limb[n];
    goFight.fist.age+=dt;
    for(int n=0;n<goFight.fist.count;n++){
        int v=go_find(goFight.fist.id[n],-1);if(v<0||go_held(v)||voxels[v].wasTethered)continue;
        if(goFight.fist.age<3)go_sweep_hits(goFight.fist.previous[n],voxels[v].pos,.3f,goAttacks[2].damage,goFight.fist.hit);
        goFight.fist.previous[n]=voxels[v].pos;
    }
    go_combat_sync();
}
static const char *go_combat_status(void){
    const char *state[]={"seeking","approaching","WINDUP","ATTACK","recovering","KNEEL / CORE OPEN","DEFEATED"};
    return state[goFight.state];
}
static void go_combat_draw(void){
    if(!goCombatEnabled)return;go_arsenal_draw();Vector3 core=go_core_position();
    if(goCombatDead){if(goFight.broken>0)DrawSphereWires(core,.4f+(1.2f-goFight.broken)*1.5f,6,8,Fade(ORANGE,goFight.broken*.4f));return;}
    float radius=.52f+.015f*sinf(goFight.time*4);Color color=goFight.flash>0?WHITE:goFight.time-goFight.lastHit>=goShieldDelay&&goFight.suppressed<=0?SKYBLUE:BLUE;
    rlDrawRenderBatchActive();rlDisableDepthMask();
    if(goFight.shield>0){DrawSphereEx(core,radius,8,12,Fade(color,.08f+.12f*goFight.shield/goShieldMax));DrawSphereWires(core,radius,6,12,Fade(color,.55f));}
    if(goFight.shield>0&&goFight.shield<goShieldMax&&goFight.time-goFight.lastHit>=goShieldDelay&&goFight.suppressed<=0){
        Vector3 ring=core;ring.y+=sinf(goFight.time*4)*.38f;
        DrawCircle3D(ring,.4f,(Vector3){1,0,0},90,SKYBLUE);
    }
    if(goFight.flash>0)DrawSphereWires(core,radius+(1-goFight.flash/.25f)*.25f,8,12,Fade(WHITE,goFight.flash*3));
    if(goFight.broken>0)DrawSphereWires(core,radius+(1-goFight.broken/.8f)*.9f,8,12,Fade(SKYBLUE,goFight.broken));
    if(goFight.state==GC_WINDUP){
        Vector3 mark=goFight.attack==2?goFight.aim:goPose.hand[goFight.attack==1?0:1];mark.y=.05f;
        DrawCircle3D(mark,goFight.attack==2?.6f:1.2f,(Vector3){1,0,0},90,ORANGE);
        DrawLine3D(goPose.hand[goFight.attack==1?0:1],goFight.aim,Fade(ORANGE,.55f));
    }
    rlDrawRenderBatchActive();rlEnableDepthMask();
}
static void go_remote_draw(void){
    if(!goRemoteActive || goRemoteDead)return;
    rlDrawRenderBatchActive();rlDisableDepthMask();
    for(int p=1;p<GM_COUNT;p++)if(goRemoteParts[p].alive && goRemoteParts[p].shield>0){
        BoundingBox b=goRemoteParts[p].box;Vector3 at=v_mul(v_add(b.min,b.max),.5f),size=v_sub(b.max,b.min);
        DrawCubeWires(at,size.x+.05f,size.y+.05f,size.z+.05f,Fade(goRemoteParts[p].flash>0?WHITE:SKYBLUE,.6f));
    }
    if(goRemoteShield>0)DrawSphereWires(goRemoteCore,.52f,6,12,SKYBLUE);
    rlDrawRenderBatchActive();rlEnableDepthMask();
    for(int side=0;side<2;side++)if(goRemoteParts[GM_TURRET_L+side].alive){
        BoundingBox b=goRemoteParts[GM_TURRET_L+side].box;Vector3 at=v_mul(v_add(b.min,b.max),.5f),dir=v_norm(v_sub(goRemoteAim,at));
        Matrix m=MatrixMultiply(MatrixRotateX(asinf(Clamp(dir.y,-1,1))),MatrixRotateY(atan2f(-dir.x,-dir.z)));m=MatrixMultiply(m,MatrixTranslate(at.x,at.y,at.z));
        WeaponVisual visual={.activity=1};Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};
        weapon_draw((WeaponPose){.kind=WEAPON_STANDARD,.transform=m},&visual,eye,(float)GetTime(),false,WEAPON_COLOR);
    }
}

