/* Enemy Lab presentation and attack state; never advanced inside a solver retry. */
static bool first_voxel_hit_filtered(Ray ray,float t_max,int ignore_id,bool static_only,int ignore_goliath_owner,VoxelHit *out_hit,const uint64_t *ignored,int ignore_count);
static bool go_part_alive(int part);
static Vector3 go_core_position(void);
static bool go_line_clear(Vector3 a,Vector3 b);
static Vector3 go_target_center(int target);
typedef struct {float value,lastHit,suppressed,flash;bool alive;int brokenFrame;} GoPartShield;
static GoPartShield goShields[GM_COUNT];
static int goDrainPart[MAX_PLAYERS];
static struct {float gunTimer,rocketTimer,rocketWindup,shot[2];int side,burst,bullets,rockets,throws;Vector3 rocketAim;} goArsenal;
static const char *go_part_name(int p){
    const char *names[]={"Root","Pelvis","Torso","Head","Left upper arm","Left forearm","Left hand","Right upper arm","Right forearm","Right hand","Left thigh","Left shin","Left foot","Right thigh","Right shin","Right foot","Chest","Reactor","Eye","Left shoulder gun","Right shoulder gun"};
    return p>=0&&p<GM_COUNT?names[p]:"Body";
}
static float go_shield_max(int p){
    if(p==GOLEM_HAND_L||p==GOLEM_HAND_R||p==GOLEM_FOOT_L||p==GOLEM_FOOT_R)return 50;
    if(p==GOLEM_HEAD||p==GOLEM_PELVIS)return 100;
    if(p==GOLEM_TORSO||p==GM_CHEST||p==GM_CORE)return 150;
    return 75;
}
static int go_shield_part(int v){
    if(!goCombatEnabled||goCombatDead||v<0||v>=voxel_count||!voxels[v].golemOwned||!voxels[v].persistentConstruction)return -1;
    for(int c=0;c<goCount;c++)if(goCells[c].id==voxels[v].identity)return goCells[c].part;
    return -1;
}
static bool go_shielded(int v){int p=go_shield_part(v);return p>=0&&goShields[p].value>0;}
static bool go_shield_hit(int v,float damage){
    int p=go_shield_part(v);if(p<0)return false;
    GoPartShield *s=&goShields[p];
    if(s->value<=0)return s->brokenFrame==recycleFrameCounter;
    s->lastHit=goFight.time;s->flash=.25f;s->value=fmaxf(0,s->value-fmaxf(0,damage));
    if(s->value<=0){s->brokenFrame=recycleFrameCounter;play_sfx_at(SFX_SHIELD_BREAK,voxels[v].pos,MAX_PLAYERS,2,.8f);}
    return true;
}
static void go_shields_tick(float dt){
    for(int p=1;p<GM_COUNT;p++){
        GoPartShield *s=&goShields[p];bool alive=go_part_alive(p);
        if(!alive||!s->alive){s->value=0;s->lastHit=goFight.time;s->brokenFrame=-1;}
        s->alive=alive;s->flash=fmaxf(0,s->flash-dt);s->suppressed=fmaxf(0,s->suppressed-dt);
        if(alive&&s->suppressed<=0&&goFight.time-s->lastHit>=8)s->value=fminf(go_shield_max(p),s->value+5*dt);
    }
}
static void go_shield_blast(Vector3 origin,float radius,float damage){
    if(!goCombatEnabled||goCombatDead)return;
    float amounts[GM_COUNT]={0};int indices[GM_COUNT];for(int p=0;p<GM_COUNT;p++)indices[p]=-1;
    for(int c=0;c<goCount;c++){
        int v=go_find(goCells[c].id,goCells[c].index);if(v<0)continue;int p=goCells[c].part;
        float d=Vector3Distance(origin,voxels[v].pos);if(d>=radius)continue;
        VoxelHit hit;Vector3 dir=v_norm(v_sub(voxels[v].pos,origin));
        if(first_voxel_hit_detailed((Ray){origin,dir},fmaxf(0,d-.2f),-1,false,&hit)&&go_shield_part(hit.id)!=p)continue;
        float amount=damage*(1-clampf((d-1)/(radius-1),0,1));
        if(amount>amounts[p]){amounts[p]=amount;indices[p]=v;}
    }
    for(int p=0;p<GM_COUNT;p++)if(indices[p]>=0)go_shield_hit(indices[p],amounts[p]);
}
static void go_drain_cancel(int i){if(i>=0&&i<MAX_PLAYERS)goDrainPart[i]=-1;}
static bool go_drain_start(int i,int v){
    if(i<0||i>=MAX_PLAYERS||!go_shielded(v))return false;
    goDrainPart[i]=go_shield_part(v);players[i].drainTarget=MAX_COMBATANTS+1;return true;
}
static Vector3 go_part_center(int part){
    Vector3 at={0};int count=0;for(int c=0;c<goCount;c++)if(goCells[c].part==part){
        int v=go_find(goCells[c].id,goCells[c].index);if(v>=0){at=v_add(at,voxels[v].pos);count++;}}
    return count?v_mul(at,1.f/count):go_core_position();
}
static void go_drain_tick(int i,float dt){
    if(!goCombatEnabled||i<0||i>=MAX_PLAYERS||goDrainPart[i]<0)return;
    int part=goDrainPart[i];Player *p=&players[i];Vector3 to=v_sub(go_part_center(part),p->pos);VoxelHit h;
    bool clear=first_voxel_hit_detailed((Ray){p->pos,v_norm(to)},v_length(to)+.5f,-1,false,&h)&&go_shield_part(h.id)==part;
    if(!goCombatEnabled||goCombatDead||p->respawn_timer>0||!p->drainTarget||!go_part_alive(part)||v_length(to)>TETHER_RANGE||v_dot(player_forward(p),v_norm(to))<cosf(65*DEG2RAD)||!clear){
        go_drain_cancel(i);p->drainTarget=0;return;}
    go_shield_hit(h.id,2*player_bullet_damage(i)*dt);
    // Keep the empty lock until release; a new press can physically grab the exposed part.
    if(goShields[part].value<=0)goShields[part].lastHit=goFight.time;
}
static Vector3 go_tether_anchor(void){
    Vector3 hand=go_part_center(GOLEM_HAND_R),dir=v_norm(v_sub(goFight.aim,hand));
    return v_add(hand,v_mul(dir,GOLD_TETHER_HOLD_FORWARD));
}
static Vector3 go_tether_aim(void){return goFight.aim;}
static void go_arsenal_drop(void){
    if(players[MAX_PLAYERS].tetherHolding){Vector3 zero={0};release_tether_velocity(MAX_PLAYERS,&zero);}
}
static void go_arsenal_reset(void){
    go_arsenal_drop();memset(goShields,0,sizeof(goShields));memset(&goArsenal,0,sizeof(goArsenal));
    for(int i=0;i<MAX_PLAYERS;i++){if(goDrainPart[i]>=0&&players[i].drainTarget==MAX_COMBATANTS+1)players[i].drainTarget=0;goDrainPart[i]=-1;}
}
static void go_arsenal_start(void){
    for(int p=1;p<GM_COUNT;p++)goShields[p]=(GoPartShield){.value=go_shield_max(p),.lastHit=-8,.alive=go_part_alive(p),.brokenFrame=-1};
    goArsenal.gunTimer=2;goArsenal.rocketTimer=15;goArsenal.shot[0]=goArsenal.shot[1]=-1000;
}
static bool go_world_ray(Ray ray,float distance,VoxelHit *hit){
    uint64_t ignored[goCount>0?goCount:1];int count=0;
    for(int c=0;c<goCount;c++)if(goCells[c].id)ignored[count++]=goCells[c].id;
    return first_voxel_hit_filtered(ray,distance,-1,false,-1,hit,ignored,count);
}
static bool go_arsenal_grab(void){
    if(players[MAX_PLAYERS].tetherHolding)return true;
    Vector3 origin=go_part_center(GOLEM_HAND_R);int best=-1;float score=-FLT_MAX;Vector3 normal={0,1,0};
    for(int v=0;v<voxel_count;v++){
        Voxel *x=&voxels[v];if(x->golemOwned||x->isBullet||x->type||voxel_is_fluid(x)||go_held(v))continue;
        float distance=Vector3Distance(origin,x->pos);if(distance>GOLD_TETHER_RANGE||distance<1)continue;
        VoxelHit hit;if(!go_world_ray((Ray){origin,v_norm(v_sub(x->pos,origin))},distance+.3f,&hit)||hit.id!=v)continue;
        float value=x->rest_volume-distance*.01f;
        if(value>score){score=value;best=v;normal=hit.normal;}
    }
    if(best<0)return false;
    players[MAX_PLAYERS].goldTetherCharged=true;acquire_tether_hit(MAX_PLAYERS,best,normal);
    return players[MAX_PLAYERS].tetherHolding;
}
static void go_arsenal_throw(void){
    if(players[MAX_PLAYERS].tetherHolding){release_tether(MAX_PLAYERS);goArsenal.throws++;goFight.cooldown[2]=10;}
}
static Vector3 go_gun_origin(int side){
    return v_add(go_part_center(GM_TURRET_L+side),v_mul(v_norm(v_sub(goFight.aim,go_part_center(GM_TURRET_L+side))),.9f));
}
static bool go_arsenal_tick(float dt,float distance,bool busy){
    goArsenal.gunTimer=fmaxf(0,goArsenal.gunTimer-dt);goArsenal.rocketTimer=fmaxf(0,goArsenal.rocketTimer-dt);
    if(busy||goFight.state==GC_RECOVER){goArsenal.rocketWindup=0;return false;}
    Vector3 head=go_part_center(GOLEM_HEAD),target=go_target_center(goFight.target);
    if(!go_part_alive(GOLEM_HEAD))goArsenal.rocketWindup=0;
    if(goArsenal.burst==0&&go_part_alive(GOLEM_HEAD)&&distance>=8&&distance<=35&&go_line_clear(head,target)&&(goArsenal.rocketTimer<=0||goArsenal.rocketWindup>0)){
        if(goArsenal.rocketWindup<.7f)goArsenal.rocketAim=target;
        goArsenal.rocketWindup+=dt;
        if(goArsenal.rocketWindup>=1){
            Vector3 dir=v_norm(v_sub(goArsenal.rocketAim,head));bool fired=false;
            for(int n=0;n<32;n++)if(!eqRockets[n].active){
                eqRockets[n]=(EquipmentRocket){true,MAX_PLAYERS,v_add(head,v_mul(dir,1)),v_mul(dir,22),5};fired=true;break;}
            goArsenal.rocketWindup=0;goArsenal.rocketTimer=fired?15:1;
            if(fired){goArsenal.rockets++;play_sfx_at(SFX_FIRE,head,MAX_PLAYERS,2,.55f);}
        }return true;
    }
    if(goArsenal.rocketWindup>0){goArsenal.rocketWindup=0;goArsenal.rocketTimer=1;}
    if(distance>30){goArsenal.burst=0;return false;}
    if(goArsenal.gunTimer>0)return goArsenal.burst>0;
    int side=goArsenal.side,arm=side?GOLEM_ARM_R:GOLEM_ARM_L;
    if(!go_part_alive(GM_TURRET_L+side)||!go_part_alive(arm)){goArsenal.side^=1;goArsenal.burst=0;goArsenal.gunTimer=.3f;return false;}
    Vector3 at=go_gun_origin(side);if(!go_line_clear(at,target)){goArsenal.burst=0;goArsenal.gunTimer=.25f;return false;}
    Vector3 dir=v_norm(v_sub(target,at)),right=v_norm(v_cross(dir,(Vector3){0,1,0})),up=v_cross(right,dir);
    float phase=(goArsenal.bullets*2.399963f),radius=tanf(4*DEG2RAD)*sqrtf(((goArsenal.bullets*37)%101)/100.f);
    dir=v_norm(v_add(dir,v_add(v_mul(right,cosf(phase)*radius),v_mul(up,sinf(phase)*radius))));
    if(spawn_default_projectile(MAX_PLAYERS,at,dir,false)){
        goArsenal.shot[side]=goFight.time;goArsenal.bullets++;goArsenal.burst++;play_sfx_at(SFX_FIRE,at,MAX_PLAYERS,1,1);
        if(goArsenal.burst==3){goArsenal.burst=0;goArsenal.side^=1;goArsenal.gunTimer=2;}else goArsenal.gunTimer=1.f/3;
    }else goArsenal.gunTimer=.2f;
    return true;
}
static void go_arsenal_draw(void){
    if(goCombatDead)return;
    Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};
    for(int side=0;side<2;side++)if(go_part_alive(GM_TURRET_L+side)&&go_part_alive(side?GOLEM_ARM_R:GOLEM_ARM_L)){
        Vector3 at=go_part_center(GM_TURRET_L+side),dir=v_norm(v_sub(goFight.aim,at));
        float yaw=atan2f(-dir.x,-dir.z),pitch=asinf(Clamp(dir.y,-1,1));
        Matrix m=MatrixMultiply(MatrixRotateX(pitch),MatrixRotateY(yaw));m=MatrixMultiply(m,MatrixTranslate(at.x,at.y,at.z));
        WeaponVisual visual={.shot_time=goArsenal.shot[side],.activity=1};
        weapon_draw((WeaponPose){.kind=WEAPON_STANDARD,.transform=m},&visual,eye,goFight.time,false,WEAPON_COLOR);
        if(goFight.time-goArsenal.shot[side]<.1f)DrawSphere(go_gun_origin(side),.15f,YELLOW);
    }
    if(go_part_alive(GOLEM_HEAD)){
        Vector3 head=go_part_center(GOLEM_HEAD),dir=v_norm(v_sub(goArsenal.rocketWindup>0?goArsenal.rocketAim:goFight.aim,head));
        Vector3 mouth=v_add(head,v_mul(dir,.5f));
        DrawCylinderEx(v_add(head,v_mul(dir,.3f)),mouth,.22f,.22f,8,DARKGRAY);
        DrawCylinderEx(mouth,v_add(mouth,v_mul(dir,.03f)),.16f,.16f,8,goArsenal.rocketWindup>0?ORANGE:BLACK);
    }
    rlDrawRenderBatchActive();rlDisableDepthMask();
    for(int p=1;p<GM_COUNT;p++)if(p!=GM_CORE&&(goShields[p].value>0||goShields[p].flash>0)&&go_part_alive(p)){
        BoundingBox box={{FLT_MAX,FLT_MAX,FLT_MAX},{-FLT_MAX,-FLT_MAX,-FLT_MAX}};
        for(int c=0;c<goCount;c++)if(goCells[c].part==p){int v=go_find(goCells[c].id,goCells[c].index);if(v>=0)
            for(int k=0;k<8;k++){box.min=Vector3Min(box.min,voxels[v].particles[k]->pos);box.max=Vector3Max(box.max,voxels[v].particles[k]->pos);}}
        Vector3 at=v_mul(v_add(box.min,box.max),.5f),size=v_add(v_sub(box.max,box.min),(Vector3){.05f,.05f,.05f});
        Color color=goShields[p].flash>0?WHITE:SKYBLUE;
        float recharge=goShields[p].value<go_shield_max(p)&&goFight.time-goShields[p].lastHit>=8&&goShields[p].suppressed<=0?.025f*(1+sinf(goFight.time*4)):0;
        if(goShields[p].value>0)DrawCube(at,size.x,size.y,size.z,Fade(color,.055f+recharge));
        float burst=goShields[p].value<=0?(1-goShields[p].flash/.25f)*.4f:0;
        DrawCubeWires(at,size.x+burst,size.y+burst,size.z+burst,Fade(color,goShields[p].value<=0?goShields[p].flash*3:.25f+.4f*goShields[p].value/go_shield_max(p)));
    }
    if(goArsenal.rocketWindup>0)DrawSphereWires(go_part_center(GOLEM_HEAD),.6f,6,10,ORANGE);
    if(players[MAX_PLAYERS].tetherHolding){
        Vector3 hand=go_part_center(GOLEM_HAND_R),end=tetherComByPlayer[MAX_PLAYERS];
        DrawCylinderEx(hand,end,.08f,.08f,6,(Color){255,185,35,100});
        DrawCylinderEx(hand,end,.05f,.05f,6,(Color){255,220,85,240});
        DrawSphere(end,.08f,GOLD);
    }

    rlDrawRenderBatchActive();rlEnableDepthMask();
}
