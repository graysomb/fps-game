/* Offline Firefight fixtures. Angles and contacts advance on the gameplay fixed
 * tick; capture cameras never drive a mechanism. Released masonry becomes PBD debris. */
typedef struct {float angle,speed;bool broken;} LabTower;
typedef struct {Vector3 pos;float angle,speed;uint32_t victims;} LabLog;
static struct {
    LabTower towers[9];
    LabLog logs[3];
    float elapsed,swingAngle,swingSpeed;
    bool triggered,open;
    int owner;
    uint64_t release;
    uint32_t falling,landed;
} lab;
static const Color labSteel={72,109,113,255},labWood={135,79,48,255},labEnd={234,177,103,255};
static const Vector3 swingPivot={0,5.4f,0};
static Vector3 lab_swing_center(void){
    return (Vector3){0,swingPivot.y-4.6f*cosf(lab.swingAngle),4.6f*sinf(lab.swingAngle)};
}
static float lab_slope(float z){return fmaxf(0,(5-z)*.25f);}
static Vector3 lab_tower_pivot(const FoundryTrap *t){
    return (Vector3){(t->x0+t->x1+1)*.25f,t->y0*.5f,(t->z1+1)*.5f};
}
static Vector3 lab_rotate(Vector3 p,Vector3 pivot,float angle){
    Vector3 d=v_sub(p,pivot);
    return v_add(pivot,(Vector3){d.x,d.y*cosf(angle)-d.z*sinf(angle),d.y*sinf(angle)+d.z*cosf(angle)});
}
static void lab_smash(uint32_t *victims,int owner,Vector3 a,Vector3 b,Vector3 half){
    Vector3 movement=v_sub(b,a);float distance=v_length(movement);
    for(int i=0;i<activePlayers;i++){
        if((*victims&(1u<<i)) || players[i].respawn_timer>0 || players[i].invuln_timer>0 || !enemy_core_open(i))continue;
        BoundingBox box=combatant_bounds(i,PLAYER_SIZE*.5f);box.min=v_sub(box.min,half);box.max=v_add(box.max,half);
        if(!CheckCollisionBoxes(box,(BoundingBox){a,a})){
            if(distance<.0001f)continue;
            RayCollision hit=GetRayCollisionBox((Ray){a,v_mul(movement,1/distance)},box);
            if(!hit.hit || hit.distance>distance)continue;
        }
        *victims|=1u<<i;
        foundry_smash_actor(i,owner);
    }
}
static bool trap_lab_switch(int id,int owner){
    if(!lab.release || voxels[id].identity!=lab.release)return false;
    if(!lab.triggered){lab.triggered=true;lab.owner=owner;play_sfx_at(SFX_GLUE_BREAK,voxels[id].pos,owner,2,1);}
    return true;
}
static void lab_panel(bool build){
    for(int z=-4;z<4;z++)for(int x=-4;x<4;x++){
        int id=table_get_static_only(x,5,z);
        if(build && id<0)add_static_voxel_at_grid(x,5,z,(x==z || x==-z-1)?ORANGE:labSteel,0);
        if(!build && id>=0)remove_voxel_index(id);
    }
}
static void lab_pit_tick(float dt){
    bool occupied=false;
    for(int i=0;i<activePlayers;i++){
        if(players[i].respawn_timer>0)continue;
        BoundingBox bounds=combatant_bounds(i,PLAYER_SIZE*.5f);
        if(CheckCollisionBoxes(bounds,(BoundingBox){{-2,0,-2},{2,3.2f,2}}))occupied=true;
        if(fabsf(players[i].pos.x)<2 && fabsf(players[i].pos.z)<2){
            if(bounds.min.y>=2.9f && bounds.min.y<3.15f && !lab.triggered){lab.triggered=true;lab.owner=i;}
            if(lab.open && players[i].pos.y<3.15f)lab.falling|=1u<<i;
        }
        if((lab.falling&(1u<<i)) && !(lab.landed&(1u<<i)) && players[i].onGround && players[i].pos.y<1.5f){
            lab.landed|=1u<<i;apply_matter_damage(i,-1,25);
        }
    }
    if(!lab.triggered)return;
    lab.elapsed+=dt;
    if(!lab.open && lab.elapsed>=.5f){lab_panel(false);lab.open=true;play_sfx_at(SFX_GLUE_BREAK,(Vector3){0,3,0},-1,2,.7f);}
    if(lab.open && lab.elapsed>4 && !occupied){lab_panel(true);lab.triggered=lab.open=false;lab.elapsed=0;lab.falling=lab.landed=0;}
}
static void trap_lab_ramp(Player *p,float oldZ){
    bool pit=trapLabCase==LAB_PIT;
    float lo=pit?2:-11,hi=pit?8:5;
    if(fabsf(p->pos.x)>(pit?2:2.7f) || p->pos.z<lo || p->pos.z>hi)return;
    float height=pit?fminf(3,(p->pos.z-2)*.6f):lab_slope(p->pos.z);
    float previous=pit?fminf(3,(oldZ-2)*.6f):lab_slope(oldZ);
    float foot=p->pos.y-BASE_EYE_HEIGHT;
    if(p->vel.y<=0 && foot<=height+.08f && foot>=previous-.18f){
        p->pos.y=height+BASE_EYE_HEIGHT;p->vel.y=0;p->onGround=true;
    }
}
static void lab_tower_step(int n,float dt){
    FoundryTrap *t=&foundryTraps[n];LabTower *body=&lab.towers[n];
    if(body->broken)return;
    if(!t->active){
        if(!t->pending)return;
        t->active=true;t->pending=false;t->activation++;body->speed=.42f;
        for(int k=0;k<t->originalCount;k++){
            Voxel *v=&t->original[k];int id=table_get_static_only(v->gx,v->gy,v->gz);
            if(id>=0)remove_voxel_index(id);
        }
        play_sfx_at(SFX_GLUE_BREAK,lab_tower_pivot(t),t->owner,2,.8f);
    }
    Vector3 pivot=lab_tower_pivot(t);
    float old=body->angle;
    body->speed+=dt*2.9f*sinf(old+.08f);
    body->angle=fminf(PI*.5f,old+dt*body->speed);
    for(int k=0;k<t->originalCount-1;k++){
        Vector3 a=lab_rotate(t->original[k].pos,pivot,old),b=lab_rotate(t->original[k].pos,pivot,body->angle);
        if(v_length(v_sub(b,a))>dt*1.25f){
            lab_smash(&t->victims,t->owner,a,b,(Vector3){.25f,.25f,.25f});
            for(int j=0;j<foundryTrapCount;j++){
                FoundryTrap *other=&foundryTraps[j];
                if(j==n || other->active || other->pending)continue;
                BoundingBox wall={{other->x0*.5f,other->y0*.5f,other->z0*.5f},{(other->x1+1)*.5f,(other->y1+1)*.5f,(other->z1+1)*.5f}};
                if(CheckCollisionBoxSphere(wall,b,.25f)){other->pending=true;other->owner=t->owner;}
            }
        }
    }
    if(body->angle<PI*.5f)return;
    body->broken=true;
    for(int k=0;k<t->originalCount-1;k++){
        Voxel *v=&t->original[k];Vector3 at=lab_rotate(v->pos,pivot,body->angle);at.y-=pivot.y;
        int id=addVoxel(at.x,at.y,at.z,false,true,v->color,0);
        if(id>=0){voxels[id].owner=voxels[id].activator=t->owner;set_voxel_velocity(&voxels[id],(Vector3){(v->gx-t->x0-.5f)*.9f,.5f,.5f+(v->gy-t->y0)*.13f});}
    }
    play_sfx_at(SFX_IMPACT,pivot,t->owner,2,.7f);
}
static void trap_lab_tick(float dt){
    if(trapLabCase==LAB_PIT){lab_pit_tick(dt);return;}
    if(trapLabCase>=LAB_TOWER){for(int i=0;i<foundryTrapCount;i++)lab_tower_step(i,dt);return;}
    if(!lab.triggered)return;
    if(lab.elapsed==0){
        int id=enemy_voxel(lab.release);
        if(id>=0)remove_voxel_index(id);
    }
    lab.elapsed+=dt;
    int count=trapLabCase==LAB_SWING?1:3;
    for(int i=0;i<count;i++){
        LabLog *l=&lab.logs[i];Vector3 old=l->pos;
        if(trapLabCase==LAB_SWING){
            lab.swingSpeed+=(-GRAVITY/4.6f*sinf(lab.swingAngle)-.22f*lab.swingSpeed)*dt;
            lab.swingAngle+=lab.swingSpeed*dt;l->pos=lab_swing_center();l->angle=lab.swingAngle;
        }else{
            if(lab.elapsed<i*.3f)continue;
            l->speed=fmaxf(0,l->speed+dt*(old.z<5?1.65f:-1.5f));
            l->pos.z+=l->speed*dt;l->pos.y=lab_slope(l->pos.z)+.72f;
            l->angle+=Vector3Distance(l->pos,old)/.72f;
        }
        if(Vector3Distance(old,l->pos)>1.4f*dt)
            lab_smash(&l->victims,lab.owner,old,l->pos,(Vector3){1.9f,.65f,.65f});
    }
}
static void trap_lab_start(TrapLabCase which){
    enemy_lab_start(ENEMY_TYPE_STANDARD);clear_world_voxels();clear_pickups();
    memset(&lab,0,sizeof(lab));lab.owner=-1;trapLabCase=which;
    enemyLabGod=false;enemyLabAll=enemyLabPaused=false;worldPeriodic=false;
    Vector3 spawn={5,BASE_EYE_HEIGHT,-6},victim={0,BASE_EYE_HEIGHT,1};
    if(which==LAB_PIT){
        add_static_box_at_grid(-4,3,5,5,-20,-5,labSteel);
        add_static_box_at_grid(-4,3,5,5,16,20,labSteel);
        lab_panel(true);spawn=(Vector3){0,3+BASE_EYE_HEIGHT,-4};victim=(Vector3){1,3.26f,-9};
    }else if(which>=LAB_TOWER){
        int count=which==LAB_FOREST?9:1;
        for(int i=0;i<count;i++){
            int x=(i/3)*12,z=(i%3)*8;
            foundry_trap_define(FOUNDRY_TIP,x,x+1,0,11,z,z+1,x-1,0,z);
            FoundryTrap *t=&foundryTraps[i];
            for(int k=0;k<t->originalCount-1;k++){
                Voxel *v=&t->original[k];int id=table_get_static_only(v->gx,v->gy,v->gz);
                v->color=voxels[id].color=v->gy==11?labEnd:v->gy%4==0?labSteel:labWood;
            }
        }
        victim=(Vector3){.5f,BASE_EYE_HEIGHT,5.5f};
        spawn=(Vector3){-4,BASE_EYE_HEIGHT,-4};
    }else{
        lab.swingAngle=-1.08f;
        for(int i=0;i<3;i++)lab.logs[i].pos=which==LAB_SWING?lab_swing_center():(Vector3){0,lab_slope(-4-i*2)+.72f,-4-i*2};
        int id=add_static_voxel_at_grid(-6,2,-5,ORANGE,0);if(id>=0)lab.release=voxels[id].identity;
        if(which==LAB_LOGS){victim.z=5;add_static_box_at_grid(5,5,0,2,7,13,labSteel);}
    }
    players[0].pos=players[0].death_pos=spawn;players[0].invuln_timer=0;
    Vector3 aim=which==LAB_PIT?(Vector3){0,3,0}:which>=LAB_TOWER?(Vector3){-.25f,.25f,.25f}:(Vector3){-2.75f,1.25f,-2.25f};
    dir_to_yaw_pitch(v_norm(v_sub(aim,spawn)),&players[0].yaw,&players[0].pitch);
    enemy_lab_spawn(4,ENEMY_TYPE_STANDARD,victim);
    enemy_lab_spawn(5,ENEMY_TYPE_STANDARD,(Vector3){3.6f,BASE_EYE_HEIGHT,victim.z});
    init_static_hash();rebuild_voxel_hash();rebuild_all_voxel_surfaces();meshDirty=true;
}
static void lab_ramp_mesh(float width,float from,float y0,float to,float y1){
    Vector3 a={-width,y0,from},b={width,y0,from},c={width,y1,to},d={-width,y1,to};
    DrawTriangle3D(a,d,b,labSteel);DrawTriangle3D(b,d,c,labSteel);
    for(int i=0;i<=12;i++){
        float t=i/12.f;Vector3 l=Vector3Lerp(a,d,t),r=Vector3Lerp(b,c,t);l.y+=.015f;r.y+=.015f;
        DrawLine3D(l,r,(Color){148,166,152,255});
    }
    DrawCylinderEx(a,d,.045f,.045f,8,labEnd);DrawCylinderEx(b,c,.045f,.045f,8,labEnd);
}
static void lab_draw_log(LabLog l){
    Vector3 a=l.pos,b=l.pos;a.x-=1.9f;b.x+=1.9f;
    DrawCylinderEx(a,b,.72f,.72f,16,labWood);
    for(int side=-1;side<=1;side+=2){
        Vector3 end=side<0?a:b,cap=end;cap.x+=side*.025f;
        DrawCylinderEx(end,cap,.70f,.70f,16,labEnd);cap.x+=side*.012f;
        for(int k=0;k<3;k++){
            float angle=l.angle+k*2.09f;Vector3 tip=cap;tip.y+=cosf(angle)*.56f;tip.z+=sinf(angle)*.56f;
            DrawCylinderEx(cap,tip,.034f,.025f,6,labWood);
        }
    }
    for(int k=0;k<7;k++){float angle=l.angle+k*2*PI/7;Vector3 u=a,v=b;
        u.y+=.706f*cosf(angle);u.z+=.706f*sinf(angle);v.y=u.y;v.z=u.z;
        DrawCylinderEx(u,v,.027f,.027f,5,(Color){87,56,42,255});
    }
}
static void trap_lab_draw(void){
    if(trapLabCase==LAB_PIT){
        lab_ramp_mesh(2,2,0,7,3);lab_ramp_mesh(2,7,3,8,3);
        if(!lab.open){
            float wobble=lab.triggered?.035f*sinf(lab.elapsed*65):0;
            for(int i=-3;i<4;i++){Vector3 a={i*.5f,3.025f+wobble,i*.5f},b={i*.5f+.36f,3.025f+wobble,i*.5f+.5f};DrawCylinderEx(a,b,.038f,.03f,6,labEnd);}
        }
    }else if(trapLabCase==LAB_SWING){
        for(int s=-1;s<=1;s+=2){
            Vector3 top={s*2.5f,swingPivot.y,0},foot={top.x,0,0};DrawCylinderEx(foot,top,.18f,.18f,10,labSteel);
            Vector3 joint=top;joint.x=s*1.5f;Vector3 end=lab.logs[0].pos;end.x=joint.x;
            DrawCylinderEx(joint,end,.065f,.065f,8,labEnd);
        }
        DrawCylinderEx((Vector3){-2.7f,swingPivot.y,0},(Vector3){2.7f,swingPivot.y,0},.14f,.14f,10,labSteel);
        lab_draw_log(lab.logs[0]);
    }else if(trapLabCase==LAB_LOGS){
        lab_ramp_mesh(2.7f,-11,4,5,0);for(int i=0;i<3;i++)lab_draw_log(lab.logs[i]);
    }else for(int n=0;n<foundryTrapCount;n++){
        FoundryTrap *t=&foundryTraps[n];if(!t->active || lab.towers[n].broken)continue;
        Vector3 p=lab_tower_pivot(t);rlPushMatrix();rlTranslatef(p.x,p.y,p.z);rlRotatef(lab.towers[n].angle*RAD2DEG,1,0,0);
        for(int k=0;k<t->originalCount-1;k++){Vector3 at=v_sub(t->original[k].pos,p);DrawCube(at,.495f,.495f,.495f,t->original[k].color);}
        rlPopMatrix();
    }
}
