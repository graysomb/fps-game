/* Gameplay feedback uses stable voxel IDs; no per-view simulation. */
typedef struct {uint64_t id; float born;} EnemyRemains;
static EnemyRemains enemyRemains[32];
typedef struct {Vector3 pos,vel;float edge;Color color;} PendingRemains;
static PendingRemains pendingRemains[32];
static int pendingRemainsCount;
static float footTravel[MAX_COMBATANTS];
static unsigned footVariant[MAX_COMBATANTS];
static Shader blastShader;
static bool blastShaderAttempted;
static int blastAgeLocation,blastCenterLocation,blastEyeLocation,blastSmokeLocation;
static bool clarity_is_remains(uint64_t id) {
    if(!id)return false;
    for(int n=0;n<32;n++)if(enemyRemains[n].id==id)return true;
    return false;
}
static bool remains_held(int v) {
    if(tetherTag[v] || tetherHoldClusterTag[v])return true;
    for(int i=0;i<activePlayers;i++)if(players[i].tetherHolding && players[i].tetherVoxelIdentity==voxels[v].identity)return true;
    return false;
}
static void feedback_body(int slot) {
    if(slot<0 || slot>=activePlayers || !is_player_bot(slot) || pendingRemainsCount>=32)return;
    Player *p=&players[slot];float half=PLAYER_SIZE*(p->enemyType==ENEMY_TYPE_SWARMER?.25f:.5f);
    BoundingBox body=combatant_bounds(slot,half);
    float edge=clampf(body.max.x-body.min.x,.5f,1.5f);
    Vector3 at=v_mul(v_add(body.min,body.max),.5f);at.y=fmaxf(edge*.5f+.05f,at.y);
    pendingRemains[pendingRemainsCount++]=(PendingRemains){at,v_mul(p->vel,.25f),edge,p->enemyType==ENEMY_TYPE_SWARMER?(Color){235,70,45,255}:p->enemyType==ENEMY_TYPE_GOLIATH?(Color){35,38,48,255}:enemy_is_special(slot)?enemy_color(p->enemyType):player_palette_color(slot)};
}
static void feedback_spawn_body(PendingRemains body) {
    int free=-1;float oldest=1e30f;
    for(int n=0;n<32;n++) {
        int v=enemy_voxel(enemyRemains[n].id);
        if(v<0){free=n;break;}
        if(!remains_held(v) && enemyRemains[n].born<oldest){oldest=enemyRemains[n].born;free=n;}
    }
    if(free<0)return;
    int old=enemy_voxel(enemyRemains[free].id);if(old>=0)remove_voxel_index(old);
    float edge=body.edge;Vector3 at=body.pos;
    int v=addVoxel(at.x,at.y,at.z,false,true,body.color,0);
    if(v<0){enemyRemains[free]=(EnemyRemains){0};return;}
    Voxel *x=&voxels[v];x->owner=ENEMY_CONSTRUCTION_OWNER;x->activator=-1;x->glueEligible=false;
    x->rest_edge=edge;x->rest_volume=edge*edge*edge;x->particle_radius=edge*.5f;
    x->vel=body.vel;x->activationBelief=1;
    for(int c=0;c<8;c++) {
        Particle *part=x->particles[c];
        part->pos=v_add(at,(Vector3){corner_signs[c][0]*edge*.5f,corner_signs[c][1]*edge*.5f,corner_signs[c][2]*edge*.5f});
        part->prev_pos=part->predicted_pos=part->pos;part->vel=x->vel;part->radius=voxel_particle_radius(x);
    }
    enemyRemains[free]=(EnemyRemains){x->identity,clarityClock};
}
static void feedback_tick(float dt) {
    if(netTransport.role!=NET_ROLE_CLIENT) {
        for(int n=0;n<pendingRemainsCount;n++)feedback_spawn_body(pendingRemains[n]);
        pendingRemainsCount=0;
        for(int n=0;n<32;n++) {
            int v=enemy_voxel(enemyRemains[n].id);
            if(v<0){enemyRemains[n].id=0;continue;}
            voxels[v].lifeFrames=0; /* Remains use seconds, not generic frame recycling. */
            if(clarityClock-enemyRemains[n].born>=45 && !remains_held(v)) {
                remove_voxel_index(v);enemyRemains[n].id=0;
            }
        }
        for(int i=MAX_PLAYERS;i<activePlayers;i++) {
            EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];
            s->healingPulse=fmaxf(0,s->healingPulse-dt);
            for(int g=0;g<4;g++)s->shieldFlash[g]=fmaxf(0,s->shieldFlash[g]-dt);
            if(p->respawn_timer>0){s->exposedAge=0;continue;}
            s->exposedAge=p->isExposed?s->exposedAge+dt:0;
            if(p->isExposed && s->exposedAge<2)continue;
            if(p->matter>=p->matterMax)continue;
            for(int a=MAX_PLAYERS;a<activePlayers;a++)if(enemy_braces(a,i)) {
                if(ex_suppressed(i))continue;
                add_player_matter(p,8*dt);
                if(p->enemyType==ENEMY_TYPE_FOUNDRY)s->coreHealth=p->matter;
                s->healingPulse=.2f;break;
            }
        }
    }
    for(int i=MAX_PLAYERS;i<activePlayers;i++) {
        Player *p=&players[i];float speed=sqrtf(p->vel.x*p->vel.x+p->vel.z*p->vel.z);
        if(p->respawn_timer>0 || !p->onGround || speed<.3f){footTravel[i]=0;continue;}
        float stride=p->enemyType==ENEMY_TYPE_SWARMER?.65f:p->enemyType==ENEMY_TYPE_GOLIATH?2.0f:1.2f;
        footTravel[i]+=speed*dt;
        if(footTravel[i]>=stride) {
            footTravel[i]=fmodf(footTravel[i],stride);
            float pitch=p->enemyType==ENEMY_TYPE_SWARMER?1.35f:p->enemyType==ENEMY_TYPE_GOLIATH?.65f:1;
            play_sfx_at(SFX_STEP_1+(footVariant[i]++%3),p->pos,i,0,pitch);
        }
    }
}
static bool clarity_shield_projectile(int v) {
    if(!clarity_shield_voxel(v))return false;
    int owner=voxels[v].owner;EnemySpecialState *s=&enemySpecial[owner];
    for(int n=0;n<s->equipmentCount;n++)if(s->equipment[n].id==voxels[v].identity) {
        int g=s->equipment[n].group;s->shieldFlash[g]=.12f;
        if(--s->shieldHits[g]<=0)clarity_shield_melee(v,-1);
        return true;
    }
    return false;
}
static void feedback_shutdown(void) {
    combat_ui_shutdown();
    if(blastShader.id && blastShader.id!=rlGetShaderIdDefault())UnloadShader(blastShader);
    blastShader=(Shader){0};blastShaderAttempted=false;
}
static void feedback_blast(ClarityEffect *e) {
    if(!blastShaderAttempted) {
        blastShaderAttempted=true;
        blastShader=LoadShaderFromMemory(
            "#version 330\n"
            "in vec3 vertexPosition; in vec4 vertexColor; uniform mat4 mvp; out vec3 p; out vec4 color;"
            "void main(){p=vertexPosition;color=vertexColor;gl_Position=mvp*vec4(vertexPosition,1.);}",
            "#version 330\n"
            "in vec3 p; in vec4 color; out vec4 finalColor; uniform float age,smoke; uniform vec3 center,eye;"
            "void main(){vec3 q=p-center;float n=.65+.18*sin(q.x*4.+q.z*2.+age*9.)+.17*sin(q.y*5.-q.x*3.-age*7.);"
            "float rim=pow(1.-abs(dot(normalize(q),normalize(eye-p))),2.);"
            "float a=mix(.15+2.5*rim,n,smoke);finalColor=vec4(color.rgb,color.a*a);}");
        blastAgeLocation=GetShaderLocation(blastShader,"age");
        blastCenterLocation=GetShaderLocation(blastShader,"center");blastEyeLocation=GetShaderLocation(blastShader,"eye");blastSmokeLocation=GetShaderLocation(blastShader,"smoke");
    }
    float t=e->age;rlDrawRenderBatchActive();rlDisableDepthMask();
    bool shader=blastShader.id && blastShader.id!=rlGetShaderIdDefault();
    if(shader){
        Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:players[0].pos;float smoke=0;
        SetShaderValue(blastShader,blastAgeLocation,&t,SHADER_UNIFORM_FLOAT);
        SetShaderValue(blastShader,blastCenterLocation,&e->pos,SHADER_UNIFORM_VEC3);SetShaderValue(blastShader,blastEyeLocation,&eye,SHADER_UNIFORM_VEC3);
        SetShaderValue(blastShader,blastSmokeLocation,&smoke,SHADER_UNIFORM_FLOAT);BeginShaderMode(blastShader);
    }
    if(t<.4f)DrawSphereEx(e->pos,4*fminf(1,t/.15f),12,20,Fade(ORANGE,.18f*(1-t/.4f)));
    if(shader){rlDrawRenderBatchActive();float smoke=1;SetShaderValue(blastShader,blastSmokeLocation,&smoke,SHADER_UNIFORM_FLOAT);}
    for(int n=0;n<5;n++) {
        float a=n*2.39996f;Vector3 p=v_add(e->pos,(Vector3){cosf(a)*t*1.5f,.3f+t,sinf(a)*t*1.5f});
        DrawSphereEx(p,.25f+t*.8f,6,8,Fade(LIGHTGRAY,.2f*fmaxf(0,1-t/.9f)));
    }
    if(shader)EndShaderMode();rlDrawRenderBatchActive();rlEnableDepthMask();
}
static bool clarity_crosscut(void) {
    return currentWorldType==WORLD_TYPE_CROSSCUT && !(gameMode==GAME_MODE_FIREFIGHT && firefightFoundry);
}
static bool clarity_supply_clear(Vector3 at) {
    /* Include bobbing pickup and a full-height approach pocket, using true bounds. */
    BoundingBox clearance={{at.x-.8f,at.y-.55f,at.z-.8f},{at.x+.8f,at.y+1.5f,at.z+.8f}};
    for(int v=0;v<voxel_count;v++) {
        if(voxels[v].isBullet)continue;
        VoxelWorldBounds b;if(voxels[v].simulate)voxel_particle_world_bounds(&voxels[v],&b);else voxel_world_bounds(&voxels[v],&b);
        if(CheckCollisionBoxes(clearance,(BoundingBox){{b.minx,b.miny,b.minz},{b.maxx,b.maxy,b.maxz}}))return false;
    }
    return true;
}
