/* One bounded wave per human slot. Simulation never runs inside a viewport. */
static Vector3 sword_lane_direction(const SwordWave *w,int lane){
    Vector3 f=w->forward,r=v_norm(v_cross(f,(Vector3){0,1,0})),u=v_cross(r,f);
    float birth=SWORD_EMIT_SECONDS+(SWORD_SLASH_END-SWORD_EMIT_SECONDS)*lane/(SWORD_LANES-1.f);
    Vector3 local=weapon_sword_sample(birth,w->direction).direction;
    return v_add(v_add(v_mul(r,local.x),v_mul(u,local.y)),v_mul(f,-local.z));
}
static void eq_sword(int i){
    if(i<0 || i>=MAX_PLAYERS || swordWaves[i].active)return;
    SwordWave *w=&swordWaves[i];
    *w=(SwordWave){.active=true,.sequence=++swordSequences[i],.direction=swordNextDirection[i]};
    swordDirection[i]=w->direction;swordNextDirection[i]^=1;
    for(int n=0;n<SWORD_LANES;n++)w->stop[n]=SWORD_REACH;
    equipment[i].cooldown[EQ_BLADE]=eqSwordCycle;
}
static bool sword_cut_hit(int owner,VoxelHit hit){
    int v=hit.id;if(v<0 || v>=voxel_count || voxel_is_fluid(&voxels[v]))return false;
    uint64_t id=voxels[v].identity;Vector3 at=voxels[v].pos;
    if(clarity_shield_voxel(v)){clarity_shield_melee(v,owner);return enemy_voxel(id)<0;}
    if(!voxels[v].simulate){
        if(foundry_trap_hit(v,owner,false))return true;
        remove_static_voxels_in_region_recycle(hit.x,hit.x,hit.y,hit.y,hit.z,hit.z);
    }else {enemy_detach(v);detach_goliath_armor_voxel(v);remove_voxel_index(v);}
    mark_surface_neighbors(at);
    return enemy_voxel(id)<0;
}
static void sword_step(int owner,SwordWave *w,float previous){
    Vector3 r=v_norm(v_cross(w->forward,(Vector3){0,1,0})),u=v_cross(r,w->forward);
    Vector3 diagonal=v_norm(v_add(v_mul(r,w->direction?-1:1),u)),normal=v_norm(v_cross(w->forward,diagonal));
    for(int lane=0;lane<SWORD_LANES;lane++){
        float t=lane/(SWORD_LANES-1.f);
        float end=weapon_sword_radius(w->age,t);if(end<0)continue;
        end=fminf(end,w->stop[lane]);
        float begin=fmaxf(0,weapon_sword_radius(previous,t));if(end<=begin)continue;
        Vector3 d=sword_lane_direction(w,lane);
        // Three overlapping rays cover the half-metre thick cutting band.
        for(int band=-1;band<=1;band++){
            Vector3 origin=v_add(w->origin,v_mul(normal,band*.20f));
            float cursor=begin;
            for(int attempt=0;attempt<=EQ_SWORD_CUTS && cursor<end;attempt++){
                Ray ray={world_wrap(v_add(origin,v_mul(d,cursor))),d};VoxelHit hit;
                if(!first_voxel_hit_detailed(ray,end-cursor,-1,false,&hit))break;
                float distance=cursor+fmaxf(0,hit.t);
                if(w->cuts>=EQ_SWORD_CUTS || !sword_cut_hit(owner,hit)){
                    end=fmaxf(begin,distance-.002f);w->stop[lane]=fminf(w->stop[lane],end);break;
                }
                w->cuts++;
                if(!w->rewarded){if(players[owner].respawn_timer<=0)add_player_matter(&players[owner],ff_harvest_amount(owner));w->rewarded=true;}
                // Re-query from this point after compaction; overlapping equipment can share a cell.
                cursor=distance;
            }
        }
        // Terrain is removed before checking actors behind it, but only as the
        // front reaches that terrain. A blocked section never damages through cover.
        if(end<=begin)continue;
        Ray ray={world_wrap(v_add(w->origin,v_mul(d,begin))),d};
        for(int j=0;j<activePlayers;j++)if(!w->hit[j] && eq_hostile(owner,j)){
            BoundingBox box=combatant_bounds_at(j,world_image(players[j].pos,ray.position),PLAYER_SIZE*.5f);
            Vector3 pad={.25f,.25f,.25f};box.min=v_sub(box.min,pad);box.max=v_add(box.max,pad);
            RayCollision hit=GetRayCollisionBox(ray,box);
            if(hit.hit && hit.distance>=0 && hit.distance<=end-begin){
                VoxelHit shield;float centerDistance=v_length(world_delta(players[j].pos,ray.position));
                if(first_voxel_hit_detailed(ray,centerDistance,-1,false,&shield))continue;
                w->hit[j]=true;eq_damage(owner,j,eqSwordDamage,true,DEATH_SWORD,w->origin);
            }
        }
    }
}
static void eq_sword_tick(float dt){
    for(int i=0;i<MAX_PLAYERS;i++){
        SwordWave *w=&swordWaves[i];if(!w->active)continue;
        if(!w->emitted && (players[i].respawn_timer>0 || ff_input_blocked(i))){w->active=false;continue;}
        float left=fminf(dt,SWORD_ANIM_SECONDS-w->age);
        while(left>0 && w->active){
            float step=fminf(left,1.f/120),previous=w->age;left-=step;w->age=fminf(SWORD_ANIM_SECONDS,w->age+step);
            if(!w->emitted && w->age>=SWORD_EMIT_SECONDS){
                w->emitted=true;w->origin=world_wrap(players[i].pos);w->forward=player_forward(&players[i]);
                play_sfx_at(SFX_MELEE,w->origin,i,1,.7f);
            }
            if(w->emitted && netTransport.role!=NET_ROLE_CLIENT)sword_step(i,w,previous);
            if(w->age>=SWORD_ANIM_SECONDS-1e-6f){w->age=SWORD_ANIM_SECONDS;w->active=false;players[i].meleeSwingActive=false;}
        }
    }
}
static void eq_sword_draw(void){
    for(int i=0;i<MAX_PLAYERS;i++){
        SwordWave *w=&swordWaves[i];if(!w->active || !w->emitted)continue;
        if(!eqSlashTried){eqSlashTried=true;eqSlashShader=LoadShader(NULL,"shaders/sword_slash.frag");}
        SwordSwingSample sample=weapon_sword_sample(w->age,w->direction);
        Vector3 r=v_norm(v_cross(w->forward,(Vector3){0,1,0})),u=v_cross(r,w->forward);
        Vector3 diagonal=v_norm(v_add(v_mul(r,w->direction?-1:1),u)),normal=v_norm(v_cross(w->forward,diagonal));
        bool shader=eqSlashShader.id && eqSlashShader.id!=rlGetShaderIdDefault();
        if(shader){float time=w->age;SetShaderValue(eqSlashShader,GetShaderLocation(eqSlashShader,"effectTime"),&time,SHADER_UNIFORM_FLOAT);BeginShaderMode(eqSlashShader);}
        BeginBlendMode(BLEND_ADDITIVE);rlDrawRenderBatchActive();rlDisableDepthMask();rlDisableBackfaceCulling();
        rlBegin(RL_TRIANGLES);
        for(int lane=0;lane<SWORD_LANES-1;lane++){
            Vector3 edge[2][2];bool visible=true;
            for(int k=0;k<2;k++){
                int n=lane+k;float t=n/(SWORD_LANES-1.f),radius=weapon_sword_radius(w->age,t);
                if(radius<0){visible=false;break;}
                radius=fminf(radius,w->stop[n]);if(radius<.01f){visible=false;break;}
                Vector3 d=sword_lane_direction(w,n);
                // Attach newborn energy to the rendered blade, then shed the
                // viewmodel offset as it expands to the authoritative world front.
                float grow=clampf((radius-1)/(SWORD_REACH-1),0,1);
                WeaponMeleePose melee={.progress=(SWORD_EMIT_SECONDS+(SWORD_SLASH_END-SWORD_EMIT_SECONDS)*t)/SWORD_ANIM_SECONDS,.sword_direction=w->direction};
                WeaponPose blade=weapon_pose(WEAPON_SWORD,w->origin,w->forward,r,u,i==drawViewPlayerIndex,weaponViewAspect,&melee,0);
                if(i!=drawViewPlayerIndex && player_uses_rat(i))blade.muzzle=v_add(w->origin,v_mul(v_sub(blade.muzzle,w->origin),.42f*CHARACTER_MODEL_SCALE));
                Vector3 offset=v_mul(v_sub(blade.muzzle,v_add(w->origin,d)),1-grow);
                Vector3 front=v_add(v_add(w->origin,v_mul(d,radius)),offset);
                float width=.25f*clampf((w->age-SWORD_EMIT_SECONDS)*20,0,1);
                edge[k][0]=v_add(front,v_mul(normal,-width));
                edge[k][1]=v_add(front,v_mul(normal,width));
            }
            if(!visible)continue;
            Vector3 v[6]={edge[0][0],edge[0][1],edge[1][1],edge[0][0],edge[1][1],edge[1][0]};
            float a=lane/(SWORD_LANES-1.f),b=(lane+1)/(SWORD_LANES-1.f);
            Vector2 uv[6]={{a,0},{a,1},{b,1},{a,0},{b,1},{b,0}};
            for(int k=0;k<6;k++){rlColor4ub(35,135,255,(unsigned char)(220*sample.fade));rlTexCoord2f(uv[k].x,uv[k].y);rlVertex3f(v[k].x,v[k].y,v[k].z);}
        }
        rlEnd();rlDrawRenderBatchActive();rlEnableBackfaceCulling();rlEnableDepthMask();EndBlendMode();if(shader)EndShaderMode();
    }
}
