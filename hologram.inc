/* Host simulation and presentation, separate from real combatant slots. */
static float hologram_max_range(void) {return 2.f*sqrtf(2.f)*FLOOR_SIZE;}
static float hologram_max_life(void) {return 8.f+hologram_max_range()/MOVE_SPEED;}
static bool hologram_deploy(int owner) {
    aim_cancel(owner);
    if(netTransport.role==NET_ROLE_CLIENT || owner<0 || owner>=MAX_PLAYERS || owner>=activePlayers ||
       !player_uses_rat(owner) || debugConsoleOpen || gameState!=GAME_STATE_PLAYING || creativeModeActive ||
       players[owner].respawn_timer>0 || combatantLife[owner]!=LIFE_ALIVE || ff_input_blocked(owner) ||
       (gameMode==GAME_MODE_FIREFIGHT && firefightWaveState==FIREFIGHT_STATE_INTERMISSION && !enemyLab) ||
       (enemyLab && enemyLabPaused) || holograms[owner].cooldown>0)return false;
    Hologram *h=&holograms[owner];*h=(Hologram){0};
    h->body=players[owner];h->body.vel=(Vector3){0};h->body.pitch=0;
    h->body.isExposed=false;h->body.invuln_timer=0;h->body.respawn_timer=0;
    h->health=50;h->life=8;h->cooldown=15;h->generation=++hologramGeneration;
    // Sample the owner's full aim once, before flattening the walking direction.
    // Empty sky has a finite destination; floor and voxel surfaces shorten it.
    Ray aim={players[owner].pos,player_forward(&players[owner])};
    float distance=hologram_max_range();
    if(!worldPeriodic) {
        float limit=fmaxf(0,FLOOR_SIZE-.4f);
        if(fabsf(aim.direction.x)>.0001f)
            distance=fminf(distance,fmaxf(0,(copysignf(limit,aim.direction.x)-aim.position.x)/aim.direction.x));
        if(fabsf(aim.direction.z)>.0001f)
            distance=fminf(distance,fmaxf(0,(copysignf(limit,aim.direction.z)-aim.position.z)/aim.direction.z));
    }
    if(aim.direction.y<-.0001f)distance=fminf(distance,fmaxf(0,-aim.position.y/aim.direction.y));
    VoxelHit surface;
    if(first_voxel_hit_detailed(aim,distance,-1,false,&surface))distance=surface.t;
    h->direction=aim.direction;h->direction.y=0;
    h->travelRemaining=distance*v_length(h->direction);
    h->life=8.f+h->travelRemaining/MOVE_SPEED;
    h->direction=v_norm(h->direction);
    h->color=player_palette_color(owner);h->weapon=player_weapon_kind(&players[owner]);h->armed=player_has_gun(owner);
    h->visualTarget=h->body.pos;
    h->pose=character_pose(h->body.pos,h->body.yaw,0,0,h->body.onGround,false,h->color,(float)GetTime());
    for(int v=0;v<voxel_count;v++)voxels[v].hologramHitMask&=~(1u<<owner);
    return true;
}
static void hologram_damage(int target,float damage) {
    if(netTransport.role==NET_ROLE_CLIENT || !target_is_hologram(target) || !isfinite(damage) || damage<=0)return;
    Hologram *h=&holograms[target-MAX_COMBATANTS];if(h->life<=0)return;
    h->damageFlicker=.4f;
    h->health=fmaxf(0,h->health-damage);if(h->health<=0)h->life=0;
}
static BoundingBox hologram_bounds(int n,Vector3 origin) {
    return combatant_bounds_at(n,world_image(holograms[n].body.pos,origin),PLAYER_SIZE*.5f);
}
static void hologram_projectile(int index,Vector3 start,Vector3 end,float damage) {
    if(netTransport.role==NET_ROLE_CLIENT)return;
    Voxel *v=&voxels[index];Vector3 path=v_sub(end,start);float distance=v_length(path);
    if(distance<.00001f)return;
    VoxelHit wall;Ray ray={start,v_mul(path,1/distance)};
    bool checkedWorld=false;
    for(int n=0;n<MAX_PLAYERS;n++)if(holograms[n].life>0 && !(v->hologramHitMask&(1u<<n)) && v->owner!=n) {
        RayCollision hit=GetRayCollisionBox(ray,expand_bounds(hologram_bounds(n,start),v->rest_edge*.5f));
        if(!hit.hit || hit.distance>distance)continue;
        if(!checkedWorld){if(first_voxel_hit_detailed(ray,distance,index,false,&wall))distance=fminf(distance,wall.t);checkedWorld=true;}
        if(hit.distance<=distance){if(!v->isBullet && !v->foundryTrapHeavy)damage=VOXEL_DAMAGE*fmaxf(1,build_glue_cluster_indices(index,glueClusterIndices));v->hologramHitMask|=1u<<n;hologram_damage(MAX_COMBATANTS+n,damage);}
    }
}
static void hologram_melee(int owner,Vector3 start,Vector3 end) {
    for(int n=0;n<MAX_PLAYERS;n++)if(n!=owner && holograms[n].life>0 && ai_hostile(owner,MAX_COMBATANTS+n)) {
        BoundingBox box=hologram_bounds(n,start);
        if(segment_intersects_aabb(start,end,box.min,box.max) && !is_view_occluded_by_voxels(start,holograms[n].body.pos))
            hologram_damage(MAX_COMBATANTS+n,MATTER_MELEE_DAMAGE);
    }
}
static void hologram_blast(Vector3 pos,float radius,float damage) {
    for(int n=0;n<MAX_PLAYERS;n++)if(holograms[n].life>0) {
        float distance=v_length(world_delta(holograms[n].body.pos,pos));
        if(distance<radius && !is_view_occluded_by_voxels(pos,holograms[n].body.pos))hologram_damage(MAX_COMBATANTS+n,damage*(1-distance/radius));
    }
}
static void hologram_tick(float dt) {
    if(enemyLab && enemyLabPaused && !enemyLabStep)return;
    dt=clampf(dt,0,.25f);
    for(int n=0;n<MAX_PLAYERS;n++) {
        Hologram *h=&holograms[n];h->cooldown=fmaxf(0,h->cooldown-dt);
        h->damageFlicker=fmaxf(0,h->damageFlicker-dt);
        if(h->life<=0)continue;
        if(netTransport.role==NET_ROLE_CLIENT){
            h->life=fmaxf(0,h->life-dt);
            h->body.pos=v_add(h->body.pos,v_mul(world_delta(h->visualTarget,h->body.pos),1-expf(-dt*18)));
        } else {
            if(n>=activePlayers || players[n].respawn_timer>0 || combatantLife[n]!=LIFE_ALIVE ||
               (netTransport.role==NET_ROLE_HOST && n>0 && !netPlayerPresent[n])){h->life=0;continue;}
            h->life=fmaxf(0,h->life-dt);if(!h->life)continue;
            float remaining=dt;
            while(remaining>0){
                float step=fminf(remaining,1.f/60);remaining-=step;
                float travel=fminf(h->travelRemaining,MOVE_SPEED*step);
                h->body.vel.x=h->direction.x*travel/step;h->body.vel.z=h->direction.z*travel/step;
                h->body.vel.y-=GRAVITY*step;
                bool blocked=human_body_sweep(&h->body,n,0,h->body.vel.x*step);
                blocked|=human_body_sweep(&h->body,n,2,h->body.vel.z*step);
                human_body_sweep(&h->body,n,1,h->body.vel.y*step);
                h->travelRemaining=fmaxf(0,h->travelRemaining-travel);
                if(blocked || h->travelRemaining<=.0001f){
                    h->travelRemaining=0;h->direction=(Vector3){0};
                    h->body.vel.x=h->body.vel.z=0;
                    h->life=fminf(h->life,8.f);
                }
                h->body.onGround=h->body.vel.y==0;
                h->body.pos=world_wrap(h->body.pos);
                if(!worldPeriodic){h->body.pos.x=clampf(h->body.pos.x,-FLOOR_SIZE+.4f,FLOOR_SIZE-.4f);h->body.pos.z=clampf(h->body.pos.z,-FLOOR_SIZE+.4f,FLOOR_SIZE-.4f);}
            }
            for(int i=0;i<activePlayers;i++)if(is_player_bot(i) && players[i].respawn_timer<=0 && ai_hostile(i,MAX_COMBATANTS+n) && players[i].contactDamageTimer<=0) {
                bool swarm=players[i].enemyType==ENEMY_TYPE_SWARMER,goliath=players[i].enemyType==ENEMY_TYPE_GOLIATH;
                if((swarm || goliath) && v_length(world_delta(players[i].pos,h->body.pos))<(swarm?.75f:1.6f)) {
                    hologram_damage(MAX_COMBATANTS+n,swarm?10:30);players[i].contactDamageTimer=swarm?.5f:.8f;
                }
            }
        }
        float speed=sqrtf(h->body.vel.x*h->body.vel.x+h->body.vel.z*h->body.vel.z);
        h->pose=character_pose(h->body.pos,h->body.yaw,0,speed,h->body.onGround,false,h->color,(float)GetTime());
        if(h->animation.initialized && v_length(v_sub(h->animation.previous,h->body.pos))>10)memset(&h->animation,0,sizeof(h->animation));
        character_animate(&h->animation,&h->pose,dt,character_ground_query,(void*)(intptr_t)n);
    }
}
static void hologram_draw(bool emission) {
    for(int n=0;n<MAX_PLAYERS;n++) {
        Hologram *h=&holograms[n];if(h->life<=0 || !periodic_visible(h->body.pos,3))continue;
        // Drop the complete projection in brief pulses, including its weapon and bloom.
        // The simulation timer keeps every viewport and render pass in agreement.
        if(h->damageFlicker>0 && ((int)((.4f-h->damageFlicker)*20)%2)==0)continue;
        if(drawViewPlayerIndex==n && v_length(world_delta(h->body.pos,players[n].pos))<.75f)continue;
        character_draw(h->pose,emission);
        if(h->armed && drawWeaponModels){
            Vector3 f=player_forward(&h->body),r=v_norm(v_cross(f,(Vector3){0,1,0})),u=v_cross(r,f);
            WeaponMeleePose melee={.progress=-1};WeaponVisual visual={0};
            WeaponPose wp=weapon_pose(h->weapon,h->body.pos,f,r,u,false,1,&melee,0);
            Matrix scale=MatrixMultiply(MatrixTranslate(-h->body.pos.x,-h->body.pos.y,-h->body.pos.z),MatrixScale(.42f*CHARACTER_MODEL_SCALE,.42f*CHARACTER_MODEL_SCALE,.42f*CHARACTER_MODEL_SCALE));
            scale=MatrixMultiply(scale,MatrixTranslate(h->body.pos.x,h->body.pos.y,h->body.pos.z));
            wp.transform=MatrixMultiply(wp.transform,scale);
            weapon_draw(wp,&visual,drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:h->body.pos,h->pose.time,false,emission?WEAPON_EMISSION:WEAPON_COLOR);
            if(!emission)weapon_draw(wp,&visual,drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:h->body.pos,h->pose.time,false,WEAPON_GLASS);
            character_draw_paws(wp,h->color,emission);
            if(!emission)for(int side=-1;side<=1;side+=2){
                Vector3 shoulder=Vector3Transform((Vector3){side*.17f,.11f,0},character_body_transform(h->pose));
                Vector3 hand=Vector3Transform(side>0?(Vector3){.17f,-.67f,.40f}:(Vector3){-.36f,-.33f,h->weapon==WEAPON_DYNAMIC_LAUNCHER?-.55f:-.30f},wp.transform);
                character_draw_arm(shoulder,hand,h->color);
            }
        }else if(!emission)character_draw_idle_arms(h->pose);
        if(!emission)character_draw_visor(h->pose);
    }
}
static void hologram_hud(int owner,int w,int h) {
    if(owner<0 || owner>=MAX_PLAYERS)return;
    Hologram *d=&holograms[owner];const char *key=playerInput[owner]==INPUT_TYPE_GAMEPAD?"LB":owner==0?"V":"/";
    const char *label=d->life>0?TextFormat("Decoy %.0f/50 | %.0fs",d->health,d->life):d->cooldown>0?TextFormat("Decoy %s: %.0fs",key,ceilf(d->cooldown)):TextFormat("Decoy %s: READY",key);
    DrawText(label,12,h-180,w<700?12:16,d->life>0?SKYBLUE:LIGHTGRAY);
}
static void hologram_net_snapshot(int peer) {
    if(netTransport.role!=NET_ROLE_HOST)return;
    uint8_t packet[512];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
    net_write_header(&w,NET_MSG_HOLOGRAMS,0,netTransport.session_id,netServerTick);
    for(int n=0;n<3;n++)net_write_u8(&w,teamUpgradeLevels[n]);
    for(int n=0;n<MAX_PLAYERS;n++){
        Hologram *h=&holograms[n];net_write_u32(&w,h->generation);
        net_write_f32(&w,h->life);net_write_f32(&w,h->health);net_write_f32(&w,h->cooldown);
        net_write_f32(&w,h->body.pos.x);net_write_f32(&w,h->body.pos.y);net_write_f32(&w,h->body.pos.z);
        net_write_f32(&w,h->body.vel.x);net_write_f32(&w,h->body.vel.y);net_write_f32(&w,h->body.vel.z);net_write_f32(&w,h->body.yaw);
        net_write_u8(&w,h->color.r);net_write_u8(&w,h->color.g);net_write_u8(&w,h->color.b);
        net_write_u8(&w,h->weapon);net_write_u8(&w,h->armed);net_write_u8(&w,h->body.onGround);
    }
    if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);
    else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
}
static bool hologram_net_receive(NetMessageType type,NetReader *r,uint32_t tick) {
    if(type!=NET_MSG_HOLOGRAMS)return false;
    if(netTransport.role!=NET_ROLE_CLIENT || tick<hologramReceiveTick)return true;
    Hologram incoming[MAX_PLAYERS]={0};unsigned levels[3];
    for(int n=0;n<3;n++){levels[n]=net_read_u8(r);if(levels[n]>teamUpgradeCaps[n])return true;}
    for(int n=0;n<MAX_PLAYERS;n++){
        Hologram *h=&incoming[n];h->generation=net_read_u32(r);h->life=net_read_f32(r);h->health=net_read_f32(r);h->cooldown=net_read_f32(r);
        h->body.pos=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        h->body.vel=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};h->body.yaw=net_read_f32(r);
        h->color=(Color){net_read_u8(r),net_read_u8(r),net_read_u8(r),255};
        h->weapon=net_read_u8(r);h->armed=net_read_u8(r)!=0;h->body.onGround=net_read_u8(r)!=0;
        if(r->failed || !v_isfinite(h->body.pos) || !v_isfinite(h->body.vel) || !isfinite(h->body.yaw) ||
           !isfinite(h->life) || h->life<0 || h->life>hologram_max_life()+.01f || !isfinite(h->health) || h->health<0 || h->health>50 ||
           !isfinite(h->cooldown) || h->cooldown<0 || h->cooldown>15 || h->weapon>WEAPON_DYNAMIC_LAUNCHER)return true;
    }
    memcpy(teamUpgradeLevels,levels,sizeof(levels));
    for(int n=0;n<MAX_PLAYERS;n++){
        Hologram old=holograms[n];holograms[n]=incoming[n];holograms[n].visualTarget=incoming[n].body.pos;
        holograms[n].pose=character_pose(incoming[n].body.pos,incoming[n].body.yaw,0,0,incoming[n].body.onGround,false,incoming[n].color,(float)GetTime());
        if(old.generation==incoming[n].generation && old.life>0){
            holograms[n].body.pos=old.body.pos;holograms[n].animation=old.animation;holograms[n].pose=old.pose;
            holograms[n].damageFlicker=incoming[n].health<old.health?.4f:old.damageFlicker;
        }
    }
    hologramReceiveTick=tick;return true;
}
