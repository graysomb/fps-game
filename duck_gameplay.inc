/* Ambient actors never enter combatant, voxel, pickup, score, or AI pools. */
static bool ducks_enabled(void) {
    return !useCustomMap && !creativeModeActive && !enemyLab && !firefightPractice &&
        ((gameMode==GAME_MODE_FIREFIGHT && firefightFoundry) ||
         (gameMode!=GAME_MODE_FIREFIGHT && currentWorldType==WORLD_TYPE_CROSSCUT));
}
static BoundingBox duck_bounds(Vector3 p,float pad) {
    return (BoundingBox){v_add(p,(Vector3){-.43f-pad,.025f-pad,-.60f-pad}),
                         v_add(p,(Vector3){.43f+pad,.68f+pad,.60f+pad})};
}
static bool duck_clear(Vector3 p,bool actors) {
    BoundingBox box=duck_bounds(p,.06f);box.min.y=.025f;
    for(int z=(int)floorf(box.min.z/VOXEL_SIZE);z<=(int)floorf(box.max.z/VOXEL_SIZE);z++)
    for(int y=0;y<=(int)floorf(box.max.y/VOXEL_SIZE);y++)
    for(int x=(int)floorf(box.min.x/VOXEL_SIZE);x<=(int)floorf(box.max.x/VOXEL_SIZE);x++)
        if(table_get(x,y,z)>=0)return false;
    if(actors)for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0 && combatantLife[i]==LIFE_ALIVE &&
        CheckCollisionBoxes(box,combatant_bounds_at(i,world_image(players[i].pos,p),PLAYER_SIZE*.5f)))return false;
    return true;
}
static void ducks_reset(void) {
    memset(mapDucks,0,sizeof(mapDucks));ducksInitialized=false;duckClock=0;duckReceiveTick=0;
}
static Vector3 duck_route(MapDuck *d,float phase) {
    return world_wrap(v_add(d->home,(Vector3){DUCK_ROUTE_RADIUS*cosf(phase),0,DUCK_ROUTE_RADIUS*sinf(phase)}));
}
static void ducks_initialize(void) {
    if(ducksInitialized || !ducks_enabled())return;
    for(int n=0;n<MAP_DUCKS;n++) {
        // Four distinct open stretches between repeated arenas, away from the
        // centerline machinery. Swim across seams; keep bounded maps inside walls.
        float edge=worldFloorHalfSize-(worldPeriodic?1.0f:DUCK_ROUTE_RADIUS+2.0f);
        Vector3 home=crosscut_rotate((Vector3){edge,0,worldFloorHalfSize*.4f},n);
        MapDuck *d=&mapDucks[n];*d=(MapDuck){.home=home,.phase=n*PI*.5f};
        d->pos=duck_route(d,d->phase);d->yaw=180-d->phase*RAD2DEG;
        d->alive=duck_clear(d->pos,true);
        if(!d->alive)d->respawn=1; // Wait for a safe space instead of spawning in cover.
    }
    ducksInitialized=true;
}
static void duck_pop(int n) {
    if(netTransport.role==NET_ROLE_CLIENT || !mapDucks[n].alive)return;
    MapDuck *d=&mapDucks[n];d->alive=false;d->respawn=DUCK_RESPAWN_SECONDS;d->puff=d->pos;++d->bursts;
}
static bool ducks_ray(Ray ray,float distance,float radius,int ignore,bool penetrate) {
    if(netTransport.role==NET_ROLE_CLIENT || !ducksInitialized || !ducks_enabled() || distance<=0)return false;
    bool any=false;float period=worldFloorHalfSize*2;
    int images=worldPeriodic && distance>period*.5f?(int)ceilf(distance/period):0;
    for(int n=0;n<MAP_DUCKS;n++)if(mapDucks[n].alive) {
        Vector3 nearest=world_image(mapDucks[n].pos,ray.position);
        for(int x=-images;x<=images && mapDucks[n].alive;x++)for(int z=-images;z<=images && mapDucks[n].alive;z++) {
            Vector3 at=v_add(nearest,(Vector3){x*period,0,z*period});
            BoundingBox box=duck_bounds(at,radius);
            Vector3 end=v_add(ray.position,v_mul(ray.direction,distance));
            BoundingBox sweep={{fminf(end.x,ray.position.x),fminf(end.y,ray.position.y),fminf(end.z,ray.position.z)},
                               {fmaxf(end.x,ray.position.x),fmaxf(end.y,ray.position.y),fmaxf(end.z,ray.position.z)}};
            if(!CheckCollisionBoxes(sweep,box))continue;
            bool inside=CheckCollisionBoxSphere(box,ray.position,0);
            RayCollision hit=GetRayCollisionBox(ray,box);float reach=inside?0:hit.distance;
            if(!inside && (!hit.hit || reach<0 || reach>distance))continue;
            VoxelHit wall;
            if(!penetrate && reach>.001f && first_voxel_hit_detailed(ray,reach,ignore,false,&wall) && wall.t<reach-.001f)continue;
            duck_pop(n);any=true;
        }
    }
    return any;
}
static void ducks_segment(Vector3 start,Vector3 end,float radius,int ignore) {
    Vector3 delta=world_delta(end,start);float length=v_length(delta);
    if(length>1e-6f)ducks_ray((Ray){start,v_mul(delta,1/length)},length,radius,ignore,false);
}
static void ducks_blast(Vector3 pos,float radius) {
    if(netTransport.role==NET_ROLE_CLIENT || !ducksInitialized || !ducks_enabled())return;
    for(int n=0;n<MAP_DUCKS;n++)if(mapDucks[n].alive) {
        Vector3 delta=world_delta(v_add(mapDucks[n].pos,(Vector3){0,.3f,0}),pos);float distance=v_length(delta);VoxelHit hit;
        if(distance<radius && (distance<.01f || !first_voxel_hit_detailed((Ray){pos,v_norm(delta)},fmaxf(0,distance-.2f),-1,false,&hit)))duck_pop(n);
    }
}
static void ducks_tick(float dt) {
    if(!ducks_enabled() || gameState!=GAME_STATE_PLAYING || dt<=0)return;
    duckClock+=dt;
    if(netTransport.role==NET_ROLE_CLIENT){
        for(int n=0;n<MAP_DUCKS;n++)if(!mapDucks[n].alive)mapDucks[n].respawn=fmaxf(0,mapDucks[n].respawn-dt);
        return; // The host alone decides hits and reappearance.
    }
    ducks_initialize();
    for(int n=0;n<MAP_DUCKS;n++) {
        MapDuck *d=&mapDucks[n];
        if(!d->alive) {
            d->respawn=fmaxf(0,d->respawn-dt);
            if(d->respawn<=0) {
                Vector3 at=duck_route(d,d->phase);
                if(duck_clear(at,true)){d->pos=at;d->alive=true;}
                else d->respawn=1; // A player-built wall cannot trap a respawning duck.
            }
            continue;
        }
        for(int i=0;i<activePlayers;i++)if(players[i].respawn_timer<=0 && combatantLife[i]==LIFE_ALIVE &&
            CheckCollisionBoxes(duck_bounds(d->pos,0),combatant_bounds_at(i,world_image(players[i].pos,d->pos),PLAYER_SIZE*.5f)))duck_pop(n);
        if(!d->alive)continue;
        if(!duck_clear(d->pos,false)){duck_pop(n);continue;}
        float next=d->phase+dt*(DUCK_SWIM_SPEED/DUCK_ROUTE_RADIUS);Vector3 at=duck_route(d,next);
        if(duck_clear(at,false)){d->pos=at;d->phase=fmodf(next,2*PI);d->yaw=180-d->phase*RAD2DEG;}
    }
}
static void duck_feathers(MapDuck *d,int index) {
    float age=DUCK_RESPAWN_SECONDS-d->respawn;
    if(!d->bursts || age<0 || age>2.6f)return;
    for(int k=0;k<28;k++) {
        float seed=k*2.399963f+index*1.7f+d->bursts*.71f;
        float spread=.5f+(k%7)*.14f,flight=1-expf(-age*2.5f);
        Vector3 p=v_add(d->puff,(Vector3){cosf(seed)*spread*flight,.3f+(1.2f+(k%5)*.16f)*age-1.1f*age*age,sinf(seed)*spread*flight});
        p.y=fmaxf(.025f,p.y);p.x+=sinf(age*9+seed)*.04f*age;
        float size=.09f*(1-smoothstepf(clampf((age-1.7f)/.9f,0,1)));
        Vector3 axis={cosf(seed+age*4)*size,sinf(seed+age*4)*size,.03f*sinf(seed)},side={-axis.y*.35f,axis.x*.35f,size*.2f};
        Vector3 a=v_add(p,axis),b=v_add(p,side),c=v_sub(p,axis),e=v_sub(p,side);
        Color color=k%6==0?(Color){43,118,74,255}:k%5==0?(Color){91,79,58,255}:(Color){231,228,208,255};
        DrawTriangle3D(a,b,c,color);DrawTriangle3D(a,c,e,color);
        DrawTriangle3D(c,b,a,color);DrawTriangle3D(e,c,a,color);
    }
}
static void ducks_draw(void) {
    if(!ducksInitialized || !ducks_enabled())return;
    duck_renderer_init(NULL);
    for(int n=0;n<MAP_DUCKS;n++) {
        MapDuck *d=&mapDucks[n];if(!periodic_visible(d->alive?d->pos:d->puff,3))continue;
        if(d->alive){DuckPose p={.position=d->pos,.yaw=d->yaw,.scale=.75f,.time=duckClock,.speed=DUCK_SWIM_SPEED,.phase=n*1.7f};duck_draw_floor(p,true);duck_draw(p);}
        else duck_feathers(d,n);
    }
}
static void ducks_write(NetWriter *w) {
    net_write_u8(w,ducksInitialized && ducks_enabled());net_write_f32(w,duckClock);
    for(int n=0;n<MAP_DUCKS;n++) {
        MapDuck *d=&mapDucks[n];net_write_u8(w,d->alive);net_write_u32(w,d->bursts);
        net_write_f32(w,d->pos.x);net_write_f32(w,d->pos.y);net_write_f32(w,d->pos.z);
        net_write_f32(w,d->puff.x);net_write_f32(w,d->puff.y);net_write_f32(w,d->puff.z);
        net_write_f32(w,d->yaw);net_write_f32(w,d->respawn);
    }
}
static void ducks_snapshot(int peer) {
    if(netTransport.role!=NET_ROLE_HOST)return;
    uint8_t packet[256];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
    net_write_header(&w,NET_MSG_DUCKS,0,netTransport.session_id,netServerTick);ducks_write(&w);
    if(w.failed)return;
    if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);
    else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
}
static bool ducks_receive(NetMessageType type,NetReader *r,uint32_t tick) {
    if(type!=NET_MSG_DUCKS)return false;
    if(netTransport.role!=NET_ROLE_CLIENT || tick<duckReceiveTick)return true;
    int initialized=net_read_u8(r);float clock=net_read_f32(r);MapDuck incoming[MAP_DUCKS]={0};
    for(int n=0;n<MAP_DUCKS;n++) {
        MapDuck *d=&incoming[n];int alive=net_read_u8(r);d->alive=alive;d->bursts=net_read_u32(r);
        d->pos=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        d->puff=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        d->yaw=net_read_f32(r);d->respawn=net_read_f32(r);
        if(alive>1 || !v_isfinite(d->pos)||!v_isfinite(d->puff)||!isfinite(d->yaw)||!isfinite(d->respawn)||d->respawn<0||d->respawn>DUCK_RESPAWN_SECONDS)return true;
    }
    if(r->failed || initialized>1 || !isfinite(clock)||clock<0)return true;
    memcpy(mapDucks,incoming,sizeof(mapDucks));ducksInitialized=initialized;duckClock=clock;duckReceiveTick=tick;return true;
}
