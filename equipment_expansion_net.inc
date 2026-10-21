/* Independently versioned chunks: at most 16 effects per datagram. Empty slots
 * still transmit their ID so removals converge after packet loss and late joins. */
static void ex_snapshot(int peer){
    if(netTransport.role!=NET_ROLE_HOST)return;
    for(int part=0;part<13;part++){
        uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_EXPANSION,0,netTransport.session_id,netServerTick);net_write_u8(&w,part);
        if(part==0){
            for(int i=0;i<MAX_PLAYERS;i++){net_write_u8(&w,exSelectedGrenade[i]);net_write_u8(&w,exGrenades[i][0]);net_write_u8(&w,exGrenades[i][1]);net_write_f32(&w,exFailedBuild[i]);}
            for(int i=0;i<MAX_COMBATANTS;i++)net_write_f32(&w,exSuppression[i]);
            for(int n=0;n<EX_FIELDS;n++){GravityField *f=&exFields[n];net_write_u32(&w,f->id);if(f->id){net_write_u8(&w,f->owner<0?255:f->owner);eq_write_vec(&w,f->pos);net_write_f32(&w,f->age);}}
            for(int n=0;n<8;n++){IglooReveal *r=&exReveals[n];net_write_u32(&w,r->id);if(r->id){eq_write_vec(&w,r->origin);net_write_f32(&w,r->age);net_write_u8(&w,r->rotation);}}
        }else if(part<=8){for(int n=(part-1)*16;n<part*16;n++){GlobProjectile *g=&exGlobs[n];net_write_u32(&w,g->id);if(g->id){net_write_u8(&w,g->owner);eq_write_vec(&w,g->pos);eq_write_vec(&w,g->velocity);net_write_f32(&w,g->age);}}}
        else {for(int n=(part-9)*16;n<(part-8)*16;n++){ExpansionSplat *s=&exSplats[n];net_write_u32(&w,s->id);if(s->id){net_write_u8(&w,s->actor<0?255:s->actor);net_write_u64(&w,s->support);eq_write_vec(&w,s->pos);eq_write_vec(&w,s->local);net_write_f32(&w,s->age);net_write_u8(&w,s->impact);}}}
        if(w.failed){TraceLog(LOG_ERROR,"Expansion packet overflow part %d",part);continue;}
        if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    }
}
static bool ex_receive(NetMessageType type,NetReader *r,uint32_t tick){
    if(type!=NET_MSG_EXPANSION)return false;if(netTransport.role!=NET_ROLE_CLIENT)return true;
    int part=net_read_u8(r);if(part>12||tick<exPartTick[part])return true;
    if(part==0){
        int selection[4],counts[4][2];float failed[4],suppression[MAX_COMBATANTS];GravityField fields[EX_FIELDS]={0};IglooReveal reveals[8]={0};
        for(int i=0;i<MAX_PLAYERS;i++){selection[i]=net_read_u8(r);counts[i][0]=net_read_u8(r);counts[i][1]=net_read_u8(r);failed[i]=net_read_f32(r);if(selection[i]>=3||counts[i][0]>2||counts[i][1]>2||!isfinite(failed[i])||failed[i]<0||failed[i]>2)return true;}
        for(int i=0;i<MAX_COMBATANTS;i++){suppression[i]=net_read_f32(r);if(!isfinite(suppression[i])||suppression[i]<0||suppression[i]>5)return true;}
        for(int n=0;n<EX_FIELDS;n++){GravityField *f=&fields[n];f->id=net_read_u32(r);if(f->id){int owner=net_read_u8(r);f->owner=owner==255?-1:owner;f->pos=eq_read_vec(r);f->age=net_read_f32(r);if(f->owner>=MAX_COMBATANTS||!v_isfinite(f->pos)||!isfinite(f->age)||f->age<0||f->age>5)return true;}}
        for(int n=0;n<8;n++){IglooReveal *v=&reveals[n];v->id=net_read_u32(r);if(v->id){v->origin=eq_read_vec(r);v->age=net_read_f32(r);v->rotation=net_read_u8(r);if(!v_isfinite(v->origin)||!isfinite(v->age)||v->age<0||v->age>1.2f||v->rotation>3)return true;}}
        if(r->failed)return true;memcpy(exSelectedGrenade,selection,sizeof(selection));memcpy(exGrenades,counts,sizeof(counts));memcpy(exFailedBuild,failed,sizeof(failed));memcpy(exSuppression,suppression,sizeof(suppression));memcpy(exFields,fields,sizeof(fields));memcpy(exReveals,reveals,sizeof(reveals));
    }else if(part<=8){GlobProjectile incoming[16]={0};for(int n=0;n<16;n++){GlobProjectile *g=&incoming[n];g->id=net_read_u32(r);if(g->id){g->owner=net_read_u8(r);g->pos=eq_read_vec(r);g->velocity=eq_read_vec(r);g->age=net_read_f32(r);if(g->owner>=MAX_COMBATANTS||!v_isfinite(g->pos)||!v_isfinite(g->velocity)||v_length(g->velocity)>61||!isfinite(g->age)||g->age<0||g->age>1)return true;}}if(r->failed)return true;memcpy(exGlobs+(part-1)*16,incoming,sizeof(incoming));}
    else {ExpansionSplat incoming[16]={0};for(int n=0;n<16;n++){ExpansionSplat *s=&incoming[n];s->id=net_read_u32(r);if(s->id){int actor=net_read_u8(r);s->actor=actor==255?-1:actor;s->support=net_read_u64(r);s->pos=eq_read_vec(r);s->local=eq_read_vec(r);s->age=net_read_f32(r);int impact=net_read_u8(r);s->impact=impact;if(s->actor>=MAX_COMBATANTS||!v_isfinite(s->pos)||!v_isfinite(s->local)||!isfinite(s->age)||s->age<0||s->age>5||impact>1)return true;}}if(r->failed)return true;memcpy(exSplats+(part-9)*16,incoming,sizeof(incoming));}
    exPartTick[part]=tick;return true;
}
static void ex_pickup_snapshot(int peer){
    for(int start=0;start<MAX_PICKUPS;start+=48){uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_WORLD_PICKUPS,0,netTransport.session_id,netServerTick);int count=MAX_PICKUPS-start;if(count>48)count=48;net_write_u8(&w,start);net_write_u8(&w,count);
        for(int i=start;i<start+count;i++){eq_write_vec(&w,pickups[i].pos);net_write_f32(&w,pickups[i].respawnTimer);net_write_u8(&w,pickups[i].type);net_write_u8(&w,pickups[i].active);}
        if(w.failed){TraceLog(LOG_ERROR,"Pickup packet overflow");continue;}if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    }
}
