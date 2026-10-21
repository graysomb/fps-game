/* Fixed two-part snapshots fit the existing 1100-byte transport budget. */
static uint32_t claritySnapshotSequence,clarityReceiveSequence,clarityPendingSequence,clarityLastEvent;
static unsigned clarityReceiveMask;
static GrenadeState clarityPending[CLARITY_GRENADES];
static void clarity_net_reset(void){clarityReceiveSequence=clarityPendingSequence=clarityLastEvent=0;clarityReceiveMask=0;}
static void clarity_net_snapshot(int peer) {
    if(netTransport.role!=NET_ROLE_HOST)return;
    uint32_t sequence=++claritySnapshotSequence;
    for(int part=0;part<2;part++) {
        uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_GRENADES,0,netTransport.session_id,netServerTick);
        net_write_u32(&w,sequence);net_write_u8(&w,part);
        for(int n=part*16;n<(part+1)*16;n++){GrenadeState empty={0},*g=grenades[n].id?&grenades[n]:&empty;net_write_u32(&w,g->id);net_write_u8(&w,g->owner<0?255:g->owner);net_write_f32(&w,g->pos.x);net_write_f32(&w,g->pos.y);net_write_f32(&w,g->pos.z);net_write_f32(&w,g->vel.x);net_write_f32(&w,g->vel.y);net_write_f32(&w,g->vel.z);net_write_f32(&w,g->fuse);net_write_u8(&w,g->kind);net_write_u8(&w,g->adhesive);net_write_u8(&w,g->attached);net_write_u64(&w,g->support);net_write_f32(&w,g->attachment.x);net_write_f32(&w,g->attachment.y);net_write_f32(&w,g->attachment.z);}
        if(w.failed){TraceLog(LOG_ERROR,"Grenade snapshot overflow");continue;}
        if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    }
}
static void clarity_net_event(ClarityEffect e) {
    if(netTransport.role!=NET_ROLE_HOST)return;
    uint8_t packet[80];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_GAME_EVENT,0,netTransport.session_id,netServerTick);
    net_write_u32(&w,e.id);net_write_u8(&w,e.kind);net_write_u8(&w,e.type);net_write_u8(&w,e.attacker<0?255:e.attacker);
    net_write_f32(&w,e.pos.x);net_write_f32(&w,e.pos.y);net_write_f32(&w,e.pos.z);net_write_f32(&w,e.yaw);
    net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_CONTROL,packet,w.length,true);
}
static bool clarity_net_receive(NetMessageType type,NetReader *r,uint32_t tick) {
    (void)tick;
    if(type==NET_MSG_GAME_EVENT){ClarityEffect e={0};e.id=net_read_u32(r);e.kind=net_read_u8(r);e.type=net_read_u8(r);int owner=net_read_u8(r);e.attacker=owner==255?-1:owner;
        e.pos=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};e.yaw=net_read_f32(r);
        if(!r->failed && e.id>clarityLastEvent && (e.kind==1 || e.kind==2) && e.type<ENEMY_TYPE_COUNT && v_isfinite(e.pos) && isfinite(e.yaw) && (owner==255 || owner<MAX_COMBATANTS)){clarityLastEvent=e.id;clarity_effect(e,false);}return true;}
    if(type!=NET_MSG_GRENADES)return false;
    uint32_t sequence=net_read_u32(r);int part=net_read_u8(r);GrenadeState incoming[16]={0};
    if(part>1 || sequence<=clarityReceiveSequence || sequence<clarityPendingSequence)return true;
    for(int n=0;n<16;n++){GrenadeState *g=&incoming[n];g->id=net_read_u32(r);int owner=net_read_u8(r);g->owner=owner==255?-1:owner;g->pos=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};g->vel=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};g->fuse=net_read_f32(r);g->kind=net_read_u8(r);g->adhesive=net_read_u8(r)!=0;g->attached=net_read_u8(r)!=0;g->support=net_read_u64(r);g->attachment=(Vector3){net_read_f32(r),net_read_f32(r),net_read_f32(r)};
        if(r->failed || g->kind>=GRENADE_KIND_COUNT || !v_isfinite(g->attachment) || !v_isfinite(g->pos) || !v_isfinite(g->vel) || !isfinite(g->fuse) || (g->fuse<0 && g->fuse!=-1) || g->fuse>1.5f || (g->id && (g->owner< -1 || g->owner>=MAX_PLAYERS)))return true;}
    if(sequence!=clarityPendingSequence){clarityPendingSequence=sequence;clarityReceiveMask=0;}
    memcpy(&clarityPending[part*16],incoming,sizeof(incoming));clarityReceiveMask|=1u<<part;
    if(clarityReceiveMask==3){for(int n=0;n<CLARITY_GRENADES;n++){GrenadeState old=grenades[n];grenades[n]=clarityPending[n];grenadeVisualTarget[n]=grenades[n].pos;if(old.id && old.id==grenades[n].id)grenades[n].pos=old.pos;}clarityReceiveSequence=sequence;}
    return true;
}
