static void eq_write_vec(NetWriter *w,Vector3 v){net_write_f32(w,v.x);net_write_f32(w,v.y);net_write_f32(w,v.z);}
static Vector3 eq_read_vec(NetReader *r){Vector3 v={net_read_f32(r),net_read_f32(r),net_read_f32(r)};return v;}
static void eq_snapshot(int peer){
    if(netTransport.role!=NET_ROLE_HOST)return;
    for(int part=0;part<4;part++){
        uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_header(&w,NET_MSG_EQUIPMENT,0,netTransport.session_id,netServerTick);net_write_u8(&w,part);
        if(part==0)for(int i=0;i<MAX_PLAYERS;i++){Equipment *e=&equipment[i];net_write_u8(&w,e->gun[0]);net_write_u8(&w,e->gun[1]);net_write_u8(&w,e->selected);net_write_u8(&w,e->upgrades);net_write_u8(&w,(e->held?1:0)|(e->buildHeld?2:0));net_write_f32(&w,e->switchTime);net_write_f32(&w,e->charge);net_write_f32(&w,e->buildTime);for(int g=0;g<EQ_GUN_COUNT;g++)net_write_f32(&w,e->cooldown[g]);eq_write_vec(&w,e->beamEnd);net_write_u8(&w,aimStates[i].active);net_write_u32(&w,netLastProcessedInput[i]);}
        if(part==1)for(int n=0;n<32;n++){EquipmentRocket *r=&eqRockets[n];net_write_u8(&w,r->active);net_write_u8(&w,r->owner);eq_write_vec(&w,r->pos);eq_write_vec(&w,r->vel);net_write_f32(&w,r->life);}
        if(part==2)for(int n=0;n<32;n++){EquipmentTrace *t=&eqTraces[n];eq_write_vec(&w,t->start);eq_write_vec(&w,t->end);net_write_f32(&w,t->life);net_write_u8(&w,t->gun);net_write_u8(&w,t->owner);}
        if(part==3)for(int i=0;i<MAX_PLAYERS;i++){
            SwordWave *s=&swordWaves[i];net_write_u8(&w,s->active);net_write_u8(&w,s->emitted);net_write_u32(&w,s->sequence);
            net_write_u8(&w,s->direction);net_write_u8(&w,swordNextDirection[i]);net_write_f32(&w,s->age);
            eq_write_vec(&w,s->origin);eq_write_vec(&w,s->forward);
            for(int n=0;n<SWORD_LANES;n++)net_write_f32(&w,s->stop[n]);
        }
        if(w.failed){TraceLog(LOG_ERROR,"Equipment packet overflow");continue;}
        if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,w.length,true);else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_SNAPSHOT,packet,w.length,false);
    }
}

static bool eq_receive(NetMessageType type,NetReader *r,uint32_t tick){
    if(type!=NET_MSG_EQUIPMENT)return false;if(netTransport.role!=NET_ROLE_CLIENT)return true;
    int part=net_read_u8(r);if(part>3 || tick<eqPartTick[part])return true;
    if(part==0){Equipment temp[MAX_PLAYERS]={0};uint8_t aimActive[MAX_PLAYERS]={0};uint32_t aimAck[MAX_PLAYERS]={0};for(int i=0;i<MAX_PLAYERS;i++){
        Equipment *e=&temp[i];e->gun[0]=net_read_u8(r);e->gun[1]=net_read_u8(r);e->selected=net_read_u8(r);e->upgrades=net_read_u8(r);int held=net_read_u8(r);e->held=(held&1)!=0;e->buildHeld=(held&2)!=0;e->switchTime=net_read_f32(r);e->charge=net_read_f32(r);e->buildTime=net_read_f32(r);for(int g=0;g<EQ_GUN_COUNT;g++){e->cooldown[g]=net_read_f32(r);if(!isfinite(e->cooldown[g])||e->cooldown[g]<0||e->cooldown[g]>3)return true;}e->beamEnd=eq_read_vec(r);e->initialized=true;aimActive[i]=net_read_u8(r);aimAck[i]=net_read_u32(r);
        if(r->failed || aimActive[i]>1 || e->selected>1 || held>3 || e->upgrades>7 || e->gun[0]>=EQ_GUN_COUNT || (e->gun[1]>=EQ_GUN_COUNT&&e->gun[1]!=EQ_EMPTY) || e->gun[e->selected]==EQ_EMPTY || !v_isfinite(e->beamEnd) || !isfinite(e->charge)||e->charge<0||e->charge>.66f || !isfinite(e->switchTime)||e->switchTime<0||e->switchTime>.26f || !isfinite(e->buildTime)||e->buildTime<0||e->buildTime>.61f)return true;
    }for(int i=0;i<MAX_PLAYERS;i++){temp[i].depleted=equipment[i].depleted;temp[i].depletionPulse=equipment[i].depletionPulse;temp[i].wasEmpty=equipment[i].wasEmpty;if(net_player_is_local(i))temp[i].buildHeld=equipment[i].buildHeld;}memcpy(equipment,temp,sizeof(temp));for(int i=0;i<MAX_PLAYERS;i++){players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;if(!aimStates[i].pending || aimAck[i]>=aimStates[i].pending){aimStates[i].active=aimActive[i]!=0 && !(net_player_is_local(i)&&(gameState!=GAME_STATE_PLAYING || !IsWindowFocused() || debugConsoleOpen));aimStates[i].pending=0;}}}
    if(part==1){EquipmentRocket temp[32]={0};for(int n=0;n<32;n++){EquipmentRocket *a=&temp[n];int active=net_read_u8(r);a->active=active!=0;a->owner=net_read_u8(r);a->pos=eq_read_vec(r);a->vel=eq_read_vec(r);a->life=net_read_f32(r);if(r->failed||active>1||a->owner>=MAX_PLAYERS||!v_isfinite(a->pos)||!v_isfinite(a->vel)||!isfinite(a->life)||a->life<-.1f||a->life>5)return true;}memcpy(eqRockets,temp,sizeof(temp));}
    if(part==2){EquipmentTrace temp[32]={0};for(int n=0;n<32;n++){EquipmentTrace *t=&temp[n];t->start=eq_read_vec(r);t->end=eq_read_vec(r);t->life=net_read_f32(r);t->gun=net_read_u8(r);t->owner=net_read_u8(r);if(r->failed||t->owner>=MAX_PLAYERS||t->gun>=EQ_GUN_COUNT||!v_isfinite(t->start)||!v_isfinite(t->end)||!isfinite(t->life)||t->life<0||t->life>.23f)return true;}memcpy(eqTraces,temp,sizeof(temp));}
    if(part==3){
        SwordWave temp[MAX_PLAYERS]={0};int next[MAX_PLAYERS]={0};
        for(int i=0;i<MAX_PLAYERS;i++){
            SwordWave *s=&temp[i];int active=net_read_u8(r),emitted=net_read_u8(r);s->active=active;s->emitted=emitted;
            s->sequence=net_read_u32(r);s->direction=net_read_u8(r);next[i]=net_read_u8(r);s->age=net_read_f32(r);
            s->origin=eq_read_vec(r);s->forward=eq_read_vec(r);
            for(int n=0;n<SWORD_LANES;n++){s->stop[n]=net_read_f32(r);if(!isfinite(s->stop[n]) || s->stop[n]<0 || s->stop[n]>SWORD_REACH)return true;}
            if(r->failed || active>1 || emitted>1 || s->direction>1 || next[i]>1 || !isfinite(s->age) || s->age<0 || s->age>SWORD_ANIM_SECONDS ||
                !v_isfinite(s->origin) || !v_isfinite(s->forward) || (emitted && fabsf(v_length(s->forward)-1)>.01f))return true;
        }
        memcpy(swordWaves,temp,sizeof(temp));
        for(int i=0;i<MAX_PLAYERS;i++){
            swordNextDirection[i]=next[i];swordDirection[i]=temp[i].direction;swordSequences[i]=temp[i].sequence;
            if(temp[i].active && players[i].respawn_timer<=0 && eq_gun(i)==EQ_BLADE){
                players[i].swordSwing=true;players[i].meleeSwingActive=true;players[i].meleeSwingStartTime=(float)GetTime()-temp[i].age;
            }
        }
    }
    eqPartTick[part]=tick;return true;
}
