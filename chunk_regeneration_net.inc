/* Preview packets share the reliable world channel with authoritative commits. */
static void regeneration_send_packet(int peer,const uint8_t *packet,size_t size) {
    if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,packet,size,true);
    else net_transport_broadcast(&netTransport,FPS_NET_CHANNEL_WORLD,packet,size,true);
}
static void regeneration_net_send(int peer,bool active) {
    if(netTransport.role!=NET_ROLE_HOST || (active && !regeneration.active))return;
    int offset=0,cursor=0;
    do {
        uint8_t packet[FPS_NET_MAX_PACKET];NetWriter w;net_writer_init(&w,packet,sizeof(packet));
        net_write_header(&w,NET_MSG_REGENERATION,0,netTransport.session_id,netServerTick);
        net_write_u32(&w,regeneration.serial);net_write_u8(&w,active?1:0);
        if(active) {
            net_write_i32(&w,regeneration.cx);net_write_i32(&w,regeneration.cy);net_write_i32(&w,regeneration.cz);
            net_write_f32(&w,Clamp(regeneration.elapsed/.6f,0,1));
            net_write_u16(&w,regeneration.count);net_write_u16(&w,offset);
            int count=regen_min(128,regeneration.count-offset);net_write_u16(&w,count);
            for(int n=0;n<count;n++) {
                while(cursor<REGEN_CELL_COUNT && !regeneration.cells[cursor].present)cursor++;
                Color c=regeneration.cells[cursor].color;
                net_write_u16(&w,cursor++);net_write_u8(&w,c.r);net_write_u8(&w,c.g);net_write_u8(&w,c.b);
            }
            offset+=count;
        }
        regeneration_send_packet(peer,packet,w.length);
    }while(active && offset<regeneration.count);
}
static bool regeneration_net_receive(NetMessageType type,NetReader *r) {
    if(type!=NET_MSG_REGENERATION)return false;
    if(netTransport.role!=NET_ROLE_CLIENT)return true;
    uint32_t serial=net_read_u32(r);int active=net_read_u8(r);
    if(r->failed || active>1 || serial<regeneration.serial)return true;
    if(!active){regeneration_reset();regeneration.serial=serial;return true;}
    int cx=net_read_i32(r),cy=net_read_i32(r),cz=net_read_i32(r);
    float progress=net_read_f32(r);
    int total=net_read_u16(r),offset=net_read_u16(r),count=net_read_u16(r);
    if(r->failed || cx< -64 || cx>64 || cy< -64 || cy>64 || cz< -64 || cz>64 ||
       !isfinite(progress) || progress<0 || progress>1 || total<1 || total>REGEN_CELL_COUNT ||
       count<1 || count>128 || offset+count>total)return true;
    int indices[128];Color colors[128];
    for(int n=0;n<count;n++) {
        indices[n]=net_read_u16(r);
        colors[n]=(Color){net_read_u8(r),net_read_u8(r),net_read_u8(r),255};
        if(indices[n]>=REGEN_CELL_COUNT || (n && indices[n]<=indices[n-1]))return true;
    }
    if(r->failed)return true;
    if(offset==0) {
        regeneration_reset();regeneration.serial=serial;
        regeneration.cx=cx;regeneration.cy=cy;regeneration.cz=cz;
        regeneration.receiveTotal=total;regeneration.remoteStart=GetTime()-progress*.6;
    }
    if(serial!=regeneration.serial || offset!=regeneration.count || total!=regeneration.receiveTotal ||
       cx!=regeneration.cx || cy!=regeneration.cy || cz!=regeneration.cz)return true;
    for(int n=0;n<count;n++)if(regeneration.cells[indices[n]].present)return true;
    for(int n=0;n<count;n++)regeneration.cells[indices[n]]=(RegenerationCell){.present=true,.color=colors[n]};
    regeneration.count+=count;
    if(regeneration.count==total)regeneration.active=regeneration.meshDirty=true;
    return true;
}
