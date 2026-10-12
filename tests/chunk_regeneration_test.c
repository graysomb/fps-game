#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void fresh(void) {
    netTransport.role=NET_ROLE_OFFLINE;activePlayers=0;clear_world_voxels();meshDirty=false;
}
static void advance(int count) {for(int i=0;i<count;i++){recycleFrameCounter++;regeneration_step(1.0f/60);}}
static void queue_cell(int x,int y,int z) {
    int i=add_static_voxel_at_grid(x,y,z,GRAY,0);assert(i>=0);
    assert(recycle_queue_push(&voxels[i]));remove_voxel_index(i);
}
static void due(void) {recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;advance(1);}
static void preview_packet(uint32_t serial,int total,int offset,int count,int first,bool truncate) {
    uint8_t data[1024];NetWriter w;net_writer_init(&w,data,sizeof(data));
    net_write_u32(&w,serial);net_write_u8(&w,1);
    net_write_i32(&w,-1);net_write_i32(&w,0);net_write_i32(&w,0);net_write_f32(&w,.25f);
    net_write_u16(&w,total);net_write_u16(&w,offset);net_write_u16(&w,count);
    for(int i=0;i<count;i++){net_write_u16(&w,first+i);net_write_u8(&w,100);net_write_u8(&w,170);net_write_u8(&w,200);}
    NetReader r;net_reader_init(&r,data,w.length-(truncate?1:0));assert(regeneration_net_receive(NET_MSG_REGENERATION,&r));
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(800,600,"Chunk regeneration regression");
    if(!DirectoryExists("artifacts"))assert(MakeDirectory("artifacts")==0);
    fresh();queue_cell(1,1,1);queue_cell(2,1,1);due();assert(regeneration.active && regeneration.count==2 && voxel_count==0);
    uint64_t generation=staticHashGeneration;advance(REGEN_PREVIEW_FRAMES-1);assert(voxel_count==0);
    advance(1);assert(!regeneration.active && voxel_count==2 && recycleQueueCount==0 && staticHashGeneration==generation+1);
    /* Later deadlines in the same chunk do not restore early. */
    fresh();queue_cell(1,1,1);recycleFrameCounter+=100;queue_cell(2,1,1);
    recycleFrameCounter=RECYCLE_STATIC_RESTORE_DELAY;advance(1);assert(regeneration.count==1);
    advance(REGEN_PREVIEW_FRAMES);assert(table_get_static_only(1,1,1)>=0 && table_get_static_only(2,1,1)<0 && recycleQueueCount==1);
    advance(120);assert(recycleQueueCount==0);
    /* Bounds crossing zero and a chunk edge remain queued until both portions restore. */
    fresh();queue_cell(-1,1,0);Voxel *snapshot=&recycleQueue[recycleQueueHead];snapshot->orig_max_gx=1;
    due();advance(REGEN_PREVIEW_FRAMES);assert(recycleQueueCount==1);
    advance(REGEN_PREVIEW_FRAMES+2);assert(recycleQueueCount==0);
    for(int x=-1;x<=1;x++)assert(table_get_static_only(x,1,0)>=0);
    /* An occupied older chunk does not starve a free chunk. */
    fresh();queue_cell(0,0,0);queue_cell(40,0,0);activePlayers=1;
    players[0].respawn_timer=0;players[0].pos=(Vector3){.25f,BASE_EYE_HEIGHT,.25f};
    due();assert(regeneration.active && regeneration.cx==1);advance(REGEN_PREVIEW_FRAMES);
    assert(table_get_static_only(40,0,0)>=0 && table_get_static_only(0,0,0)<0);
    players[0].pos.x=10;advance(REGEN_PREVIEW_FRAMES+2);assert(recycleQueueCount==0);
    /* Entering the preview cancels commit; later work can retry. */
    fresh();queue_cell(0,0,0);due();activePlayers=1;players[0].respawn_timer=0;players[0].pos=(Vector3){.25f,BASE_EYE_HEIGHT,.25f};
    advance(REGEN_PREVIEW_FRAMES);assert(!regeneration.active && voxel_count==0 && recycleQueueCount==1);
    /* Construction and queue cancellation invalidate prepared data. */
    fresh();queue_cell(1,1,1);queue_cell(2,1,1);due();add_static_voxel_at_grid(1,1,1,RED,0);
    advance(REGEN_PREVIEW_FRAMES);assert(!regeneration.active && voxel_count==1 && voxels[0].color.r==RED.r);
    fresh();queue_cell(1,1,1);due();recycle_queue_remove_region(1,1,1,1,1,1);advance(REGEN_PREVIEW_FRAMES);assert(voxel_count==0 && !regeneration.active);
    /* Capacity failure cannot publish part of a chunk. */
    fresh();queue_cell(1,1,1);due();int pool=particle_pool_count,freeCount=free_particle_count;
    particle_pool_count=MAX_PARTICLES;free_particle_count=0;advance(REGEN_PREVIEW_FRAMES);
    assert(voxel_count==0 && recycleQueueCount==1);particle_pool_count=pool;free_particle_count=freeCount;
    /* Dynamic debris blocks exact candidate cells, not the whole chunk AABB. */
    fresh();queue_cell(1,1,1);addVoxel(.75f,.75f,.75f,false,true,BLUE,0);due();assert(!regeneration.active);
    remove_voxel_index(0);advance(16);assert(regeneration.active);
    /* Preview captures, including four viewport sizes, share the same pose. */
    fresh();for(int x=2;x<12;x++)for(int y=0;y<8;y++)queue_cell(x,y,3);
    due();advance(18);Camera3D cam={{10,7,12},{3,2,1.5},{0,1,0},55,CAMERA_PERSPECTIVE};
    RenderTexture2D full=LoadRenderTexture(800,600);
    BeginTextureMode(full);ClearBackground((Color){20,25,35,255});BeginMode3D(cam);DrawGrid(24,1);regeneration_draw();EndMode3D();EndTextureMode();
    Image capture=LoadImageFromTexture(full.texture);ImageFlipVertical(&capture);ExportImage(capture,"artifacts/regeneration-preview.png");UnloadImage(capture);assert(regenerationMeshCount>0);
    RenderTexture2D small=LoadRenderTexture(400,300);
    BeginTextureMode(small);ClearBackground((Color){20,25,35,255});BeginMode3D(cam);DrawGrid(24,1);regeneration_draw();EndMode3D();EndTextureMode();
    BeginTextureMode(full);ClearBackground(BLACK);for(int q=0;q<4;q++)DrawTextureRec(small.texture,(Rectangle){0,0,400,-300},(Vector2){(q%2)*400,(q/2)*300},WHITE);EndTextureMode();
    capture=LoadImageFromTexture(full.texture);ImageFlipVertical(&capture);ExportImage(capture,"artifacts/regeneration-four-view.png");UnloadImage(capture);UnloadRenderTexture(small);
    cam.position=(Vector3){3,2,10};cam.target=(Vector3){3,2,1.5f};
    BeginTextureMode(full);ClearBackground(BLACK);BeginMode3D(cam);DrawCube((Vector3){3,2,5},2,6,1,RED);regeneration_draw();EndMode3D();EndTextureMode();
    capture=LoadImageFromTexture(full.texture);Color center=GetImageColor(capture,400,300);
    assert(center.r>200 && center.g<60 && center.b<70); /* Preview respects opaque depth. */
    UnloadImage(capture);UnloadRenderTexture(full);
    fresh();assert(regenerationMeshCount==0 && !regeneration.active && !recycleQueueCount);
    /* Reliable multipart previews publish only after the complete packet set. */
    netTransport.role=NET_ROLE_CLIENT;regeneration.serial=0;
    preview_packet(1,2,0,1,5,false);assert(!regeneration.active && regeneration.count==1);
    preview_packet(1,2,1,1,6,true);assert(!regeneration.active && regeneration.count==1);
    preview_packet(1,2,1,1,6,false);assert(regeneration.active && regeneration.count==2);
    uint8_t packet[8];NetWriter w;net_writer_init(&w,packet,sizeof(packet));net_write_u32(&w,1);net_write_u8(&w,0);
    NetReader r;net_reader_init(&r,packet,w.length);regeneration_net_receive(NET_MSG_REGENERATION,&r);assert(!regeneration.active);
    preview_packet(2,1,0,1,REGEN_CELL_COUNT,false);assert(!regeneration.active);
    regeneration_shutdown();clear_world_voxels();CloseWindow();puts("Chunk regeneration regression passed");return 0;
}
