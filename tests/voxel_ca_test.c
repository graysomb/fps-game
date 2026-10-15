#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static int compare_ms(const void *a,const void *b) {
    double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);
}
static void check_ca_rule(void) {
    assert(voxel_ca_init());
    Image initial=LoadImageFromTexture(voxel_ca_texture());Color *state=LoadImageColors(initial);UnloadImage(initial);
    Color *next=MemAlloc(256*256*sizeof(Color));
    for(unsigned gen=0;gen<12;gen++) {
        for(int y=0;y<256;y++)for(int x=0;x<256;x++) {
            int baseX=x/128*128,baseY=y/128*128;
            int current=state[y*256+x].r/255,previous=state[y*256+x].g/255,count=0;
            for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++)if(dx || dy) {
                int xx=baseX+(x%128+dx+128)%128,yy=baseY+(y%128+dy+128)%128;
                count+=state[yy*256+xx].r/255;
            }
            // Independent expression of the source rule, rather than GPU table lookup.
            bool alive=previous ? (current ? count==1 || count==2 : count==6 || count==8)
                                : !current && (count==2 || count==7);
            int trail=alive!=current?255:state[y*256+x].b>16?state[y*256+x].b-16:0;
            next[y*256+x]=(Color){alive*255,current*255,trail,255};
        }
        voxel_ca_step();Image image=LoadImageFromTexture(voxel_ca_texture());Color *gpu=LoadImageColors(image);
        for(int i=0;i<256*256;i++) {
            if(gpu[i].r!=next[i].r || gpu[i].g!=next[i].g)fprintf(stderr,"CA mismatch generation=%u cell=%d GPU=%u/%u CPU=%u/%u\n",gen,i,gpu[i].r,gpu[i].g,next[i].r,next[i].g);
            assert(gpu[i].r==next[i].r && gpu[i].g==next[i].g && gpu[i].b==next[i].b);
        }
        memcpy(state,next,256*256*sizeof(Color));UnloadImageColors(gpu);UnloadImage(image);
    }
    MemFree(next);UnloadImageColors(state);
    unsigned before=voxelCaGeneration;voxel_ca_update(voxelCaLastTime+10);
    assert(voxelCaGeneration==before+2 && voxelCaAccumulator<VOXEL_CA_STEP);
    voxelCaEnabled=false;before=voxelCaGeneration;voxel_ca_update(voxelCaLastTime+1);assert(voxelCaGeneration==before);voxelCaEnabled=true;
    puts("CA CPU/GPU rules, independent tile wrapping, second-order memory and trails and update cap passed");
}
static void check_sustained_activity(void) {
    unsigned long activity[4]={0};
    for(int step=0;step<1200;step++) {
        voxel_ca_step();
        if(step<1000)continue;
        Image image=LoadImageFromTexture(voxel_ca_texture());Color *cells=LoadImageColors(image);
        for(int y=0;y<256;y++)for(int x=0;x<256;x++) {
            assert((cells[y*256+x].r==0 || cells[y*256+x].r==255) && (cells[y*256+x].g==0 || cells[y*256+x].g==255));
            activity[(y/128)*2+x/128]+=cells[y*256+x].r!=cells[y*256+x].g;
        }
        UnloadImageColors(cells);UnloadImage(image);
    }
    for(int tile=0;tile<4;tile++) {
        double rate=(double)activity[tile]/(200*128*128);
        fprintf(stderr,"Tile %d last-ten-second activity: %.2f%% cells/step\n",tile,rate*100);
        assert(rate>.01);
    }
}
static void capture_animation(void) {
    clear_all_chunks();clear_world_voxels();
    for(int x=-1;x<=1;x++)for(int y=0;y<2;y++)add_static_voxel_at_grid(x,y,0,(Color){130,167,177,255},0);
    Camera3D cam={{2,1.4f,2.7f},{.25f,.5f,.25f},{0,1,0},1.8f,CAMERA_ORTHOGRAPHIC};
    RenderTexture2D rt=LoadRenderTexture(640,480);
    for(int f=0;f<192;f++) {
        voxel_ca_update(voxelCaLastTime+1.0/24);voxel_render_frame++;
        BeginTextureMode(rt);ClearBackground((Color){30,35,40,255});BeginMode3D(cam);DrawVoxels(cam);EndMode3D();EndTextureMode();
        Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);
        ExportImage(image,TextFormat("artifacts/voxel-second-order-%03d.png",f));UnloadImage(image);
    }
    UnloadRenderTexture(rt);
}
static void profile_maps(void) {
    for(int map=0;map<2;map++) {
        gameMode=map?GAME_MODE_DEATHMATCH:GAME_MODE_FIREFIGHT;firefightFoundry=true;currentWorldType=WORLD_TYPE_CROSSCUT;multiplayerActivePlayers=1;ResetGame();
        for(int views=1;views<=4;views+=3) {
            RenderTexture2D rt=LoadRenderTexture(views==1?1280:640,views==1?720:360);
            for(int enabled=0;enabled<2;enabled++) {
                voxelCaEnabled=enabled;double times[100],updates[100];int updateCount=0;
                for(int f=0;f<110;f++) {
                    double start=GetTime();unsigned before=voxelCaGeneration;
                    voxel_ca_update(voxelCaLastTime+1.0/60);voxel_render_frame++;
                    for(int v=0;v<views;v++) {
                        Camera3D cam={crosscut_rotate((Vector3){12,1.2f,21},v),{0,2,0},{0,1,0},65,CAMERA_PERSPECTIVE};
                        BeginTextureMode(rt);ClearBackground(SKYBLUE);BeginMode3D(cam);draw_world_surfaces();DrawVoxels(cam);EndMode3D();EndTextureMode();
                    }
                    BeginDrawing();ClearBackground(BLACK);EndDrawing();
                    if(f>=10){double ms=(GetTime()-start)*1000;times[f-10]=ms;if(voxelCaGeneration!=before)updates[updateCount++]=ms;}
                }
                qsort(times,100,sizeof(double),compare_ms);qsort(updates,updateCount,sizeof(double),compare_ms);
                fprintf(stderr,"map=%d views=%d CA=%d p50=%.3f p95=%.3f updateMax=%.3f ms\n",map,views,enabled,times[50],times[95],updateCount?updates[updateCount-1]:0);
                Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);
                ExportImage(image,TextFormat("artifacts/voxel-second-order-map%d-views%d-enabled%d.png",map,views,enabled));UnloadImage(image);
            }
            UnloadRenderTexture(rt);
        }
    }
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Voxel CA validation");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    check_ca_rule();check_sustained_activity();capture_animation();profile_maps();
    voxel_ca_shutdown();assert(!voxelCaReady && !voxelCaTargets[0].id && !voxelCaTargets[1].id);
    assert(voxel_ca_init());voxel_ca_shutdown();
    clear_all_chunks();clear_world_voxels();shutdown_world_visuals();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
