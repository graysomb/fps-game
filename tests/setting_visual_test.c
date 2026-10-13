#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static int ms_compare(const void *a,const void *b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static uint64_t layout_hash(void) {
    uint64_t h=1469598103934665603ull;
    for(int i=0;i<voxel_count;i++) {
        Voxel *v=&voxels[i];int data[]={v->min_gx,v->max_gx,v->min_gy,v->max_gy,v->min_gz,v->max_gz,v->simulate,v->type};
        for(size_t j=0;j<sizeof(data);j++){h^=((unsigned char*)data)[j];h*=1099511628211ull;}
    }
    return h;
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,640,"Setting visual regression");
    init_pbd_thread_pool();assert(initialize_physics_backend());
    for(int map=1;map<=2;map++) {
        gameMode=map==1?GAME_MODE_FIREFIGHT:GAME_MODE_DEATHMATCH;firefightFoundry=true;
        currentWorldType=WORLD_TYPE_CROSSCUT;multiplayerActivePlayers=1;multiplayerPlayerInput[0]=INPUT_TYPE_KEYBOARD;
        ResetGame();activePlayers=1;
        fprintf(stderr,"map=%d voxels=%d layout=%016llx\n",map,voxel_count,(unsigned long long)layout_hash());
        assert(layout_hash()==(map==1?0xc33455752e176ef3ull:0x2cddfd9bb3219393ull));
        assert(settingMap==map && settingDetailCount>0);
        uint64_t unchanged=layout_hash();int particles=particle_pool_count,details=settingDetailCount;
        setting_receive_layout(map);assert(settingDetailCount==details && layout_hash()==unchanged && particle_pool_count==particles);
        RenderTexture2D target=LoadRenderTexture(960,640);
        for(int view=0;view<5;view++) {
            Camera3D cam={{37,32,37},{0,0,0},{0,1,0},65,CAMERA_PERSPECTIVE};
            if(view){cam.position=crosscut_rotate((Vector3){12,1.1f,21},view-1);cam.target=crosscut_rotate((Vector3){0,2,0},view-1);}
            voxel_render_frame++;BeginTextureMode(target);ClearBackground((Color){33,40,43,255});BeginMode3D(cam);draw_world_surfaces();DrawVoxels(cam);EndMode3D();EndTextureMode();
            Image img=LoadImageFromTexture(target.texture);ImageFlipVertical(&img);
            ExportImage(img,TextFormat("artifacts/setting-%s-%d-%d.png",getenv("SETTING_BASELINE")?"before":"after",map,view));UnloadImage(img);
        }
        UnloadRenderTexture(target);
        // Same scene, targets and camera with and without the decorative pass.
        for(int damaged=0;damaged<2;damaged++) {
        if(damaged) {
            if(map==1)assert(activate_static_voxels_near_region(39,46,0,9,-10,10,0));
            else assert(activate_static_voxels_near_region(21,35,0,9,6,16,0));
            simulate_voxel_pbd_steps(1.0f/120,2);
        }
        for(int views=1;views<=4;views+=3) {
            RenderTexture2D rt=LoadRenderTexture(views==1?960:480,views==1?640:320);
            for(int enabled=0;enabled<2;enabled++) {
                settingVisualsEnabled=enabled;double times[100];
                for(int frame=0;frame<110;frame++) {
                    double start=GetTime();voxel_render_frame++;
                    for(int v=0;v<views;v++) {
                        Camera3D cam={crosscut_rotate((Vector3){12,1.1f,21},v),{0,2,0},{0,1,0},65,CAMERA_PERSPECTIVE};
                        BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(cam);draw_world_surfaces();DrawVoxels(cam);
                        if(map==1)for(int type=0;type<ENEMY_TYPE_COUNT;type++) {
                            Vector3 p={(type%5-2)*2.0f,1,(type/5)*3.0f};
                            if(enabled)robot_draw(type,p,0,.5f,SKYBLUE,1.25f,1,0,false);
                            else DrawCube(p,.5f,.5f,.5f,GRAY);
                        }
                        EndMode3D();EndTextureMode();
                    }
                    BeginDrawing();ClearBackground(BLACK);EndDrawing();
                    if(frame>=10)times[frame-10]=(GetTime()-start)*1000;
                }
                qsort(times,100,sizeof(double),ms_compare);
                fprintf(stderr,"render map=%d damaged=%d views=%d setting=%d p50=%.3f p95=%.3f ms\n",map,damaged,views,enabled,times[50],times[95]);
            }
            UnloadRenderTexture(rt);
        }
        }
        settingVisualsEnabled=true;setting_update();
        int chosen=-1;for(int n=0;n<settingDetailCount;n++)if(settingDetails[n].visible){chosen=n;break;}assert(chosen>=0);
        SettingDetail *d=&settingDetails[chosen];int *cell=d->cells[2];
        int victim=table_get_static_only(cell[0],cell[1],cell[2]);assert(victim>=0);
        assert(recycle_queue_push(&voxels[victim]));remove_voxel_index(victim);setting_update();assert(!d->visible);
        recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;
        // Destruction can queue other chunks ahead of this panel's backing cell.
        for(int f=0;f<2000 && !d->visible;f++){recycle_dead_voxels(1.0f/60);setting_update();}
        if(!d->visible)fprintf(stderr,"panel %d map %d did not restore; pending=%d\n",chosen,map,recycleQueueCount);
        assert(d->visible);
        unchanged=layout_hash();setting_receive_layout(map);assert(layout_hash()==unchanged);setting_update();assert(settingDetails[chosen].visible);
        useCustomMap=true;assert(!setting_active());useCustomMap=false;
        enemyLab=true;assert(!setting_active());enemyLab=false;
        clear_world_voxels();assert(settingMap==0 && settingDetailCount==0);
    }
    RenderTexture2D lineup=LoadRenderTexture(1200,600);
    for(int type=0;type<ENEMY_TYPE_COUNT;type++) {
        clarity_menu_render(type,210,1.25f);
        BeginTextureMode(lineup);
        if(type==0)ClearBackground((Color){15,25,30,255});
        DrawTexturePro(clarityShowcase.texture,(Rectangle){0,0,720,-720},(Rectangle){(type%5)*240,(type/5)*300,240,260},(Vector2){0},0,WHITE);
        DrawText(robot_roles[type],(type%5)*240+8,(type/5)*300+268,12,RAYWHITE);
        EndTextureMode();
    }
    Image image=LoadImageFromTexture(lineup.texture);ImageFlipVertical(&image);ExportImage(image,"artifacts/setting-robot-lineup.png");UnloadImage(image);UnloadRenderTexture(lineup);
    clarity_menu_shutdown();
    puts("Setting geometry, decoration lifecycle, robot previews and rendering checks passed");
    shutdown_pbd_thread_pool();gpu_physics_shutdown();shutdown_world_visuals();CloseWindow();return 0;
}
