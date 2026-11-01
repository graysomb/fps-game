#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
#ifndef GRAPHICS_API_OPENGL_43
#include <GL/gl.h>
#endif
bool ffTestReviveHeld[4];
#include "../victory_town.inc"
static void town_render(RenderTexture2D rt,WeaponBloom *b,TownSample s,bool bloom){
 BeginTextureMode(rt);ClearBackground((Color){19,31,43,255});
 /* Existing generated sky and meteors, with dusk horizon for this preview. */
 town_sky(s.camera,s.time);
 BeginMode3D(s.camera);victory_town_draw(s,false);EndMode3D();EndTextureMode();
 if(bloom){weapon_bloom_begin(b,rt);BeginMode3D(s.camera);victory_town_draw(s,true);EndMode3D();weapon_bloom_composite(b,rt);EndTextureMode();}
}
static int compare_times(const void *a,const void *b){double x=*(double*)a,y=*(double*)b;return(x>y)-(x<y);}
int main(int argc,char **argv){
 bool qa=false,individual=false,capture=false,bench=false,fallback=false,bloom=true;int frames=640,width=1280,height=720;
 for(int i=1;i<argc;i++){if(!strcmp(argv[i],"--capture"))capture=true;if(!strcmp(argv[i],"--scenes")){capture=true;individual=true;}if(!strcmp(argv[i],"--benchmark"))bench=true;if(!strcmp(argv[i],"--fallback"))fallback=true;if(!strcmp(argv[i],"--no-bloom"))bloom=false;if(!strcmp(argv[i],"--qa")){capture=true;qa=true;frames=12;}if(!strcmp(argv[i],"--stills")){capture=true;frames=4;}}
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);if(capture||bench)SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(width,height,"Reclaimed Foundry - standalone victory preview");
 SetRandomSeed(73021);srand(73021);init_pbd_thread_pool();assert(initialize_physics_backend());character_renderer_init();weapon_renderer_init();duck_renderer_init(NULL);camp_init();
 clear_all_chunks();clear_world_voxels();activePlayers=0;worldPeriodic=false;worldFloorHalfSize=40;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;foundry_selection_seed(73);buildFoundryHoldout();init_static_hash();victory_town_init();victory_town_validate();clear_all_chunks();clear_world_voxels();
 if(fallback){Shader sh={.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};for(int i=0;i<T_COUNT;i++)town.models[i].materials[0].shader=sh;}
 RenderTexture2D a=LoadRenderTexture(width,height),b=LoadRenderTexture(width,height),out=LoadRenderTexture(width,height);WeaponBloom glow={0};assert(weapon_bloom_resize(&glow,width,height));
 double times[640]={0};float time=0;bool pause=false;int fixed=-1,frame=0;const char *names[4]={"Supper on the platform","Homes in the Foundry","Garden and workshop","A place to call home"};
 while(!WindowShouldClose()){
 if(IsKeyPressed(KEY_SPACE))pause=!pause;if(IsKeyPressed(KEY_B))bloom=!bloom;if(IsKeyPressed(KEY_ZERO))fixed=-1;for(int k=0;k<4;k++)if(IsKeyPressed(KEY_ONE+k))fixed=k;
 if(IsKeyDown(KEY_RIGHT))time+=GetFrameTime()*3;if(IsKeyDown(KEY_LEFT))time=fmaxf(0,time-GetFrameTime()*3);
 if(capture||bench)time=qa?(frame/3)*8+(frame%3==0?0:frame%3==1?4:7.95f):frames==4?frame*8+3:frame/20.f;else if(!pause)time+=GetFrameTime();
 int scene=fixed>=0?fixed:((int)(time/8)%4);float local=fmodf(time,8),blend=fixed<0&&!individual&&!qa?fmaxf(0,local-7):0;TownSample s=victory_town_sample(scene,local);
 assert(voxel_count==0);assert(isfinite(s.camera.position.x)&&isfinite(s.camera.position.y));
 double start=GetTime();town_render(a,&glow,s,bloom);if(blend>0)town_render(b,&glow,victory_town_sample((scene+1)%4,local-8),bloom);
 BeginTextureMode(out);ClearBackground(BLACK);DrawTextureRec(a.texture,(Rectangle){0,0,width,-height},(Vector2){0},WHITE);if(blend>0)DrawTextureRec(b.texture,(Rectangle){0,0,width,-height},(Vector2){0},Fade(WHITE,blend));EndTextureMode();
 glFinish();if(frame<640)times[frame]=(GetTime()-start)*1000;
 if(capture){Image im=LoadImageFromTexture(out.texture);ImageFlipVertical(&im);char path[256];snprintf(path,sizeof(path),qa?"artifacts/victory-town/qa/%02d.png":frames==4?"artifacts/victory-town/scene-%d.png":individual?"artifacts/victory-town/scenes/%04d.png":"artifacts/victory-town/frames/%04d.png",frame);assert(ExportImage(im,path));UnloadImage(im);}
 BeginDrawing();ClearBackground(BLACK);DrawTextureRec(out.texture,(Rectangle){0,0,width,-height},(Vector2){0},WHITE);if(!capture&&!bench){DrawText(names[scene],26,24,24,RAYWHITE);DrawText("1-4 scenes | 0 montage | Space pause | arrows scrub | B bloom",26,height-34,18,RAYWHITE);}EndDrawing();
 frame++;if((capture&&frame==frames)||(bench&&frame==640))break;
 }
 if(bench){qsort(times+40,600,sizeof(double),compare_times);printf("1280x720 render + GPU completion, bloom=%d: median %.3f ms; p95 %.3f ms (600 frames, 40 warmup)\n",bloom,times[340],times[610]);}
 weapon_bloom_unload(&glow);UnloadRenderTexture(a);UnloadRenderTexture(b);UnloadRenderTexture(out);victory_town_shutdown();camp_shutdown();duck_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();shutdown_world_visuals();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
