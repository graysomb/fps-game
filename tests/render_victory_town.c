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
#include "../victory_town_credits.inc"
static void town_render(RenderTexture2D rt,WeaponBloom *b,TownSample s,bool bloom){
 BeginTextureMode(rt);ClearBackground((Color){19,31,43,255});
 /* Existing generated sky and meteors, with dusk horizon for this preview. */
 town_sky(s);
 BeginMode3D(s.camera);victory_town_draw(s,false);EndMode3D();EndTextureMode();
 if(bloom){weapon_bloom_begin(b,rt);BeginMode3D(s.camera);victory_town_draw(s,true);EndMode3D();weapon_bloom_composite(b,rt);EndTextureMode();}
}
static int compare_times(const void *a,const void *b){double x=*(double*)a,y=*(double*)b;return(x>y)-(x<y);}
int main(int argc,char **argv){
 bool creditStills=false,video=false,audioTest=false,qa=false,individual=false,capture=false,bench=false,fallback=false,bloom=true;int frames=(int)(town_duration()*20),width=1280,height=720;
 for(int i=1;i<argc;i++){if(!strcmp(argv[i],"--credits-stills")){creditStills=true;capture=true;frames=5;}if(!strcmp(argv[i],"--video"))video=true;if(!strcmp(argv[i],"--audio-test"))audioTest=true;if(!strcmp(argv[i],"--capture"))capture=true;if(!strcmp(argv[i],"--scenes")){capture=true;individual=true;}if(!strcmp(argv[i],"--benchmark"))bench=true;if(!strcmp(argv[i],"--fallback"))fallback=true;if(!strcmp(argv[i],"--no-bloom"))bloom=false;if(!strcmp(argv[i],"--qa")){capture=true;qa=true;frames=32;}if(!strcmp(argv[i],"--stills")){capture=true;frames=8;}}
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);if(capture||bench||video||audioTest)SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(width,height,"Reclaimed Foundry - standalone victory preview");
 if(video)frames=(int)ceilf(TOWN_MUSIC_SECONDS*20);
 Music soundtrack={0};bool audioReady=false;
 if((!capture&&!bench&&!video)||audioTest){
  InitAudioDevice();
  if(IsAudioDeviceReady()){soundtrack=LoadMusicStream(TOWN_MUSIC_PATH);audioReady=IsMusicValid(soundtrack);
   if(audioReady){soundtrack.looping=false;SetMusicVolume(soundtrack,.75f);if(audioTest)PlayMusicStream(soundtrack);}}
  if(!audioReady)fprintf(stderr,"Victory soundtrack unavailable; using the timed silent preview.\n");
 }
 if(audioTest){assert(audioReady);assert(fabsf(GetMusicTimeLength(soundtrack)-TOWN_MUSIC_SECONDS)<.1f);
  for(int i=0;i<30;i++){UpdateMusicStream(soundtrack);WaitTime(.01);}
  PauseMusicStream(soundtrack);SeekMusicStream(soundtrack,30);assert(!IsMusicStreamPlaying(soundtrack));
  ResumeMusicStream(soundtrack);assert(IsMusicStreamPlaying(soundtrack));
  SeekMusicStream(soundtrack,TOWN_MUSIC_SECONDS-.1f);UpdateMusicStream(soundtrack);
  StopMusicStream(soundtrack);UnloadMusicStream(soundtrack);CloseAudioDevice();CloseWindow();puts("Victory music duration, pause, seek, resume and cleanup passed");return 0;}
 FILE *movie=NULL;
 if(video){
  const char *ffmpeg=getenv("TOWN_FFMPEG");if(!ffmpeg)ffmpeg="ffmpeg";
  assert(!strchr(ffmpeg,'"'));char command[2048];
  snprintf(command,sizeof(command),"\"%s\" -y -loglevel error -f rawvideo -pixel_format rgb24 -video_size 1280x720 -framerate 20 -i - -i \"" TOWN_MUSIC_PATH "\" -map 0:v:0 -map 1:a:0 -c:v libx264 -preset fast -crf 20 -pix_fmt yuv420p -c:a aac -b:a 192k -t %.6f -movflags +faststart artifacts/victory-town/reclaimed-foundry-credits.mp4",ffmpeg,(double)TOWN_MUSIC_SECONDS);
#ifdef _WIN32
  char windowsCommand[2052];snprintf(windowsCommand,sizeof(windowsCommand),"\"%s\"",command);movie=_popen(windowsCommand,"wb");
#else
  movie=popen(command,"w");
#endif
  assert(movie);
 }
 SetRandomSeed(73021);srand(73021);init_pbd_thread_pool();assert(initialize_physics_backend());character_renderer_init();weapon_renderer_init();duck_renderer_init(NULL);camp_init();
 clear_all_chunks();clear_world_voxels();activePlayers=0;worldPeriodic=false;worldFloorHalfSize=40;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;foundry_selection_seed(73);buildFoundryHoldout();init_static_hash();victory_town_init();victory_town_validate();clear_all_chunks();clear_world_voxels();
 if(fallback){Shader sh={.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};for(int i=0;i<T_COUNT;i++)town.models[i].materials[0].shader=sh;for(int i=0;i<SP_COUNT;i++)storyModels[i].materials[0].shader=sh;}
 RenderTexture2D a=LoadRenderTexture(width,height),b=LoadRenderTexture(width,height),out=LoadRenderTexture(width,height);WeaponBloom glow={0};assert(weapon_bloom_resize(&glow,width,height));
 FILE *manifest=fopen("artifacts/victory-town/scenes.json","w");assert(manifest);fprintf(manifest,"[");
 for(int i=0;i<SC_COUNT;i++)fprintf(manifest,"%s{\"id\":\"%s\",\"title\":\"%s\",\"frames\":%d,\"playback_seconds\":%.6f}",i?",":"",townScenes[i].id,townScenes[i].title,(int)(townScenes[i].duration*20),(double)(townScenes[i].duration*town_playback_scale()));
 fprintf(manifest,"]");fclose(manifest);
 double times[1280]={0};float time=0;bool pause=false;int fixed=-1,frame=0;
 if(audioReady)PlayMusicStream(soundtrack);
 while(!WindowShouldClose()){
 bool seek=false;
 if(IsKeyPressed(KEY_SPACE)){pause=!pause;if(audioReady){if(pause)PauseMusicStream(soundtrack);else ResumeMusicStream(soundtrack);}}
 if(IsKeyPressed(KEY_B))bloom=!bloom;
 if(IsKeyPressed(KEY_ZERO)){fixed=-1;time=0;seek=true;}
 for(int k=0;k<SC_COUNT;k++)if(IsKeyPressed(KEY_ONE+k)){fixed=k;time=k*8;seek=true;}
 float dt=GetFrameTime();
 if(IsKeyDown(KEY_RIGHT)){time+=dt*3/town_playback_scale();seek=true;}
 if(IsKeyDown(KEY_LEFT)){time=fmaxf(0,time-dt*3/town_playback_scale());seek=true;}
 if(creditStills)time=(frame+.5f)*town_duration()/5;
 else if(capture||bench)time=qa?(frame/4)*8+(frame%4==0?0:frame%4==1?3:frame%4==2?5.6f:7.95f):frames==8?frame*8+3:frame/20.f;
 else if(video)time=frame/20.f/town_playback_scale();
 else {
  if(audioReady){UpdateMusicStream(soundtrack);
   if(seek){time=fmodf(time,town_duration());SeekMusicStream(soundtrack,time*town_playback_scale());}
   else if(!pause){if(!IsMusicStreamPlaying(soundtrack)){PlayMusicStream(soundtrack);time=0;}else time=GetMusicTimePlayed(soundtrack)/town_playback_scale();}
  }else if(!pause&&!seek)time+=dt/town_playback_scale();
  if(fixed>=0&&(time<fixed*8||time>=(fixed+1)*8)){time=fixed*8;if(audioReady)SeekMusicStream(soundtrack,time*town_playback_scale());}
 }
 float local;int timelineScene=town_timeline(time,&local);int scene=fixed>=0?fixed:timelineScene;float blend=fixed<0&&!qa&&frames!=8&&(!individual||townScenes[scene].oneWay)?fmaxf(0,local-(townScenes[scene].duration-1)):0;TownSample s=victory_town_sample(scene,local);
 assert(voxel_count==0);assert(isfinite(s.camera.position.x)&&isfinite(s.camera.position.y));
 double start=GetTime();town_render(a,&glow,s,bloom);if(blend>0)town_render(b,&glow,victory_town_sample(individual?scene:(scene+1)%SC_COUNT,0),bloom);
 BeginTextureMode(out);ClearBackground(BLACK);DrawTextureRec(a.texture,(Rectangle){0,0,width,-height},(Vector2){0},WHITE);if(blend>0)DrawTextureRec(b.texture,(Rectangle){0,0,width,-height},(Vector2){0},Fade(WHITE,blend));town_credits_draw(time,width,height);EndTextureMode();
 glFinish();if(frame<1280)times[frame]=(GetTime()-start)*1000;
 if(video){Image im=LoadImageFromTexture(out.texture);ImageFlipVertical(&im);ImageFormat(&im,PIXELFORMAT_UNCOMPRESSED_R8G8B8);size_t bytes=(size_t)width*height*3;assert(fwrite(im.data,1,bytes,movie)==bytes);UnloadImage(im);}
 if(capture){Image im=LoadImageFromTexture(out.texture);ImageFlipVertical(&im);char path[256];snprintf(path,sizeof(path),creditStills?"artifacts/victory-town/credit-%d.png":qa?"artifacts/victory-town/qa/%02d.png":frames==8?"artifacts/victory-town/scene-%d.png":individual?"artifacts/victory-town/scenes/%04d.png":"artifacts/victory-town/frames/%04d.png",frame);assert(ExportImage(im,path));UnloadImage(im);}
 BeginDrawing();ClearBackground(BLACK);DrawTextureRec(out.texture,(Rectangle){0,0,width,-height},(Vector2){0},WHITE);if(!capture&&!bench&&!video){DrawText(townScenes[scene].title,26,24,24,RAYWHITE);DrawText("1-8 scenes | 0 montage | Space pause | arrows scrub | B bloom",26,height-34,18,RAYWHITE);}EndDrawing();
 frame++;if(((capture||video)&&frame==frames)||(bench&&frame==1280))break;
 }
 if(bench){qsort(times+40,1240,sizeof(double),compare_times);printf("1280x720 render + GPU completion, bloom=%d: median %.3f ms; p95 %.3f ms (1240 frames, 40 warmup)\n",bloom,times[660],times[1218]);}
 if(movie){
#ifdef _WIN32
  int result=_pclose(movie);
#else
  int result=pclose(movie);
#endif
  assert(result==0);puts("Saved soundtrack-synchronized MP4");
 }
 if(audioReady){StopMusicStream(soundtrack);UnloadMusicStream(soundtrack);}
 if(IsAudioDeviceReady())CloseAudioDevice();
 weapon_bloom_unload(&glow);UnloadRenderTexture(a);UnloadRenderTexture(b);UnloadRenderTexture(out);victory_town_shutdown();camp_shutdown();duck_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();shutdown_world_visuals();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
