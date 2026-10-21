#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static int compare(const void *a,const void *b){double x=*(double*)a,y=*(double*)b;return (x>y)-(x<y);}
int main(int argc,char**argv){assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1600,900,"Progression preview");init_pbd_thread_pool();assert(initialize_physics_backend());weapon_renderer_init();character_renderer_init();
 gameMode=GAME_MODE_FIREFIGHT;firefightPractice=false;firefightFoundry=true;firefightCheckpoints=true;multiplayerActivePlayers=4;ResetGame();gameState=GAME_STATE_PLAYING;firefightWaveBannerTimer=0;for(int i=0;i<4;i++){players[i].pos=(Vector3){(i%2?3:-3),3,4+(i/2)*2};players[i].yaw=0;players[i].pitch=-5;}
 RenderTexture2D views[4]={0};int count=0,w=0,h=0;double times[40];
 for(int phase=0;phase<2;phase++){
  if(phase){assert(fp_reserve(FP_RULE_COUNT));fpCount=FP_RULE_COUNT;for(int n=0;n<FP_RULE_COUNT;n++){fpItems[n]=(FfItem){.id=fpNextId++,.type=fpRules[n].type,.owner=-1,.slot=-1,.pos={(n%5-2)*2.f,3,(n/5-1)*3.f-4}};fp_ground(&fpItems[n]);}
   fp_queue_unlock(PICKUP_RAILGUN);fpBannerRead=1;fpBannerTimer=3.5f;
  }
  for(int frame=0;frame<45;frame++){double before=GetTime();render_gameplay_view(views,&count,&w,&h,false);rlDrawRenderBatchActive();glFinish();if(frame>=5)times[frame-5]=(GetTime()-before)*1000;assert(count==4);}
  qsort(times,40,sizeof(double),compare);fprintf(stderr,"Four-view progression phase%d median %.3f p95 %.3f ms\n",phase,times[20],times[38]);
  if(phase){for(int n=0;n<count;n++){Image im=LoadImageFromTexture(views[n].texture);ImageFlipVertical(&im);ExportImage(im,TextFormat(".build/playtest-run/firefight-earned-%d.png",n));UnloadImage(im);}}
 }
 // A long run must exceed the old fixed pickup pool without losing IDs.
 assert(fp_reserve(2048));fpCount=2048;for(size_t n=0;n<fpCount;n++)fpItems[n]=(FfItem){.id=(uint32_t)n+1,.type=fpRules[n%FP_RULE_COUNT].type,.owner=-1,.slot=-1,.pos={25,.8f,25}};
 double begin=GetTime();for(int n=0;n<240;n++)fp_tick(1.f/60);fprintf(stderr,"2048 persistent drops: %.3f ms average simulation, %zu bytes registry capacity\n",(GetTime()-begin)*1000/240,fpCapacity*sizeof(FfItem));assert(fpCount==2048);for(int n=0;n<2048;n++)assert(fpItems[n].id==(uint32_t)n+1&&v_isfinite(fpItems[n].pos));
 for(int frame=0;frame<45;frame++){double before=GetTime();render_gameplay_view(views,&count,&w,&h,false);rlDrawRenderBatchActive();glFinish();if(frame>=5)times[frame-5]=(GetTime()-before)*1000;}
 qsort(times,40,sizeof(double),compare);fprintf(stderr,"2048-drop pile four-view median %.3f p95 %.3f ms\n",times[20],times[38]);
 for(int i=0;i<count;i++){UnloadRenderTexture(views[i]);weapon_bloom_unload(&weaponBlooms[i]);}fp_shutdown();weapon_renderer_shutdown();character_renderer_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Progression four-player preview passed");return 0;}
