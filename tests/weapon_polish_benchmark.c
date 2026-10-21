#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
#if defined(__APPLE__)
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif
bool ffTestReviveHeld[4];
static int cmp(const void *a,const void *b){double x=*(double*)a,y=*(double*)b;return (x>y)-(x<y);}
int main(int argc,char **argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1280,720,"Weapon render benchmark");SetTargetFPS(0);init_pbd_thread_pool();assert(initialize_physics_backend());
 gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;multiplayerActivePlayers=4;firefightPractice=false;SetRandomSeed(1234);ResetGame();gameState=GAME_STATE_PLAYING;start_firefight_wave(5);for(int n=0;n<30;n++)update_firefight_logic(2.1f);
 weaponRenderingReady=weapon_renderer_init();character_renderer_init();drawWeaponModels=true;weaponBloomEnabled=true;
 int guns[]={EQ_SHOTGUN,EQ_CUTTER,EQ_RAIL,EQ_ROCKET};
 for(int i=0;i<4;i++){equipment[i].gun[0]=guns[i];players[i].pos=(Vector3){-3+i*2.f,4,8};players[i].yaw=0;players[i].pitch=-4;players[i].respawn_timer=0;players[i].matter=100;combatantLife[i]=LIFE_ALIVE;weaponVisuals[i].charge=.7f;weaponVisuals[i].activity=1;}
 RenderTexture2D views[4]={0};int count=0,w=0,h=0;double samples[180];
 for(int frame=0;frame<240;frame++){double t=GetTime();render_gameplay_view(views,&count,&w,&h,false);glFinish();if(frame>=60)samples[frame-60]=(GetTime()-t)*1000;}
 qsort(samples,180,sizeof(double),cmp);printf("WEAPON_BENCH four-player Foundry 1280x720 bloom: median %.3f p95 %.3f ms\n",samples[90],samples[171]);
 for(int i=0;i<count;i++)UnloadRenderTexture(views[i]);for(int i=0;i<4;i++)weapon_bloom_unload(&weaponBlooms[i]);shutdown_world_visuals();weapon_renderer_shutdown();character_renderer_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
