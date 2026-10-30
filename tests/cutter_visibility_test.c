#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
int main(void){
 SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Incoming cutter visibility");
 activePlayers=2;worldPeriodic=false;eq_reset_all();weapon_renderer_init();
 players[0].pos=(Vector3){0,2,0};players[1].pos=(Vector3){0,2,-6};players[1].yaw=180;
 equipment[1].gun[0]=EQ_CUTTER;equipment[1].held=true;equipment[1].beamEnd=(Vector3){0,2,-.1f};
 drawViewPlayerIndex=0;
 assert(eq_cutter_contact_visibility(1,equipment[1].beamEnd)==0);
 assert(eq_cutter_contact_visibility(1,(Vector3){0,2,-4})==1);
 assert(eq_cutter_contact_visibility(0,equipment[1].beamEnd)==1);
 float middle=eq_cutter_contact_visibility(1,(Vector3){0,2,-2.25f});assert(middle>.49f && middle<.51f);
 for(int small=0;small<2;small++)for(int emission=0;emission<2;emission++)for(int fallback=0;fallback<2;fallback++){
  int w=small?480:960,h=small?270:540;RenderTexture2D rt=LoadRenderTexture(w,h);
  eq_beam_shutdown();if(fallback)eqBeamTried=true;
  equipmentEmission=emission;Camera3D cam={players[0].pos,{0,2,-10},{0,1,0},60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){24,32,45,255});BeginMode3D(cam);eq_world();EndMode3D();EndTextureMode();
  Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);
  Color *pixels=LoadImageColors(im);int bright=0;
  for(int y=h/2-h/10;y<h/2+h/10;y++)for(int x=w/2-w/10;x<w/2+w/10;x++){Color c=pixels[y*w+x];bright+=c.r>180 && c.g>180 && c.b>140;}
  assert(bright<w*h*.005f);UnloadImageColors(pixels);
  ExportImage(im,TextFormat("artifacts/cutter-incoming-%d-%d-%d.png",small,emission,fallback));UnloadImage(im);UnloadRenderTexture(rt);
 }
 equipmentEmission=false;eq_beam_shutdown();weapon_renderer_shutdown();CloseWindow();puts("Cutter camera visibility, color/emission, fallback and viewport captures passed");return 0;
}
