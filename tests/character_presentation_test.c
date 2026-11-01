/* Presentation lighting must be opt-in and restore the existing character look. */
#include "../character_renderer.c"
#include <assert.h>
#include <stdio.h>
static Image sample(RenderTexture2D target,Vector3 tint,float strength,bool emission){
 character_presentation_light(tint,strength);
 BeginTextureMode(target);ClearBackground(BLACK);BeginMode3D((Camera3D){{2,1.3f,-3},{0,.5f,0},{0,1,0},40,CAMERA_PERSPECTIVE});
 CharacterPose p=character_pose((Vector3){0,1,0},0,0,0,true,false,SKYBLUE,0);character_draw(p,emission);if(!emission)character_draw_idle_arms(p);EndMode3D();EndTextureMode();return LoadImageFromTexture(target.texture);
}
int main(void){SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Character presentation isolation");assert(character_renderer_init());RenderTexture2D rt=LoadRenderTexture(320,240);
 Image base=sample(rt,(Vector3){1,1,1},0,false),warm=sample(rt,(Vector3){1.3f,.8f,.6f},1,false),restored=sample(rt,(Vector3){1,1,1},0,false);
 assert(base.format==PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);assert(memcmp(base.data,restored.data,320*240*4)==0);assert(memcmp(base.data,warm.data,320*240*4)!=0);
 Image a=sample(rt,(Vector3){1,1,1},0,true),b=sample(rt,(Vector3){1.3f,.8f,.6f},1,true);assert(memcmp(a.data,b.data,320*240*4)==0);
 UnloadImage(base);UnloadImage(warm);UnloadImage(restored);UnloadImage(a);UnloadImage(b);UnloadRenderTexture(rt);character_renderer_shutdown();CloseWindow();puts("Presentation lighting isolation passed");return 0;}
