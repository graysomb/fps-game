/* Deterministic production-mesh validation and animation contact sheets.
 * Includes renderer internals to check geometry budgets, packed parts and cleanup. */
#include "../weapon_renderer.c"
#include "../character_renderer.h"
#define main original_preview_main
#include "render_weapon_preview.c"
#undef main
static const int kinds[]={WEAPON_IMPULSE_SHOTGUN,WEAPON_CUTTER,WEAPON_ADHESIVE_GRENADE,WEAPON_RAILGUN,WEAPON_ROCKET_LAUNCHER,WEAPON_REINFORCED_BUILDER};
static const char *names[]={"shotgun","cutter","adhesive","railgun","rocket","builder"};
static WeaponVisual sample(int kind,float t){
 WeaponVisual v;weapon_visual_reset(&v);
 if(kind==WEAPON_IMPULSE_SHOTGUN){float shot=t<2?.6f:2.f;if(t>=shot){v.shot_time=shot;v.cooldown=fmaxf(0,1-(t-shot)/.8f);}}
 if(kind==WEAPON_CUTTER)v.activity=t>=.6f && t<2.4f;
 if(kind==WEAPON_ADHESIVE_GRENADE){v.attached=t>=.8f;v.activity=v.attached;v.charge=Clamp((t-.8f)/2,0,1);}
 if(kind==WEAPON_RAILGUN){v.charge=t<1.25f?Clamp((t-.6f)/.65f,0,1):0;if(t>=1.25f){v.shot_time=1.25f;v.cooldown=fmaxf(0,1-(t-1.25f)/1.8f);}}
 if(kind==WEAPON_ROCKET_LAUNCHER && t>=.6f){v.shot_time=.6f;v.cooldown=fmaxf(0,1-(t-.6f)/1.5f);}
 if(kind==WEAPON_REINFORCED_BUILDER){v.deployment=Clamp((t-.5f)/.15f,0,1)*Clamp((2.4f-t)/.15f,0,1);v.build_pulse=t>=2?fmaxf(0,1-(t-2)/.6f):0;}
 return v;
}
static void meshes(void){
 for(int k=0;k<6;k++){
  Mesh *m=&wr_parts[wr_models[kinds[k]].part].meshes[0];assert(m->triangleCount<=2000 && m->vertexCount>0);
  for(int v=0;v<m->vertexCount;v++){
   Vector3 p,n;memcpy(&p,m->vertices+v*3,sizeof(p));memcpy(&n,m->normals+v*3,sizeof(n));
   assert(isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z)&&fabsf(p.x)<2&&fabsf(p.y)<2&&fabsf(p.z)<3);
   assert(fabsf(Vector3Length(n)-1)<.001f);assert(m->texcoords[v*2]>=0 && m->texcoords[v*2]<104);
  }
  printf("%s: %d triangles, one cached model\n",names[k],m->triangleCount);
 }
}
static void dark_emission(RenderTexture2D scene,Camera3D cam,int kind){
 WeaponVisual v;weapon_visual_reset(&v);v.depleted=1;
 WeaponPose p={.kind=kind,.transform=MatrixIdentity()};
 BeginTextureMode(scene);ClearBackground(BLACK);BeginMode3D(cam);weapon_draw(p,&v,cam.position,10,false,WEAPON_EMISSION);EndMode3D();EndTextureMode();
 Image image=LoadImageFromTexture(scene.texture);Color *c=LoadImageColors(image);
 for(int n=0;n<image.width*image.height;n++)assert(c[n].r==0&&c[n].g==0&&c[n].b==0);
 UnloadImageColors(c);UnloadImage(image);
}
static void occlusion(RenderTexture2D scene,int kind){
 Camera3D cam={{0,0,4},{0,0,0},{0,1,0},42,CAMERA_PERSPECTIVE};
 WeaponVisual v=sample(kind,1);WeaponPose pose={.kind=kind,.transform=MatrixIdentity()};
 assert(weapon_bloom_resize(&bloom,scene.texture.width,scene.texture.height));
 BeginTextureMode(scene);ClearBackground(BLACK);BeginMode3D(cam);DrawCube((Vector3){0,0,2},5,5,.3f,GRAY);EndMode3D();EndTextureMode();
 weapon_bloom_begin(&bloom,scene);BeginMode3D(cam);weapon_draw(pose,&v,cam.position,1,false,WEAPON_EMISSION);EndMode3D();EndTextureMode();
 Image image=LoadImageFromTexture(bloom.emission.texture);Color *c=LoadImageColors(image);
 for(int n=0;n<image.width*image.height;n++)assert(c[n].r==0&&c[n].g==0&&c[n].b==0);
 UnloadImageColors(c);UnloadImage(image);
}
int main(int argc,char **argv){
 const char *out=argc>1?argv[1]:"artifacts/weapon-polish/animation";
 SetConfigFlags(FLAG_WINDOW_HIDDEN);SetTraceLogLevel(LOG_WARNING);InitWindow(960,540,"Weapon polish");
 assert(weapon_renderer_init() && !wr_fallback);assert(character_renderer_init());meshes();test_visual_state();test_melee_pose();test_cube_motion();
 Camera3D cam={{-2.8f,1.5f,-3.5f},{0,-.1f,-.4f},{0,1,0},42,CAMERA_PERSPECTIVE};
 Camera3D fp={{0},{0,0,-1},{0,1,0},60,CAMERA_PERSPECTIVE};
 RenderTexture2D scene=LoadRenderTexture(480,480),first=LoadRenderTexture(640,360),sheet=LoadRenderTexture(960,480);
 use_bloom=true;
 for(int k=0;k<6;k++){
  WeaponPose pose={.kind=kinds[k],.transform=MatrixIdentity()};
  WeaponPose hand=weapon_pose(kinds[k],(Vector3){0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},true,16.f/9,NULL,0);
  dark_emission(scene,cam,kinds[k]);occlusion(scene,kinds[k]);
  for(int n=0;n<(argc>2?1:72);n++){
   render_time=n/20.f;WeaponVisual visual=sample(kinds[k],render_time);
   render(scene,cam,pose,&visual,false,false,false);render(first,fp,hand,&visual,false,false,true);
   BeginTextureMode(first);BeginMode3D(fp);character_draw_paws(hand,BLUE,false);EndMode3D();EndTextureMode();
   BeginTextureMode(sheet);ClearBackground((Color){20,31,49,255});
   DrawTexturePro(scene.texture,(Rectangle){0,0,480,-480},(Rectangle){0,0,480,480},(Vector2){0},0,WHITE);
   DrawTexturePro(first.texture,(Rectangle){0,0,640,-360},(Rectangle){480,105,480,270},(Vector2){0},0,WHITE);
   DrawText(names[k],24,22,25,WHITE);DrawText("First-person / paws",510,68,20,LIGHTGRAY);
   DrawText(render_time<.5f?"IDLE":render_time<2.4f?"ACTION":"RECOVERY",24,442,20,LIGHTGRAY);
   EndTextureMode();char file[80];snprintf(file,sizeof(file),"%s-%03d",names[k],n);save(sheet,out,file);
  }
 }
 UnloadRenderTexture(sheet);UnloadRenderTexture(first);UnloadRenderTexture(scene);weapon_bloom_unload(&bloom);character_renderer_shutdown();weapon_renderer_shutdown();assert(!wr_ready);assert(weapon_renderer_init());weapon_renderer_shutdown();CloseWindow();puts("Weapon polish geometry, depletion, animation and cleanup passed");
}
