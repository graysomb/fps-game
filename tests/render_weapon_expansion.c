/* Standalone presentation sandbox. No fps_ray.c, equipment, physics or networking.
 * Includes private mesh data solely for resource/budget assertions. */
#include "../weapon_renderer.c"
#include "../character_renderer.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>
#ifdef _WIN32
__declspec(dllimport) void __stdcall glFinish(void);
#endif
static const WeaponKind items[]={WEAPON_GLOB_GUN,WEAPON_BLACK_HOLE_BOMB,WEAPON_GLOCK,WEAPON_IGLOO_GRENADE,WEAPON_DMR};
static const char *slugs[]={"glob","black-hole","glock","igloo","dmr"};
static WeaponBloom preview_blooms[4];
static WeaponBloom *preview_bloom(RenderTexture2D rt){
 for(int i=0;i<4;i++)if(preview_blooms[i].width==rt.texture.width && preview_blooms[i].height==rt.texture.height)return &preview_blooms[i];
 for(int i=0;i<4;i++)if(!preview_blooms[i].emission.id){assert(weapon_bloom_resize(&preview_blooms[i],rt.texture.width,rt.texture.height));return &preview_blooms[i];}
 assert(!"Too many preview target sizes");return &preview_blooms[0];
}
static bool bloom_on=true,effects_on=true,depleted=false,wall=false;
static const Color backdrop={16,27,43,255};
static Matrix at(Vector3 p){return MatrixTranslate(p.x,p.y,p.z);}
static WeaponVisual sample(int k,float t){
 WeaponVisual v;weapon_visual_reset(&v);v.depleted=depleted;
 if(k==WEAPON_GLOB_GUN && t>=.7f&&t<2.4f){v.shot_time=.7f+floorf((t-.7f)/.18f)*.18f;v.activity=1;}
 if(k==WEAPON_GLOCK||k==WEAPON_DMR){float shot=t<3.7f?1.4f:3.7f;if(t>=shot)v.shot_time=shot;}
 if(k==WEAPON_BLACK_HOLE_BOMB)v.deployment=Clamp((t-.6f)/.8f,0,1)*Clamp((6-t)/.5f,0,1);
 if(k==WEAPON_IGLOO_GRENADE)v.deployment=Clamp((t-.7f)/.7f,0,1);
 return v;
}
static void robot(Matrix root,WeaponRenderPass pass){
 // Scripted target geometry, not an AI actor; reused for both effect demos.
 Color body=pass==WEAPON_EMISSION?BLACK:(Color){165,176,188,255};
 Color dark=pass==WEAPON_EMISSION?BLACK:(Color){45,54,65,255};
 Vector3 positions[]={{0,.7f,0},{0,1.18f,0},{-.28f,.7f,0},{.28f,.7f,0},{-.14f,.23f,0},{.14f,.23f,0}};
 Vector3 sizes[]={{.43f,.50f,.27f},{.34f,.32f,.30f},{.14f,.48f,.15f},{.14f,.48f,.15f},{.16f,.46f,.19f},{.16f,.46f,.19f}};
 if(pass==WEAPON_GLASS)return;
 for(int i=0;i<6;i++)DrawCubeV(Vector3Transform(positions[i],root),sizes[i],i<2?body:dark);
 DrawCube(Vector3Transform((Vector3){0,1.2f,-.156f},root),.19f,.07f,.02f,(Color){245,155,40,255});
}
static void demo(int kind,float t,Vector3 eye,WeaponRenderPass pass,bool cut){
 WeaponEffect e={.transform=MatrixIdentity(),.age=fmaxf(0,t-1),.lifetime=20,.progress=Clamp((t-1)/3,0,1),.seed=17};
 if(kind==WEAPON_GLOB_GUN){
  Matrix target=at((Vector3){sinf(t*.8f)*.35f,0,0});robot(target,pass);
  if(!effects_on)return;
  if(t>=.7f&&t<1.7f){e.age=t-.7f;e.lifetime=1;e.transform=at((Vector3){-.8f+.8f*e.age,.85f,-2+1.82f*e.age});weapon_draw_glob(&e,eye,t,pass);}
  if(t>=1.7f){e.attached=true;e.age=t-1.7f;e.lifetime=5;e.transform=MatrixMultiply(at((Vector3){-.06f,.85f,-.16f}),target);weapon_draw_glob(&e,eye,t,pass);}
 }else if(kind==WEAPON_BLACK_HOLE_BOMB){
  WeaponVisual v=sample(kind,t);WeaponPose p={.kind=kind,.transform=at((Vector3){0,.9f,0})};weapon_draw(p,&v,eye,t,false,pass);
  if(!effects_on)return;
  e.age=t-1;e.lifetime=5;e.transform=at((Vector3){0,.9f,0});weapon_draw_gravity(&e,eye,t,pass);
  float pull=Clamp((t-1)/4,0,1);float r=2.4f*(1-pull);
  if(t<5){robot(at((Vector3){r,.4f*pull,.25f}),pass);if(pass==WEAPON_COLOR)for(int i=0;i<9;i++){
   float a=i*2.4f+t*pull;Vector3 p={cosf(a)*r,.12f+pull*.8f,sinf(a)*r};DrawCube(p,.15f,.15f,.15f,(Color){130,141,160,255});}}
 }else if(kind==WEAPON_IGLOO_GRENADE){
  CharacterPose rat=character_pose((Vector3){0,1.30f,0},0,0,0,true,false,SKYBLUE,t);
  if(pass!=WEAPON_GLASS){character_draw(rat,pass==WEAPON_EMISSION);if(pass==WEAPON_COLOR)character_draw_idle_arms(rat);}else character_draw_visor(rat);
  WeaponVisual v=sample(kind,t);WeaponPose p={.kind=kind,.transform=MatrixMultiply(MatrixScale(.45f,.45f,.45f),at((Vector3){0,.48f,-1.4f}))};weapon_draw(p,&v,eye,t,false,pass);
  if(!effects_on)return;
  e.cutaway=cut;weapon_draw_construction(&e,eye,t,pass);
 }else if(kind==WEAPON_DMR){
  float shot=t<3.7f?1.4f:3.7f,age=t-shot;
  if(pass!=WEAPON_GLASS)for(int y=0;y<3;y++)for(int x=-1;x<=1;x++){
   if(x==0&&y==1&&age>=0&&effects_on)continue;
   float impulse=effects_on&&age>=0?expf(-age*2)*age*.35f:0;
   Vector3 p={x*(.4f+impulse),.2f+y*.4f,0};DrawCube(p,.39f,.39f,.39f,pass==WEAPON_EMISSION?BLACK:(Color){145,157,173,255});
  }
  if(!effects_on)return;
  if(age>=0){e.age=age;e.lifetime=1.2f;e.transform=at((Vector3){0,.6f,-.21f});weapon_draw_dmr_impact(&e,eye,t,pass);}
  if(age>=0&&age<.11f&&pass!=WEAPON_GLASS)DrawCylinderEx((Vector3){0,.6f,-2},(Vector3){0,.6f,-.2f},.009f,.009f,5,SKYBLUE);
 }else if(kind==WEAPON_GLOCK){
  if(pass==WEAPON_COLOR){DrawCube((Vector3){0,.7f,0},.9f,1.4f,.15f,(Color){62,75,89,255});DrawCube((Vector3){0,.8f,-.09f},.25f,.25f,.02f,LIGHTGRAY);}
 }
}
static void world(RenderTexture2D rt,int kind,float t,bool demonstrate,int view,bool cut){
 WeaponBloom *fx_bloom=preview_bloom(rt);
 Camera3D cam={{-2.8f,1.5f,-3.5f},{0,-.12f,-.3f},{0,1,0},42,CAMERA_PERSPECTIVE};
 if(view==1)cam.position=(Vector3){0,.15f,-4.5f};
 if(view==2)cam.position=(Vector3){-4,.4f,0};
 if(view==3)cam.position=(Vector3){0,.6f,4.3f};
 if(demonstrate){cam=(Camera3D){{-3.8f,2.5f,-5.0f},{0,.65f,0},{0,1,0},45,CAMERA_PERSPECTIVE};if(kind==WEAPON_IGLOO_GRENADE)cam=(Camera3D){{-6,5,-7},{0,1.1f,0},{0,1,0},46,CAMERA_PERSPECTIVE};}
 WeaponVisual v=sample(kind,t);WeaponPose p={.kind=kind,.transform=MatrixIdentity()};
 BeginTextureMode(rt);ClearBackground(backdrop);BeginMode3D(cam);
 DrawPlane((Vector3){0,demonstrate?-.02f:-1.12f,0},(Vector2){30,30},(Color){25,38,54,255});
 if(wall)DrawCube((Vector3){0,0,-2},10,10,.2f,GRAY);
 if(demonstrate)demo(kind,t,cam.position,WEAPON_COLOR,cut);else weapon_draw(p,&v,cam.position,t,false,WEAPON_COLOR);
 EndMode3D();EndTextureMode();
 if(bloom_on && !wr_fallback){
  assert(weapon_bloom_resize(fx_bloom,rt.texture.width,rt.texture.height));weapon_bloom_begin(fx_bloom,rt);BeginMode3D(cam);
  if(demonstrate)demo(kind,t,cam.position,WEAPON_EMISSION,cut);else weapon_draw(p,&v,cam.position,t,false,WEAPON_EMISSION);
  EndMode3D();weapon_bloom_composite(fx_bloom,rt);EndTextureMode();
 }
 BeginTextureMode(rt);BeginMode3D(cam);if(demonstrate)demo(kind,t,cam.position,WEAPON_GLASS,cut);else weapon_draw(p,&v,cam.position,t,false,WEAPON_GLASS);EndMode3D();EndTextureMode();
}
static void first_person(RenderTexture2D rt,int kind,float t){
 WeaponBloom *fx_bloom=preview_bloom(rt);
 Camera3D cam={{0,1.2f,-4},{0,1.2f,0},{0,1,0},60,CAMERA_PERSPECTIVE};
 bool zoom=kind==WEAPON_DMR&&t>=2&&t<3.5f;
 if(zoom)cam.fovy=30;
 BeginTextureMode(rt);ClearBackground(backdrop);BeginMode3D(cam);DrawGrid(20,1);
 DrawCube((Vector3){0,1.2f,5},1.2f,1.8f,.7f,(Color){110,125,145,255});EndMode3D();
 // Separate fixed-FOV hands, exactly as the game's clear-view aiming path.
 weapon_clear_depth();cam=(Camera3D){{0},{0,0,-1},{0,1,0},60,CAMERA_PERSPECTIVE};
 WeaponVisual v=sample(kind,t);WeaponPose p=weapon_pose(kind,(Vector3){0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},true,(float)rt.texture.width/rt.texture.height,NULL,weapon_shot_amount(&v,t));
 BeginMode3D(cam);weapon_draw(p,&v,cam.position,t,false,WEAPON_COLOR);character_draw_paws(p,SKYBLUE,false);EndMode3D();EndTextureMode();
 if(bloom_on&&!wr_fallback){weapon_bloom_resize(fx_bloom,rt.texture.width,rt.texture.height);weapon_bloom_begin(fx_bloom,rt);BeginMode3D(cam);weapon_draw(p,&v,cam.position,t,false,WEAPON_EMISSION);character_draw_paws(p,SKYBLUE,true);EndMode3D();weapon_bloom_composite(fx_bloom,rt);EndTextureMode();}
 BeginTextureMode(rt);BeginMode3D(cam);weapon_draw(p,&v,cam.position,t,false,WEAPON_GLASS);EndMode3D();
 if(kind==WEAPON_DMR)weapon_draw_dmr_reticle(rt.texture.width,rt.texture.height,1);else DrawCircle(rt.texture.width/2,rt.texture.height/2,2,WHITE);
 // Pixel-sized HUD reference stays fixed at all zoom levels.
 DrawRectangle(18,rt.texture.height-32,28,16,SKYBLUE);DrawRectangleLines(52,rt.texture.height-32,28,16,GRAY);
 DrawText(zoom?"DMR / 2x preview":"FIRST PERSON",18,18,16,LIGHTGRAY);EndTextureMode();
}
static void save_png(RenderTexture2D rt,const char *dir,const char *name){
 char file[1024];snprintf(file,sizeof(file),"%s/%s.png",dir,name);Image img=LoadImageFromTexture(rt.texture);ImageFlipVertical(&img);assert(ExportImage(img,file));UnloadImage(img);
}
static void blit(RenderTexture2D rt,Rectangle dest){DrawTexturePro(rt.texture,(Rectangle){0,0,(float)rt.texture.width,-(float)rt.texture.height},dest,(Vector2){0},0,WHITE);}
static void geometry_checks(void){
 assert(WEAPON_REINFORCED_BUILDER==8 && WEAPON_GLOB_GUN==9);
 for(int i=0;i<5;i++){
  Mesh m=wr_parts[wr_models[items[i]].part].meshes[0];int max=items[i]==WEAPON_GLOCK?4000:2000;assert(m.triangleCount<=max);
  int triangles=m.triangleCount;if(items[i]==WEAPON_GLOB_GUN)triangles+=wr_parts[WR_GEL_GLASS].meshes[0].triangleCount;
  if(items[i]==WEAPON_DMR)triangles+=wr_parts[WR_SCOPE_GLASS].meshes[0].triangleCount;
  assert(triangles<=max);printf("GEOMETRY %s %d triangles\n",slugs[i],triangles);
  WeaponAnchors a=weapon_anchors(items[i]);assert(isfinite(a.muzzle.z));
 }
 for(int i=WR_GLOB;i<WR_PART_COUNT;i++)for(int v=0;v<wr_parts[i].meshes[0].vertexCount;v++){
  Mesh m=wr_parts[i].meshes[0];Vector3 p,n;memcpy(&p,m.vertices+v*3,sizeof(p));memcpy(&n,m.normals+v*3,sizeof(n));
  assert(isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z));if(fabsf(Vector3Length(n)-1)>=.005f){fprintf(stderr,"Bad normal part=%d vertex=%d norm=%f\n",i,v,Vector3Length(n));abort();}
 }
 // Flood fill the full grid to catch unsupported/disconnected courses and flue closures.
 bool visited[WX_HIGH][WX_GRID][WX_GRID]={0};int queue[WX_HIGH*WX_GRID*WX_GRID][3],head=0,tail=0,total=0;
 for(int y=0;y<WX_HIGH;y++)for(int z=0;z<WX_GRID;z++)for(int x=0;x<WX_GRID;x++)if(wx_cell(x,y,z,false)){
  total++;if(tail==0){queue[tail][0]=x;queue[tail][1]=y;queue[tail++][2]=z;visited[y][z][x]=true;}}
 int d[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
 while(head<tail){int x=queue[head][0],y=queue[head][1],z=queue[head++][2];for(int i=0;i<6;i++){
  int a=x+d[i][0],b=y+d[i][1],c=z+d[i][2];if(wx_cell(a,b,c,false)&&!visited[b][c][a]){visited[b][c][a]=true;queue[tail][0]=a;queue[tail][1]=b;queue[tail++][2]=c;}}}
 assert(tail==total);for(int y=1;y<=10;y++)assert(!wx_cell(9,y,9,false));
 printf("IGLOO %d connected cells, %d cached surface triangles, open chimney\n",total,wr_parts[WR_FX_IGLOO].meshes[0].triangleCount);
}
static void emission_checks(void){
 RenderTexture2D rt=LoadRenderTexture(320,240);Camera3D cam={{0,0,-4},{0,0,0},{0,1,0},45,CAMERA_PERSPECTIVE};
 for(int k=0;k<5;k++){
  WeaponVisual v;weapon_visual_reset(&v);v.depleted=1;WeaponPose p={.kind=items[k],.transform=MatrixIdentity()};
  for(int occluded=0;occluded<2;occluded++){
   BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(cam);
   if(occluded){DrawCube((Vector3){0,0,-2},8,8,.2f,BLACK);v.depleted=0;}
   weapon_draw(p,&v,cam.position,4,false,WEAPON_EMISSION);EndMode3D();EndTextureMode();
   Image im=LoadImageFromTexture(rt.texture);Color *px=LoadImageColors(im);for(int j=0;j<im.width*im.height;j++)assert(px[j].r==0&&px[j].g==0&&px[j].b==0);UnloadImageColors(px);UnloadImage(im);
  }
 }
 // Completed glob effects must leave no pixels, independent of the source's power.
 WeaponEffect e={.transform=MatrixIdentity(),.age=5,.lifetime=5,.attached=true};
 BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(cam);weapon_draw_glob(&e,cam.position,7,WEAPON_COLOR);EndMode3D();EndTextureMode();
 Image im=LoadImageFromTexture(rt.texture);Color *px=LoadImageColors(im);for(int j=0;j<im.width*im.height;j++)assert(px[j].r==0&&px[j].g==0&&px[j].b==0);UnloadImageColors(px);UnloadImage(im);UnloadRenderTexture(rt);
}
static int compare_double(const void *a,const void *b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static void benchmark(const char *out){
 RenderTexture2D views[4];for(int i=0;i<4;i++)views[i]=LoadRenderTexture(640,360);
 char path[1024];snprintf(path,sizeof(path),"%s/benchmark.csv",out);FILE *file=fopen(path,"w");assert(file);fprintf(file,"run,effects,frame,ms\n");
 for(int run=0;run<3;run++)for(int mode=0;mode<2;mode++){
  double times[180];effects_on=mode;for(int frame=0;frame<240;frame++){
   double start=GetTime();float t=1.3f+(frame%120)/60.f;
   for(int i=0;i<4;i++)world(views[i],items[i==2?4:i],t,true,0,false);
   rlDrawRenderBatchActive();glFinish();double ms=(GetTime()-start)*1000;
   if(frame>=60){times[frame-60]=ms;fprintf(file,"%d,%d,%d,%.6f\n",run,mode,frame-60,ms);}
  }
  qsort(times,180,sizeof(double),compare_double);printf("BENCH run=%d effects=%d median=%.3f p95=%.3f ms\n",run,mode,times[90],times[171]);
 }
 fclose(file);for(int i=0;i<4;i++)UnloadRenderTexture(views[i]);
}
int main(int argc,char **argv){
 const char *out=argc>2?argv[2]:"artifacts/weapon-expansion";
 bool fallback=argc>1&&strcmp(argv[1],"--fallback")==0;
 bool capture=argc>1&&strcmp(argv[1],"--capture")==0,bench=argc>1&&strcmp(argv[1],"--benchmark")==0,checks=argc>1&&strcmp(argv[1],"--checks")==0;
 SetTraceLogLevel(LOG_WARNING);if(capture||bench||checks||fallback)SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(1200,720,"Weapon visual lab");
 char original[1024];snprintf(original,sizeof(original),"%s",GetWorkingDirectory());
 if(fallback){MakeDirectory(".build/weapon-empty-assets");assert(ChangeDirectory(".build/weapon-empty-assets"));}
 assert(weapon_renderer_init());if(fallback){assert(wr_fallback);assert(ChangeDirectory(original));}else{assert(!wr_fallback&&wx_shader(wr_effect)&&wx_shader(wr_glock));}
 assert(character_renderer_init());geometry_checks();emission_checks();
 if(bench)benchmark(out);
 else{
  RenderTexture2D model=LoadRenderTexture(600,480),effect=LoadRenderTexture(600,480),fp=LoadRenderTexture(640,360),sheet=LoadRenderTexture(1200,720);
  if(fallback){for(int k=0;k<5;k++){world(model,items[k],0,false,0,false);char name[80];snprintf(name,sizeof(name),"%s-fallback",slugs[k]);save_png(model,out,name);}}
  else if(capture){
   for(int k=0;k<5;k++){
    if(argc>3&&k!=atoi(argv[3]))continue;
    char name[128];for(int view=0;view<4;view++){world(model,items[k],0,false,view,false);snprintf(name,sizeof(name),"%s-view-%d",slugs[k],view);save_png(model,out,name);}
    for(int small=0;small<2;small++){RenderTexture2D full=LoadRenderTexture(small?640:1280,small?360:720);first_person(full,items[k],2.5f);snprintf(name,sizeof(name),"%s-first-person-%s",slugs[k],small?"quarter":"full");save_png(full,out,name);UnloadRenderTexture(full);}
    RenderTexture2D four=LoadRenderTexture(1280,720);
    first_person(fp,items[k],2.5f);BeginTextureMode(four);ClearBackground(backdrop);
    for(int q=0;q<4;q++)blit(fp,(Rectangle){(q%2)*640,(q/2)*360,640,360});EndTextureMode();
    snprintf(name,sizeof(name),"%s-four-player",slugs[k]);save_png(four,out,name);UnloadRenderTexture(four);
    bloom_on=false;world(model,items[k],0,false,0,false);snprintf(name,sizeof(name),"%s-no-bloom",slugs[k]);save_png(model,out,name);bloom_on=true;
    depleted=true;world(model,items[k],0,false,0,false);snprintf(name,sizeof(name),"%s-depleted",slugs[k]);save_png(model,out,name);depleted=false;
    for(int frame=0;frame<160;frame++){
     float t=frame/20.f;world(model,items[k],t,false,0,false);world(effect,items[k],t,true,0,items[k]==WEAPON_IGLOO_GRENADE&&t>=5.5f);first_person(fp,items[k],t);
     BeginTextureMode(sheet);ClearBackground(backdrop);blit(model,(Rectangle){0,40,600,480});blit(effect,(Rectangle){600,40,600,480});blit(fp,(Rectangle){600,470,440,247.5f});
     DrawText(weapon_visual_name(items[k]),24,16,25,WHITE);DrawText("SCRIPTED EFFECT / NO GAMEPLAY",620,20,16,LIGHTGRAY);
     DrawText(t<.7f?"IDLE":t<5.5f?"ACTION / DEPLOYMENT":"RECOVERY / DISAPPEARANCE",30,550,20,SKYBLUE);
     DrawText(TextFormat("%.2fs   |   visual prototype",t),30,582,18,LIGHTGRAY);
     if(items[k]==WEAPON_IGLOO_GRENADE)DrawText(t>=5.5f?"CUTAWAY / OPEN CHIMNEY":"GROUNDED COURSE REVEAL",30,620,18,LIGHTGRAY);
     EndTextureMode();snprintf(name,sizeof(name),"%s-%03d",slugs[k],frame);save_png(sheet,out,name);
    }
   }
   world(effect,WEAPON_IGLOO_GRENADE,5,true,0,false);save_png(effect,out,"igloo-exterior");world(effect,WEAPON_IGLOO_GRENADE,6,true,0,true);save_png(effect,out,"igloo-cutaway");
  }else if(checks){
   // Exercise color, emission, glass and all effect programs on either GL backend.
   for(int k=0;k<5;k++){
    world(model,items[k],2.5f,false,0,false);world(effect,items[k],2.5f,true,0,false);first_person(fp,items[k],2.5f);
    Image a=LoadImageFromTexture(effect.texture);Color *pa=LoadImageColors(a);
    // Other views may change uniforms, but cannot advance or alter this sample.
    world(model,items[(k+1)%5],1.4f,false,1,false);world(effect,items[k],2.5f,true,0,false);
    Image b=LoadImageFromTexture(effect.texture);Color *pb=LoadImageColors(b);
    assert(memcmp(pa,pb,(size_t)a.width*a.height*sizeof(Color))==0);
    UnloadImageColors(pa);UnloadImageColors(pb);UnloadImage(a);UnloadImage(b);
   }
   first_person(fp,WEAPON_DMR,0);Image a=LoadImageFromTexture(fp.texture);Color *pa=LoadImageColors(a);
   first_person(fp,WEAPON_DMR,2.5f);Image b=LoadImageFromTexture(fp.texture);Color *pb=LoadImageColors(b);
   // Texture readback is upside-down: the bottom-left HUD occupies rows 16..32.
   for(int y=16;y<32;y++)for(int x=18;x<46;x++)assert(memcmp(&pa[y*a.width+x],&pb[y*a.width+x],sizeof(Color))==0);
   UnloadImageColors(pa);UnloadImageColors(pb);UnloadImage(a);UnloadImage(b);
  }else if(!checks){
   SetTargetFPS(60);int k=0;float t=0;bool paused=false,cut=false;int view=0;
   while(!WindowShouldClose()){
    if(IsKeyPressed(KEY_SPACE))paused=!paused;if(IsKeyPressed(KEY_TAB))k=(k+1)%5;
    for(int i=0;i<5;i++)if(IsKeyPressed(KEY_ONE+i)){k=i;t=0;}
    if(IsKeyDown(KEY_RIGHT))t+=GetFrameTime()*2;if(IsKeyDown(KEY_LEFT))t=fmaxf(0,t-GetFrameTime()*2);
    if(IsKeyPressed(KEY_R))t=0;if(IsKeyPressed(KEY_B))bloom_on=!bloom_on;if(IsKeyPressed(KEY_P))depleted=!depleted;
    if(IsKeyPressed(KEY_C))cut=!cut;if(IsKeyPressed(KEY_V))view=(view+1)%4;if(IsKeyPressed(KEY_O))wall=!wall;
    if(!paused)t=fmodf(t+GetFrameTime(),8);world(model,items[k],t,false,view,false);world(effect,items[k],t,true,0,cut);first_person(fp,items[k],t);
    BeginDrawing();ClearBackground(backdrop);blit(model,(Rectangle){0,40,600,480});blit(effect,(Rectangle){600,40,600,480});blit(fp,(Rectangle){600,470,440,247.5f});
    DrawText(weapon_visual_name(items[k]),24,12,25,WHITE);DrawText(TextFormat("%.2fs %s",t,paused?"PAUSED":""),24,540,24,SKYBLUE);
    DrawText("1-5 / Tab: item   Space: pause   Arrows: scrub   R: reset",24,580,18,LIGHTGRAY);
    DrawText("V: view   B: bloom   P: depleted   C: cutaway   O: wall",24,608,18,LIGHTGRAY);EndDrawing();
   }
  }
  UnloadRenderTexture(model);UnloadRenderTexture(effect);UnloadRenderTexture(fp);UnloadRenderTexture(sheet);
 }
 for(int i=0;i<4;i++)weapon_bloom_unload(&preview_blooms[i]);character_renderer_shutdown();weapon_renderer_shutdown();assert(!wr_ready);assert(weapon_renderer_init());weapon_renderer_shutdown();CloseWindow();puts("Expansion geometry, depletion, occlusion, lifetime and cleanup passed");return 0;
}
