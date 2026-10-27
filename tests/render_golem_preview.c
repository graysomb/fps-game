/* Standalone presentation lab: no gameplay, physics or network dependencies.
 * Include the renderer to validate cached geometry and resource ownership. */
#include "../golem_renderer.c"
#include "../character_renderer.h"
#include <stdlib.h>
#ifdef _WIN32
__declspec(dllimport) void __stdcall glFinish(void);
#endif
static const char *names[]={"idle","walk","slam","sweep","throw","arm-loss","leg-loss","chest-loss","both-legs","hit"};
static bool bloom=true,effects=true,skeleton=false,wall=false,showModel=true;
static float orbit=0;
static WeaponBloom fx;
static Vector3 path_at(float t){return (Vector3){3*sinf(.35f*t),0,3*cosf(.35f*t)-3};}
static float path_yaw(float t){return (.35f*t-PI/2)*RAD2DEG;}
static Vector3 contact_at(float t,int side){return Vector3Add(path_at(t),Vector3RotateByAxisAngle((Vector3){side?.63f:-.63f,.33f,0},(Vector3){0,1,0},path_yaw(t)*DEG2RAD));}
static GolemInput base_input(int scenario,float t){
 GolemInput in=golem_input();in.time=t;in.side=1;in.scale=1.5f;
 if(scenario==1){
  in.action=GOLEM_WALK;in.position=path_at(t);in.yaw=path_yaw(t);in.travel=.7f*t;in.speed=.7f;in.feetResolved=true;
  float period=1.2f/.7f;
  for(int s=0;s<2;s++){
   float u=t/period+s*.5f,step=floorf(u),phase=u-step;
   float center=(step-s*.5f+.3f)*period;
   in.feet[s]=contact_at(center,s);in.footYaw[s]=path_yaw(center);
   if(phase>=.6f){float f=gr_ease(.6f,1,phase);in.feet[s]=Vector3Lerp(in.feet[s],contact_at(center+period,s),f);in.feet[s].y+=.18f*sinf((phase-.6f)/.4f*PI);in.footYaw[s]=Lerp(path_yaw(center),path_yaw(center+period),f);}
  }
 }
 if(scenario>=2&&scenario<=4){in.action=(GolemAction)scenario;in.actionProgress=Clamp((t-.7f)/golem_action_seconds(in.action),0,1);}
 if(scenario==5){in.action=GOLEM_SLAM;in.actionProgress=Clamp((t-.4f)/2.2f,0,1);}
 return in;
}
static void scripted_loss(GolemInput *in,int section,bool forearm,float eventTime,int scenario){
 if(in->time<eventTime)return;
 GolemInput original=base_input(scenario,eventTime);GolemPose initial=golem_sample(&original);
 GolemDamage *d=&in->damage[section];golem_capture_damage(&initial,(GolemSection)section,forearm,d);
 d->age=in->time-eventTime;d->regeneration=Clamp((d->age-1.25f)/2,0,1);d->recovery=gr_ease(3.25f,4.05f,d->age);
}
static GolemPose sample(int scenario,float t,unsigned manual){
 GolemInput in=base_input(scenario,t);
 if(scenario==4)scripted_loss(&in,GOLEM_RIGHT_ARM,true,.7f+2.4f*.60f,scenario);
 if(scenario>=5&&scenario<=8){
  int section=scenario==5?GOLEM_RIGHT_ARM:scenario==7?GOLEM_CHEST:GOLEM_LEFT_LEG;
  scripted_loss(&in,section,false,1,scenario);
  if(scenario==8)scripted_loss(&in,GOLEM_RIGHT_LEG,false,1,scenario);
 }
 if(scenario==9){in.hitJoint=GOLEM_TORSO;in.hitAge=fmodf(t,1.5f);}
 for(int s=0;s<5;s++)if(manual&(1u<<s)){scripted_loss(&in,s,false,0,scenario);in.damage[s].regeneration=0;in.damage[s].recovery=0;}
 return golem_sample(&in);
}
static Camera3D camera_for(GolemPose *p,int view){
 Vector3 target=Vector3Add(p->input.position,(Vector3){0,2.7f,0});
 Vector3 offset=view==1?(Vector3){0,.6f,-8}:view==2?(Vector3){-8,.6f,0}:view==3?(Vector3){0,.6f,8}:(Vector3){-5,2.0f,-7};
 offset=Vector3RotateByAxisAngle(Vector3Scale(offset,1.5f),(Vector3){0,1,0},orbit);
 return (Camera3D){Vector3Add(target,offset),target,{0,1,0},43,CAMERA_PERSPECTIVE};
}
static void scene(GolemPose *p,Camera3D cam,WeaponRenderPass pass){
 if(showModel)golem_draw(p,cam.position,pass,effects);
 if(!wall){
  CharacterPose rat=character_pose(Vector3Add(p->input.position,(Vector3){-3.0f,1,-1.05f}),25,0,0,true,false,SKYBLUE,p->input.time);
  if(pass==WEAPON_GLASS)character_draw_visor(rat);else{character_draw(rat,pass==WEAPON_EMISSION);if(pass==WEAPON_COLOR)character_draw_idle_arms(rat);}
 }
 if(skeleton&&pass==WEAPON_GLASS)golem_draw_skeleton(p);
}
static void render(RenderTexture2D rt,GolemPose *p,int view){
 Camera3D cam=camera_for(p,view);
 BeginTextureMode(rt);ClearBackground((Color){16,24,34,255});BeginMode3D(cam);
 DrawPlane((Vector3){0,-.014f,0},(Vector2){80,80},(Color){31,42,51,255});
 for(int i=-20;i<=20;i++){DrawLine3D((Vector3){i,0,-20},(Vector3){i,0,20},(Color){50,64,75,255});DrawLine3D((Vector3){-20,0,i},(Vector3){20,0,i},(Color){50,64,75,255});}
 if(wall)DrawCube(Vector3Lerp(cam.position,cam.target,.5f),12,12,.7f,GRAY);
 scene(p,cam,WEAPON_COLOR);EndMode3D();EndTextureMode();
 if(bloom&&grShaded){
  assert(weapon_bloom_resize(&fx,rt.texture.width,rt.texture.height));weapon_bloom_begin(&fx,rt);BeginMode3D(cam);
  scene(p,cam,WEAPON_EMISSION);EndMode3D();weapon_bloom_composite(&fx,rt);EndTextureMode();
 }
 BeginTextureMode(rt);BeginMode3D(cam);scene(p,cam,WEAPON_GLASS);EndMode3D();EndTextureMode();
}
static void save(RenderTexture2D rt,const char *out,const char *name){
 char file[1024];snprintf(file,sizeof(file),"%s/%s.png",out,name);Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);assert(ExportImage(image,file));UnloadImage(image);
}
static void geometry_checks(void){
 GolemStats stats=golem_renderer_stats();assert(stats.cubes<=384&&stats.triangles<=5000&&stats.models==18);
 unsigned vao=gm[GOLEM_TORSO].meshes[0].vaoId;assert(golem_renderer_init(NULL));assert(gm[GOLEM_TORSO].meshes[0].vaoId==vao);
 for(int k=1;k<GM_COUNT;k++){
  if(!gm[k].meshCount)continue; /* Physics-only shoulder mounts have no cached mesh. */
  Mesh m=gm[k].meshes[0];for(int i=0;i<m.vertexCount;i++){
   Vector3 n={m.normals[i*3],m.normals[i*3+1],m.normals[i*3+2]};assert(fabsf(Vector3Length(n)-1)<.0001f);
   for(int a=0;a<3;a++)assert(isfinite(m.vertices[i*3+a]));
  }
 }
 for(int scenario=0;scenario<10;scenario++)for(int f=0;f<600;f++){
  GolemPose p=sample(scenario,f/100.0f,0),same=sample(scenario,f/100.0f,0);assert(memcmp(&p,&same,sizeof(p))==0);
  for(int j=0;j<GOLEM_JOINT_COUNT;j++){float16 values=MatrixToFloatV(p.joint[j]);for(int a=0;a<16;a++)assert(isfinite(values.v[a]));}
  for(int s=0;s<2;s++){
   int arm=s?GOLEM_ARM_R:GOLEM_ARM_L;
   Vector3 end=Vector3Transform((Vector3){0,-.72f,0},p.joint[arm]);assert(Vector3Distance(end,gr_origin(p.joint[arm+1]))<.005f);
   end=Vector3Transform((Vector3){0,-.68f,0},p.joint[arm+1]);assert(Vector3Distance(end,gr_origin(p.joint[arm+2]))<.005f);
   int thigh=s?GOLEM_THIGH_R:GOLEM_THIGH_L;
   end=Vector3Transform((Vector3){0,-.82f,0},p.joint[thigh]);assert(Vector3Distance(end,gr_origin(p.joint[thigh+1]))<.005f);
   end=Vector3Transform((Vector3){0,-.85f,0},p.joint[thigh+1]);assert(Vector3Distance(end,gr_origin(p.joint[thigh+2]))<.005f);
   if(scenario==1){float phase=fmodf(p.input.travel/1.2f+s*.5f,1);if(phase<.6f){Vector3 ankle=gr_origin(p.joint[s?GOLEM_FOOT_R:GOLEM_FOOT_L]);assert(Vector3Distance(ankle,p.input.feet[s])<.001f);assert(fabsf(ankle.y-.33f)<.001f);}}
  }
 }
 for(int mask=0;mask<32;mask++)for(int action=0;action<=4;action++){
  GolemPose p=sample(action,2,mask);assert(p.kneel==((mask&12)?1:0));
  if(mask&12)assert(p.impact==0);
  if(mask&2)assert(p.impact==0);
  for(int s=0;s<2;s++)if(mask&(1u<<(s+2)))assert(p.knee[s].y>.30f && p.knee[s].y<.45f);
 }
 GolemPose slam=sample(2,.7f+2.2f*.63f,0);printf("Slam hand above floor: %.3f m\n",slam.hand[1].y);assert(slam.hand[1].y<.08f);
 GolemPose kneel=sample(6,2,0);printf("Kneel brace above floor: %.3f m\n",kneel.hand[0].y);assert(kneel.hand[0].y<.08f);
 assert(sample(5,2,0).visibility[GOLEM_ARM_R]==0);assert(sample(5,5.5f,0).visibility[GOLEM_ARM_R]==1);assert(sample(6,5.5f,0).kneel==0);
 printf("PASS: %d cubes, %d triangles, %d cached models; deterministic poses, joints, contacts, damage and recovery\n",stats.cubes,stats.triangles,stats.models);
}
static void occlusion_check(RenderTexture2D rt){
 GolemPose p=sample(5,3.2f,0);wall=true;showModel=false;render(rt,&p,1);Image a=LoadImageFromTexture(rt.texture);
 showModel=true;render(rt,&p,1);Image b=LoadImageFromTexture(rt.texture);
 Color *ca=LoadImageColors(a),*cb=LoadImageColors(b);int differing=0;
 for(int i=0;i<a.width*a.height;i++)if(abs(ca[i].r-cb[i].r)+abs(ca[i].g-cb[i].g)+abs(ca[i].b-cb[i].b)>3)differing++;
 printf("Occluded color/emission/ghost differing pixels: %d\n",differing);assert(differing==0);
 UnloadImageColors(ca);UnloadImageColors(cb);UnloadImage(a);UnloadImage(b);wall=false;
}
static int compare_double(const void *a,const void *b){double aa=*(const double*)a,bb=*(const double*)b;return (aa>bb)-(aa<bb);}
static void benchmark(const char *out){
 RenderTexture2D rt=LoadRenderTexture(640,360);char filename[1024];snprintf(filename,sizeof(filename),"%s/benchmark.csv",out);FILE *f=fopen(filename,"w");assert(f);fprintf(f,"run,effects,frame,ms\n");
 for(int run=0;run<3;run++)for(int active=0;active<2;active++){
  double timings[180];effects=active;
  for(int frame=-60;frame<180;frame++){
   GolemPose p=sample(8,3.0f+fmodf((frame+60)/60.0f,1),0);
   double start=GetTime();for(int view=0;view<4;view++)render(rt,&p,view);glFinish();double ms=(GetTime()-start)*1000;
   if(frame>=0){timings[frame]=ms;fprintf(f,"%d,%d,%d,%.5f\n",run,active,frame,ms);}
  }
  qsort(timings,180,sizeof(double),compare_double);printf("Four views run %d effects %d: median %.3f ms p95 %.3f ms\n",run,active,timings[90],timings[171]);
 }
 fclose(f);UnloadRenderTexture(rt);
}
int main(int argc,char **argv){
 setvbuf(stdout,NULL,_IONBF,0);
 bool automated=argc>1,capture=automated&&!strcmp(argv[1],"--capture"),bench=automated&&!strcmp(argv[1],"--benchmark"),fallback=automated&&!strcmp(argv[1],"--fallback");
 const char *out=argc>2?argv[2]:"artifacts/golem";SetTraceLogLevel(LOG_WARNING);if(automated)SetConfigFlags(FLAG_WINDOW_HIDDEN);
 InitWindow(1100,800,"Foundry golem - visual prototype");assert(weapon_renderer_init());assert(character_renderer_init());assert(golem_renderer_init(fallback?"missing-golem-shaders":NULL));assert(grShaded!=fallback);
 geometry_checks();RenderTexture2D rt=LoadRenderTexture(900,680);occlusion_check(rt);
 if(bench)benchmark(out);
 else if(capture){
  for(int view=0;view<4;view++){GolemPose p=sample(0,0,0);render(rt,&p,view);char name[64];snprintf(name,sizeof(name),"view-%d",view);save(rt,out,name);}
  skeleton=true;showModel=false;GolemPose rig=sample(0,0,0);render(rt,&rig,0);save(rt,out,"skeleton");skeleton=false;showModel=true;
  bloom=false;render(rt,&rig,0);save(rt,out,"no-bloom");bloom=true;
  RenderTexture2D four=LoadRenderTexture(1280,720),tile=LoadRenderTexture(640,360);
  for(int view=0;view<4;view++){
   GolemPose p=sample(view+5,3.1f,0);render(tile,&p,view);
   BeginTextureMode(four);DrawTextureRec(tile.texture,(Rectangle){0,0,640,-360},(Vector2){(view%2)*640,(view/2)*360},WHITE);EndTextureMode();
  }
  save(four,out,"four-view");UnloadRenderTexture(four);UnloadRenderTexture(tile);
  RenderTexture2D frames=LoadRenderTexture(720,544);
  for(int s=0;s<10;s++)for(int frame=0;frame<120;frame++){
   GolemPose p=sample(s,frame/20.0f,0);render(frames,&p,0);
   BeginTextureMode(frames);DrawText("FOUNDRY / VOXEL GOLEM",20,18,22,RAYWHITE);DrawText(TextFormat("%s   %.2fs",names[s],frame/20.0f),20,48,18,SKYBLUE);DrawText("PRESENTATION PROTOTYPE - RAT SHOWN FOR SCALE",20,frames.texture.height-27,14,LIGHTGRAY);EndTextureMode();
   char name[80];snprintf(name,sizeof(name),"%s-%03d",names[s],frame);save(frames,out,name);
  }
  UnloadRenderTexture(frames);
 }else if(automated){GolemPose p=sample(6,3.3f,0);render(rt,&p,0);save(rt,out,fallback?"fallback":"checked");}
 else{
  SetTargetFPS(60);float t=0;int scenario=0;bool paused=false;unsigned missing=0;
  while(!WindowShouldClose()){
   if(IsKeyPressed(KEY_SPACE))paused=!paused;
   if(IsKeyPressed(KEY_TAB)){scenario=(scenario+1)%10;t=0;}
   if(IsKeyPressed(KEY_R)){t=0;missing=0;}if(IsKeyPressed(KEY_S)){skeleton=!skeleton;showModel=!skeleton;}
   if(IsKeyPressed(KEY_B))bloom=!bloom;
   if(IsKeyPressed(KEY_E))effects=!effects;
   for(int i=0;i<5;i++)if(IsKeyPressed(KEY_ONE+i))missing^=1u<<i;
   if(IsKeyDown(KEY_LEFT))t=fmaxf(0,t-GetFrameTime());
   if(IsKeyDown(KEY_RIGHT))t=fminf(5.99f,t+GetFrameTime());
   if(IsKeyPressed(KEY_PERIOD)){paused=true;t=fminf(5.99f,t+.05f);}if(IsKeyPressed(KEY_COMMA)){paused=true;t=fmaxf(0,t-.05f);}
   if(IsMouseButtonDown(MOUSE_BUTTON_RIGHT))orbit+=GetMouseDelta().x*.006f;
   if(!paused)t=fmodf(t+GetFrameTime(),6);
   GolemPose p=sample(scenario,t,missing);render(rt,&p,0);
   BeginDrawing();ClearBackground((Color){12,18,25,255});DrawTextureRec(rt.texture,(Rectangle){0,0,900,-680},(Vector2){100,0},WHITE);
   DrawText(TextFormat("FOUNDRY GOLEM | %s | %.2fs %s",names[scenario],t,paused?"PAUSED":""),20,690,24,RAYWHITE);
   DrawText("Tab: scenario   Space: pause   arrows: scrub   , / .: step   RMB drag: orbit",20,728,18,LIGHTGRAY);
   DrawText("1/2: arm L/R   3/4: lower leg L/R   5: chest   R: reset   S: skeleton   B: bloom   E: effects",20,759,17,SKYBLUE);EndDrawing();
  }
 }
 UnloadRenderTexture(rt);weapon_bloom_unload(&fx);golem_renderer_shutdown();assert(!grReady&&grStats.models==0);golem_renderer_shutdown();assert(golem_renderer_init(NULL));golem_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();puts("PASS: occlusion and resource lifecycle");return 0;
}
