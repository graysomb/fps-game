/* Cached presentation geometry; story sampling never reads gameplay state. */
enum { SP_WHEEL,SP_SEAT,SP_SCOPE,SP_BREAD,SP_SCRAP,SP_COUNT };
static Model storyModels[SP_COUNT];
static float story_ease(float x){x=Clamp(x,0,1);return x*x*(3-2*x);}
static void story_box(TownMesh *b,int m,Vector3 p,Vector3 size,Color c){town_box(&b[m],p,size,c,NULL);}
static void town_story_geometry(TownMesh *b){
 Color steel={112,137,145,255},cream={219,203,161,255},blue={83,185,220,255};
 story_box(b,T_METAL,(Vector3){-22,.65f,-27.6f},(Vector3){1.8f,1.3f,.8f},steel);
 for(int k=0;k<6;k++)town_ring(&b[T_GLOW],(Vector3){-22.5f+k*.2f,.85f,-28.04f},.23f,.035f,blue);
 for(int k=0;k<3;k++)story_box(b,T_METAL,(Vector3){-21.45f,.5f+k*.25f,-28.03f},(Vector3){.12f,.12f,.06f},cream);
 town_tube(&b[T_METAL],(Vector3){-21.4f,.55f,-28},(Vector3){-20,.08f,-28.5f},.035f,.035f,steel);
 story_box(b,T_WOOD,(Vector3){-19.5f,1,-27.5f},(Vector3){1.6f,2,.25f},townWood);
 for(int r=0;r<3;r++){story_box(b,T_WOOD,(Vector3){-19.5f,.3f+r*.55f,-27.8f},(Vector3){1.7f,.08f,.6f},townWood);for(int k=0;k<9;k++)story_box(b,T_CLOTH,(Vector3){-20.15f+k*.16f,.52f+r*.55f,-27.8f},(Vector3){.12f,.35f+(k%3)*.04f,.3f},k%2?cream:steel);}
 story_box(b,T_WOOD,(Vector3){-18.5f,.35f,-29},(Vector3){.8f,.15f,.8f},townWood);story_box(b,T_WOOD,(Vector3){-18.5f,.7f,-28.65f},(Vector3){.8f,.8f,.12f},townWood);
 town_lantern(b,(Vector3){-19.5f,2.2f,-27.7f});
 for(int k=0;k<5;k++)story_box(b,T_GLOW,(Vector3){-22+k*.7f,2.2f,-27.9f},(Vector3){.10f,.15f,.08f},k<2?blue:townGold);
 for(int stall=0;stall<2;stall++){float z=-1.8f+stall*3.7f;
 story_box(b,T_WOOD,(Vector3){-22,.8f,z},(Vector3){2.5f,.12f,1.1f},townWood);
 for(int k=0;k<4;k++)story_box(b,T_WOOD,(Vector3){-23+(k%2)*2,1.2f,z+(k/2-.5f)*.9f},(Vector3){.09f,2.4f,.09f},townWood);
 for(int k=0;k<8;k++)story_box(b,T_CLOTH,(Vector3){-23.2f+k*.34f,2.35f,z},(Vector3){.34f,.06f,1.7f},k%2?cream:(Color){97,139,132,255});
 for(int k=0;k<5;k++)town_bush(&b[T_WOOD],(Vector3){-22.85f+k*.4f,.98f,z},(Vector3){.3f,.22f,.45f},(Color){190,133,64,255});
 town_crate(b,(Vector3){-23,.25f,z+1},.7f);town_lantern(b,(Vector3){-23,1.7f,z});}
 town_tube(&b[T_METAL],(Vector3){-21.3f,.85f,-1.7f},(Vector3){-21.3f,1.4f,-1.7f},.04f,.04f,steel);
 town_tube(&b[T_METAL],(Vector3){-21.7f,1.4f,-1.7f},(Vector3){-20.9f,1.4f,-1.7f},.025f,.025f,steel);
 for(int k=0;k<2;k++){float x=-21.7f+k*.8f;town_tube(&b[T_METAL],(Vector3){x,1.4f,-1.7f},(Vector3){x,1.05f,-1.7f},.01f,.01f,cream);story_box(b,T_METAL,(Vector3){x,1.05f,-1.7f},(Vector3){.35f,.03f,.35f},steel);}
 /* Hollow robot head, steps and a salvaged arm supporting the swing. */
 story_box(b,T_METAL,(Vector3){33,1.4f,-32},(Vector3){3,2.8f,.3f},steel);
 for(int side=-1;side<=1;side+=2)story_box(b,T_METAL,(Vector3){33+side*1.35f,1.2f,-31},(Vector3){.3f,2.4f,2},steel);
 story_box(b,T_METAL,(Vector3){33,2.5f,-31},(Vector3){3,.35f,2.4f},steel);
 for(int k=0;k<4;k++)story_box(b,T_WOOD,(Vector3){34.8f,(.3f+k*.28f)/2,-29.2f-k*.4f},(Vector3){.8f,.3f+k*.28f,.4f},townWood);
 /* Recognizable face: one hollow eye for entry, the other a bolted lens. */
 for(int k=0;k<2;k++){float x=32.3f+k*1.4f;town_ring(&b[T_METAL],(Vector3){x,1.5f,-29.76f},.57f,.12f,cream);}
 story_box(b,T_METAL,(Vector3){33,2.1f,-29.8f},(Vector3){2.8f,.3f,.28f},steel);
 story_box(b,T_METAL,(Vector3){33,.35f,-29.8f},(Vector3){2.8f,.25f,.28f},steel);
 for(int k=0;k<6;k++)story_box(b,T_METAL,(Vector3){32.2f+k*.32f,.4f,-29.62f},(Vector3){.16f,.12f,.08f},cream);
 town_disc(&b[T_GLASS],(Vector3){33.7f,1.5f,-29.78f},.46f,.02f,(Color){83,130,143,140});
 town_tube(&b[T_METAL],(Vector3){33,2.6f,-31},(Vector3){33.2f,3.1f,-31},.04f,.04f,steel);
 town_bush(&b[T_METAL],(Vector3){33.2f,3.15f,-31},(Vector3){.18f,.18f,.18f},cream);
 town_tube(&b[T_METAL],(Vector3){30,0,-31},(Vector3){30,3,-31},.17f,.17f,steel);town_tube(&b[T_METAL],(Vector3){30,3,-31},(Vector3){32,3,-31},.18f,.12f,steel);
 for(int k=0;k<18;k++)town_plant(&b[T_LEAF],(Vector3){29+k*.36f,.03f,-32.7f},k,.65f,true);
 for(int k=0;k<3;k++)town_tube(&b[T_METAL],(Vector3){21,6.9f,21},(Vector3){21+cosf(k*2*PI/3)*.7f,5.75f,21+sinf(k*2*PI/3)*.7f},.055f,.08f,steel);
 story_box(b,T_WOOD,(Vector3){23,6.35f,21},(Vector3){1.4f,.12f,1},townWood);story_box(b,T_CLOTH,(Vector3){23,6.42f,21},(Vector3){1.2f,.015f,.8f},cream);
 for(int k=0;k<10;k++)story_box(b,T_METAL,(Vector3){22.5f+(k%4)*.28f,6.435f,20.7f+(k/4)*.25f},(Vector3){.035f,.01f,.035f},blue);
 for(int k=0;k<5;k++)town_tube(&b[T_WOOD],(Vector3){18+k*1.4f,5.75f,24},(Vector3){18+k*1.4f,6.6f,24},.055f,.055f,townWood);
 town_tube(&b[T_WOOD],(Vector3){18,6.6f,24},(Vector3){24,6.6f,24},.055f,.055f,townWood);
 for(int k=0;k<2;k++)story_box(b,T_CLOTH,(Vector3){20+k*2,5.9f,22.7f},(Vector3){.8f,.25f,.65f},k?cream:blue);
 town_lantern(b,(Vector3){23.5f,6.55f,21});
 TownMesh p[SP_COUNT]={0};
 town_ring(&p[SP_WHEEL],(Vector3){0},.42f,.055f,cream);for(int k=0;k<6;k++)town_tube(&p[SP_WHEEL],(Vector3){0},(Vector3){cosf(k*PI/3)*.4f,sinf(k*PI/3)*.4f,0},.025f,.025f,steel);
 town_box(&p[SP_SEAT],(Vector3){0},(Vector3){.75f,.09f,.42f},townWood,NULL);
 town_tube(&p[SP_SCOPE],(Vector3){0,0,-.55f},(Vector3){0,.4f,.75f},.15f,.3f,steel);town_tube(&p[SP_SCOPE],(Vector3){0,.4f,.75f},(Vector3){0,.42f,.80f},.29f,.29f,blue);
 for(int k=0;k<3;k++)town_tube(&p[SP_SCOPE],(Vector3){0,.12f*k,-.15f+.39f*k},(Vector3){0,.12f*k+.04f,-.02f+.39f*k},.2f+k*.035f,.2f+k*.035f,cream);
 town_tube(&p[SP_SCOPE],(Vector3){-.24f,.05f,0},(Vector3){-.33f,.05f,0},.08f,.08f,townWood);
 for(int k=0;k<4;k++)story_box(b,T_WOOD,(Vector3){22.5f+(k%2),6.03f,20.65f+(k/2)*.7f},(Vector3){.07f,.56f,.07f},townWood);
 town_bush(&p[SP_BREAD],(Vector3){0},(Vector3){.35f,.22f,.5f},(Color){202,145,70,255});town_box(&p[SP_SCRAP],(Vector3){0},(Vector3){.18f,.14f,.18f},steel,NULL);
 for(int k=0;k<SP_COUNT;k++)storyModels[k]=town_upload(&p[k]);
}
static StoryActor *story_actor(TownStoryFrame *s,int role,Vector3 floor,float yaw,float scale,float t){
 assert(s->count<4);StoryActor *a=&s->actors[s->count++];*a=(StoryActor){.role=role,.floor=floor,.scale=scale};
 a->pose=character_pose((Vector3){0,1,0},yaw,0,0,true,false,role%2?GOLD:SKYBLUE,t);
 a->pose.breath=sinf(t*PI/2)*.004f;a->pose.time=sin(t*PI/4);
 a->pose.feetResolved=true;
 Matrix rot=MatrixRotateY(yaw*DEG2RAD);
 for(int j=0;j<2;j++)a->pose.soles[j]=Vector3Transform((Vector3){j?.19f:-.19f,0,-.08f},rot);
 Matrix body=character_body_transform(a->pose);
 for(int j=0;j<2;j++){float side=j?1:-1;a->shoulders[j]=Vector3Transform((Vector3){side*.17f,.11f,0},body);a->hands[j]=Vector3Transform((Vector3){side*.16f,.08f,-.3f},body);}return a;
}
static void story_hand(StoryActor *a,int j,Vector3 world){a->hands[j]=Vector3Scale(Vector3Subtract(world,a->floor),1/a->scale);}
static Vector3 story_world_hand(StoryActor *a,int j){return Vector3Add(a->floor,Vector3Scale(a->hands[j],a->scale));}
static TownStoryFrame town_story_sample(int scene,float t,Camera3D camera){
 TownStoryFrame s={.power=1};float phase=fmodf(t,8)*PI/4;
 if(scene==SC_POWER){s.power=story_ease((t-2.5f)/2);s.wheel=360*story_ease(t/3);
 StoryActor *a=story_actor(&s,R_WHEEL,(Vector3){-22.8f,0,-28.65f},180,1,t);float angle=s.wheel*DEG2RAD;story_hand(a,0,(Vector3){-22.7f+.4f*cosf(angle),.9f+.4f*sinf(angle),-28.1f});
 a=story_actor(&s,R_CABLE,(Vector3){-20.8f,0,-28.65f},180,1,t);story_hand(a,1,Vector3Lerp((Vector3){-20.9f,.7f,-28.7f},(Vector3){-21.45f,.75f,-28.1f},story_ease(t/2.5f)));
 a=story_actor(&s,R_READER,(Vector3){-18.5f,.42f,-29},0,1,t);a->pose.pitch=12;

 }else if(scene==SC_MARKET){story_actor(&s,R_VENDOR,(Vector3){-22,0,-.6f},0,1,t);
 float x=-19.5f-story_ease(t/2)*1.3f+story_ease((t-5)/2)*2;StoryActor *a=story_actor(&s,R_BUYER,(Vector3){x,0,-2.9f},180,1,t);
 s.scrapOwner=t<2?R_BUYER:-1;s.breadOwner=t<4?R_VENDOR:R_BUYER;Vector3 tray={-20.9f,1.16f,-1.7f},exchange={-21.4f,1.1f,-1.9f};
 s.scrap=t<2?Vector3Lerp(story_world_hand(a,0),tray,story_ease((t-1.3f)/.7f)):tray;
 if(t<2)story_hand(a,0,s.scrap);
 Vector3 buyerHand=story_world_hand(a,1);
 story_hand(a,1,Vector3Lerp(buyerHand,exchange,story_ease((t-3)/1)*(1-story_ease((t-4)/1))));
 s.bread=t<4?Vector3Lerp(story_world_hand(&s.actors[0],1),(Vector3){-21.4f,1.1f,-1.9f},story_ease((t-2)/2)):story_world_hand(a,1);if(t>=2&&t<4)story_hand(&s.actors[0],1,s.bread);
 }else if(scene==SC_PLAYGROUND){s.swingAngle=.35f*sinf(phase);s.swingSeat=(Vector3){31,3-2.1f*cosf(s.swingAngle),-31+2.1f*sinf(s.swingAngle)};
 StoryActor *a=story_actor(&s,R_SWING,Vector3Add(s.swingSeat,(Vector3){0,.04f,0}),180,.65f,t);for(int j=0;j<2;j++)story_hand(a,j,Vector3Add(s.swingSeat,(Vector3){j?.31f:-.31f,.7f,0}));
 float step=(.5f-.5f*cosf(phase))*3;
 int stair=(int)floorf(step);float f=step-stair;
 a=story_actor(&s,R_CLIMB,(Vector3){34.8f,.3f+.28f*(stair+story_ease(f)),-29.2f-step*.4f},180,.65f,t);
 for(int j=0;j<2;j++){float y=(j==0?ceilf(step):floorf(step))*.28f+.3f;a->pose.soles[j].y=(y-a->floor.y)/a->scale;}
 story_actor(&s,R_PEEK,(Vector3){32.3f,.75f,-30.1f+.2f*sinf(phase)},180,.65f,t);
 }else if(scene==SC_OBSERVATORY){s.telescopeAngle=8*sinf(phase);StoryActor *a=story_actor(&s,R_ASTRONOMER,(Vector3){21,5.75f,20},180,1,t);story_hand(a,1,(Vector3){21.15f,6.9f,20.6f});story_actor(&s,R_CHART,(Vector3){23,5.75f,20},180,1,t);
 float meteor=Clamp((t-5)/1.2f,0,1);s.meteorGlow=t>5&&t<6.2f?sinf(meteor*PI):0;
 Vector3 f=Vector3Normalize(Vector3Subtract(camera.target,camera.position)),r=Vector3Normalize(Vector3CrossProduct(f,camera.up)),u=Vector3CrossProduct(r,f);
 Vector3 head=Vector3Normalize(Vector3Add(f,Vector3Add(Vector3Scale(u,.35f),Vector3Scale(r,-.3f+meteor*.6f))));Vector3 tail=Vector3Normalize(Vector3Add(head,Vector3Scale(r,-.10f)));
 for(int k=0;k<s.count;k++)s.actors[k].pose.pitch=25*story_ease((t-4.8f)/.5f)*(1-story_ease((t-6.5f)/1.5f));
 s.meteor=(Vector4){head.x/(1+fmaxf(head.y,0)),head.z/(1+fmaxf(head.y,0)),tail.x/(1+fmaxf(tail.y,0)),tail.z/(1+fmaxf(tail.y,0))};
 }return s;
}
static void town_story_draw(TownStoryFrame s,int scene,float t,bool emission){
 if(!emission){int mode=T_METAL;SetShaderValue(town.shader,GetShaderLocation(town.shader,"townMode"),&mode,SHADER_UNIFORM_INT);if(scene==SC_POWER)DrawModelEx(storyModels[SP_WHEEL],(Vector3){-22.7f,.9f,-28.1f},(Vector3){0,0,1},s.wheel,(Vector3){1,1,1},WHITE);
 if(scene==SC_MARKET){DrawModel(storyModels[SP_BREAD],s.bread,1,WHITE);DrawModel(storyModels[SP_SCRAP],s.scrap,1,WHITE);}
 if(scene==SC_PLAYGROUND){DrawModel(storyModels[SP_SEAT],s.swingSeat,1,WHITE);for(int k=0;k<2;k++)DrawCylinderEx((Vector3){31+(k?.32f:-.32f),3,-31},Vector3Add(s.swingSeat,(Vector3){k?.32f:-.32f,0,0}),.018f,.018f,6,(Color){184,165,119,255});}
 if(scene==SC_OBSERVATORY)DrawModelEx(storyModels[SP_SCOPE],(Vector3){21,6.9f,21},(Vector3){0,1,0},s.telescopeAngle,(Vector3){1,1,1},WHITE);}
 for(int k=0;k<s.count;k++){StoryActor a=s.actors[k];if(!emission){rlDrawRenderBatchActive();rlDisableDepthMask();DrawCylinder(Vector3Add(a.floor,(Vector3){0,.01f,0}),.32f*a.scale,.32f*a.scale,.001f,16,(Color){12,17,20,35});rlDrawRenderBatchActive();rlEnableDepthMask();}
 rlPushMatrix();rlTranslatef(a.floor.x,a.floor.y,a.floor.z);rlScalef(a.scale,a.scale,a.scale);character_presentation_light((Vector3){1,.8f,.6f},.3f);character_draw(a.pose,emission);
 if(!emission){for(int j=0;j<2;j++)character_draw_arm(a.shoulders[j],a.hands[j],a.role%2?GOLD:SKYBLUE);if(a.role==R_READER){Vector3 p=Vector3Scale(Vector3Add(a.hands[0],a.hands[1]),.5f);DrawCube(p,.38f,.04f,.28f,(Color){217,202,162,255});}character_draw_visor(a.pose);}rlPopMatrix();}
 character_presentation_light((Vector3){1,1,1},0);
}
static void town_story_shutdown(void){for(int k=0;k<SP_COUNT;k++){storyModels[k].materials[0].shader=(Shader){.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};UnloadModel(storyModels[k]);storyModels[k]=(Model){0};}}
