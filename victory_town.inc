/* Presentation-only Foundry town. Included by the standalone preview harness.
 * All geometry is baked at init; samples never inspect or mutate simulation. */
typedef struct {float *v,*n;unsigned char *c;int count,capacity;} TownMesh;
enum TownMaterial { T_STONE,T_LEAF,T_GLOW,T_WATER,T_WOOD,T_METAL,T_CLOTH,T_GLASS,T_SOIL,T_COUNT };
enum TownActivity { A_EAT,A_IDLE,A_READ,A_WORK,A_CARRY,A_WATER,A_LIFT };
typedef struct {Model models[T_COUNT];Shader shader,depthShader;RenderTexture2D shadow;Matrix lightMatrix;bool ready;int triangles;} VictoryTown;
typedef struct {float time,phase;int scene;Camera3D camera;} TownSample;
static VictoryTown town;
static Vector3 townHands[16][2];
static CharacterPose townVisors[16];static int townVisorCount;
static const Color townWood={99,63,38,255},townLeaf={76,131,70,255},townGold={255,189,76,255};
static void town_tri(TownMesh *b,Vector3 a,Vector3 c,Vector3 d,Color color){
 if(b->count+3>b->capacity){b->capacity=b->capacity?b->capacity*2:4096;b->v=realloc(b->v,b->capacity*3*sizeof(float));b->n=realloc(b->n,b->capacity*3*sizeof(float));b->c=realloc(b->c,b->capacity*4);assert(b->v&&b->n&&b->c);}
 Vector3 normal=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(c,a),Vector3Subtract(d,a))),vs[3]={a,c,d};
 for(int j=0;j<3;j++){int i=b->count++;memcpy(b->v+i*3,&vs[j],12);memcpy(b->n+i*3,&normal,12);memcpy(b->c+i*4,&color,4);}
}
static void town_box(TownMesh *b,Vector3 p,Vector3 size,Color color,const bool *faces){
 Vector3 v[8];for(int i=0;i<8;i++)v[i]=(Vector3){p.x+((i&1)?1:-1)*size.x/2,p.y+((i&2)?1:-1)*size.y/2,p.z+((i&4)?1:-1)*size.z/2};
 const int q[6][4]={{1,3,7,5},{4,6,2,0},{2,6,7,3},{4,0,1,5},{5,7,6,4},{0,2,3,1}};
 for(int f=0;f<6;f++)if(!faces||faces[f]){Color c=color;float shade=.94f+.02f*(f%4);c.r*=shade;c.g*=shade;c.b*=shade;town_tri(b,v[q[f][0]],v[q[f][1]],v[q[f][2]],c);town_tri(b,v[q[f][0]],v[q[f][2]],v[q[f][3]],c);}
}
static void town_bush(TownMesh *b,Vector3 p,Vector3 size,Color c){
 Vector3 top={p.x,p.y+size.y*.6f,p.z},bottom={p.x,p.y-size.y*.3f,p.z};
 for(int k=0;k<7;k++){float a=k*2*PI/7,z=(k+1)*2*PI/7;Vector3 v={p.x+cosf(a)*size.x*.5f,p.y,p.z+sinf(a)*size.z*.5f},w={p.x+cosf(z)*size.x*.5f,p.y,p.z+sinf(z)*size.z*.5f};Color shade=c;shade.g-=k*3;town_tri(b,top,w,v,shade);town_tri(b,bottom,v,w,shade);}
}
static void town_disc(TownMesh *b,Vector3 p,float r,float depth,Color c){
 for(int i=0;i<16;i++){float a=i*2*PI/16,d=(i+1)*2*PI/16;Vector3 x={p.x+r*cosf(a),p.y+r*sinf(a),p.z},y={p.x+r*cosf(d),p.y+r*sinf(d),p.z};town_tri(b,p,y,x,c);town_box(b,(Vector3){x.x,x.y,p.z+depth/2},(Vector3){.10f,.10f,depth},townWood,NULL);}
}
static Model town_upload(TownMesh *b){Mesh mesh={0};mesh.vertices=b->v;mesh.normals=b->n;mesh.colors=b->c;mesh.vertexCount=b->count;mesh.triangleCount=b->count/3;UploadMesh(&mesh,false);town.triangles+=mesh.triangleCount;return LoadModelFromMesh(mesh);}
#include "victory_town_detail.inc"
static void victory_town_init(void){
 assert(!town.ready);TownMesh b[T_COUNT]={0};
 /* Snapshot exposed authored surfaces, then release the simulation world. */
 rebuild_all_voxel_surfaces();
 for(int i=0;i<voxel_count;i++){Voxel *v=&voxels[i];if(v->simulate)continue;Color c=v->color;c.a=254;if(c.r>c.g*1.25f)c=(Color){126,119,96,254};
 Vector3 size={(v->max_gx-v->min_gx+1)*VOXEL_SIZE,(v->max_gy-v->min_gy+1)*VOXEL_SIZE,(v->max_gz-v->min_gz+1)*VOXEL_SIZE};
 town_box(&b[0],v->pos,size,c,v->surface);}
 town_box(&b[0],(Vector3){0,-.12f,0},(Vector3){600,.2f,600},(Color){58,76,67,255},NULL);
 for(int side=-1;side<=1;side+=2)for(int row=-1;row<=1;row+=2)town_home(b,side*21.5f,row*20.6f-5.35f,side+row*3);
 /* Low dining tables and stools remain on top of the original platform. */
 for(int k=0;k<4;k++){float a=k*PI/2;Vector3 p={sinf(a)*2.1f,2.0f,cosf(a)*2.1f};town_box(&b[0],p,(Vector3){.75f,.5f,.75f},townWood,NULL);}
 town_box(&b[0],(Vector3){0,2.15f,0},(Vector3){2.2f,.16f,1.1f},(Color){130,125,112,255},NULL);
 /* Garden beds and a stone-edged irrigation pool beside the east route. */
 for(int row=0;row<3;row++){float z=7+row*1.3f;town_box(&b[0],(Vector3){15,.18f,z},(Vector3){4,.35f,.9f},townWood,NULL);town_box(&b[0],(Vector3){15,.38f,z},(Vector3){3.8f,.06f,.7f},(Color){48,39,30,255},NULL);
 for(int k=0;k<10;k++)for(int l=0;l<3;l++)town_leaf(&b[T_LEAF],(Vector3){13.4f+k*.35f,.54f+l*.07f,z},k+l*2.4f,.22f,townLeaf);}
 town_box(&b[3],(Vector3){19,.025f,9},(Vector3){2.8f,.02f,4},(Color){40,103,112,255},NULL);
 for(int k=0;k<10;k++)for(int s=-1;s<=1;s+=2)town_box(&b[0],(Vector3){19+s*1.6f,.10f,7+k*.44f},(Vector3){.4f,.25f,.45f},(Color){112,126,123,255},NULL);
 town_box(&b[0],(Vector3){12,.8f,12},(Vector3){2.4f,.16f,1.1f},townWood,NULL);
 for(int s=-1;s<=1;s+=2)town_box(&b[0],(Vector3){12+s,.4f,12},(Vector3){.14f,.8f,.8f},townWood,NULL);
 /* Laundry and planted paths add domestic scale without changing routes. */
 for(int k=0;k<6;k++)town_box(&b[T_CLOTH],(Vector3){17+k*.65f,2.4f,14},(Vector3){.5f,.7f,.04f},k%2?(Color){143,181,188,255}:(Color){219,202,167,255},NULL);
 for(int k=0;k<32;k++){float z=-14+k;town_box(&b[0],(Vector3){11,.025f,z},(Vector3){.9f,.06f,.7f},(Color){127,137,124,255},NULL);if(k%3==0){town_plant(&b[T_LEAF],(Vector3){12,.12f,z},k,1,true);town_box(&b[2],(Vector3){10,.28f,z},(Vector3){.12f,.12f,.12f},townGold,NULL);}}
 town_details(b);
 town.shader=LoadShader("shaders/victory_town.vert","shaders/victory_town.frag");
 for(int i=0;i<T_COUNT;i++){town.models[i]=town_upload(&b[i]);if(town.shader.id!=rlGetShaderIdDefault())town.models[i].materials[0].shader=town.shader;}
 town_bake_shadow();
 town.ready=true;printf("Town cached triangles: %d\n",town.triangles);
}
static TownSample victory_town_sample(int scene,float time){
 TownSample s={.scene=scene,.time=time,.phase=fmodf(time,8)/8};float a=sinf(s.phase*2*PI);
 Camera3D cams[4]={{{4.1f,3.6f,-4.8f},{0,2.65f,0},{0,1,0},48,0},{{23.5f,3.5f,10.8f},{21.2f,1.25f,14.8f},{0,1,0},63,0},{{21.5f,5.4f,13.3f},{15.7f,.8f,9.5f},{0,1,0},60,0},{{39,17,-38},{0,6,0},{0,1,0},57,0}};
 s.camera=cams[scene];s.camera.position.x+=a*(scene==3?3:.16f);s.camera.position.y+=(scene==3?1.2f:0)*(1-cosf(s.phase*2*PI));return s;
}
static Vector3 town_carry_floor(float time){float a=time*PI/4;return (Vector3){19+1.2f*sinf(a),0,12.7f+.2f*cosf(a)};}
static float town_carry_yaw(float time){float a=time*PI/4;return atan2f(-1.2f*cosf(a),.2f*sinf(a))*RAD2DEG;}
static Vector3 town_carry_sole(float time,int foot){float start=floorf(time+foot*.5f)-foot*.5f,phase=time-start;float blend=Clamp((phase-.6f)/.4f,0,1);blend=blend*blend*(3-2*blend);Vector3 points[2];for(int k=0;k<2;k++){float at=start+.25f+k;points[k]=Vector3Add(town_carry_floor(at),Vector3Transform((Vector3){(foot?1:-1)*.19f,0,0},MatrixRotateY(town_carry_yaw(at)*DEG2RAD)));}Vector3 sole=Vector3Lerp(points[0],points[1],blend);sole.y=.065f*sinf(blend*PI);return sole;}
static void town_mouse(Vector3 floor,float yaw,Color gear,int activity,float t,bool emission){
 float phase=t*PI/4,offset=floor.x*.63f+floor.z*.42f;float cycle=sinf(phase*2+offset),walk=activity==A_CARRY?1:0;
 character_presentation_light((Vector3){1.12f,1.02f,.9f},1);CharacterPose p=character_pose(Vector3Add(floor,(Vector3){0,activity==0?.82f:1,0}),yaw,activity==2?18:0,walk,true,false,gear,t);
 p.breath=.008f*sinf(phase+offset);p.stride=walk*sinf(phase*8);p.time=sinf(phase)*.3f;
 p.feetResolved=true;Matrix rotation=MatrixRotateY(yaw*DEG2RAD);
 for(int j=0;j<2;j++){Vector3 local={(j?1:-1)*.19f,0,activity==0?-.3f:0};local.z+=walk*.12f*sinf(phase*8+j*PI);p.soles[j]=activity==A_CARRY?town_carry_sole(t,j):Vector3Add(floor,Vector3Transform(local,rotation));}
 if(!emission){rlDrawRenderBatchActive();rlDisableDepthMask();for(int k=0;k<4;k++)DrawCylinder(Vector3Add(floor,(Vector3){0,.008f+k*.001f,0}),.38f-k*.045f,.38f-k*.045f,.001f,20,(Color){12,17,20,13});rlDrawRenderBatchActive();rlEnableDepthMask();}
 character_draw(p,emission);if(emission)return;
 Matrix body=character_body_transform(p);
 for(int j=0;j<2;j++){float side=j?1:-1;Vector3 shoulder=Vector3Transform((Vector3){side*.17f,.11f,0},body),hand={side*.2f,.04f,-.25f};
 if(activity==A_EAT){hand.y=j?.10f+.09f*(.5f+.5f*cycle):.07f;hand.z=j?-.24f+.08f*cycle:-.26f;}
 if(activity==A_READ){hand.y=.10f;hand.x=side*.12f;}
 if(activity==A_WORK){hand.y+=.09f*cycle;hand.z-=.05f*cycle;}
 if(activity==A_WATER){hand.y=.04f;hand.z=-.31f;}
 if(activity==A_LIFT){hand.y=.09f;hand.z=-.32f;}
 if(activity==A_CARRY){hand.x=side*.13f;hand.y=.07f;hand.z=-.28f;}
 hand=Vector3Transform(hand,body);townHands[townVisorCount][j]=hand;character_draw_arm(shoulder,hand,gear);
 if(j==0&&activity==A_EAT){DrawCylinder(hand,.12f,.08f,.09f,10,(Color){169,125,75,255});DrawCylinder(Vector3Add(hand,(Vector3){0,.05f,0}),.105f,.105f,.008f,10,(Color){98,134,43,255});}
 if(j==1&&activity==A_EAT){Vector3 tip=Vector3Add(hand,Vector3Transform((Vector3){0,0,-.16f},rotation));DrawCylinderEx(hand,tip,.012f,.012f,5,(Color){170,150,105,255});DrawSphere(tip,.035f,(Color){181,166,128,255});}
 if(j==1&&activity==A_CARRY){Vector3 center=Vector3Transform((Vector3){0,.04f,-.3f},body);rlPushMatrix();rlTranslatef(center.x,center.y,center.z);rlRotatef(yaw,0,1,0);DrawCube((Vector3){0},.43f,.24f,.3f,townWood);for(int q=0;q<4;q++)DrawSphere((Vector3){(q%2-.5f)*.16f,.14f,(q/2-.5f)*.12f},.075f,q%2?ORANGE:GREEN);rlPopMatrix();}
 if(j==1&&activity==A_WORK){Vector3 tip=Vector3Add(hand,(Vector3){0,-.22f,0});DrawCylinderEx(hand,tip,.02f,.02f,6,townWood);DrawCube(tip,.16f,.06f,.08f,GRAY);}
 if(j==1&&activity==A_WATER){Vector3 p=Vector3Add(hand,(Vector3){0,-.07f,0});DrawCylinder(p,.10f,.12f,.19f,10,(Color){98,139,128,255});Vector3 spout=Vector3Add(p,Vector3Transform((Vector3){0,.04f,-.27f},rotation));DrawCylinderEx(p,spout,.03f,.018f,6,GRAY);for(int q=0;q<5;q++)DrawSphere(Vector3Add(spout,(Vector3){0,-.08f-q*.06f+.025f*cycle,0}),.012f,(Color){117,188,209,180});}
 if(j==1&&activity==A_READ){Vector3 center=Vector3Transform((Vector3){0,.09f,-.27f},body);rlPushMatrix();rlTranslatef(center.x,center.y,center.z);rlRotatef(yaw,0,1,0);DrawCube((Vector3){0},.38f,.035f,.26f,(Color){106,69,42,255});DrawCube((Vector3){0,.025f,0},.35f,.015f,.23f,(Color){216,204,166,255});rlPushMatrix();rlRotatef(70*sinf(phase),0,0,1);DrawCube((Vector3){.085f,.04f,0},.17f,.006f,.23f,(Color){225,217,189,255});rlPopMatrix();rlPopMatrix();}

 }
 assert(townVisorCount<16);townVisors[townVisorCount++]=p;
}
static void victory_town_draw(TownSample s,bool emission){
 if(!emission)townVisorCount=0;
 int pass=emission;bool shader=town.shader.id&&town.shader.id!=rlGetShaderIdDefault();
 if(shader){SetShaderValueMatrix(town.shader,GetShaderLocation(town.shader,"lightMatrix"),town.lightMatrix);SetShaderValueTexture(town.shader,GetShaderLocation(town.shader,"shadowMap"),town.shadow.texture);SetShaderValue(town.shader,GetShaderLocation(town.shader,"townTime"),&s.time,SHADER_UNIFORM_FLOAT);SetShaderValue(town.shader,GetShaderLocation(town.shader,"eye"),&s.camera.position,SHADER_UNIFORM_VEC3);SetShaderValue(town.shader,GetShaderLocation(town.shader,"emissionPass"),&pass,SHADER_UNIFORM_INT);}
 for(int mode=0;mode<T_COUNT;mode++){if(emission&&mode!=2)continue;if(shader)SetShaderValue(town.shader,GetShaderLocation(town.shader,"townMode"),&mode,SHADER_UNIFORM_INT);if(mode==T_GLASS){rlDrawRenderBatchActive();rlDisableDepthMask();}DrawModel(town.models[mode],(Vector3){0},1,WHITE);if(mode==T_GLASS){rlDrawRenderBatchActive();rlEnableDepthMask();}}
 Color colors[4]={SKYBLUE,RED,GREEN,GOLD};
 for(int k=0;k<4;k++){float a=k*PI/2;town_mouse((Vector3){sinf(a)*2.1f,2.25f,cosf(a)*2.1f},k*90,colors[k],A_EAT,s.time,emission);}
 town_mouse((Vector3){23,.19f,14},-25,SKYBLUE,A_READ,s.time,emission);
 town_mouse(town_carry_floor(s.time),town_carry_yaw(s.time),GOLD,A_CARRY,s.time,emission);
 town_mouse((Vector3){14,0,6.4f},180,GREEN,A_WATER,s.time,emission);
 town_mouse((Vector3){12,0,11},180,RED,A_WORK,s.time,emission);
 for(int k=0;k<2;k++)town_mouse((Vector3){8+k*3,0,12},k?90:-90,k?SKYBLUE:GOLD,A_LIFT,s.time,emission);
 town_mouse((Vector3){-21,5.75f,-20},0,SKYBLUE,A_READ,s.time,emission);
 if(emission)camp_draw_at(s.camera,(Vector3){0,2.24f,0},s.time,true);
 if(!emission){DuckPose duck={.position={19+.55f*cosf(s.time*PI/4),.04f,9+sinf(s.time*PI/4)},.yaw=-s.time*45,.scale=.65f,.time=s.time,.speed=.5f};duck_draw_floor(duck,true);duck_draw(duck);
 Vector3 block={9.5f,1.4f+.2f*sinf(s.time*PI/2),12};DrawCube(block,.6f,.6f,.6f,(Color){140,159,165,255});weapon_draw_drain_state(townHands[8][1],block,s.time,true);weapon_draw_drain_state(townHands[9][1],block,s.time,true);
 for(int i=0;i<townVisorCount;i++)for(int j=i+1;j<townVisorCount;j++)if(Vector3DistanceSqr(townVisors[j].position,s.camera.position)>Vector3DistanceSqr(townVisors[i].position,s.camera.position)){CharacterPose p=townVisors[i];townVisors[i]=townVisors[j];townVisors[j]=p;}
 for(int i=0;i<townVisorCount;i++)character_draw_visor(townVisors[i]);
 character_presentation_light((Vector3){1,1,1},0);
 camp_draw_at(s.camera,(Vector3){0,2.24f,0},s.time,false);
 }
 character_presentation_light((Vector3){1,1,1},0);
}
static void victory_town_shutdown(void){
 if(!town.ready)return;Shader fallback={.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};for(int i=0;i<T_COUNT;i++){town.models[i].materials[0].shader=fallback;UnloadModel(town.models[i]);}if(town.shadow.id)UnloadRenderTexture(town.shadow);if(town.depthShader.id&&town.depthShader.id!=fallback.id)UnloadShader(town.depthShader);if(town.shader.id&&town.shader.id!=fallback.id)UnloadShader(town.shader);town=(VictoryTown){0};
}

static void town_sky(Camera3D camera,float time) {
 const Color horizon={19,31,43,255};
    // One sky per viewport, before periodic world images. Direction only: no
    // parallax, tile repetition or discontinuity when the player wraps.
    if(!nightSkyShader.id)nightSkyShader=LoadShader(NULL,"shaders/victory_sky.frag");
    if(nightSkyShader.id!=rlGetShaderIdDefault() && nightSkyShader.id) {
        int w=rlGetFramebufferWidth(),h=rlGetFramebufferHeight();
        Vector3 forward=Vector3Normalize(Vector3Subtract(camera.target,camera.position));
        Vector3 right=Vector3Normalize(Vector3CrossProduct(forward,camera.up));
        Vector3 up=Vector3CrossProduct(right,forward);
        Vector4 viewport={(float)w,(float)h,tanf(camera.fovy*DEG2RAD*.5f),(float)w/h};
        Vector3 haze={horizon.r/255.f,horizon.g/255.f,horizon.b/255.f};
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewForward"),&forward,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewRight"),&right,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewUp"),&up,SHADER_UNIFORM_VEC3);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"viewport"),&viewport,SHADER_UNIFORM_VEC4);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"hazeColor"),&haze,SHADER_UNIFORM_VEC3);
        Vector4 trail;float glow;shooting_star_sample(time,&trail,&glow);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"meteorTrail"),&trail,SHADER_UNIFORM_VEC4);
        SetShaderValue(nightSkyShader,GetShaderLocation(nightSkyShader,"meteorGlow"),&glow,SHADER_UNIFORM_FLOAT);
        Texture2D space=generated_space_texture();
        rlDisableDepthMask();BeginShaderMode(nightSkyShader);
        SetShaderValueTexture(nightSkyShader,GetShaderLocation(nightSkyShader,"spaceTexture"),space);
        DrawRectangle(0,0,w,h,WHITE);EndShaderMode();
        rlDrawRenderBatchActive();rlEnableDepthMask();
    }

    // Sky geometry follows translation but retains world orientation, like a skybox.
    BeginMode3D(camera);rlDisableDepthMask();
    for(int n=0;n<12;n++) {
        float angle=n*(2*PI/12)+.13f,r=130,height=30+(n*17%43);
        Vector3 p=Vector3Add(camera.position,(Vector3){sinf(angle)*r,height*.5f-8-camera.position.y,cosf(angle)*r});
        Color c={28+(n%3)*3,42+(n%3)*3,53+(n%3)*3,255};
        DrawCube(p,10+(n%3)*5,height,12,c);
        DrawCube(Vector3Add(p,(Vector3){2,height*.45f,0}),7,height*.25f,8,c);
        DrawCube(Vector3Add(p,(Vector3){-9,height*.15f,0}),9,height*.6f,9,c);
        for(int k=0;k<3;k++)DrawCube(Vector3Add(p,(Vector3){0,k*4-height*.1f,-6.1f}),.6f,1.6f,.1f,(Color){47,76,91,255});
        if(n%3==0) {
            DrawCube(Vector3Add(p,(Vector3){9,height*.25f,0}),30,3,5,c);
            DrawCube(Vector3Add(p,(Vector3){-3,height*.5f+5,0}),3,10,3,c);
        }
    }
    rlDrawRenderBatchActive();rlEnableDepthMask();EndMode3D();
}

static void victory_town_validate(void){
 for(int m=0;m<T_COUNT;m++)for(int i=0;i<town.models[m].meshes[0].vertexCount*3;i++){assert(isfinite(town.models[m].meshes[0].vertices[i]));assert(isfinite(town.models[m].meshes[0].normals[i]));}
 for(int scene=0;scene<4;scene++)for(int k=0;k<9;k++){TownSample s=victory_town_sample(scene,k);Vector3 delta=Vector3Subtract(s.camera.target,s.camera.position);RayCollision hit=GetRayCollisionMesh((Ray){s.camera.position,Vector3Normalize(delta)},town.models[T_STONE].meshes[0],MatrixIdentity());if(hit.hit&&hit.distance<.5f){fprintf(stderr,"Camera intersects scenery: scene %d time %d distance %.3f\n",scene,k,hit.distance);assert(false);}}
 puts("Town finite geometry and camera clearance passed");
}
