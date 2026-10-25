#include "golem_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "golem_body.h"
static Model gm[GM_COUNT];
static BoundingBox grBounds[GM_COUNT];
static Shader gs;
static bool grReady,grShaded;
static GolemStats grStats;
static int grMode,grRebuild,grGlow,grTime,grHit,grOpacity;
static bool grUniformValid;
static struct { int mode,glow;float rebuild,time,hit,alpha; } grUniform;
static const int grParent[GOLEM_JOINT_COUNT]={-1,0,1,2,2,4,5,2,7,8,1,10,11,1,13,14};
static const Color grSteel={119,137,150,255},grDark={47,57,66,255},grOrange={204,121,44,255};
static float gr_ease(float a,float b,float x){float t=Clamp((x-a)/(b-a),0,1);return t*t*(3-2*t);}
static Vector3 gr_add(Vector3 a,Vector3 b){return Vector3Add(a,b);}
static Matrix gr_at(Vector3 p){return MatrixTranslate(p.x,p.y,p.z);}
static Matrix gr_local(Vector3 p,Vector3 r,Matrix parent){return MatrixMultiply(MatrixMultiply(MatrixRotateXYZ(r),gr_at(p)),parent);}
static Vector3 gr_origin(Matrix m){return Vector3Transform((Vector3){0},m);}
static Matrix gr_aim(Vector3 a,Vector3 b,Matrix root){
    Vector3 d=Vector3Normalize(Vector3Subtract(b,a));
    Matrix r=QuaternionToMatrix(QuaternionFromVector3ToVector3((Vector3){0,-1,0},d));
    return MatrixMultiply(MatrixMultiply(r,gr_at(a)),root);
}
/* Two-bone analytic IK, with an explicit bend direction and clamped reach. */
static Vector3 gr_elbow(Vector3 a,Vector3 target,float l1,float l2,Vector3 bend){
    Vector3 delta=Vector3Subtract(target,a);float length=Vector3Length(delta);
    Vector3 axis=length>.00001f?Vector3Scale(delta,1/length):(Vector3){0,-1,0};
    float d=Clamp(length,fabsf(l1-l2)+.001f,l1+l2-.001f);
    Vector3 side=Vector3Subtract(bend,Vector3Scale(axis,Vector3DotProduct(bend,axis)));
    if(Vector3Length(side)<.001f)side=Vector3CrossProduct(axis,(Vector3){1,0,0});
    side=Vector3Normalize(side);
    float along=(l1*l1-l2*l2+d*d)/(2*d),height=sqrtf(fmaxf(0,l1*l1-along*along));
    return gr_add(a,gr_add(Vector3Scale(axis,along),Vector3Scale(side,height)));
}
GolemInput golem_input(void){GolemInput p={0};p.scale=1;p.hitJoint=-1;p.hitAge=100;p.seed=17;return p;}
float golem_action_seconds(GolemAction a){return a==GOLEM_SLAM?2.2f:a==GOLEM_SWEEP?1.8f:a==GOLEM_THROW?2.4f:4;}
static float gr_visible(GolemDamage d){return d.removed?Clamp(d.regeneration,0,1):1;}
GolemPose golem_sample(const GolemInput *source){
    GolemPose p={0};p.input=*source;GolemInput *in=&p.input;
    if(!isfinite(in->scale)||in->scale<=0)in->scale=1;
    float time=in->time,progress=Clamp(in->actionProgress,0,1);
    float left=in->damage[GOLEM_LEFT_LEG].removed?1-Clamp(in->damage[GOLEM_LEFT_LEG].recovery,0,1):0;
    float right=in->damage[GOLEM_RIGHT_LEG].removed?1-Clamp(in->damage[GOLEM_RIGHT_LEG].recovery,0,1):0;
    left*=gr_ease(0,.6f,in->damage[GOLEM_LEFT_LEG].age);right*=gr_ease(0,.6f,in->damage[GOLEM_RIGHT_LEG].age);
    p.kneel=fmaxf(left,right);
    float stagger=0;
    for(int s=0;s<2;s++)if(in->damage[s].removed)stagger+=(s?1:-1)*sinf(PI*Clamp(in->damage[s].age/.6f,0,1))*.14f;
    bool damagedAction=in->damage[in->side?GOLEM_RIGHT_ARM:GOLEM_LEFT_ARM].removed && gr_visible(in->damage[in->side?1:0])<1;
    GolemAction action=p.kneel>.01f||damagedAction||fabsf(stagger)>.01f?GOLEM_IDLE:in->action;
    float slam=action==GOLEM_SLAM?gr_ease(.40f,.60f,progress)*(1-gr_ease(.70f,1,progress)):0;
    Matrix root=MatrixMultiply(MatrixMultiply(MatrixScale(in->scale,in->scale,in->scale),MatrixRotateY(in->yaw*DEG2RAD)),gr_at(in->position));
    p.joint[GOLEM_ROOT]=root;
    float gait=in->travel/1.2f;
    bool walking=in->action==GOLEM_WALK && p.kneel<.001f;
    float bob=walking?.025f*cosf(gait*4*PI):.006f*sinf(time*1.6f);
    float pelvisY=Lerp(1.82f,.78f,p.kneel)+bob*(1-p.kneel)-slam*.95f;
    float lowSweep=in->combatReach&&action==GOLEM_SWEEP?gr_ease(0,.35f,progress)*(1-gr_ease(.7f,1,progress)):0;
    pelvisY-=.70f*lowSweep;
    Matrix pelvis=gr_local((Vector3){0,pelvisY,0},(Vector3){0,0,stagger},MatrixIdentity());
    p.joint[GOLEM_PELVIS]=MatrixMultiply(pelvis,root);
    float recoil=sinf(PI*Clamp(in->hitAge/.6f,0,1))*.10f;
    float torsion=0,lean=-.72f*p.kneel;
    if(in->hitJoint==GOLEM_TORSO)lean+=recoil;
    if(action==GOLEM_SWEEP)torsion=(in->side?1:-1)*(-.55f*gr_ease(0,.35f,progress)+1.2f*gr_ease(.35f,.62f,progress))*(1-gr_ease(.72f,1,progress));
    if(action==GOLEM_SLAM)lean-=.90f*slam;
    lean-=.65f*lowSweep;
    Matrix torso=gr_local((Vector3){0,.28f,0},(Vector3){lean,torsion,0},pelvis);
    p.joint[GOLEM_TORSO]=MatrixMultiply(torso,root);
    p.joint[GOLEM_HEAD]=gr_local((Vector3){0,1.28f,0},(Vector3){.15f*p.kneel,.08f*sinf(time*.6f),0},p.joint[GOLEM_TORSO]);
    p.chest=p.core=p.joint[GOLEM_TORSO];p.eye=p.joint[GOLEM_HEAD];
    p.chestVisibility=gr_visible(in->damage[GOLEM_CHEST]);
    for(int i=0;i<GOLEM_JOINT_COUNT;i++)p.visibility[i]=1;
    for(int s=0;s<2;s++){
        int thigh=s?GOLEM_THIGH_R:GOLEM_THIGH_L,shin=thigh+1,foot=thigh+2;
        float side=s?1:-1,loss=s?right:left;
        Vector3 hip=Vector3Transform((Vector3){side*.42f,-.03f,0},pelvis);
        float phase=gait+s*.5f;phase-=floorf(phase);
        float step=phase<.6f?-.36f+1.2f*phase:Lerp(.36f,-.36f,gr_ease(.6f,1,phase));
        Vector3 ankle={side*.42f,.22f,walking?step:0};float footYaw=in->yaw;
        if(walking&&phase>=.6f)ankle.y+=.18f*sinf((phase-.6f)/.4f*PI);
        if(walking&&in->feetResolved){ankle=Vector3Transform(in->feet[s],MatrixInvert(root));footYaw=in->footYaw[s];}
        ankle.z=Lerp(ankle.z,-.65f,p.kneel);
        /* Fold the absent shin behind the knee. Its regeneration target must
         * stay above the floor, and both bones retain their actual lengths. */
        ankle.z=Lerp(ankle.z,.18f,loss);
        Vector3 knee=gr_elbow(hip,ankle,.82f,.85f,(Vector3){0,0,-1});
        p.joint[thigh]=gr_aim(hip,knee,root);p.joint[shin]=gr_aim(knee,ankle,root);
        p.joint[foot]=gr_local(ankle,(Vector3){0,(footYaw-in->yaw)*DEG2RAD,0},root);
        p.knee[s]=Vector3Transform(knee,root);
        p.visibility[shin]=p.visibility[foot]=gr_visible(in->damage[GOLEM_LEFT_LEG+s]);
        int arm=s?GOLEM_ARM_R:GOLEM_ARM_L,fore=arm+1,hand=arm+2;
        float swing=walking?.23f*sinf(gait*2*PI+s*PI):0,spread=-side*.06f;
        float elbow=.13f,turn=0;
        if(in->hitJoint>=arm && in->hitJoint<=hand)swing-=recoil;
        if(s==in->side){
            if(action==GOLEM_SLAM){float raise=gr_ease(0,.4f,progress),strike=gr_ease(.45f,.60f,progress),recover=gr_ease(.7f,1,progress);swing=Lerp(2.7f*raise-2.7f*strike,0,recover);elbow=.5f*(1-strike);}
            if(action==GOLEM_SWEEP){swing=1.25f*(gr_ease(0,.3f,progress)-gr_ease(.7f,1,progress));turn=side*(-.85f+1.7f*gr_ease(.35f,.65f,progress))*(1-gr_ease(.72f,1,progress));elbow=.25f;}
            if(action==GOLEM_THROW){swing=(2.4f*gr_ease(0,.4f,progress)-1.0f*gr_ease(.45f,.62f,progress))*(1-gr_ease(.72f,1,progress));elbow=1.0f*(1-gr_ease(.4f,.6f,progress));}
        }
        Matrix upper=gr_local((Vector3){side*1.20f,1,0},(Vector3){swing,turn,spread},torso);
        Matrix lower=gr_local((Vector3){0,-.72f,0},(Vector3){elbow,0,0},upper);
        Matrix palm=gr_local((Vector3){0,-.68f,0},(Vector3){0},lower);
        bool brace=p.kneel>.01f && !in->damage[s].removed && ((left>.01f&&s==0) || (right>.01f&&s==1) || in->damage[1-s].removed);
        bool combatSweep=in->combatReach&&action==GOLEM_SWEEP&&s==in->side;
        if(brace || combatSweep || (action==GOLEM_SLAM && s==in->side && slam>0)){
            float weight=brace?p.kneel:combatSweep?lowSweep:slam;
            Vector3 shoulder=gr_origin(upper),wrist={side*1.20f,.39f,shoulder.z-.10f};
            if(in->combatReach&&!brace){
                wrist=Vector3Transform(in->attackTarget,MatrixInvert(root));
                wrist.y=combatSweep?1.0f:.42f;
                if(combatSweep)wrist.x+=Lerp(-1.2f,1.2f,gr_ease(.35f,.7f,progress));
            }
            wrist=Vector3Lerp(gr_origin(palm),wrist,weight);
            Vector3 reach=Vector3Subtract(wrist,shoulder);
            if(Vector3Length(reach)>1.399f)wrist=gr_add(shoulder,Vector3Scale(Vector3Normalize(reach),1.399f));
            Vector3 joint=gr_elbow(shoulder,wrist,.72f,.68f,(Vector3){side,0,0});
            upper=gr_aim(shoulder,joint,MatrixIdentity());lower=gr_aim(joint,wrist,MatrixIdentity());
            palm=gr_at(wrist);
        }
        p.joint[arm]=MatrixMultiply(upper,root);p.joint[fore]=MatrixMultiply(lower,root);p.joint[hand]=MatrixMultiply(palm,root);
        p.hand[s]=Vector3Transform((Vector3){0,-.39f,0},p.joint[hand]);
        if(!in->damage[s].forearmOnly)p.visibility[arm]=gr_visible(in->damage[s]);
        p.visibility[fore]=p.visibility[hand]=gr_visible(in->damage[s]);
    }
    if(action==GOLEM_SLAM){p.impact=(gr_ease(.60f,.62f,progress)-gr_ease(.62f,.85f,progress));p.contact=p.hand[in->side?1:0];p.contact.y=in->position.y+.015f;}
    if(action==GOLEM_SWEEP)p.sweep=gr_ease(.35f,.4f,progress)-gr_ease(.6f,.7f,progress);
    return p;
}

/* Rigid-group vertices carry square face UVs, cube centers and assembly order. */
typedef struct GrBuilder {float v[36*128*3],n[36*128*3],uv[36*128*2],cell[36*128*4];unsigned char c[36*128*4];int vertices,cubes;} GrBuilder;
static GrBuilder gb;
static void gr_cube(Vector3 center,Vector3 size,Color col){
    assert(gb.cubes<128);
    const int face[6][4]={{0,2,3,1},{4,5,7,6},{0,1,5,4},{2,6,7,3},{0,4,6,2},{1,3,7,5}};
    const int tri[6]={0,1,2,0,2,3};const Vector2 uv[4]={{0,0},{0,1},{1,1},{1,0}};
    Vector3 corner[8];for(int k=0;k<8;k++)corner[k]=gr_add(center,(Vector3){(k&1?.5f:-.5f)*size.x,(k&2?.5f:-.5f)*size.y,(k&4?.5f:-.5f)*size.z});
    for(int f=0;f<6;f++){
        Vector3 a=corner[face[f][0]],b=corner[face[f][1]],c=corner[face[f][2]];
        Vector3 normal=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b,a),Vector3Subtract(c,a)));
        float shade=.93f+.012f*((gb.cubes*17+f*3)%6);
        for(int k=0;k<6;k++){
            int i=gb.vertices++,j=tri[k];memcpy(gb.v+i*3,&corner[face[f][j]],sizeof(Vector3));memcpy(gb.n+i*3,&normal,sizeof(Vector3));
            gb.uv[i*2]=uv[j].x;gb.uv[i*2+1]=uv[j].y;
            gb.cell[i*4]=center.x;gb.cell[i*4+1]=center.y;gb.cell[i*4+2]=center.z;gb.cell[i*4+3]=gb.cubes;
            gb.c[i*4]=col.r*shade;gb.c[i*4+1]=col.g*shade;gb.c[i*4+2]=col.b*shade;gb.c[i*4+3]=255;
        }
    }gb.cubes++;
}
static void gr_finish(int kind){
    if(!gb.vertices)return;
    Mesh m={0};m.vertexCount=gb.vertices;m.triangleCount=m.vertexCount/3;
    m.vertices=MemAlloc(m.vertexCount*3*sizeof(float));m.normals=MemAlloc(m.vertexCount*3*sizeof(float));
    m.texcoords=MemAlloc(m.vertexCount*2*sizeof(float));m.tangents=MemAlloc(m.vertexCount*4*sizeof(float));m.colors=MemAlloc(m.vertexCount*4);
    for(int i=0;i<m.vertexCount;i++)gb.cell[i*4+3]/=fmaxf(1,gb.cubes-1);
    memcpy(m.vertices,gb.v,m.vertexCount*3*sizeof(float));memcpy(m.normals,gb.n,m.vertexCount*3*sizeof(float));
    memcpy(m.texcoords,gb.uv,m.vertexCount*2*sizeof(float));memcpy(m.tangents,gb.cell,m.vertexCount*4*sizeof(float));memcpy(m.colors,gb.c,m.vertexCount*4);
    UploadMesh(&m,false);gm[kind]=LoadModelFromMesh(m);grBounds[kind]=GetMeshBoundingBox(m);
    grStats.cubes+=gb.cubes;grStats.triangles+=m.triangleCount;grStats.models++;memset(&gb,0,sizeof(gb));
}
bool golem_renderer_init(const char *dir){
    if(grReady)return true;
    grStats=(GolemStats){0};grUniformValid=false;memset(&gb,0,sizeof(gb));
    golem_body_generate(gr_cube,gr_finish);
    char v[512],f[512];if(!dir)dir="shaders";snprintf(v,sizeof(v),"%s/golem.vert",dir);snprintf(f,sizeof(f),"%s/golem.frag",dir);
    gs=LoadShader(v,f);grShaded=gs.id && gs.id!=rlGetShaderIdDefault();
    if(grShaded){
        gs.locs[SHADER_LOC_VERTEX_TANGENT]=GetShaderLocationAttrib(gs,"vertexTangent");
        grMode=GetShaderLocation(gs,"renderMode");grRebuild=GetShaderLocation(gs,"rebuild");grGlow=GetShaderLocation(gs,"luminous");
        grTime=GetShaderLocation(gs,"clockTime");grHit=GetShaderLocation(gs,"hitFlash");grOpacity=GetShaderLocation(gs,"effectOpacity");
        for(int k=0;k<GM_COUNT;k++)if(gm[k].meshCount)gm[k].materials[0].shader=gs;
    }
    grStats.shader=grShaded;grReady=true;return true;
}
static void gr_draw_part(int model,Matrix transform,float rebuild,int mode,int glow,float time,float flash,float alpha){
    if(!gm[model].meshCount)return;
    if(grShaded){
        /* The shader belongs exclusively to this renderer. Most groups share
         * every uniform, and split-screen reuses the same sampled time. */
        if(!grUniformValid||grUniform.mode!=mode)SetShaderValue(gs,grMode,&mode,SHADER_UNIFORM_INT);
        if(!grUniformValid||grUniform.rebuild!=rebuild)SetShaderValue(gs,grRebuild,&rebuild,SHADER_UNIFORM_FLOAT);
        if(!grUniformValid||grUniform.glow!=glow)SetShaderValue(gs,grGlow,&glow,SHADER_UNIFORM_INT);
        if(!grUniformValid||grUniform.time!=time)SetShaderValue(gs,grTime,&time,SHADER_UNIFORM_FLOAT);
        if(!grUniformValid||grUniform.hit!=flash)SetShaderValue(gs,grHit,&flash,SHADER_UNIFORM_FLOAT);
        if(!grUniformValid||grUniform.alpha!=alpha)SetShaderValue(gs,grOpacity,&alpha,SHADER_UNIFORM_FLOAT);
        grUniform.mode=mode;grUniform.rebuild=rebuild;grUniform.glow=glow;grUniform.time=time;grUniform.hit=flash;grUniform.alpha=alpha;grUniformValid=true;
    }else{
        if(mode>=2||rebuild<.5f)return;
        gm[model].materials[0].maps[MATERIAL_MAP_DIFFUSE].color=mode==1&&!glow?BLACK:WHITE;
    }
    DrawMesh(gm[model].meshes[0],gm[model].materials[0],transform);
}
static int gr_section_parts(int section,bool forearm,int indices[3]){
    if(section<2){int first=section?GOLEM_ARM_R:GOLEM_ARM_L;indices[0]=first+(forearm?1:0);indices[1]=first+2;if(forearm)return 2;indices[1]=first+1;indices[2]=first+2;return 3;}
    if(section<4){indices[0]=section==GOLEM_LEFT_LEG?GOLEM_SHIN_L:GOLEM_SHIN_R;indices[1]=indices[0]+1;return 2;}
    indices[0]=GM_CHEST;return 1;
}
void golem_capture_damage(const GolemPose *p,GolemSection section,bool forearm,GolemDamage *event){
    *event=(GolemDamage){.removed=true,.captured=true,.forearmOnly=forearm};
    int parts[3],count=gr_section_parts(section,forearm,parts);
    for(int i=0;i<count;i++)event->detached[i]=parts[i]==GM_CHEST?p->chest:p->joint[parts[i]];
    event->origin=gr_origin(event->detached[0]);
    event->yaw=p->input.yaw;event->floorY=p->input.position.y;
}
static void gr_detached(const GolemPose *p,int section,WeaponRenderPass pass){
    /* Opaque world depth is already copied into the emission target. Detached
     * unlit armor needs no duplicate black draw in that pass. */
    GolemDamage d=p->input.damage[section];if(!d.removed||!d.captured||d.age<0||d.age>4||pass!=WEAPON_COLOR)return;
    int parts[3],count=gr_section_parts(section,d.forearmOnly,parts);Matrix transforms[3];
    float t=fminf(d.age,1.3f),side=section==GOLEM_LEFT_ARM||section==GOLEM_LEFT_LEG?-1:1;
    Vector3 travel={side*(1-expf(-t*2))*.9f,2*t-4.9f*t*t,-.9f*t};
    if(d.forearmOnly){travel.x*=.2f;travel.z=-5*t;travel.y=2.7f*t-4.9f*t*t;}
    travel=Vector3RotateByAxisAngle(travel,(Vector3){0,1,0},d.yaw*DEG2RAD);
    Matrix move=MatrixMultiply(gr_at(Vector3Negate(d.origin)),MatrixRotateXYZ((Vector3){t*.7f,t*.5f,side*t*.8f}));
    move=MatrixMultiply(move,gr_at(gr_add(d.origin,travel)));
    float bottom=10000;
    for(int i=0;i<count;i++){
        transforms[i]=MatrixMultiply(d.detached[i],move);BoundingBox b=grBounds[parts[i]];
        for(int k=0;k<8;k++){Vector3 v={k&1?b.max.x:b.min.x,k&2?b.max.y:b.min.y,k&4?b.max.z:b.min.z};bottom=fminf(bottom,Vector3Transform(v,transforms[i]).y);}
    }
    float lift=fmaxf(0,d.floorY+.015f-bottom),alpha=1-gr_ease(3,4,d.age);
    for(int i=0;i<count;i++)gr_draw_part(parts[i],MatrixMultiply(transforms[i],MatrixTranslate(0,lift,0)),1,pass==WEAPON_EMISSION?1:0,0,p->input.time,0,alpha);
}
void golem_draw(const GolemPose *p,Vector3 eye,WeaponRenderPass pass,bool effects){
    (void)eye;if(!grReady)return;
    int mode=pass==WEAPON_EMISSION?1:pass==WEAPON_GLASS?2:0;
    if(pass==WEAPON_GLASS){rlDrawRenderBatchActive();rlDisableDepthMask();}
    for(int j=1;j<GOLEM_JOINT_COUNT;j++){
        float progress=p->visibility[j];if(pass==WEAPON_GLASS && (progress<=0||progress>=1||!effects))continue;
        if(progress<=0 || (pass==WEAPON_EMISSION && progress>=1))continue;
        float flash=j==p->input.hitJoint?1-gr_ease(0,.22f,p->input.hitAge):0;
        gr_draw_part(j,p->joint[j],progress,mode,0,p->input.time,flash,1);
    }
    if(p->chestVisibility>0 && !(pass==WEAPON_EMISSION && p->chestVisibility>=1) &&
       (pass!=WEAPON_GLASS || (effects && p->chestVisibility<1)))
        gr_draw_part(GM_CHEST,p->chest,p->chestVisibility,mode,0,p->input.time,p->input.hitJoint==GOLEM_TORSO?1-gr_ease(0,.22f,p->input.hitAge):0,1);
    if(pass!=WEAPON_GLASS){
        gr_draw_part(GM_CORE,p->core,1,mode,1,p->input.time,p->input.hitJoint==GM_CORE?1-gr_ease(0,.22f,p->input.hitAge):0,1);
        gr_draw_part(GM_EYE,p->eye,1,mode,1,p->input.time,0,1);
        if(effects)for(int s=0;s<GOLEM_SECTION_COUNT;s++)gr_detached(p,s,pass);
    }
    if(pass==WEAPON_GLASS){rlDrawRenderBatchActive();rlEnableDepthMask();}
    if(effects && pass==WEAPON_EMISSION && grShaded){
        rlDrawRenderBatchActive();rlDisableDepthMask();
        for(int j=1;j<GOLEM_JOINT_COUNT;j++)if(p->visibility[j]>0 && p->visibility[j]<1)gr_draw_part(j,p->joint[j],p->visibility[j],3,0,p->input.time,0,1);
        if(p->chestVisibility>0 && p->chestVisibility<1)gr_draw_part(GM_CHEST,p->chest,p->chestVisibility,3,0,p->input.time,0,1);
        rlDrawRenderBatchActive();rlEnableDepthMask();
    }
    if(effects && pass==WEAPON_COLOR && p->impact>0){
        float age=1-p->impact;
        DrawCircle3D(p->contact,.25f+age*1.6f,(Vector3){1,0,0},90,Fade(ORANGE,p->impact));
        for(int n=0;n<12;n++){float a=n*2.39996f;Vector3 at=gr_add(p->contact,(Vector3){cosf(a)*age,sinf(age*PI)*.22f,sinf(a)*age});DrawCube(at,.06f,.06f,.06f,Fade(grSteel,p->impact));}
    }
}
void golem_draw_skeleton(const GolemPose *p){
    for(int j=1;j<GOLEM_JOINT_COUNT;j++){
        Vector3 a=gr_origin(p->joint[grParent[j]]),b=gr_origin(p->joint[j]);
        DrawLine3D(a,b,YELLOW);DrawSphere(b,.045f*p->input.scale,j>=GOLEM_THIGH_L?SKYBLUE:ORANGE);
    }
}
GolemStats golem_renderer_stats(void){return grStats;}
void golem_renderer_shutdown(void){
    if(!grReady)return;
    for(int i=0;i<GM_COUNT;i++)if(gm[i].meshCount)UnloadModel(gm[i]);
    if(grShaded)UnloadShader(gs);
    memset(gm,0,sizeof(gm));memset(grBounds,0,sizeof(grBounds));grStats=(GolemStats){0};grReady=grShaded=false;
}
