#include "character_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <string.h>

/* Model coordinates: feet at -0.25, forward -Z. All triangles have flat normals. */
enum { CR_BODY, CR_HEAD, CR_EAR, CR_INNER, CR_PAW, CR_LEG, CR_TAIL,
       CR_LENS, CR_HOUSING, CR_EYE, CR_NOSE, CR_BADGE, CR_COUNT };
static Model parts[CR_COUNT];
static Shader shader;
static bool ready, shaded;
static int maskLoc, glowLoc, tintLoc;
typedef struct Builder { float v[6144], n[6144]; unsigned char c[8192]; int count; } Builder;
static Builder builder;
static const Color fur={83,84,86,255}, white={224,226,224,255};

static void triangle(Vector3 a,Vector3 b,Vector3 c,Color color) {
    Vector3 normal=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b,a),Vector3Subtract(c,a)));
    Vector3 points[]={a,b,c};
    float shade=.85f+.025f*((builder.count/3*17)%7);
    for(int k=0;k<3;k++) {
        int i=builder.count++;
        memcpy(builder.v+3*i,&points[k],sizeof(Vector3));
        memcpy(builder.n+3*i,&normal,sizeof(Vector3));
        builder.c[4*i]=(unsigned char)(color.r*shade);
        builder.c[4*i+1]=(unsigned char)(color.g*shade);
        builder.c[4*i+2]=(unsigned char)(color.b*shade);
        builder.c[4*i+3]=255;
    }
}
/* Faceted ellipsoid with staggered rings, including a white front on the body. */
static void ellipsoid(int part,Vector3 size,Color color,bool belly) {
    memset(&builder,0,sizeof(builder));
    const int sides=10,rings=6;
    for(int r=0;r<rings;r++) for(int j=0;j<sides;j++) {
        Vector3 p[4];
        for(int k=0;k<4;k++) {
            int row=r+(k>=2),col=j+(k==1||k==2);
            float lat=-PI/2+PI*row/rings,angle=2*PI*col/sides;
            p[k]=(Vector3){size.x*cosf(lat)*cosf(angle),size.y*sinf(lat),size.z*cosf(lat)*sinf(angle)};
        }
        Color face=belly && (j==6||j==7||j==8) ? white:color;
        // Winding faces outward; poles emit only one nondegenerate triangle.
        if(r>0) triangle(p[0],p[2],p[1],face);
        if(r<rings-1) triangle(p[0],p[3],p[2],face);
    }
    Mesh mesh={0};mesh.vertexCount=builder.count;mesh.triangleCount=builder.count/3;
    mesh.vertices=MemAlloc(builder.count*3*sizeof(float));mesh.normals=MemAlloc(builder.count*3*sizeof(float));
    mesh.colors=MemAlloc(builder.count*4);
    memcpy(mesh.vertices,builder.v,builder.count*3*sizeof(float));
    memcpy(mesh.normals,builder.n,builder.count*3*sizeof(float));memcpy(mesh.colors,builder.c,builder.count*4);
    UploadMesh(&mesh,false);parts[part]=LoadModelFromMesh(mesh);
}
bool character_renderer_init(void) {
    if(ready) return true;
    ellipsoid(CR_BODY,(Vector3){.205f,.235f,.16f},fur,true);
    ellipsoid(CR_HEAD,(Vector3){.135f,.105f,.19f},fur,false);
    // Pull the forward tip into a long tapered snout.
    Mesh *head=&parts[CR_HEAD].meshes[0];
    for(int i=0;i<head->vertexCount;i++) if(head->vertices[3*i+2]<-.08f) {
        head->vertices[3*i]*=.65f;head->vertices[3*i+1]*=.65f;
        head->vertices[3*i+2]-=.035f;
    }
    for(int i=0;i<head->vertexCount;i+=3) {
        Vector3 *a=(Vector3*)(head->vertices+i*3);
        Vector3 n=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(a[1],a[0]),Vector3Subtract(a[2],a[0])));
        for(int j=0;j<3;j++) memcpy(head->normals+(i+j)*3,&n,sizeof(n));
    }
    UpdateMeshBuffer(*head,0,head->vertices,head->vertexCount*3*sizeof(float),0);
    UpdateMeshBuffer(*head,2,head->normals,head->vertexCount*3*sizeof(float),0);
    ellipsoid(CR_EAR,(Vector3){.080f,.102f,.025f},fur,false);
    ellipsoid(CR_INNER,(Vector3){.061f,.080f,.012f},(Color){164,131,138,255},false);
    ellipsoid(CR_PAW,(Vector3){.040f,.045f,.060f},white,false);
    ellipsoid(CR_LEG,(Vector3){.060f,.085f,.065f},fur,false);
    ellipsoid(CR_TAIL,(Vector3){.023f,.023f,.12f},fur,false);
    parts[CR_LENS]=LoadModelFromMesh(GenMeshPlane(.066f,.066f,1,1));
    parts[CR_HOUSING]=LoadModelFromMesh(GenMeshCube(.1f,.1f,.1f));
    ellipsoid(CR_EYE,(Vector3){.012f,.015f,.009f},(Color){6,7,8,255},false);
    ellipsoid(CR_NOSE,(Vector3){.032f,.030f,.033f},(Color){28,29,31,255},false);
    ellipsoid(CR_BADGE,(Vector3){.055f,.043f,.012f},WHITE,false);
    shader=LoadShader("shaders/character.vert","shaders/character.frag");
    shaded=shader.id!=0 && shader.id!=rlGetShaderIdDefault();
    if(shaded) {
        maskLoc=GetShaderLocation(shader,"emissionOnly");glowLoc=GetShaderLocation(shader,"luminous");
        tintLoc=GetShaderLocation(shader,"statusTint");
        for(int i=0;i<CR_COUNT;i++) parts[i].materials[0].shader=shader;
    }
    ready=true;return true;
}
void character_renderer_shutdown(void) {
    if(!ready) return;
    for(int i=0;i<CR_COUNT;i++) UnloadModel(parts[i]);
    if(shaded) UnloadShader(shader);
    memset(parts,0,sizeof(parts));ready=shaded=false;
}
CharacterPose character_pose(Vector3 pos,float yaw,float pitch,float speed,bool grounded,bool downed,Color gear,float time) {
    return (CharacterPose){.position=pos,.yaw=yaw,.pitch=Clamp(pitch,-25,25),
        .stride=sinf(time*10)*Clamp(speed/3,0,1),.breath=sinf(time*2)*.003f,
        .time=time,.airborne=!grounded,.downed=downed,.gear=gear};
}
/* Body and view share an eye anchor; only the torso gets breathing/bobbing. */
Matrix character_body_transform(CharacterPose p) {
    Matrix root=MatrixScale(CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE);
    if(p.downed) root=MatrixMultiply(root,MatrixRotateZ(PI/2));
    root=MatrixMultiply(root,MatrixRotateY(p.yaw*DEG2RAD));
    float lift=p.downed?.205f*CHARACTER_MODEL_SCALE:.25f*CHARACTER_MODEL_SCALE+p.breath;
    return MatrixMultiply(root,MatrixTranslate(p.position.x,p.position.y-1+lift,p.position.z));
}
static Vector3 nominal_sole(CharacterPose p,int foot) {
    float side=foot?1:-1;
    Vector3 local={side*.10f*CHARACTER_MODEL_SCALE,0,-.045f*CHARACTER_MODEL_SCALE};
    local=Vector3Transform(local,MatrixRotateY(p.yaw*DEG2RAD));
    return Vector3Add((Vector3){p.position.x,p.position.y-1,p.position.z},local);
}
void character_animate(CharacterAnimation *s,CharacterPose *p,float dt,CharacterGroundQuery query,void *context) {
    dt=Clamp(dt,0,.1f);
    Vector3 movement=Vector3Subtract(p->position,s->previous);
    bool reset=!s->initialized || Vector3Length(movement)>1.5f || fabsf(movement.y)>.6f;
    if(reset) {
        *s=(CharacterAnimation){.initialized=true,.previous=p->position,.yaw=p->yaw,.grounded=!p->airborne,.landing=.1f};
        for(int f=0;f<2;f++)s->soles[f]=nominal_sole(*p,f);
        movement=(Vector3){0};
    }
    if(p->downed){s->initialized=false;return;}
    float travel=sqrtf(movement.x*movement.x+movement.z*movement.z);
    // Subtract motion of a shared support before advancing the walk cycle.
    Vector3 platform={0};int supports=0;
    for(int f=0;f<2;f++) if(!s->swinging[f] && s->support[f].valid && s->support[f].identity) {
        CharacterSupport now=query(s->soles[f],context);
        if(now.valid && now.identity==s->support[f].identity) {
            Vector3 delta=Vector3Subtract(now.origin,s->support[f].origin);
            s->soles[f]=Vector3Add(s->soles[f],delta);platform=Vector3Add(platform,delta);supports++;
        }
    }
    if(supports){movement=Vector3Subtract(movement,Vector3Scale(platform,1.0f/supports));travel=sqrtf(movement.x*movement.x+movement.z*movement.z);}
    float turn=fabsf(remainderf(p->yaw-s->yaw,360));
    bool overextended=false;
    for(int f=0;f<2;f++)if(Vector3Distance(s->soles[f],nominal_sole(*p,f))>.26f)overextended=true;
    bool walking=travel>.00001f;
    if(!p->airborne) {
        if(walking)s->phase+=travel/.60f; // short strides suit the rat's legs
        else if(turn>20 || overextended || s->swinging[0] || s->swinging[1])s->phase+=dt*1.8f;
    }
    if(!p->airborne && !s->grounded) {
        s->landing=0;
        for(int f=0;f<2;f++){s->landingStart[f]=s->soles[f];s->swinging[f]=false;}
    }
    s->landing=fminf(.1f,s->landing+dt);
    for(int f=0;f<2;f++) {
        Vector3 desired=nominal_sole(*p,f);
        if(p->airborne) {
            desired.y+=.065f*CHARACTER_MODEL_SCALE;desired.z+=.02f;
            s->soles[f]=desired;s->support[f]=(CharacterSupport){0};s->swinging[f]=false;
        } else {
            // Land ahead of the hips, then let the body pass the planted paw.
            if(walking) {
                desired.x+=movement.x/travel*.15f;
                desired.z+=movement.z/travel*.15f;
            }
            CharacterSupport ground=query(desired,context);
            if(ground.valid && fabsf(ground.height-desired.y)<=.12f)desired.y=ground.height;
            float phase=s->phase+f*.5f;phase-=floorf(phase);
            bool swing=phase>=.6f && (walking || overextended || s->swinging[f] || turn>20);
            CharacterSupport anchor=query(s->soles[f],context);
            if(!anchor.valid || (s->support[f].identity && anchor.identity!=s->support[f].identity)) {
                s->soles[f]=desired;s->support[f]=ground;
            }
            if(swing) {
                if(!s->swinging[f])s->swingStart[f]=s->soles[f];
                float t=Clamp((phase-.6f)/.4f,0,1),ease=t*t*(3-2*t);
                s->soles[f]=Vector3Lerp(s->swingStart[f],desired,ease);
                s->soles[f].y+=sinf(PI*t)*.10f;s->support[f]=ground;
            } else if(s->swinging[f] || reset || !s->grounded) {
                s->soles[f]=desired;s->support[f]=ground;
            } else if(anchor.valid) {s->soles[f].y=anchor.height;s->support[f]=anchor;}
            s->swinging[f]=swing;
            if(s->landing<.1f)s->soles[f]=Vector3Lerp(s->landingStart[f],s->soles[f],s->landing/.1f);
        }
        p->soles[f]=s->soles[f];
    }
    p->feetResolved=true;
    p->breath+=p->airborne?0:sinf(s->phase*4*PI)*.008f*fminf(1,travel/fmaxf(dt,.001f));
    s->grounded=!p->airborne;s->previous=p->position;
    if(walking || turn>20)s->yaw=p->yaw;
}
static void part(int id,Vector3 pos,Vector3 scale,Vector3 rotation,Matrix root,Color color,bool glow,bool emission) {
    if(!ready) return;
    if(!shaded && emission) return;
    Matrix m=MatrixMultiply(MatrixScale(scale.x,scale.y,scale.z),MatrixRotateXYZ(rotation));
    m=MatrixMultiply(m,MatrixTranslate(pos.x,pos.y,pos.z));m=MatrixMultiply(m,root);
    parts[id].transform=m;
    if(shaded) {int g=glow,e=emission;SetShaderValue(shader,glowLoc,&g,SHADER_UNIFORM_INT);SetShaderValue(shader,maskLoc,&e,SHADER_UNIFORM_INT);}
    DrawModel(parts[id],(Vector3){0},1,color);
}
#define ONE (Vector3){1,1,1}
#define ZERO (Vector3){0}
static void status(CharacterPose p) {
    if(shaded) {float t[]={p.effect.r/255.0f,p.effect.g/255.0f,p.effect.b/255.0f,p.effectAmount};SetShaderValue(shader,tintLoc,t,SHADER_UNIFORM_VEC4);}
}
static Matrix scouter_frame(Matrix head) {
    Matrix lens=MatrixMultiply(MatrixRotateY(35*DEG2RAD),MatrixTranslate(-.074f,.042f,-.139f));
    return MatrixMultiply(lens,head);
}
void character_draw(CharacterPose p,bool emission) {
    status(p);
    Matrix root=character_body_transform(p);
    part(CR_BODY,(Vector3){0,.015f,.035f},ONE,ZERO,root,WHITE,false,emission);
    part(CR_BADGE,(Vector3){0,.11f,.184f},ONE,ZERO,root,p.gear,true,emission);
    Matrix head=MatrixMultiply(MatrixRotateX(p.pitch*DEG2RAD),MatrixTranslate(0,.245f,-.055f));
    head=MatrixMultiply(head,root);
    part(CR_HEAD,ZERO,ONE,ZERO,head,WHITE,false,emission);
    part(CR_NOSE,(Vector3){0,-.018f,-.225f},ONE,ZERO,head,WHITE,false,emission);
    for(int s=-1;s<=1;s+=2) {
        Vector3 tilt={0,0,-s*.2f};
        part(CR_EAR,(Vector3){s*.105f,.105f,.045f},ONE,tilt,head,WHITE,false,emission);
        part(CR_INNER,(Vector3){s*.105f,.105f,.024f},ONE,tilt,head,WHITE,false,emission);
        part(CR_EYE,(Vector3){s*.059f,.042f,-.118f},ONE,ZERO,head,WHITE,false,emission);
        if(p.feetResolved && !p.downed) {
            Vector3 sole=p.soles[s>0?1:0];
            Vector3 center=Vector3Add(sole,(Vector3){0,.045f*CHARACTER_MODEL_SCALE,0});
            Vector3 localCenter=Vector3Transform(center,MatrixInvert(root));
            Vector3 delta=Vector3Subtract(localCenter,(Vector3){s*.10f,-.105f,.015f});
            float length=Vector3Length(delta);
            Vector3 middle=Vector3Lerp(localCenter,(Vector3){s*.10f,-.105f,.015f},.5f);
            Vector3 rotation={atan2f(delta.z,-delta.y),0,atan2f(delta.x,sqrtf(delta.y*delta.y+delta.z*delta.z))};
            part(CR_LEG,middle,(Vector3){.8f,fmaxf(.3f,length/.17f),.8f},rotation,root,WHITE,false,emission);
            Matrix footRoot=MatrixMultiply(MatrixScale(CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE),MatrixRotateY(p.yaw*DEG2RAD));
            footRoot=MatrixMultiply(footRoot,MatrixTranslate(center.x,center.y,center.z));
            part(CR_PAW,ZERO,(Vector3){1,1,1.35f},ZERO,footRoot,WHITE,false,emission);
        } else {
            bool tuck=p.airborne && !p.downed;
            float step=p.downed?0:tuck?.055f:p.stride*s*.035f;
            part(CR_LEG,(Vector3){s*.10f,tuck?-.095f:-.125f,step},ONE,(Vector3){tuck?-.5f:p.stride*s*.3f,0,0},root,WHITE,false,emission);
            part(CR_PAW,(Vector3){s*.10f,tuck?-.14f:-.205f,-.045f+step},(Vector3){1,1,1.35f},ZERO,root,WHITE,false,emission);
        }
        if(p.downed) {
            part(CR_LEG,(Vector3){s*.18f,.09f,-.07f},(Vector3){.45f,.7f,.5f},(Vector3){.7f,0,0},root,WHITE,false,emission);
            part(CR_PAW,(Vector3){s*.18f,.035f,-.12f},ONE,ZERO,root,WHITE,false,emission);
        }
    }
    // One temple-mounted scouter, separate from the rat's natural eyes.
    Color metal={145,153,159,255};
    part(CR_HOUSING,(Vector3){-.122f,.043f,-.055f},(Vector3){.28f,.40f,.35f},ZERO,head,metal,false,emission);
    part(CR_HOUSING,(Vector3){-.105f,.045f,-.107f},(Vector3){.10f,.09f,.82f},ZERO,head,metal,false,emission);
    part(CR_HOUSING,(Vector3){-.138f,.043f,-.060f},(Vector3){.035f,.14f,.16f},ZERO,head,p.gear,true,emission);
    // Fine luminous square border; the panel itself belongs to the glass pass.
    Matrix lens=scouter_frame(head);
    for(int side=-1;side<=1;side+=2) {
        part(CR_HOUSING,(Vector3){side*.034f,0,0},(Vector3){.016f,.696f,.016f},ZERO,lens,p.gear,true,emission);
        part(CR_HOUSING,(Vector3){0,side*.034f,0},(Vector3){.664f,.016f,.016f},ZERO,lens,p.gear,true,emission);
    }
    Vector3 tail={0,-.10f,.18f};
    for(int i=0;i<3;i++) {
        float angle=.15f+i*.30f,sway=sinf(p.time*2+i*.4f)*.13f;
        Vector3 delta={sway*.12f,sinf(angle)*.16f,cosf(angle)*.16f};
        part(CR_TAIL,Vector3Add(tail,Vector3Scale(delta,.5f)),(Vector3){1-i*.24f,1-i*.24f,.8f},(Vector3){-angle,sway,0},root,WHITE,false,emission);
        tail=Vector3Add(tail,delta);
    }
}
void character_draw_visor(CharacterPose p) {
    if(!ready) return;
    Matrix root=character_body_transform(p);
    Matrix head=MatrixMultiply(MatrixRotateX(p.pitch*DEG2RAD),MatrixTranslate(0,.245f,-.055f));
    head=MatrixMultiply(head,root);
    status((CharacterPose){0});
    Color glass=p.gear;glass.a=48;
    // Render after opaque actors, preserve depth testing but never write glass depth.
    rlDrawRenderBatchActive();rlDisableDepthMask();rlDisableBackfaceCulling();
    part(CR_LENS,ZERO,ONE,(Vector3){PI/2,0,0},scouter_frame(head),glass,true,false);
    rlDrawRenderBatchActive();rlEnableBackfaceCulling();rlEnableDepthMask();
}
void character_draw_paws(WeaponPose weapon,Color gear,bool emission) {
    status((CharacterPose){0});
    // Positions are in the existing weapon's authored local coordinate system.
    part(CR_PAW,(Vector3){.17f,-.67f,.40f},(Vector3){3.2f,3.2f,2.7f},(Vector3){-.25f,0,0},weapon.transform,WHITE,false,emission);
    part(CR_PAW,(Vector3){-.36f,-.33f,weapon.kind==WEAPON_DYNAMIC_LAUNCHER?-.55f:-.30f},(Vector3){3,2.5f,3},ZERO,weapon.transform,WHITE,false,emission);
    part(CR_LEG,(Vector3){.25f,-.92f,.65f},(Vector3){2,3,3},(Vector3){-.6f,0,0},weapon.transform,WHITE,false,emission);
    part(CR_LEG,(Vector3){-.42f,-.55f,weapon.kind==WEAPON_DYNAMIC_LAUNCHER?-.32f:-.08f},(Vector3){1.6f,2.5f,2.5f},(Vector3){-.65f,0,0},weapon.transform,WHITE,false,emission);
    part(CR_BADGE,(Vector3){.28f,-.76f,.43f},(Vector3){1.6f,1.3f,1},(Vector3){0,PI/2,0},weapon.transform,gear,false,emission);
}
void character_draw_arm(Vector3 shoulder,Vector3 hand,Color gear) {
    status((CharacterPose){0});
    Vector3 d=Vector3Subtract(hand,shoulder);float length=Vector3Length(d);
    Matrix root=MatrixTranslate(hand.x,hand.y,hand.z);
    part(CR_PAW,ZERO,(Vector3){CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE,CHARACTER_MODEL_SCALE},ZERO,root,WHITE,false,false);
    if(length>.001f) DrawCylinderEx(shoulder,hand,.022f*CHARACTER_MODEL_SCALE,.028f*CHARACTER_MODEL_SCALE,6,fur);
    (void)gear;
}
