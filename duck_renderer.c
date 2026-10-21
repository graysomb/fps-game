#include "duck_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

enum { BODY, BREAST, NECK, COLLAR, HEAD, BILL, BILL_LINE, NOSTRIL,
       EYE, GLINT, WING, WING_BAR, SPECULUM, TAIL, TAIL_CURL, PART_COUNT };
static Model duckParts[PART_COUNT], floorModel;
static Shader duckShader, floorShader;
static bool ready, shaded, floorShaded;
static int materialLoc, clockLoc, speedLoc, wakeLoc;
static int triangles;

typedef struct DuckMeshBuilder {
    float position[4096*3], normal[4096*3];
    unsigned char color[4096*4];
    int count;
} DuckMeshBuilder;
static DuckMeshBuilder meshBuilder;

static void face(Vector3 a,Vector3 b,Vector3 c,Color color) {
    Vector3 normal=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b,a),Vector3Subtract(c,a)));
    Vector3 points[]={a,b,c};
    float variation=.94f+.012f*((meshBuilder.count/3*37)%6);
    assert(meshBuilder.count+3<=4096);
    for(int k=0;k<3;k++) {
        int i=meshBuilder.count++;
        memcpy(meshBuilder.position+i*3,&points[k],sizeof(Vector3));
        memcpy(meshBuilder.normal+i*3,&normal,sizeof(Vector3));
        meshBuilder.color[i*4]=(unsigned char)(color.r*variation);
        meshBuilder.color[i*4+1]=(unsigned char)(color.g*variation);
        meshBuilder.color[i*4+2]=(unsigned char)(color.b*variation);
        meshBuilder.color[i*4+3]=255;
    }
}

/* Flat triangles, broad planar facets and a flattened underside. */
static void ellipsoid(int part,Vector3 radii,Color color,int sides,int rings,float bottom) {
    meshBuilder.count=0;
    for(int row=0;row<rings;row++) for(int side=0;side<sides;side++) {
        Vector3 p[4];
        for(int k=0;k<4;k++) {
            int r=row+(k>=2),s=side+(k==1||k==2);
            float latitude=-PI/2+PI*r/rings,angle=2*PI*s/sides;
            p[k]=(Vector3){radii.x*cosf(latitude)*cosf(angle),
                fmaxf(bottom,radii.y*sinf(latitude)),radii.z*cosf(latitude)*sinf(angle)};
        }
        if(row>0) face(p[0],p[2],p[1],color);
        if(row<rings-1) face(p[0],p[3],p[2],color);
    }
    Mesh mesh={0};
    mesh.vertexCount=meshBuilder.count;mesh.triangleCount=mesh.vertexCount/3;
    mesh.vertices=MemAlloc(mesh.vertexCount*3*sizeof(float));
    mesh.normals=MemAlloc(mesh.vertexCount*3*sizeof(float));mesh.colors=MemAlloc(mesh.vertexCount*4);
    memcpy(mesh.vertices,meshBuilder.position,mesh.vertexCount*3*sizeof(float));
    memcpy(mesh.normals,meshBuilder.normal,mesh.vertexCount*3*sizeof(float));
    memcpy(mesh.colors,meshBuilder.color,mesh.vertexCount*4);
    UploadMesh(&mesh,false);duckParts[part]=LoadModelFromMesh(mesh);triangles+=mesh.triangleCount;
}

bool duck_renderer_init(const char *directory) {
    if(ready)return true;
    triangles=0;
    ellipsoid(BODY,(Vector3){.34f,.29f,.55f},(Color){171,169,150,255},12,6,-.18f);
    ellipsoid(BREAST,(Vector3){.29f,.235f,.27f},(Color){103,53,37,255},12,6,-.18f);
    ellipsoid(NECK,(Vector3){.135f,.21f,.14f},(Color){17,101,65,255},10,5,-.21f);
    ellipsoid(COLLAR,(Vector3){.148f,.028f,.153f},(Color){240,233,209,255},10,4,-.028f);
    ellipsoid(HEAD,(Vector3){.187f,.195f,.215f},(Color){18,123,71,255},12,6,-.195f);
    ellipsoid(BILL,(Vector3){.115f,.043f,.205f},(Color){227,177,47,255},8,4,-.032f);
    ellipsoid(BILL_LINE,(Vector3){.11f,.005f,.17f},(Color){87,62,24,255},8,4,-.005f);
    ellipsoid(NOSTRIL,(Vector3){.009f,.004f,.015f},(Color){62,53,29,255},6,4,-.004f);
    ellipsoid(EYE,(Vector3){.018f,.022f,.025f},(Color){7,12,9,255},8,4,-.022f);
    ellipsoid(GLINT,(Vector3){.004f,.005f,.006f},(Color){236,239,216,255},6,4,-.005f);
    ellipsoid(WING,(Vector3){.045f,.17f,.365f},(Color){126,129,120,255},10,5,-.17f);
    ellipsoid(WING_BAR,(Vector3){.009f,.067f,.17f},(Color){230,225,207,255},8,4,-.067f);
    ellipsoid(SPECULUM,(Vector3){.012f,.044f,.135f},(Color){45,67,161,255},8,4,-.044f);
    ellipsoid(TAIL,(Vector3){.17f,.045f,.21f},(Color){218,210,183,255},8,4,-.045f);
    ellipsoid(TAIL_CURL,(Vector3){.025f,.025f,.085f},(Color){29,37,33,255},6,4,-.025f);
    char vertex[512],fragment[512];if(!directory)directory="shaders";
    snprintf(vertex,sizeof(vertex),"%s/duck.vert",directory);snprintf(fragment,sizeof(fragment),"%s/duck.frag",directory);
    duckShader=LoadShader(vertex,fragment);
    shaded=duckShader.id && duckShader.id!=rlGetShaderIdDefault();
    if(shaded) {
        materialLoc=GetShaderLocation(duckShader,"duckMaterial");
        for(int i=0;i<PART_COUNT;i++)duckParts[i].materials[0].shader=duckShader;
    }
    floorModel=LoadModelFromMesh(GenMeshPlane(2.4f,3.5f,1,1));
    snprintf(fragment,sizeof(fragment),"%s/duck_floor.frag",directory);
    floorShader=LoadShader(vertex,fragment);
    floorShaded=floorShader.id && floorShader.id!=rlGetShaderIdDefault();
    if(floorShaded) {
        floorModel.materials[0].shader=floorShader;
        clockLoc=GetShaderLocation(floorShader,"duckTime");speedLoc=GetShaderLocation(floorShader,"swimSpeed");
        wakeLoc=GetShaderLocation(floorShader,"showWake");
    }
    ready=true;return true;
}

static bool valid_pose(DuckPose p) {
    return isfinite(p.position.x)&&isfinite(p.position.y)&&isfinite(p.position.z)&&
        isfinite(p.yaw)&&isfinite(p.time)&&isfinite(p.speed)&&isfinite(p.phase)&&isfinite(p.scale)&&p.scale>0;
}
static Matrix root_transform(DuckPose p) {
    Matrix m=MatrixScale(p.scale,p.scale,p.scale);
    m=MatrixMultiply(m,MatrixRotateY(p.yaw*DEG2RAD));
    return MatrixMultiply(m,MatrixTranslate(p.position.x,p.position.y,p.position.z));
}
static void part(int kind,Vector3 position,Vector3 rotation,Matrix parent,int material) {
    Matrix m=MatrixRotateXYZ(rotation);
    m=MatrixMultiply(m,MatrixTranslate(position.x,position.y,position.z));
    m=MatrixMultiply(m,parent);
    if(shaded)SetShaderValue(duckShader,materialLoc,&material,SHADER_UNIFORM_INT);
    DrawMesh(duckParts[kind].meshes[0],duckParts[kind].materials[0],m);
}

void duck_draw(DuckPose p) {
    if(!ready||!valid_pose(p))return;
    float swim=Clamp(p.speed/.8f,0,1),cycle=p.time*(2*PI*.75f)+p.phase;
    float roll=sinf(cycle)*.025f*swim;
    /* Bottom stays within 4 mm of the floor; no legs or visible foot skating. */
    float bob=.002f+.002f*sinf(cycle*2);
    Matrix root=root_transform(p);
    Matrix body=MatrixMultiply(MatrixRotateZ(roll),MatrixTranslate(0,bob,0));
    body=MatrixMultiply(body,root);
    part(BODY,(Vector3){0,.18f,0},(Vector3){0},body,0);
    part(BREAST,(Vector3){0,.235f,-.31f},(Vector3){0},body,0);
    part(NECK,(Vector3){0,.48f,-.36f},(Vector3){0},body,1);
    part(COLLAR,(Vector3){0,.475f,-.36f},(Vector3){0},body,0);
    Matrix head=MatrixMultiply(MatrixRotateZ(-roll*.8f),MatrixTranslate(0,.67f+.004f*sinf(cycle),-.39f));
    head=MatrixMultiply(head,body);
    part(HEAD,(Vector3){0},(Vector3){0},head,1);
    part(BILL,(Vector3){0,-.055f,-.29f},(Vector3){.04f,0,0},head,2);
    part(BILL_LINE,(Vector3){0,-.071f,-.31f},(Vector3){.04f,0,0},head,0);
    for(int s=-1;s<=1;s+=2) {
        part(NOSTRIL,(Vector3){s*.046f,-.019f,-.28f},(Vector3){0},head,0);
        part(EYE,(Vector3){s*.158f,.018f,-.096f},(Vector3){0,s*-.38f,0},head,2);
        part(GLINT,(Vector3){s*.173f,.025f,-.106f},(Vector3){0},head,0);
        float flex=.008f*sinf(cycle+s*.6f)*swim;
        Matrix wing=MatrixMultiply(MatrixRotateZ(s*flex),MatrixTranslate(s*.296f,.25f,.075f));
        wing=MatrixMultiply(wing,body);
        part(WING,(Vector3){0},(Vector3){-.12f,0,0},wing,0);
        part(WING_BAR,(Vector3){s*.039f,-.01f,.12f},(Vector3){-.22f,0,0},wing,0);
        part(SPECULUM,(Vector3){s*.047f,-.01f,.12f},(Vector3){-.22f,0,0},wing,1);
    }
    Matrix tail=MatrixMultiply(MatrixRotateY(sinf(cycle-.7f)*.06f*swim),MatrixTranslate(0,.24f,.53f));
    tail=MatrixMultiply(tail,body);
    part(TAIL,(Vector3){0},(Vector3){-.24f,0,0},tail,0);
    for(int s=-1;s<=1;s+=2) {
        part(TAIL_CURL,(Vector3){s*.029f,.04f,.055f},(Vector3){-.48f,0,0},tail,0);
        part(TAIL_CURL,(Vector3){s*.029f,.085f,.09f},(Vector3){.48f,0,0},tail,0);
    }
}

void duck_draw_floor(DuckPose p,bool wake) {
    if(!ready||!floorShaded||!valid_pose(p))return;
    float speed=Clamp(p.speed/.8f,0,1);int enabled=wake;
    float clock=p.time+p.phase/(2*PI*.75f);
    SetShaderValue(floorShader,clockLoc,&clock,SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader,speedLoc,&speed,SHADER_UNIFORM_FLOAT);
    SetShaderValue(floorShader,wakeLoc,&enabled,SHADER_UNIFORM_INT);
    Matrix m=MatrixMultiply(MatrixTranslate(0,.006f,.38f),root_transform(p));
    rlDrawRenderBatchActive();rlDisableDepthMask();
    DrawMesh(floorModel.meshes[0],floorModel.materials[0],m);
    rlDrawRenderBatchActive();rlEnableDepthMask();
}

int duck_renderer_triangle_count(void){return triangles;}
bool duck_renderer_has_shader(void){return shaded;}
void duck_renderer_shutdown(void) {
    if(!ready)return;
    for(int i=0;i<PART_COUNT;i++)UnloadModel(duckParts[i]);
    UnloadModel(floorModel);
    if(shaded)UnloadShader(duckShader);
    if(floorShaded)UnloadShader(floorShader);
    memset(duckParts,0,sizeof(duckParts));floorModel=(Model){0};
    ready=shaded=floorShaded=false;triangles=0;
}
