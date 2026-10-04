#include "weapon_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <string.h>

enum { WR_BODY, WR_UPPER, WR_LOWER, WR_FLASH, WR_PART_COUNT };
enum { WR_SILVER, WR_DARK, WR_ENERGY };
static const Color wr_colors[] = { {184,190,201,255}, {37,43,58,255}, {255,255,255,255} };
typedef struct WrBuilder {
    float vertices[4096 * 3], normals[4096 * 3], uv[4096 * 2];
    float edge_uv[4096 * 2];
    unsigned char colors[4096 * 4];
    int count;
} WrBuilder;
static Model wr_parts[WR_PART_COUNT];
static Shader wr_shader, wr_blur;
static int wr_eye, wr_energy, wr_strength, wr_mask, wr_direction;
static bool wr_ready;

static Vector3 wr_rotate(Vector3 p, float angle) {
    return (Vector3){p.x, p.y*cosf(angle)-p.z*sinf(angle), p.y*sinf(angle)+p.z*cosf(angle)};
}

static void wr_triangle(WrBuilder *b, Vector3 a, Vector3 c, Vector3 d, int material, float tint, int edge_mask) {
    if (b->count + 3 > 4096) return;
    Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(c,a), Vector3Subtract(d,a)));
    Vector3 points[3] = {a,c,d};
    Color color = wr_colors[material];
    for (int i=0; i<3; ++i) {
        int j = b->count++;
        memcpy(b->vertices+j*3, &points[i], sizeof(Vector3));
        memcpy(b->normals+j*3, &n, sizeof(Vector3));
        b->uv[j*2] = material == WR_ENERGY ? 1.0f : 0.0f;
        b->uv[j*2+1] = (float)edge_mask;
        b->edge_uv[j*2] = i==0 ? 1.0f : 0.0f;
        b->edge_uv[j*2+1] = i==1 ? 1.0f : 0.0f;
        b->colors[j*4] = (unsigned char)(color.r*tint);
        b->colors[j*4+1] = (unsigned char)(color.g*tint);
        b->colors[j*4+2] = (unsigned char)(color.b*tint);
        b->colors[j*4+3] = 255;
    }
}

/* Octagonal cross sections, with an additional bevel around both end caps. */
static void wr_prism(WrBuilder *b, Vector3 center, float w, float h, float length,
                     float bevel, int material, float angle) {
    Vector3 rings[4][8];
    float zs[4] = {-length*.5f, -length*.5f+bevel, length*.5f-bevel, length*.5f};
    for (int r=0; r<4; ++r) {
        float inset = (r==0 || r==3) ? bevel*.45f : 0;
        float x=w*.5f-inset, y=h*.5f-inset, k=bevel*.65f;
        Vector3 points[8] = {{-x+k,-y,0},{x-k,-y,0},{x,-y+k,0},{x,y-k,0},
                              {x-k,y,0},{-x+k,y,0},{-x,y-k,0},{-x,-y+k,0}};
        for (int i=0;i<8;++i) {
            points[i].z=zs[r];
            rings[r][i]=Vector3Add(center,wr_rotate(points[i],angle));
        }
    }
    for (int r=0;r<3;++r) for (int i=0;i<8;++i) {
        int j=(i+1)%8;
        float tint=material==WR_ENERGY ? (0.84f+0.02f*i) : 1.0f;
        // Mark actual face boundaries, excluding the quad's triangulation diagonal.
        wr_triangle(b,rings[r][i],rings[r][j],rings[r+1][j],material,tint,5);
        wr_triangle(b,rings[r][i],rings[r+1][j],rings[r+1][i],material,tint,3);
    }
    Vector3 front=Vector3Add(center,wr_rotate((Vector3){0,0,zs[0]},angle));
    Vector3 back=Vector3Add(center,wr_rotate((Vector3){0,0,zs[3]},angle));
    for (int i=0;i<8;++i) {
        int j=(i+1)%8;
        // End caps have a perimeter only; hide the fan's internal spokes.
        wr_triangle(b,front,rings[0][j],rings[0][i],material,1,1);
        wr_triangle(b,back,rings[3][i],rings[3][j],material,1,1);
    }
}

static void wr_build_part(WrBuilder *b, int part) {
    memset(b,0,sizeof(*b));
    if (part==WR_BODY) {
        wr_prism(b,(Vector3){0,0,0},.65f,.65f,.55f,.10f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,0,-.32f},.54f,.54f,.14f,.055f,WR_DARK,0);
        wr_prism(b,(Vector3){0,0,-.39f},.64f,.64f,.10f,.055f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,0,-.49f},.40f,.40f,.17f,.045f,WR_ENERGY,0);
        wr_prism(b,(Vector3){0,0,-.60f},.35f,.35f,.14f,.06f,WR_DARK,0);
        wr_prism(b,(Vector3){0,0,-.72f},.17f,.17f,.15f,.035f,WR_DARK,0);
        wr_prism(b,(Vector3){0,0,-.81f},.17f,.17f,.06f,.025f,WR_ENERGY,0);
        // Rear energy cartridge: exposed sides bounded by real silver rails.
        wr_prism(b,(Vector3){0,0,.37f},.38f,.32f,.34f,.065f,WR_ENERGY,0);
        wr_prism(b,(Vector3){0,.24f,.36f},.58f,.13f,.44f,.045f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,-.24f,.36f},.58f,.13f,.44f,.045f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,0,.59f},.58f,.56f,.13f,.065f,WR_SILVER,0);
        // Carry handle (three bars), grip, and trigger guard (three bars).
        wr_prism(b,(Vector3){0,.48f,.12f},.15f,.13f,.38f,.025f,WR_DARK,0);
        wr_prism(b,(Vector3){0,.39f,-.045f},.15f,.24f,.10f,.025f,WR_DARK,-.3f);
        wr_prism(b,(Vector3){0,.39f,.285f},.15f,.24f,.10f,.025f,WR_DARK,.3f);
        wr_prism(b,(Vector3){0,-.65f,.30f},.30f,.74f,.28f,.05f,WR_DARK,-.25f);
        wr_prism(b,(Vector3){0,-1.01f,.39f},.35f,.17f,.33f,.05f,WR_SILVER,-.25f);
        wr_prism(b,(Vector3){0,-.50f,-.19f},.12f,.32f,.11f,.025f,WR_DARK,0);
        wr_prism(b,(Vector3){0,-.67f,-.025f},.12f,.10f,.41f,.025f,WR_DARK,0);
        wr_prism(b,(Vector3){0,-.48f,-.02f},.09f,.22f,.09f,.022f,WR_ENERGY,-.2f);
    } else if (part==WR_FLASH) {
        wr_prism(b,(Vector3){0},.25f,.25f,.12f,.035f,WR_ENERGY,0);
    } else {
        float sign=part==WR_UPPER ? 1.0f : -1.0f;
        // Coordinates relative to the hinge on the front housing.
        wr_prism(b,(Vector3){0,sign*.10f,-.15f},.21f,.18f,.37f,.035f,WR_DARK,sign*.42f);
        wr_prism(b,(Vector3){0,sign*.24f,-.44f},.25f,.18f,.37f,.045f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,sign*.18f,-.65f},.25f,.18f,.22f,.04f,WR_DARK,-sign*.65f);
        wr_prism(b,(Vector3){0,sign*.13f,-.75f},.23f,.16f,.065f,.025f,WR_ENERGY,-sign*.65f);
    }
}

bool weapon_renderer_init(void) {
    if (wr_ready) return true;
    wr_shader=LoadShader("shaders/weapon.vert","shaders/weapon.frag");
    wr_blur=LoadShader(NULL,"shaders/weapon_blur.frag");
    unsigned int default_shader=rlGetShaderIdDefault();
    if (!wr_shader.id || !wr_blur.id || wr_shader.id==default_shader || wr_blur.id==default_shader) {
        if (wr_shader.id && wr_shader.id!=default_shader) UnloadShader(wr_shader);
        if (wr_blur.id && wr_blur.id!=default_shader) UnloadShader(wr_blur);
        wr_shader=(Shader){0}; wr_blur=(Shader){0};
        TraceLog(LOG_WARNING,"WEAPON: shaders unavailable; skipping weapon rendering");
        return false;
    }
    WrBuilder *builder=MemAlloc(sizeof(WrBuilder));
    if (!builder) { UnloadShader(wr_shader); UnloadShader(wr_blur); return false; }
    for (int p=0;p<WR_PART_COUNT;++p) {
        wr_build_part(builder,p);
        Mesh mesh={0}; mesh.vertexCount=builder->count; mesh.triangleCount=builder->count/3;
        mesh.vertices=MemAlloc(mesh.vertexCount*3*sizeof(float));
        mesh.normals=MemAlloc(mesh.vertexCount*3*sizeof(float));
        mesh.texcoords=MemAlloc(mesh.vertexCount*2*sizeof(float));
        mesh.texcoords2=MemAlloc(mesh.vertexCount*2*sizeof(float));
        mesh.colors=MemAlloc(mesh.vertexCount*4);
        if (!mesh.vertices || !mesh.normals || !mesh.texcoords || !mesh.texcoords2 || !mesh.colors) {
            UnloadMesh(mesh); MemFree(builder);
            for (int i=0;i<p;++i) UnloadModel(wr_parts[i]);
            UnloadShader(wr_shader); UnloadShader(wr_blur); return false;
        }
        memcpy(mesh.vertices,builder->vertices,mesh.vertexCount*3*sizeof(float));
        memcpy(mesh.normals,builder->normals,mesh.vertexCount*3*sizeof(float));
        memcpy(mesh.texcoords,builder->uv,mesh.vertexCount*2*sizeof(float));
        memcpy(mesh.texcoords2,builder->edge_uv,mesh.vertexCount*2*sizeof(float));
        memcpy(mesh.colors,builder->colors,mesh.vertexCount*4);
        UploadMesh(&mesh,false);
        wr_parts[p]=LoadModelFromMesh(mesh);
        wr_parts[p].materials[0].shader=wr_shader;
    }
    MemFree(builder);
    wr_eye=GetShaderLocation(wr_shader,"eyePosition");
    wr_energy=GetShaderLocation(wr_shader,"energyColor");
    wr_strength=GetShaderLocation(wr_shader,"energyStrength");
    wr_mask=GetShaderLocation(wr_shader,"emissionOnly");
    wr_direction=GetShaderLocation(wr_blur,"direction");
    wr_ready=true; return true;
}

void weapon_renderer_shutdown(void) {
    if (!wr_ready) return;
    for (int i=0;i<WR_PART_COUNT;++i) UnloadModel(wr_parts[i]);
    UnloadShader(wr_shader); UnloadShader(wr_blur);
    memset(wr_parts,0,sizeof(wr_parts)); wr_ready=false;
}

void weapon_visual_reset(WeaponVisual *v) { *v=(WeaponVisual){.shot_time=-1000}; }
void weapon_visual_shot(WeaponVisual *v,float now) { ++v->shot_sequence; v->shot_time=now; }
void weapon_visual_receive(WeaponVisual *v,uint32_t sequence,uint8_t age_ms,float now) {
    // Serial arithmetic also handles sequence wrap. Ignore old/duplicate snapshots.
    uint32_t advance=sequence-v->shot_sequence;
    if (!advance || advance>=0x80000000u) return;
    v->shot_sequence=sequence;
    v->shot_time=age_ms<120 ? now-(float)age_ms*.001f : -1000;
}
void weapon_visual_update(WeaponVisual *v,bool tether,float dt) {
    float step=fmaxf(0,fminf(dt, .15f))/.15f;
    v->claw_open=tether ? fminf(1,v->claw_open+step) : fmaxf(0,v->claw_open-step);
}
float weapon_shot_amount(const WeaponVisual *v,float now) {
    float t=now-v->shot_time;
    return t>=0 && t<WEAPON_SHOT_SECONDS ? 1-t/WEAPON_SHOT_SECONDS : 0;
}

WeaponPose weapon_pose(Vector3 pos,Vector3 forward,Vector3 right,Vector3 up,
                       bool first_person,float aspect,float melee,float recoil) {
    float scale=first_person ? .36f*fminf(1.0f,aspect) : .56f;
    float distance=(first_person ? 1.0f : .62f)-.045f*recoil;
    float right_offset=.48f;
    if (first_person) {
        // Anchor the rear cap at 97% of the viewport's right half-width.
        // Its outer corner includes the existing 22-degree inward gun rotation.
        float angle=22*DEG2RAD;
        float rear_x=.58f*.5f-.065f*.45f, rear_z=.59f+.13f*.5f;
        float corner_x=scale*(rear_x*cosf(angle)+rear_z*sinf(angle));
        float corner_z=scale*(-rear_x*sinf(angle)+rear_z*cosf(angle));
        right_offset=.97f*tanf(30*DEG2RAD)*aspect*(distance-corner_z)-corner_x;
    }
    Vector3 offset=Vector3Add(Vector3Scale(right,right_offset),
                    Vector3Add(Vector3Scale(up,-.25f-.23f*melee),
                               Vector3Scale(forward,distance)));
    pos=Vector3Add(pos,offset);
    // Mesh forward is -Z. right/up/forward form an orthonormal camera basis.
    Matrix m={right.x*scale,up.x*scale,-forward.x*scale,pos.x,
              right.y*scale,up.y*scale,-forward.y*scale,pos.y,
              right.z*scale,up.z*scale,-forward.z*scale,pos.z,0,0,0,1};
    if (first_person) m=MatrixMultiply(MatrixRotateY(22*DEG2RAD),m);
    return (WeaponPose){m,Vector3Transform((Vector3){0,0,-.84f},m)};
}

void weapon_draw(WeaponPose pose,const WeaponVisual *v,Vector3 eye,float now,
                 bool gold,bool powered_shot,bool emission_only) {
    if (!wr_ready) return;
    float shot=weapon_shot_amount(v,now);
    Vector3 energy=(gold || (powered_shot && shot>0)) ? (Vector3){1,.67f,.10f} : (Vector3){0,.95f,1};
    float strength=.82f+v->claw_open*(.34f+.12f*sinf(now*9)) + shot*.70f;
    int mask=emission_only ? 1 : 0;
    SetShaderValue(wr_shader,wr_eye,&eye,SHADER_UNIFORM_VEC3);
    SetShaderValue(wr_shader,wr_energy,&energy,SHADER_UNIFORM_VEC3);
    SetShaderValue(wr_shader,wr_strength,&strength,SHADER_UNIFORM_FLOAT);
    SetShaderValue(wr_shader,wr_mask,&mask,SHADER_UNIFORM_INT);
    for (int i=0;i<WR_FLASH;++i) {
        Matrix transform=pose.transform;
        if (i!=WR_BODY) {
            float sign=i==WR_UPPER ? 1 : -1;
            Matrix local=MatrixMultiply(MatrixRotateX(sign*v->claw_open*10*DEG2RAD),
                                         MatrixTranslate(0,sign*.24f,-.32f));
            transform=MatrixMultiply(local,transform);
        }
        DrawMesh(wr_parts[i].meshes[0],wr_parts[i].materials[0],transform);
    }
    if (shot>.65f) {
        // Small faceted muzzle burst shares the emissive shader and depth behavior.
        Matrix tip=MatrixMultiply(MatrixScale(shot,shot,shot),MatrixTranslate(0,0,-.89f));
        tip=MatrixMultiply(tip,pose.transform);
        DrawMesh(wr_parts[WR_FLASH].meshes[0],wr_parts[WR_FLASH].materials[0],tip);
    }
}

void weapon_bloom_unload(WeaponBloom *b) {
    if (b->emission.id) UnloadRenderTexture(b->emission);
    for (int i=0;i<2;++i) if (b->blur[i].id) UnloadRenderTexture(b->blur[i]);
    *b=(WeaponBloom){0};
}
bool weapon_bloom_resize(WeaponBloom *b,int width,int height) {
    if (!wr_ready) return false;
    if (b->width==width && b->height==height && b->emission.id) return true;
    weapon_bloom_unload(b);
    b->emission=LoadRenderTexture(width,height);
    int w=width/2>0 ? width/2 : 1, h=height/2>0 ? height/2 : 1;
    for (int i=0;i<2;++i) { b->blur[i]=LoadRenderTexture(w,h); SetTextureFilter(b->blur[i].texture,TEXTURE_FILTER_BILINEAR); }
    SetTextureFilter(b->emission.texture,TEXTURE_FILTER_BILINEAR);
    if (!b->emission.id || !b->blur[0].id || !b->blur[1].id) { weapon_bloom_unload(b); return false; }
    b->width=width; b->height=height; return true;
}
void weapon_clear_depth(void) {
    rlDrawRenderBatchActive();
    rlEnableDepthMask(); rlColorMask(false,false,false,false);
    rlClearScreenBuffers(); rlColorMask(true,true,true,true);
}
void weapon_bloom_begin(WeaponBloom *b,RenderTexture2D world) {
    BeginTextureMode(b->emission); ClearBackground(BLACK);
    rlDrawRenderBatchActive();
    rlBindFramebuffer(RL_READ_FRAMEBUFFER,world.id);
    rlBindFramebuffer(RL_DRAW_FRAMEBUFFER,b->emission.id);
    rlBlitFramebuffer(0,0,b->width,b->height,0,0,b->width,b->height,0x00000100); // depth
    rlEnableFramebuffer(b->emission.id);
}
static void wr_fullscreen(Texture2D texture,int width,int height) {
    DrawTexturePro(texture,(Rectangle){0,0,(float)texture.width,-(float)texture.height},
                    (Rectangle){0,0,(float)width,(float)height},(Vector2){0},0,WHITE);
}
void weapon_bloom_composite(WeaponBloom *b,RenderTexture2D scene) {
    EndTextureMode();
    for (int pass=0;pass<2;++pass) {
        BeginTextureMode(b->blur[pass]); ClearBackground(BLACK);
        Vector2 direction=pass==0 ? (Vector2){6.0f/b->width,0} : (Vector2){0,3.0f/b->blur[0].texture.height};
        SetShaderValue(wr_blur,wr_direction,&direction,SHADER_UNIFORM_VEC2);
        BeginShaderMode(wr_blur);
        wr_fullscreen(pass==0 ? b->emission.texture : b->blur[0].texture,
                       b->blur[pass].texture.width,b->blur[pass].texture.height);
        EndShaderMode(); EndTextureMode();
    }
    BeginTextureMode(scene);
    BeginBlendMode(BLEND_ADDITIVE);
    DrawTexturePro(b->blur[1].texture,(Rectangle){0,0,(float)b->blur[1].texture.width,-(float)b->blur[1].texture.height},
                   (Rectangle){0,0,(float)b->width,(float)b->height},(Vector2){0},0,(Color){110,110,110,255});
    EndBlendMode();
}
