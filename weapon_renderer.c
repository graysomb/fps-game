#include "weapon_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <string.h>

enum { WR_BODY, WR_UPPER, WR_LOWER, WR_FLASH, WR_LAUNCH_BARREL, WR_LAUNCH_FRAME,
       WR_LAUNCH_GRIP, WR_LAUNCH_STOCK, WR_CUBE, WR_GLASS_LEFT, WR_GLASS_RIGHT,
       WR_GLASS_FRONT, WR_GLASS_BACK, WR_AURA, WR_PART_COUNT };
enum { WR_SILVER, WR_DARK, WR_ENERGY, WR_CHARCOAL, WR_PARTICLE, WR_GLASS, WR_GOLD_AURA };
static const Color wr_colors[] = { {184,190,201,255}, {37,43,58,255}, {255,255,255,255},
                                 {47,49,56,255}, {255,255,255,255}, {255,235,170,24}, {255,255,255,255} };
typedef struct WrBuilder {
    float vertices[4096 * 3], normals[4096 * 3], uv[4096 * 2];
    float edge_uv[4096 * 2];
    unsigned char colors[4096 * 4];
    int count;
    bool overflow;
} WrBuilder;
static Model wr_parts[WR_PART_COUNT];
static Shader wr_shader, wr_blur;
static int wr_eye, wr_energy, wr_strength, wr_mask, wr_direction, wr_time, wr_launcher;
static bool wr_ready;

static Vector3 wr_rotate(Vector3 p, float angle) {
    return (Vector3){p.x, p.y*cosf(angle)-p.z*sinf(angle), p.y*sinf(angle)+p.z*cosf(angle)};
}

static void wr_triangle(WrBuilder *b, Vector3 a, Vector3 c, Vector3 d, int material, float tint, int edge_mask) {
    if (b->count + 3 > 4096) { b->overflow=true; return; }
    Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(c,a), Vector3Subtract(d,a)));
    Vector3 points[3] = {a,c,d};
    Color color = wr_colors[material];
    for (int i=0; i<3; ++i) {
        int j = b->count++;
        memcpy(b->vertices+j*3, &points[i], sizeof(Vector3));
        memcpy(b->normals+j*3, &n, sizeof(Vector3));
        b->uv[j*2] = material == WR_ENERGY ? 1.0f : material == WR_PARTICLE ? 2.0f :
                         material == WR_GOLD_AURA ? 3.0f : 0.0f;
        b->uv[j*2+1] = (float)edge_mask;
        b->edge_uv[j*2] = i==0 ? 1.0f : 0.0f;
        b->edge_uv[j*2+1] = i==1 ? 1.0f : 0.0f;
        b->colors[j*4] = (unsigned char)(color.r*tint);
        b->colors[j*4+1] = (unsigned char)(color.g*tint);
        b->colors[j*4+2] = (unsigned char)(color.b*tint);
        b->colors[j*4+3] = color.a;
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

static void wr_quad(WrBuilder *b,Vector3 a,Vector3 c,Vector3 d,Vector3 e,int material,float tint) {
    wr_triangle(b,a,c,d,material,tint,5);
    wr_triangle(b,a,d,e,material,tint,3);
}

static Vector3 wr_octagon(int i,float w,float h,float bevel,float z) {
    float x=w*.5f,y=h*.5f;
    Vector3 points[8]={{-x+bevel,-y,z},{x-bevel,-y,z},{x,-y+bevel,z},{x,y-bevel,z},
                       {x-bevel,y,z},{-x+bevel,y,z},{-x,y-bevel,z},{-x,-y+bevel,z}};
    return points[i];
}

/* Real hollow geometry: end annuli, inner walls and beveled outer walls. */
static void wr_tube(WrBuilder *b,Vector3 center,float w,float h,float length,
                    float thickness,float bevel,int material) {
    Vector3 outer[4][8],inner[2][8];
    float zs[4]={-length*.5f,-length*.5f+bevel,length*.5f-bevel,length*.5f};
    for (int r=0;r<4;++r) for (int i=0;i<8;++i) {
        float inset=(r==0 || r==3) ? bevel*.45f : 0;
        outer[r][i]=Vector3Add(center,wr_octagon(i,w-2*inset,h-2*inset,bevel*.65f,zs[r]));
    }
    for (int r=0;r<2;++r) for (int i=0;i<8;++i)
        inner[r][i]=Vector3Add(center,wr_octagon(i,w-2*thickness,h-2*thickness,
                                                  bevel*.65f,zs[r==0 ? 0 : 3]));
    for (int i=0;i<8;++i) {
        int j=(i+1)%8;
        for (int r=0;r<3;++r)
            wr_quad(b,outer[r][i],outer[r][j],outer[r+1][j],outer[r+1][i],material,1);
        wr_quad(b,inner[0][j],inner[0][i],inner[1][i],inner[1][j],
                 material==WR_ENERGY ? WR_ENERGY : WR_CHARCOAL,1);
        wr_quad(b,outer[0][j],outer[0][i],inner[0][i],inner[0][j],material,1);
        wr_quad(b,outer[3][i],outer[3][j],inner[1][j],inner[1][i],material,1);
    }
}

static void wr_launcher_part(WrBuilder *b,int part) {
    if (part==WR_LAUNCH_BARREL) {
        wr_tube(b,(Vector3){0,0,-.83f},.82f,.82f,1.12f,.15f,.10f,WR_CHARCOAL);
        wr_tube(b,(Vector3){0,0,-1.37f},.86f,.86f,.18f,.12f,.065f,WR_SILVER);
        wr_tube(b,(Vector3){0,0,-1.30f},.61f,.61f,.06f,.055f,.02f,WR_ENERGY);
        wr_prism(b,(Vector3){0,0,-.30f},.84f,.84f,.16f,.06f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,.48f,-.83f},.45f,.18f,.56f,.045f,WR_CHARCOAL,-.05f);
        wr_prism(b,(Vector3){-.228f,.48f,-.99f},.018f,.075f,.19f,.009f,WR_ENERGY,0);
        wr_prism(b,(Vector3){0,-.49f,-.87f},.52f,.24f,.60f,.055f,WR_CHARCOAL,0);
    } else if (part==WR_LAUNCH_FRAME) {
        // Thin side frames joined by rails; no opaque wall behind the particles.
        wr_tube(b,(Vector3){0,0,-.37f},.98f,1.18f,.08f,.13f,.018f,WR_SILVER);
        wr_tube(b,(Vector3){0,0,.37f},.98f,1.18f,.08f,.13f,.018f,WR_SILVER);
        // Turn the chamber tube sideways: its openings face either side of the gun.
        for (int i=0;i<b->count;++i) {
            float x=b->vertices[i*3],z=b->vertices[i*3+2];
            b->vertices[i*3]=z; b->vertices[i*3+1]+=.23f; b->vertices[i*3+2]=.20f-x;
            x=b->normals[i*3]; z=b->normals[i*3+2];
            b->normals[i*3]=z; b->normals[i*3+2]=-x;
        }
        wr_prism(b,(Vector3){0,.775f,.20f},.66f,.09f,.80f,.025f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,-.315f,.20f},.66f,.09f,.80f,.025f,WR_CHARCOAL,0);
        for (int y=0;y<2;++y) for (int z=0;z<2;++z)
            wr_prism(b,(Vector3){0,.23f+(y==0 ? -.51f : .51f),.20f+(z==0 ? -.43f : .43f)},
                      .66f,.075f,.075f,.015f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,.64f,-.18f},.56f,.22f,.19f,.045f,WR_CHARCOAL,-.20f);
        wr_prism(b,(Vector3){-.285f,.64f,-.18f},.025f,.085f,.085f,.012f,WR_ENERGY,-.20f);
        wr_prism(b,(Vector3){0,-.41f,-.13f},.53f,.19f,.27f,.035f,WR_CHARCOAL,.10f);
        wr_prism(b,(Vector3){-.269f,-.41f,-.13f},.02f,.085f,.10f,.01f,WR_ENERGY,.10f);
    } else if (part==WR_LAUNCH_GRIP) {
        wr_prism(b,(Vector3){0,-.65f,.30f},.34f,.74f,.31f,.05f,WR_CHARCOAL,-.25f);
        wr_prism(b,(Vector3){0,-1.01f,.39f},.38f,.17f,.36f,.05f,WR_CHARCOAL,-.25f);
        wr_prism(b,(Vector3){0,-.50f,-.19f},.13f,.32f,.11f,.025f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,-.67f,-.025f},.13f,.10f,.41f,.025f,WR_SILVER,0);
        wr_prism(b,(Vector3){0,-.48f,-.02f},.09f,.22f,.09f,.022f,WR_ENERGY,-.2f);
    } else if (part==WR_LAUNCH_STOCK) {
        wr_prism(b,(Vector3){0,.06f,.65f},.36f,.34f,.34f,.05f,WR_CHARCOAL,0);
        wr_prism(b,(Vector3){0,-.12f,.82f},.48f,.77f,.16f,.075f,WR_CHARCOAL,0);
        wr_prism(b,(Vector3){0,.10f,.775f},.49f,.11f,.06f,.025f,WR_SILVER,0);
    } else if (part==WR_CUBE) {
        for (int face=0;face<6;++face) {
            Vector3 points[4]={{-.5f,-.5f,.5f},{.5f,-.5f,.5f},{.5f,.5f,.5f},{-.5f,.5f,.5f}};
            Matrix rot=face<4 ? MatrixRotateY(face*PI*.5f) : MatrixRotateX((face==4 ? 1 : -1)*PI*.5f);
            for (int i=0;i<4;++i) points[i]=Vector3Transform(points[i],rot);
            wr_quad(b,points[0],points[1],points[2],points[3],WR_PARTICLE,.84f+face*.025f);
        }
    } else {
        bool side_pane=part==WR_GLASS_LEFT || part==WR_GLASS_RIGHT;
        float side=(part==WR_GLASS_LEFT || part==WR_GLASS_FRONT) ? -1 : 1;
        Vector3 center=side_pane ? (Vector3){side*.411f,.23f,.20f} : (Vector3){0,.23f,.20f+side*.36f};
        for (int i=0;i<8;++i) {
            Vector3 a=wr_octagon(i,.72f,.92f,.049f,0),c=wr_octagon((i+1)%8,.72f,.92f,.049f,0);
            a=side_pane ? (Vector3){center.x,center.y+a.y,center.z-a.x} : Vector3Add(center,a);
            c=side_pane ? (Vector3){center.x,center.y+c.y,center.z-c.x} : Vector3Add(center,c);
            // Two-sided panes without changing the caller's culling state.
            wr_triangle(b,center,a,c,WR_GLASS,1,1);
            wr_triangle(b,center,c,a,WR_GLASS,1,1);
        }
    }
}

static void wr_build_part(WrBuilder *b, int part) {
    memset(b,0,sizeof(*b));
    if (part==WR_AURA) {
        // Parameter-space ribbons: the vertex shader orbits and ripples them.
        for (int strand=0;strand<2;++strand) for (int step=0;step<64;++step) {
            float t=step/64.0f,next=(step+1)/64.0f;
            float a=strand*PI+t*7.5f,c=strand*PI+next*7.5f;
            wr_quad(b,(Vector3){a,-1,t},(Vector3){a,1,t},
                      (Vector3){c,1,next},(Vector3){c,-1,next},WR_GOLD_AURA,1);
        }
        return;
    }
    if (part>=WR_LAUNCH_BARREL) { wr_launcher_part(b,part); return; }
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
        if (builder->overflow || !builder->count) {
            TraceLog(LOG_WARNING,"WEAPON: invalid procedural mesh %d",p);
            MemFree(builder);
            for (int i=0;i<p;++i) UnloadModel(wr_parts[i]);
            UnloadShader(wr_shader); UnloadShader(wr_blur); return false;
        }
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
    wr_time=GetShaderLocation(wr_shader,"effectTime");
    wr_launcher=GetShaderLocation(wr_shader,"launcher");
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

static float wr_ease(float t) {
    t=fmaxf(0,fminf(1,t));
    return t*t*(3-2*t);
}

WeaponPose weapon_pose(WeaponKind kind,Vector3 pos,Vector3 forward,Vector3 right,Vector3 up,
                       bool first_person,float aspect,const WeaponMeleePose *melee,float recoil) {
    bool swinging=melee && melee->progress>=0 && melee->progress<=1;
    float windup=0,extend=0,recover=0,anchor=0;
    if (swinging) {
        windup=wr_ease(melee->progress/WEAPON_MELEE_ACTIVE_START);
        extend=wr_ease((melee->progress-WEAPON_MELEE_ACTIVE_START)/
                       (WEAPON_MELEE_PEAK-WEAPON_MELEE_ACTIVE_START));
        recover=wr_ease((melee->progress-WEAPON_MELEE_ACTIVE_END)/
                        (1-WEAPON_MELEE_ACTIVE_END));
        anchor=windup*(1-recover);
        recoil*=1-anchor;
    }
    float scale=first_person ? .36f*fminf(1.0f,aspect) : .56f;
    float distance=(first_person ? 1.0f : .62f)-.045f*recoil;
    float right_offset=.48f;
    if (first_person) {
        // Keep the standard rear cap near the right edge. For the launcher,
        // place the stock's innermost corner beyond it, clipping the whole stock.
        // Both include the existing 22-degree inward gun rotation.
        float angle=22*DEG2RAD;
        bool launcher=kind==WEAPON_DYNAMIC_LAUNCHER;
        float rear_x=launcher ? -.48f*.5f : .58f*.5f-.065f*.45f;
        float rear_z=launcher ? .82f-.16f*.5f : .59f+.13f*.5f;
        float corner_x=scale*(rear_x*cosf(angle)+rear_z*sinf(angle));
        float corner_z=scale*(-rear_x*sinf(angle)+rear_z*cosf(angle));
        float edge=launcher ? 1.02f : .97f;
        right_offset=edge*tanf(30*DEG2RAD)*aspect*(distance-corner_z)-corner_x;
    }
    Vector3 offset=Vector3Add(Vector3Scale(right,right_offset),
                    Vector3Add(Vector3Scale(up,-.25f),
                               Vector3Scale(forward,distance)));
    pos=Vector3Add(pos,offset);
    // Mesh forward is -Z. right/up/forward form an orthonormal camera basis.
    Matrix m={right.x*scale,up.x*scale,-forward.x*scale,pos.x,
              right.y*scale,up.y*scale,-forward.y*scale,pos.y,
              right.z*scale,up.z*scale,-forward.z*scale,pos.z,0,0,0,1};
    if (first_person) m=MatrixMultiply(MatrixRotateY(22*(1-anchor)*DEG2RAD),m);
    // Center of the underside of the silver grip cap (which is tilted -0.25 rad).
    const Vector3 heel={0,-1.01f-.085f*cosf(.25f),.39f+.085f*sinf(.25f)};
    if (swinging) {
        // Rigid rotation around the grip; the heel turns forward for the strike.
        const Vector3 pivot={0,-.55f,.28f};
        // Tip past vertical in windup so the claws clear the top of the view.
        float windup_angle=kind==WEAPON_DYNAMIC_LAUNCHER ? 150 : 130;
        float pitch=(windup_angle*windup+(100-windup_angle)*extend)*(1-recover)*DEG2RAD;
        float roll=(-12*windup+6*extend)*(1-recover)*DEG2RAD;
        Matrix swing=MatrixMultiply(MatrixTranslate(-pivot.x,-pivot.y,-pivot.z),
                                     MatrixRotateX(pitch));
        swing=MatrixMultiply(swing,MatrixRotateZ(roll));
        swing=MatrixMultiply(swing,MatrixTranslate(pivot.x,pivot.y,pivot.z));
        m=MatrixMultiply(swing,m);
        Vector3 delta=Vector3Scale(Vector3Subtract(melee->strike_endpoint,
                                    Vector3Transform(heel,m)),anchor);
        m.m12+=delta.x; m.m13+=delta.y; m.m14+=delta.z;
        // Dip during the flip, then meet the unmodified sweep at active-start.
        if (windup<1) {
            Vector3 dip=Vector3Scale(up,-.50f*sinf(windup*PI));
            m.m12+=dip.x; m.m13+=dip.y; m.m14+=dip.z;
        }
    }
    Vector3 muzzle={0,0,kind==WEAPON_DYNAMIC_LAUNCHER ? -1.46f : -.84f};
    return (WeaponPose){m,Vector3Transform(muzzle,m),Vector3Transform(heel,m),kind};
}

Matrix weapon_cube_transform(int cube,float now) {
    if (cube<0 || cube>=WEAPON_LAUNCHER_CUBES) return MatrixIdentity();
    const float zs[3]={-.08f,.10f,-.05f};
    float phase=cube*2*PI/3,speed=.65f+.09f*cube;
    Vector3 pos={-.20f+.02f*sinf(now*speed+phase),
                 .23f+(1-cube)*.20f+.045f*sinf(now*(speed+.13f)+phase),
                 .20f+zs[cube]+.04f*cosf(now*(speed+.19f)+phase)};
    Matrix m=MatrixMultiply(MatrixScale(.19f,.19f,.19f),
        MatrixRotateXYZ((Vector3){now*(.25f+.06f*cube),now*(.31f+.04f*cube),phase+now*.18f}));
    return MatrixMultiply(m,MatrixTranslate(pos.x,pos.y,pos.z));
}

static void wr_draw_part(int part,Matrix transform) {
    DrawMesh(wr_parts[part].meshes[0],wr_parts[part].materials[0],transform);
}

void weapon_draw(WeaponPose pose,const WeaponVisual *v,Vector3 eye,float now,
                 bool gold,WeaponRenderPass pass) {
    if (!wr_ready) return;
    bool launcher=pose.kind==WEAPON_DYNAMIC_LAUNCHER;
    if (pass==WEAPON_GLASS && !launcher) return;
    float shot=weapon_shot_amount(v,now);
    Vector3 energy=launcher ? (Vector3){1,.85f,.035f} : gold ? (Vector3){1,.67f,.10f} : (Vector3){0,.95f,1};
    float strength=.82f+v->claw_open*(.34f+.12f*sinf(now*9)) + shot*.70f;
    int mask=pass==WEAPON_EMISSION ? 1 : 0;
    SetShaderValue(wr_shader,wr_eye,&eye,SHADER_UNIFORM_VEC3);
    SetShaderValue(wr_shader,wr_energy,&energy,SHADER_UNIFORM_VEC3);
    SetShaderValue(wr_shader,wr_strength,&strength,SHADER_UNIFORM_FLOAT);
    SetShaderValue(wr_shader,wr_mask,&mask,SHADER_UNIFORM_INT);
    SetShaderValue(wr_shader,wr_time,&now,SHADER_UNIFORM_FLOAT);
    int variant=launcher ? 1 : 0;
    SetShaderValue(wr_shader,wr_launcher,&variant,SHADER_UNIFORM_INT);
    if (pass==WEAPON_GLASS) {
        const Vector3 centers[4]={{-.411f,.23f,.20f},{.411f,.23f,.20f},{0,.23f,-.16f},{0,.23f,.56f}};
        int order[4]; float distances[4];
        for (int pane=0;pane<4;++pane) {
            float distance=Vector3DistanceSqr(eye,Vector3Transform(centers[pane],pose.transform));
            int at=pane;
            while (at>0 && distances[at-1]<distance) {
                order[at]=order[at-1]; distances[at]=distances[at-1]; --at;
            }
            order[at]=WR_GLASS_LEFT+pane; distances[at]=distance;
        }
        rlDrawRenderBatchActive();
        BeginBlendMode(BLEND_ALPHA); rlDisableDepthMask();
        for (int pane=0;pane<4;++pane) wr_draw_part(order[pane],pose.transform);
        rlDrawRenderBatchActive(); rlEnableDepthMask(); EndBlendMode();
        return;
    }
    if (launcher) {
        for (int part=WR_LAUNCH_BARREL;part<=WR_LAUNCH_STOCK;++part) wr_draw_part(part,pose.transform);
        for (int cube=0;cube<WEAPON_LAUNCHER_CUBES;++cube)
            wr_draw_part(WR_CUBE,MatrixMultiply(weapon_cube_transform(cube,now),pose.transform));
    } else {
        for (int i=0;i<WR_FLASH;++i) {
            Matrix transform=pose.transform;
            if (i!=WR_BODY) {
                float sign=i==WR_UPPER ? 1 : -1;
                Matrix local=MatrixMultiply(MatrixRotateX(sign*v->claw_open*10*DEG2RAD),
                                             MatrixTranslate(0,sign*.24f,-.32f));
                transform=MatrixMultiply(local,transform);
            }
            wr_draw_part(i,transform);
        }
    }
    if (shot>.65f) {
        // Small faceted muzzle burst shares the emissive shader and depth behavior.
        float size=shot*(launcher ? 1.55f : 1);
        Matrix tip=MatrixMultiply(MatrixScale(size,size,size),
                                  MatrixTranslate(0,0,launcher ? -1.49f : -.89f));
        tip=MatrixMultiply(tip,pose.transform);
        wr_draw_part(WR_FLASH,tip);
    }
    if (gold) {
        // Keep world occlusion, but do not write translucent ribbons into depth.
        rlDrawRenderBatchActive();
        BeginBlendMode(BLEND_ADDITIVE); rlDisableDepthMask(); rlDisableBackfaceCulling();
        wr_draw_part(WR_AURA,pose.transform);
        rlDrawRenderBatchActive();
        rlEnableBackfaceCulling(); rlEnableDepthMask(); EndBlendMode();
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
