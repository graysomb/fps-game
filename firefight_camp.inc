/* Cosmetic intermission camp: no voxels, collision bodies, or gameplay effects. */
typedef struct { bool active; Vector3 pos; float age; } CampState;
static CampState camp;
static Model campSolid, campFlames, campSmoke;
static Shader campShader;
static bool campResources;
static int campTimeLoc, campModeLoc, campEmissionLoc;
/* Keep the stew below a standing rat's eyes. */
static const float campHeightScale=.75f;
static const float campScale=.5f;
typedef struct { Mesh mesh; int count; } CampMesh;
#define CAMP_VERTICES 18000
static CampMesh camp_mesh_begin(void) {
    CampMesh b={0};
    b.mesh.vertices=MemAlloc(CAMP_VERTICES*3*sizeof(float));
    b.mesh.normals=MemAlloc(CAMP_VERTICES*3*sizeof(float));
    b.mesh.texcoords=MemAlloc(CAMP_VERTICES*2*sizeof(float));
    b.mesh.colors=MemAlloc(CAMP_VERTICES*4);
    return b;
}
static void camp_triangle(CampMesh *b,Vector3 a,Vector3 c,Vector3 d,Color color) {
    if(b->count+3>CAMP_VERTICES)return;
    Vector3 normal=Vector3Normalize(Vector3CrossProduct(v_sub(c,a),v_sub(d,a)));
    Vector3 p[3]={a,c,d};float shade=.90f+.025f*((b->count/3)%5);
    for(int k=0;k<3;k++) {
        int n=b->count++;memcpy(&b->mesh.vertices[n*3],&p[k],sizeof(Vector3));
        memcpy(&b->mesh.normals[n*3],&normal,sizeof(Vector3));
        b->mesh.texcoords[n*2]=p[k].x*.5f+.5f;b->mesh.texcoords[n*2+1]=p[k].y;
        b->mesh.colors[n*4]=(unsigned char)(color.r*shade);b->mesh.colors[n*4+1]=(unsigned char)(color.g*shade);
        b->mesh.colors[n*4+2]=(unsigned char)(color.b*shade);b->mesh.colors[n*4+3]=color.a;
    }
}
static void camp_tube(CampMesh *b,Vector3 a,Vector3 c,float radius,Color color) {
    Vector3 axis=v_norm(v_sub(c,a)),right=v_norm(v_cross(axis,fabsf(axis.y)>.9f?(Vector3){1,0,0}:(Vector3){0,1,0})),up=v_cross(axis,right);
    for(int n=0;n<7;n++){
        float t=n*2*PI/7,u=(n+1)*2*PI/7;
        Vector3 r=v_mul(v_add(v_mul(right,cosf(t)),v_mul(up,sinf(t))),radius),s=v_mul(v_add(v_mul(right,cosf(u)),v_mul(up,sinf(u))),radius);
        camp_triangle(b,v_add(a,r),v_add(c,r),v_add(c,s),color);camp_triangle(b,v_add(a,r),v_add(c,s),v_add(a,s),color);
        camp_triangle(b,a,v_add(a,r),v_add(a,s),color);camp_triangle(b,c,v_add(c,s),v_add(c,r),color);
    }
}
static void camp_band(CampMesh *b,float y0,float r0,float y1,float r1,Color color){
    for(int n=0;n<16;n++){
        float a=n*2*PI/16,c=(n+1)*2*PI/16;
        Vector3 p={cosf(a)*r0,y0,sinf(a)*r0},q={cosf(c)*r0,y0,sinf(c)*r0};
        Vector3 r={cosf(a)*r1,y1,sinf(a)*r1},s={cosf(c)*r1,y1,sinf(c)*r1};
        camp_triangle(b,p,r,s,color);camp_triangle(b,p,s,q,color);
    }
}
static void camp_stone(CampMesh *b,Vector3 p,Vector3 scale,Color color){
    for(int y=0;y<4;y++)for(int n=0;n<7;n++){
        Vector3 v[4];for(int k=0;k<4;k++){
            float lat=-PI/2+(y+(k/2))*PI/4,angle=(n+k%2)*2*PI/7;
            v[k]=v_add(p,(Vector3){cosf(lat)*cosf(angle)*scale.x,sinf(lat)*scale.y,cosf(lat)*sinf(angle)*scale.z});
        }
        camp_triangle(b,v[0],v[2],v[3],color);camp_triangle(b,v[0],v[3],v[1],color);
    }
}
static Model camp_mesh_end(CampMesh b){b.mesh.vertexCount=b.count;b.mesh.triangleCount=b.count/3;UploadMesh(&b.mesh,false);return LoadModelFromMesh(b.mesh);}
static void camp_init(void){
    if(campResources)return;campResources=true;
    CampMesh b=camp_mesh_begin();Color iron={70,73,77,255};
    camp_band(&b,.53f,0,.57f,.35f,iron);camp_band(&b,.57f,.35f,.78f,.53f,iron);
    camp_band(&b,.78f,.53f,1.16f,.60f,iron);camp_band(&b,1.16f,.60f,1.23f,.64f,iron);
    camp_band(&b,1.23f,.64f,1.23f,.55f,(Color){95,99,103,255});camp_band(&b,1.23f,.55f,1.12f,.53f,iron);
    camp_band(&b,1.19f,.55f,1.19f,0,(Color){176,91,22,255});
    for(int side=-1;side<=1;side+=2){
        Vector3 a={side*.58f,1.14f,-.13f},c={side*.80f,1.13f,-.13f},d={side*.80f,1.13f,.13f},e={side*.58f,1.14f,.13f};
        camp_tube(&b,a,c,.045f,iron);camp_tube(&b,c,d,.045f,iron);camp_tube(&b,d,e,.045f,iron);
    }
    camp_band(&b,.50f,.44f,.56f,.44f,iron);
    for(int n=0;n<3;n++){
        float a=n*2*PI/3;camp_tube(&b,(Vector3){cosf(a)*.65f,.08f,sinf(a)*.65f},(Vector3){cosf(a)*.44f,.55f,sinf(a)*.44f},.045f,iron);
        Vector3 end={cosf(a)*.53f,.14f,sinf(a)*.53f};camp_tube(&b,end,(Vector3){-end.x,.12f,-end.z},.085f,(Color){65,36,20,255});
    }
    for(int n=0;n<11;n++){float a=n*2*PI/11;camp_stone(&b,(Vector3){cosf(a)*.78f,.12f,sinf(a)*.78f},(Vector3){.15f,.12f,.14f},(Color){109,111,117,255});}
    for(int n=0;n<20;n++){
        float a=n*2.39996f,r=.46f*sqrtf((n+.5f)/20);Vector3 p={cosf(a)*r,1.195f,sinf(a)*r};
        if(n%3==0)camp_stone(&b,p,(Vector3){.065f,.035f,.06f},(Color){238,131,30,255});
        else{Color green=n%2?(Color){65,120,21,255}:(Color){116,166,34,255};
            Vector3 tip=v_add(p,(Vector3){cosf(a)*.15f,.025f,sinf(a)*.15f});
            Vector3 left=v_add(p,(Vector3){-sinf(a)*.085f,0,cosf(a)*.085f}),right=v_add(p,(Vector3){sinf(a)*.085f,0,-cosf(a)*.085f});
            Vector3 mid=v_add(p,(Vector3){cosf(a)*.04f,.045f,sinf(a)*.04f});
            camp_triangle(&b,left,tip,mid,green);camp_triangle(&b,mid,tip,right,green);camp_triangle(&b,left,mid,right,green);
        }
    }
    campSolid=camp_mesh_end(b);b=camp_mesh_begin();
    for(int n=0;n<9;n++){
        float a=n*2.4f,r=n?.27f:0;Vector3 p={cosf(a)*r,.15f,sinf(a)*r};
        for(int k=0;k<5;k++){
            float t=k*2*PI/5,u=(k+1)*2*PI/5;
            camp_triangle(&b,v_add(p,(Vector3){cosf(t)*.14f,0,sinf(t)*.14f}),v_add(p,(Vector3){.04f,.35f+.07f*(n%3),0}),v_add(p,(Vector3){cosf(u)*.14f,0,sinf(u)*.14f}),(Color){255,163,28,255});
        }
    }
    campFlames=camp_mesh_end(b);b=camp_mesh_begin();
    for(int n=0;n<12;n++){float y=n/12.f,z=(n+1)/12.f;camp_triangle(&b,(Vector3){-1,y,0},(Vector3){1,y,0},(Vector3){1,z,0},WHITE);camp_triangle(&b,(Vector3){-1,y,0},(Vector3){1,z,0},(Vector3){-1,z,0},WHITE);}
    campSmoke=camp_mesh_end(b);
    campShader=LoadShader("shaders/camp.vert","shaders/camp.frag");
    if(campShader.id && campShader.id!=rlGetShaderIdDefault()){
        campTimeLoc=GetShaderLocation(campShader,"campTime");campModeLoc=GetShaderLocation(campShader,"campMode");campEmissionLoc=GetShaderLocation(campShader,"emissionPass");
        campSolid.materials[0].shader=campFlames.materials[0].shader=campSmoke.materials[0].shader=campShader;
    }
}
static void camp_end(void){camp=(CampState){0};}
static bool camp_ground(Vector3 p,float *height){
    p=world_wrap(p);if(!worldPeriodic && (fabsf(p.x)>FLOOR_SIZE-1 || fabsf(p.z)>FLOOR_SIZE-1))return false;
    Ray ray={{p.x,32,p.z},{0,-1,0}};VoxelHit hit;*height=0;
    if(first_voxel_hit_detailed(ray,32,-1,true,&hit))*height=32-hit.t;
    return true;
}
static bool camp_site(Vector3 center,Vector3 *standing){
    float ground;if(!camp_ground(center,&ground))return false;center.y=ground;
    for(int n=0;n<9;n++){
        Vector3 at=v_add(center,(Vector3){n?cosf(n*PI/4)*.9f:0,0,n?sinf(n*PI/4)*.9f:0});float h;
        if(!camp_ground(at,&h)||fabsf(h-ground)>.08f)return false;
        at.y=ground+BASE_EYE_HEIGHT+.04f;if(!spawn_position_clear_ex(at,.15f,0,0))return false;
    }
    for(int i=0;i<firefightHumans;i++){
        float angle=i*2*PI/firefightHumans;
        Vector3 at=world_wrap(v_add(center,(Vector3){sinf(angle)*2.5f,0,cosf(angle)*2.5f}));float h;
        if(!camp_ground(at,&h)||fabsf(h-ground)>.25f)return false;
        at.y=h+BASE_EYE_HEIGHT+.04f;
        if(!spawn_position_clear_ex(at,spawn_clear_half(i),0,i))return false;
        standing[i]=at;
    }
    camp.pos=center;return true;
}
static void camp_begin(void){
    camp_end();if(firefightPractice || enemyLab || netTransport.role!=NET_ROLE_OFFLINE)return;
    Vector3 center={0};if(!firefightFoundry){center=players[0].pos;Vector3 sum={0};for(int i=0;i<firefightHumans;i++)sum=v_add(sum,world_delta(players[i].pos,center));center=world_wrap(v_add(center,v_mul(sum,1.f/firefightHumans)));}
    Vector3 standing[MAX_PLAYERS];
    for(int n=0;n<=160;n++){
        float angle=(n?((n-1)%16):0)*2*PI/16,radius=n?2*(1+(n-1)/16):0;
        Vector3 at=world_wrap(v_add(center,(Vector3){cosf(angle)*radius,0,sinf(angle)*radius}));
        if(!camp_site(at,standing))continue;
        camp.active=true;
        for(int i=0;i<firefightHumans;i++){
            players[i].pos=players[i].death_pos=standing[i];players[i].vel=(Vector3){0};
            Vector3 dir=world_delta(camp.pos,standing[i]);players[i].yaw=atan2f(-dir.x,-dir.z)*RAD2DEG;players[i].pitch=0;
            characterAnimations[i].initialized=false;aim_cancel(i);
        }
        return;
    }
}
static bool camp_visible(void){return camp.active && gameMode==GAME_MODE_FIREFIGHT && firefightWaveState==FIREFIGHT_STATE_INTERMISSION && (gameState==GAME_STATE_PLAYING || gameState==GAME_STATE_PAUSED) && !enemyLab && !firefightPractice;}
static void camp_update(float dt){if(camp_visible()){if(gameState==GAME_STATE_PLAYING)camp.age+=fmaxf(0,dt);}else if(camp.active)camp_end();}
static void camp_draw(Camera3D camera,bool emission){
    if(!camp_visible() || !periodic_visible_box(camera,v_add(camp.pos,(Vector3){0,1.4f*campScale,0}),v_mul((Vector3){1,1.5f,1},campScale)))return;
    camp_init();bool shader=campShader.id && campShader.id!=rlGetShaderIdDefault();int pass=emission?1:0,mode=0;
    if(shader){SetShaderValue(campShader,campTimeLoc,&camp.age,SHADER_UNIFORM_FLOAT);SetShaderValue(campShader,campEmissionLoc,&pass,SHADER_UNIFORM_INT);SetShaderValue(campShader,campModeLoc,&mode,SHADER_UNIFORM_INT);}
    if(!emission)DrawModelEx(campSolid,camp.pos,(Vector3){0,1,0},0,v_mul((Vector3){1,campHeightScale,1},campScale),WHITE);
    mode=1;if(shader)SetShaderValue(campShader,campModeLoc,&mode,SHADER_UNIFORM_INT);
    rlDisableBackfaceCulling();DrawModelEx(campFlames,camp.pos,(Vector3){0,1,0},0,v_mul((Vector3){1,campHeightScale*(shader?1:.93f+.12f*sinf(camp.age*8)),1},campScale),WHITE);
    if(!emission){
        mode=2;if(shader)SetShaderValue(campShader,campModeLoc,&mode,SHADER_UNIFORM_INT);
        Vector3 origins[3];int order[3]={0,1,2};float distances[3];
        for(int n=0;n<3;n++){origins[n]=v_add(camp.pos,v_mul((Vector3){(n-1)*.13f,1.16f*campHeightScale,.08f*sinf(n*2.f)},campScale));distances[n]=Vector3DistanceSqr(origins[n],camera.position);}
        for(int a=0;a<3;a++)for(int b=a+1;b<3;b++)if(distances[order[a]]<distances[order[b]]){int t=order[a];order[a]=order[b];order[b]=t;}
        rlDrawRenderBatchActive();rlDisableDepthMask();
        for(int k=0;k<3;k++){int n=order[k];Vector3 p=origins[n];float angle=atan2f(camera.position.x-p.x,camera.position.z-p.z)*RAD2DEG;
            if(!shader){p.x+=.06f*campScale*sinf(camp.age+n);DrawModelEx(campSmoke,p,(Vector3){0,1,0},angle,v_mul((Vector3){.06f,1.15f,1},campScale),(Color){185,185,180,15});}
            else DrawModelEx(campSmoke,p,(Vector3){0,1,0},angle,v_mul((Vector3){.20f+n*.05f,(1.25f+n*.18f)*campHeightScale,1},campScale),(Color){220,223,218,255});
        }
        rlDrawRenderBatchActive();rlEnableDepthMask();
    }
    rlEnableBackfaceCulling();
}
static void camp_shutdown(void){
    camp_end();if(!campResources)return;
    Shader fallback={.id=rlGetShaderIdDefault(),.locs=rlGetShaderLocsDefault()};
    campSolid.materials[0].shader=campFlames.materials[0].shader=campSmoke.materials[0].shader=fallback;
    UnloadModel(campSolid);UnloadModel(campFlames);UnloadModel(campSmoke);
    if(campShader.id && campShader.id!=fallback.id)UnloadShader(campShader);
    campSolid=campFlames=campSmoke=(Model){0};campShader=(Shader){0};campResources=false;
}
