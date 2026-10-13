/* Batched, flat-shaded procedural geometry shared by factories and robots.
 * No physics objects or per-frame allocations are created here. */
static const Color factoryShell={178,181,168,255},factoryDark={39,48,51,255},factoryRust={143,83,47,255},factoryLight={71,198,201,255};
static Texture2D settingAtlas;
static float settingTime;
static uint64_t settingSampleFrame;
static Mesh robotMeshes[ENEMY_TYPE_COUNT];
static Color robotMeshColors[ENEMY_TYPE_COUNT];
static Mesh settingWorldMesh,settingMeshes[SETTING_DETAILS_MAX],settingEmissions[SETTING_DETAILS_MAX];
static Material settingMaterial;
static uint64_t settingMeshLayout;
static bool factoryRecording;
static Mesh factoryCapture;
static int factoryCaptureCapacity;
static void factory_vertex(Vector3 p,Vector2 uv,Color c) {
    int n=factoryCapture.vertexCount;
    if(n>=factoryCaptureCapacity) {
        int cap=factoryCaptureCapacity?factoryCaptureCapacity*2:1024;
        factoryCapture.vertices=MemRealloc(factoryCapture.vertices,cap*3*sizeof(float));
        factoryCapture.texcoords=MemRealloc(factoryCapture.texcoords,cap*2*sizeof(float));
        factoryCapture.colors=MemRealloc(factoryCapture.colors,cap*4);
        factoryCaptureCapacity=cap;
    }
    p=Vector3Transform(p,rlGetMatrixTransform());
    memcpy(factoryCapture.vertices+n*3,&p,sizeof(Vector3));memcpy(factoryCapture.texcoords+n*2,&uv,sizeof(Vector2));memcpy(factoryCapture.colors+n*4,&c,4);
    factoryCapture.vertexCount++;
}
static void factory_capture_begin(void) {factoryRecording=true;factoryCapture=(Mesh){0};factoryCaptureCapacity=0;}
static Mesh factory_capture_end(void) {
    factoryRecording=false;factoryCapture.triangleCount=factoryCapture.vertexCount/3;
    if(factoryCapture.vertexCount)UploadMesh(&factoryCapture,false);
    return factoryCapture;
}
static void setting_meshes_clear(void) {
    if(settingWorldMesh.vaoId)UnloadMesh(settingWorldMesh);settingWorldMesh=(Mesh){0};
    for(int n=0;n<SETTING_DETAILS_MAX;n++) {
        if(settingMeshes[n].vaoId)UnloadMesh(settingMeshes[n]);settingMeshes[n]=(Mesh){0};
        if(settingEmissions[n].vaoId)UnloadMesh(settingEmissions[n]);settingEmissions[n]=(Mesh){0};
    }
    settingMeshLayout=0;
}
static const char *settingLabels[]={"01 / CHASSIS ASSEMBLY","02 / ACTUATOR FITTING","03 / QUALITY INSPECTION","04 / PACKING + DISPATCH","AWAITING OPERATOR","PERSONNEL: 0","AUTONOMOUS PRODUCTION","CAUTION / LIVE MACHINERY","MOUSE ROUTE >>","HUMAN ACCESS / 2.2m","MAINTENANCE OVERDUE","CELL 07 / OFFLINE"};
static void factory_box(Vector3 p,Vector3 size,Color color) {
    static const int quads[6][4]={{4,5,7,6},{1,0,2,3},{5,1,3,7},{0,4,6,2},{6,7,3,2},{0,1,5,4}};
    static const float light[6]={.78f,.65f,.87f,.64f,1,.48f};
    rlSetTexture(rlGetTextureIdDefault());
    rlPushMatrix();rlTranslatef(p.x,p.y,p.z);rlScalef(size.x,size.y,size.z);
    if(factoryRecording) {
        static const int triangles[6]={0,1,2,0,2,3};
        for(int f=0;f<6;f++)for(int k=0;k<6;k++) {
            int c=quads[f][triangles[k]];Color shaded={(unsigned char)(color.r*light[f]),(unsigned char)(color.g*light[f]),(unsigned char)(color.b*light[f]),color.a};
            factory_vertex((Vector3){c&1?.5f:-.5f,c&2?.5f:-.5f,c&4?.5f:-.5f},(Vector2){.99f,.99f},shaded);
        }
        rlPopMatrix();return;
    }
    rlBegin(RL_QUADS);
    for(int f=0;f<6;f++) {
        rlColor4ub((unsigned char)(color.r*light[f]),(unsigned char)(color.g*light[f]),(unsigned char)(color.b*light[f]),color.a);
        for(int k=0;k<4;k++){int c=quads[f][k];rlVertex3f(c&1?.5f:-.5f,c&2?.5f:-.5f,c&4?.5f:-.5f);}
    }
    rlEnd();rlPopMatrix();
}
static void setting_init(void) {
    if(settingAtlas.id)return;
    Image atlas=GenImageColor(512,512,(Color){20,34,38,255});
    for(int i=0;i<12;i++)ImageDrawText(&atlas,settingLabels[i],8,i*32+7,18,i==8?(Color){239,177,98,255}:factoryLight);
    for(int n=0;n<10;n++){char digit[2]={(char)('0'+n),0};ImageDrawText(&atlas,digit,n*32+6,388,24,factoryLight);}
    ImageDrawRectangle(&atlas,500,500,12,12,WHITE);
    settingMaterial=LoadMaterialDefault();
    settingAtlas=LoadTextureFromImage(atlas);UnloadImage(atlas);SetTextureFilter(settingAtlas,TEXTURE_FILTER_BILINEAR);settingMaterial.maps[MATERIAL_MAP_DIFFUSE].texture=settingAtlas;
}
static void setting_sign(int row,Vector3 p,float width,float height) {
    float top=row/16.0f,bottom=(row+1)/16.0f;
    if(factoryRecording) {
        Vector3 corners[4]={{p.x-width*.5f,p.y-height*.5f,p.z},{p.x+width*.5f,p.y-height*.5f,p.z},{p.x+width*.5f,p.y+height*.5f,p.z},{p.x-width*.5f,p.y+height*.5f,p.z}};
        Vector2 uv[4]={{0,bottom},{1,bottom},{1,top},{0,top}};int tri[6]={0,1,2,0,2,3};
        for(int k=0;k<6;k++)factory_vertex(corners[tri[k]],uv[tri[k]],WHITE);return;
    }
    rlSetTexture(settingAtlas.id);rlBegin(RL_QUADS);rlColor4ub(255,255,255,255);
    rlTexCoord2f(0,bottom);rlVertex3f(p.x-width*.5f,p.y-height*.5f,p.z);
    rlTexCoord2f(1,bottom);rlVertex3f(p.x+width*.5f,p.y-height*.5f,p.z);
    rlTexCoord2f(1,top);rlVertex3f(p.x+width*.5f,p.y+height*.5f,p.z);
    rlTexCoord2f(0,top);rlVertex3f(p.x-width*.5f,p.y+height*.5f,p.z);
    rlEnd();rlSetTexture(rlGetTextureIdDefault());
}
static void setting_counter(Vector3 p) {
    if(factoryRecording)return;
    unsigned count=104287+(unsigned)(settingTime*.25f);
    for(int n=0;n<6;n++) {
        int digit=count%10;count/=10;float x=p.x+(2.5f-n)*.22f;
        rlSetTexture(settingAtlas.id);rlBegin(RL_QUADS);rlColor4ub(255,255,255,255);
        float u=digit/16.0f,v=(digit+1)/16.0f;
        rlTexCoord2f(u,416/512.0f);rlVertex3f(x-.11f,p.y-.16f,p.z);
        rlTexCoord2f(v,416/512.0f);rlVertex3f(x+.11f,p.y-.16f,p.z);
        rlTexCoord2f(v,384/512.0f);rlVertex3f(x+.11f,p.y+.16f,p.z);
        rlTexCoord2f(u,384/512.0f);rlVertex3f(x-.11f,p.y+.16f,p.z);
        rlEnd();rlSetTexture(rlGetTextureIdDefault());
    }
}
static void setting_update(void) {
    if(settingSampleFrame!=voxel_render_frame){settingSampleFrame=voxel_render_frame;settingTime=(float)GetTime();}
    if(settingCheckedGeneration==staticHashGeneration)return;
    settingCheckedGeneration=staticHashGeneration;
    for(int d=0;d<settingDetailCount;d++) {
        SettingDetail *p=&settingDetails[d];p->visible=true;
        for(int n=0;n<6;n++) {
            int *cell=p->cells[n],i=table_get_static_only(cell[0],cell[1],cell[2]);
            if(i<0 || memcmp(&voxels[i].color,&p->colors[n],sizeof(Color))) {p->visible=false;break;}
            int dx=(int)roundf(sinf(p->yaw*DEG2RAD)),dz=(int)roundf(cosf(p->yaw*DEG2RAD));
            if(table_get_static_only(cell[0]+dx,cell[1],cell[2]+dz)>=0){p->visible=false;break;}
        }
        // Check the complete backing rectangle, not just its corner samples:
        // a hole through its middle must not leave an intact floating panel.
        int lo[3],hi[3];for(int a=0;a<3;a++){lo[a]=hi[a]=p->cells[0][a];for(int n=1;n<6;n++){if(p->cells[n][a]<lo[a])lo[a]=p->cells[n][a];if(p->cells[n][a]>hi[a])hi[a]=p->cells[n][a];}}
        for(int z=lo[2];p->visible && z<=hi[2];z++)for(int y=lo[1];p->visible && y<=hi[1];y++)for(int x=lo[0];x<=hi[0];x++)if(table_get_static_only(x,y,z)<0){p->visible=false;break;}
    }
}
static bool setting_active(void) {return settingVisualsEnabled && settingMap && (netTransport.role==NET_ROLE_CLIENT || !useCustomMap) && !enemyLab;}
static void setting_detail_draw(SettingDetail *d,bool emission) {
    rlPushMatrix();rlTranslatef(d->center.x,d->center.y,d->center.z);rlRotatef(d->yaw,0,1,0);
    float w=d->width,h=d->height;
    if(!emission) {
        factory_box((Vector3){0,0,.004f},(Vector3){w,h,.014f},factoryDark);
        factory_box((Vector3){0,-h*.14f,.014f},(Vector3){w*.96f,h*.66f,.015f},factoryShell);
        if(d->style==0)for(int n=0;n<5;n++)factory_box((Vector3){-w*.3f,-h*.35f+n*h*.095f,.027f},(Vector3){w*.25f,.026f,.01f},factoryDark);
        if(d->style==1) {
            // An inset assembly schematic, including unattached limbs.
            factory_box((Vector3){-w*.22f,-h*.1f,.03f},(Vector3){w*.24f,h*.32f,.007f},factoryDark);
            factory_box((Vector3){-w*.22f,h*.12f,.03f},(Vector3){w*.12f,h*.09f,.007f},factoryLight);
            for(int side=-1;side<=1;side+=2){factory_box((Vector3){-w*.22f+side*w*.19f,-h*.06f,.03f},(Vector3){w*.06f,h*.25f,.007f},factoryDark);factory_box((Vector3){-w*.22f+side*w*.08f,-h*.33f,.03f},(Vector3){w*.06f,h*.15f,.007f},factoryDark);}
        }
        if(d->style==2)for(int n=0;n<8;n++)factory_box((Vector3){-w*.42f+n*w*.075f,-h*.13f,.03f},(Vector3){.035f,h*.57f,.013f},factoryDark);
        if(d->style==3) {
            factory_box((Vector3){-w*.21f,-h*.11f,.03f},(Vector3){w*.5f,h*.5f,.007f},factoryDark);
            for(int n=0;n<3;n++)factory_box((Vector3){-w*.28f,-h*.24f+n*h*.12f,.036f},(Vector3){w*(.22f+n*.05f),.024f,.004f},(Color){53,109,108,255});
            rlPushMatrix();rlTranslatef(-w*.23f,-h*.1f,.04f);rlRotatef(18,0,0,1);factory_box((Vector3){0},(Vector3){.012f,h*.48f,.004f},factoryShell);rlPopMatrix();
        }
        for(int n=0;n<4;n++)factory_box((Vector3){w*.2f+n*.065f,-h*.25f,.028f},(Vector3){.025f,h*.26f,.014f},factoryRust);
        factory_box((Vector3){0,-h*.46f,.023f},(Vector3){w*.95f,.06f,.012f},factoryRust);
        for(int n=0;n<5;n++)factory_box((Vector3){-w*.4f+n*w*.2f,-h*.46f,.032f},(Vector3){w*.06f,.06f,.01f},factoryDark);
        // Quiet decay and a warm mouse-scale wayfinding mark.
        factory_box((Vector3){w*.43f,-h*.33f,.03f},(Vector3){.08f,h*.2f,.006f},(Color){74,91,65,255});
        if(d->style==0)setting_sign(8,(Vector3){w*.15f,-h*.32f,.037f},w*.4f,.12f);
    }
    setting_sign(d->label,(Vector3){0,h*.28f,.037f},w*.9f,fminf(.28f,h*.22f));
    factory_box((Vector3){w*.4f,h*.04f,.032f},(Vector3){.055f,.055f,.01f},d->style==3?factoryRust:factoryLight);
    rlPopMatrix();
}
static void setting_meshes_prepare(void);
static void setting_draw(bool emission) {
    if(!setting_active())return;setting_init();setting_update();setting_meshes_prepare();
    Matrix inverse=MatrixInvert(rlGetMatrixModelview());Vector3 eye={inverse.m12,inverse.m13,inverse.m14};
    for(int d=0;d<settingDetailCount;d++) {
        SettingDetail *detail=&settingDetails[d];Vector3 delta=v_sub(eye,detail->center);
        Vector3 normal={sinf(detail->yaw*DEG2RAD),0,cosf(detail->yaw*DEG2RAD)};
        if(detail->visible && Vector3LengthSqr(delta)<55*55 && Vector3DotProduct(delta,normal)>0)
            DrawMesh(emission?settingEmissions[d]:settingMeshes[d],settingMaterial,MatrixIdentity());
    }
}
static void setting_world_geometry(void) {
    setting_init();setting_update();
    factory_box((Vector3){0,-.08f,0},(Vector3){FLOOR_SIZE*2,.15f,FLOOR_SIZE*2},(Color){63,70,69,255});
    for(int q=0;q<4;q++) {
        rlPushMatrix();rlRotatef(q*90,0,1,0);
        // Backdrop sits outside the unchanged arena collision boundary.
        float edge=FLOOR_SIZE;
        factory_box((Vector3){0,7,edge+10},(Vector3){edge*2+20,14,1},factoryDark);
        factory_box((Vector3){0,.3f,edge+.15f},(Vector3){edge*2,.6f,.25f},factoryRust);
        for(int x=-30;x<=30;x+=10) {
            factory_box((Vector3){x,5,edge+5},(Vector3){.8f,10,.8f},factoryShell);
            factory_box((Vector3){x,9,edge+4},(Vector3){.5f,.5f,6},factoryRust);
            factory_box((Vector3){x+1,9,edge+1},(Vector3){2.5f,.4f,.4f},factoryDark);
            factory_box((Vector3){x+2,8.3f,edge+1},(Vector3){.35f,1.4f,.35f},factoryShell);
            // Dormant robot shells above a conveyor; deliberately unreachable.
            factory_box((Vector3){x,1,edge+3},(Vector3){7,.6f,2},factoryDark);
            for(int n=0;n<3;n++) {
                factory_box((Vector3){x-2+n*2,1.9f,edge+3},(Vector3){.75f,1.2f,.65f},factoryShell);
                factory_box((Vector3){x-2+n*2,2.7f,edge+3},(Vector3){.45f,.35f,.4f},factoryDark);
            }
        }
        // Human door and abandoned control station establish the scale.
        factory_box((Vector3){-5,1.2f,edge+1},(Vector3){1.2f,2.4f,.2f},factoryShell);
        factory_box((Vector3){-4.65f,1.1f,edge+.88f},(Vector3){.18f,.06f,.06f},factoryRust);
        factory_box((Vector3){-2,1.2f,edge+1},(Vector3){1.5f,.15f,.6f},factoryShell);
        factory_box((Vector3){-2,1.6f,edge+1.1f},(Vector3){.8f,.55f,.1f},factoryDark);
        factory_box((Vector3){-2,.65f,edge+2},(Vector3){.65f,.15f,.65f},factoryRust);
        factory_box((Vector3){-2,1.05f,edge+2.28f},(Vector3){.65f,.7f,.12f},factoryRust);
        for(int leg=-1;leg<=1;leg+=2)factory_box((Vector3){-2+leg*.23f,.3f,edge+2},(Vector3){.08f,.6f,.5f},factoryDark);
        // A mouse-sized shelter facade cut into an abandoned equipment cabinet.
        factory_box((Vector3){5,.3f,edge+.9f},(Vector3){.7f,.6f,.08f},factoryRust);
        factory_box((Vector3){5,.24f,edge+.84f},(Vector3){.3f,.45f,.02f},factoryDark);
        rlPushMatrix();rlTranslatef(0,0,edge+.65f);rlRotatef(180,0,1,0);
        setting_sign(q%2?5:4,(Vector3){0,3,0},5,.4f);
        setting_sign(settingMap==1?q:3,(Vector3){0,4,0},7,.5f);
        setting_counter((Vector3){0,2.45f,0});
        setting_sign(8,(Vector3){-5,.8f,0},1.1f,.14f);
        setting_sign(9,(Vector3){5,2.6f,0},2,.18f);
        rlPopMatrix();
        for(int z=-30;z<30;z+=4)factory_box((Vector3){-edge+.8f,.008f,z},(Vector3){.08f,.012f,2.5f},factoryRust);
        rlPopMatrix();
    }
}
static void setting_meshes_prepare(void) {
    if(settingMeshLayout==settingLayoutVersion)return;
    setting_meshes_clear();
    factory_capture_begin();setting_world_geometry();settingWorldMesh=factory_capture_end();
    for(int n=0;n<settingDetailCount;n++) {
        factory_capture_begin();setting_detail_draw(&settingDetails[n],false);settingMeshes[n]=factory_capture_end();
        factory_capture_begin();setting_detail_draw(&settingDetails[n],true);settingEmissions[n]=factory_capture_end();
    }
    settingMeshLayout=settingLayoutVersion;
}
static void setting_world(void) {
    setting_init();setting_update();setting_meshes_prepare();
    DrawMesh(settingWorldMesh,settingMaterial,MatrixIdentity());setting_draw(false);
    for(int q=0;q<4;q++){rlPushMatrix();rlRotatef(q*90,0,1,0);rlTranslatef(0,0,FLOOR_SIZE+.65f);rlRotatef(180,0,1,0);setting_counter((Vector3){0,2.45f,0});rlPopMatrix();}
}
static void setting_shutdown(void) {
    setting_meshes_clear();
    for(int n=0;n<ENEMY_TYPE_COUNT;n++){if(robotMeshes[n].vaoId)UnloadMesh(robotMeshes[n]);robotMeshes[n]=(Mesh){0};}
    if(settingAtlas.id){settingMaterial.maps[MATERIAL_MAP_DIFFUSE].texture=(Texture2D){0};UnloadMaterial(settingMaterial);UnloadTexture(settingAtlas);}
    settingAtlas=(Texture2D){0};
}

static const char *robot_roles[]={"SECURITY / S-01","INSPECTION / I-02","MATERIAL HANDLER / H-03","SALVAGE / R-04","SAFETY / S-05","DEMOLITION / D-06","CARGO / C-07","MAINTENANCE / M-08","REPAIR POWER / P-09","MOBILE PRODUCTION / F-10"};
static void robot_geometry(int type,Vector3 pos,float yaw,float size,Color accent,float time,float speed,float action,bool coreOpen) {
    float stride=sinf(time*9)*fminf(speed*.1f,.12f);
    rlPushMatrix();rlTranslatef(pos.x,pos.y,pos.z);rlRotatef(yaw,0,1,0);rlScalef(size,size,size);
    factory_box((Vector3){0,.03f,0},(Vector3){.82f,.66f,.76f},factoryDark);
    factory_box((Vector3){0,.06f,-.36f},(Vector3){.73f,.52f,.08f},factoryShell);
    factory_box((Vector3){0,.34f,0},(Vector3){.62f,.18f,.6f},factoryShell);
    factory_box((Vector3){0,.21f,-.413f},(Vector3){.46f,.10f,.015f},factoryDark);
    factory_box((Vector3){0,.21f,-.425f},(Vector3){.22f,.04f,.01f},accent);
    for(int s=-1;s<=1;s+=2) {
        if(!factoryRecording)factory_box((Vector3){s*.26f,-.36f,stride*s},(Vector3){.21f,.25f,.45f},factoryDark);
        if(!factoryRecording)factory_box((Vector3){s*.26f,-.45f,-.08f+stride*s},(Vector3){.24f,.1f,.5f},factoryShell);
        factory_box((Vector3){s*.44f,0,-.02f},(Vector3){.12f,.45f,.24f},accent);
    }
    for(int n=0;n<4;n++)factory_box((Vector3){-.22f+n*.14f,-.08f,-.412f},(Vector3){.05f,.16f,.01f},factoryDark);
    factory_box((Vector3){.2f,.03f,-.419f},(Vector3){.12f,.025f,.015f},factoryRust);
    switch(type) {
    case ENEMY_TYPE_SWARMER:
        for(int s=-1;s<=1;s+=2)for(int n=0;n<3;n++)factory_box((Vector3){s*.49f,-.33f,(n-1)*.26f+stride*s},(Vector3){.24f,.12f,.09f},factoryRust);
        break;
    case ENEMY_TYPE_STANDARD:
        factory_box((Vector3){.27f,.02f,-.48f},(Vector3){.16f,.18f,.2f},factoryDark);break;
    case ENEMY_TYPE_GOLIATH:
        factory_box((Vector3){0,0,-.43f},(Vector3){.25f,.22f,.06f},factoryRust);break;
    case ENEMY_TYPE_SCRAPPER:case ENEMY_TYPE_LOBBER:
        for(int s=-1;s<=1;s+=2){factory_box((Vector3){s*.4f,-.2f,-.38f-action*.06f},(Vector3){.13f,.15f,.28f},factoryRust);factory_box((Vector3){s*.32f,-.15f,-.5f},(Vector3){.2f,.16f,.08f},factoryShell);}
        if(type==ENEMY_TYPE_LOBBER)factory_box((Vector3){0,.37f,.1f},(Vector3){.72f,.12f,.45f},accent);break;
    case ENEMY_TYPE_SHIELDBEARER:
        factory_box((Vector3){0,.4f,0},(Vector3){.82f,.1f,.7f},accent);break;
    case ENEMY_TYPE_SAPPER:
        for(int n=0;n<3;n++)factory_box((Vector3){(n-1)*.22f,.1f,.4f},(Vector3){.17f,.48f,.16f},factoryRust);
        factory_box((Vector3){0,.47f,0},(Vector3){.13f,.1f,.13f},sinf(time*7)>0?accent:factoryRust);break;
    case ENEMY_TYPE_BUILDER:
        factory_box((Vector3){0,.43f,0},(Vector3){.85f,.12f,.7f},factoryRust);
        factory_box((Vector3){.38f,-.17f,-.4f},(Vector3){.12f,.12f,.26f},action>0?factoryLight:factoryDark);break;
    case ENEMY_TYPE_ANCHOR:
        factory_box((Vector3){0,.47f,.15f},(Vector3){.1f,.18f,.1f},factoryShell);
        factory_box((Vector3){0,.53f,.15f},(Vector3){.45f,.06f,.06f},accent);break;
    case ENEMY_TYPE_FOUNDRY:
        factory_box((Vector3){0,.02f,-.43f},(Vector3){.35f,.35f,.08f},coreOpen?GOLD:factoryDark);
        for(int s=-1;s<=1;s+=2)factory_box((Vector3){s*.25f,.3f,.1f},(Vector3){.16f,.3f,.18f},factoryRust);break;
    }
    rlPopMatrix();
}

static void robot_draw(int type,Vector3 pos,float yaw,float size,Color accent,float time,float speed,float action,bool coreOpen) {
    if(type<0 || type>=ENEMY_TYPE_COUNT)return;
    setting_init();
    if(!robotMeshes[type].vaoId || memcmp(&robotMeshColors[type],&accent,sizeof(Color))) {
        if(robotMeshes[type].vaoId)UnloadMesh(robotMeshes[type]);
        factory_capture_begin();robot_geometry(type,(Vector3){0},0,1,accent,0,0,0,false);robotMeshes[type]=factory_capture_end();robotMeshColors[type]=accent;
    }
    Matrix transform=MatrixMultiply(MatrixMultiply(MatrixScale(size,size,size),MatrixRotateY(yaw*DEG2RAD)),MatrixTranslate(pos.x,pos.y,pos.z));
    DrawMesh(robotMeshes[type],settingMaterial,transform);
    float stride=sinf(time*9)*fminf(speed*.1f,.12f);
    rlPushMatrix();rlTranslatef(pos.x,pos.y,pos.z);rlRotatef(yaw,0,1,0);rlScalef(size,size,size);
    for(int side=-1;side<=1;side+=2) {
        factory_box((Vector3){side*.26f,-.36f,stride*side},(Vector3){.21f,.25f,.45f},factoryDark);
        factory_box((Vector3){side*.26f,-.45f,-.08f+stride*side},(Vector3){.24f,.1f,.5f},factoryShell);
    }
    if(type==ENEMY_TYPE_FOUNDRY && coreOpen)factory_box((Vector3){0,.02f,-.476f},(Vector3){.35f,.35f,.012f},GOLD);
    if(type==ENEMY_TYPE_SAPPER)factory_box((Vector3){0,.474f,0},(Vector3){.135f,.105f,.135f},sinf(time*7)>0?accent:factoryRust);
    if(type==ENEMY_TYPE_BUILDER && action>0)factory_box((Vector3){.38f,-.17f,-.54f},(Vector3){.09f,.09f,.025f},factoryLight);
    if(type==ENEMY_TYPE_LOBBER || type==ENEMY_TYPE_SCRAPPER)for(int side=-1;side<=1;side+=2)factory_box((Vector3){side*(.32f+.025f*sinf(action*8)),-.15f,-.5f},(Vector3){.16f,.17f,.09f},factoryShell);
    rlPopMatrix();
}

static void robot_emission(void) {
    if(gameMode!=GAME_MODE_FIREFIGHT)return;
    for(int i=MAX_PLAYERS;i<activePlayers;i++)if(players[i].respawn_timer<=0) {
        Player *p=&players[i];float size=p->enemyType==ENEMY_TYPE_FOUNDRY?1.5f:p->enemyType==ENEMY_TYPE_SWARMER?PLAYER_SIZE*.5f:PLAYER_SIZE;
        Color color=enemy_is_special(i)?enemy_color(p->enemyType):p->enemyType==ENEMY_TYPE_SWARMER?RED:SKYBLUE;
        rlPushMatrix();rlTranslatef(p->pos.x,p->pos.y,p->pos.z);rlRotatef(p->yaw,0,1,0);rlScalef(size,size,size);
        factory_box((Vector3){0,.21f,-.426f},(Vector3){.22f,.04f,.01f},color);rlPopMatrix();
    }
}
