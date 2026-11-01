/* Cached domestic scenery and material-specific procedural geometry. */
static float town_noise(int n){unsigned v=(unsigned)n*747796405u+2891336453u;v=((v>>((v>>28)+4))^v)*277803737u;return ((v>>22)^v)%10000/10000.f;}
static void town_tube(TownMesh *b,Vector3 a,Vector3 z,float ra,float rz,Color color){
 Vector3 axis=Vector3Normalize(Vector3Subtract(z,a)),u=Vector3Normalize(Vector3CrossProduct(axis,fabsf(axis.y)>.9f?(Vector3){1,0,0}:(Vector3){0,1,0})),v=Vector3CrossProduct(axis,u);
 for(int i=0;i<8;i++){float t=i*PI/4,q=(i+1)*PI/4;Vector3 d=Vector3Add(Vector3Scale(u,cosf(t)),Vector3Scale(v,sinf(t))),e=Vector3Add(Vector3Scale(u,cosf(q)),Vector3Scale(v,sinf(q)));Vector3 p=Vector3Add(a,Vector3Scale(d,ra)),r=Vector3Add(a,Vector3Scale(e,ra)),s=Vector3Add(z,Vector3Scale(d,rz)),w=Vector3Add(z,Vector3Scale(e,rz));town_tri(b,p,r,w,color);town_tri(b,p,w,s,color);town_tri(b,z,s,w,color);town_tri(b,a,r,p,color);}
}
static void town_ring(TownMesh *b,Vector3 p,float r,float width,Color color){
 for(int k=0;k<24;k++){float a=k*PI/12,z=(k+1)*PI/12;Vector3 v={p.x+r*cosf(a),p.y+r*sinf(a),p.z},w={p.x+r*cosf(z),p.y+r*sinf(z),p.z};town_tube(b,v,w,width,width,color);}
}
static void town_leaf(TownMesh *b,Vector3 p,float a,float length,Color c){
 Vector3 tip={p.x+cosf(a)*length,p.y+.11f,p.z+sinf(a)*length},left={p.x-sinf(a)*length*.3f,p.y,p.z+cosf(a)*length*.3f},right={p.x+sinf(a)*length*.3f,p.y,p.z-cosf(a)*length*.3f},ridge={p.x+cosf(a)*length*.35f,p.y+.075f,p.z+sinf(a)*length*.35f};
 town_tri(b,left,ridge,tip,c);town_tri(b,ridge,right,tip,c);town_tri(b,tip,ridge,left,c);town_tri(b,tip,right,ridge,c);
}
static void town_plant(TownMesh *b,Vector3 p,int seed,float scale,bool flower){
 for(int k=0;k<7;k++){float a=k*2.399f+seed,dy=k*.045f*scale;town_leaf(b,Vector3Add(p,(Vector3){0,dy,0}),a,(.15f+town_noise(seed+k)*.2f)*scale,(Color){(unsigned char)(45+k*5),(unsigned char)(95+k*8),48,255});}
 if(flower)for(int k=0;k<5;k++)town_leaf(b,Vector3Add(p,(Vector3){0,.4f*scale,0}),k*2*PI/5,.095f*scale,seed%2?(Color){233,177,90,255}:(Color){177,132,180,255});
}
static void town_crate(TownMesh *b,Vector3 p,float scale){
 town_box(&b[T_WOOD],p,(Vector3){scale,scale*.8f,scale*.8f},(Color){93,65,40,255},NULL);
 for(int k=-1;k<=1;k+=2){town_box(&b[T_WOOD],Vector3Add(p,(Vector3){k*scale*.42f,0,-scale*.42f}),(Vector3){scale*.1f,scale*.85f,.06f},(Color){137,98,56,255},NULL);town_box(&b[T_WOOD],Vector3Add(p,(Vector3){0,k*scale*.3f,-scale*.43f}),(Vector3){scale,.06f,.07f},(Color){137,98,56,255},NULL);}
}
static void town_lantern(TownMesh *b,Vector3 p){
 town_box(&b[T_GLOW],p,(Vector3){.17f,.25f,.17f},townGold,NULL);
 for(int s=-1;s<=1;s+=2){town_box(&b[T_METAL],Vector3Add(p,(Vector3){0,s*.15f,0}),(Vector3){.25f,.055f,.25f},(Color){49,55,58,255},NULL);for(int q=-1;q<=1;q+=2)town_tube(&b[T_METAL],Vector3Add(p,(Vector3){s*.1f,-.14f,q*.1f}),Vector3Add(p,(Vector3){s*.1f,.14f,q*.1f}),.016f,.016f,(Color){58,61,60,255});}
 town_ring(&b[T_METAL],Vector3Add(p,(Vector3){0,.23f,0}),.08f,.014f,(Color){58,61,60,255});
}
static void town_home(TownMesh *b,float x,float z,int seed){
 Color timber={117,79,43,255},stone={151,143,119,255};
 /* Deep circular door, planked inset and chunky masonry surround. */
 town_disc(&b[T_WOOD],(Vector3){x,1.25f,z-.11f},1.16f,.08f,timber);
 town_ring(&b[T_STONE],(Vector3){x,1.25f,z-.17f},1.27f,.17f,stone);
 town_ring(&b[T_WOOD],(Vector3){x,1.25f,z-.21f},1.12f,.075f,(Color){73,46,29,255});
 for(int k=-4;k<=4;k++){float xx=k*.23f,h=sqrtf(fmaxf(0,1.1f*1.1f-xx*xx));town_box(&b[T_WOOD],(Vector3){x+xx,1.25f,z-.145f},(Vector3){.013f,h*2,.018f},(Color){54,35,24,255},NULL);}
 for(int k=-1;k<=1;k+=2){town_box(&b[T_METAL],(Vector3){x-.70f,1.25f+k*.42f,z-.23f},(Vector3){.42f,.085f,.05f},(Color){51,57,59,255},NULL);}
 town_ring(&b[T_METAL],(Vector3){x+.63f,1.23f,z-.27f},.10f,.025f,(Color){181,141,68,255});
 for(int side=-1;side<=1;side+=2){float wx=x+side*3.4f;
 town_disc(&b[T_GLOW],(Vector3){wx,2.15f,z-.09f},.66f,.05f,(Color){234,166,84,255});
 town_ring(&b[T_STONE],(Vector3){wx,2.15f,z-.2f},.77f,.13f,stone);
 town_ring(&b[T_WOOD],(Vector3){wx,2.15f,z-.22f},.65f,.055f,timber);
 town_box(&b[T_WOOD],(Vector3){wx,2.15f,z-.24f},(Vector3){1.3f,.07f,.1f},timber,NULL);town_box(&b[T_WOOD],(Vector3){wx,2.15f,z-.25f},(Vector3){.07f,1.3f,.1f},timber,NULL);
 town_box(&b[T_WOOD],(Vector3){wx,1.32f,z-.36f},(Vector3){1.8f,.14f,.65f},timber,NULL);
 for(int k=-1;k<=1;k+=2)town_box(&b[T_CLOTH],(Vector3){wx+k*.45f,2.16f,z-.19f},(Vector3){.23f,1.05f,.035f},seed%2?(Color){143,104,72,255}:(Color){97,137,129,255},NULL);
 town_box(&b[T_GLASS],(Vector3){wx,2.15f,z-.26f},(Vector3){1.1f,1.1f,.01f},(Color){139,183,194,38},NULL);
 town_box(&b[T_WOOD],(Vector3){wx,1.48f,z-.55f},(Vector3){1.4f,.22f,.36f},timber,NULL);
 for(int k=0;k<5;k++)town_plant(&b[T_LEAF],(Vector3){wx-.55f+k*.27f,1.59f,z-.55f},seed+k,.7f,true);
 }
 /* Plank porch and small awning on braced salvaged supports. */
 for(int k=0;k<12;k++)town_box(&b[T_WOOD],(Vector3){x-1.8f+k*.33f,.1f,z-1.1f},(Vector3){.31f,.16f,1.8f},timber,NULL);
 for(int s=-1;s<=1;s+=2){Vector3 a={x+s*2,0,z-1.8f},v={a.x,3.05f,a.z};town_tube(&b[T_WOOD],a,v,.075f,.06f,timber);town_tube(&b[T_WOOD],(Vector3){a.x,2.3f,a.z},(Vector3){a.x-s*.5f,3.05f,a.z},.04f,.04f,timber);town_lantern(b,(Vector3){a.x,2.5f,a.z+.18f});}
 for(int k=0;k<12;k++){float xx=x-2+k*.35f;Color c=k%3?(Color){145,102,65,255}:(Color){188,156,107,255};town_tri(&b[T_CLOTH],(Vector3){xx,3.5f,z-.15f},(Vector3){xx+.35f,3.5f,z-.15f},(Vector3){xx+.35f,3.08f,z-2},c);town_tri(&b[T_CLOTH],(Vector3){xx,3.5f,z-.15f},(Vector3){xx+.35f,3.08f,z-2},(Vector3){xx,3.08f,z-2},c);}
 for(int side=-1;side<=1;side+=2){for(int k=0;k<3;k++)town_tube(&b[T_WOOD],(Vector3){x+side*1.9f,.17f,z-.25f-k*.45f},(Vector3){x+side*1.9f,.7f,z-.25f-k*.45f},.035f,.035f,timber);town_tube(&b[T_WOOD],(Vector3){x+side*1.9f,.72f,z-.25f},(Vector3){x+side*1.9f,.72f,z-1.2f},.045f,.045f,timber);}
 town_crate(b,(Vector3){x-2.6f,.4f,z-1},.65f);town_crate(b,(Vector3){x-2.6f,1,z-.9f},.5f);
 for(int k=0;k<7;k++)town_tube(&b[T_WOOD],(Vector3){x+2.7f+(k%3)*.19f,.12f+(k/3)*.18f,z-1.3f},(Vector3){x+2.7f+(k%3)*.19f,.12f+(k/3)*.18f,z-.75f},.10f,.08f,(Color){92,60,37,255});
 for(int k=0;k<12;k++){float xx=x-5+k*.9f;town_tube(&b[T_METAL],(Vector3){xx,4.15f+.2f*cosf(k),z-.5f},(Vector3){xx+.9f,4.15f+.2f*cosf(k+1),z-.5f},.012f,.012f,(Color){43,50,49,255});if(k%2==0)town_lantern(b,(Vector3){xx,3.92f+.2f*cosf(k),z-.5f});}
 for(int vine=0;vine<9;vine++){float xx=x-5.4f+vine*1.35f,len=.7f+town_noise(seed+vine)*2.1f;for(int k=0;k<12;k++){float y=5.5f-k*len/12;town_tube(&b[T_LEAF],(Vector3){xx+.08f*sinf(k),y,z-.3f},(Vector3){xx+.08f*sinf(k+1),y-len/12,z-.3f},.02f,.018f,(Color){65,83,42,255});town_plant(&b[T_LEAF],(Vector3){xx+.10f*sinf(k),y,z-.35f},vine*13+k,.7f,false);
 for(int side=-1;side<=1;side+=2){Vector3 stem={xx,y,z-.37f},tip={xx+side*.36f,y-.25f,z-.40f},left={xx+side*.10f,y-.20f,z-.42f},right={xx+side*.29f,y+.025f,z-.40f};town_tri(&b[T_LEAF],stem,left,tip,(Color){61,106,49,255});town_tri(&b[T_LEAF],stem,tip,right,(Color){85,128,57,255});town_tri(&b[T_LEAF],tip,left,stem,(Color){61,106,49,255});town_tri(&b[T_LEAF],right,tip,stem,(Color){85,128,57,255});}}}
}
static void town_details(TownMesh *b){
 /* Each material has geometry even when a distant feature is culled. */
 town_box(&b[T_SOIL],(Vector3){15,.39f,8.3f},(Vector3){4,.02f,3.6f},(Color){48,36,26,255},NULL);
 for(int row=0;row<3;row++)for(int k=0;k<9;k++){Vector3 p={13.5f+k*.36f,.65f,7+row*1.3f};town_plant(&b[T_LEAF],p,k+row*11,.8f,row==2);if(row==1)town_bush(&b[T_LEAF],Vector3Add(p,(Vector3){0,.1f,0}),(Vector3){.16f,.18f,.16f},(Color){205,88,42,255});}
 for(int k=0;k<42;k++){float a=k*2.399f;Vector3 p={11.8f+cosf(a)*.35f,.04f,-13+k*.65f};town_plant(&b[T_LEAF],p,k,.75f,k%3==0);}
 for(int k=0;k<12;k++)town_plant(&b[T_LEAF],(Vector3){20.5f,.12f,7+k*.33f},k,1,k%4==0);
 town_tube(&b[T_METAL],(Vector3){20.5f,.5f,11.5f},(Vector3){20.5f,.5f,9.3f},.09f,.09f,(Color){91,124,126,255});town_tube(&b[T_METAL],(Vector3){20.5f,.5f,9.3f},(Vector3){20.2f,.5f,9.3f},.09f,.09f,(Color){91,124,126,255});town_tube(&b[T_WATER],(Vector3){20.2f,.46f,9.3f},(Vector3){20.2f,.05f,9.3f},.035f,.05f,(Color){103,167,177,255});
 /* Workshop backboard, hung tools, salvaged machine and loose parts. */
 town_box(&b[T_WOOD],(Vector3){12,1.55f,12.48f},(Vector3){2.5f,1.2f,.09f},townWood,NULL);
 for(int k=0;k<7;k++){float x=11.05f+k*.3f;town_tube(&b[T_METAL],(Vector3){x,1.1f,12.37f},(Vector3){x,1.8f,12.37f},.025f,.025f,(Color){137,150,153,255});town_box(&b[T_METAL],(Vector3){x,1.75f,12.37f},(Vector3){.14f,.07f,.09f},(Color){145,156,155,255},NULL);}
 town_box(&b[T_METAL],(Vector3){12.5f,1.03f,12},(Vector3){.55f,.4f,.46f},(Color){115,135,141,255},NULL);town_ring(&b[T_METAL],(Vector3){12.5f,1.06f,11.74f},.16f,.035f,(Color){53,61,65,255});town_box(&b[T_GLOW],(Vector3){12.5f,1.06f,11.72f},(Vector3){.08f,.08f,.015f},(Color){83,192,211,255},NULL);
 for(int k=0;k<4;k++)town_crate(b,(Vector3){10.7f+k*.5f,.24f,12.9f},.45f);
 town_lantern(b,(Vector3){11.2f,1.02f,12});
 /* Lantern canopy over the feast, made from reclaimed posts and sagging wires. */
 for(int side=-1;side<=1;side+=2){town_tube(&b[T_WOOD],(Vector3){side*3.5f,1.8f,2.8f},(Vector3){side*3.5f,4.6f,2.8f},.065f,.05f,townWood);town_crate(b,(Vector3){side*3.5f,2.02f,2.8f},.45f);}
 for(int k=0;k<12;k++){float x=-3.5f+k*7.f/12,y=4.6f-.65f*sinf(k*PI/12);town_tube(&b[T_METAL],(Vector3){x,y,2.8f},(Vector3){x+7.f/12,4.6f-.65f*sinf((k+1)*PI/12),2.8f},.012f,.012f,(Color){48,55,52,255});if(k%2==0)town_lantern(b,(Vector3){x,y-.2f,2.8f});}
 town_box(&b[T_CLOTH],(Vector3){0,1.78f,0},(Vector3){4.7f,.015f,4.7f},(Color){69,100,102,255},NULL);
 for(int k=-2;k<=2;k++)town_box(&b[T_CLOTH],(Vector3){k*.85f,1.792f,0},(Vector3){.04f,.008f,4.65f},(Color){176,149,103,255},NULL);
 /* Feast: planked seats, bowls on side tables, rugs, cups and provisions. */
 for(int k=0;k<4;k++){float a=k*PI/2;Vector3 p={sinf(a)*2.1f,2.25f,cosf(a)*2.1f};town_box(&b[T_CLOTH],p,(Vector3){.7f,.055f,.7f},k%2?(Color){129,82,65,255}:(Color){72,124,124,255},NULL);}
 town_crate(b,(Vector3){-2.4f,2.05f,-1.6f},.6f);town_lantern(b,(Vector3){-2.4f,2.5f,-1.6f});
 for(int k=0;k<6;k++)town_plant(&b[T_LEAF],(Vector3){-3.7f+k*.4f,1.8f,2.9f},k,.8f,true);
 town_box(&b[T_WOOD],(Vector3){-1.6f,2.1f,1.3f},(Vector3){.9f,.12f,.6f},townWood,NULL);
 for(int k=0;k<3;k++)town_tube(&b[T_METAL],(Vector3){-1.9f+k*.3f,2.18f,1.3f},(Vector3){-1.9f+k*.3f,2.33f,1.3f},.075f,.08f,(Color){156,139,104,255});
}
static void town_bake_shadow(void){
 town.depthShader=LoadShader("shaders/victory_shadow.vert","shaders/victory_shadow.frag");town.shadow=LoadRenderTexture(2048,2048);SetTextureFilter(town.shadow.texture,TEXTURE_FILTER_POINT);
 Camera3D sun={{-45,80,-35},{0,0,0},{0,1,0},100,CAMERA_ORTHOGRAPHIC};
 BeginTextureMode(town.shadow);ClearBackground(WHITE);BeginMode3D(sun);town.lightMatrix=MatrixMultiply(rlGetMatrixModelview(),rlGetMatrixProjection());
 for(int i=0;i<T_COUNT;i++)if(i!=T_GLASS && i!=T_WATER){Shader old=town.models[i].materials[0].shader;town.models[i].materials[0].shader=town.depthShader;DrawModel(town.models[i],(Vector3){0},1,WHITE);town.models[i].materials[0].shader=old;}
 EndMode3D();EndTextureMode();
}
