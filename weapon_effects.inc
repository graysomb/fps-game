/* Stateless, bounded presentation effects. All geometry is initialized once. */
static void wx_float(Shader s,const char *name,float x){SetShaderValue(s,wx_location(s,name),&x,SHADER_UNIFORM_FLOAT);}
static void wx_int(Shader s,const char *name,int x){SetShaderValue(s,wx_location(s,name),&x,SHADER_UNIFORM_INT);}
static void wx_vec(Shader s,const char *name,Vector3 x){SetShaderValue(s,wx_location(s,name),&x,SHADER_UNIFORM_VEC3);}
static bool wx_shader(Shader s){return s.id && s.id!=rlGetShaderIdDefault();}
static void wx_draw_glock(WeaponPose p,const WeaponVisual *v,Vector3 eye,float time,WeaponRenderPass pass){
 if(pass==WEAPON_GLASS)return;
 if(wx_shader(wr_glock)){
  wx_vec(wr_glock,"eyePosition",eye);wx_int(wr_glock,"emissionOnly",pass==WEAPON_EMISSION);
  wx_int(wr_glock,"weaponKind",WEAPON_GLOCK);wx_float(wr_glock,"effectTime",time);
  wx_float(wr_glock,"shotAge",fmaxf(0,time-v->shot_time));
 }else if(pass!=WEAPON_COLOR)return;
 wr_fallback_dark=0;wr_draw_part(WR_GLOCK,p.transform);
}
static bool wx_setup(const WeaponEffect *e,int kind,Vector3 tint,Vector3 eye,float time,WeaponRenderPass pass){
 if(!wr_ready||!e||pass==WEAPON_GLASS||!isfinite(e->age)||e->age<0||!isfinite(e->progress))return false;
 if(e->lifetime<=0 || e->age>=e->lifetime)return false;
 const float *m=(const float *)&e->transform;
 for(int i=0;i<16;i++)if(!isfinite(m[i]))return false;
 if(!wx_shader(wr_effect))return pass==WEAPON_COLOR;
 wx_int(wr_effect,"effectKind",kind);wx_int(wr_effect,"emissionOnly",pass==WEAPON_EMISSION);
 wx_vec(wr_effect,"eyePosition",eye);wx_vec(wr_effect,"tint",tint);
 wx_float(wr_effect,"age",e->age);wx_float(wr_effect,"effectTime",time);
 wx_float(wr_effect,"lifetime",e->lifetime);wx_float(wr_effect,"progress",Clamp(e->progress,0,1));
 wx_float(wr_effect,"seed",(e->seed%4096)*.123f);wr_fallback_dark=0;return true;
}
static void wx_part(int part,Matrix m,Vector3 tint){
 Material mat=wr_parts[part].materials[0];
 if(!wx_shader(wr_effect))mat.maps[MATERIAL_MAP_DIFFUSE].color=(Color){(unsigned char)(tint.x*255),(unsigned char)(tint.y*255),(unsigned char)(tint.z*255),255};
 DrawMesh(wr_parts[part].meshes[0],mat,m);
}
void weapon_draw_glob(const WeaponEffect *e,Vector3 eye,float time,WeaponRenderPass pass){
 if(!wx_setup(e,0,(Vector3){.4f,1,.035f},eye,time,pass))return;
 float wobble=sinf(e->age*12+e->seed)*.035f*expf(-e->age*.7f);
 float squash=e->attached?.4f+.3f*expf(-e->age*12):1;
 Matrix m=MatrixMultiply(MatrixScale(.14f*(1+wobble),.14f*(1-wobble),.14f*squash),e->transform);
 wx_part(e->attached?WR_FX_ATTACHED:WR_FX_GLOB,m,(Vector3){.4f,1,.035f});
}
void weapon_draw_gravity(const WeaponEffect *e,Vector3 eye,float time,WeaponRenderPass pass){
 if(!wx_setup(e,1,(Vector3){.65f,.16f,1},eye,time,pass))return;
 float size=Clamp(e->age*2,0,1)*Clamp((e->lifetime-e->age)*2,0,1);
 for(int i=0;i<3;i++){
  Matrix m=MatrixMultiply(MatrixScale(size*(.5f+i*.22f),size*(.5f+i*.22f),size),MatrixRotateXYZ((Vector3){PI/2+i*.4f,e->age*(1+i*.3f),i*.8f}));
  wx_part(WR_FX_RING,MatrixMultiply(m,e->transform),(Vector3){.65f,.16f,1});
 }
 if(wx_shader(wr_effect)){
  wx_int(wr_effect,"effectKind",2);
  rlDrawRenderBatchActive();BeginBlendMode(BLEND_ADDITIVE);rlDisableDepthMask();rlDisableBackfaceCulling();
  Matrix m=MatrixMultiply(MatrixScale(size,size,size),e->transform);
  wx_part(WR_FX_RIBBON,m,(Vector3){.65f,.16f,1});
  rlDrawRenderBatchActive();rlEnableBackfaceCulling();rlEnableDepthMask();EndBlendMode();
 }
}
void weapon_draw_construction(const WeaponEffect *e,Vector3 eye,float time,WeaponRenderPass pass){
 if(!wx_setup(e,3,(Vector3){1,1,1},eye,time,pass))return;
 // Fallback cannot discard unrevealed courses: show the completed shelter.
 wx_part(e->cutaway?WR_FX_CUTAWAY:WR_FX_IGLOO,e->transform,(Vector3){1,1,1});
 if(e->progress>0 && e->progress<1){
  wx_int(wr_effect,"effectKind",1);wx_vec(wr_effect,"tint",(Vector3){.15f,.65f,1});
  float r=3*e->progress;
  Matrix m=MatrixMultiply(MatrixScale(r,r,r),MatrixRotateX(PI/2));m=MatrixMultiply(m,MatrixTranslate(0,.04f,0));
  wx_part(WR_FX_RING,MatrixMultiply(m,e->transform),(Vector3){.15f,.65f,1});
 }
}
void weapon_draw_dmr_impact(const WeaponEffect *e,Vector3 eye,float time,WeaponRenderPass pass){
 if(!wx_setup(e,4,(Vector3){.65f,.8f,1},eye,time,pass))return;
 if(pass!=WEAPON_COLOR)return;
 for(int i=0;i<12;i++){
  float a=i*2.39996f+e->seed, speed=.3f+(i%4)*.18f;
  Vector3 p={cosf(a)*speed*e->age,sinf(a)*speed*e->age-e->age*e->age*.6f,-(.4f+(i%3)*.2f)*e->age};
  float s=.07f*Clamp(1-e->age/e->lifetime,0,1);
  Matrix m=MatrixMultiply(MatrixScale(s,s,s),MatrixRotateXYZ((Vector3){a,e->age*3,a*.7f}));
  m=MatrixMultiply(m,MatrixTranslate(p.x,p.y,p.z));wx_part(WR_FX_CUBE,MatrixMultiply(m,e->transform),(Vector3){.65f,.8f,1});
 }
}
void weapon_draw_dmr_reticle(int w,int h,float opacity){
 Color c=Fade((Color){180,220,245,255},Clamp(opacity,0,1));int x=w/2,y=h/2;
 DrawLine(x-26,y,x-5,y,c);DrawLine(x+5,y,x+26,y,c);DrawLine(x,y-22,x,y-5,c);DrawLine(x,y+5,x,y+28,c);
 for(int i=1;i<=3;i++)DrawLine(x-3,y+i*9,x+3,y+i*9,c);
 DrawPixel(x,y,c);
}
