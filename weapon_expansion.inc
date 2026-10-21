/* Authored visual catalog. Mesh identity never depends on enum arithmetic. */
typedef struct WrModel { int part; const char *name; WeaponAnchors anchors; } WrModel;
static const WrModel wr_models[WEAPON_KIND_COUNT]={
 [WEAPON_STANDARD]={WR_BODY,"Standard"}, [WEAPON_DYNAMIC_LAUNCHER]={WR_LAUNCH_BARREL,"Voxel launcher"},
 [WEAPON_IMPULSE_SHOTGUN]={WR_SHOTGUN,"Shotgun"}, [WEAPON_CUTTER]={WR_CUTTER,"Cutter"},
 [WEAPON_ADHESIVE_GRENADE]={WR_ADHESIVE,"Adhesive"}, [WEAPON_RAILGUN]={WR_RAIL,"Railgun"},
 [WEAPON_ROCKET_LAUNCHER]={WR_ROCKET_GUN,"Rocket launcher"}, [WEAPON_SWORD]={WR_SWORD,"Sword"},
 [WEAPON_REINFORCED_BUILDER]={WR_BUILDER,"Builder"},
 [WEAPON_GLOB_GUN]={WR_GLOB,"Glob Gun",{{.13f,-.67f,.35f},{-.31f,-.26f,-.48f},{0,0,-1.03f}}},
 [WEAPON_BLACK_HOLE_BOMB]={WR_BLACK_HOLE,"Black Hole Bomb",{{.24f,-.28f,.1f},{-.24f,-.28f,.1f},{0,0,-.47f}}},
 [WEAPON_GLOCK]={WR_GLOCK,"Glock",{{.10f,-.57f,.34f},{-.17f,-.53f,.23f},{0,.06f,-.84f}}},
 [WEAPON_IGLOO_GRENADE]={WR_IGLOO,"Igloo Grenade",{{.24f,-.28f,.1f},{-.24f,-.28f,.1f},{0,0,-.42f}}},
 [WEAPON_DMR]={WR_DMR,"DMR",{{.13f,-.67f,.35f},{-.22f,-.25f,-.74f},{0,.03f,-1.60f}}}
};
const char *weapon_visual_name(WeaponKind k){return k>=0&&k<WEAPON_KIND_COUNT?wr_models[k].name:"Unknown";}
WeaponAnchors weapon_anchors(WeaponKind k){return k>=0&&k<WEAPON_KIND_COUNT?wr_models[k].anchors:(WeaponAnchors){0};}
static int wr_part_kind(int part){
 for(int k=0;k<WEAPON_KIND_COUNT;k++)if(wr_models[k].part==part)return k;
 return part==WR_GEL_GLASS||part==WR_FX_GLOB?WEAPON_GLOB_GUN:part==WR_SCOPE_GLASS?WEAPON_DMR:WEAPON_ROCKET_LAUNCHER;
}
static Vector3 wx_energy(int kind){
 switch(kind){case WEAPON_GLOB_GUN:return (Vector3){.38f,1,.025f};case WEAPON_BLACK_HOLE_BOMB:return (Vector3){.65f,.18f,1};
 case WEAPON_IGLOO_GRENADE:return (Vector3){.35f,.75f,1};case WEAPON_DMR:return (Vector3){.08f,.50f,1};default:return wc_energy(kind);}
}
/* Low-poly ellipsoid with flat normals; poles contain no degenerate triangles. */
static void wx_ball(WrBuilder *b,Vector3 center,Vector3 size,int material){
 for(int y=0;y<6;y++)for(int x=0;x<12;x++){
  Vector3 p[4];int ix[4]={x,x+1,x+1,x},iy[4]={y,y,y+1,y+1};
  for(int i=0;i<4;i++){float a=ix[i]*2*PI/12,t=-PI/2+iy[i]*PI/6;
   p[i]=Vector3Add(center,(Vector3){size.x*cosf(t)*cosf(a),size.y*sinf(t),size.z*cosf(t)*sinf(a)});}
  if(y>0)wr_triangle(b,p[0],p[2],p[1],material,1,0);
  if(y<5)wr_triangle(b,p[0],p[3],p[2],material,1,0);
 }
}
static void wx_ring(WrBuilder *b,float radius,float tube,int material){
 for(int i=0;i<32;i++)for(int j=0;j<4;j++){
  Vector3 p[4];int u[4]={i,i+1,i+1,i},v[4]={j,j,j+1,j+1};
  for(int k=0;k<4;k++){float a=u[k]*2*PI/32,c=v[k]*2*PI/4;
   p[k]=(Vector3){cosf(a)*(radius+tube*cosf(c)),sinf(a)*(radius+tube*cosf(c)),tube*sinf(c)};}
  wr_quad(b,p[0],p[1],p[2],p[3],material,1);
 }
}
static void wx_flash(WrBuilder *b,float z,float y){
 if(wr_fallback)return;
 b->part=11;
 for(int i=0;i<4;i++){
  float a=i*PI/2;Vector3 c={0,y,z-.26f},p={cosf(a)*.09f,y+sinf(a)*.09f,z},q={-p.x,2*y-p.y,z};
  wr_triangle(b,p,q,c,WR_ENERGY,1,0);wr_triangle(b,q,p,c,WR_ENERGY,1,0);
 }b->part=0;
}
#include "igloo_blueprint.h"
enum {WX_GRID=IGLOO_SIDE,WX_HIGH=IGLOO_HEIGHT};
static bool wx_cell(int x,int y,int z,bool cut){return igloo_cell(x,y,z,cut);}
static void wx_igloo(WrBuilder *b,bool cut){
 const int dirs[6][3]={{0,0,-1},{0,0,1},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0}};
 const int f[6][4]={{0,3,2,1},{4,5,6,7},{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2}};
 for(int y=0;y<WX_HIGH;y++)for(int z=0;z<WX_GRID;z++)for(int x=0;x<WX_GRID;x++)if(wx_cell(x,y,z,cut)){
  Vector3 p[8]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
  for(int i=0;i<8;i++)p[i]=(Vector3){(x-5.5f+p[i].x)*IGLOO_CELL_SIZE,(y+p[i].y)*IGLOO_CELL_SIZE,(z-5.5f+p[i].z)*IGLOO_CELL_SIZE};
  b->part=y;
  for(int i=0;i<6;i++)if(!wx_cell(x+dirs[i][0],y+dirs[i][1],z+dirs[i][2],cut))
   wr_quad(b,p[f[i][0]],p[f[i][1]],p[f[i][2]],p[f[i][3]],y==0?WR_CHARCOAL:WR_SILVER,.88f+.02f*((x*7+z*3+y)%6));
 }
 b->part=0;
}
static void wx_build(WrBuilder *b,int part){
 if(part==WR_FX_IGLOO||part==WR_FX_CUTAWAY){wx_igloo(b,part==WR_FX_CUTAWAY);return;}
 if(part==WR_FX_ATTACHED){
  wx_ball(b,(Vector3){0},(Vector3){1,1,1},WR_ENERGY);
  for(int i=0;i<5;i++){float a=i*2*PI/5;wx_ball(b,(Vector3){cosf(a)*.93f,sinf(a)*.93f,0},(Vector3){.39f,.43f,.5f},WR_ENERGY);}return;
 }
 if(part==WR_FX_GLOB){wx_ball(b,(Vector3){0},(Vector3){1,1,1},WR_ENERGY);return;}
 if(part==WR_FX_RING){wx_ring(b,1,.022f,WR_ENERGY);return;}
 if(part==WR_FX_CUBE){wc_detail(b,0,0,0,1,1,1,WR_SILVER);return;}
 if(part==WR_FX_RIBBON){
  for(int strand=0;strand<6;strand++)for(int i=0;i<48;i++){float a=i/48.f,c=(i+1)/48.f;wr_quad(b,(Vector3){a,-1,strand},(Vector3){a,1,strand},(Vector3){c,1,strand},(Vector3){c,-1,strand},WR_ENERGY,1);}return;
 }
 if(part==WR_GEL_GLASS){
  b->part=15;
  for(int side=-1;side<=1;side+=2)wr_quad(b,(Vector3){side*.275f,-.16f,-.39f},(Vector3){side*.275f,.20f,-.39f},(Vector3){side*.275f,.20f,.27f},(Vector3){side*.275f,-.16f,.27f},WR_GLASS,1);
  return;
 }
 if(part==WR_SCOPE_GLASS){
  b->part=16;wc_detail(b,0,.51f,-.29f,.22f,.22f,.008f,WR_GLASS);wc_detail(b,0,.51f,.49f,.22f,.22f,.008f,WR_GLASS);return;
 }
 if(part==WR_GLOB){
  wc_grip(b);wc_box(b,0,.02f,-.02f,.48f,.43f,.85f,WR_CHARCOAL);
  b->part=14;wc_box(b,0,.02f,-.06f,.535f,.34f,.64f,WR_ENERGY);b->part=0;
  for(int y=-1;y<=1;y+=2)wc_box(b,0,y*.27f,-.04f,.65f,.14f,.93f,WR_SILVER);
  wc_box(b,0,0,.48f,.65f,.59f,.25f,WR_SILVER);wc_detail(b,0,0,.61f,.20f,.22f,.02f,WR_DARK);
  wc_box(b,0,0,-.57f,.67f,.67f,.24f,WR_SILVER);
  for(int side=-1;side<=1;side+=2){
   wc_detail(b,side*.34f,0,-.58f,.012f,.48f,.15f,WR_DARK);
   for(int i=0;i<5;i++)wc_detail(b,side*.35f,-.18f+i*.09f,-.58f,.013f,.045f,.10f,WR_ENERGY);
   wc_box(b,0,side*.31f,-.92f,.58f,.16f,.48f,WR_SILVER);
   wc_detail(b,0,side*.315f,-1.02f,.59f,.165f,.04f,WR_AMBER);
  }
  b->part=13;wr_tube(b,(Vector3){0,0,-.83f},.39f,.39f,.33f,.07f,.03f,WR_DARK);
  wc_box(b,0,0,-.95f,.23f,.23f,.10f,WR_ENERGY);b->part=0;
  wc_detail(b,0,-.44f,0,.09f,.20f,.08f,WR_AMBER);wx_flash(b,-1.05f,0);return;
 }
 if(part==WR_BLACK_HOLE){
  b->part=32;wr_tube(b,(Vector3){0},.72f,.72f,.54f,.15f,.10f,WR_CHARCOAL);
  wx_ball(b,(Vector3){0,0,.13f},(Vector3){.23f,.23f,.18f},WR_CHARCOAL);
  b->part=0;wr_tube(b,(Vector3){0,0,-.28f},.49f,.49f,.04f,.025f,.015f,WR_ENERGY);
  for(int y=-1;y<=1;y+=2)wc_box(b,0,y*.42f,0,.43f,.12f,.44f,WR_SILVER);
  b->part=20;wx_ring(b,.365f,.025f,WR_ENERGY);b->part=21;wx_ring(b,.42f,.018f,WR_ENERGY);
  b->part=0;
  for(int i=0;i<4;i++){
   b->part=22+i;float x=(i&1)?1:-1,y=(i&2)?1:-1;
   wc_box(b,x*.32f,y*.32f,0,.23f,.23f,.75f,WR_SILVER);
   wc_box(b,x*.27f,y*.27f,-.35f,.32f,.32f,.17f,WR_SILVER);
   wc_box(b,x*.27f,y*.27f,.35f,.32f,.32f,.17f,WR_SILVER);
   wc_detail(b,x*.27f,y*.27f,-.44f,.12f,.055f,.02f,WR_DARK);
   wc_detail(b,x*.32f,y*.32f,-.45f,.10f,.10f,.025f,WR_ENERGY);
  }b->part=0;return;
 }
 if(part==WR_IGLOO){
  wc_box(b,0,0,0,.64f,.64f,.64f,WR_CHARCOAL);
  for(int i=0;i<4;i++){
   float x=(i&1)?1:-1,z=(i&2)?1:-1;b->part=26+i;
   wc_box(b,x*.22f,0,z*.22f,.38f,.75f,.38f,WR_SILVER);
   wc_detail(b,x*.415f,0,z*.22f,.014f,.34f,.20f,WR_ICE);
   wc_detail(b,x*.22f,0,z*.415f,.20f,.34f,.014f,WR_ICE);
   wc_detail(b,x*.22f,.25f,z*.42f,.23f,.055f,.018f,WR_AMBER);
  }b->part=0;wc_box(b,0,.37f,0,.29f,.06f,.29f,WR_ENERGY);return;
 }
 if(part==WR_GLOCK){
  // Fixed barrel stays in line while slide reciprocates independently.
  wr_tube(b,(Vector3){0,.06f,-.71f},.14f,.14f,.26f,.035f,.008f,WR_STEEL);
  wc_box(b,0,-.16f,-.13f,.29f,.15f,1.26f,WR_POLYMER);
  b->part=30;wc_box(b,0,.045f,-.13f,.31f,.28f,1.32f,WR_STEEL);
  wr_tube(b,(Vector3){0,.045f,-.805f},.31f,.28f,.02f,.065f,.012f,WR_STEEL);
  wc_detail(b,0,.20f,-.68f,.035f,.055f,.08f,WR_STEEL);
  for(int side=-1;side<=1;side+=2){
   wc_detail(b,side*.105f,.21f,.44f,.055f,.065f,.08f,WR_STEEL);
   for(int i=0;i<8;i++)wc_detail(b,side*.159f,.025f,.12f+i*.043f,.012f,.17f,.013f,WR_POLYMER);
  }
  wc_detail(b,.10f,.189f,-.04f,.09f,.008f,.22f,WR_POLYMER);
  b->part=31;wr_prism(b,(Vector3){0,-.53f,.29f},.29f,.70f,.34f,.025f,WR_POLYMER,-.21f);
  wc_box(b,0,-.88f,.37f,.31f,.065f,.37f,WR_POLYMER);
  wc_detail(b,0,-.39f,-.27f,.09f,.37f,.07f,WR_POLYMER);
  wc_detail(b,0,-.56f,-.065f,.09f,.07f,.46f,WR_POLYMER);
  wr_prism(b,(Vector3){0,-.35f,-.08f},.07f,.21f,.06f,.01f,WR_POLYMER,-.3f);
  for(int i=0;i<5;i++)wc_detail(b,0,-.33f-i*.10f,.11f+i*.021f,.285f,.025f,.035f,WR_POLYMER);
  wc_detail(b,-.165f,-.17f,.20f,.022f,.035f,.16f,WR_STEEL);
  b->part=0;if(wr_glock.id!=rlGetShaderIdDefault())wx_flash(b,-.85f,.06f);return;
 }
 if(part==WR_DMR){
  wc_grip(b);wc_box(b,0,0,-.25f,.43f,.42f,1.21f,WR_CHARCOAL);
  wc_box(b,0,.16f,-.23f,.52f,.15f,1.15f,WR_SILVER);
  wc_box(b,0,-.18f,-.23f,.47f,.09f,1.13f,WR_SILVER);
  wc_box(b,0,.03f,-1.09f,.32f,.32f,.61f,WR_CHARCOAL);
  for(int side=-1;side<=1;side+=2){
   wc_box(b,side*.19f,.08f,-1.06f,.10f,.24f,.62f,WR_SILVER);
   wc_detail(b,side*.245f,.08f,-1.06f,.012f,.045f,.40f,WR_ENERGY);
  }
  wc_box(b,0,.24f,-1.06f,.40f,.10f,.60f,WR_SILVER);
  wr_tube(b,(Vector3){0,.03f,-1.47f},.38f,.38f,.26f,.075f,.025f,WR_SILVER);
  for(int side=-1;side<=1;side+=2){
   wc_detail(b,side*.265f,.07f,-.15f,.018f,.07f,.53f,WR_ENERGY);
   for(int i=0;i<5;i++)wc_detail(b,side*.222f,-.025f,-.54f+i*.085f,.025f,.08f,.035f,WR_DARK);
   for(int i=0;i<4;i++)wc_detail(b,side*.17f,.04f,-1.25f+i*.10f,.025f,.14f,.045f,WR_SILVER);
  }
  wc_box(b,0,-.38f,-.06f,.29f,.39f,.34f,WR_CHARCOAL);
  wc_box(b,0,-.04f,.63f,.34f,.30f,.62f,WR_CHARCOAL);
  wc_box(b,0,-.12f,.94f,.39f,.65f,.15f,WR_SILVER);
  wc_detail(b,0,-.12f,1.03f,.40f,.65f,.045f,WR_DARK);
  for(int z=0;z<2;z++)wc_detail(b,0,.30f,-.15f+z*.44f,.13f,.24f,.12f,WR_DARK);
  wr_tube(b,(Vector3){0,.51f,.10f},.31f,.31f,.73f,.045f,.09f,WR_CHARCOAL);
  wr_tube(b,(Vector3){0,.51f,-.12f},.37f,.37f,.30f,.055f,.095f,WR_SILVER);
  for(int z=0;z<2;z++)wr_tube(b,(Vector3){0,.51f,-.30f+z*.80f},.36f,.36f,.16f,.055f,.07f,WR_SILVER);
  wc_box(b,.20f,.51f,.18f,.13f,.17f,.17f,WR_DARK);
  wc_detail(b,0,.70f,.18f,.16f,.06f,.16f,WR_DARK);
  wx_flash(b,-1.61f,.03f);return;
 }
}
