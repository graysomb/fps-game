/* Enemy-earned equipment: a copy is minted only when a carrier is assigned. */
static const FfEquipmentRule fpRules[]={
 {PICKUP_GLOB,2,ENEMY_TYPE_STANDARD,'E'},{PICKUP_REINFORCED,2,ENEMY_TYPE_BUILDER,'E'},
 {PICKUP_GLOCK,3,ENEMY_TYPE_STANDARD,'C'},{PICKUP_GRENADE,3,ENEMY_TYPE_SAPPER,'C'},
 {PICKUP_ADHESIVE,4,ENEMY_TYPE_SAPPER,'C'},
 {PICKUP_SHOTGUN,5,ENEMY_TYPE_STANDARD,'B'},{PICKUP_IGLOO,5,ENEMY_TYPE_SAPPER,'B'},
 {PICKUP_CUTTER,6,ENEMY_TYPE_STANDARD,'B'},{PICKUP_DMR,6,ENEMY_TYPE_STANDARD,'B'},
 {PICKUP_DYNAMIC_SHOT,7,ENEMY_TYPE_STANDARD,'A'},{PICKUP_SWORD,7,ENEMY_TYPE_STANDARD,'A'},
 {PICKUP_BLACK_HOLE,8,ENEMY_TYPE_SAPPER,'A'},{PICKUP_GOLD_TETHER,9,ENEMY_TYPE_LOBBER,'S'},
 {PICKUP_RAILGUN,10,ENEMY_TYPE_STANDARD,'S'},{PICKUP_ROCKET,11,ENEMY_TYPE_STANDARD,'S'}
};
#define FP_RULE_COUNT ((int)(sizeof(fpRules)/sizeof(fpRules[0])))
static bool fp_active(void){return gameMode==GAME_MODE_FIREFIGHT&&!firefightPractice&&!enemyLab&&!equipmentLab;}
static int fp_team_size(void){if(netTransport.role==NET_ROLE_OFFLINE)return firefightHumans;int count=0;for(int i=0;i<MAX_PLAYERS;i++)count+=netPlayerPresent[i];return count>0?count:1;}
static bool fp_carrier(int i){return fp_active()&&i>=MAX_PLAYERS&&i<MAX_COMBATANTS&&fpCarrierType[i]>0;}
static int fp_rule(int type){for(int n=0;n<FP_RULE_COUNT;n++)if(fpRules[n].type==type)return n;return -1;}
static bool fp_reserve(size_t count){
 if(count<=fpCapacity)return true;size_t cap=fpCapacity?fpCapacity*2:128;while(cap<count)cap*=2;
 if(cap>SIZE_MAX/sizeof(FfItem))return false;FfItem *items=realloc(fpItems,cap*sizeof(*items));if(!items)return false;
 fpItems=items;fpCapacity=cap;return true;
}
static void fp_new_run(void){
 memset(fpCarrierType,0,sizeof(fpCarrierType));memset(fpCarrierClock,0,sizeof(fpCarrierClock));memset(fpCarrierPanel,0,sizeof(fpCarrierPanel));
 fpPendingCount=0;fpRestartTimer=0;fpSupplyTimer=0;fpBannerRead=fpBannerWrite=0;fpBannerTimer=0;
 if(fpRestoring)return;
 fpCount=0;fpNextId=1;fpUnlocked=0;fpRepeatCursor=0;++fpRevision;++fpGeneration;
 free(fpCheckpoint.items);fpCheckpoint=(FfCheckpoint){0};
}
static void fp_shutdown(void){free(fpIncoming);fpIncoming=NULL;free(fpItems);fpItems=NULL;fpCount=fpCapacity=0;free(fpCheckpoint.items);fpCheckpoint=(FfCheckpoint){0};}
static void fp_ground(FfItem *item){
 Vector3 origin=world_wrap(item->pos),safe;
 for(int n=0;n<49;n++){
  float radius=n?((n-1)/8+1)*.75f:0,angle=(n-1)*PI/4;
  Vector3 at=world_wrap(v_add(origin,(Vector3){cosf(angle)*radius,0,sinf(angle)*radius}));
  if(clarity_supply_position(at,0,&safe)){item->pos=safe;return;}
 }
 item->pos=origin;item->pos.y=fmaxf(.8f,item->pos.y);
}
static void fp_drop(FfItem *item,Vector3 at){item->owner=-1;item->slot=-1;item->pos=at;fp_ground(item);++fpRevision;}
static void fp_death(int i,Vector3 at){
 if(!fp_active()||netTransport.role==NET_ROLE_CLIENT)return;
 for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner==i)fp_drop(&fpItems[n],at);
 if(i>=MAX_PLAYERS&&i<MAX_COMBATANTS)fpCarrierType[i]=0;
}
static void fp_queue_unlock(int type){
 int rule=fp_rule(type);if(rule<0 || (fpUnlocked&(1u<<rule)))return;fpUnlocked|=1u<<rule;++fpRevision;
 if(fpBannerWrite-fpBannerRead<32)fpBannerQueue[fpBannerWrite++%32]=type;
}
static bool fp_collect(int i,size_t n,bool replace){
 if(!fp_active()||netTransport.role==NET_ROLE_CLIENT||i<0||i>=firefightHumans||players[i].respawn_timer>0||n>=fpCount)return false;
 FfItem *item=&fpItems[n];if(item->owner!=-1||combatantLife[i]!=LIFE_ALIVE)return false;int type=item->type,gun=eq_pickup_gun(type);Equipment *e=&equipment[i];
 if(gun>=0){
  if(e->gun[0]==gun||e->gun[1]==gun)return false;
  int slot=e->gun[0]==EQ_EMPTY?0:e->gun[1]==EQ_EMPTY?1:replace?e->selected:-1;if(slot<0)return false;
  int old=e->gun[slot];if(!eq_pickup(i,type,replace))return false;
  for(size_t k=0;k<fpCount;k++)if(fpItems[k].owner==i&&fpItems[k].slot==slot&&eq_pickup_gun(fpItems[k].type)==old)fp_drop(&fpItems[k],players[i].pos);
  item->owner=i;item->slot=slot;
 }else if(type==PICKUP_GOLD_TETHER){if(players[i].goldTetherCharged)return false;players[i].goldTetherCharged=true;item->owner=i;}
 else if(type==PICKUP_ADHESIVE || type==PICKUP_REINFORCED){if(!eq_pickup(i,type,false))return false;item->owner=i;}
 else {
  int kind=type==PICKUP_BLACK_HOLE?GRENADE_GRAVITY:type==PICKUP_IGLOO?GRENADE_IGLOO:GRENADE_REGULAR;
  if(kind){if(!ex_grenade_pickup(i,kind))return false;}
  else {if(grenadeInventory[i]>=2)return false;if(!ex_grenade_count(i,exSelectedGrenade[i]))exSelectedGrenade[i]=GRENADE_REGULAR;grenadeInventory[i]++;}
  item->owner=-2;
 }
 fp_queue_unlock(type);++fpRevision;play_sfx_at(SFX_SHIELD,players[i].pos,i,1,1);return true;
}
static int fp_near(int i,bool automatic){
 if(!fp_active()||i<0||i>=firefightHumans||players[i].respawn_timer>0)return -1;
 int best=-1;float distance=1.4f;
 for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner==-1){
  int gun=eq_pickup_gun(fpItems[n].type);Equipment *e=&equipment[i];
  if(gun>=0&&(e->gun[0]==gun||e->gun[1]==gun))continue;
  if(automatic&&gun>=0&&e->gun[0]!=EQ_EMPTY&&e->gun[1]!=EQ_EMPTY)continue;
  if(fpItems[n].type==PICKUP_GOLD_TETHER&&players[i].goldTetherCharged)continue;
  if((fpItems[n].type==PICKUP_ADHESIVE&&(e->upgrades&EQ_ADHESIVE))||(fpItems[n].type==PICKUP_REINFORCED&&(e->upgrades&EQ_BUILD)))continue;
  int grenade=fpItems[n].type==PICKUP_GRENADE?GRENADE_REGULAR:fpItems[n].type==PICKUP_BLACK_HOLE?GRENADE_GRAVITY:fpItems[n].type==PICKUP_IGLOO?GRENADE_IGLOO:-1;
  if(grenade>=0&&ex_grenade_count(i,grenade)>=2)continue;
  Vector3 d=world_delta(fpItems[n].pos,players[i].pos);float length=v_length(d);VoxelHit hit;
  if(length<distance && (length<.15f || !first_voxel_hit_detailed((Ray){players[i].pos,v_norm(d)},fmaxf(0,length-.15f),-1,false,&hit))){best=(int)n;distance=length;}
 }return best;
}
static bool fp_interact(int i){int n=fp_near(i,false);return n>=0 && fp_collect(i,(size_t)n,true);}
static void fp_prepare_wave(int wave){
 fpPendingCount=0;memset(fpCarrierType,0,sizeof(fpCarrierType));if(!fp_active())return;
 size_t kept=0;for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner!=-2)fpItems[kept++]=fpItems[n];
 if(kept!=fpCount){fpCount=kept;++fpRevision;}
 int wanted[8],count=0;
 for(int n=0;n<FP_RULE_COUNT;n++)if(fpRules[n].wave==wave)wanted[count++]=n;
 if(wave>=3){int repeats=(fp_team_size()+1)/2;for(int k=0;k<repeats;k++)for(int n=0;n<FP_RULE_COUNT;n++){
   int rule=fpRepeatCursor++%FP_RULE_COUNT;if(fpRules[rule].wave<wave){wanted[count++]=rule;break;}
 }}
 for(int n=0;n<count;n++){
  int replace=ffRemaining[ENEMY_TYPE_STANDARD]>0?ENEMY_TYPE_STANDARD:ffRemaining[ENEMY_TYPE_SWARMER]>0?ENEMY_TYPE_SWARMER:-1;
  // Reserve replacements before adding carriers, so later assignments cannot
  // accidentally consume an earlier carrier's quota.
  if(replace<0)break;ffRemaining[replace]--;fpPending[fpPendingCount++]=wanted[n];
 }
 for(int n=0;n<fpPendingCount;n++)ffRemaining[fpRules[fpPending[n]].chassis]++;
}
static void fp_spawn(int i,int chassis){
 eq_reset(i);ex_reset_player(i);fpCarrierType[i]=0;
 if(!fp_active())return;
 for(int n=0;n<fpPendingCount;n++)if(fpPending[n]>=0 && fpRules[fpPending[n]].chassis==chassis){
  int type=fpRules[fpPending[n]].type;fpPending[n]=-1;
  if(!fp_reserve(fpCount+1)){TraceLog(LOG_ERROR,"Could not allocate carrier item");return;}
  fpItems[fpCount++]=(FfItem){.id=fpNextId++,.type=type,.owner=i,.slot=0,.pos=players[i].pos};++fpRevision;
  fpCarrierType[i]=type;fpCarrierClock[i]=0;fpCarrierPanel[i]=0;
  int gun=eq_pickup_gun(type);if(gun>=0)equipment[i].gun[0]=gun;
  equipment[i].upgrades=type==PICKUP_ADHESIVE?EQ_ADHESIVE:type==PICKUP_REINFORCED?EQ_BUILD:0;
  players[i].goldTetherCharged=type==PICKUP_GOLD_TETHER;players[i].dynamicShotActive=gun==EQ_VOXEL;
  grenadeInventory[i]=type==PICKUP_GRENADE||type==PICKUP_ADHESIVE?2:0;
  if(type==PICKUP_BLACK_HOLE||type==PICKUP_IGLOO){int kind=type==PICKUP_BLACK_HOLE?GRENADE_GRAVITY:GRENADE_IGLOO;exGrenades[i][kind-1]=2;exSelectedGrenade[i]=kind;}
  return;
 }
}
static int fp_carrier_visual(int i){
 int type=fpCarrierType[i];if(type==PICKUP_REINFORCED)return WEAPON_REINFORCED_BUILDER;
 if(type==PICKUP_GRENADE||type==PICKUP_ADHESIVE)return WEAPON_ADHESIVE_GRENADE;
 if(type==PICKUP_BLACK_HOLE)return WEAPON_BLACK_HOLE_BOMB;if(type==PICKUP_IGLOO)return WEAPON_IGLOO_GRENADE;
 return eqGuns[eq_gun(i)].visual;
}
static void fp_tick(float dt){
 if(!fp_active())return;
 if(fpBannerTimer>0)fpBannerTimer=fmaxf(0,fpBannerTimer-dt);
 if(fpBannerTimer<=0 && fpBannerRead<fpBannerWrite){fpBannerTimer=3.5f;fpBannerRead++;play_sfx(SFX_WIN);}
 if(netTransport.role==NET_ROLE_CLIENT)return;
 if(fpRestartTimer>0){fpRestartTimer=fmaxf(0,fpRestartTimer-dt);if(fpRestartTimer<=0)fp_restore();return;}
 for(int i=0;i<firefightHumans;i++)if(combatantLife[i]==LIFE_ALIVE){int n=fp_near(i,true);if(n>=0)fp_collect(i,(size_t)n,false);}
 fpSupplyTimer-=dt;if(fpSupplyTimer<=0){
  fpSupplyTimer=.25f;static size_t cursor=0;
  for(size_t n=0;n<16 && n<fpCount;n++){if(cursor>=fpCount)cursor=0;FfItem *item=&fpItems[cursor++];if(item->owner==-1){Vector3 old=item->pos;fp_ground(item);if(v_length(world_delta(old,item->pos))>.01f)++fpRevision;}}
 }
}
static void fp_draw(WeaponRenderPass pass){
 // Identical copies in a pile share a display model; ownership is never merged.
 // Keep late endless runs bounded even when many unclaimed weapons accumulate.
 struct DisplayCell {int x,y,z,type;bool used;} shown[256]={0};int detailed=0;
 if(!fp_active())return;Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};float now=weaponFrameSampled?weaponFrameTime:(float)GetTime();
 for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner==-1 && periodic_visible(fpItems[n].pos,2) && v_length(world_delta(fpItems[n].pos,eye))<60){
  FfItem *item=&fpItems[n];Vector3 p=item->pos;int x=(int)floorf(p.x*2),y=(int)floorf(p.y*2),z=(int)floorf(p.z*2);
  unsigned hash=((uint32_t)x*73856093u^(uint32_t)y*19349663u^(uint32_t)z*83492791u^(uint32_t)item->type*977u)&255;bool duplicate=false;
  for(int probe=0;probe<256;probe++){struct DisplayCell *cell=&shown[(hash+probe)&255];if(!cell->used){*cell=(struct DisplayCell){x,y,z,item->type,true};break;}if(cell->x==x&&cell->y==y&&cell->z==z&&cell->type==item->type){duplicate=true;break;}}
  if(duplicate)continue;p.y+=sinf(now*2+item->id)*.12f;
  if(v_length(world_delta(p,eye))>18 || detailed>=48){if(pass!=WEAPON_GLASS)DrawCube(p,.18f,.18f,.18f,GOLD);continue;}
  detailed++;int kind=eq_pickup_kind(item->type);
  if(item->type==PICKUP_GOLD_TETHER)weapon_draw_pickup(false,p,eye,now,pass);
  else {if(kind<0)kind=WEAPON_ADHESIVE_GRENADE;WeaponVisual visual={.shot_time=-1000,.charge=.7f,.activity=1};Matrix m=MatrixMultiply(MatrixScale(.4f,.4f,.4f),MatrixRotateY(now*.6f));m=MatrixMultiply(m,MatrixTranslate(p.x,p.y,p.z));weapon_draw((WeaponPose){.kind=kind,.transform=m},&visual,eye,now,false,pass);}
 }
}
static void fp_icon(int type,int x,int y,Color c){
 int gun=eq_pickup_gun(type);if(gun>=0){eq_icon(gun,x,y,c);return;}
 if(type==PICKUP_REINFORCED){for(int a=0;a<3;a++)for(int b=0;b<3;b++)DrawRectangleLines(x+a*8,y+b*8,7,7,c);return;}
 if(type==PICKUP_GOLD_TETHER){DrawLineEx((Vector2){x+12,y+24},(Vector2){x+12,y+12},3,c);DrawLineEx((Vector2){x+12,y+12},(Vector2){x,y+5},3,c);DrawLineEx((Vector2){x,y+5},(Vector2){x+5,y},3,c);DrawLineEx((Vector2){x+12,y+12},(Vector2){x+24,y+5},3,c);DrawLineEx((Vector2){x+24,y+5},(Vector2){x+19,y},3,c);return;}
 if(type==PICKUP_IGLOO){DrawCircleSector((Vector2){x+12,y+23},12,180,360,8,c);DrawRectangle(x+9,y+17,6,7,BLACK);DrawRectangle(x+15,y+5,5,9,c);return;}
 if(type==PICKUP_BLACK_HOLE){DrawCircleLines(x+12,y+12,10,c);DrawEllipseLines(x+12,y+12,17,6,c);return;}
 DrawCircle(x+12,y+15,9,c);DrawRectangle(x+8,y+3,8,5,c);DrawCircleLines(x+18,y+3,4,c);
}
static void fp_hud(int i,int w,int h){
 if(!fp_active())return;
 if(fpBannerTimer>0 && fpBannerRead>0){
  int type=fpBannerQueue[(fpBannerRead-1)%32],r=fp_rule(type),size=h<500?18:25;
  const char *text=TextFormat("WEAPON UNLOCKED!  %s  [%c]",pickup_type_label(type),r>=0?fpRules[r].tier:'?');
  while(size>12&&MeasureText(text,size)>w-90)size--;
  int width=MeasureText(text,size)+48,x=(w-width)/2,y=46;
  DrawRectangle(x-8,y-6,width+16,38,Fade(BLACK,.85f));fp_icon(type,x+2,y,GOLD);DrawText(text,x+42,y+3,size,GOLD);
 }
 if(fpRestartTimer>0){const char *text=TextFormat("Checkpoint: restarting wave %d in %.0f",firefightWave,ceilf(fpRestartTimer));int size=w<700?17:28;DrawText(text,(w-MeasureText(text,size))/2,h/2,size,GOLD);}
 int n=fp_near(i,false);if(n>=0 && eq_pickup_gun(fpItems[n].type)>=0){Equipment *e=&equipment[i];int gun=eq_pickup_gun(fpItems[n].type);if(e->gun[0]!=gun&&e->gun[1]!=gun&&e->gun[0]!=EQ_EMPTY&&e->gun[1]!=EQ_EMPTY){
  const char *key=playerInput[i]==INPUT_TYPE_GAMEPAD?"Hold X":i==0?"B":"Y";int size=w<700?12:18;const char *text=TextFormat("%s: replace with %s",key,pickup_type_label(fpItems[n].type));DrawText(text,(w-MeasureText(text,size))/2,h/2+55,size,WHITE);
 }}
}
static void fp_radar(int i,int cx,int cy,int radius,float range){
 if(!fp_active())return;Vector3 f=player_forward(&players[i]);f.y=0;f=v_norm(f);Vector3 right={-f.z,0,f.x};
 for(size_t n=0;n<fpCount;n++)if(fpItems[n].owner==-1){Vector3 d=world_delta(fpItems[n].pos,players[i].pos);d.y=0;if(v_length(d)>range)continue;
  int x=cx+(int)(v_dot(d,right)*radius/range),y=cy-(int)(v_dot(d,f)*radius/range);DrawRectangle(x-2,y-2,4,4,GOLD);
 }
}
#include "firefight_carriers.inc"
#include "firefight_checkpoint.inc"
#include "firefight_progression_net.inc"
