#include <ctype.h>

static void enemy_lab_spawn(int slot,EnemyType type,Vector3 position);
/* Local arena editor. All edits affect cmDocument, never a live playtest. */
typedef struct {int tool,item,brush[3],color,material,panel,row,group,trap,objectKind,objectIndex;bool reinforced,selected,corner,stroke,latched;Vector3 a,b;float distance,yaw,repeat;Vector3 preview,lastPos;float lastYaw,lastPitch;bool previewValid;} CeSeat;
static CeSeat ceSeats[4];
static char ceFileSummary[280][160];
static RenderTexture2D ceThumb[4];
static int ceThumbKind[4]={-99,-99,-99,-99};
static struct {bool active,dirty,testing,lobby,typing,sandbox,physicsPreview;float physicsTime;int menu,row,pending,textAction,textKey,sandboxEnemy,sandboxGun;char text[CM_PATH],path[CM_PATH],notice[256];CmMap *history[128];int historyCount,cursor;size_t bytes;Player cameras[4];int seats;FilePathList files;} ce;
static const char *ceTools[]={"Terrain","Weapons","Upgrades","Supplies","Spawns","Traps","Ducks","Selection"};
static const char *ceSpawns[]={"Player spawn","Enemy entrance","Camp","Boss spawn"};
static const int ceWeapons[]={PICKUP_DYNAMIC_SHOT,PICKUP_SHOTGUN,PICKUP_CUTTER,PICKUP_RAILGUN,PICKUP_ROCKET,PICKUP_SWORD,PICKUP_GLOB,PICKUP_GLOCK,PICKUP_DMR};
static const int ceUpgrades[]={PICKUP_GOLD_TETHER,PICKUP_ADHESIVE,PICKUP_REINFORCED,PICKUP_BLACK_HOLE,PICKUP_IGLOO,PICKUP_SHIELD_REGEN};
static const int ceSupplies[]={PICKUP_AMMO,PICKUP_HEALTH,PICKUP_GRENADE,PICKUP_VOID};
static bool ce_mouse_blocked(void){return ce.active&&!ce.testing&&(ce.menu|| (mouseSlot>=0&&mouseSlot<4&&ceSeats[mouseSlot].panel));}
static void ce_notice(const char *s){snprintf(ce.notice,sizeof(ce.notice),"%s",s);}
static size_t ce_size(const CmMap *m){return sizeof(*m)+(size_t)m->capacity*sizeof(CmCell);}
static void ce_clear_history(void){for(int n=0;n<ce.historyCount;n++)cm_free(ce.history[n]);ce.historyCount=ce.cursor=0;ce.bytes=0;}
static void ce_history(void){
 while(ce.historyCount>ce.cursor){CmMap *m=ce.history[--ce.historyCount];ce.bytes-=ce_size(m);cm_free(m);}
 CmMap *m=cm_clone(cmDocument);if(!m){ce_notice("Cannot allocate undo snapshot");return;}
 while(ce.historyCount&&(ce.historyCount>=128||ce.bytes+ce_size(m)>64u*1024u*1024u)){ce.bytes-=ce_size(ce.history[0]);cm_free(ce.history[0]);memmove(ce.history,ce.history+1,(--ce.historyCount)*sizeof(CmMap*));if(ce.cursor)ce.cursor--;}
 ce.history[ce.historyCount++]=m;ce.cursor=ce.historyCount;ce.bytes+=ce_size(m);
}
static void ce_normalize_traps(void){for(int t=cmDocument->trapCount;t>=1;t--){bool exists=false;for(int n=0;n<cmDocument->count;n++)exists|=cmDocument->cells[n].trap==t;if(exists)continue;for(int n=0;n<cmDocument->count;n++)if(cmDocument->cells[n].trap>t)cmDocument->cells[n].trap--;memmove(cmDocument->seeds+t-1,cmDocument->seeds+t,(cmDocument->trapCount-t)*sizeof(unsigned));cmDocument->trapCount--;}}
static void ce_commit(void){ce_normalize_traps();qsort(cmDocument->cells,cmDocument->count,sizeof(CmCell),cm_compare);ce.dirty=true;ce_history();}
static void ce_undo(int seat,bool redo){for(int i=0;i<4;i++)if(ceSeats[i].stroke){ce_commit();ceSeats[i].stroke=false;}int next=ce.cursor+(redo?1:-1);if(next<1||next>ce.historyCount)return;CmMap *m=cm_clone(ce.history[next-1]);if(!m)return;cm_free(cmDocument);cmDocument=m;ce.cursor=next;for(int n=0;n<4;n++)ceSeats[n].objectKind=0;cm_apply(m,true);ce.dirty=true;ce_notice(TextFormat("P%d %s shared edit history",seat+1,redo?"redid":"undid"));}
static int ce_item_count(int tool){const int c[]={2,9,6,4,4,7,1,5};return c[clampi(tool,0,7)];}
static int ce_pickup(CeSeat *s){return s->tool==1?ceWeapons[s->item%9]:s->tool==2?ceUpgrades[s->item%6]:ceSupplies[s->item%4];}
static const char *ce_item_name(CeSeat *s){if(s->tool==0)return s->material?"Fluid":"Solid";if(s->tool<4)return pickup_type_label(ce_pickup(s));if(s->tool==4)return ceSpawns[s->item%4];if(s->tool==5){static const char *n[]={"Select trap region","Generate branching crack","Paint removed + trigger cells","Paint body cells","Drop slab preset","Cracked pillar preset","Cracked bay preset"};return n[s->item%7];}if(s->tool==6)return "Duck home";static const char *n[]={"Select two corners","Move selection here","Copy selection here","Rotate 90 degrees","Delete selection"};return n[s->item%5];}
static Vector3 ce_target(int i){CeSeat *s=&ceSeats[i];Player *p=&players[i];if(s->previewValid&&Vector3Equals(s->lastPos,p->pos)&&s->lastYaw==p->yaw&&s->lastPitch==p->pitch)return s->preview;
 Ray ray={p->pos,player_forward(p)};VoxelHit h;Vector3 at=v_add(p->pos,v_mul(ray.direction,s->distance));if(first_voxel_hit_detailed(ray,s->distance,-1,true,&h))at=v_add((Vector3){(h.x+.5f)*VOXEL_SIZE,(h.y+.5f)*VOXEL_SIZE,(h.z+.5f)*VOXEL_SIZE},v_mul(h.normal,s->tool==0?VOXEL_SIZE:0));else if(ray.direction.y<-.001f){float t=-p->pos.y/ray.direction.y;if(t>=0&&t<s->distance)at=v_add(p->pos,v_mul(ray.direction,t));}
 at=world_wrap(at);at.x=(floorf(at.x/VOXEL_SIZE)+.5f)*VOXEL_SIZE;at.y=fmaxf(.25f,(floorf(at.y/VOXEL_SIZE)+.5f)*VOXEL_SIZE);at.z=(floorf(at.z/VOXEL_SIZE)+.5f)*VOXEL_SIZE;s->lastPos=p->pos;s->lastYaw=p->yaw;s->lastPitch=p->pitch;s->preview=at;s->previewValid=true;return at;
}
static bool ce_in_box(CmCell c,Vector3 a,Vector3 b){Vector3 p={(c.x+.5f)*VOXEL_SIZE,(c.y+.5f)*VOXEL_SIZE,(c.z+.5f)*VOXEL_SIZE};return p.x>=fminf(a.x,b.x)-.01f&&p.x<=fmaxf(a.x,b.x)+.01f&&p.y>=fminf(a.y,b.y)-.01f&&p.y<=fmaxf(a.y,b.y)+.01f&&p.z>=fminf(a.z,b.z)-.01f&&p.z<=fmaxf(a.z,b.z)+.01f;}
static void ce_add_cell(CmMap *m,int x,int y,int z,Color c,int material,bool reinforced,int trap,int flags){m->cells[m->count++]=(CmCell){x,y,z,c,material,reinforced,trap,flags};}
static void ce_refresh(void){ce_normalize_traps();cm_apply(cmDocument,true);for(int i=0;i<4;i++)ceSeats[i].previewValid=false;}
static void ce_reset_physics(void){
 if(!ce.physicsPreview)return;
 ce_refresh();
 ce.physicsPreview=false;
 ce.physicsTime=0;
 ce_notice("Structure restored");
}

static bool ce_release_structure(int first){
 // Creative maps use unit cells. Follow shared faces, leaving authored traps independent.
 static UnitVoxelBuffer cells;
 static bool seen[MAX_VOXELS];
 memset(seen,0,(size_t)voxel_count*sizeof(*seen));
 for(int n=0;n<cmDocument->count;n++){
  CmCell *cell=&cmDocument->cells[n];
  if(!cell->trap)continue;
  int id=table_get_static_only(cell->x,cell->y,cell->z);
  if(id>=0)seen[id]=true;
 }
 static const int offsets[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
 cells.count=0;
 Voxel *v=&voxels[first];
 unit_voxel_buffer_push(&cells,v->gx,v->gy,v->gz,v->color,v->type,false,first,0,-1);
 seen[first]=true;
 for(int head=0;head<cells.count;head++){
  UnitVoxelSeed *seed=&cells.voxels[head];
  for(int d=0;d<6;d++){
   int id=table_get_static_only(seed->gx+offsets[d][0],seed->gy+offsets[d][1],seed->gz+offsets[d][2]);
   if(id<0||seen[id]||voxels[id].type!=0)continue;
   seen[id]=true;
   v=&voxels[id];
   unit_voxel_buffer_push(&cells,v->gx,v->gy,v->gz,v->color,v->type,false,id,0,-1);
  }
 }
 // Descending removal keeps lower indices valid as the voxel array compacts.
 ++staticBeliefBatchDepth;
 static int indices[MAX_VOXELS];
 for(int n=0;n<cells.count;n++)indices[n]=cells.voxels[n].voxelIndex;
 qsort(indices,cells.count,sizeof(*indices),activation_index_descending);
 for(int n=0;n<cells.count;n++)remove_voxel_index(indices[n]);
 --staticBeliefBatchDepth;
 int spawned=emit_unit_voxels_from_units(&cells,false,true,-1);
 if(spawned!=cells.count){ce_refresh();return false;}
 rebuild_voxel_hash();
 refresh_activation_surfaces(&cells);
 refresh_static_voxel_beliefs();
 return true;
}

static void ce_toggle_physics(int seat){
 if(ce.physicsPreview){ce_reset_physics();return;}
 Player *p=&players[seat];
 VoxelHit hit;
 if(!first_voxel_hit_detailed((Ray){p->pos,player_forward(p)},cmDocument->width*2,-1,true,&hit)||voxels[hit.id].type!=0){
  ce_notice("Aim at a solid structure");return;
 }
 bool released=false,trap=false;
 for(int n=0;n<foundryTrapCount;n++){
  if(!foundry_explicit_flags(&foundryTraps[n],hit.x,hit.y,hit.z))continue;
  foundry_trap_release(n);
  released=foundryTraps[n].active;
  trap=true;
  break;
 }
 if(!trap)released=ce_release_structure(hit.id);
 ce.physicsPreview=released;
 ce.physicsTime=0;
 ce_notice(released?"Physics running":"Insufficient physics capacity");
}

static void ce_step_physics(float dt){
 if(!ce.physicsPreview)return;
 ce.physicsTime=fminf(ce.physicsTime+dt,PBD_MAX_STEP_DT*PBD_MAX_ACCUM_STEPS);
 int steps=(int)(ce.physicsTime/PBD_MAX_STEP_DT);
 if(!steps)return;
 ce.physicsTime-=steps*PBD_MAX_STEP_DT;
 simulate_voxel_pbd_steps(PBD_MAX_STEP_DT,steps);
 for(int i=0;i<4;i++)ceSeats[i].previewValid=false;
}
static bool ce_edit(int i,bool erase){
 CeSeat *s=&ceSeats[i];CmMap *m=cmDocument;Vector3 at=ce_target(i);if(erase){for(int n=0;n<4;n++)ceSeats[n].objectKind=0;}if(erase&&s->tool==0){VoxelHit hit;if(first_voxel_hit_detailed((Ray){players[i].pos,player_forward(&players[i])},s->distance,-1,true,&hit))at=(Vector3){(hit.x+.5f)*VOXEL_SIZE,(hit.y+.5f)*VOXEL_SIZE,(hit.z+.5f)*VOXEL_SIZE};}if(!cm_inside(m,at)){ce_notice("Outside world bounds");return false;}
 if(s->tool==7 || (s->tool==5&&s->item==0)){
  if(s->tool==5||s->item==0){if(!s->corner){s->a=at;s->corner=true;s->selected=false;}else{s->b=at;s->corner=false;s->selected=true;}return false;}
  if(!s->selected){ce_notice("Select two corners first");return false;}
  int count=0;for(int n=0;n<m->count;n++)count+=ce_in_box(m->cells[n],s->a,s->b);if(!count)return false;
  if(s->item==4||erase){int out=0;for(int n=0;n<m->count;n++)if(!ce_in_box(m->cells[n],s->a,s->b))m->cells[out++]=m->cells[n];m->count=out;}
  else {CmMap *copy=cm_clone(m);if(!copy)return false;int original=m->count;if(!cm_reserve(copy,original+(s->item==2?count:0))){cm_free(copy);ce_notice("Voxel budget exceeded");return false;}Vector3 center=v_mul(v_add(s->a,s->b),.5f),delta=v_sub(at,center);int out=0;
   int remap[FOUNDRY_TRAPS+1]={0};for(int t=1;t<=m->trapCount;t++){int total=0,selected=0;for(int n=0;n<original;n++)if(m->cells[n].trap==t){total++;selected+=ce_in_box(m->cells[n],s->a,s->b);}if(selected&&selected!=total){cm_free(copy);ce_notice("Select the entire linked trap to transform it");return false;}if(selected&&s->item==2){if(copy->trapCount>=FOUNDRY_TRAPS){cm_free(copy);ce_notice("Trap capacity reached");return false;}remap[t]=++copy->trapCount;copy->seeds[remap[t]-1]=m->seeds[t-1];}}
   for(int n=0;n<original;n++){CmCell c=m->cells[n];if(!ce_in_box(c,s->a,s->b)){copy->cells[out++]=c;continue;}if(s->item==2)copy->cells[out++]=c;
    if(c.trap&&remap[c.trap])c.trap=remap[c.trap];
    if(s->item==3){int cx=(int)floorf(center.x/VOXEL_SIZE),cz=(int)floorf(center.z/VOXEL_SIZE),x=c.x-cx;c.x=cx-(c.z-cz);c.z=cz+x;}else{c.x+=(int)roundf(delta.x/VOXEL_SIZE);c.y+=(int)roundf(delta.y/VOXEL_SIZE);c.z+=(int)roundf(delta.z/VOXEL_SIZE);}copy->cells[out++]=c;
   }copy->count=out;if(!cm_validate(copy)){cm_free(copy);ce_notice(cmError);return false;}cm_free(cmDocument);cmDocument=m=copy;
  }ce_refresh();return true;
 }
 if(s->tool>=1&&s->tool<=3){if(!erase&&s->objectKind==1&&s->objectIndex<m->pickupCount){Pickup *p=&m->supplies[s->objectIndex];p->pos=v_add(at,(Vector3){0,.5f,0});p->type=ce_pickup(s);ce_refresh();return true;}if(erase){int nearest=-1;float d=2;for(int n=0;n<m->pickupCount;n++){float x=v_length(world_delta(m->supplies[n].pos,at));if(x<d)d=x,nearest=n;}if(nearest<0)return false;memmove(m->supplies+nearest,m->supplies+nearest+1,(--m->pickupCount-nearest)*sizeof(Pickup));}else{if(m->pickupCount==MAX_PICKUPS){ce_notice("Pickup capacity reached");return false;}m->supplies[m->pickupCount++]=(Pickup){.pos=v_add(at,(Vector3){0,.5f,0}),.active=true,.type=ce_pickup(s)};}ce_refresh();return true;}
 if(s->tool==4){if(!erase&&s->objectKind==2&&s->objectIndex<m->spawnCount){CmSpawn *p=&m->spawns[s->objectIndex];p->pos=at;if(s->item==CM_PLAYER||s->item==CM_ENTRANCE)p->pos.y+=BASE_EYE_HEIGHT;p->kind=s->item;p->group=s->group;p->yaw=s->yaw;return true;}if(erase){int nearest=-1;float d=3;for(int n=0;n<m->spawnCount;n++){float x=v_length(world_delta(m->spawns[n].pos,at));if(x<d)d=x,nearest=n;}if(nearest<0)return false;memmove(m->spawns+nearest,m->spawns+nearest+1,(--m->spawnCount-nearest)*sizeof(CmSpawn));}else{if(m->spawnCount==CM_SPAWNS)return false;Vector3 pos=at;if(s->item==CM_PLAYER||s->item==CM_ENTRANCE)pos.y+=BASE_EYE_HEIGHT;m->spawns[m->spawnCount++]=(CmSpawn){pos,s->yaw,s->item,s->group};}return true;}
 if(s->tool==6){if(!erase&&s->objectKind==3&&s->objectIndex<m->duckCount){at.y=0;m->ducks[s->objectIndex]=at;return true;}if(erase){int nearest=-1;float d=4;for(int n=0;n<m->duckCount;n++){float x=v_length(world_delta(m->ducks[n],at));if(x<d)d=x,nearest=n;}if(nearest<0)return false;memmove(m->ducks+nearest,m->ducks+nearest+1,(--m->duckCount-nearest)*sizeof(Vector3));}else{if(m->duckCount==MAP_DUCKS){ce_notice("Four duck homes maximum");return false;}at.y=0;m->ducks[m->duckCount++]=at;}return true;}
 if(s->tool==5){
  if(s->item>=4){int sx=s->item==4?8:s->item==5?8:16,sy=s->item==4?8:s->item==5?20:12,sz=s->item==4?8:s->item==5?8:16;int need=sx*sy*sz;if(m->count+need>MAX_VOXELS||!cm_reserve(m,m->count+need)){ce_notice("Voxel capacity reached");return false;}
   int x0=floorf(at.x/VOXEL_SIZE),y0=floorf(at.y/VOXEL_SIZE),z0=floorf(at.z/VOXEL_SIZE);s->a=(Vector3){(x0+.5f)*VOXEL_SIZE,(y0+.5f)*VOXEL_SIZE,(z0+.5f)*VOXEL_SIZE};s->b=v_add(s->a,(Vector3){(sx-1)*VOXEL_SIZE,(sy-1)*VOXEL_SIZE,(sz-1)*VOXEL_SIZE});if(!cm_inside(m,s->b)){ce_notice("Preset outside world");return false;}
   for(int x=0;x<sx;x++)for(int y=0;y<sy;y++)for(int z=0;z<sz;z++){bool put=s->item==4?(y==sy-1||((x==0||x==sx-1)&&(z==0||z==sz-1))):s->item==5?true:(x==0||x==sx-1||z==0||z==sz-1||y==sy-1);if(put&&table_get_static_only(x0+x,y0+y,z0+z)<0)ce_add_cell(m,x0+x,y0+y,z0+z,creative_block_palette(s->color),0,false,0,0);}
   s->selected=true;ce_refresh();
  }
  if(s->item==1||s->item>=4){if(!s->selected){ce_notice("Select trap region first");return false;}int count=0;for(int n=0;n<m->count;n++)if(ce_in_box(m->cells[n],s->a,s->b))count++;if(count>FOUNDRY_TRAP_PAYLOAD||count<2||m->trapCount>=FOUNDRY_TRAPS){ce_notice("Trap needs 2–2048 cells; maximum 32 traps");return s->item>=4;}
   for(int n=0;n<m->count;n++)if(ce_in_box(m->cells[n],s->a,s->b)&&m->cells[n].trap){ce_notice("Region already belongs to a trap");return s->item>=4;}
   int t=++m->trapCount;s->trap=t;unsigned seed=0x12345u+t*9781u;m->seeds[t-1]=seed;int first=-1;
   bool crackMask[64][64]={0};foundry_crack_path(crackMask,0,30,32,32,&seed);foundry_crack_path(crackMask,32,32,63,30,&seed);foundry_crack_path(crackMask,32,32,12,63,&seed);foundry_crack_path(crackMask,32,32,52,0,&seed);
   float minx=fminf(s->a.x,s->b.x),miny=fminf(s->a.y,s->b.y),spanx=fmaxf(VOXEL_SIZE,fabsf(s->b.x-s->a.x)),spany=fmaxf(VOXEL_SIZE,fabsf(s->b.y-s->a.y));
   for(int n=0;n<m->count;n++){CmCell *c=&m->cells[n];if(!ce_in_box(*c,s->a,s->b))continue;if(first<0)first=n;c->trap=t;c->flags=CM_BODY;c->material=0;int u=clampi((int)roundf(((c->x+.5f)*VOXEL_SIZE-minx)/spanx*63),0,63),v=clampi((int)roundf(((c->y+.5f)*VOXEL_SIZE-miny)/spany*63),0,63);bool crack=false;int pad=(int)ceilf(32*VOXEL_SIZE/fminf(spanx,spany));for(int x=clampi(u-pad,0,63);x<=clampi(u+pad,0,63);x++)for(int y=clampi(v-pad,0,63);y<=clampi(v+pad,0,63);y++)crack|=crackMask[x][y];if(crack){c->flags=CM_REMOVE|CM_TRIGGER;c->color=ORANGE;}}
   if(first>=0){m->cells[first].flags=CM_REMOVE|CM_TRIGGER;m->cells[first].color=ORANGE;for(int n=m->count-1;n>=0;n--)if(n!=first&&m->cells[n].trap==t){m->cells[n].flags=CM_BODY;m->cells[n].color=creative_block_palette(s->color);break;}}
  }else {if(s->trap<1||s->trap>m->trapCount){ce_notice("Choose a trap in Properties");return false;}CmCell key={.x=floorf(at.x/VOXEL_SIZE),.y=floorf(at.y/VOXEL_SIZE),.z=floorf(at.z/VOXEL_SIZE)};CmCell *c=bsearch(&key,m->cells,m->count,sizeof(CmCell),cm_compare);if(!c)return false;if(c->trap&&c->trap!=s->trap)return false;c->trap=s->trap;c->flags=s->item==2?CM_REMOVE|CM_TRIGGER:CM_BODY;c->color=s->item==2?ORANGE:creative_block_palette(s->color);}
  ce_refresh();return true;
 }
 // Terrain brush: batch the document change, then update only affected runtime cells.
 int sx=s->brush[0],sy=s->brush[1],sz=s->brush[2],x0=floorf(at.x/VOXEL_SIZE)-sx/2,y0=floorf(at.y/VOXEL_SIZE)-sy/2,z0=floorf(at.z/VOXEL_SIZE)-sz/2;
 Vector3 a={(x0+.5f)*VOXEL_SIZE,(y0+.5f)*VOXEL_SIZE,(z0+.5f)*VOXEL_SIZE},b=v_add(a,(Vector3){(sx-1)*VOXEL_SIZE,(sy-1)*VOXEL_SIZE,(sz-1)*VOXEL_SIZE});
 if(!cm_inside(m,a)||!cm_inside(m,b)){ce_notice("Brush outside world");return false;}int existing=0;for(int n=0;n<m->count;n++)existing+=ce_in_box(m->cells[n],a,b);int final=m->count-existing+(erase?0:sx*sy*sz);if(!cm_reserve(m,final)){ce_notice("Voxel budget exceeded");return false;}
 if(erase&&!existing)return false;int out=0;for(int n=0;n<m->count;n++)if(!ce_in_box(m->cells[n],a,b))m->cells[out++]=m->cells[n];m->count=out;
 for(int x=x0;x<x0+sx;x++)for(int y=y0;y<y0+sy;y++)for(int z=z0;z<z0+sz;z++){int id=table_get_static_only(x,y,z);if(erase){if(id>=0)remove_voxel_index(id);}else{ce_add_cell(m,x,y,z,creative_block_palette(s->color),s->material,s->reinforced,0,0);int v=add_static_voxel_at_grid(x,y,z,creative_block_palette(s->color),s->material);if(v>=0)voxels[v].reinforced=s->reinforced;}}
 qsort(m->cells,m->count,sizeof(CmCell),cm_compare);rebuild_static_hash_if_dirty();rebuild_voxel_hash();meshDirty=true;for(int k=0;k<4;k++)ceSeats[k].previewValid=false;return true;
}
static void ce_defaults(CmMap *m){for(int i=0;i<4;i++)m->spawns[m->spawnCount++]=(CmSpawn){{(i-1.5f)*3,BASE_EYE_HEIGHT,5},180,CM_PLAYER,0};m->spawns[m->spawnCount++]=(CmSpawn){{0,0,0},0,CM_CAMP,0};m->spawns[m->spawnCount++]=(CmSpawn){{0,0,-12},0,CM_BOSS,0};for(int i=0;i<2;i++)m->spawns[m->spawnCount++]=(CmSpawn){{(i?1:-1)*(m->width*.5f-4),BASE_EYE_HEIGHT,0},0,CM_ENTRANCE,i};}
static void ce_accept(CmMap *m,const char *path){if(!cm_apply(m,true)){cm_free(m);ce_notice(cmError);return;}gameMode=GAME_MODE_DEATHMATCH;cm_free(cmDocument);cmDocument=m;snprintf(ce.path,CM_PATH,"%s",path?path:"");snprintf(cmSelectedPath,CM_PATH,"%s",ce.path);ce_clear_history();ce_history();ce.dirty=false;for(int i=0;i<4;i++){ceSeats[i].previewValid=false;ceSeats[i].objectKind=0;}ce.menu=0;ce.row=0;ce_notice("Map ready. Y / Tab opens the palette.");}
static void ce_begin(void){
 ce.active=true;ce.testing=cmTesting=false;eq_reset_all();ce.menu=1;ce.row=0;ce.seats=1;netRequestedCreative=false;enemyLab=false;useCustomMap=false;
 int device=couch_menu_device();multiplayerActivePlayers=activePlayers=1;playerInput[0]=multiplayerPlayerInput[0]=device>=0?INPUT_TYPE_GAMEPAD:INPUT_TYPE_KEYBOARD;if(device>=0)localSeats[0].device=device;
 for(int i=0;i<4;i++){ceSeats[i]=(CeSeat){.brush={3,3,3},.distance=8};players[i].matterMax=players[i].matter=100;players[i].respawn_timer=0;}
 if(!cmDocument){cmDocument=cm_new();ce_defaults(cmDocument);ce_history();}ce_refresh();ce_notice("Choose a map. Other controllers: Start to join. F9 adds keyboard P2.");
}
static int ce_path_compare(const void *a,const void *b){return strcmp(*(const char *const*)a,*(const char *const*)b);}
static void ce_scan(void){if(ce.files.paths)UnloadDirectoryFiles(ce.files);MakeDirectory("maps");ce.files=LoadDirectoryFilesEx("maps",".map",false);qsort(ce.files.paths,ce.files.count,sizeof(char*),ce_path_compare);ce.row=0;for(unsigned n=0;n<ce.files.count&&n<280;n++){CmMap *m=cm_read(ce.files.paths[n]);int players=0,entrances=0,camps=0,boss=0;if(m)for(int k=0;k<m->spawnCount;k++){players+=m->spawns[k].kind==CM_PLAYER;entrances+=m->spawns[k].kind==CM_ENTRANCE;camps+=m->spawns[k].kind==CM_CAMP;boss+=m->spawns[k].kind==CM_BOSS;}snprintf(ceFileSummary[n],sizeof(ceFileSummary[n]),"%s | %s",GetFileName(ce.files.paths[n]),!m?"INVALID":players&&entrances&&camps==1&&boss==1?"Versus / Firefight markers":players?"Versus markers":"Needs spawns");cm_free(m);}}
static void ce_library(bool lobby){ce.lobby=lobby;ce.menu=3;ce_scan();}
static bool ce_ready(bool firefight){if(!cm_validate(cmDocument)){ce_notice(cmError);return false;}int playersCount=0,entrances=0,camps=0,bosses=0;for(int n=0;n<cmDocument->spawnCount;n++){CmSpawn *s=&cmDocument->spawns[n];playersCount+=s->kind==CM_PLAYER;entrances+=s->kind==CM_ENTRANCE;camps+=s->kind==CM_CAMP;bosses+=s->kind==CM_BOSS;if(s->kind==CM_PLAYER||s->kind==CM_ENTRANCE){if(player_aabb_hits_world_ex(s->pos,PLAYER_SIZE*.5f,-1)){ce_notice("Spawn marker intersects terrain");return false;}}if(s->kind==CM_BOSS&&firefight){BoundingBox b={v_add(s->pos,(Vector3){-3,.1f,-3}),v_add(s->pos,(Vector3){3,10,3})};for(int k=0;k<cmDocument->count;k++){CmCell *c=&cmDocument->cells[k];Vector3 p={(c->x+.5f)*VOXEL_SIZE,(c->y+.5f)*VOXEL_SIZE,(c->z+.5f)*VOXEL_SIZE};if(CheckCollisionBoxSphere(b,p,.25f)){ce_notice("Boss marker needs clear 6 x 10 x 6 metre space");return false;}}}}
 if(playersCount<ce.seats||(firefight&&(!entrances||camps!=1||bosses!=1))){ce_notice("Need player spawns, plus entrances, one camp and one boss for Firefight");return false;}if(firefight){Vector3 stand[4];int old=firefightHumans;firefightHumans=4;bool valid=false;for(int n=0;n<cmDocument->spawnCount;n++)if(cmDocument->spawns[n].kind==CM_CAMP)valid=camp_site(cmDocument->spawns[n].pos,stand);firefightHumans=old;camp_end();if(!valid){ce_notice("Camp needs clear flat space for four players");return false;}}return true;}
static bool ce_validate_selected(bool firefight){CmMap *m=cm_read(cmSelectedPath);if(!m){ce_notice(cmError);return false;}CmMap *saved=cmDocument;int seats=ce.seats;cmDocument=m;ce.seats=multiplayerActivePlayers;bool okay=cm_apply(m,true)&&ce_ready(firefight);cmDocument=saved;ce.seats=seats;if(!okay)cm_fail(ce.notice);else cmError[0]=0;if(okay){cm_free(cmRuntime);cmRuntime=m;}else cm_free(m);return okay;}
static void ce_play(int mode){if(!ce_ready(mode==2))return;for(int i=0;i<4;i++)ce.cameras[i]=players[i];ce.seats=activePlayers;ce.testing=cmTesting=true;useCustomMap=true;gameMode=mode==2?GAME_MODE_FIREFIGHT:GAME_MODE_DEATHMATCH;firefightFoundry=mode==2;firefightPractice=false;enemyLab=false;multiplayerActivePlayers=ce.seats;ce.menu=0;ResetGame();gameState=GAME_STATE_PLAYING;ce.sandbox=mode==0;if(ce.sandbox){gameMode=GAME_MODE_FIREFIGHT;enemyLab=true;enemyLabPaused=false;enemyLabCase=ENEMY_TYPE_STANDARD;firefightHumans=ce.seats;activePlayers=MAX_COMBATANTS;for(int n=0;n<MAX_COMBATANTS;n++){combatantTeam[n]=n<MAX_PLAYERS?0:1;combatantLife[n]=n<ce.seats?LIFE_ALIVE:LIFE_UNUSED;if(n>=ce.seats)players[n].respawn_timer=9999;}ce.menu=7;ce.row=0;}ce_notice("Playtest: F10 / Back+Start returns to editor. Console available for sandbox spawning.");}
static void ce_return(void){enemy_cleanup(true);camp_end();ffend_reset();cmTesting=ce.testing=ce.sandbox=false;creativeModeActive=true;gameMode=GAME_MODE_DEATHMATCH;useCustomMap=false;enemyLab=false;activePlayers=multiplayerActivePlayers=ce.seats;for(int i=0;i<4;i++){players[i]=ce.cameras[i];playerInput[i]=multiplayerPlayerInput[i];}gameState=GAME_STATE_CREATIVE;ce_refresh();ce.menu=0;}
static void ce_text(int action,const char *initial){ce.textAction=action;ce.typing=true;ce.textKey=0;snprintf(ce.text,CM_PATH,"%s",initial?initial:"");}
static bool ce_save(const char *path){if(!cm_validate(cmDocument)){ce_notice(cmError);return false;}if(!cm_write(cmDocument,path)){ce_notice(cmError);return false;}snprintf(ce.path,CM_PATH,"%s",path);snprintf(cmSelectedPath,CM_PATH,"%s",path);ce.dirty=false;ce_notice("Map saved");return true;}
static void ce_action(int action){
 if(action==70){ct_open();return;}
 if((action==10||action==11||action==12||action==20||action==30||action==42)&&ce.dirty){ce.pending=action;ce.menu=6;ce.row=0;return;}
 if(action>=10&&action<=12){CmMap *m=NULL;if(action==10){m=cm_new();ce_defaults(m);}else{clear_world_voxels();clear_pickups();worldFloorHalfSize=action==11?40:32;worldPeriodic=true;useCustomMap=false;currentWorldType=WORLD_TYPE_CROSSCUT;gameMode=action==11?GAME_MODE_FIREFIGHT:GAME_MODE_DEATHMATCH;firefightFoundry=action==11;if(action==11)buildFoundryHoldout();else buildCrosscut();init_pickups();m=cm_capture(action==11?"Foundry copy":"Crosscut copy");if(m){ce_defaults(m);for(int n=0;n<4;n++)m->spawns[n].pos=pick_player_spawn(n);
   Vector3 standing[4];int old=firefightHumans;firefightHumans=4;for(int n=0;n<161;n++){float angle=n*2*PI/16,radius=n?2*(1+(n-1)/16):0;Vector3 pos={cosf(angle)*radius,0,sinf(angle)*radius};if(camp_site(pos,standing)){m->spawns[4].pos=camp.pos;break;}}firefightHumans=old;camp_end();
   for(int n=0;n<64;n++){Vector3 pos={(n&1?1:-1)*(m->width*.5f-6),.04f,(n/2-16)*1.f};if(go_move_clear(pos)){m->spawns[5].pos=pos;break;}}
   for(int n=0;n<4;n++)m->ducks[m->duckCount++]=(Vector3){(n&1?1:-1)*(m->width*.5f-5),0,(n&2?1:-1)*10};}}if(m)ce_accept(m,NULL);}
 else if(action==20)ce_library(false);else if(action==30){ce.active=false;creativeModeActive=false;gameState=GAME_STATE_MENU;ce.menu=0;}else if(action==40){if(ce.path[0])ce_save(ce.path);else ce_text(1,cmDocument->name);}else if(action==41)ce_text(1,cmDocument->name);else if(action==42)ce_text(2,"");else if(action==43)ce_text(3,ce.path);else if(action>=50&&action<=52)ce_play(action-50);
}
static bool ce_pad_press(int i,int button){return playerInput[i]==INPUT_TYPE_GAMEPAD&&IsGamepadAvailable(seat_device(i))&&IsGamepadButtonPressed(seat_device(i),button);}
static bool ce_pad_down(int i,int button){return playerInput[i]==INPUT_TYPE_GAMEPAD&&IsGamepadAvailable(seat_device(i))&&IsGamepadButtonDown(seat_device(i),button);}
static void ce_property(CeSeat *s,int row,int delta){switch(row){case 0:s->objectKind=0;s->item=(s->item+delta+ce_item_count(s->tool))%ce_item_count(s->tool);if(s->tool==0)s->material=s->item;break;case 1:case 2:case 3:s->brush[row-1]=clampi(s->brush[row-1]+delta,1,20);break;case 4:s->distance=clampf(s->distance+delta,1,CM_MAX_WIDTH);break;case 5:s->color=(s->color+delta+6)%6;break;case 6:s->material=!s->material;break;case 7:s->reinforced=!s->reinforced;break;case 8:s->yaw=remainderf(s->yaw+delta*90,360);break;case 9:s->group=clampi(s->group+delta,0,31);break;case 10:s->trap=clampi(s->trap+delta,1,cmDocument->trapCount?cmDocument->trapCount:1);break;}s->previewValid=false;}
static void ce_update(float dt){
 if(ce.menu||ce.testing)return;
 if(!IsWindowFocused()){for(int n=0;n<4;n++){if(ceSeats[n].stroke)ce_commit();ceSeats[n].stroke=ceSeats[n].latched=false;}return;}
 for(int d=0;d<4;d++)if(IsGamepadAvailable(d)&&IsGamepadButtonPressed(d,GAMEPAD_BUTTON_MIDDLE_RIGHT)){bool used=false;for(int i=0;i<activePlayers;i++)used|=playerInput[i]==INPUT_TYPE_GAMEPAD&&seat_device(i)==d;if(!used&&activePlayers<4){int i=activePlayers++;ce.seats=multiplayerActivePlayers=activePlayers;localSeats[i].device=d;playerInput[i]=multiplayerPlayerInput[i]=INPUT_TYPE_GAMEPAD;players[i].pos=(Vector3){0,5,0};return;}}
 int keyboards=0;for(int n=0;n<activePlayers;n++)keyboards+=playerInput[n]==INPUT_TYPE_KEYBOARD;
 if(ct_pressed(0,CT_JOIN_EDITOR)&&activePlayers<4&&keyboards<2){int i=activePlayers++;ce.seats=multiplayerActivePlayers=activePlayers;playerInput[i]=multiplayerPlayerInput[i]=INPUT_TYPE_KEYBOARD;players[i].pos=(Vector3){0,5,0};}
 for(int i=0;i<activePlayers;i++){CeSeat *s=&ceSeats[i];bool pad=playerInput[i]==INPUT_TYPE_GAMEPAD;bool secondary=false;for(int n=0;n<i;n++)secondary|=playerInput[n]==INPUT_TYPE_KEYBOARD;int device=seat_device(i);if(pad&&!IsGamepadAvailable(device)){if(s->stroke)ce_commit();s->stroke=s->latched=false;continue;}
  if(ct_pressed(i,CT_PAUSE)||ce_pad_press(i,GAMEPAD_BUTTON_MIDDLE_RIGHT)||(!pad&&IsKeyPressed(KEY_ESCAPE)&&!s->panel)){if(s->stroke){ce_commit();s->stroke=false;}ce_reset_physics();ce.menu=2;ce.row=0;return;}
  bool palette=ct_pressed(i,CT_PALETTE);bool props=ct_pressed(i,CT_PROPERTIES);if(palette||props){ce_reset_physics();s->panel=s->panel?0:palette?1:2;s->row=0;}
  int vertical=pad?(ce_pad_press(i,GAMEPAD_BUTTON_LEFT_FACE_DOWN)-ce_pad_press(i,GAMEPAD_BUTTON_LEFT_FACE_UP)):(IsKeyPressed(KEY_DOWN)-IsKeyPressed(KEY_UP));int horizontal=pad?(ce_pad_press(i,GAMEPAD_BUTTON_LEFT_FACE_RIGHT)-ce_pad_press(i,GAMEPAD_BUTTON_LEFT_FACE_LEFT)):(IsKeyPressed(KEY_RIGHT)-IsKeyPressed(KEY_LEFT));
  if(s->panel){if(!pad&&mouseSlot==i){int vx,vy,vw,vh;get_viewport(i,activePlayers,&vx,&vy,&vw,&vh);Vector2 p=GetMousePosition();int fs=vw<700?13:18;int row=(int)((p.y-vy-75)/(fs+5));if(p.x>=vx+6&&p.x<vx+vw-6&&p.y>=vy+75&&row>=0&&row<(s->panel==1?8:11)&&IsMouseButtonPressed(MOUSE_BUTTON_LEFT)){s->row=row;if(s->panel==1){s->tool=row;s->item=0;s->objectKind=0;s->previewValid=false;s->panel=2;s->row=0;}else ce_property(s,row,1);}if(GetMouseWheelMove())ce_property(s,s->row,GetMouseWheelMove()>0?1:-1);}if(s->stroke){ce_commit();s->stroke=false;}if(ce_pad_press(i,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT)||(!pad&&IsKeyPressed(KEY_ESCAPE))){s->panel=0;continue;}s->row=(s->row+vertical+(s->panel==1?8:11))%(s->panel==1?8:11);bool confirm=ce_pad_press(i,GAMEPAD_BUTTON_RIGHT_FACE_DOWN)||(!pad&&IsKeyPressed(KEY_ENTER));if(s->panel==1&&confirm){s->tool=s->row;s->item=0;s->objectKind=0;s->previewValid=false;s->panel=2;s->row=0;}else if(s->panel==2&&(horizontal||confirm))ce_property(s,s->row,horizontal?horizontal:1);continue;}
  if(ct_pressed(i,CT_UNDO)||ct_pressed(i,CT_REDO))ce_reset_physics();
  if(ct_pressed(i,CT_UNDO))ce_undo(i,false);if(ct_pressed(i,CT_REDO))ce_undo(i,true);
  int size=ct_pressed(i,CT_SIZE_UP)-ct_pressed(i,CT_SIZE_DOWN);if(size)for(int a=0;a<3;a++)s->brush[a]=clampi(s->brush[a]+size,1,20);
  float distance=ct_pressed(i,CT_DIST_UP)-ct_pressed(i,CT_DIST_DOWN);if(distance){s->distance=clampf(s->distance+distance,1,256);s->previewValid=false;}
  Player *p=&players[i];Vector2 axes=ct_move_axes(i);float mx=axes.x,mz=axes.y,my=ct_down(i,CT_UP)-ct_down(i,CT_DOWN);
  ct_edit_look(i,dt);
  Vector3 f=player_forward(p),r=v_norm(v_cross(f,(Vector3){0,1,0})),move=v_add(v_add(v_mul(f,mz),v_mul(r,mx)),(Vector3){0,my,0});float speed=ct_down(i,CT_FAST)?22:10;if(v_length(move)>1)move=v_norm(move);p->pos=world_wrap(v_add(p->pos,v_mul(move,speed*dt)));p->pos.y=clampf(p->pos.y,.5f,cmDocument->width+10);
  if(ct_pressed(i,CT_ACTIVATE_STRUCTURE)){if(s->stroke){ce_commit();s->stroke=false;}ce_toggle_physics(i);continue;}
  bool sample=ct_pressed(i,CT_SAMPLE);if(sample){
   Ray ray={p->pos,player_forward(p)};float best=s->distance;int kind=0,index=-1;
   for(int category=1;category<=3;category++){int count=category==1?cmDocument->pickupCount:category==2?cmDocument->spawnCount:cmDocument->duckCount;for(int n=0;n<count;n++){Vector3 pos=category==1?cmDocument->supplies[n].pos:category==2?cmDocument->spawns[n].pos:cmDocument->ducks[n];pos=world_image(pos,p->pos);RayCollision h=GetRayCollisionSphere(ray,pos,.8f);if(h.hit&&h.distance<best){best=h.distance;kind=category;index=n;}}}
   s->objectKind=kind;s->objectIndex=index;
   if(kind==1){int type=cmDocument->supplies[index].type;for(int tool=1;tool<=3;tool++)for(int item=0;item<ce_item_count(tool);item++){CeSeat probe={.tool=tool,.item=item};if(ce_pickup(&probe)==type){s->tool=tool;s->item=item;}}}
   else if(kind==2){CmSpawn *spawn=&cmDocument->spawns[index];s->tool=4;s->item=spawn->kind;s->group=spawn->group;s->yaw=spawn->yaw;}
   else if(kind==3){s->tool=6;s->item=0;}
   else {VoxelHit hit;if(first_voxel_hit_detailed(ray,s->distance,-1,true,&hit)){Voxel *v=&voxels[hit.id];s->material=v->type;s->reinforced=v->reinforced;for(int n=0;n<6;n++)if(ColorToInt(creative_block_palette(n))==ColorToInt(v->color))s->color=n;}}
   if(kind)ce_notice("Object selected: aim and apply to move it. Palette selects a new object.");
  }
  bool apply=ct_down(i,CT_PLACE);bool erase=ct_down(i,CT_ERASE);
  s->repeat-=dt;bool discrete=s->tool!=0;bool busy=false;for(int n=0;n<activePlayers;n++)if(n!=i&&ceSeats[n].stroke)busy=true;if(busy){apply=erase=false;}if((apply||erase)&&s->repeat<=0&&(!discrete||!s->latched)){ce_reset_physics();bool changed=ce_edit(i,erase);s->latched=true;s->stroke|=changed;s->repeat=discrete?.25f:.12f;}
  if(!(apply||erase)){if(s->stroke)ce_commit();s->stroke=false;s->latched=false;s->repeat=0;}
 }
 ce_step_physics(dt);
}
static void ce_world(int i){if(!ce.active||ce.testing||!cmDocument||i<0||i>=4)return;CeSeat *s=&ceSeats[i];Vector3 at=ce_target(i);Color c=cm_inside(cmDocument,at)?SKYBLUE:RED;DrawCubeWires(at,s->brush[0]*VOXEL_SIZE,s->brush[1]*VOXEL_SIZE,s->brush[2]*VOXEL_SIZE,c);
 if(s->selected||s->corner){Vector3 a=s->a,b=s->corner?at:s->b;DrawBoundingBox((BoundingBox){v_sub(Vector3Min(a,b),(Vector3){.25f,.25f,.25f}),v_add(Vector3Max(a,b),(Vector3){.25f,.25f,.25f})},GOLD);}
 DrawCubeWires((Vector3){0,cmDocument->width*.5f,0},cmDocument->width,cmDocument->width,cmDocument->width,Fade(SKYBLUE,.3f));
 for(int n=0;n<cmDocument->spawnCount;n++){CmSpawn *p=&cmDocument->spawns[n];Color col=p->kind==CM_PLAYER?GREEN:p->kind==CM_ENTRANCE?RED:p->kind==CM_CAMP?GOLD:PURPLE;DrawCubeWires(p->pos,.7f,1.5f,.7f,col);float a=p->yaw*DEG2RAD;DrawLine3D(p->pos,v_add(p->pos,(Vector3){-sinf(a)*2,0,-cosf(a)*2}),col);}
 for(int n=0;n<cmDocument->duckCount;n++){duck_renderer_init(NULL);DuckPose pose={.position=cmDocument->ducks[n],.scale=.75f,.time=0};duck_draw(pose);DrawSphereWires(v_add(cmDocument->ducks[n],(Vector3){0,.3f,0}),.6f,4,8,GREEN);}
}
static void ce_hud(int i,int w,int h){if(!ce.active||!cmDocument)return;CeSeat *s=&ceSeats[i];int font=w<700?13:18;DrawRectangle(0,0,w,58,Fade(BLACK,.85f));DrawText(TextFormat("P%d %s%s | %d/%d voxels",i+1,cmDocument->name,ce.dirty?" *":"",cmDocument->count,MAX_VOXELS),8,6,font,WHITE);DrawText(TextFormat("%s: %s | %.0fm | %dx%dx%d",ceTools[s->tool],ce_item_name(s),s->distance,s->brush[0],s->brush[1],s->brush[2]),8,28,font,SKYBLUE);
 DrawRectangle(0,h-64,w,64,Fade(BLACK,.8f));DrawText(TextFormat("%s  %s",ct_hint(i,CT_ACTIVATE_STRUCTURE),ce.physicsPreview?"Reset physics":"Activate structure physics"),8,h-60,font,SKYBLUE);DrawText(TextFormat("%s apply | %s erase | %s palette | %s properties",ct_hint(i,CT_PLACE),ct_hint(i,CT_ERASE),ct_hint(i,CT_PALETTE),ct_hint(i,CT_PROPERTIES)),8,h-40,font,WHITE);DrawText(TextFormat("%s",ce.notice),8,h-21,font-2,GOLD);
 if(s->panel){int rows=s->panel==1?8:11;DrawRectangle(6,65,w-12,rows*(font+5)+36,Fade((Color){15,24,32,255},.97f));for(int n=0;n<rows;n++){const char *label;if(s->panel==1)label=ceTools[n];else {switch(n){case 0:label=ce_item_name(s);break;case 1:case 2:case 3:label=TextFormat("Brush %s: %d",n==1?"width":n==2?"height":"depth",s->brush[n-1]);break;case 4:label=TextFormat("Distance: %.0f m",s->distance);break;case 5:label=TextFormat("Color: %d",s->color+1);break;case 6:label=s->material?"Material: fluid":"Material: solid";break;case 7:label=s->reinforced?"Reinforced: yes":"Reinforced: no";break;case 8:label=TextFormat("Rotation: %.0f",s->yaw);break;case 9:label=TextFormat("Entrance group: %d",s->group);break;default:label=TextFormat("Linked trap: %d",s->trap);}}
 DrawText(label,16,75+n*(font+5),font,n==s->row?GOLD:WHITE);}DrawText("Up/down select | Left/right change | B / Tab close",16,80+rows*(font+5),font-2,SKYBLUE);if(s->panel==2&&ceThumb[i].id&&s->tool>=1&&s->tool<=3)DrawTexturePro(ceThumb[i].texture,(Rectangle){0,0,240,-140},(Rectangle){w-220,90,210,123},(Vector2){0},0,WHITE);}
}
static void ce_overlay(void){if(ce.testing)DrawText(ce.sandbox?"SANDBOX | F5 / Start: tools | F10 / Back+Start: editor":"PLAYTEST | F10 / Back+Start: return to editor",16,GetScreenHeight()-55,18,GOLD);}
static void ce_shutdown(void){for(int i=0;i<4;i++){if(ceThumb[i].id)UnloadRenderTexture(ceThumb[i]);ceThumb[i]=(RenderTexture2D){0};ceThumbKind[i]=-99;}ce_clear_history();cm_free(cmDocument);cmDocument=NULL;cm_free(cmRuntime);cmRuntime=NULL;if(ce.files.paths)UnloadDirectoryFiles(ce.files);ce.files=(FilePathList){0};}
static void ce_thumbnails(void){
 if(!ce.active||ce.testing||ce.menu)return;
 for(int i=0;i<activePlayers&&i<4;i++){CeSeat *s=&ceSeats[i];if(s->panel!=2||s->tool<1||s->tool>3)continue;int type=ce_pickup(s);if(ceThumbKind[i]==type&&ceThumb[i].id)continue;if(!ceThumb[i].id)ceThumb[i]=LoadRenderTexture(240,140);ceThumbKind[i]=type;
 BeginTextureMode(ceThumb[i]);ClearBackground((Color){15,24,32,255});Camera3D cam={.position={1.5f,.9f,1.8f},.target={0,.1f,0},.up={0,1,0},.fovy=35,.projection=CAMERA_PERSPECTIVE};BeginMode3D(cam);int kind=eq_pickup_kind(type);weapon_renderer_init();WeaponVisual vis={.shot_time=-100,.charge=.8f,.activity=1};if(kind>=0){WeaponPose p={.kind=kind,.transform=MatrixIdentity()};weapon_draw(p,&vis,cam.position,0,false,WEAPON_COLOR);weapon_draw(p,&vis,cam.position,0,false,WEAPON_GLASS);}else if(type==PICKUP_GOLD_TETHER)weapon_draw_pickup(false,(Vector3){0},cam.position,0,WEAPON_COLOR);else DrawCube((Vector3){0},.4f,.4f,.4f,type==PICKUP_HEALTH?GREEN:type==PICKUP_AMMO?SKYBLUE:GOLD);EndMode3D();EndTextureMode();}
}
static bool ce_frame(void){
 ce_thumbnails();
 if(ce.active&&ce.menu&&!ce.testing){for(int d=0;d<4;d++)if(IsGamepadAvailable(d)&&IsGamepadButtonPressed(d,GAMEPAD_BUTTON_MIDDLE_RIGHT)){bool used=false;for(int i=0;i<activePlayers;i++)used|=playerInput[i]==INPUT_TYPE_GAMEPAD&&seat_device(i)==d;if(!used&&activePlayers<4){int i=activePlayers++;ce.seats=multiplayerActivePlayers=activePlayers;localSeats[i].device=d;playerInput[i]=multiplayerPlayerInput[i]=INPUT_TYPE_GAMEPAD;players[i].pos=(Vector3){0,5,0};ce_notice(TextFormat("P%d joined the editor",i+1));}}}
 if(ce.testing){bool back=ct_pressed(0,CT_TEST_RETURN);for(int i=0;i<ce.seats;i++)back|=ce_pad_down(i,GAMEPAD_BUTTON_MIDDLE_LEFT)&&ce_pad_press(i,GAMEPAD_BUTTON_MIDDLE_RIGHT);if(back){ce_return();return false;}if(ce.sandbox&&(ct_pressed(0,CT_TEST_TOOLS)||ce_pad_press(0,GAMEPAD_BUTTON_MIDDLE_RIGHT))){ce.menu=7;ce.row=0;}if(!ce.menu)return false;}
 if(!ce.menu)return false;bool confirm=couch_confirm(),back=couch_back();int nav=couch_axis(false),change=couch_axis(true);if(mouseCaptured){EnableCursor();mouseCaptured=false;mouseButtonsArmed=false;}
 BeginDrawing();ClearBackground((Color){15,24,32,255});int w=GetScreenWidth(),h=GetScreenHeight();int font=h<650?18:23;DrawText(ce.lobby?"MAP LIBRARY":"CREATIVE ARENA EDITOR",32,24,30,RAYWHITE);DrawText(ce.notice,32,h-40,16,GOLD);
 if(ce.typing){const char *alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 _-./:";int k;while((k=GetCharPressed()))if(k>=32&&k<127&&strlen(ce.text)<CM_PATH-2){size_t n=strlen(ce.text);ce.text[n]=k;ce.text[n+1]=0;}if(IsKeyPressed(KEY_BACKSPACE)&&ce.text[0])ce.text[strlen(ce.text)-1]=0;ce.textKey=(ce.textKey+change+nav*10+(int)strlen(alphabet)*10)%(int)strlen(alphabet);if(couch_menu_device()>=0&&confirm&&strlen(ce.text)<CM_PATH-2){int n=strlen(ce.text);ce.text[n]=alphabet[ce.textKey];ce.text[n+1]=0;}if(couch_menu_button(GAMEPAD_BUTTON_RIGHT_FACE_LEFT)&&ce.text[0])ce.text[strlen(ce.text)-1]=0;
 DrawText(ce.textAction==1?"Map name":"File path",32,90,font,WHITE);DrawText(ce.text,32,125,font,SKYBLUE);DrawText(TextFormat("Character: %c | A add | X erase | Start accept | B cancel",alphabet[ce.textKey]),32,170,18,WHITE);DrawText("Keyboard: type, Enter accepts",32,205,18,WHITE);
 bool done=IsKeyPressed(KEY_ENTER)||couch_menu_button(GAMEPAD_BUTTON_MIDDLE_RIGHT);if(back)ce.typing=false;else if(done&&ce.text[0]){if(ce.textAction==1){bool valid=strlen(ce.text)<CM_NAME;for(const char *p=ce.text;*p;p++)if(!isalnum((unsigned char)*p)&&*p!=' '&&*p!='_'&&*p!='-')valid=false;if(!valid)ce_notice("Use letters, numbers, spaces, - or _ (under 64 characters)");else{char path[CM_PATH];snprintf(path,sizeof(path),"maps/%s.map",ce.text);MakeDirectory("maps");if(FileExists(path)&&strcmp(path,ce.path)){ce_notice("Name exists; open it first or choose another name");}else {char oldName[CM_NAME];snprintf(oldName,CM_NAME,"%s",cmDocument->name);snprintf(cmDocument->name,CM_NAME,"%s",ce.text);if(ce_save(path)){for(int n=0;n<ce.historyCount;n++)snprintf(ce.history[n]->name,CM_NAME,"%s",cmDocument->name);ce.typing=false;if(ce.menu==6){int action=ce.pending;ce_action(action);}}else snprintf(cmDocument->name,CM_NAME,"%s",oldName);}}}
 else if(ce.textAction==2){CmMap *m=cm_read(ce.text);if(m){if(ce.lobby){snprintf(cmSelectedPath,CM_PATH,"%s",ce.text);cm_free(cmRuntime);cmRuntime=m;useCustomMap=true;ce.lobby=false;ce.menu=0;}else{ce_accept(m,NULL);ce.dirty=true;}ce.typing=false;}else ce_notice(cmError);}else if(ce_save(ce.text))ce.typing=false;}EndDrawing();return true;}
 const char *rows[300];int actions[300]={0},count=0;
 if(ce.menu==1){rows[count]="New blank map";actions[count++]=10;rows[count]="Open map";actions[count++]=20;rows[count]="Edit Foundry copy";actions[count++]=11;rows[count]="Edit Crosscut copy";actions[count++]=12;rows[count]="Continue editing";actions[count++]=0;rows[count]="Exit";actions[count++]=30;}
 else if(ce.menu==2){const char *r[]={"Resume","Save","Save As","Open","Import path","Export path","World settings","Test sandbox","Test versus","Test Firefight","New / arena templates","Controls & Settings","Exit"};int a[]={0,40,41,20,42,43,60,50,51,52,61,70,30};for(int n=0;n<13;n++){rows[count]=r[n];actions[count++]=a[n];}}
 else if(ce.menu==3){if(ce.lobby){rows[count]="Use built-in arena";actions[count++]=90;}for(unsigned n=0;n<ce.files.count&&count<280;n++){rows[count]=ceFileSummary[n];actions[count++]=100+n;}rows[count]="Import external path";actions[count++]=42;rows[count]="Back";actions[count++]=0;}
 else if(ce.menu==4){rows[count++]=TextFormat("Width: %.0f metres (40–256)",cmDocument->width);rows[count++]=cmDocument->periodic?"Boundaries: periodic":"Boundaries: bounded";rows[count++]="Back";}
 else if(ce.menu==7){rows[count++]=TextFormat("Enemy: %s",ce.sandboxEnemy==ENEMY_TYPE_COUNT?"Golem":enemy_names[ce.sandboxEnemy]);rows[count++]="Spawn selected enemy ahead";rows[count++]=TextFormat("Weapon: %s",eqGuns[ce.sandboxGun].name);rows[count++]="Equip selected weapon";rows[count++]="Restore Matter and grenades";rows[count++]="Resume sandbox";rows[count++]="Return to editor";}
 else if(ce.menu==6){rows[count++]="Save changes";rows[count++]="Discard changes and continue";rows[count++]="Cancel";}
 ce.row=(ce.row+nav+count)%count;int visible=(h-150)/(font+13),start=ce.row>=visible?ce.row-visible+1:0;
 for(int n=start;n<count&&n<start+visible;n++){Rectangle r={24,85+(n-start)*(font+13),w-48,font+9};if(CheckCollisionPointRec(GetMousePosition(),r)&&IsMouseButtonPressed(MOUSE_BUTTON_LEFT)){ce.row=n;confirm=true;}if(n==ce.row)DrawRectangleRec(r,(Color){35,65,80,255});DrawText(rows[n],34,r.y+3,font,n==ce.row?GOLD:RAYWHITE);}
 EndDrawing();
 if(back){if(ce.menu==7){ce.menu=0;return true;}if(ce.menu==1)ce_action(30);else if(ce.lobby){ce.menu=0;ce.lobby=false;}else ce.menu=2;return true;}
 if(ce.menu==4){if(ce.row<2&&(change||confirm)){CmMap *m=cm_clone(cmDocument);if(ce.row==0)m->width=clampf(m->width+(change?change:1)*8,40,256);else m->periodic=!m->periodic;if(cm_validate(m)){cm_free(cmDocument);cmDocument=m;ce_refresh();ce_commit();}else{for(int n=0;n<m->count;n++){CmCell *c=&m->cells[n];Vector3 p={(c->x+.5f)*VOXEL_SIZE,(c->y+.5f)*VOXEL_SIZE,(c->z+.5f)*VOXEL_SIZE};if(!cm_inside(m,p)){ceSeats[0].a=ceSeats[0].b=p;ceSeats[0].selected=true;break;}}cm_free(m);ce_notice(cmError);}}else if(confirm){ce.menu=2;ce.row=0;}return true;}
 if(ce.menu==7){if(ce.row==0&&(change||confirm))ce.sandboxEnemy=(ce.sandboxEnemy+(change?change:1)+ENEMY_TYPE_COUNT+1)%(ENEMY_TYPE_COUNT+1);else if(ce.row==2&&(change||confirm))ce.sandboxGun=(ce.sandboxGun+(change?change:1)+EQ_GUN_COUNT)%EQ_GUN_COUNT;else if(confirm){if(ce.row==1){Vector3 at=world_wrap(v_add(players[0].pos,v_mul(player_forward(&players[0]),10)));at.y=BASE_EYE_HEIGHT;if(ce.sandboxEnemy==ENEMY_TYPE_COUNT){at.y=.04f;if(golem_physics_start_at(at,180))go_combat_start();}else for(int n=MAX_PLAYERS;n<MAX_COMBATANTS;n++)if(players[n].respawn_timer>0){enemy_lab_spawn(n,ce.sandboxEnemy,at);break;}}else if(ce.row==3){equipment[0].gun[equipment[0].selected]=ce.sandboxGun;equipment[0].switchTime=0;}else if(ce.row==4){players[0].matter=players[0].matterMax;players[0].isExposed=false;grenadeInventory[0]=2;exGrenades[0][0]=exGrenades[0][1]=2;}else if(ce.row==5)ce.menu=0;else if(ce.row==6)ce_return();}return true;}
 if(!confirm)return true;
 if(ce.menu==6){if(ce.row==0){if(ce.path[0]&&ce_save(ce.path)){int a=ce.pending;ce_action(a);}else if(!ce.path[0])ce_text(1,cmDocument->name);return true;}if(ce.row==1){int a=ce.pending;ce.dirty=false;ce_action(a);}else ce.menu=2;return true;}
 int action=actions[ce.row];if(action==90){useCustomMap=false;cmSelectedPath[0]=0;firefightFoundry=true;if(gameMode==GAME_MODE_DEATHMATCH)currentWorldType=WORLD_TYPE_CROSSCUT;ce.menu=0;ce.lobby=false;return true;}if(action>=100){CmMap *m=cm_read(ce.files.paths[action-100]);if(!m){ce_notice(cmError);return true;}if(ce.lobby){snprintf(cmSelectedPath,CM_PATH,"%s",ce.files.paths[action-100]);cm_free(cmRuntime);cmRuntime=m;useCustomMap=true;ce.lobby=false;ce.menu=0;}else ce_accept(m,ce.files.paths[action-100]);}
 else if(action==60){ce.menu=4;ce.row=0;}else if(action==61){ce.menu=1;ce.row=0;}else if(!action){ce.menu=0;ce.lobby=false;}else ce_action(action);return true;
}
