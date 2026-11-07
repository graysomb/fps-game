static const char *cmPickupNames[]={"matter","health","void","voxel_launcher","upgraded_tether","grenade","shotgun","cutter","railgun","rocket","sword","adhesive","reinforced_build","glob","glock","dmr","black_hole","igloo","shield_regenerator"};
static void cm_fail(const char *s){snprintf(cmError,sizeof(cmError),"%s",s);}
static CmMap *cm_new(void){CmMap *m=calloc(1,sizeof(*m));if(m){m->width=80;m->periodic=true;strcpy(m->name,"Untitled");}return m;}
static void cm_free(CmMap *m){if(m){free(m->cells);free(m);}}
static bool cm_reserve(CmMap *m,int n){if(n<0||n>MAX_VOXELS)return false;if(n<=m->capacity)return true;int cap=clampi(n+1024,1,MAX_VOXELS);CmCell *p=realloc(m->cells,(size_t)cap*sizeof(*p));if(!p)return false;m->cells=p;m->capacity=cap;return true;}
static CmMap *cm_clone(const CmMap *m){if(!m)return NULL;CmMap *p=malloc(sizeof(*p));if(!p)return NULL;*p=*m;p->cells=NULL;p->capacity=0;if(!cm_reserve(p,m->count)){free(p);return NULL;}if(m->count)memcpy(p->cells,m->cells,m->count*sizeof(CmCell));return p;}
static int cm_compare(const void *a,const void *b){const CmCell *x=a,*y=b;if(x->x!=y->x)return x->x<y->x?-1:1;if(x->y!=y->y)return x->y<y->y?-1:1;return (x->z>y->z)-(x->z<y->z);}
static bool cm_inside(const CmMap *m,Vector3 p){float r=m->width*.5f;return v_isfinite(p)&&p.x>=-r&&p.x<r&&p.z>=-r&&p.z<r&&p.y>=0&&p.y<=m->width;}
static bool cm_validate(CmMap *m){
 if(!m||!isfinite(m->width)||m->width<40||m->width>256||m->count<0||m->count>MAX_VOXELS||m->pickupCount<0||m->pickupCount>MAX_PICKUPS||m->spawnCount<0||m->spawnCount>CM_SPAWNS||m->duckCount<0||m->duckCount>MAP_DUCKS||m->trapCount<0||m->trapCount>FOUNDRY_TRAPS){cm_fail("Invalid map dimensions or capacity");return false;}
 qsort(m->cells,m->count,sizeof(CmCell),cm_compare);int members[FOUNDRY_TRAPS]={0},triggers[FOUNDRY_TRAPS]={0},bodies[FOUNDRY_TRAPS]={0};
 for(int n=0;n<m->count;n++){CmCell *c=&m->cells[n];if(c->x< -512||c->x>512||c->y<0||c->y>512||c->z< -512||c->z>512||!cm_inside(m,(Vector3){(c->x+.5f)*VOXEL_SIZE,(c->y+.5f)*VOXEL_SIZE,(c->z+.5f)*VOXEL_SIZE})||c->material>1||c->reinforced>1||c->trap>m->trapCount||c->flags>7||(!c->trap&&c->flags)||(n&&cm_compare(c,c-1)==0)){cm_fail("Invalid, duplicate or out-of-bounds voxel");return false;}if(c->trap){int t=c->trap-1;members[t]++;triggers[t]+=(c->flags&CM_TRIGGER)!=0;bodies[t]+=(c->flags&CM_BODY)!=0;if(c->material||!c->flags||(c->flags&CM_BODY&&c->flags&CM_REMOVE)){cm_fail("Invalid trap membership");return false;}}}
 for(int t=0;t<m->trapCount;t++)if(members[t]>FOUNDRY_TRAP_PAYLOAD||!triggers[t]||!bodies[t]){cm_fail("Trap needs trigger and body cells (maximum 2048 cells)");return false;}
 for(int n=0;n<m->pickupCount;n++)if(!cm_inside(m,m->supplies[n].pos)||m->supplies[n].type<0||m->supplies[n].type>=PICKUP_TYPE_COUNT){cm_fail("Invalid pickup");return false;}
 for(int n=0;n<m->spawnCount;n++)if(!cm_inside(m,m->spawns[n].pos)||!isfinite(m->spawns[n].yaw)||m->spawns[n].kind<0||m->spawns[n].kind>CM_BOSS||m->spawns[n].group<0||m->spawns[n].group>31){cm_fail("Invalid spawn marker");return false;}
 for(int n=0;n<m->duckCount;n++)if(!cm_inside(m,m->ducks[n])||(!m->periodic&&(fabsf(m->ducks[n].x)+DUCK_ROUTE_RADIUS>m->width*.5f||fabsf(m->ducks[n].z)+DUCK_ROUTE_RADIUS>m->width*.5f))){cm_fail("Duck home outside world");return false;}
 return true;
}
static CmMap *cm_capture(const char *name){
 CmMap *m=cm_new();if(!m)return NULL;snprintf(m->name,CM_NAME,"%s",name);m->width=FLOOR_SIZE*2;m->periodic=worldPeriodic;
 if(!cm_reserve(m,voxel_count)){cm_free(m);return NULL;}
 for(int n=0;n<voxel_count;n++){Voxel *v=&voxels[n];if(v->simulate||v->isBullet)continue;CmCell c={v->gx,v->gy,v->gz,v->color,v->type,v->reinforced,0,0};
  for(int t=0;t<foundryTrapCount;t++){FoundryTrap *f=&foundryTraps[t];bool member=false;for(int k=0;k<f->originalCount;k++)if(f->original[k].gx==c.x&&f->original[k].gy==c.y&&f->original[k].gz==c.z){member=true;break;}if(!member)continue;c.trap=t+1;bool trigger=foundry_trap_trigger(f,c.x,c.y,c.z);bool removed=trigger||(f->fourPosts&&c.y>=f->ry&&c.y<f->y0);c.flags=removed?CM_REMOVE:CM_BODY;if(trigger)c.flags|=CM_TRIGGER;break;}
  m->cells[m->count++]=c;
 }
 m->trapCount=foundryTrapCount;for(int t=0;t<m->trapCount;t++)m->seeds[t]=t+1;
 for(int n=0;n<MAX_PICKUPS;n++)if(pickups[n].active||pickups[n].respawnTimer>0){Pickup p=pickups[n];p.active=true;p.respawnTimer=0;m->supplies[m->pickupCount++]=p;}
 qsort(m->cells,m->count,sizeof(CmCell),cm_compare);return m;
}
static bool cm_apply(const CmMap *m,bool editing){
 if(!m)return false;
 Voxel *originals[FOUNDRY_TRAPS]={0};for(int n=0;n<m->trapCount;n++){originals[n]=calloc(FOUNDRY_TRAP_PAYLOAD,sizeof(Voxel));if(!originals[n]){for(int k=0;k<n;k++)free(originals[k]);cm_fail("Cannot allocate trap state; map preserved");return false;}}
 clear_world_voxels();clear_pickups();worldFloorHalfSize=m->width*.5f;worldPeriodic=m->periodic;
 for(int n=0;n<m->count;n++){const CmCell *c=&m->cells[n];int v=add_static_voxel_at_grid(c->x,c->y,c->z,c->color,c->material);if(v<0){cm_fail("Insufficient voxel capacity");return false;}voxels[v].reinforced=c->reinforced;}
 for(int n=0;n<m->pickupCount;n++){Pickup p=m->supplies[n];p.active=true;p.respawnTimer=0;p.bobTimer=n*.37f;if(!editing&&gameMode==GAME_MODE_FIREFIGHT&&p.type!=PICKUP_HEALTH&&p.type!=PICKUP_AMMO&&p.type!=PICKUP_GRENADE)p.active=false;pickups[n]=p;}
 for(int n=0;n<m->trapCount;n++){FoundryTrap *t=&foundryTraps[n];*t=(FoundryTrap){.kind=FOUNDRY_BAY,.owner=-1,.explicitCells=true,.x0=INT_MAX,.y0=INT_MAX,.z0=INT_MAX,.x1=INT_MIN,.y1=INT_MIN,.z1=INT_MIN};t->original=originals[n];foundryTrapCount=n+1;
  for(int k=0;k<m->count;k++){const CmCell *c=&m->cells[k];if(c->trap!=n+1)continue;int id=table_get_static_only(c->x,c->y,c->z);if(id<0)continue;int j=t->originalCount++;t->original[j]=voxels[id];t->cellFlags[j]=c->flags;t->x0=fminf(t->x0,c->x);t->x1=fmaxf(t->x1,c->x);t->y0=fminf(t->y0,c->y);t->y1=fmaxf(t->y1,c->y);t->z0=fminf(t->z0,c->z);t->z1=fmaxf(t->z1,c->z);}
 }
 rebuild_all_voxel_surfaces();init_static_hash();rebuild_voxel_hash();meshDirty=true;ducks_reset();return true;
}
static bool cm_spawn(int kind,int ordinal,Vector3 *pos,float *yaw){if(!useCustomMap||!cmRuntime)return false;int count=0;for(int n=0;n<cmRuntime->spawnCount;n++)count+=cmRuntime->spawns[n].kind==kind;if(!count)return false;int which=(ordinal%count+count)%count;for(int n=0;n<cmRuntime->spawnCount;n++){CmSpawn *s=&cmRuntime->spawns[n];if(s->kind==kind&&!which--){*pos=s->pos;if(yaw)*yaw=s->yaw;return true;}}return false;}
static bool cm_write(const CmMap *m,const char *path){
 if(!m||!path||strlen(path)>CM_PATH-8){cm_fail("Invalid path");return false;}char temp[CM_PATH];snprintf(temp,sizeof(temp),"%s.tmp",path);FILE *f=fopen(temp,"wb");if(!f){cm_fail("Cannot create file; check folder and permissions");return false;}
 fprintf(f,"FPSMAP2\nNAME %s\nWORLD %.0f %d\nVOXELS %d\n",m->name,m->width,m->periodic,m->count);
 for(int n=0;n<m->count;n++){const CmCell *c=&m->cells[n];fprintf(f,"%d %d %d %u %u %u %u %u %u %u %u\n",c->x,c->y,c->z,c->color.r,c->color.g,c->color.b,c->color.a,c->material,c->reinforced,c->trap,c->flags);}
 fprintf(f,"PICKUPS %d\n",m->pickupCount);for(int n=0;n<m->pickupCount;n++){const Pickup *p=&m->supplies[n];fprintf(f,"%s %.9g %.9g %.9g\n",cmPickupNames[p->type],p->pos.x,p->pos.y,p->pos.z);}
 fprintf(f,"SPAWNS %d\n",m->spawnCount);for(int n=0;n<m->spawnCount;n++){const CmSpawn *s=&m->spawns[n];fprintf(f,"%d %d %.9g %.9g %.9g %.9g\n",s->kind,s->group,s->pos.x,s->pos.y,s->pos.z,s->yaw);}
 fprintf(f,"DUCKS %d\n",m->duckCount);for(int n=0;n<m->duckCount;n++)fprintf(f,"%.9g %.9g %.9g\n",m->ducks[n].x,m->ducks[n].y,m->ducks[n].z);
 fprintf(f,"TRAPS %d\n",m->trapCount);for(int n=0;n<m->trapCount;n++)fprintf(f,"%u\n",m->seeds[n]);fprintf(f,"END\n");bool ok=!ferror(f);if(fclose(f))ok=false;
#ifdef _WIN32
 if(ok)ok=MoveFileExA(temp,path,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0;
#else
 if(ok)ok=rename(temp,path)==0;
#endif
 if(!ok){remove(temp);cm_fail("Save failed; existing map preserved");}return ok;
}
static bool cm_section(FILE *f,const char *expected,int *count,int max){char key[32];return fscanf(f,"%31s %d",key,count)==2&&!strcmp(key,expected)&&*count>=0&&*count<=max;}
static CmMap *cm_read(const char *path){
 FILE *f=fopen(path,"rb");if(!f){cm_fail("Cannot open map file");return NULL;}CmMap *m=cm_new();if(!m){fclose(f);return NULL;}char key[64];bool legacy=false;
 if(fscanf(f,"%63s",key)!=1)goto bad;legacy=!strcmp(key,"FPSMAP1");if(!legacy&&strcmp(key,"FPSMAP2"))goto bad;
 if(!legacy){int wrap;if(fscanf(f," %63s",key)!=1||strcmp(key,"NAME")||fscanf(f," %63[^\r\n]",m->name)!=1)goto bad;if(fscanf(f," %63s %f %d",key,&m->width,&wrap)!=3||strcmp(key,"WORLD")||(wrap!=0&&wrap!=1))goto bad;m->periodic=wrap;}
 else {m->width=40;m->periodic=false;snprintf(m->name,CM_NAME,"Imported legacy map");}
 int count;if(!cm_section(f,"VOXELS",&count,MAX_VOXELS)||!cm_reserve(m,count))goto bad;m->count=count;
 for(int n=0;n<count;n++){CmCell *c=&m->cells[n];int r,g,b,a,mat,rein=0,trap=0,flags=0;int got=legacy?fscanf(f,"%d %d %d %d %d %d %d %d",&c->x,&c->y,&c->z,&r,&g,&b,&a,&mat):fscanf(f,"%d %d %d %d %d %d %d %d %d %d %d",&c->x,&c->y,&c->z,&r,&g,&b,&a,&mat,&rein,&trap,&flags);if(got!=(legacy?8:11)||r<0||r>255||g<0||g>255||b<0||b>255||a<0||a>255||mat<0||mat>1||rein<0||rein>1||trap<0||trap>32||flags<0||flags>7)goto bad;c->color=(Color){r,g,b,a};c->material=mat;c->reinforced=rein;c->trap=trap;c->flags=flags;}
 if(!cm_section(f,"PICKUPS",&m->pickupCount,MAX_PICKUPS))goto bad;
 for(int n=0;n<m->pickupCount;n++){Pickup *p=&m->supplies[n];int type=-1;if(legacy){int active;float timer;if(fscanf(f,"%f %f %f %d %d %f",&p->pos.x,&p->pos.y,&p->pos.z,&type,&active,&timer)!=6||!isfinite(timer))goto bad;}
 else {if(fscanf(f,"%63s %f %f %f",key,&p->pos.x,&p->pos.y,&p->pos.z)!=4)goto bad;for(int k=0;k<PICKUP_TYPE_COUNT;k++)if(!strcmp(key,cmPickupNames[k]))type=k;}p->type=type;p->active=true;}
 if(!legacy){if(!cm_section(f,"SPAWNS",&m->spawnCount,CM_SPAWNS))goto bad;for(int n=0;n<m->spawnCount;n++){CmSpawn *s=&m->spawns[n];if(fscanf(f,"%d %d %f %f %f %f",&s->kind,&s->group,&s->pos.x,&s->pos.y,&s->pos.z,&s->yaw)!=6)goto bad;}
 if(!cm_section(f,"DUCKS",&m->duckCount,MAP_DUCKS))goto bad;for(int n=0;n<m->duckCount;n++)if(fscanf(f,"%f %f %f",&m->ducks[n].x,&m->ducks[n].y,&m->ducks[n].z)!=3)goto bad;
 if(!cm_section(f,"TRAPS",&m->trapCount,FOUNDRY_TRAPS))goto bad;for(int n=0;n<m->trapCount;n++)if(fscanf(f,"%u",&m->seeds[n])!=1)goto bad;if(fscanf(f,"%63s",key)!=1||strcmp(key,"END"))goto bad;}
 if(fscanf(f,"%63s",key)==1)goto bad;fclose(f);if(!cm_validate(m)){cm_free(m);return NULL;}return m;
 bad:fclose(f);cm_free(m);cm_fail("Malformed or unsupported map; current map preserved");return NULL;
}
static bool cm_load_selected(void){CmMap *m=cmTesting?cm_clone(cmDocument):cm_read(cmSelectedPath);if(!m)return false;if(!cm_apply(m,false)){cm_free(m);return false;}cm_free(cmRuntime);cmRuntime=m;return true;}
