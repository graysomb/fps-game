/* Cover and ballistic decisions are authoritative; grenades use the shared pool,
   bounce fuse, damage, avoidance, visuals and network replication. */
#define SAPPER_THROW_RANGE 28.0f
#define SAPPER_SAFE_RANGE 7.0f
#define SAPPER_THROW_COOLDOWN 6.0f
static Vector3 sapper_origin(Vector3 eye){return v_add(eye,(Vector3){0,.35f,0});}
static bool sapper_arc(Vector3 eye,Vector3 target,Vector3 *velocity){
 Vector3 start=sapper_origin(eye),end=world_image(target,start);end.y-=BASE_EYE_HEIGHT-.2f;
 Vector3 delta=v_sub(end,start);float distance=sqrtf(delta.x*delta.x+delta.z*delta.z);
 if(distance<SAPPER_SAFE_RANGE || distance>SAPPER_THROW_RANGE)return false;
 const Vector3 offsets[]={{0,0,0},{.17f,0,0},{-.17f,0,0},{0,.17f,0},{0,-.17f,0},{0,0,.17f},{0,0,-.17f}};
 for(int attempt=0;attempt<5;attempt++){
  float flight=1.2f+attempt*.3f;
  Vector3 v=v_mul(delta,1/flight);v.y+=.5f*GRAVITY*(flight+1.f/120);
  if(v_length(v)>24)continue;
  bool clear=true;Vector3 previous=start;
  for(int step=1;step<=48 && clear;step++){
   float t=flight*step/48;Vector3 at=v_add(start,v_mul(v,t));at.y-=.5f*GRAVITY*t*t;
   Vector3 segment=v_sub(at,previous);float length=v_length(segment);
   if(at.y<.15f){clear=false;break;}
   for(int n=0;n<7;n++)if(first_voxel_hit((Ray){v_add(previous,offsets[n]),v_mul(segment,1/length)},length,-1)>=0){clear=false;break;}
   previous=at;
  }
  if(clear){*velocity=v;return true;}
 }
 return false;
}
static bool sapper_hidden(int i,Vector3 eye,Vector3 threat,uint64_t *cover){
 Vector3 d=world_delta(eye,threat);float length=v_length(d);VoxelHit hit;
 if(length<.1f || !first_voxel_hit_detailed((Ray){threat,v_mul(d,1/length)},length,-1,false,&hit) || voxels[hit.id].simulate)return false;
 if(!is_view_occluded_by_voxels(threat,v_add(eye,(Vector3){0,.2f,0})))return false;
 if(cover)*cover=voxels[hit.id].identity;return true;
}
static bool sapper_ground(int i,Vector3 *eye){
 // Sample the local floor, including platforms, without choosing unreachable roofs.
 Vector3 origin=v_add(*eye,(Vector3){0,1.5f,0});VoxelHit hit;
 float ground=0;
 if(first_voxel_hit_detailed((Ray){origin,{0,-1,0}},5,-1,false,&hit)){
  if(voxels[hit.id].simulate)return false;ground=origin.y-hit.t;
 }
 float height=ground+BASE_EYE_HEIGHT+.02f;
 if(fabsf(height-players[i].pos.y)>1)return false;
 eye->y=height;
 return (worldPeriodic || (fabsf(eye->x)<FLOOR_SIZE-1 && fabsf(eye->z)<FLOOR_SIZE-1)) && !player_aabb_hits_world_ex(*eye,player_collision_half(&players[i]),i);
}
static bool sapper_find_cover(int i){
 EnemySpecialState *s=&enemySpecial[i];Vector3 threat=aiCombat[i].known,origin=players[i].pos,best={0};
 float scores[4]={FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};Vector3 candidates[4];uint64_t covers[4]={0};
 for(int ring=0;ring<=4;ring++)for(int n=0;n<(ring?16:1);n++){
  float angle=n*2*PI/16;Vector3 at=world_wrap(v_add(origin,(Vector3){cosf(angle)*ring*2,0,sinf(angle)*ring*2}));
  Vector3 delta=world_delta(at,threat);float range=v_length(delta);
  if(range<SAPPER_SAFE_RANGE || range>SAPPER_THROW_RANGE || !sapper_ground(i,&at))continue;
  uint64_t cover;
  if(!sapper_hidden(i,at,threat,&cover))continue;
  float cost=v_length(world_delta(at,origin))+.15f*fabsf(range-14);
  for(int k=0;k<4;k++)if(cost<scores[k]){
   for(int j=3;j>k;j--){scores[j]=scores[j-1];candidates[j]=candidates[j-1];covers[j]=covers[j-1];}
   scores[k]=cost;candidates[k]=at;covers[k]=cover;break;
  }
 }
 uint64_t support=0;
 // Only the four best cover positions need the more expensive swept arc test.
 for(int k=0;k<4;k++)if(covers[k]){Vector3 velocity;if(sapper_arc(candidates[k],threat,&velocity)){best=candidates[k];support=covers[k];break;}}
 s->selected=support;if(support)s->buildOrigin=best;return support!=0;
}
static bool sapper_throw(int i,Vector3 velocity){
 if(netTransport.role==NET_ROLE_CLIENT || players[i].respawn_timer>0)return false;
 int free=-1,owned=0;
 for(int n=0;n<CLARITY_GRENADES;n++){if(!grenades[n].id)free=n;else if(grenades[n].owner==i)owned++;}
 if(free<0 || owned>=2)return false;
 // Carriers use their finite specialty supply, then fall back to the Sapper's
 // innate regular bombs. This does not refill the player's collectible supply.
 int kind=exSelectedGrenade[i];bool supply=ex_grenade_count(i,kind)>0;
 if(!supply)kind=GRENADE_REGULAR;
 grenades[free]=(GrenadeState){.id=clarityNextId++,.owner=i,.pos=sapper_origin(players[i].pos),.vel=velocity,.fuse=-1,
  .kind=kind,.incarnation=exIncarnation[i],.launchOrigin=players[i].pos,.adhesive=kind==GRENADE_REGULAR && eq_upgrade(i,EQ_ADHESIVE)};
 if(supply){if(kind)exGrenades[i][kind-1]--;else grenadeInventory[i]--;if(!ex_grenade_count(i,kind))ex_cycle(i);}
 grenadeCooldown[i]=.4f;play_sfx_at(SFX_MELEE,players[i].pos,i,2,1.4f);return true;
}
static void sapper_update(int i,float dt){
 if(netTransport.role==NET_ROLE_CLIENT || (enemyLab && enemyLabPaused && !enemyLabStep))return;
 Player *p=&players[i];EnemySpecialState *s=&enemySpecial[i];
 ai_tick(i,dt);s->target=aiCombat[i].target;s->cooldown=fmaxf(0,s->cooldown-dt);s->searchTimer-=dt;s->melee-=dt;
 if(clarity_evade(i,dt)){s->windup=0;s->state=0;s->reason="avoiding blast";return;}
 if(s->target<0){
  s->selected=0;s->windup=0;s->state=0;s->reason="searching for players";
  BotState *b=&botStates[i];b->pathTimer-=dt;
  if(b->pathTimer<=0){float angle=ai_random(i)*2*PI;b->moveTarget=world_wrap(v_add(p->pos,(Vector3){cosf(angle)*8,0,sinf(angle)*8}));b->pathTimer=3;}
  ai_aim(i,b->moveTarget,dt);bot_apply_move(p,b,world_delta(b->moveTarget,p->pos),MOVE_SPEED*.6f,dt);return;
 }
 Vector3 to=world_delta(aiCombat[i].known,p->pos);float distance=v_length(to);ai_aim(i,aiCombat[i].known,dt);
 if(s->selected){int v=enemy_voxel(s->selected);if(v<0 || voxels[v].simulate || !sapper_hidden(i,s->buildOrigin,aiCombat[i].known,NULL))s->selected=0;}
 if(s->searchTimer<=0){sapper_find_cover(i);s->searchTimer=1+(i%4)*.13f;}
 Vector3 move={0};
 if(s->selected){move=world_delta(s->buildOrigin,p->pos);if(v_length(move)<.45f)move=(Vector3){0};s->reason="moving to cover";}
 else {if(distance<10)move=v_mul(to,-1);else if(distance>22)move=to;else move=(Vector3){to.z,0,-to.x};s->reason="seeking cover";}
 if(p->isExposed){s->windup=0;s->state=0;bot_apply_move(p,&botStates[i],move,MOVE_SPEED*.75f,dt);return;}
 Vector3 velocity;bool arrived=s->selected && v_length(world_delta(s->buildOrigin,p->pos))<.7f;
 bool canThrow=distance>=SAPPER_SAFE_RANGE && distance<=SAPPER_THROW_RANGE && (arrived || !s->selected) && aiCombat[i].memory>0;
 if(s->windup>0){
  p->vel.x=p->vel.z=0;s->windup-=dt;s->reason="lob wind-up";s->work=2;
  if(!canThrow){s->windup=0;s->state=0;s->cooldown=1;return;}
  if(s->windup<=0){
   bool fired=sapper_arc(p->pos,s->aim,&velocity) && sapper_throw(i,velocity);
   s->actions+=fired;s->cooldown=fired?SAPPER_THROW_COOLDOWN:1;s->state=0;s->work=3;
   s->reason=fired?"reloading behind cover":"arc blocked; repositioning";
   if(!fired){s->selected=0;s->searchTimer=1;}
  }return;
 }
 if(canThrow && s->cooldown<=0 && (!aiCombat[i].visible || aiCombat[i].reaction<=0) && sapper_arc(p->pos,aiCombat[i].known,&velocity)){
  s->aim=aiCombat[i].known;s->state=6;s->windup=.8f;s->work=2;enemy_cue(i);p->vel.x=p->vel.z=0;return;
 }
 if(arrived){s->reason="holding cover";p->vel.x=p->vel.z=0;s->work=0;}
 else bot_apply_move(p,&botStates[i],move,MOVE_SPEED*.75f,dt);
}
