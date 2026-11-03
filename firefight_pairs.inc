/* Host-only rifle teams. No inventory or protocol state is stored here. */
static bool fp_pair_eligible(int i){
 return i>=MAX_PLAYERS && i<activePlayers && players[i].respawn_timer<=0 &&
  players[i].enemyType==ENEMY_TYPE_STANDARD && fp_carrier(i) &&
  (eq_gun(i)==EQ_GLOCK || eq_gun(i)==EQ_GLOB);
}
static void fp_pair_clear(int i){
 if(i<0 || i>=MAX_COMBATANTS)return;
 int j=fpPairs[i].mate-1;
 if(j>=0 && j<MAX_COMBATANTS && fpPairs[j].mate==i+1)fpPairs[j]=(FfRiflePair){0};
 fpPairs[i]=(FfRiflePair){0};
}
static bool fp_pair_valid(int i){
 if(!fp_pair_eligible(i))return false;
 int j=fpPairs[i].mate-1;
 return fp_pair_eligible(j) && fpPairs[j].mate==i+1 && eq_gun(i)!=eq_gun(j) && combatantTeam[i]==combatantTeam[j];
}
static void fp_pair_link(int a,int b){
 fp_pair_clear(a);fp_pair_clear(b);
 int glob=eq_gun(a)==EQ_GLOB?a:b,other=glob==a?b:a;
 int source=ai_hostile(glob,aiCombat[glob].target) && aiCombat[glob].memory>0?glob:other;
 AiCombat *sense=&aiCombat[source];
 FfRiflePair state={.target=ai_hostile(source,sense->target)&&sense->memory>0?sense->target:-1,
  .known=sense->known,.memory=sense->memory};
 fpPairs[a]=fpPairs[b]=state;fpPairs[a].mate=b+1;fpPairs[b].mate=a+1;
 if(towerRifles[a].tag)tower_release(a);if(towerRifles[b].tag)tower_release(b);
}
static bool fp_pair_sees(int i,int target){return ai_hostile(i,target) && check_ai_awareness_of_player(i,target,NULL);}
static void fp_pairs_tick(float dt){
 if(!fp_active() || netTransport.role==NET_ROLE_CLIENT)return;
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(fpPairs[i].mate && !fp_pair_valid(i))fp_pair_clear(i);
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(fp_pair_eligible(i) && !fpPairs[i].mate){
  int best=-1;float distance=15.001f;
  for(int j=MAX_PLAYERS;j<activePlayers;j++)if(j!=i && fp_pair_eligible(j) && !fpPairs[j].mate && eq_gun(i)!=eq_gun(j) && combatantTeam[i]==combatantTeam[j]){
   float d=v_length(world_delta(players[j].pos,players[i].pos));if(d<distance){distance=d;best=j;}
  }
  if(best>=0)fp_pair_link(i,best);
 }
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(fp_pair_valid(i)){
  int j=fpPairs[i].mate-1;if(i>j)continue;
  FfRiflePair s=fpPairs[i];s.memory=fmaxf(0,s.memory-dt);
  if(!ai_hostile(i,s.target))s.target=-1;
  bool seen=s.target>=0 && (fp_pair_sees(i,s.target)||fp_pair_sees(j,s.target));
  // Like unpaired riflemen, a nearby decoy can draw a team off a distant
  // player. Don't retain a visible real player indefinitely.
  if(seen && !target_is_hologram(s.target)){
   float nearest=.7f*fminf(v_length(world_delta(s.known,players[i].pos)),v_length(world_delta(s.known,players[j].pos)));
   for(int h=0;h<MAX_PLAYERS;h++){
    int target=MAX_COMBATANTS+h;
    if(!fp_pair_sees(i,target)&&!fp_pair_sees(j,target))continue;
    Vector3 at=ai_target_body(target)->pos;
    float distance=fminf(v_length(world_delta(at,players[i].pos)),v_length(world_delta(at,players[j].pos)));
    if(distance<nearest){nearest=distance;s.target=target;}
   }
  }
  if(s.target<0 || (!seen && s.memory<=0)){
   s.target=-1;float nearest=FLT_MAX;
   for(int t=0;t<MAX_COMBATANTS+MAX_PLAYERS;t++)if(fp_pair_sees(i,t)||fp_pair_sees(j,t)){
    float d=fminf(v_length(world_delta(ai_target_body(t)->pos,players[i].pos)),v_length(world_delta(ai_target_body(t)->pos,players[j].pos)));
    if(d<nearest){nearest=d;s.target=t;}
   }
   seen=s.target>=0;
  }
  if(seen){s.known=ai_target_body(s.target)->pos;s.memory=4;}
  s.patrolTime=fmaxf(0,s.patrolTime-dt);
  if(s.target<0 && s.patrolTime<=0){
   Vector3 center=v_add(players[i].pos,v_mul(world_delta(players[j].pos,players[i].pos),.5f));
   float angle=ai_random(i)*2*PI;
   Vector3 direction=v_length((Vector3){center.x,0,center.z})>20?v_norm((Vector3){-center.x,0,-center.z}):(Vector3){cosf(angle),0,sinf(angle)};
   s.patrol=world_wrap(v_add(center,v_mul(direction,8)));s.patrolTime=3;
  }
  float separation=v_length(world_delta(players[j].pos,players[i].pos));
  if(separation>6)s.regroup=true;else if(separation<4)s.regroup=false;
  fpPairs[i]=fpPairs[j]=s;fpPairs[i].mate=j+1;fpPairs[j].mate=i+1;
 }
}
static void fp_pair_perception(int i){
 FfRiflePair *pair=&fpPairs[i];AiCombat *a=&aiCombat[i];
 if(pair->target<0){a->target=-1;a->visible=a->wasVisible=false;a->memory=0;return;}
 bool visible=fp_pair_sees(i,pair->target);
 if(a->target!=pair->target)ai_select(i,pair->target,pair->known,visible?AI_ALERT_SEEN:AI_ALERT_HEARD);
 a->known=pair->known;a->memory=pair->memory;a->visible=visible;
 if(visible && !a->wasVisible)a->reaction=fmaxf(a->reaction,ai_reaction(i));
 a->wasVisible=visible;
 if(!visible){ai_cancel_burst(i);a->tracking=0;}
 botStates[i].targetIndex=a->target;botStates[i].hasTarget=true;botStates[i].lastKnownTargetPos=a->known;
}
static void fp_pair_move(int i,float dt){
 FfRiflePair *pair=&fpPairs[i];int j=pair->mate-1;Player *p=&players[i];
 Vector3 other=world_image(players[j].pos,p->pos),center=v_mul(v_add(p->pos,other),.5f),goal=p->pos;
 if(pair->regroup)goal=center;
 else if(pair->target<0){goal=world_image(pair->patrol,p->pos);goal.x+=eq_gun(i)==EQ_GLOB?-1.5f:1.5f;ai_aim(i,goal,dt);}
 else if(pair->target>=0){
  Vector3 target=world_image(pair->known,center),forward=v_sub(target,center);forward.y=0;
  if(v_length(forward)<.01f)forward=(Vector3){0,0,-1};else forward=v_norm(forward);
  Vector3 side={forward.z,0,-forward.x};
  goal=v_add(v_sub(target,v_mul(forward,10)),v_mul(side,eq_gun(i)==EQ_GLOB?-1.5f:1.5f));goal.y=p->pos.y;
 }
 Vector3 move=world_delta(goal,p->pos);move.y=0;
 if(v_length(move)<.5f){p->vel.x=p->vel.z=0;return;}
 bot_apply_move(p,&botStates[i],move,fminf(MOVE_SPEED*(pair->target<0?.45f:.8f),v_length(move)/fmaxf(dt,.001f)),dt);
 // Keep collision-safe steering from cutting across a partner's active firing lane.
 Vector3 next=v_add(p->pos,v_mul(p->vel,fmaxf(dt,1.f/60)));
 Vector3 from=world_delta(next,players[j].pos),facing=player_forward(&players[j]);float along=v_dot(from,facing);
 Vector3 current=world_delta(p->pos,players[j].pos);float oldAlong=v_dot(current,facing);
 float lane=v_length(v_sub(from,v_mul(facing,along))),oldLane=v_length(v_sub(current,v_mul(facing,oldAlong)));
 // An actor already inside a lane must be allowed to step out of it.
 bool crowd=v_length(from)<1.2f && v_length(from)<=v_length(current);
 bool crossing=aiCombat[j].visible && along>0 && along<10 && lane<.8f && !(oldAlong>0 && oldAlong<10 && lane>oldLane);
 if(crowd || crossing)p->vel.x=p->vel.z=0;
}
static bool fp_pair_blocks_fire(int i){
 int j=fpPairs[i].mate-1;Vector3 to=world_delta(players[j].pos,players[i].pos),forward=player_forward(&players[i]);float along=v_dot(to,forward);
 return along>0 && along<v_length(world_delta(aiCombat[i].known,players[i].pos)) && v_length(v_sub(to,v_mul(forward,along)))<player_collision_half(&players[j])+.15f;
}
static bool fp_pair_pending(void){
 if(!fp_active())return false;bool glob=false,glock=false;
 for(int n=0;n<fpPendingCount;n++)if(fpPending[n]>=0){glob|=fpRules[fpPending[n]].type==PICKUP_GLOB;glock|=fpRules[fpPending[n]].type==PICKUP_GLOCK;}
 return glob && glock;
}
static bool fp_spawn_pair(int living,int cap){
 if(!fp_pair_pending() || ffRemaining[ENEMY_TYPE_STANDARD]<2 || living+2>cap)return false;
 int slots[2],count=0;
 for(int i=MAX_PLAYERS;i<MAX_COMBATANTS && count<2;i++)if(players[i].respawn_timer>0 && !go_combat_actor(i))slots[count++]=i;
 if(count<2)return false;
 int entrance=ffEntrance,positions[2]={ffEntrancePosition[0],ffEntrancePosition[1]};
 EnemyType old=players[slots[0]].enemyType;players[slots[0]].enemyType=ENEMY_TYPE_STANDARD;
 Vector3 at[2];bool clear=ff_spawn_position(slots[0],&at[0]);players[slots[0]].enemyType=old;
 bool second=false;
 if(clear)for(int n=0;n<8;n++){
  float angle=n*PI/4;at[1]=world_wrap(v_add(at[0],(Vector3){cosf(angle)*3,0,sinf(angle)*3}));
  if(spawn_position_clear_ex(at[1],PLAYER_SIZE*.5f,4,slots[1])){second=true;break;}
 }
 if(!clear || !second || !fp_reserve(fpCount+2)){
  ffEntrance=entrance;ffEntrancePosition[0]=positions[0];ffEntrancePosition[1]=positions[1];return false;
 }
 // Move one rule of each kind to the front; leave the rest of the queue intact.
 for(int k=0;k<2;k++)for(int n=k;n<fpPendingCount;n++)if(fpPending[n]>=0 && fpRules[fpPending[n]].type==(k?PICKUP_GLOB:PICKUP_GLOCK)){
  int swap=fpPending[k];fpPending[k]=fpPending[n];fpPending[n]=swap;break;
 }
 int level=ff_ai_level(firefightWave,ENEMY_TYPE_STANDARD);
 for(int k=0;k<2;k++){
  int i=slots[k];init_firefight_bot_entity(i,ENEMY_TYPE_STANDARD,level==0?INPUT_TYPE_BOT_EASY:level==1?INPUT_TYPE_BOT_MEDIUM:INPUT_TYPE_BOT_HARD);
  players[i].pos=players[i].death_pos=at[k];clarity_spawn(i);fp_spawn(i,ENEMY_TYPE_STANDARD);
 }
 fp_pair_link(slots[0],slots[1]);ffRemaining[ENEMY_TYPE_STANDARD]-=2;firefightEnemiesSpawned+=2;ffSpawnOrdinal+=2;
 ffEntrance=1-ffEntrance;ffSpawnTimer=2;return true;
}
