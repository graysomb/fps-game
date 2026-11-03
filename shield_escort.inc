/* Support formations are host-side; positions/shots use existing replication.
   shieldEscort stores either a shield bearer or an anchor, plus one. */
static bool shield_rifle_eligible(int i){
 return gameMode==GAME_MODE_FIREFIGHT && i>=MAX_PLAYERS && i<activePlayers &&
  is_player_bot(i) && players[i].respawn_timer<=0 && players[i].enemyType==ENEMY_TYPE_STANDARD && eq_gun(i)!=EQ_BLADE;
}
static bool shield_bearer_eligible(int i){
 return i>=MAX_PLAYERS && i<activePlayers && is_player_bot(i) && players[i].respawn_timer<=0 &&
  players[i].enemyType==ENEMY_TYPE_SHIELDBEARER && enemySpecial[i].equipmentCount>0;
}
static bool anchor_escort_eligible(int i){
 return i>=MAX_PLAYERS && i<activePlayers && is_player_bot(i) && players[i].respawn_timer<=0 &&
  players[i].enemyType==ENEMY_TYPE_ANCHOR;
}
static bool rifle_friendly_shield(int shooter,int v){
 if(!clarity_shield_voxel(v))return false;
 int owner=voxels[v].owner;
 if(shooter==owner)return true; // The bearer's own shield cannot blind it.
 if(!shield_rifle_eligible(shooter))return false;
 return shield_bearer_eligible(owner) && combatantTeam[shooter]==combatantTeam[owner];
}
static bool rifle_view_blocked(int i,Vector3 target){
 for(int pass=0;pass<2;pass++){
  Vector3 point=target;if(pass)point.y-=.4f;
  Vector3 d=world_delta(point,players[i].pos);float length=v_length(d);
  if(length<=.4f)return false;d=v_mul(d,1/length);VoxelHit hit;
  if(!rifle_shield_ray((Ray){v_add(players[i].pos,v_mul(d,.2f)),d},length-.4f,-1,i,&hit))return false;
 }
 return true;
}
static bool shield_escort_valid(int i){
 if(!shield_rifle_eligible(i))return false;
 int b=shieldEscort[i]-1;
 return (shield_bearer_eligible(b) || anchor_escort_eligible(b)) && combatantTeam[i]==combatantTeam[b] && v_length(world_delta(players[i].pos,players[b].pos))<18;
}
static bool shield_has_escort(int b){
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(shieldEscort[i]==b+1 && shield_escort_valid(i))return true;
 return false;
}
static void shield_escort_tick(void){
 if(netTransport.role==NET_ROLE_CLIENT)return;
 int assigned[MAX_COMBATANTS]={0};
 for(int i=MAX_PLAYERS;i<activePlayers;i++){
  if(!shield_escort_valid(i))shieldEscort[i]=0;
  else assigned[shieldEscort[i]-1]++;
 }
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(shield_rifle_eligible(i) && !shieldEscort[i]){
  int best=-1;float nearest=12.001f;
  int pair=fp_pair_valid(i)?fpPairs[i].mate-1:-1;
  for(int b=MAX_PLAYERS;b<activePlayers;b++)if((shield_bearer_eligible(b) || anchor_escort_eligible(b)) && assigned[b]<(anchor_escort_eligible(b)?3:2) && combatantTeam[b]==combatantTeam[i]){
   float distance=v_length(world_delta(players[b].pos,players[i].pos));
   if(distance<=12 && pair>=0 && shieldEscort[pair]==b+1){best=b;break;}
   if(distance<nearest){nearest=distance;best=b;}
  }
  if(best>=0){shieldEscort[i]=best+1;assigned[best]++;if(towerRifles[i].tag)tower_release(i);}
 }
 // Share only actual sightings, never another AI's unverified current target position.
 for(int i=MAX_PLAYERS;i<activePlayers;i++)if(shield_escort_valid(i)){
  int b=shieldEscort[i]-1;
  if(aiCombat[i].visible && ai_hostile(b,aiCombat[i].target))ai_alert(b,aiCombat[i].target,aiCombat[i].known,AI_ALERT_HEARD);
  else if(aiCombat[b].visible && ai_hostile(i,aiCombat[b].target))ai_alert(i,aiCombat[b].target,aiCombat[b].known,AI_ALERT_HEARD);
 }
}
static void shield_escort_move(int i,float dt){
 int b=shieldEscort[i]-1;Player *p=&players[i];
 Vector3 forward=player_forward(&players[b]);forward.y=0;forward=v_norm(forward);
 Vector3 right={-forward.z,0,forward.x};int rank=0;
 for(int j=MAX_PLAYERS;j<i;j++)if(shieldEscort[j]==b+1 && shield_escort_valid(j))rank++;
 bool anchor=anchor_escort_eligible(b);
 if(anchor && aiCombat[b].target>=0){forward=world_delta(aiCombat[b].known,players[b].pos);forward.y=0;forward=v_norm(forward);right=(Vector3){-forward.z,0,forward.x};}
 Vector3 goal=world_wrap(v_add(v_add(players[b].pos,v_mul(forward,anchor?3:-1.4f)),v_mul(right,anchor?(rank-1)*1.8f:(rank?.9f:-.9f))));goal.y=p->pos.y;
 Vector3 move=world_delta(goal,p->pos);move.y=0;
 if(v_length(move)<.2f){p->vel.x=p->vel.z=0;return;}
 bot_apply_move(p,&botStates[i],move,fminf(MOVE_SPEED*.8f,v_length(move)/fmaxf(dt,.001f)),dt);
}
static bool shield_escort_clear_shot(int i){
 Vector3 direction=player_forward(&players[i]);Ray ray={players[i].pos,direction};
 float range=v_length(world_delta(aiCombat[i].known,players[i].pos));
 for(int j=0;j<activePlayers;j++)if(j!=i && players[j].respawn_timer<=0 && combatantTeam[j]==combatantTeam[i]){
  RayCollision hit=GetRayCollisionBox(ray,expand_bounds(combatant_bounds_at(j,world_image(players[j].pos,ray.position),PLAYER_SIZE*.5f),.25f));
  if(hit.hit && hit.distance>=0 && hit.distance<range)return false;
 }
 return true;
}
static bool shield_escort_default(int i,float dt){
 if(!shield_escort_valid(i) || fp_carrier(i))return false;
 bool anchor=anchor_escort_eligible(shieldEscort[i]-1);
 if(!anchor && (players[i].isExposed || players[i].matter<MATTER_SHOT_COST))return false;
 if(towerRifles[i].tag)tower_release(i);
 ai_tick(i,dt);if(clarity_evade(i,dt))return true;
 shield_escort_move(i,dt);
 if(aiCombat[i].target>=0){
  ai_aim(i,aiCombat[i].known,dt);
  if(!players[i].isExposed && players[i].matter>=MATTER_SHOT_COST && ai_ready_to_fire(i) && shield_escort_clear_shot(i))FireVoxel(i);
 }else {int b=shieldEscort[i]-1;ai_aim(i,v_add(players[b].pos,v_mul(player_forward(&players[b]),10)),dt);}
 return true;
}
