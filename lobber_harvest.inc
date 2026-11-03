/* Native lobbers recover through the same timed melee/terrain path as players. */
static bool sapper_hidden(int i,Vector3 eye,Vector3 threat,uint64_t *cover);
static bool sapper_ground(int i,Vector3 *eye);
static bool lobber_cover(int i,float dt,Vector3 *move){
 EnemySpecialState *s=&enemySpecial[i];Player *p=&players[i];
 s->coverSearch-=dt;
 if(aiCombat[i].target<0){s->coverSupport=0;return false;}
 int support=enemy_voxel(s->coverSupport);
 if(support<0 || voxels[support].simulate || !sapper_hidden(i,s->coverPosition,aiCombat[i].known,NULL))s->coverSupport=0;
 if(s->coverSearch<=0){
  s->coverSearch=1.0f+.05f*(i%4);s->coverSupport=0;float best=FLT_MAX;
  for(int ring=0;ring<=4;ring++)for(int n=0;n<(ring?16:1);n++){
   float angle=n*2*PI/16;Vector3 at=world_wrap(v_add(p->pos,(Vector3){cosf(angle)*ring*2,0,sinf(angle)*ring*2}));
   float range=v_length(world_delta(at,aiCombat[i].known));uint64_t id;
   if(range<7 || !sapper_ground(i,&at) || !sapper_hidden(i,at,aiCombat[i].known,&id))continue;
   float score=v_length(world_delta(at,p->pos))+.1f*fabsf(range-16);
   if(score<best){best=score;s->coverPosition=at;s->coverSupport=id;}
  }
 }
 if(!s->coverSupport)return false;
 *move=world_delta(s->coverPosition,p->pos);move->y=0;
 if(v_length(*move)<.5f)*move=(Vector3){0};
 s->reason="lobbing from cover";return true;
}
static bool lobber_mineable(int i,int v){
 if(v<0 || v>=voxel_count)return false;
 Voxel *cell=&voxels[v];Player *p=&players[i];
 return !cell->simulate && !cell->isBullet && !cell->pendingActivation && !cell->towerTag && cell->type==0 &&
  cell->pos.y+VOXEL_SIZE*.5f>=p->pos.y-.25f && cell->pos.y-p->pos.y<=MELEE_RANGE*.8f &&
  v_length(world_delta(cell->pos,p->pos))<=lobberMineRange;
}
static int lobber_find_mine(int i){
 int best=-1;float nearest=lobberMineRange;EnemySpecialState *s=&enemySpecial[i];
 for(int v=0;v<voxel_count;v++)if(lobber_mineable(i,v) && voxels[v].identity!=s->harvestRejected){
  Vector3 delta=world_delta(voxels[v].pos,players[i].pos);float distance=v_length(delta);
  if(distance>=nearest)continue;
  VoxelHit obstruction;
  if(distance>.3f && first_voxel_hit_detailed((Ray){players[i].pos,v_norm(delta)},distance-.26f,-1,false,&obstruction) && obstruction.id!=v)continue;
  nearest=distance;best=v;
 }
 return best;
}
static bool lobber_harvest(int i,float dt){
 Player *p=&players[i];EnemySpecialState *s=&enemySpecial[i];
 if(p->enemyType!=ENEMY_TYPE_LOBBER && !(p->enemyType==ENEMY_TYPE_STANDARD && fp_carrier(i)))return false;
 if(!s->harvesting && (p->isExposed || p->matter<player_max_matter(p)*lobberMineStart)){
  s->harvesting=true;s->harvestSearch=s->harvestStuck=0;s->harvestTarget=s->harvestRejected=0;s->harvestLast=p->pos;
  enemy_drop_payload(i);s->windup=0;s->state=0;s->selected=0;
 }
 if(!s->harvesting)return false;
 if(p->matter>=player_max_matter(p)*lobberMineFinish){
  if(p->meleeSwingActive){p->vel.x=p->vel.z=0;return true;}
  s->harvesting=false;s->harvestTarget=0;s->searchTimer=0;s->reason="recovered";return false;
 }
 s->harvestSearch=fmaxf(0,s->harvestSearch-dt);
 int v=enemy_voxel(s->harvestTarget);
 if(!lobber_mineable(i,v)){
  s->harvestTarget=0;v=-1;
  if(s->harvestSearch<=0){v=lobber_find_mine(i);s->harvestTarget=v>=0?voxels[v].identity:0;s->harvestSearch=.75f;s->harvestStuck=0;s->harvestLast=p->pos;if(v<0)s->harvestRejected=0;}
 }
 if(v<0){s->reason="seeking harvestable terrain";return false;}
 Vector3 to=world_delta(voxels[v].pos,p->pos);float distance=v_length(to),yaw,pitch;
 dir_to_yaw_pitch(to,&yaw,&pitch);
 p->yaw+=angle_diff_deg(p->yaw,yaw)*fminf(1,12*dt);p->pitch+=(pitch-p->pitch)*fminf(1,12*dt);
 s->work=1;s->reason="harvesting Matter";
 if(distance<MELEE_RANGE*.8f){
  p->vel.x=p->vel.z=0;s->harvestStuck=0;
  if(v_dot(player_forward(p),v_norm(to))>.94f)perform_melee(i);
 }else{
  Vector3 progress=world_delta(p->pos,s->harvestLast);progress.y=0;
  if(v_length(progress)>.15f){s->harvestStuck=0;s->harvestLast=p->pos;}else s->harvestStuck+=dt;
  if(s->harvestStuck>3){s->harvestRejected=s->harvestTarget;s->harvestTarget=0;s->harvestSearch=0;}
  bot_apply_move(p,&botStates[i],to,MOVE_SPEED*.7f,dt);
 }
 return true;
}
