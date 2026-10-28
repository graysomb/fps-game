static void fp_update_carrier(int i,float dt){
 if(players[i].enemyType==ENEMY_TYPE_BUILDER){tower_builder(i,dt);return;}
 if(tower_rifle(i,dt))return;
 Player *p=&players[i];Equipment *e=&equipment[i];int type=fpCarrierType[i];
 ai_tick(i,dt);if(clarity_evade(i,dt)){eq_cancel_charge(i);e->held=false;return;}
 AiCombat *a=&aiCombat[i];fpCarrierClock[i]=fmaxf(0,fpCarrierClock[i]-dt);
 if(a->target<0){
  eq_cancel_charge(i);BotState *b=&botStates[i];b->pathTimer-=dt;
  if(b->pathTimer<=0){b->pathTimer=3; a->known=world_wrap(v_add(p->pos,(Vector3){cosf(ai_random(i)*2*PI)*8,0,sinf(ai_random(i)*2*PI)*8}));}
  Vector3 wander=world_delta(a->known,p->pos);wander.y=0;bot_apply_move(p,b,wander,MOVE_SPEED*.45f,dt);ai_aim(i,a->known,dt);return;
 }
 Vector3 delta=world_delta(a->known,p->pos);float distance=v_length(delta);ai_aim(i,a->known,dt);
 int gun=eq_gun(i);float desired=gun==EQ_BLADE?SWORD_REACH*.65f:gun==EQ_SHOTGUN?6:gun==EQ_CUTTER?12:gun==EQ_RAIL||gun==EQ_DMR?22:10;
 Vector3 move=delta;move.y=0;
 if(distance>desired)bot_apply_move(p,&botStates[i],move,MOVE_SPEED*.8f,dt);else p->vel.x=p->vel.z=0;
 if(!a->visible || a->reaction>0){eq_cancel_charge(i);return;}
 bool supply=type==PICKUP_GRENADE||type==PICKUP_ADHESIVE||type==PICKUP_BLACK_HOLE||type==PICKUP_IGLOO;
 float needed=type==PICKUP_GOLD_TETHER?eqGoldGrabCost:eqGuns[gun].cost;
 if(!supply && p->matter<needed && e->charge<=0){
  float d2;int v=find_nearest_static_voxel(&p->pos,&d2);
  if(v>=0){Vector3 to=world_delta(voxels[v].pos,p->pos);ai_aim(i,world_image(voxels[v].pos,p->pos),dt);bot_apply_move(p,&botStates[i],to,MOVE_SPEED*.6f,dt);if(d2<MELEE_RANGE*MELEE_RANGE)perform_melee(i);}
  else if(distance<MELEE_RANGE)perform_melee(i);return;
 }
 if(type==PICKUP_GOLD_TETHER){
  if(!p->tetherHolding && fpCarrierClock[i]<=0 && p->matter>=eqGoldGrabCost){
   fpCarrierClock[i]=.5f;int v=enemy_voxel(enemy_find_payload(i));if(v>=0 && v_length(world_delta(voxels[v].pos,p->pos))<=GOLD_TETHER_RANGE){
    acquire_tether_hit(i,v,(Vector3){0,1,0});if(p->tetherHolding)fpCarrierClock[i]=.7f;
   }
  }else if(p->tetherHolding&&fpCarrierClock[i]<=0){release_tether(i);fpCarrierClock[i]=3;}
  return;
 }
 if(supply){
  if(!ex_grenade_count(i,exSelectedGrenade[i])){bot_apply_move(p,&botStates[i],move,MOVE_SPEED*.8f,dt);if(distance<MELEE_RANGE)perform_melee(i);return;}
  if(fpCarrierClock[i]<=0 && distance<24){float pitch=p->pitch;
   if(type==PICKUP_IGLOO)p->pitch=-55;
   else {float flight=fmaxf(.2f,sqrtf(delta.x*delta.x+delta.z*delta.z)/14);p->pitch=clampf(atan2f(delta.y/flight+4.9f*flight-2,14)*RAD2DEG,-40,60);}
   if(clarity_throw(i))fpCarrierClock[i]=3;p->pitch=pitch;
  }return;
 }
 if(p->isExposed && e->charge<=0){
  if(distance<MELEE_RANGE)perform_melee(i);return;
 }
 if(gun==EQ_BLADE){if(distance<SWORD_REACH*.9f)eq_fire(i);return;}
 if(gun==EQ_ROCKET && distance<4.5f)return;
 if(distance>eqGuns[gun].range)return;
 if(v_dot(player_forward(p),v_norm(delta))<cosf(8*DEG2RAD))return;
 if(gun==EQ_RAIL){NetInputCommand c={0};if(e->charge<eqRailCharge)c.held=NET_INPUT_FIRE;eq_command(i,&c,dt);return;}
 if(a->recovery>0)return;
 if(gun==EQ_CUTTER||gun==EQ_GLOB){NetInputCommand c={.held=NET_INPUT_FIRE};eq_command(i,&c,dt);return;}
 if(gun==EQ_VOXEL){FireVoxel(i);return;}
 if(eq_fire(i))ai_shot_fired(i);
}
