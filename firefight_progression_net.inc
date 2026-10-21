static uint32_t fpIncomingRevision,fpIncomingCount,fpIncomingNext;
static uint32_t fpSentRevision,fpReceivedRevision,fpActorTick[MAX_COMBATANTS],fpStateTick;
static void fp_send(int peer,NetWriter *w,bool reliable){
 if(w->failed){TraceLog(LOG_ERROR,"Firefight packet overflow");return;}
 if(peer>=0)net_transport_send(&netTransport,peer,FPS_NET_CHANNEL_WORLD,w->data,w->length,true);
 else net_transport_broadcast(&netTransport,reliable?FPS_NET_CHANNEL_WORLD:FPS_NET_CHANNEL_SNAPSHOT,w->data,w->length,reliable);
}
static void fp_snapshot(int peer){
 if(netTransport.role!=NET_ROLE_HOST || !fp_active())return;
 uint8_t data[FPS_NET_MAX_PACKET];NetWriter w;
 net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_FIREFIGHT,0,netTransport.session_id,netServerTick);net_write_u8(&w,0);
 net_write_u32(&w,fpGeneration);net_write_u8(&w,firefightCheckpoints);net_write_u8(&w,firefightEndless);net_write_u8(&w,firefightFoundry);net_write_u8(&w,firefightHumans);
 net_write_i32(&w,firefightWave);net_write_u8(&w,firefightWaveState);net_write_i32(&w,firefightTeamLives);
 net_write_i32(&w,firefightScore);net_write_i32(&w,firefightKillsTotal);net_write_i32(&w,firefightWaveEnemiesTotal);net_write_i32(&w,firefightWaveEnemiesRemaining);net_write_i32(&w,firefightEnemiesSpawned);
 net_write_f32(&w,firefightIntermissionTimer);net_write_f32(&w,fpRestartTimer);net_write_u32(&w,fpUnlocked);
 for(int k=0;k<3;k++)net_write_u8(&w,teamUpgradeLevels[k]);
 for(int i=0;i<MAX_PLAYERS;i++){net_write_u8(&w,combatantLife[i]);net_write_u8(&w,prepReady[i]);net_write_u8(&w,prepVotes[i]);net_write_u8(&w,prepMenuOpen[i]);net_write_f32(&w,reviveProgress[i]);net_write_u8(&w,reviveOwner[i]+1);}
 fp_send(peer,&w,false);
 if(peer>=0 || fpSentRevision!=fpRevision){
  for(size_t start=0;start<fpCount || start==0;start+=32){
   int count=(int)fmin(32,fpCount-start);net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_FIREFIGHT,0,netTransport.session_id,netServerTick);net_write_u8(&w,1);
   net_write_u32(&w,fpRevision);net_write_u32(&w,(uint32_t)start);net_write_u32(&w,(uint32_t)fpCount);net_write_u8(&w,count);
   for(int n=0;n<count;n++){FfItem *item=&fpItems[start+n];net_write_u32(&w,item->id);net_write_u8(&w,item->type);net_write_u8(&w,item->owner+2);net_write_u8(&w,item->slot+1);eq_write_vec(&w,item->pos);}
   fp_send(peer,&w,true);if(!fpCount)break;
  }if(peer<0)fpSentRevision=fpRevision;
 }
 for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++){
  Player *p=&players[i];Equipment *e=&equipment[i];SwordWave *s=&swordWaves[i];
  net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_FIREFIGHT,0,netTransport.session_id,netServerTick);net_write_u8(&w,2);net_write_u8(&w,i);
  net_write_u8(&w,combatantLife[i]);net_write_u8(&w,p->enemyType);net_write_u8(&w,fpCarrierType[i]);
  eq_write_vec(&w,p->pos);eq_write_vec(&w,p->vel);net_write_f32(&w,p->yaw);net_write_f32(&w,p->pitch);net_write_f32(&w,p->matter);net_write_f32(&w,p->respawn_timer);
  net_write_u8(&w,p->goldTetherCharged);net_write_u8(&w,p->goldTetherHolding);Vector3 target={0};bool tether=player_tether_visual_target(i,&target);net_write_u8(&w,tether);eq_write_vec(&w,target);
  net_write_u8(&w,e->gun[0]);net_write_u8(&w,e->upgrades);net_write_f32(&w,e->charge);net_write_u8(&w,e->held);net_write_f32(&w,e->buildTime);eq_write_vec(&w,e->beamEnd);
  net_write_u32(&w,weaponVisuals[i].shot_sequence);float age=(float)GetTime()-weaponVisuals[i].shot_time;net_write_f32(&w,clampf(age,0,1000));
  net_write_u8(&w,s->active);net_write_u8(&w,s->emitted);net_write_u32(&w,s->sequence);net_write_u8(&w,s->direction);net_write_f32(&w,s->age);eq_write_vec(&w,s->origin);eq_write_vec(&w,s->forward);for(int k=0;k<SWORD_LANES;k++)net_write_f32(&w,s->stop[k]);
  EnemySpecialState *special=&enemySpecial[i];net_write_u8(&w,special->coreOpen);net_write_u8(&w,special->lostMask);net_write_f32(&w,special->coreHealth);net_write_f32(&w,special->coreMax);
  net_write_f32(&w,special->windup);net_write_f32(&w,special->disabled);net_write_f32(&w,special->healingPulse);net_write_u8(&w,special->state);eq_write_vec(&w,special->aim);eq_write_vec(&w,special->buildOrigin);net_write_u8(&w,special->buildKind);
  net_write_u8(&w,special->equipmentCount);for(int n=0;n<special->equipmentCount;n++){eq_write_vec(&w,special->equipment[n].offset);net_write_u8(&w,special->equipment[n].group);}
  fp_send(peer,&w,false);
 }
}
static bool fp_receive(NetMessageType type,NetReader *r,uint32_t tick){
 if(type==NET_MSG_FIREFIGHT_RESET){
  if(netTransport.role!=NET_ROLE_CLIENT)return true;uint32_t generation=net_read_u32(r),session=net_read_u32(r);
  if(r->failed || generation<=fpGeneration || !session)return true;
  netTransport.session_id=session;fpGeneration=generation;fpCount=0;fpUnlocked=0;fpReceivedRevision=fpStateTick=0;memset(fpActorTick,0,sizeof(fpActorTick));
  free(fpIncoming);fpIncoming=NULL;fpIncomingNext=fpIncomingCount=0;fpBannerRead=fpBannerWrite=0;fpBannerTimer=0;
  memset(netPredictionCount,0,sizeof(netPredictionCount));memset(netPredictionHead,0,sizeof(netPredictionHead));
  clarity_reset();eq_reset_all();memset(holograms,0,sizeof(holograms));netWorldReady=false;return true;
 }
 if(type!=NET_MSG_FIREFIGHT)return false;if(netTransport.role!=NET_ROLE_CLIENT)return true;
 int part=net_read_u8(r);
 if(part==0){
  if(tick<fpStateTick)return true;uint32_t generation=net_read_u32(r);bool checkpoints=net_read_u8(r),endless=net_read_u8(r),foundry=net_read_u8(r);int humans=net_read_u8(r);
  int wave=net_read_i32(r),state=net_read_u8(r),lives=net_read_i32(r),score=net_read_i32(r),kills=net_read_i32(r),total=net_read_i32(r),remaining=net_read_i32(r),spawned=net_read_i32(r);
  float prep=net_read_f32(r),restart=net_read_f32(r);uint32_t unlocked=net_read_u32(r);unsigned upgrades[3];for(int k=0;k<3;k++)upgrades[k]=net_read_u8(r);
  int life[4],ready[4],votes[4],menus[4],owners[4];float progress[4];for(int i=0;i<4;i++){life[i]=net_read_u8(r);ready[i]=net_read_u8(r);votes[i]=net_read_u8(r);menus[i]=net_read_u8(r);progress[i]=net_read_f32(r);owners[i]=net_read_u8(r)-1;}
  for(int n=0;n<4;n++)if(life[n]>LIFE_UNUSED||ready[n]>1||menus[n]>1||(votes[n]>2&&votes[n]!=255)||!isfinite(progress[n])||progress[n]<0||owners[n]>=4)return true;
  if(r->failed || humans<1||humans>4||wave<1||state>FIREFIGHT_STATE_RESTART||!isfinite(prep)||!isfinite(restart))return true;
  gameMode=GAME_MODE_FIREFIGHT;firefightPractice=false;activePlayers=MAX_COMBATANTS;firefightHumans=humans;firefightCheckpoints=checkpoints;firefightEndless=endless;firefightFoundry=foundry;
  firefightWave=wave;firefightWaveState=state;firefightTeamLives=lives;firefightScore=score;firefightKillsTotal=kills;firefightWaveEnemiesTotal=total;firefightWaveEnemiesRemaining=remaining;firefightEnemiesSpawned=spawned;
  firefightIntermissionTimer=prep;fpRestartTimer=restart;memcpy(teamUpgradeLevels,upgrades,sizeof(upgrades));
  if(fpStateTick)for(int n=0;n<FP_RULE_COUNT;n++)if((unlocked&(1u<<n)) && !(fpUnlocked&(1u<<n)))fp_queue_unlock(fpRules[n].type);fpUnlocked=unlocked;fpGeneration=generation;fpStateTick=tick;
  for(int i=0;i<MAX_COMBATANTS;i++){combatantTeam[i]=i<MAX_PLAYERS?0:1;if(i<4){combatantLife[i]=life[i];prepReady[i]=ready[i];prepVotes[i]=votes[i]==255?-1:votes[i];prepMenuOpen[i]=menus[i];reviveProgress[i]=progress[i];reviveOwner[i]=owners[i];}else playerInput[i]=INPUT_TYPE_BOT_EASY;}
  if(state==FIREFIGHT_STATE_DEFEAT||state==FIREFIGHT_STATE_VICTORY)gameState=GAME_STATE_GAMEOVER;else if(netWorldReady)gameState=GAME_STATE_PLAYING;
 }else if(part==1){
  uint32_t revision=net_read_u32(r),start=net_read_u32(r),total=net_read_u32(r);int count=net_read_u8(r);
  if(r->failed||revision<fpReceivedRevision||total>1000000||count>32||start>total||(uint32_t)count>total-start)return true;
  FfItem items[32];for(int n=0;n<count;n++){items[n].id=net_read_u32(r);items[n].type=net_read_u8(r);items[n].owner=net_read_u8(r)-2;items[n].slot=net_read_u8(r)-1;items[n].pos=eq_read_vec(r);
   if(!items[n].id||fp_rule(items[n].type)<0||items[n].owner>=MAX_COMBATANTS||items[n].slot>1||!v_isfinite(items[n].pos))return true;
  }if(r->failed)return true;
  if(start==0){FfItem *next=total?malloc(total*sizeof(FfItem)):NULL;if(total&&!next)return true;free(fpIncoming);fpIncoming=next;fpIncomingCount=total;fpIncomingRevision=revision;fpIncomingNext=0;}
  if(revision!=fpIncomingRevision||total!=fpIncomingCount||start!=fpIncomingNext)return true;
  if(count)memcpy(fpIncoming+start,items,count*sizeof(FfItem));fpIncomingNext+=count;
  if(fpIncomingNext==total){free(fpItems);fpItems=fpIncoming;fpIncoming=NULL;fpCount=fpCapacity=total;fpReceivedRevision=revision;}
 }else if(part==2){
  int i=net_read_u8(r);if(i<MAX_PLAYERS||i>=MAX_COMBATANTS||tick<fpActorTick[i])return true;
  int life=net_read_u8(r),enemy=net_read_u8(r),carrier=net_read_u8(r);Vector3 pos=eq_read_vec(r),vel=eq_read_vec(r);
  float yaw=net_read_f32(r),pitch=net_read_f32(r),matter=net_read_f32(r),respawn=net_read_f32(r);bool gold=net_read_u8(r),holding=net_read_u8(r),tether=net_read_u8(r);Vector3 target=eq_read_vec(r);
  Equipment e={.initialized=true,.gun={EQ_DEFAULT,EQ_EMPTY}};e.gun[0]=net_read_u8(r);e.upgrades=net_read_u8(r);e.charge=net_read_f32(r);e.held=net_read_u8(r);e.buildTime=net_read_f32(r);e.beamEnd=eq_read_vec(r);
  uint32_t shot=net_read_u32(r);float shotAge=net_read_f32(r);SwordWave s={0};s.active=net_read_u8(r);s.emitted=net_read_u8(r);s.sequence=net_read_u32(r);s.direction=net_read_u8(r);s.age=net_read_f32(r);s.origin=eq_read_vec(r);s.forward=eq_read_vec(r);for(int k=0;k<SWORD_LANES;k++)s.stop[k]=net_read_f32(r);
  EnemySpecialState special={0};special.coreOpen=net_read_u8(r);special.lostMask=net_read_u8(r);special.coreHealth=net_read_f32(r);special.coreMax=net_read_f32(r);
  special.windup=net_read_f32(r);special.disabled=net_read_f32(r);special.healingPulse=net_read_f32(r);special.state=net_read_u8(r);special.aim=eq_read_vec(r);special.buildOrigin=eq_read_vec(r);special.buildKind=net_read_u8(r);special.equipmentCount=net_read_u8(r);
  if(special.equipmentCount>ENEMY_EQUIPMENT_MAX)return true;for(int n=0;n<special.equipmentCount;n++){special.equipment[n].offset=eq_read_vec(r);special.equipment[n].group=net_read_u8(r);if(!v_isfinite(special.equipment[n].offset))return true;}
  if(!isfinite(pitch)||!isfinite(e.charge)||e.charge<0||e.charge>eqRailCharge+.01f||!isfinite(e.buildTime)||!v_isfinite(e.beamEnd)||!v_isfinite(target)||!isfinite(shotAge)||!isfinite(s.age)||!v_isfinite(s.origin)||!v_isfinite(s.forward)||s.direction>1)return true;
  if(!isfinite(special.coreHealth)||!isfinite(special.coreMax)||!isfinite(special.windup)||!isfinite(special.disabled)||!isfinite(special.healingPulse)||!v_isfinite(special.aim)||!v_isfinite(special.buildOrigin))return true;
  for(int n=0;n<SWORD_LANES;n++)if(!isfinite(s.stop[n]))return true;
  if(r->failed||enemy>=ENEMY_TYPE_COUNT||life>LIFE_UNUSED||(carrier&&fp_rule(carrier)<0)||e.gun[0]>=EQ_GUN_COUNT||!v_isfinite(pos)||!v_isfinite(vel)||!isfinite(matter)||!isfinite(yaw)||!isfinite(respawn))return true;
  enemySpecial[i]=special;Player *p=&players[i];p->pos=pos;p->vel=vel;p->yaw=yaw;p->pitch=pitch;p->matter=matter;p->isExposed=matter<=0;p->respawn_timer=respawn;p->enemyType=enemy;
  p->goldTetherCharged=gold;p->goldTetherHolding=holding;p->netTetherVisualActive=tether;p->netGoldTetherVisualActive=holding;p->netTetherVisualTarget=target;
  combatantLife[i]=life;combatantTeam[i]=1;playerInput[i]=INPUT_TYPE_BOT_EASY;fpCarrierType[i]=carrier;equipment[i]=e;swordWaves[i]=s;swordDirection[i]=s.direction;
  if(s.active){p->swordSwing=true;p->meleeSwingActive=true;p->meleeSwingStartTime=(float)GetTime()-s.age;}else p->meleeSwingActive=false;
  weapon_visual_receive(&weaponVisuals[i],shot,shotAge<.12f?(uint8_t)(shotAge*1000):255,(float)GetTime());fpActorTick[i]=tick;
 }return true;
}
static void fp_network_reset(void){
 if(netTransport.role!=NET_ROLE_HOST)return;uint8_t data[40];NetWriter w;uint32_t next=netTransport.session_id+1;if(!next)next=1;
 net_writer_init(&w,data,sizeof(data));net_write_header(&w,NET_MSG_FIREFIGHT_RESET,0,netTransport.session_id,netServerTick);net_write_u32(&w,fpGeneration);net_write_u32(&w,next);fp_send(-1,&w,true);
 netTransport.session_id=next;memset(netHasPendingInput,0,sizeof(netHasPendingInput));memset(netPreviousHeld,0,sizeof(netPreviousHeld));memset(netLastProcessedInput,0,sizeof(netLastProcessedInput));fpSentRevision=0;
 for(int peer=1;peer<MAX_PLAYERS;peer++)if(netTransport.peers[peer].welcomed)net_send_world_to(peer);
 netLastStaticGenerationSent=staticHashGeneration;
}
