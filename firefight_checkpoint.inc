static void fp_capture(void){
 if(!fp_active()||!firefightCheckpoints||fpRestoring||netTransport.role==NET_ROLE_CLIENT)return;
 FfItem *copy=fpCount?malloc(fpCount*sizeof(*copy)):NULL;
 if(fpCount&&!copy){TraceLog(LOG_ERROR,"Checkpoint allocation failed");fpCheckpoint.valid=false;return;}
 if(fpCount)memcpy(copy,fpItems,fpCount*sizeof(*copy));free(fpCheckpoint.items);
 fpCheckpoint=(FfCheckpoint){.valid=true,.wave=firefightWave,.lives=firefightTeamLives,.score=firefightScore,.kills=firefightKillsTotal,
  .repeatCursor=fpRepeatCursor,.nextId=fpNextId,.unlocked=fpUnlocked,.items=copy,.count=fpCount};
 memcpy(fpCheckpoint.remaining,ffRemaining,sizeof(ffRemaining));memcpy(fpCheckpoint.pending,fpPending,sizeof(fpPending));fpCheckpoint.pendingCount=fpPendingCount;fpCheckpoint.total=firefightWaveEnemiesTotal;
 for(int k=0;k<3;k++)fpCheckpoint.upgrades[k]=teamUpgradeLevels[k];
 for(int i=0;i<firefightHumans;i++){
  eq_cancel_charge(i);fpCheckpoint.humans[i]=players[i];fpCheckpoint.equipment[i]=equipment[i];
  fpCheckpoint.grenades[i]=grenadeInventory[i];fpCheckpoint.selected[i]=exSelectedGrenade[i];fpCheckpoint.revives[i]=reviveCount[i];
  memcpy(fpCheckpoint.special[i],exGrenades[i],sizeof(fpCheckpoint.special[i]));
 }
}
static bool fp_retry(void){
 if(!fp_active()||!firefightCheckpoints||!fpCheckpoint.valid)return false;
 if(fpRestartTimer<=0)fpRestartTimer=3;
 firefightWaveState=FIREFIGHT_STATE_RESTART;return true;
}
static void fp_restore(void){
 if(!fpCheckpoint.valid || !fp_reserve(fpCheckpoint.count))return;
 fpRestoring=true;ResetGame();gameState=GAME_STATE_PLAYING;
 fpCount=fpCheckpoint.count;if(fpCount)memcpy(fpItems,fpCheckpoint.items,fpCount*sizeof(*fpItems));
 fpNextId=fpCheckpoint.nextId;fpUnlocked=fpCheckpoint.unlocked;fpRepeatCursor=fpCheckpoint.repeatCursor;
 firefightScore=fpCheckpoint.score;firefightKillsTotal=fpCheckpoint.kills;firefightTeamLives=fpCheckpoint.lives;
 for(int k=0;k<3;k++)teamUpgradeLevels[k]=fpCheckpoint.upgrades[k];
 for(int i=0;i<firefightHumans;i++){
  Player saved=fpCheckpoint.humans[i];players[i].kills=saved.kills;players[i].deaths=saved.deaths;players[i].debrisKills=saved.debrisKills;players[i].kd_ratio=saved.kd_ratio;
  equipment[i]=fpCheckpoint.equipment[i];equipment[i].charge=equipment[i].switchTime=equipment[i].buildTime=0;equipment[i].held=equipment[i].buildHeld=false;
  memset(equipment[i].cooldown,0,sizeof(equipment[i].cooldown));ff_restore_human(i,MATTER_MAX_DEFAULT);
  players[i].goldTetherCharged=saved.goldTetherCharged;players[i].dynamicShotActive=eq_gun(i)==EQ_VOXEL;
  grenadeInventory[i]=fpCheckpoint.grenades[i];exSelectedGrenade[i]=fpCheckpoint.selected[i];reviveCount[i]=fpCheckpoint.revives[i];
  memcpy(exGrenades[i],fpCheckpoint.special[i],sizeof(fpCheckpoint.special[i]));
 }
 for(size_t n=0;n<fpCount;n++){if(netTransport.role==NET_ROLE_HOST&&fpItems[n].owner>=0&&fpItems[n].owner<MAX_PLAYERS&&!netPlayerPresent[fpItems[n].owner])fp_drop(&fpItems[n],fpCheckpoint.humans[fpItems[n].owner].pos);else if(fpItems[n].owner==-1)fp_ground(&fpItems[n]);}
 // Keep the original wave roster even if a LAN seat joined or disconnected.
 start_firefight_wave(fpCheckpoint.wave);
 memcpy(ffRemaining,fpCheckpoint.remaining,sizeof(ffRemaining));memcpy(fpPending,fpCheckpoint.pending,sizeof(fpPending));fpPendingCount=fpCheckpoint.pendingCount;
 firefightWaveEnemiesTotal=firefightWaveEnemiesRemaining=fpCheckpoint.total;fpRepeatCursor=fpCheckpoint.repeatCursor;
 fpRestoring=false;fpRestartTimer=0;++fpGeneration;++fpRevision;fp_network_reset();
}
