/* Persistent supply stands use the run's shared unlock mask, never automatic pickup. */
static void fp_cache_place(void){
 fpCacheReady=false;if(!fp_active() || netTransport.role==NET_ROLE_CLIENT)return;
 fpCacheOrigin=players[0].pos;
 for(int r=0;r<FP_RULE_COUNT;r++){
  FfItem stand={.pos=world_wrap(v_add(fpCacheOrigin,(Vector3){((r%4)-1.5f)*2.8f,0,4+(r/4)*2.8f}))};
  fp_ground(&stand);fpCachePositions[r]=stand.pos;
 }
 fpCacheReady=true;
}
static bool fp_cache_eligible(int i,int type){
 Equipment *e=&equipment[i];int gun=eq_pickup_gun(type);
 if(gun>=0)return e->gun[0]!=gun && e->gun[1]!=gun;
 if(type==PICKUP_ADHESIVE)return !(e->upgrades&EQ_ADHESIVE);
 if(type==PICKUP_REINFORCED)return !(e->upgrades&EQ_BUILD);
 if(type==PICKUP_SHIELD_REGEN)return !(e->upgrades&EQ_SHIELD);
 if(type==PICKUP_GOLD_TETHER)return !players[i].goldTetherCharged;
 int kind=type==PICKUP_BLACK_HOLE?GRENADE_GRAVITY:type==PICKUP_IGLOO?GRENADE_IGLOO:GRENADE_REGULAR;
 return ex_grenade_count(i,kind)<2;
}
static int fp_cache_near(int i){
 if(!fp_active() || !fpCacheReady || i<0 || i>=firefightHumans || players[i].respawn_timer>0 || combatantLife[i]!=LIFE_ALIVE)return -1;
 int best=-1;float nearest=1.8f;
 for(int r=0;r<FP_RULE_COUNT;r++)if((fpUnlocked&(1u<<r)) && fp_cache_eligible(i,fpRules[r].type)){
  Vector3 d=world_delta(fpCachePositions[r],players[i].pos);float distance=v_length(d);VoxelHit hit;
  if(distance<nearest && (distance<.15f || !first_voxel_hit_detailed((Ray){players[i].pos,v_norm(d)},distance-.15f,-1,false,&hit))){best=r;nearest=distance;}
 }
 return best;
}
static bool fp_cache_collect(int i){
 if(netTransport.role==NET_ROLE_CLIENT)return false;
 int r=fp_cache_near(i);if(r<0 || !fp_reserve(fpCount+1))return false;
 size_t n=fpCount++;fpItems[n]=(FfItem){.id=fpNextId++,.type=fpRules[r].type,.owner=-1,.slot=-1,.pos=fpCachePositions[r]};
 if(!fp_collect(i,n,true)){fpCount--;return false;}
 return true;
}
