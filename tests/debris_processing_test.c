#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fresh(void){
 clear_world_voxels();activePlayers=0;worldPeriodic=false;memset(players,0,sizeof(players));init_static_hash();
}
static void candidates(void){
 fresh();for(int x=0;x<8;x++)add_static_voxel_at_grid(x,0,0,GRAY,0);
 rebuild_static_hash_if_dirty();for(int v=0;v<voxel_count;v++){voxels[v].freezeBelief=.9f;atomic_store(&activationClaims[v],0);}
 activation_candidates_reset();ActivationCandidates *c=&activationCandidates[0];
 assert(activation_candidate_rejected(c,0,3,0,0,100,0,VOXEL_ACTIVATION_UNIT_BUDGET));
 unsigned long long before=c->searches;
 for(int v=0;v<8;v++)assert(activation_candidate_rejected(c,v,3,0,0,100,0,VOXEL_ACTIVATION_UNIT_BUDGET));
 assert(c->searches==before);
 assert(!activation_candidate_rejected(c,0,3,0,0,100,0,4));
 assert(!activation_candidate_rejected(c,0,3,0,0,100,100,VOXEL_ACTIVATION_UNIT_BUDGET));
 static UnitVoxelBuffer b;unit_voxel_buffer_clear(&b);
 int n=collect_static_activation_cluster(0,0,3,0,0,100,&b);assert(n==8);
 float belief=compute_cluster_freeze_belief(&b,0);rollback_activation_buffer(&b,0);
 assert(fabsf(belief-c->resistance[0])<.00001f);
 for(int v=0;v<8;v++)for(int f=0;f<10;f++){
  float force=f*.2f;
  activation_candidates_reset();
  bool cached=activation_candidate_rejected(c,v,3,0,0,100,force,VOXEL_ACTIVATION_UNIT_BUDGET);
  activationCandidateCacheDisabled=true;
  bool uncached=activation_candidate_rejected(c,v,3,0,0,100,force,VOXEL_ACTIVATION_UNIT_BUDGET);
  activationCandidateCacheDisabled=false;
  assert(cached==uncached);
  assert(cached==!dynamic_belief_overcomes_static(force,belief));
 }
 atomic_store(&activationClaims[0],1);atomic_fetch_add(&activationClaimRevision,1);
 assert(!activation_candidate_rejected(c,0,3,0,0,100,0,VOXEL_ACTIVATION_UNIT_BUDGET));
 atomic_store(&activationClaims[0],0);atomic_fetch_add(&activationClaimRevision,1);
 fresh();for(int x=0;x<VOXEL_ACTIVATION_UNIT_BUDGET+1;x++)add_static_voxel_at_grid(x,0,0,GRAY,0);
 rebuild_static_hash_if_dirty();activation_candidates_reset();
 assert(!activation_candidate_rejected(c,0,0,0,0,-1,0,VOXEL_ACTIVATION_UNIT_BUDGET)); /* incomplete: no caching */
 puts("Candidate equivalence, stronger hits, claims, truncation: PASS");
}
static void islands(void){
 for(int held=0;held<4;held++){
  fresh();int a=addVoxel(.25f,.25f,.25f,false,true,GRAY,0),b=addVoxel(.75f,.25f,.25f,false,true,GRAY,0);
  assert(a>=0&&b>=0);voxels[a].owner=voxels[b].owner=ENEMY_SCRAP_OWNER;
  /* Match a shared corner while deliberately leaving glue edges absent. */
  Particle *saved=NULL;int bc=-1;
  for(int x=0;x<8&&bc<0;x++)for(int y=0;y<8;y++)
   if(v_length(v_sub(voxels[a].particles[x]->pos,voxels[b].particles[y]->pos))<.001f){
    bc=y;saved=voxels[b].particles[y];voxels[b].particles[y]=voxels[a].particles[x];break;
   }
  assert(bc>=0);
  voxels[a].sleepFrames=voxels[b].sleepFrames=VOXEL_DEACTIVATION_FRAMES;
  voxels[a].wake_timer=voxels[b].wake_timer=0;
  if(held==2)voxels[b].sleeping=true;
  if(held==3)voxels[b].persistentConstruction=true;
  if(held==1){activePlayers=1;players[0].tetherHolding=true;players[0].tetherVoxel=a;players[0].tetherVoxelIdentity=voxels[a].identity;tetherTag[a]=1;}
  deactivate_sleeping_voxels();assert(sleepGroupsChecked==1);
  if(held==2){assert(!voxels[a].sleeping&&voxels[b].sleeping);assert(sleepMixed==1);}
  else if(held){assert(!voxels[a].sleeping&&!voxels[b].sleeping);assert(sleepGroupsBlocked==1);}
  else {assert(voxels[a].sleeping&&voxels[b].sleeping);assert(sleepGroupsFrozen==1);}
  voxels[b].particles[bc]=saved;
 }
 puts("Shared-particle island sleeps once; held, mixed and persistent islands protected: PASS");
}
static void wakes(void){
 fresh();int a=addVoxel(.25f,2.25f,.25f,false,true,GRAY,0),b=addVoxel(.75f,2.25f,.25f,false,true,GRAY,0);
 voxels[a].sleeping=true;voxels[b].wasTethered=true;rebuild_voxel_hash();
 voxels[b].vgs_active=false;assert(!wake_sleeping_voxels_near_awake());assert(voxels[a].sleeping);
 voxels[b].vgs_active=true;assert(wake_sleeping_voxels_near_awake());assert(!voxels[a].sleeping);
 fresh();a=addVoxel(.25f,2.25f,.25f,false,true,GRAY,0);b=addVoxel(10.25f,2.25f,.25f,false,true,GRAY,0);
 voxels[a].sleeping=true;voxels[b].wasTethered=true;rebuild_voxel_hash();
 assert(!wake_sleeping_voxels_near_awake());assert(voxels[a].sleeping);
 fresh();worldPeriodic=true;worldFloorHalfSize=20;
 a=addVoxel(-19.75f,2.25f,.25f,false,true,GRAY,0);b=addVoxel(19.75f,2.25f,.25f,false,true,GRAY,0);
 voxels[a].sleeping=true;voxels[b].isBullet=true;rebuild_voxel_hash();
 assert(wake_sleeping_voxels_near_awake());assert(!voxels[a].sleeping);
 puts("Contact wake, distant exclusion, periodic bullet wake: PASS");
}
int main(void){
 SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(320,240,"Debris processing");
 candidates();islands();wakes();fresh();CloseWindow();return 0;
}
