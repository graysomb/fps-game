#ifndef FIREFIGHT_RULES_H
#define FIREFIGHT_RULES_H
#include <stdbool.h>
#define FF_SEATS 4
#define FF_ENEMIES 12
#define FF_ACTORS 16
#define FF_BLEED_SECONDS 15.0f
#define FF_REVIVE_SECONDS 3.0f
#define FF_REVIVE_COST 25.0f
#define FF_REVIVE_MATTER 35.0f
#define FF_PREP_SECONDS 20.0f
typedef enum { LIFE_ALIVE, LIFE_DOWNED, LIFE_WAITING, LIFE_SPECTATING, LIFE_UNUSED } CombatantLifeState;
typedef enum { UPGRADE_BUILD, UPGRADE_HARVEST, UPGRADE_THROW, UPGRADE_NONE } TeamUpgrade;
typedef struct { const char *name; int count[3]; } WaveDefinition;
static const WaveDefinition ff_waves[5] = {
    /* Storage order is normal, swarmer, Goliath; difficulty starts with swarmers. */
    {"First Contact", {0,3,0}}, {"Swarm Patrol", {0,6,0}},
    {"Heavy Contact", {1,6,0}}, {"First Goliath", {2,6,1}},
    {"Last Stand", {3,9,2}}
};
static inline int ff_count(int wave, int humans, int type) {
    int base = ff_waves[(wave-1)%5].count[type];
    int cycle = (wave-1)/5;
    return (base * (humans+1) * (4+cycle) + 7)/8;
}
static inline int ff_cap(int humans) { int n=4+2*humans; return n>12?12:n; }
static inline int ff_wave_cap(int wave,int humans) {
    if(wave==1) return humans+1;
    if(wave==2) return humans+2;
    if(wave==3) return humans+3;
    return ff_cap(humans);
}
static inline TeamUpgrade ff_vote(const int *votes, int humans) {
    int counts[3]={0};
    for(int i=0;i<humans;i++) if(votes[i]>=0 && votes[i]<3) counts[votes[i]]++;
    int best=0;
    for(int i=1;i<3;i++) if(counts[i]>counts[best]) best=i;
    for(int i=0;i<3;i++) if(i!=best && counts[i]==counts[best]) return UPGRADE_BUILD;
    return (TeamUpgrade)best;
}
#endif
