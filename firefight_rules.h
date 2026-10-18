#ifndef FIREFIGHT_RULES_H
#define FIREFIGHT_RULES_H
#include <stdbool.h>
#define FF_SEATS 4
#define FF_ENEMIES 12
#define FF_ACTORS 16
#define FF_REVIVE_SECONDS 3.0f
#define FF_REVIVE_COST 0.0f
#define FF_REVIVE_MATTER 35.0f
#define FF_PREP_SECONDS 20.0f
typedef enum { LIFE_ALIVE, LIFE_DOWNED, LIFE_WAITING, LIFE_SPECTATING, LIFE_UNUSED } CombatantLifeState;
typedef enum { UPGRADE_BUILD, UPGRADE_HARVEST, UPGRADE_THROW, UPGRADE_NONE } TeamUpgrade;
#define FF_TYPE_COUNT 10
#define FF_WAVE_COUNT 12
typedef struct { const char *name; int count[FF_TYPE_COUNT]; } WaveDefinition;
static const WaveDefinition ff_waves[FF_WAVE_COUNT] = {
    /* normal, swarmer, Goliath, scrapper, shield, sapper, lobber, builder, anchor, foundry */
    {"First Contact", {0,3}}, {"Heavy Patrol", {1,4,1}},
    {"Demolition", {0,4,0,0,0,1}}, {"Field Works", {1,3,0,0,1,0,0,1}},
    {"Incoming Scrap", {1,4,0,0,0,0,1}}, {"Anchored Assault", {0,3,1,0,0,1,0,0,1}},
    {"Fortified Patrol", {1,4,0,0,1,0,0,1}}, {"Siege Support", {1,4,0,0,0,0,1,0,1}},
    {"Heavy Works", {1,4,1,0,0,1,0,1}}, {"Shielded Siege", {1,4,0,0,1,0,1,0,1}},
    {"Combined Assault", {2,9,2,0,1,1,1,1,1}},
    {"Walking Foundry", {0,6,0,0,0,0,0,0,0,1}}
};
/* AI tiers: 0 easy, 1 medium, 2 hard. Riflemen ramp more slowly. */
static inline int ff_ai_level(int wave, int type) {
    if(type==0) return wave<5?0:wave<10?1:2;
    return wave==1?0:wave<=3?1:2;
}
static inline int ff_count(int wave, int humans, int type) {
    if(wave<1 || type<0 || type>=FF_TYPE_COUNT) return 0;
    int base = ff_waves[(wave-1)%FF_WAVE_COUNT].count[type];
    if(type==0 && (wave-1)%FF_WAVE_COUNT<10) return base;
    if(type==9 || (wave-1)%FF_WAVE_COUNT==FF_WAVE_COUNT-1) return base;
    /* Specialists scale with team size, retaining their existing cycle counts. */
    if(type>=3) return (base * (humans+1) + 1)/2;
    int cycle = (wave-1)/FF_WAVE_COUNT;
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
