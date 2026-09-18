#include "../greek_grammar.h"
#include "../megalith_grammar.h"
#include "../hyperborean_grammar.h"
#include "../forerunner_grammar.h"
#include "../unified_sanctum_grammar.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static uint64_t add(uint64_t h,uint64_t v){h^=v;return h*UINT64_C(1099511628211);}
static uint64_t greek_hash(const TemplePlan *p){uint64_t h=add(UINT64_C(1469598103934665603),p->layout_family);h=add(h,p->optional_feature_mask);for(int i=0;i<p->node_count;i++){const ArchNode*n=&p->nodes[i];h=add(h,n->type);h=add(h,(uint32_t)n->box.min_x);h=add(h,(uint32_t)n->box.min_y);h=add(h,(uint32_t)n->box.min_z);h=add(h,(uint32_t)n->box.max_x);h=add(h,(uint32_t)n->box.max_y);h=add(h,(uint32_t)n->box.max_z);}return h;}
static uint64_t mega_hash(const MegalithPlan*p){uint64_t h=add(UINT64_C(1469598103934665603),p->layout_family);for(int i=0;i<p->node_count;i++){const MegalithNode*n=&p->nodes[i];h=add(h,n->type);h=add(h,(uint32_t)n->box.min_x);h=add(h,(uint32_t)n->box.min_z);h=add(h,(uint32_t)n->box.max_x);h=add(h,(uint32_t)n->box.max_z);}return h;}
static uint64_t hyper_hash(const HyperPlan*p){uint64_t h=add(UINT64_C(1469598103934665603),p->layout_family);for(int i=0;i<p->node_count;i++){const HyperNode*n=&p->nodes[i];h=add(h,n->type);h=add(h,(uint32_t)n->box.min_x);h=add(h,(uint32_t)n->box.min_z);h=add(h,(uint32_t)n->box.max_x);h=add(h,(uint32_t)n->box.max_z);}return h;}
static uint64_t forerunner_hash(const ForerunnerPlan*p){uint64_t h=add(UINT64_C(1469598103934665603),p->layout_family);for(int i=0;i<p->node_count;i++){const ForerunnerNode*n=&p->nodes[i];h=add(h,n->type);h=add(h,(uint32_t)n->box.min_x);h=add(h,(uint32_t)n->box.min_y);h=add(h,(uint32_t)n->box.min_z);h=add(h,(uint32_t)n->box.max_x);h=add(h,(uint32_t)n->box.max_y);h=add(h,(uint32_t)n->box.max_z);}return h;}
static uint64_t sanctum_hash(const SanctumCitadelPlan*p){uint64_t h=add(UINT64_C(1469598103934665603),p->layout_family);for(int i=0;i<p->node_count;i++){const SanctumNode*n=&p->nodes[i];h=add(h,n->motif);h=add(h,(uint32_t)n->x);h=add(h,(uint32_t)n->y);h=add(h,(uint32_t)n->z);h=add(h,(uint32_t)n->w);h=add(h,(uint32_t)n->h);h=add(h,(uint32_t)n->d);}return h;}
static void check(uint64_t hashes[256],int families[4],const char*label){int unique=0;for(int i=0;i<256;i++){bool first=true;for(int j=0;j<i;j++)if(hashes[i]==hashes[j]){first=false;break;}unique+=first;}for(int i=0;i<4;i++)assert(families[i]>=25);if(unique<230){fprintf(stderr,"%s unique=%d/256\n",label,unique);assert(unique>=230);}}
int main(void){uint64_t hashes[256];int family[4]={0};
    for(uint32_t s=1;s<=256;s++){TemplePlan a=generate_greek_temple(s,TEMPLE_STAGE_SANCTUARY,2),b=generate_greek_temple(s,TEMPLE_STAGE_SANCTUARY,2);assert(memcmp(&a,&b,sizeof(a))==0);family[a.layout_family]++;hashes[s-1]=greek_hash(&a);}check(hashes,family,"greek");
    memset(family,0,sizeof(family));for(uint32_t s=1;s<=256;s++){MegalithPlan p=generate_megalith_structure(s,MEGALITH_ARCHETYPE_PASSAGE_GRAVE,MEGALITH_STAGE_TUMULUS);family[p.layout_family]++;hashes[s-1]=mega_hash(&p);}check(hashes,family,"megalith");
    memset(family,0,sizeof(family));for(uint32_t s=1;s<=256;s++){HyperPlan p=generate_hyperborean_structure(s,HYPER_STAGE_FULL_SANCTUM);family[p.layout_family]++;hashes[s-1]=hyper_hash(&p);}check(hashes,family,"hyperborean");
    memset(family,0,sizeof(family));for(uint32_t s=1;s<=256;s++){ForerunnerPlan p=generate_forerunner_structure(s,FORERUNNER_ARCHETYPE_CARTOGRAPHER,FORERUNNER_STAGE_APEX);family[p.layout_family]++;hashes[s-1]=forerunner_hash(&p);}check(hashes,family,"forerunner");
    memset(family,0,sizeof(family));for(uint32_t s=1;s<=256;s++){SanctumCitadelPlan p=generate_unified_sanctum(s,20);family[p.layout_family]++;hashes[s-1]=sanctum_hash(&p);}check(hashes,family,"citadel");
    puts("building diversity: five styles pass family distribution and 90% fingerprint uniqueness");return 0;}
