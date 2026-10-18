#include <assert.h>
#include <stdio.h>
#include "../firefight_rules.h"
int main(void) {
    assert(ff_count(1,1,0)==0 && ff_count(1,1,1)==3);
    assert(ff_count(2,4,1)==10 && ff_count(2,1,2)==1);
    assert(ff_count(3,1,5)==1 && ff_count(4,1,7)==1 && ff_count(6,1,8)==1);
    for(int h=1;h<=4;h++)for(int wave=1;wave<=10;wave++)
        assert(ff_count(wave,h,0)==(wave==2 || wave==4 || wave==5 || wave>=7));
    for(int wave=1;wave<=50;wave++) {
        assert(ff_ai_level(wave,0)==(wave<5?0:wave<10?1:2));
        for(int type=1;type<FF_TYPE_COUNT;type++)
            assert(ff_ai_level(wave,type)==(wave==1?0:wave<=3?1:2));
    }
    assert(ff_count(13,1,1)==4 && ff_count(11,4,2)==5);
    assert(ff_count(12,4,1)==6 && ff_count(12,4,9)==1);
    assert(ff_wave_cap(1,1)==2 && ff_wave_cap(2,1)==3);
    assert(ff_cap(1)==6 && ff_cap(4)==12);
    int votes[4]={-1,-1,-1,-1}; assert(ff_vote(votes,4)==UPGRADE_BUILD);
    votes[0]=1; votes[1]=2; assert(ff_vote(votes,2)==UPGRADE_BUILD);
    votes[2]=1; assert(ff_vote(votes,3)==UPGRADE_HARVEST);
    for(int h=1;h<=4;h++) for(int w=1;w<=50;w++) for(int t=0;t<FF_TYPE_COUNT;t++) {
        assert(ff_count(w,h,t)>=0);
        assert(ff_count(w+FF_WAVE_COUNT,h,t)>=ff_count(w,h,t));
    }
    puts("Firefight rules passed");
}
