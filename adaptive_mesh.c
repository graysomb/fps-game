#include "adaptive_mesh.h"
#include <stdlib.h>
#include <string.h>

_Static_assert(sizeof(AdaptiveLeaf) == 32, "GPU leaf layout");
uint64_t adaptive_morton(const int32_t xyz[3]) {
    uint64_t k = 0;
    for (unsigned a=0;a<3;a++) {
        uint32_t p=(uint32_t)(xyz[a]+ADAPTIVE_COORD_LIMIT);
        for (unsigned b=0;b<21;b++) k|=(uint64_t)((p>>b)&1)<<(3*b+a);
    }
    return k;
}
void adaptive_origin(uint64_t k,int32_t xyz[3]) {
    for (unsigned a=0;a<3;a++) {
        uint32_t p=0;
        for (unsigned b=0;b<21;b++) p|=(uint32_t)((k>>(3*b+a))&1)<<b;
        xyz[a]=(int32_t)p-ADAPTIVE_COORD_LIMIT;
    }
}
static int cmp(const void *aa,const void *bb) {
    const AdaptiveLeaf *a=aa,*b=bb;
    if (a->arena!=b->arena) return a->arena<b->arena?-1:1;
    return a->key<b->key?-1:a->key>b->key;
}
static size_t find(const AdaptiveLeaf *a,size_t n,uint32_t arena,uint64_t key) {
    size_t lo=0,hi=n;
    AdaptiveLeaf needle={.key=key,.arena=arena};
    while(lo<hi) {size_t m=lo+(hi-lo)/2;if(cmp(a+m,&needle)<=0)lo=m+1;else hi=m;}
    if(!lo)return n;
    size_t i=lo-1;
    return a[i].arena==arena && key-a[i].key<(UINT64_C(1)<<(3*a[i].level))?i:n;
}
void adaptive_mesh_release(AdaptiveMesh *m) {if(m){free(m->leaves);*m=(AdaptiveMesh){0};}}
int adaptive_mesh_reference(const AdaptiveInput *in,size_t n,unsigned levels,
                            int surface,int balance,AdaptiveMesh *out) {
    if(!out||(!in&&n)||levels>20||n>UINT32_MAX||n>SIZE_MAX/sizeof(AdaptiveLeaf))return ADAPTIVE_BAD_INPUT;
    if(!n){adaptive_mesh_release(out);return ADAPTIVE_SUCCESS;}
    AdaptiveLeaf *a=calloc(n,sizeof(*a)),*b=calloc(n,sizeof(*b)),*original=calloc(n,sizeof(*original));
    unsigned char *split=calloc(n,1);
    int result=ADAPTIVE_ALLOCATION_FAILED;
    if(!a||!b||!original||!split)goto done;
    for(size_t i=0;i<n;i++) {
        for(unsigned d=0;d<3;d++)if(in[i].xyz[d]<-ADAPTIVE_COORD_LIMIT||in[i].xyz[d]>=ADAPTIVE_COORD_LIMIT){result=ADAPTIVE_BAD_INPUT;goto done;}
        a[i]=(AdaptiveLeaf){adaptive_morton(in[i].xyz),in[i].arena,in[i].material,0,in[i].flags,(uint32_t)i,0};
    }
    qsort(a,n,sizeof(*a),cmp);
    for(size_t i=1;i<n;i++)if(!cmp(a+i-1,a+i)){result=ADAPTIVE_BAD_INPUT;goto done;}
    if(surface)for(size_t i=0;i<n;i++) {
        int32_t p[3];adaptive_origin(a[i].key,p);
        for(unsigned f=0;f<6;f++) {
            int32_t q[3]={p[0],p[1],p[2]};q[f/2]+=(f&1)?1:-1;
            if(q[f/2]<-ADAPTIVE_COORD_LIMIT||q[f/2]>=ADAPTIVE_COORD_LIMIT||find(a,n,a[i].arena,adaptive_morton(q))==n)a[i].flags|=ADAPTIVE_PROTECTED;
        }
    }
    memcpy(original,a,n*sizeof(*a));
    size_t count=n;
    for(unsigned level=0;level<levels;level++) {
        size_t w=0;int merged=0;
        uint64_t step=UINT64_C(1)<<(3*level);
        for(size_t i=0;i<count;) {
            int yes=count-i>=8 && a[i].level==level && !(a[i].key&(8*step-1));
            for(unsigned j=0;yes&&j<8;j++)
                yes=a[i+j].level==level&&a[i+j].arena==a[i].arena&&a[i+j].material==a[i].material&&
                    !(a[i+j].flags&ADAPTIVE_PROTECTED)&&a[i+j].key==a[i].key+j*step;
            b[w]=a[i];if(yes){b[w].level++;i+=8;merged=1;}else i++;w++;
        }
        AdaptiveLeaf *t=a;a=b;b=t;count=w;
        if(!merged)break;
    }
    /* Independent oracle: exhaustive positive-area face adjacency, not the
     * production point-location/request algorithm. O(n^2), test use only. */
    if(balance)for(;;) {
        memset(split,0,count);int any=0;
        for(size_t i=0;i<count;i++)for(size_t j=i+1;j<count;j++) {
            if(a[i].arena!=a[j].arena)continue;
            int32_t p[3],q[3];adaptive_origin(a[i].key,p);adaptive_origin(a[j].key,q);
            int s=1<<a[i].level,t=1<<a[j].level;
            for(unsigned d=0;d<3;d++) {
                unsigned u=(d+1)%3,v=(d+2)%3;
                if((p[d]+s==q[d]||q[d]+t==p[d])&&p[u]<q[u]+t&&q[u]<p[u]+s&&p[v]<q[v]+t&&q[v]<p[v]+s) {
                    if(a[i].level>a[j].level+1){split[i]=1;any=1;}
                    if(a[j].level>a[i].level+1){split[j]=1;any=1;}
                }
            }
        }
        if(!any)break;
        size_t w=0;
        for(size_t i=0;i<count;i++) {
            unsigned parts=split[i]?8:1;
            AdaptiveLeaf c=a[i];if(split[i])c.level--;
            for(unsigned j=0;j<parts;j++){b[w]=c;b[w++].key+=j*(UINT64_C(1)<<(3*c.level));}
        }
        AdaptiveLeaf *t=a;a=b;b=t;count=w;
    }
    for(size_t i=0;i<count;i++){size_t source=find(original,n,a[i].arena,a[i].key);a[i].source=original[source].source;}
    adaptive_mesh_release(out);out->leaves=a;out->count=count;a=NULL;result=ADAPTIVE_SUCCESS;
done:free(a);free(b);free(original);free(split);return result;
}
