#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "adaptive_metal.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned rng=0x784920;
static unsigned random32(void){rng^=rng<<13;rng^=rng>>17;rng^=rng<<5;return rng;}
static void check_coverage(const AdaptiveInput *in,size_t n,const AdaptiveMesh *m) {
    unsigned char *seen=calloc(n?n:1,1);size_t volume=0;
    for(size_t i=0;i<m->count;i++) {
        const AdaptiveLeaf *c=m->leaves+i;int32_t p[3];adaptive_origin(c->key,p);int side=1<<c->level;
        assert(!(c->key&((UINT64_C(1)<<(3*c->level))-1)));
        size_t hits=0;
        for(size_t j=0;j<n;j++)if(in[j].arena==c->arena&&in[j].xyz[0]>=p[0]&&in[j].xyz[0]<p[0]+side&&in[j].xyz[1]>=p[1]&&in[j].xyz[1]<p[1]+side&&in[j].xyz[2]>=p[2]&&in[j].xyz[2]<p[2]+side){assert(in[j].material==c->material);assert(!seen[j]);seen[j]=1;hits++;}
        assert(hits==(size_t)side*side*side);volume+=hits;
    }
    assert(volume==n);for(size_t i=0;i<n;i++)assert(seen[i]);free(seen);
}
static void run(id<MTLDevice> device,id<MTLCommandQueue> queue,AdaptiveInput *in,size_t n,unsigned level,int surface,int balance) {
    AdaptiveMesh oracle={0},actual={0};
    assert(adaptive_mesh_reference(in,n,level,surface,balance,&oracle)==0);check_coverage(in,n,&oracle);
    char error[4096]={0};AdaptiveMetal *gpu=adaptive_metal_create((void *)device,"shaders/adaptive/mesh.metal",n?n:1,error,sizeof(error));
    if(!gpu){fprintf(stderr,"%s\n",error);abort();}
    assert(adaptive_metal_upload(gpu,in,n));
    id<MTLCommandBuffer> command=[queue commandBuffer];id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
    assert(adaptive_metal_encode(gpu,(void *)encoder,level,surface,balance));[encoder endEncoding];[command commit];[command waitUntilCompleted];
    if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"%s\n",command.error.localizedDescription.UTF8String);abort();}
    AdaptiveGpuHeader h;assert(adaptive_metal_read(gpu,&actual,&h));
    if(actual.count!=oracle.count){fprintf(stderr,"count mismatch input=%zu L=%u surface=%d balance=%d gpu=%zu cpu=%zu\n",n,level,surface,balance,actual.count,oracle.count);abort();}
    for(size_t i=0;i<actual.count;i++)if(actual.leaves[i].key!=oracle.leaves[i].key||actual.leaves[i].level!=oracle.leaves[i].level||actual.leaves[i].arena!=oracle.leaves[i].arena||actual.leaves[i].material!=oracle.leaves[i].material||actual.leaves[i].source!=oracle.leaves[i].source){fprintf(stderr,"leaf mismatch %zu input=%zu\n",i,n);abort();}
    check_coverage(in,n,&actual);adaptive_mesh_release(&oracle);adaptive_mesh_release(&actual);adaptive_metal_destroy(gpu);
}
static void reject_duplicate(id<MTLDevice> device,id<MTLCommandQueue> queue,const AdaptiveInput *in) {
    char error[1024]={0};
    AdaptiveMetal *gpu=adaptive_metal_create((void *)device,"shaders/adaptive/mesh.metal",2,error,sizeof(error));
    assert(gpu&&adaptive_metal_upload(gpu,in,2));
    id<MTLCommandBuffer> command=[queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
    assert(adaptive_metal_encode(gpu,(void *)encoder,20,true,true));
    [encoder endEncoding];[command commit];[command waitUntilCompleted];
    assert(command.status==MTLCommandBufferStatusCompleted);
    AdaptiveMesh output={0};AdaptiveGpuHeader header;
    assert(!adaptive_metal_read(gpu,&output,&header)&&header.error&&output.count==0);
    AdaptiveInput bad=*in;bad.xyz[0]=ADAPTIVE_COORD_LIMIT;
    assert(!adaptive_metal_upload(gpu,&bad,1));
    adaptive_metal_destroy(gpu);
}
int main(void){@autoreleasepool{
    id<MTLDevice> d=MTLCreateSystemDefaultDevice();if(!d){fputs("Metal unavailable\n",stderr);return 77;}id<MTLCommandQueue> q=[d newCommandQueue];
    AdaptiveInput in[8192];size_t n=0;
    run(d,q,NULL,0,20,1,1);
    for(int z=0;z<8;z++)for(int y=0;y<8;y++)for(int x=0;x<8;x++)in[n++]=(AdaptiveInput){{x,y,z},0,0,0};
    run(d,q,in,n,20,0,0);run(d,q,in,n,20,1,1);
    /* Tiny off-center face neighbor, forcing multiple balance rounds. */
    in[n++]=(AdaptiveInput){{8,3,5},0,0,0};run(d,q,in,n,20,0,1);
    for(size_t i=0;i<n;i++)in[i].xyz[0]-=16;
    run(d,q,in,n,20,0,1);
    for(size_t i=0;i<n;i++){in[i].material=(uint32_t)(i%5==0);in[i].flags=(i%7==0)?1:0;}
    run(d,q,in,n,20,1,1);
    for(unsigned trial=0;trial<16;trial++){
        n=0;for(int z=-4;z<5;z++)for(int y=-4;y<5;y++)for(int x=-4;x<5;x++)if(random32()%5)in[n++]=(AdaptiveInput){{x,y,z},trial%3,random32()%7==0,random32()%13==0};
        for(size_t i=n;i>1;i--){size_t j=random32()%i;AdaptiveInput tmp=in[i-1];in[i-1]=in[j];in[j]=tmp;}
        run(d,q,in,n,20,trial%2,1);
    }
    /* Hierarchical scan spans more than 128*128 request entries. */
    n=0;for(int z=0;z<16;z++)for(int y=0;y<16;y++)for(int x=0;x<16;x++)in[n++]=(AdaptiveInput){{x,y,z},2,0,0};
    run(d,q,in,n,20,1,1);
    AdaptiveMesh out={0};AdaptiveInput duplicate[2]={{{0,0,0},0,0,0},{{0,0,0},0,1,0}};
    assert(adaptive_mesh_reference(duplicate,2,20,0,0,&out)==ADAPTIVE_BAD_INPUT);
    reject_duplicate(d,q,duplicate);
    duplicate[1].arena=1;run(d,q,duplicate,2,20,1,1);
    duplicate[0].xyz[0]=ADAPTIVE_COORD_LIMIT;assert(adaptive_mesh_reference(duplicate,2,20,0,0,&out)==ADAPTIVE_BAD_INPUT);
    [q release];[d release];puts("adaptive mesh: CPU coverage + Metal oracle parity passed");return 0;
}}
