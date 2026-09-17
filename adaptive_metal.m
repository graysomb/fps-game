#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "adaptive_metal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct AdaptiveMetal {
    id<MTLComputePipelineState> pipeline;
    id<MTLBuffer> leaves[2], requests[2], work, header, fine;
    id<MTLBuffer> fracture_voxels;
    uint32_t capacity, count, arena_bits;
    unsigned leaf_side, request_side;
    size_t bytes;
    uint64_t dispatches;
};
typedef struct {uint32_t mode,count,bit,level,src,dst,aux,surface;} Params;
static void error_text(char *out,size_t n,NSString *s){if(out&&n)snprintf(out,n,"%s",s.UTF8String?:"Metal error");}
void adaptive_metal_destroy(AdaptiveMetal *m) {
    if(!m)return;
    [m->pipeline release];for(int i=0;i<2;i++){[m->leaves[i] release];[m->requests[i] release];}
    [m->fine release];[m->work release];[m->header release];free(m);
}
AdaptiveMetal *adaptive_metal_create(void *device,const char *path,size_t cap,char *error,size_t len) {
    if(!device||!cap||cap>UINT32_MAX/32u){error_text(error,len,@"invalid adaptive capacity");return NULL;}
    AdaptiveMetal *m=calloc(1,sizeof(*m));if(!m)return NULL;m->capacity=(uint32_t)cap;
    id<MTLDevice> d=(id<MTLDevice>)device;
    NSError *err=nil;
    NSString *src=[NSString stringWithContentsOfFile:[NSString stringWithUTF8String:path] encoding:NSUTF8StringEncoding error:&err];
    id<MTLLibrary> lib=src?[d newLibraryWithSource:src options:nil error:&err]:nil;
    id<MTLFunction> fn=[lib newFunctionWithName:@"adaptive_mesh"];
    MTLComputePipelineDescriptor *descriptor=[[MTLComputePipelineDescriptor alloc] init];
    descriptor.computeFunction=fn;descriptor.supportIndirectCommandBuffers=YES;
    m->pipeline=fn?[d newComputePipelineStateWithDescriptor:descriptor options:MTLPipelineOptionNone reflection:nil error:&err]:nil;
    [descriptor release];[fn release];[lib release];
    if(!m->pipeline){error_text(error,len,err.localizedDescription);adaptive_metal_destroy(m);return NULL;}
    if(m->pipeline.maxTotalThreadsPerThreadgroup<128){error_text(error,len,@"adaptive workgroup requires 128 threads");adaptive_metal_destroy(m);return NULL;}
    for(int i=0;i<2;i++){
        m->leaves[i]=[d newBufferWithLength:cap*sizeof(AdaptiveLeaf) options:MTLResourceStorageModeShared];
        m->requests[i]=[d newBufferWithLength:6*cap*sizeof(uint32_t) options:MTLResourceStorageModeShared];
    }
    /* Six requests/cell, scan input/output, and recursive block sums. */
    m->fine=[d newBufferWithLength:cap*sizeof(AdaptiveLeaf) options:MTLResourceStorageModeShared];
    m->work=[d newBufferWithLength:(16*cap+1024)*sizeof(uint32_t) options:MTLResourceStorageModeShared];
    m->header=[d newBufferWithLength:sizeof(AdaptiveGpuHeader) options:MTLResourceStorageModeShared];
    if(!m->leaves[0]||!m->leaves[1]||!m->requests[0]||!m->requests[1]||!m->fine||!m->work||!m->header){adaptive_metal_destroy(m);return NULL;}
    m->bytes=3*cap*sizeof(AdaptiveLeaf)+(28*cap+1024)*sizeof(uint32_t)+sizeof(AdaptiveGpuHeader);
    return m;
}
bool adaptive_metal_upload(AdaptiveMetal *m,const AdaptiveInput *input,size_t count) {
    if(!m||count>m->capacity||(!input&&count))return false;
    for(size_t i=0;i<count;i++)for(int d=0;d<3;d++)if(input[i].xyz[d]<-ADAPTIVE_COORD_LIMIT||input[i].xyz[d]>=ADAPTIVE_COORD_LIMIT)return false;
    AdaptiveLeaf *a=m->leaves[0].contents;uint32_t max_arena=0;
    for(size_t i=0;i<count;i++){
        a[i]=(AdaptiveLeaf){adaptive_morton(input[i].xyz),input[i].arena,input[i].material,0,input[i].flags,(uint32_t)i,0};
        if(input[i].arena>max_arena)max_arena=input[i].arena;
    }
    m->arena_bits=0;while(max_arena){m->arena_bits++;max_arena>>=1;}
    m->count=(uint32_t)count;m->leaf_side=m->request_side=0;
    AdaptiveGpuHeader h={0};h.leaves=h.input_count=(uint32_t)count;h.active=1;
    memcpy(m->header.contents,&h,sizeof(h));return true;
}
static void dispatch(AdaptiveMetal *m,id<MTLComputeCommandEncoder> e,Params p) {
    m->dispatches++;
    [e memoryBarrierWithScope:MTLBarrierScopeBuffers];
    [e setComputePipelineState:m->pipeline];
    [e setBuffer:m->leaves[m->leaf_side] offset:0 atIndex:0];
    [e setBuffer:m->leaves[1-m->leaf_side] offset:0 atIndex:1];
    [e setBuffer:m->requests[m->request_side] offset:0 atIndex:2];
    [e setBuffer:m->requests[1-m->request_side] offset:0 atIndex:3];
    [e setBuffer:m->work offset:0 atIndex:4];[e setBuffer:m->header offset:0 atIndex:5];
    [e setBytes:&p length:sizeof(p) atIndex:6];
    [e setBuffer:m->fracture_voxels?:m->work offset:0 atIndex:7];
    [e setBuffer:m->fine offset:0 atIndex:8];
    [e dispatchThreadgroups:MTLSizeMake((p.count+127)/128,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
}
static uint32_t scan(AdaptiveMetal *m,id<MTLComputeCommandEncoder> e,uint32_t n,uint32_t src,uint32_t dst,uint32_t temp) {
    uint32_t groups=(n+127)/128;
    dispatch(m,e,(Params){.mode=2,.count=n,.src=src,.dst=dst,.aux=temp});
    if(groups>1){uint32_t block_prefix=temp+groups;scan(m,e,groups,temp,block_prefix,block_prefix+groups);dispatch(m,e,(Params){.mode=3,.count=n,.src=block_prefix,.dst=dst});}
    return dst;
}
static void prefix(AdaptiveMetal *m,id<MTLComputeCommandEncoder> e,uint32_t n){scan(m,e,n,0,6*m->count,12*m->count);}
bool adaptive_metal_encode(AdaptiveMetal *m,void *encoder,unsigned levels,bool surface,bool balance) {
    if(!m||!encoder||levels>20)return false;
    id<MTLComputeCommandEncoder> e=(id<MTLComputeCommandEncoder>)encoder;uint32_t n=m->count,off=6*n;
    if(!n)return true;
    for(unsigned bit=0;bit<63+m->arena_bits;bit++){
        dispatch(m,e,(Params){.mode=1,.count=n,.bit=bit});prefix(m,e,n);
        dispatch(m,e,(Params){.mode=4,.count=n,.bit=bit,.dst=off});m->leaf_side^=1;
    }
    dispatch(m,e,(Params){.mode=27,.count=n});
    dispatch(m,e,(Params){.mode=18,.count=n});prefix(m,e,n);dispatch(m,e,(Params){.mode=19,.count=1,.dst=off,.aux=n});
    dispatch(m,e,(Params){.mode=5,.count=n,.surface=surface});dispatch(m,e,(Params){.mode=6,.count=n});
    /* Stop at the occupancy-derived maximum: no larger full cube can exist. */
    unsigned bound=0;for(size_t volume=8;volume<=n&&bound<20;volume*=8)bound++;
    if(levels>bound)levels=bound;
    for(unsigned level=0;level<levels;level++){
        dispatch(m,e,(Params){.mode=7,.count=n,.level=level});prefix(m,e,n);
        dispatch(m,e,(Params){.mode=8,.count=n,.dst=off});
        dispatch(m,e,(Params){.mode=9,.count=1,.dst=off,.aux=n});
    }
    dispatch(m,e,(Params){.mode=20,.count=1});
    if(balance)for(unsigned pass=0;pass<levels;pass++){
        dispatch(m,e,(Params){.mode=10,.count=6*n});prefix(m,e,6*n);
        dispatch(m,e,(Params){.mode=11,.count=6*n,.dst=off});m->request_side^=1;
        unsigned bits=0;for(uint32_t v=n;v;v>>=1)bits++;
        for(unsigned bit=0;bit<bits;bit++){
            dispatch(m,e,(Params){.mode=12,.count=6*n,.bit=bit});prefix(m,e,6*n);
            dispatch(m,e,(Params){.mode=13,.count=6*n,.bit=bit,.dst=off});m->request_side^=1;
        }
        dispatch(m,e,(Params){.mode=14,.count=6*n});prefix(m,e,6*n);
        dispatch(m,e,(Params){.mode=15,.count=6*n,.dst=off});m->request_side^=1;
        dispatch(m,e,(Params){.mode=16,.count=n});prefix(m,e,n);
        dispatch(m,e,(Params){.mode=21,.count=1,.dst=off,.aux=n});
        dispatch(m,e,(Params){.mode=17,.count=n,.dst=off});
        dispatch(m,e,(Params){.mode=9,.count=1,.dst=off,.aux=n});
    }
    dispatch(m,e,(Params){.mode=20,.count=1});
    if(balance){dispatch(m,e,(Params){.mode=25,.count=n});prefix(m,e,n);dispatch(m,e,(Params){.mode=19,.count=1,.dst=off,.aux=n});}
    dispatch(m,e,(Params){.mode=23,.count=n});dispatch(m,e,(Params){.mode=24,.count=1});
    return true;
}
void *adaptive_metal_fine_buffer(AdaptiveMetal *m){return m?(void *)m->fine:NULL;}
bool adaptive_metal_encode_refine(AdaptiveMetal *m,void *encoder,void *voxels,unsigned levels) {
    if(!m||!encoder||!voxels||levels>20)return false;
    m->fracture_voxels=(id<MTLBuffer>)voxels;
    id<MTLComputeCommandEncoder> e=(id<MTLComputeCommandEncoder>)encoder;
    uint32_t n=m->count,off=6*n;if(!n)return true;
    dispatch(m,e,(Params){.mode=20,.count=1});
    for(unsigned level=0;level<levels;level++) {
        dispatch(m,e,(Params){.mode=30,.count=n});prefix(m,e,n);
        dispatch(m,e,(Params){.mode=21,.count=1,.dst=off,.aux=n});
        dispatch(m,e,(Params){.mode=17,.count=n,.dst=off});
        dispatch(m,e,(Params){.mode=9,.count=1,.dst=off,.aux=n});
    }
    /* Exact owner-computes face search: inspect every fine face tile, never
     * sample just four quadrants. Bounded by occupied surface area. */
    dispatch(m,e,(Params){.mode=20,.count=1});
    for(unsigned level=0;level<levels;level++) {
        dispatch(m,e,(Params){.mode=31,.count=n});prefix(m,e,n);
        dispatch(m,e,(Params){.mode=21,.count=1,.dst=off,.aux=n});
        dispatch(m,e,(Params){.mode=17,.count=n,.dst=off});
        dispatch(m,e,(Params){.mode=9,.count=1,.dst=off,.aux=n});
    }
    dispatch(m,e,(Params){.mode=20,.count=1});
    dispatch(m,e,(Params){.mode=25,.count=n});prefix(m,e,n);dispatch(m,e,(Params){.mode=19,.count=1,.dst=off,.aux=n});
    dispatch(m,e,(Params){.mode=23,.count=n});dispatch(m,e,(Params){.mode=24,.count=1});
    return true;
}
void *adaptive_metal_leaf_buffer(AdaptiveMetal *m){return m?(void *)m->leaves[m->leaf_side]:NULL;}
void *adaptive_metal_header_buffer(AdaptiveMetal *m){return m?(void *)m->header:NULL;}
uint64_t adaptive_metal_encoded_dispatches(const AdaptiveMetal *m){return m?m->dispatches:0;}
size_t adaptive_metal_scratch_bytes(const AdaptiveMetal *m){return m?m->bytes:0;}
bool adaptive_metal_read(const AdaptiveMetal *m,AdaptiveMesh *out,AdaptiveGpuHeader *h) {
    if(!m||!out||!h)return false;memcpy(h,m->header.contents,sizeof(*h));
    if(h->error||h->leaves>m->count)return false;
    AdaptiveLeaf *a=h->leaves?malloc(h->leaves*sizeof(*a)):NULL;if(h->leaves&&!a)return false;
    if(h->leaves)memcpy(a,m->leaves[m->leaf_side].contents,h->leaves*sizeof(*a));
    adaptive_mesh_release(out);out->leaves=a;out->count=h->leaves;return true;
}
