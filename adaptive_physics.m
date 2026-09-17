#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "adaptive_physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adaptive_commands.inc"
struct AdaptivePhysics {
    AdaptiveCommandRecorder *refinement;
    id<MTLDevice> device;
    AdaptiveMetal *mesh;
    id<MTLComputePipelineState> pipeline,momentum_pipeline;
    id<MTLBuffer> info,controls,shapes,refs[2],delta,work,owner,arena_ids,arenas,momentum,slot_order,ref_ranges,fine_arena;
    uint32_t n,p,extent,levels,arena_count;
    unsigned side;
    size_t bytes;
    uint64_t dispatches;
};
typedef struct {uint32_t arena,id;} ArenaPair;
static int arena_compare(const void *aa,const void *bb){const ArenaPair *a=aa,*b=bb;return a->arena<b->arena?-1:a->arena>b->arena?1:a->id<b->id?-1:a->id>b->id;}
typedef struct {uint32_t mode,count,n,p,extent,bit,src,dst,aux;} PParams;
void adaptive_physics_destroy(AdaptivePhysics *s){if(!s)return;[s->refinement release];adaptive_metal_destroy(s->mesh);[s->pipeline release];[s->momentum_pipeline release];[s->info release];[s->controls release];[s->shapes release];[s->refs[0] release];[s->refs[1] release];[s->delta release];[s->work release];[s->owner release];[s->arena_ids release];[s->arenas release];[s->momentum release];[s->slot_order release];[s->ref_ranges release];[s->fine_arena release];free(s);}
AdaptivePhysics *adaptive_physics_create(void *device,const AdaptiveInput *in,uint32_t n,const AdaptiveParticleInput *pi,uint32_t p,char *error,size_t error_size){
    if(!device||!n||!p||n>UINT32_MAX/64||p>UINT32_MAX/64)return NULL;
    for(uint32_t i=0;i<p;i++) {
        if(pi[i].count>8)return NULL;
        for(uint32_t j=0;j<pi[i].count;j++) {
            uint32_t cell=pi[i].incident[j]/8,corner=pi[i].incident[j]%8;
            if(cell>=n||pi[i].arena!=in[cell].arena)return NULL;
            for(unsigned axis=0;axis<3;axis++)if(pi[i].xyz[axis]!=in[cell].xyz[axis]+(int)((corner>>axis)&1u))return NULL;
        }
    }
    ArenaPair *pairs=NULL;
    AdaptivePhysics *s=calloc(1,sizeof(*s));if(!s)return NULL;s->n=n;s->p=p;s->extent=8*n+9*p;
    if(s->extent>(UINT32_MAX-1024u)/3u){free(s);return NULL;}
    for(size_t volume=8;volume<=n&&s->levels<20;volume*=8)s->levels++;
    s->mesh=adaptive_metal_create(device,"shaders/adaptive/mesh.metal",n,error,error_size);if(!s->mesh)goto fail;
    if(!adaptive_metal_upload(s->mesh,in,n))goto fail;
    id<MTLDevice> d=(id<MTLDevice>)device;s->device=d;NSError *err=nil;
    NSMutableString *source=[NSMutableString string];
    for(NSString *path in @[@"shaders/pbd/pbd_pipeline.metal",@"shaders/adaptive/mesh.metal",@"shaders/adaptive/physics.metal"]){NSString *part=[NSString stringWithContentsOfFile:path encoding:NSUTF8StringEncoding error:&err];if(!part)goto shader_fail;[source appendString:part];[source appendString:@"\n"];}
    id<MTLLibrary> lib=[d newLibraryWithSource:source options:nil error:&err];
    id<MTLFunction> fn=[lib newFunctionWithName:@"adaptive_physics"];
    MTLComputePipelineDescriptor *descriptor=[[MTLComputePipelineDescriptor alloc] init];descriptor.computeFunction=fn;descriptor.supportIndirectCommandBuffers=YES;
    s->pipeline=fn?[d newComputePipelineStateWithDescriptor:descriptor options:MTLPipelineOptionNone reflection:nil error:&err]:nil;[descriptor release];[fn release];
    fn=[lib newFunctionWithName:@"adaptive_momentum"];
    descriptor=[[MTLComputePipelineDescriptor alloc] init];descriptor.computeFunction=fn;descriptor.supportIndirectCommandBuffers=YES;
    s->momentum_pipeline=fn?[d newComputePipelineStateWithDescriptor:descriptor options:MTLPipelineOptionNone reflection:nil error:&err]:nil;
    [descriptor release];[fn release];[lib release];
    if(!s->pipeline||!s->momentum_pipeline)goto shader_fail;
    if(s->pipeline.maxTotalThreadsPerThreadgroup<128||s->momentum_pipeline.maxTotalThreadsPerThreadgroup<128)goto fail;
#define ALLOC(member,nbytes) s->member=[d newBufferWithLength:(nbytes) options:MTLResourceStorageModeShared];if(!s->member)goto fail;s->bytes+=(nbytes)
    ALLOC(info,(size_t)p*sizeof(*pi));memcpy(s->info.contents,pi,(size_t)p*sizeof(*pi));
    ALLOC(controls,(size_t)p*80);ALLOC(shapes,(size_t)n*64);
    ALLOC(refs[0],(size_t)s->extent*8);ALLOC(refs[1],(size_t)s->extent*8);
    ALLOC(delta,(size_t)s->extent*16);ALLOC(work,((size_t)s->extent*3+1024)*4);ALLOC(owner,(size_t)n*4);
    pairs=malloc((size_t)p*sizeof(*pairs));if(!pairs)goto fail;
    for(uint32_t i=0;i<p;i++)pairs[i]=(ArenaPair){pi[i].arena,i};qsort(pairs,p,sizeof(*pairs),arena_compare);
    for(uint32_t i=0;i<p;i++)if(!i||pairs[i].arena!=pairs[i-1].arena)s->arena_count++;
    ALLOC(arena_ids,(size_t)p*4);ALLOC(arenas,(size_t)s->arena_count*8);ALLOC(momentum,(size_t)s->arena_count*128);
    uint32_t *ids=s->arena_ids.contents,*ranges=s->arenas.contents;AdaptiveParticleInput *mapped=s->info.contents;
    uint32_t arena=0;
    for(uint32_t i=0;i<p;i++){if(i&&pairs[i].arena!=pairs[i-1].arena)arena++;if(!i||pairs[i].arena!=pairs[i-1].arena){ranges[2*arena]=i;ranges[2*arena+1]=0;}ranges[2*arena+1]++;ids[i]=pairs[i].id;mapped[pairs[i].id].arena=arena;}
    ALLOC(slot_order,(size_t)s->extent*4);ALLOC(ref_ranges,(size_t)s->arena_count*16);ALLOC(fine_arena,(size_t)n*4);
    uint32_t *ref_ranges=s->ref_ranges.contents,*fine_arena=s->fine_arena.contents,*slot_order=s->slot_order.contents;
    arena=0;for(uint32_t i=0;i<p;i++)if(!i||pairs[i].arena!=pairs[i-1].arena){ref_ranges[4*arena+2]=pairs[i].arena;arena++;}
    for(uint32_t i=0;i<n;i++){uint32_t lo=0,hi=s->arena_count;while(lo<hi){uint32_t mid=lo+(hi-lo)/2;if(ref_ranges[4*mid+2]<in[i].arena)lo=mid+1;else hi=mid;}if(lo==s->arena_count||ref_ranges[4*lo+2]!=in[i].arena)goto fail;fine_arena[i]=lo;}
    free(pairs);pairs=malloc((size_t)s->extent*sizeof(*pairs));if(!pairs)goto fail;
    for(uint32_t slot=0;slot<s->extent;slot++){uint32_t group=slot<8*n?fine_arena[slot/8]:mapped[(slot-8*n)/9].arena;pairs[slot]=(ArenaPair){group,slot};}
    qsort(pairs,s->extent,sizeof(*pairs),arena_compare);
    for(uint32_t i=0;i<s->extent;i++){uint32_t group=pairs[i].arena;if(!i||pairs[i-1].arena!=group){ref_ranges[4*group]=i;ref_ranges[4*group+1]=0;}ref_ranges[4*group+1]++;slot_order[i]=pairs[i].id;}
    free(pairs);pairs=NULL;
#undef ALLOC
    s->bytes+=adaptive_metal_scratch_bytes(s->mesh);return s;
shader_fail: if(error&&error_size)snprintf(error,error_size,"%s",err.localizedDescription.UTF8String?:"adaptive shader error");
fail:free(pairs);adaptive_physics_destroy(s);return NULL;
}
static void dispatch(AdaptivePhysics *s,id<MTLComputeCommandEncoder> e,PParams params){
    s->dispatches++;
    params.n=s->n;params.p=s->p;params.extent=s->extent;
    [e memoryBarrierWithScope:MTLBarrierScopeBuffers];[e setComputePipelineState:(params.mode==40||params.mode==41)?s->momentum_pipeline:s->pipeline];
    [e setBuffer:s->info offset:0 atIndex:2];[e setBuffer:s->controls offset:0 atIndex:3];
    [e setBuffer:(id<MTLBuffer>)adaptive_metal_leaf_buffer(s->mesh) offset:0 atIndex:4];
    [e setBuffer:(id<MTLBuffer>)adaptive_metal_fine_buffer(s->mesh) offset:0 atIndex:5];
    [e setBuffer:(id<MTLBuffer>)adaptive_metal_header_buffer(s->mesh) offset:0 atIndex:6];
    [e setBuffer:s->shapes offset:0 atIndex:7];[e setBuffer:s->refs[s->side] offset:0 atIndex:8];[e setBuffer:s->refs[1-s->side] offset:0 atIndex:9];
    [e setBuffer:s->delta offset:0 atIndex:10];[e setBuffer:s->work offset:0 atIndex:11];[e setBuffer:s->owner offset:0 atIndex:12];
    [e setBytes:&params length:sizeof(params) atIndex:13];
    [e setBuffer:s->arena_ids offset:0 atIndex:17];[e setBuffer:s->arenas offset:0 atIndex:18];[e setBuffer:s->momentum offset:0 atIndex:19];
    [e setBuffer:s->slot_order offset:0 atIndex:20];[e setBuffer:s->ref_ranges offset:0 atIndex:21];[e setBuffer:s->fine_arena offset:0 atIndex:22];
    [e dispatchThreadgroups:MTLSizeMake((params.count+127)/128,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
}
static void scan(AdaptivePhysics *s,id<MTLComputeCommandEncoder> e,uint32_t n,uint32_t src,uint32_t dst,uint32_t temp){
    uint32_t groups=(n+127)/128;dispatch(s,e,(PParams){.mode=3,.count=n,.src=src,.dst=dst,.aux=temp});
    if(groups>1){uint32_t bp=temp+groups;scan(s,e,groups,temp,bp,bp+groups);dispatch(s,e,(PParams){.mode=4,.count=n,.src=bp,.dst=dst});}
}
static void bind(id<MTLComputeCommandEncoder> e,void *particle,void *voxel,void *sim,void *control,const FpsGpuUniforms *u){
    [e setBuffer:(id<MTLBuffer>)particle offset:0 atIndex:0];[e setBuffer:(id<MTLBuffer>)voxel offset:0 atIndex:1];
    [e setBytes:u length:sizeof(*u) atIndex:14];[e setBuffer:(id<MTLBuffer>)sim offset:0 atIndex:15];[e setBuffer:(id<MTLBuffer>)control offset:0 atIndex:16];
}
static void rebuild(AdaptivePhysics *s,id<MTLComputeCommandEncoder> e){
    dispatch(s,e,(PParams){.mode=50,.count=s->n});
    dispatch(s,e,(PParams){.mode=0,.count=s->n});dispatch(s,e,(PParams){.mode=1,.count=s->n});dispatch(s,e,(PParams){.mode=2,.count=s->p});
    s->side=0;
    dispatch(s,e,(PParams){.mode=5,.count=s->extent});
    uint32_t bits=0;for(uint32_t value=s->p;value;value>>=1)bits++;
    for(uint32_t bit=0;bit<bits;bit++){
        dispatch(s,e,(PParams){.mode=6,.count=s->extent,.bit=bit});scan(s,e,s->extent,0,s->extent,2*s->extent);
        dispatch(s,e,(PParams){.mode=7,.count=s->extent,.bit=bit,.dst=s->extent});s->side^=1;
    }
    dispatch(s,e,(PParams){.mode=12,.count=s->extent});s->side=0;
    dispatch(s,e,(PParams){.mode=8,.count=s->p});
    dispatch(s,e,(PParams){.mode=40,.count=s->arena_count*128});
    dispatch(s,e,(PParams){.mode=9,.count=s->p});
    dispatch(s,e,(PParams){.mode=41,.count=s->arena_count*128});
    dispatch(s,e,(PParams){.mode=42,.count=s->p});
    dispatch(s,e,(PParams){.mode=13,.count=s->p});
    scan(s,e,s->p,0,s->extent,2*s->extent);dispatch(s,e,(PParams){.mode=10,.count=s->p,.dst=s->extent});dispatch(s,e,(PParams){.mode=11,.count=1});
}
bool adaptive_physics_build(AdaptivePhysics *s,void *encoder,void *particle,void *voxel,void *sim,void *control,const FpsGpuUniforms *u){
    if(!s||!encoder)return false;id<MTLComputeCommandEncoder> e=(id<MTLComputeCommandEncoder>)encoder;
    if(!adaptive_metal_encode(s->mesh,encoder,s->levels,true,true))return false;
    bind(e,particle,voxel,sim,control,u);dispatch(s,e,(PParams){.mode=49,.count=s->arena_count});rebuild(s,e);
    s->refinement=[[AdaptiveCommandRecorder alloc] initWithDevice:s->device];
    if(s->refinement->failed)return false;
    if(!adaptive_metal_encode_refine(s->mesh,(void *)s->refinement,voxel,s->levels))return false;
    id<MTLComputeCommandEncoder> recorder=(id<MTLComputeCommandEncoder>)s->refinement;
    bind(recorder,particle,voxel,sim,control,u);rebuild(s,recorder);
    return !s->refinement->failed;
}
bool adaptive_physics_stage(AdaptivePhysics *s,void *encoder,void *particle,void *voxel,void *sim,void *control,const FpsGpuUniforms *u,uint32_t mode){
    if(!s||!encoder)return false;id<MTLComputeCommandEncoder> e=(id<MTLComputeCommandEncoder>)encoder;
    if(mode==ADAPTIVE_STAGE_FRACTURE){
        bind(e,particle,voxel,sim,control,u);
        dispatch(s,e,(PParams){.mode=45,.count=s->n});
        /* The scan must run for detection even when no rebuild is pending. */
        dispatch(s,e,(PParams){.mode=47,.count=1});
        scan(s,e,s->n,0,s->extent,2*s->extent);
        dispatch(s,e,(PParams){.mode=48,.count=s->arena_count,.dst=s->extent});
        dispatch(s,e,(PParams){.mode=46,.count=1,.dst=s->extent,.bit=(uint32_t)s->refinement->count});
        [e memoryBarrierWithScope:MTLBarrierScopeBuffers];[s->refinement makeResident:e];
        [e executeCommandsInBuffer:s->refinement->commands indirectBuffer:(id<MTLBuffer>)adaptive_metal_header_buffer(s->mesh) indirectBufferOffset:16*sizeof(uint32_t)];
        [e memoryBarrierWithScope:MTLBarrierScopeBuffers];return true;
    }
    bind(e,particle,voxel,sim,control,u);
    dispatch(s,e,(PParams){.mode=mode,.count=mode==ADAPTIVE_STAGE_SHAPE?s->n:s->p});
    if(mode==ADAPTIVE_STAGE_SHAPE||mode==ADAPTIVE_STAGE_ATTACH)dispatch(s,e,(PParams){.mode=mode==ADAPTIVE_STAGE_SHAPE?34:35,.count=s->p});
    return true;
}
size_t adaptive_physics_bytes(const AdaptivePhysics *s){return s?s->bytes+(s->refinement?s->refinement->constants.length+s->refinement->commands.allocatedSize:0):0;}
void adaptive_physics_diagnostics(const AdaptivePhysics *s,uint32_t *leaves,uint32_t *generation,uint32_t *error){
    if(!s)return;id<MTLBuffer> b=(id<MTLBuffer>)adaptive_metal_header_buffer(s->mesh);AdaptiveGpuHeader h;memcpy(&h,b.contents,sizeof(h));
    if(leaves)*leaves=h.leaves;if(generation)*generation=h.generation;if(error)*error=h.error;
}

bool adaptive_physics_independent(const AdaptivePhysics *s,uint32_t i){
    if(!s||i>=s->p)return false;const uint32_t *c=s->controls.contents;return c[i*20+16]!=0;
}
double adaptive_physics_mass_error(const AdaptivePhysics *s,void *particle_buffer){
    if(!s)return 0;id<MTLBuffer> b=(id<MTLBuffer>)particle_buffer;
    const float *p=b.contents;const AdaptiveParticleInput *in=s->info.contents;
    double before=0,after=0;
    for(uint32_t i=0;i<s->p;i++){if(in[i].inverse_mass>0)before+=1.0/in[i].inverse_mass;if(p[16*i+7]>0)after+=1.0/p[16*i+7];}
    return before>0?fabs(after-before)/before:0;
}

double adaptive_physics_interface_error(const AdaptivePhysics *s,void *particle_buffer){
    if(!s)return 0;id<MTLBuffer> b=(id<MTLBuffer>)particle_buffer;
    const float *p=b.contents;const uint32_t *ctrl=s->controls.contents;double maximum=0;
    for(uint32_t i=0;i<s->p;i++)if(ctrl[i*20+17]) {
        const float *weights=(const float *)(ctrl+i*20+8);double target[3]={0};
        for(int k=0;k<8;k++)for(int d=0;d<3;d++)target[d]+=weights[k]*p[16*ctrl[i*20+k]+8+d];
        double square=0;for(int d=0;d<3;d++){double diff=target[d]-p[16*i+8+d];square+=diff*diff;}
        double distance=sqrt(square);if(!isfinite(distance))return INFINITY;if(distance>maximum)maximum=distance;
    }
    return maximum;
}

uint64_t adaptive_physics_dispatches(const AdaptivePhysics *s,bool gpu){
    if(!s)return 0;
    uint64_t count=s->dispatches+adaptive_metal_encoded_dispatches(s->mesh)-(s->refinement?s->refinement->count:0);
    if(gpu){id<MTLBuffer> b=(id<MTLBuffer>)adaptive_metal_header_buffer(s->mesh);const uint32_t *h=b.contents;count+=h[10];}
    return count;
}
