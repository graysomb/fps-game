/* A tiny recorder for the subset of MTLComputeCommandEncoder used by the
 * adaptive kernels. It caches the immutable refinement/rebuild schedule in a
 * Metal ICB; the GPU chooses a zero-length range when topology is unchanged. */
@interface AdaptiveCommandRecorder : NSObject {
@public
    id<MTLIndirectCommandBuffer> commands;
    id<MTLBuffer> constants;
    NSMutableSet *resources;
    NSUInteger count, constant_bytes;
    BOOL failed;
@private
    id<MTLBuffer> bound[23];
    NSUInteger offsets[23];
    id<MTLComputePipelineState> pipeline;
}
- (id)initWithDevice:(id<MTLDevice>)device;
- (void)memoryBarrierWithScope:(MTLBarrierScope)scope;
- (void)setComputePipelineState:(id<MTLComputePipelineState>)state;
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index;
- (void)setBytes:(const void *)bytes length:(NSUInteger)length atIndex:(NSUInteger)index;
- (void)dispatchThreadgroups:(MTLSize)groups threadsPerThreadgroup:(MTLSize)threads;
- (void)makeResident:(id<MTLComputeCommandEncoder>)encoder;
@end
@implementation AdaptiveCommandRecorder
- (id)initWithDevice:(id<MTLDevice>)device {
    if((self=[super init])) {
        MTLIndirectCommandBufferDescriptor *d=[[MTLIndirectCommandBufferDescriptor alloc] init];
        d.commandTypes=MTLIndirectCommandTypeConcurrentDispatch;d.inheritPipelineState=NO;d.inheritBuffers=NO;d.maxKernelBufferBindCount=23;
        commands=[device newIndirectCommandBufferWithDescriptor:d maxCommandCount:4096 options:MTLResourceStorageModeShared];[d release];
        constants=[device newBufferWithLength:4096*1024 options:MTLResourceStorageModeShared];
        resources=[[NSMutableSet alloc] init];failed=!commands||!constants;
    }return self;
}
- (void)dealloc {[commands release];[constants release];[resources release];[super dealloc];}
- (void)memoryBarrierWithScope:(MTLBarrierScope)scope {(void)scope;}
- (void)setComputePipelineState:(id<MTLComputePipelineState>)state {pipeline=state;}
- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSUInteger)offset atIndex:(NSUInteger)index {
    if(index>=23){failed=YES;return;}bound[index]=buffer;offsets[index]=offset;if(buffer)[resources addObject:buffer];
}
- (void)setBytes:(const void *)bytes length:(NSUInteger)length atIndex:(NSUInteger)index {
    NSUInteger aligned=(length+255)&~(NSUInteger)255;
    if(constant_bytes+aligned>constants.length){failed=YES;return;}
    memcpy((char *)constants.contents+constant_bytes,bytes,length);
    [self setBuffer:constants offset:constant_bytes atIndex:index];constant_bytes+=aligned;
}
- (void)dispatchThreadgroups:(MTLSize)groups threadsPerThreadgroup:(MTLSize)threads {
    if(failed||count>=4096){failed=YES;return;}
    id<MTLIndirectComputeCommand> c=[commands indirectComputeCommandAtIndex:count++];
    [c setComputePipelineState:pipeline];
    for(NSUInteger i=0;i<23;i++)if(bound[i])[c setKernelBuffer:bound[i] offset:offsets[i] atIndex:i];
    [c concurrentDispatchThreadgroups:groups threadsPerThreadgroup:threads];[c setBarrier];
}
- (void)makeResident:(id<MTLComputeCommandEncoder>)encoder {
    for(id<MTLResource> resource in resources)[encoder useResource:resource usage:MTLResourceUsageRead|MTLResourceUsageWrite];
    [encoder useResource:commands usage:MTLResourceUsageRead];
}
@end
