#if defined(__APPLE__) && defined(FPS_GPU_METAL)

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "physics_gpu_metal.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct FpsMetalState {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLLibrary> library;
    id<MTLComputePipelineState> pipeline;
    id<MTLCommandBuffer> command_buffer;
    id<MTLComputeCommandEncoder> encoder;
    id<MTLBuffer> bound[FPS_GPU_BUFFER_COUNT];
    FpsGpuUniforms uniforms;
    MTLResourceOptions resource_options;
    double wait_start;
    bool managed;
    bool stage_encoder_mode;
    id<MTLCounterSampleBuffer> stage_samples;
    id<MTLBuffer> stage_resolved;
    unsigned stage_pairs, stage_open[FPS_METAL_STAGES];
    unsigned stage_kind[1024];
    bool stage_complete[1024];
} FpsMetalState;

static FpsMetalState metal_state;

static void fps_metal_error(char *out, size_t capacity, NSString *message) {
    if (!out || capacity == 0) return;
    const char *text = message ? [message UTF8String] : "unknown Metal error";
    snprintf(out, capacity, "%s", text ? text : "unknown Metal error");
}

bool fps_metal_initialize(const char *library_path, long long *max_buffer_size,
                          char *error, size_t error_capacity) {
    @autoreleasepool {
        memset(&metal_state, 0, sizeof(metal_state));
        metal_state.device = MTLCreateSystemDefaultDevice();
        if (!metal_state.device) {
            fps_metal_error(error, error_capacity, @"no Metal device is available");
            return false;
        }
        metal_state.queue = [metal_state.device newCommandQueue];
        if (!metal_state.queue) {
            fps_metal_error(error, error_capacity, @"unable to create the Metal command queue");
            fps_metal_shutdown();
            return false;
        }

#if defined(__arm64__)
        metal_state.managed = false;
        metal_state.resource_options = MTLResourceStorageModeShared;
#else
        metal_state.managed = true;
        metal_state.resource_options = MTLResourceStorageModeManaged;
#endif

        NSString *path = [NSString stringWithUTF8String:library_path ? library_path : ""];
        NSError *library_error = nil;
        NSString *source_path = [path stringByReplacingOccurrencesOfString:@".metallib" withString:@".metal"];
        if (![[NSFileManager defaultManager] fileExistsAtPath:source_path]) {
            source_path = @"shaders/pbd/pbd_pipeline.metal";
        }
        if ([[NSFileManager defaultManager] fileExistsAtPath:source_path]) {
            NSString *source = [NSString stringWithContentsOfFile:source_path encoding:NSUTF8StringEncoding error:nil];
            if (source) {
                MTLCompileOptions *opts = [[MTLCompileOptions alloc] init];
                metal_state.library = [metal_state.device newLibraryWithSource:source options:opts error:&library_error];
                [opts release];
            }
        }
        if (!metal_state.library) {
            metal_state.library = [metal_state.device newLibraryWithFile:path error:&library_error];
        }
        if (!metal_state.library) {
            fps_metal_error(error, error_capacity, [library_error localizedDescription]);
            fps_metal_shutdown();
            return false;
        }
        id<MTLFunction> function = [metal_state.library newFunctionWithName:@"pbd_pipeline"];
        if (!function) {
            fps_metal_error(error, error_capacity, @"pbd_pipeline is missing from the Metal library");
            fps_metal_shutdown();
            return false;
        }
        NSError *pipeline_error = nil;
        metal_state.pipeline = [metal_state.device newComputePipelineStateWithFunction:function
                                                                                   error:&pipeline_error];
        [function release];
        if (!metal_state.pipeline) {
            fps_metal_error(error, error_capacity, [pipeline_error localizedDescription]);
            fps_metal_shutdown();
            return false;
        }
        if (metal_state.pipeline.maxTotalThreadsPerThreadgroup < 128) {
            fps_metal_error(error, error_capacity,
                            @"Metal device does not support the solver's 128-thread workgroup");
            fps_metal_shutdown();
            return false;
        }
        if (max_buffer_size) {
            if ([metal_state.device respondsToSelector:@selector(maxBufferLength)])
                *max_buffer_size = (long long)metal_state.device.maxBufferLength;
            else
                *max_buffer_size = 256ll * 1024ll * 1024ll;
        }
        if(getenv("FPS_ADAPTIVE_PROFILE") &&
           ([metal_state.device supportsCounterSampling:MTLCounterSamplingPointAtDispatchBoundary] ||
            [metal_state.device supportsCounterSampling:MTLCounterSamplingPointAtStageBoundary])) {
            metal_state.stage_encoder_mode=![metal_state.device supportsCounterSampling:MTLCounterSamplingPointAtDispatchBoundary];
            for(id<MTLCounterSet> set in metal_state.device.counterSets)if([set.name isEqualToString:MTLCommonCounterSetTimestamp]) {
                MTLCounterSampleBufferDescriptor *descriptor=[[MTLCounterSampleBufferDescriptor alloc] init];
                descriptor.counterSet=set;descriptor.storageMode=MTLStorageModeShared;descriptor.sampleCount=2048;
                metal_state.stage_samples=[metal_state.device newCounterSampleBufferWithDescriptor:descriptor error:nil];[descriptor release];
                metal_state.stage_resolved=[metal_state.device newBufferWithLength:2048*sizeof(MTLCounterResultTimestamp) options:MTLResourceStorageModeShared];
                break;
            }
        }
        return true;
    }
}

void fps_metal_shutdown(void) {
    @autoreleasepool {
        if (metal_state.encoder) {
            [metal_state.encoder endEncoding];
            [metal_state.encoder release];
            metal_state.encoder = nil;
        }
        if (metal_state.command_buffer) {
            [metal_state.command_buffer waitUntilCompleted];
            [metal_state.command_buffer release];
            metal_state.command_buffer = nil;
        }
        [metal_state.stage_samples release];[metal_state.stage_resolved release];
        [metal_state.pipeline release];
        [metal_state.library release];
        [metal_state.queue release];
        [metal_state.device release];
        memset(&metal_state, 0, sizeof(metal_state));
    }
}

void *fps_metal_buffer_create(size_t size) {
    @autoreleasepool {
        if (!metal_state.device || size == 0) return NULL;
        id<MTLBuffer> buffer = [metal_state.device newBufferWithLength:size
                                                               options:metal_state.resource_options];
        return (void *)buffer;
    }
}

void fps_metal_buffer_destroy(void *handle) {
    @autoreleasepool {
        id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
        if (!buffer) return;
        for (int i = 0; i < FPS_GPU_BUFFER_COUNT; ++i) {
            if (metal_state.bound[i] == buffer) metal_state.bound[i] = nil;
        }
        [buffer release];
    }
}

bool fps_metal_buffer_update(void *handle, const void *data, size_t size, size_t offset) {
    id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
    if (!buffer || !data || offset > buffer.length || size > buffer.length - offset) return false;
    memcpy((uint8_t *)buffer.contents + offset, data, size);
    if (metal_state.managed) {
        [buffer didModifyRange:NSMakeRange(offset, size)];
    }
    return true;
}

bool fps_metal_buffer_read(void *handle, void *data, size_t size, size_t offset) {
    id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
    if (!buffer || !data || offset > buffer.length || size > buffer.length - offset) return false;
    memcpy(data, (const uint8_t *)buffer.contents + offset, size);
    return true;
}

void fps_metal_bind_buffer(int slot, void *handle) {
    if (slot < 0 || slot >= FPS_GPU_BUFFER_COUNT) return;
    metal_state.bound[slot] = (id<MTLBuffer>)handle;
    if (metal_state.encoder && metal_state.bound[slot]) {
        [metal_state.encoder setBuffer:metal_state.bound[slot] offset:0 atIndex:(NSUInteger)slot];
    }
}

void fps_metal_get_uniforms(FpsGpuUniforms *uniforms) { if(uniforms)*uniforms=metal_state.uniforms; }

void fps_metal_set_uniforms(const FpsGpuUniforms *uniforms) {
    if (uniforms) metal_state.uniforms = *uniforms;
}

void fps_metal_set_vgs_color(int color) {
    metal_state.uniforms.vgs_color = color;
}

static FpsMetalProfileInfo metal_profile = { 0 };

static inline double metal_time_now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1000000.0;
}

void fps_metal_get_profile_info(FpsMetalProfileInfo *out) {
    if (out) *out = metal_profile;
}

bool fps_metal_begin_batch(void) {
    @autoreleasepool {
        if (!metal_state.queue || metal_state.command_buffer || metal_state.encoder) return false;
        metal_profile.last_encode_ms = 0.0;
        metal_profile.last_dispatch_count = 0;
        metal_profile.stage_counters_available=metal_state.stage_samples&&metal_state.stage_resolved;
        memset(metal_profile.stage_ticks,0,sizeof(metal_profile.stage_ticks));
        metal_state.stage_pairs=0;memset(metal_state.stage_complete,0,sizeof(metal_state.stage_complete));
        for(unsigned i=0;i<FPS_METAL_STAGES;i++)metal_state.stage_open[i]=UINT32_MAX;
        metal_state.command_buffer = [[metal_state.queue commandBuffer] retain];
        if (!metal_state.command_buffer) return false;
        metal_state.encoder = [[metal_state.command_buffer computeCommandEncoder] retain];
        if (!metal_state.encoder) {
            [metal_state.command_buffer release];
            metal_state.command_buffer = nil;
            return false;
        }
        [metal_state.encoder setComputePipelineState:metal_state.pipeline];
        for (int i = 0; i < FPS_GPU_BUFFER_COUNT; ++i) {
            if (!metal_state.bound[i]) {
                [metal_state.encoder endEncoding];
                [metal_state.encoder release];
                metal_state.encoder = nil;
                [metal_state.command_buffer release];
                metal_state.command_buffer = nil;
                return false;
            }
            [metal_state.encoder setBuffer:metal_state.bound[i] offset:0 atIndex:(NSUInteger)i];
        }
        return true;
    }
}

void fps_metal_profile_mark(unsigned stage,bool end) {
    if(!metal_state.stage_samples||!metal_state.stage_resolved||!metal_state.encoder||stage>=FPS_METAL_STAGES)return;
    if(!end){
        if(metal_state.stage_pairs>=1024)return;
        unsigned pair=metal_state.stage_pairs++;metal_state.stage_open[stage]=pair;metal_state.stage_kind[pair]=stage;
        if(metal_state.stage_encoder_mode) {
            [metal_state.encoder endEncoding];[metal_state.encoder release];
            MTLComputePassDescriptor *descriptor=[MTLComputePassDescriptor computePassDescriptor];
            descriptor.sampleBufferAttachments[0].sampleBuffer=metal_state.stage_samples;
            descriptor.sampleBufferAttachments[0].startOfEncoderSampleIndex=2*pair;
            descriptor.sampleBufferAttachments[0].endOfEncoderSampleIndex=2*pair+1;
            metal_state.encoder=[[metal_state.command_buffer computeCommandEncoderWithDescriptor:descriptor] retain];
            fps_metal_restore_bindings();
        } else [metal_state.encoder sampleCountersInBuffer:metal_state.stage_samples atSampleIndex:2*pair withBarrier:YES];
    } else {
        unsigned pair=metal_state.stage_open[stage];if(pair==UINT32_MAX)return;
        if(metal_state.stage_encoder_mode) {
            [metal_state.encoder endEncoding];[metal_state.encoder release];
            metal_state.encoder=[[metal_state.command_buffer computeCommandEncoder] retain];fps_metal_restore_bindings();
        } else [metal_state.encoder sampleCountersInBuffer:metal_state.stage_samples atSampleIndex:2*pair+1 withBarrier:YES];
        metal_state.stage_complete[pair]=true;metal_state.stage_open[stage]=UINT32_MAX;
    }
}

void *fps_metal_device(void) { return (void *)metal_state.device; }
void *fps_metal_encoder(void) { return (void *)metal_state.encoder; }
void fps_metal_restore_bindings(void) {
    if (!metal_state.encoder) return;
    [metal_state.encoder memoryBarrierWithScope:MTLBarrierScopeBuffers];
    [metal_state.encoder setComputePipelineState:metal_state.pipeline];
    for (int i=0;i<FPS_GPU_BUFFER_COUNT;i++)
        [metal_state.encoder setBuffer:metal_state.bound[i] offset:0 atIndex:(NSUInteger)i];
}

static bool fps_metal_encode(int mode, int count, id<MTLBuffer> indirect, size_t offset) {
    @autoreleasepool {
        if (!metal_state.command_buffer || !metal_state.encoder || !metal_state.pipeline || count < 0) return false;
        double t0 = metal_time_now_ms();
        int profile_stage=mode==GPU_MODE_VGS?FPS_METAL_SHAPE:
            (mode==GPU_MODE_HASH_CLEAR||mode==GPU_MODE_HASH_BUILD||mode==GPU_MODE_PAIR_COLLISIONS||mode==GPU_MODE_STATIC_COLLISIONS)?FPS_METAL_CONTACT:-1;
        if(profile_stage>=0)fps_metal_profile_mark((unsigned)profile_stage,false);
        if (metal_profile.last_dispatch_count > 0) {
            [metal_state.encoder memoryBarrierWithScope:MTLBarrierScopeBuffers];
        }
        metal_state.uniforms.mode = mode;
        [metal_state.encoder setBytes:&metal_state.uniforms length:sizeof(metal_state.uniforms) atIndex:FPS_GPU_BUFFER_COUNT];
        MTLSize threads = MTLSizeMake(128, 1, 1);
        if (indirect) {
            [metal_state.encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                                   indirectBufferOffset:offset
                                                  threadsPerThreadgroup:threads];
        } else if (count > 0) {
            MTLSize groups = MTLSizeMake(((NSUInteger)count + 127u) / 128u, 1, 1);
            [metal_state.encoder dispatchThreadgroups:groups threadsPerThreadgroup:threads];
        }
        if(profile_stage>=0)fps_metal_profile_mark((unsigned)profile_stage,true);
        metal_profile.last_encode_ms += (metal_time_now_ms() - t0);
        metal_profile.last_dispatch_count++;
        return true;
    }
}

bool fps_metal_dispatch(int mode, int count) {
    return count <= 0 || fps_metal_encode(mode, count, nil, 0);
}

bool fps_metal_dispatch_indirect(int mode, void *buffer, size_t offset) {
    return fps_metal_encode(mode, 1, (id<MTLBuffer>)buffer, offset);
}

bool fps_metal_commit_batch(bool wait) {
    @autoreleasepool {
        if (!metal_state.command_buffer) return false;
        if (metal_state.encoder) {
            [metal_state.encoder endEncoding];
            [metal_state.encoder release];
            metal_state.encoder = nil;
        }
        if(metal_state.stage_samples&&metal_state.stage_resolved&&metal_state.stage_pairs) {
            id<MTLBlitCommandEncoder> counter_blit=[metal_state.command_buffer blitCommandEncoder];
            [counter_blit resolveCounters:metal_state.stage_samples inRange:NSMakeRange(0,2*metal_state.stage_pairs) destinationBuffer:metal_state.stage_resolved destinationOffset:0];
            [counter_blit endEncoding];
        }
        if (metal_state.managed) {
            static const int readback_slots[] = {
                FPS_GPU_BUFFER_PARTICLE, FPS_GPU_BUFFER_CELL, FPS_GPU_BUFFER_SIM_ID, FPS_GPU_BUFFER_VOXEL,
                FPS_GPU_BUFFER_CONTROL, FPS_GPU_BUFFER_CLONE_PARENT
            };
            id<MTLBlitCommandEncoder> blit = [metal_state.command_buffer blitCommandEncoder];
            if (!blit) {
                [metal_state.command_buffer release];
                metal_state.command_buffer = nil;
                return false;
            }
            for (size_t i = 0; i < sizeof(readback_slots)/sizeof(readback_slots[0]); ++i) {
                id<MTLBuffer> buffer = metal_state.bound[readback_slots[i]];
                if (buffer) [blit synchronizeResource:buffer];
            }
            [blit endEncoding];
        }
        metal_state.wait_start = metal_time_now_ms();
        [metal_state.command_buffer commit];
        if (wait) {
            return fps_metal_wait_batch();
        }
        return true;
    }
}

bool fps_metal_wait_batch(void) {
    @autoreleasepool {
        if (!metal_state.command_buffer) return true;
        [metal_state.command_buffer waitUntilCompleted];
        double t_wait_end = metal_time_now_ms();
        metal_profile.last_wait_ms = t_wait_end - metal_state.wait_start;
        CFTimeInterval gpu_start = metal_state.command_buffer.GPUStartTime;
        CFTimeInterval gpu_end = metal_state.command_buffer.GPUEndTime;
        metal_profile.last_gpu_exec_ms = (gpu_end > gpu_start) ? (gpu_end - gpu_start) * 1000.0 : 0.0;
        bool ok = metal_state.command_buffer.status == MTLCommandBufferStatusCompleted;
        if(ok&&metal_state.stage_resolved) {
            const MTLCounterResultTimestamp *samples=metal_state.stage_resolved.contents;
            for(unsigned i=0;i<metal_state.stage_pairs;i++)if(metal_state.stage_complete[i]) {
                uint64_t begin=samples[2*i].timestamp,end=samples[2*i+1].timestamp;
                if(begin!=MTLCounterErrorValue&&end!=MTLCounterErrorValue&&end>=begin)metal_profile.stage_ticks[metal_state.stage_kind[i]]+=end-begin;
            }
        }
        [metal_state.command_buffer release];
        metal_state.command_buffer = nil;
        return ok;
    }
}

bool fps_metal_end_batch(void) {
    return fps_metal_commit_batch(true);
}

#else
typedef int fps_metal_translation_unit_is_empty;
#endif
