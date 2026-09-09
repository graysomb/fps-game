#if defined(__APPLE__) && defined(FPS_GPU_METAL)

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "physics_gpu_metal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { FPS_METAL_PROFILE_STAGE_COUNT = 25, FPS_METAL_PROFILE_MAX_DISPATCHES = 2048 };

typedef struct FpsMetalState {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLLibrary> library;
    id<MTLComputePipelineState> pipeline;
    id<MTLCommandBuffer> command_buffer;
    id<MTLBuffer> bound[FPS_GPU_BUFFER_COUNT];
    FpsGpuUniforms uniforms;
    MTLResourceOptions resource_options;
    bool managed;
    bool timing_enabled;
    bool profile_enabled;
    id<MTLCounterSampleBuffer> profile_samples;
    NSUInteger profile_dispatch_count;
    int profile_modes[FPS_METAL_PROFILE_MAX_DISPATCHES];
    uint64_t profile_ticks[FPS_METAL_PROFILE_STAGE_COUNT];
    uint64_t profile_dispatches[FPS_METAL_PROFILE_STAGE_COUNT];
    double profile_gpu_seconds;
} FpsMetalState;

static FpsMetalState metal_state;

static const char *fps_metal_mode_name(int mode) {
    switch (mode) {
        case GPU_MODE_RESET: return "reset";
        case GPU_MODE_INTEGRATE: return "integrate";
        case GPU_MODE_HASH_CLEAR: return "hash-clear";
        case GPU_MODE_HASH_BUILD: return "hash-build";
        case GPU_MODE_SCENE_COLLISIONS: return "scene-contacts";
        case GPU_MODE_PAIR_COLLISIONS: return "pair-contacts";
        case GPU_MODE_APPLY: return "apply";
        case GPU_MODE_VGS: return "vgs";
        case GPU_MODE_STATIC_COLLISIONS: return "static-contacts";
        case GPU_MODE_BREAK_MASK: return "break-mask";
        case GPU_MODE_FINALIZE_PARTICLES: return "finalize-particles";
        case GPU_MODE_FINALIZE_VOXELS: return "finalize-voxels";
        case GPU_MODE_SPLIT_BREAKS: return "split-breaks";
        case GPU_MODE_PREPARE_INDIRECT: return "prepare-indirect";
        case GPU_MODE_WAKE_GATHER: return "wake-gather";
        case GPU_MODE_WAKE_APPLY: return "wake-apply";
        case GPU_MODE_TOPOLOGY_REBUILD_SERIAL: return "topology-rebuild";
        case GPU_MODE_REFRESH_DEPENDENCIES: return "refresh-interfaces";
        case GPU_MODE_FINALIZE_DEPENDENCIES: return "finalize-dependencies";
        case GPU_MODE_LIFT_GREEDY_FLOOR: return "floor-lift";
        case GPU_MODE_RESET_STATIC_MAPPED: return "reset-static-mapped";
        case GPU_MODE_APPLY_STATIC_MAPPED: return "apply-static-mapped";
        case GPU_MODE_BUILD_RENDER_MATRICES: return "render-matrices";
        case GPU_MODE_GREEDY_FLOOR_ISLANDS: return "floor-islands";
        case GPU_MODE_ATTACHMENTS: return "attachments";
        default:
            if (mode >= GPU_MODE_GREEDY_FLOOR_BATCH_BASE) return "floor-batch";
            if (mode >= GPU_MODE_STATIC_GREEDY_BATCH_BASE) return "static-batch";
            return "unknown";
    }
}

static int fps_metal_profile_stage(int mode) {
    if (mode >= 0 && mode < FPS_METAL_PROFILE_STAGE_COUNT) return mode;
    if (mode >= GPU_MODE_GREEDY_FLOOR_BATCH_BASE) return GPU_MODE_GREEDY_FLOOR_ISLANDS;
    if (mode >= GPU_MODE_STATIC_GREEDY_BATCH_BASE) return GPU_MODE_STATIC_COLLISIONS;
    return -1;
}

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
        metal_state.library = [metal_state.device newLibraryWithFile:path error:&library_error];
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
        metal_state.timing_enabled = getenv("FPS_METAL_GPU_TIME") != NULL ||
                                     getenv("FPS_METAL_STAGE_PROFILE") != NULL;
        if (getenv("FPS_METAL_STAGE_PROFILE") != NULL &&
            [metal_state.device supportsCounterSampling:MTLCounterSamplingPointAtStageBoundary]) {
            id<MTLCounterSet> timestamp_set = nil;
            for (id<MTLCounterSet> set in metal_state.device.counterSets) {
                if ([set.name isEqualToString:MTLCommonCounterSetTimestamp]) {
                    timestamp_set = set;
                    break;
                }
            }
            if (timestamp_set) {
                MTLCounterSampleBufferDescriptor *descriptor = [[MTLCounterSampleBufferDescriptor alloc] init];
                descriptor.counterSet = timestamp_set;
                descriptor.storageMode = MTLStorageModeShared;
                descriptor.sampleCount = FPS_METAL_PROFILE_MAX_DISPATCHES * 2u;
                descriptor.label = @"FPS PBD stage timestamps";
                NSError *sample_error = nil;
                metal_state.profile_samples = [metal_state.device
                    newCounterSampleBufferWithDescriptor:descriptor error:&sample_error];
                [descriptor release];
                metal_state.profile_enabled = metal_state.profile_samples != nil;
                if (!metal_state.profile_enabled)
                    fprintf(stderr, "metal-stage-profile unavailable: %s\n",
                            sample_error ? sample_error.localizedDescription.UTF8String : "counter buffer creation failed");
            }
        }
        return true;
    }
}

void fps_metal_shutdown(void) {
    @autoreleasepool {
        if (metal_state.command_buffer) {
            [metal_state.command_buffer waitUntilCompleted];
            [metal_state.command_buffer release];
        }
        [metal_state.pipeline release];
        [metal_state.library release];
        if (metal_state.timing_enabled && !metal_state.profile_enabled)
            fprintf(stderr, "metal-gpu-time gpuMs=%.3f\n", metal_state.profile_gpu_seconds * 1000.0);
        if (metal_state.profile_enabled) {
            uint64_t total_ticks = 0;
            for (int i = 0; i < FPS_METAL_PROFILE_STAGE_COUNT; ++i) total_ticks += metal_state.profile_ticks[i];
            fprintf(stderr, "metal-stage-profile gpuMs=%.3f sampledTicks=%llu\n",
                    metal_state.profile_gpu_seconds * 1000.0, (unsigned long long)total_ticks);
            for (int i = 0; i < FPS_METAL_PROFILE_STAGE_COUNT; ++i) {
                if (metal_state.profile_dispatches[i] == 0) continue;
                double fraction = total_ticks ? (double)metal_state.profile_ticks[i] / (double)total_ticks : 0.0;
                fprintf(stderr, "metal-stage mode=%s dispatches=%llu gpuMs=%.3f share=%.2f%%\n",
                        fps_metal_mode_name(i),
                        (unsigned long long)metal_state.profile_dispatches[i],
                        metal_state.profile_gpu_seconds * 1000.0 * fraction, fraction * 100.0);
            }
        }
        [metal_state.profile_samples release];
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
    if (!handle) return;
    [(id<MTLBuffer>)handle release];
}

bool fps_metal_buffer_update(void *handle, const void *data, size_t size, size_t offset) {
    id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
    if (!buffer || !data || offset > buffer.length || size > buffer.length - offset) return false;
    memcpy((uint8_t *)buffer.contents + offset, data, size);
    if (metal_state.managed) [buffer didModifyRange:NSMakeRange(offset, size)];
    return true;
}

bool fps_metal_buffer_read(void *handle, void *data, size_t size, size_t offset) {
    id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
    if (!buffer || !data || offset > buffer.length || size > buffer.length - offset) return false;
    memcpy(data, (const uint8_t *)buffer.contents + offset, size);
    return true;
}

const void *fps_metal_buffer_contents(void *handle) {
    id<MTLBuffer> buffer = (id<MTLBuffer>)handle;
    return buffer ? buffer.contents : NULL;
}

void fps_metal_bind_buffer(int slot, void *handle) {
    if (slot < 0 || slot >= FPS_GPU_BUFFER_COUNT) return;
    metal_state.bound[slot] = (id<MTLBuffer>)handle;
}

void fps_metal_set_uniforms(const FpsGpuUniforms *uniforms) {
    if (uniforms) metal_state.uniforms = *uniforms;
}

bool fps_metal_begin_batch(void) {
    @autoreleasepool {
        if (!metal_state.queue || metal_state.command_buffer) return false;
        metal_state.command_buffer = [[metal_state.queue commandBuffer] retain];
        metal_state.profile_dispatch_count = 0;
        return metal_state.command_buffer != nil;
    }
}

static bool fps_metal_encode(int mode, int count, id<MTLBuffer> indirect, size_t offset) {
    @autoreleasepool {
        if (!metal_state.command_buffer || !metal_state.pipeline || count < 0) return false;
        NSUInteger profile_slot = metal_state.profile_dispatch_count;
        bool profile_this = metal_state.profile_enabled &&
                            profile_slot < FPS_METAL_PROFILE_MAX_DISPATCHES;
        id<MTLComputeCommandEncoder> encoder = nil;
        if (profile_this) {
            MTLComputePassDescriptor *descriptor = [MTLComputePassDescriptor computePassDescriptor];
            MTLComputePassSampleBufferAttachmentDescriptor *attachment =
                descriptor.sampleBufferAttachments[0];
            attachment.sampleBuffer = metal_state.profile_samples;
            attachment.startOfEncoderSampleIndex = profile_slot * 2u;
            attachment.endOfEncoderSampleIndex = profile_slot * 2u + 1u;
            encoder = [metal_state.command_buffer computeCommandEncoderWithDescriptor:descriptor];
        } else {
            encoder = [metal_state.command_buffer computeCommandEncoder];
        }
        if (!encoder) return false;
        [encoder setComputePipelineState:metal_state.pipeline];
        for (int i = 0; i < FPS_GPU_BUFFER_COUNT; ++i) {
            if (!metal_state.bound[i]) { [encoder endEncoding]; return false; }
            [encoder setBuffer:metal_state.bound[i] offset:0 atIndex:(NSUInteger)i];
        }
        metal_state.uniforms.mode = mode;
        [encoder setBytes:&metal_state.uniforms length:sizeof(metal_state.uniforms) atIndex:FPS_GPU_BUFFER_COUNT];
        if (mode == GPU_MODE_GREEDY_FLOOR_ISLANDS) {
            const NSUInteger floor_bytes =
                (NSUInteger)FPS_GPU_MAX_FLOOR_ISLAND_CONTROLS * 2u * sizeof(float) * 4u;
            [encoder setThreadgroupMemoryLength:floor_bytes atIndex:0];
        }
        MTLSize threads = MTLSizeMake(128, 1, 1);
        if (indirect) {
            [encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                       indirectBufferOffset:offset
                                      threadsPerThreadgroup:threads];
        } else if (count > 0) {
            MTLSize groups = MTLSizeMake(((NSUInteger)count + 127u) / 128u, 1, 1);
            [encoder dispatchThreadgroups:groups threadsPerThreadgroup:threads];
        }
        if (profile_this) {
            metal_state.profile_modes[profile_slot] = mode;
            metal_state.profile_dispatch_count++;
        }
        [encoder endEncoding];
        return true;
    }
}

bool fps_metal_dispatch(int mode, int count) {
    return count <= 0 || fps_metal_encode(mode, count, nil, 0);
}

bool fps_metal_dispatch_indirect(int mode, void *buffer, size_t offset) {
    return fps_metal_encode(mode, 1, (id<MTLBuffer>)buffer, offset);
}

bool fps_metal_end_batch(void) {
    @autoreleasepool {
        if (!metal_state.command_buffer) return false;
        if (metal_state.managed) {
            static const int readback_slots[] = {
                FPS_GPU_BUFFER_PARTICLE, FPS_GPU_BUFFER_CELL, FPS_GPU_BUFFER_SIM_ID, FPS_GPU_BUFFER_VOXEL,
                FPS_GPU_BUFFER_CONTROL, FPS_GPU_BUFFER_CLONE_PARENT, FPS_GPU_BUFFER_RENDER_MATRIX
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
        [metal_state.command_buffer commit];
        [metal_state.command_buffer waitUntilCompleted];
        bool ok = metal_state.command_buffer.status == MTLCommandBufferStatusCompleted;
        if (ok && metal_state.timing_enabled && !metal_state.profile_enabled) {
            CFTimeInterval gpu_start = metal_state.command_buffer.GPUStartTime;
            CFTimeInterval gpu_end = metal_state.command_buffer.GPUEndTime;
            if (gpu_end >= gpu_start) metal_state.profile_gpu_seconds += gpu_end - gpu_start;
        }
        if (ok && metal_state.profile_enabled && metal_state.profile_dispatch_count > 0) {
            CFTimeInterval gpu_start = metal_state.command_buffer.GPUStartTime;
            CFTimeInterval gpu_end = metal_state.command_buffer.GPUEndTime;
            if (gpu_end >= gpu_start) metal_state.profile_gpu_seconds += gpu_end - gpu_start;
            NSUInteger count = metal_state.profile_dispatch_count * 2u;
            NSData *resolved = [metal_state.profile_samples resolveCounterRange:NSMakeRange(0, count)];
            const MTLCounterResultTimestamp *samples = resolved.bytes;
            if (samples && resolved.length >= count * sizeof(*samples)) {
                for (NSUInteger i = 0; i < metal_state.profile_dispatch_count; ++i) {
                    uint64_t begin = samples[i * 2u].timestamp;
                    uint64_t end = samples[i * 2u + 1u].timestamp;
                    int stage = fps_metal_profile_stage(metal_state.profile_modes[i]);
                    if (stage < 0 || begin == MTLCounterErrorValue || end == MTLCounterErrorValue || end < begin) continue;
                    metal_state.profile_ticks[stage] += end - begin;
                    metal_state.profile_dispatches[stage]++;
                }
            }
        }
        [metal_state.command_buffer release];
        metal_state.command_buffer = nil;
        return ok;
    }
}

#else
typedef int fps_metal_translation_unit_is_empty;
#endif
