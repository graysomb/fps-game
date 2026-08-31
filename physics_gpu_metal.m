#if defined(__APPLE__) && defined(FPS_GPU_METAL)

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "physics_gpu_metal.h"
#include <stdio.h>
#include <string.h>

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
        return metal_state.command_buffer != nil;
    }
}

static bool fps_metal_encode(int mode, int count, id<MTLBuffer> indirect, size_t offset) {
    @autoreleasepool {
        if (!metal_state.command_buffer || !metal_state.pipeline || count < 0) return false;
        id<MTLComputeCommandEncoder> encoder = [metal_state.command_buffer computeCommandEncoder];
        if (!encoder) return false;
        [encoder setComputePipelineState:metal_state.pipeline];
        for (int i = 0; i < FPS_GPU_BUFFER_COUNT; ++i) {
            if (!metal_state.bound[i]) { [encoder endEncoding]; return false; }
            [encoder setBuffer:metal_state.bound[i] offset:0 atIndex:(NSUInteger)i];
        }
        metal_state.uniforms.mode = mode;
        [encoder setBytes:&metal_state.uniforms length:sizeof(metal_state.uniforms) atIndex:FPS_GPU_BUFFER_COUNT];
        MTLSize threads = MTLSizeMake(128, 1, 1);
        if (indirect) {
            [encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                       indirectBufferOffset:offset
                                      threadsPerThreadgroup:threads];
        } else if (count > 0) {
            MTLSize groups = MTLSizeMake(((NSUInteger)count + 127u) / 128u, 1, 1);
            [encoder dispatchThreadgroups:groups threadsPerThreadgroup:threads];
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
        [metal_state.command_buffer commit];
        [metal_state.command_buffer waitUntilCompleted];
        bool ok = metal_state.command_buffer.status == MTLCommandBufferStatusCompleted;
        [metal_state.command_buffer release];
        metal_state.command_buffer = nil;
        return ok;
    }
}

#else
typedef int fps_metal_translation_unit_is_empty;
#endif
