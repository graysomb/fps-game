#ifndef FPS_ADAPTIVE_METAL_H
#define FPS_ADAPTIVE_METAL_H
#include "adaptive_mesh.h"
#include <stdbool.h>
typedef struct AdaptiveMetal AdaptiveMetal;
typedef struct {
    uint32_t leaves, requests, error, generation, new_count, changed, active, input_count;
    uint32_t reserved[8];
    uint32_t indirect[4];
} AdaptiveGpuHeader;
/* Borrow the existing device/encoder; none of these methods submits or waits. */
AdaptiveMetal *adaptive_metal_create(void *device,const char *source_path,size_t capacity,
                                    char *error,size_t error_capacity);
void adaptive_metal_destroy(AdaptiveMetal *);
bool adaptive_metal_upload(AdaptiveMetal *,const AdaptiveInput *,size_t count);
bool adaptive_metal_encode(AdaptiveMetal *,void *encoder,unsigned max_level,bool surface,bool balance);
bool adaptive_metal_encode_refine(AdaptiveMetal *,void *encoder,void *voxel_buffer,unsigned max_level);
void *adaptive_metal_fine_buffer(AdaptiveMetal *);
void *adaptive_metal_leaf_buffer(AdaptiveMetal *);
void *adaptive_metal_header_buffer(AdaptiveMetal *);
uint64_t adaptive_metal_encoded_dispatches(const AdaptiveMetal *);
size_t adaptive_metal_scratch_bytes(const AdaptiveMetal *);
/* For tests/diagnostics only, after caller's existing completion boundary. */
bool adaptive_metal_read(const AdaptiveMetal *,AdaptiveMesh *,AdaptiveGpuHeader *);
#endif
