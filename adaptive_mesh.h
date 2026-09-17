#ifndef FPS_ADAPTIVE_MESH_H
#define FPS_ADAPTIVE_MESH_H
#include <stdint.h>
#include <stddef.h>
/* All coordinates are integer fine-grid coordinates, never current positions.
 * A leaf's key remains in fine Morton units at every level. */
typedef struct {
    uint64_t key;
    uint32_t arena, material, level, flags, source, reserved;
} AdaptiveLeaf;
typedef struct { int32_t xyz[3]; uint32_t arena, material, flags; } AdaptiveInput;
typedef struct { AdaptiveLeaf *leaves; size_t count; } AdaptiveMesh;
enum { ADAPTIVE_PROTECTED = 1, ADAPTIVE_COORD_LIMIT = 1 << 20 };
enum { ADAPTIVE_SUCCESS, ADAPTIVE_BAD_INPUT, ADAPTIVE_ALLOCATION_FAILED };
/* Transactional: out must be zero-initialized or owned by this API. */
int adaptive_mesh_reference(const AdaptiveInput *, size_t, unsigned max_level,
                            int protect_surface, int balance, AdaptiveMesh *out);
void adaptive_mesh_release(AdaptiveMesh *);
uint64_t adaptive_morton(const int32_t xyz[3]);
void adaptive_origin(uint64_t key, int32_t xyz[3]);
#endif
