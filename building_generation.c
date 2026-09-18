#include "building_generation.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t mix32(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15;
    x *= 0x846ca68bu; x ^= x >> 16; return x;
}

uint32_t building_seed_stream(uint32_t seed, BuildingStyle style,
                              uint32_t family, uint32_t purpose) {
    uint32_t value = seed ^ (0x9e3779b9u * ((uint32_t)style + 1u));
    value ^= 0x85ebca6bu * (family + 1u);
    value ^= 0xc2b2ae35u * (purpose + 1u);
    value = mix32(value);
    return value ? value : 0x6d2b79f5u;
}

uint32_t building_rng_next(uint32_t *state) {
    uint32_t x = state && *state ? *state : 0x6d2b79f5u;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    if (state) *state = x; return x;
}

int building_rng_range(uint32_t *state, int lo, int hi) {
    if (hi <= lo) return lo;
    return lo + (int)(building_rng_next(state) % (uint32_t)(hi - lo + 1));
}

static int voxel_compare(const void *ap, const void *bp) {
    const BuildingVoxelSpec *a = ap, *b = bp;
    if (a->y != b->y) return a->y < b->y ? -1 : 1;
    if (a->z != b->z) return a->z < b->z ? -1 : 1;
    if (a->x != b->x) return a->x < b->x ? -1 : 1;
    return 0;
}

bool building_blueprint_init(BuildingBlueprint *b, size_t capacity,
                             BuildingStyle style, uint32_t seed, uint32_t family,
                             uint32_t structural_group) {
    if (!b || !capacity) return false;
    memset(b, 0, sizeof(*b));
    b->voxels = malloc(capacity * sizeof(*b->voxels));
    b->slot_capacity = 1; while (b->slot_capacity < capacity * 2) b->slot_capacity <<= 1;
    b->slots = malloc(b->slot_capacity * sizeof(*b->slots));
    if (!b->voxels || !b->slots) { free(b->voxels); free(b->slots); memset(b, 0, sizeof(*b)); return false; }
    for (size_t i = 0; i < b->slot_capacity; ++i) b->slots[i] = -1;
    b->capacity = capacity; b->min_y = INT_MAX; b->style = (uint32_t)style;
    b->seed = seed; b->family = family; b->structural_group = structural_group;
    return true;
}

void building_blueprint_release(BuildingBlueprint *b) {
    if (!b) return; free(b->voxels); free(b->slots); memset(b, 0, sizeof(*b));
}

static uint32_t coord_hash(int x, int y, int z) {
    return mix32((uint32_t)x * 0x8da6b343u ^ (uint32_t)y * 0xd8163841u ^
                 (uint32_t)z * 0xcb1ab31fu);
}

static ptrdiff_t find_index(const BuildingBlueprint *b, int x, int y, int z) {
    if (!b || !b->slots || !b->slot_capacity) return -1;
    size_t mask = b->slot_capacity - 1, slot = coord_hash(x, y, z) & mask;
    for (size_t probe = 0; probe < b->slot_capacity; ++probe, slot = (slot + 1) & mask) {
        int32_t index = b->slots[slot];
        if (index == -1) return -1;
        if (index >= 0) { const BuildingVoxelSpec *v = &b->voxels[index];
            if (v->x == x && v->y == y && v->z == z) return index; }
    }
    return -1;
}

static bool index_voxel(BuildingBlueprint *b, size_t index) {
    BuildingVoxelSpec *v = &b->voxels[index];
    size_t mask = b->slot_capacity - 1, slot = coord_hash(v->x, v->y, v->z) & mask, tomb = SIZE_MAX;
    for (size_t probe = 0; probe < b->slot_capacity; ++probe, slot = (slot + 1) & mask) {
        if (b->slots[slot] == -2 && tomb == SIZE_MAX) tomb = slot;
        if (b->slots[slot] == -1) { b->slots[tomb == SIZE_MAX ? slot : tomb] = (int32_t)index; return true; }
    }
    return false;
}

bool building_blueprint_put(BuildingBlueprint *b, BuildingVoxelSpec v,
                            BuildingValidationReport *report) {
    if (!b || !b->voxels) return false;
    ptrdiff_t existing = find_index(b, v.x, v.y, v.z);
    if (existing >= 0) {
        BuildingVoxelSpec *old = &b->voxels[existing];
        if (report) report->duplicate_writes++;
        if (old->material != v.material && report) report->conflict_writes++;
        /* Structural information is monotonic across overlapping primitives. */
        old->rgba = v.rgba; old->flags |= v.flags;
        if (v.role < old->role) old->role = v.role;
        if (old->source_node < 0) old->source_node = v.source_node;
        return true;
    }
    if (b->count == b->capacity) return false;
    b->voxels[b->count] = v;
    if (!index_voxel(b, b->count)) return false;
    b->count++;
    if (v.y < b->min_y) b->min_y = v.y;
    return true;
}

bool building_blueprint_remove(BuildingBlueprint *b, int x, int y, int z) {
    ptrdiff_t index = find_index(b, x, y, z);
    if (index < 0) return true;
    size_t mask = b->slot_capacity - 1, slot = coord_hash(x, y, z) & mask;
    while (b->slots[slot] != index) slot = (slot + 1) & mask;
    b->slots[slot] = -2;
    size_t last = --b->count;
    if ((size_t)index != last) {
        BuildingVoxelSpec moved = b->voxels[last];
        ptrdiff_t moved_index = find_index(b, moved.x, moved.y, moved.z);
        if (moved_index != (ptrdiff_t)last) return false;
        size_t moved_slot = coord_hash(moved.x, moved.y, moved.z) & mask;
        while (b->slots[moved_slot] != (int32_t)last) moved_slot = (moved_slot + 1) & mask;
        b->voxels[index] = moved; b->slots[moved_slot] = (int32_t)index;
    }
    return true;
}

const BuildingVoxelSpec *building_blueprint_find(const BuildingBlueprint *b,
                                                  int x, int y, int z) {
    ptrdiff_t index = find_index(b, x, y, z);
    return index < 0 ? NULL : &b->voxels[index];
}

static bool is_anchor(const BuildingVoxelSpec *v, int min_y) {
    return (v->flags & BUILDING_VOXEL_ANCHOR) || v->y == min_y;
}

static bool mark_grounded(const BuildingBlueprint *b, unsigned char *grounded) {
    size_t *queue = malloc((b->count ? b->count : 1) * sizeof(*queue));
    if (!queue) return false;
    size_t head = 0, tail = 0;
    for (size_t i = 0; i < b->count; ++i) if (is_anchor(&b->voxels[i], b->min_y)) {
        grounded[i] = 1; queue[tail++] = i;
    }
    static const int delta[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    while (head < tail) {
        BuildingVoxelSpec v = b->voxels[queue[head++]];
        for (int d = 0; d < 6; ++d) {
            ptrdiff_t j = find_index(b, v.x + delta[d][0], v.y + delta[d][1], v.z + delta[d][2]);
            if (j >= 0 && !grounded[j]) { grounded[j] = 1; queue[tail++] = (size_t)j; }
        }
    }
    free(queue); return true;
}

static bool add_pier(BuildingBlueprint *b, const BuildingVoxelSpec *from,
                     BuildingValidationReport *report) {
    for (int y = from->y - 1; y >= b->min_y; --y) {
        if (building_blueprint_find(b, from->x, y, from->z)) break;
        BuildingVoxelSpec pier = *from; pier.y = y; pier.role = BUILDING_ROLE_BUTTRESS;
        pier.flags = BUILDING_VOXEL_DAMAGEABLE | BUILDING_VOXEL_REPAIR;
        if (y == b->min_y) pier.flags |= BUILDING_VOXEL_ANCHOR;
        if (!building_blueprint_put(b, pier, report)) return false;
        if (report) report->piers_added++;
    }
    return true;
}

bool building_blueprint_validate_and_repair(BuildingBlueprint *b,
                                             BuildingValidationReport *report) {
    if (!b || !b->count) return false;
    BuildingValidationReport local = {0}; if (!report) report = &local;
    report->valid = false; report->reason[0] = 0;
    /* The lowest physical layer is the visible foundation and remains fixed. */
    for (size_t i = 0; i < b->count; ++i) if (b->voxels[i].y == b->min_y) {
        b->voxels[i].role = BUILDING_ROLE_FOUNDATION;
        b->voxels[i].flags |= BUILDING_VOXEL_ANCHOR;
    }
    for (unsigned pass = 0; pass < 8; ++pass) {
        size_t count = b->count;
        unsigned char *grounded = calloc(count, 1);
        if (!grounded || !mark_grounded(b, grounded)) { free(grounded); snprintf(report->reason, sizeof(report->reason), "allocation failure"); return false; }
        size_t first = count;
        for (size_t i = 0; i < count; ++i) if (!grounded[i] && !(b->voxels[i].flags & BUILDING_VOXEL_SEMANTIC)) { first = i; break; }
        free(grounded);
        if (first == count) {
            qsort(b->voxels, b->count, sizeof(*b->voxels), voxel_compare);
            for (size_t i = 0; i < b->slot_capacity; ++i) b->slots[i] = -1;
            for (size_t i = 0; i < b->count; ++i) if (!index_voxel(b, i)) {
                snprintf(report->reason, sizeof(report->reason), "failed to rebuild occupancy index"); return false;
            }
            report->valid = true; report->repair_passes = pass; return true;
        }
        report->disconnected_components++;
        if (!add_pier(b, &b->voxels[first], report)) { snprintf(report->reason, sizeof(report->reason), "blueprint capacity exhausted during repair"); return false; }
    }
    snprintf(report->reason, sizeof(report->reason), "structure remained disconnected after eight repairs");
    return false;
}

uint64_t building_blueprint_fingerprint(const BuildingBlueprint *b) {
    if (!b) return 0;
    uint64_t h = UINT64_C(1469598103934665603);
    h = (h ^ b->style) * UINT64_C(1099511628211);
    h = (h ^ b->family) * UINT64_C(1099511628211);
    for (size_t i = 0; i < b->count; ++i) {
        const BuildingVoxelSpec *v = &b->voxels[i];
        uint64_t q = (uint32_t)(v->x / 4) ^ ((uint64_t)(uint32_t)(v->y / 4) << 21) ^
                     ((uint64_t)(uint32_t)(v->z / 4) << 42) ^ ((uint64_t)v->role << 57);
        h = (h ^ q) * UINT64_C(1099511628211);
    }
    return h;
}
