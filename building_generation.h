#ifndef BUILDING_GENERATION_H
#define BUILDING_GENERATION_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
    BUILDING_STYLE_GREEK,
    BUILDING_STYLE_MEGALITH,
    BUILDING_STYLE_HYPERBOREAN,
    BUILDING_STYLE_FORERUNNER,
    BUILDING_STYLE_CITADEL,
    BUILDING_STYLE_COUNT
} BuildingStyle;

typedef enum {
    BUILDING_ROLE_FOUNDATION,
    BUILDING_ROLE_COLUMN,
    BUILDING_ROLE_WALL,
    BUILDING_ROLE_BEAM,
    BUILDING_ROLE_SLAB,
    BUILDING_ROLE_BUTTRESS,
    BUILDING_ROLE_DECORATION,
    BUILDING_ROLE_VOID,
    BUILDING_ROLE_SEMANTIC_ANCHOR
} BuildingStructuralRole;

enum {
    BUILDING_VOXEL_ANCHOR = 1u << 0,
    BUILDING_VOXEL_SEMANTIC = 1u << 1,
    BUILDING_VOXEL_DAMAGEABLE = 1u << 2,
    BUILDING_VOXEL_REPAIR = 1u << 3
};

typedef struct {
    int32_t x, y, z;
    uint32_t rgba;
    uint16_t material;
    uint16_t role;
    int32_t source_node;
    uint32_t flags;
} BuildingVoxelSpec;

typedef struct {
    BuildingVoxelSpec *voxels;
    size_t count, capacity;
    int32_t *slots;
    size_t slot_capacity;
    int32_t min_y;
    uint32_t seed;
    uint32_t style, family, structural_group;
} BuildingBlueprint;

typedef struct {
    bool valid;
    uint32_t repair_passes;
    uint32_t piers_added;
    uint32_t bearing_repairs;
    uint32_t unsupported_members;
    uint32_t disconnected_components;
    uint32_t duplicate_writes;
    uint32_t conflict_writes;
    char reason[160];
} BuildingValidationReport;

/* Stable integer streams. Purpose constants isolate layout, structure,
 * decoration, and gameplay so one subsystem cannot perturb another. */
uint32_t building_seed_stream(uint32_t seed, BuildingStyle style,
                              uint32_t family, uint32_t purpose);
uint32_t building_rng_next(uint32_t *state);
int building_rng_range(uint32_t *state, int lo, int hi);

bool building_blueprint_init(BuildingBlueprint *blueprint, size_t capacity,
                             BuildingStyle style, uint32_t seed, uint32_t family,
                             uint32_t structural_group);
void building_blueprint_release(BuildingBlueprint *blueprint);
bool building_blueprint_put(BuildingBlueprint *blueprint, BuildingVoxelSpec voxel,
                            BuildingValidationReport *report);
bool building_blueprint_remove(BuildingBlueprint *blueprint, int x, int y, int z);
const BuildingVoxelSpec *building_blueprint_find(const BuildingBlueprint *blueprint,
                                                  int x, int y, int z);
bool building_blueprint_validate_and_repair(BuildingBlueprint *blueprint,
                                             BuildingValidationReport *report);
uint64_t building_blueprint_fingerprint(const BuildingBlueprint *blueprint);

#endif
