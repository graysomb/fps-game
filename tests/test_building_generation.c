#include "../building_generation.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

static BuildingVoxelSpec voxel(int x,int y,int z,uint16_t role){
    return (BuildingVoxelSpec){x,y,z,0xff808080u,0,role,-1,BUILDING_VOXEL_DAMAGEABLE};
}

int main(void){
    uint32_t a=building_seed_stream(42,BUILDING_STYLE_GREEK,2,1);
    uint32_t decoration=building_seed_stream(42,BUILDING_STYLE_GREEK,2,2);
    for(int i=0;i<100;i++)building_rng_next(&decoration);
    assert(a==building_seed_stream(42,BUILDING_STYLE_GREEK,2,1));

    BuildingBlueprint b={0};BuildingValidationReport report={0};
    assert(building_blueprint_init(&b,128,BUILDING_STYLE_GREEK,42,2,7));
    assert(building_blueprint_put(&b,voxel(0,0,0,BUILDING_ROLE_FOUNDATION),&report));
    assert(building_blueprint_put(&b,voxel(0,1,0,BUILDING_ROLE_COLUMN),&report));
    assert(building_blueprint_put(&b,voxel(8,4,0,BUILDING_ROLE_BEAM),&report));
    assert(building_blueprint_put(&b,voxel(8,4,0,BUILDING_ROLE_BEAM),&report));
    assert(report.duplicate_writes==1);
    assert(building_blueprint_validate_and_repair(&b,&report));
    assert(report.valid&&report.piers_added==4);
    assert(building_blueprint_find(&b,8,0,0));
    assert(building_blueprint_fingerprint(&b)!=0);
    building_blueprint_release(&b);
    puts("building generation: deterministic streams, transactional occupancy and structural repair passed");
    return 0;
}
