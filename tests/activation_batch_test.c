#define ACTIVATION_PROFILE
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Activation regression");
    activePlayers=0;clear_world_voxels();
    static UnitVoxelBuffer buffer;
    remove_buffered_static_voxels(NULL);
    remove_buffered_static_voxels(&buffer);
    assert(staticBeliefBatchDepth==0);
    assert(!activate_static_voxel_for_tether(-1,0,ACTIVATION_TETHER_BELIEF));
    assert(!activate_static_voxels_near_region(1,0,0,0,0,0,0));
    assert(staticBeliefBatchDepth==0);

    // Removal of the first voxel compacts a distant static voxel into its slot.
    int first=add_static_voxel_at_grid(-32,10,0,GRAY,0);
    add_static_voxel_at_grid(96,10,0,BLUE,0);
    refresh_static_voxel_beliefs();
    memset(&activationProfile,0,sizeof(activationProfile));
    assert(unit_voxel_buffer_push(&buffer,-32,10,0,GRAY,0,false,first,0,0));
    ++staticBeliefBatchDepth; // An enclosing transaction must survive nesting.
    remove_buffered_static_voxels(&buffer);
    assert(staticBeliefBatchDepth==1);--staticBeliefBatchDepth;
    assert(activationProfile.columnScans==0 && staticBeliefsForceFullRefresh);
    assert(voxel_count==1 && voxels[0].gx==96);
    assert(recycleQueueCount==1 && buffer.voxels[0].restorationQueued);
    assert(buffer.voxels[0].restorationDeadlineFrame==recycleFrameCounter+RECYCLE_STATIC_RESTORE_DELAY);
    refresh_static_voxel_beliefs();
    assert(!staticBeliefsForceFullRefresh && staticBeliefDirtyCount==0);
    assert(activationProfile.supportRefreshes==1);
    // Ordinary single-voxel invalidation still takes the existing path.
    mark_static_beliefs_dirty_for_voxel(&voxels[0]);
    assert(activationProfile.columnScans==1);
    refresh_static_voxel_beliefs();
    assert(activationProfile.supportRefreshes==2);
    unit_voxel_buffer_clear(&buffer);remove_buffered_static_voxels(&buffer);
    assert(staticBeliefBatchDepth==0 && !staticBeliefsForceFullRefresh);
    clear_world_voxels();clear_all_chunks();CloseWindow();
    puts("Activation batching regression passed");return 0;
}
