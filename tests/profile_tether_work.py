"""Generate a diagnostic-only copy for wave7_tether_profile.c.
Run from the repository root, then build the harness with
-DTETHER_PROFILE_INSTRUMENTED -I. in addition to its normal GPU build flags.
No production source is modified.
"""
from pathlib import Path
s = Path("fps_ray.c").read_text()
def replace(old, new):
    global s
    assert old in s, old
    s = s.replace(old, new, 1)
replace("static int recycleFrameCounter = 0;", """static atomic_ullong profileCollect,profileVisited,profileReject,profileRejectedCells;
static unsigned long long profileWakeCalls,profileWakeCells,profileFreezeCalls,profileFreezeCells;
static int recycleFrameCounter = 0;""")
replace("    int queue[VOXEL_ACTIVATION_UNIT_BUDGET];\n    int head = 0;",
        "    atomic_fetch_add(&profileCollect,1);\n    int queue[VOXEL_ACTIVATION_UNIT_BUDGET];\n    int head = 0;")
replace("        int current_idx = queue[head++];",
        "        atomic_fetch_add(&profileVisited,1);\n        int current_idx = queue[head++];")
replace("if (!dynamic_belief_overcomes_static(dynamic->activationBelief, clusterBelief)) {",
        """if (!dynamic_belief_overcomes_static(dynamic->activationBelief, clusterBelief)) {
                        atomic_fetch_add(&profileReject,1);
                        atomic_fetch_add(&profileRejectedCells,buffer->count-previousCount);""")
replace("static bool freeze_dynamic_cluster_in_place(const int *cluster, int cluster_count)\n{",
        "static bool freeze_dynamic_cluster_in_place(const int *cluster, int cluster_count)\n{\n    profileFreezeCalls++;")
replace("    int frozen_stamp = ++particle_sync_stamp;",
        "    profileFreezeCells += cluster_count;\n    int frozen_stamp = ++particle_sync_stamp;")
replace("static bool wake_sleeping_cluster_in_place(int start_idx)\n{",
        "static bool wake_sleeping_cluster_in_place(int start_idx)\n{\n    profileWakeCalls++;")
replace("        voxel->sleeping = false;\n        voxel->sleepFrames = 0;",
        "        profileWakeCells++;\n        voxel->sleeping = false;\n        voxel->sleepFrames = 0;")
Path(".build").mkdir(exist_ok=True)
Path(".build/fps_tether_instrumented.c").write_text(s)
print("Generated .build/fps_tether_instrumented.c")
