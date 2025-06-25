## Octree Meshing Integration Tasks

| Step | Description                                                  | Status    |
|:----:|--------------------------------------------------------------|:---------:|
| 0    | Add `GreedyMeshChunk` & `OctreeNode` structs                  | ✅ Done   |
| 1    | Add `#define CHUNK_SIZE 16` & `#define MAX_LEVELS 6`         | ✅ Done   |
| 2    | Build phase in `ResetGame()` (BuildOctree/InsertVoxel/Remesh) | ✅ Done   |
| 3    | Greedy‑mesh helpers (`BuildGreedyMesh`/`GreedyRemeshLeaf`)     | ✅ Done   |
| 4    | Dirty‑flag propagation & `FlushDirty()` pass                  | ✅ Done   |
| 5    | Frustum‑culled drawing (`DrawOctreeNode`)                     | ✅ Done   |
| 6    | Camera helper (`GetCurrentCamera`)                            | ✅ Done   |
| 7    | Utility routines (`BuildOctree`, `LocateLeaf`, `InsertVoxel`) | ✅ Done   |
| 8    | Build command unchanged                                       | (n/a)     |
| 9    | Test plan & milestones                                        | (n/a)     |