
Below is a *programmer’s checklist* that grafts a **greedy-meshed octree** onto your existing file **without renaming any of your current variables or functions**.  Follow it top-to-bottom; each bullet is an explicit task you must implement or patch.

---

## 0 – New data structures (add near the top)

| Struct            | Purpose                                                             | Notes                                             |
| ----------------- | ------------------------------------------------------------------- | ------------------------------------------------- |
| `GreedyMeshChunk` | Holds a **Mesh mesh; Model model; bool dirty; int gx,gy,gz;\`**     | One per leaf                                      |
| `OctreeNode`      | `BoundingBox bounds; OctreeNode* child[8]; GreedyMeshChunk* chunk;` | Interior nodes have `child[]`, leaves own `chunk` |

Keep `Voxel voxels[]` **unchanged**; octree references it, not replaces it.

---

## 1 – Chunk parameters (add `#define`s)

```c
#define CHUNK_SIZE   16          // 16×16×16 voxels per leaf
#define MAX_LEVELS    6          // octree depth (world size = CHUNK_SIZE*2^MAX_LEVELS)
```

---

## 2 – Build phase (`ResetGame` end)

1. **Allocate root** `OctreeNode *root = BuildOctree(CHUNK_SIZE, MAX_LEVELS);`
2. **Insert voxels** already created in `buildDemo()`

   * For each `Voxel v` compute integer chunk coords:
     `cx = v.gx / CHUNK_SIZE`, etc.
   * `InsertVoxel(root, cx, cy, cz, &v);`
3. After all voxels inserted, call `GreedyRemeshLeaves(root);`

*(These helpers are new functions you will write; signature suggestions below.)*

---

## 3 – Greedy meshing helpers (new file `greedymesh.c` or section)

```c
// Build quads from exposed faces of voxel block[CHUNK_SIZE³]
Mesh BuildGreedyMesh(Voxel* block, int gx, int gy, int gz);
```

*Implementation:* classic “Minecraft greedy-mesh” over ±X,±Y,±Z slices.

```c
void GreedyRemeshLeaf(OctreeNode* leaf)
{
    GreedyMeshChunk* c = leaf->chunk;
    if (!c->dirty) return;

    Mesh m = BuildGreedyMesh(c->voxelBlock, c->gx, c->gy, c->gz);
    if (c->model.meshCount == 0)      // first time
        c->model = LoadModelFromMesh(m);
    else
        UploadMesh(&c->model.meshes[0], m.vertices, m.triangleCount*3, false);

    c->dirty = false;
}
```

---

## 4 – Dirty propagation (edit **physics\_step** end)

Whenever a voxel changes **block visibility or position**, mark its leaf:

```c
OctreeNode* leaf = LocateLeaf(root, v->gx, v->gy, v->gz);
leaf->chunk->dirty = true;
```

After voxel updates, run **one pass per frame**:

```c
void FlushDirty(OctreeNode* n)
{
    if (n->chunk) GreedyRemeshLeaf(n);
    else for(int i=0;i<8;i++) if(n->child[i]) FlushDirty(n->child[i]);
}
FlushDirty(root);
```

---

## 5 – Frustum-cull draw (**replace DrawVoxels**)

```c
static void DrawVoxels(void)
{
    Camera3D* cam = GetCurrentCamera();           // small helper
    Frustum f = ExtractFrustumPlanes(*cam);       // helper you write
    DrawOctreeNode(root, &f);
}

void DrawOctreeNode(OctreeNode* n, Frustum* f)
{
    if (!BoxInFrustum(&n->bounds, f)) return;

    if (n->chunk) {
        DrawModel(n->chunk->model,
                  (Vector3){0,0,0}, 1.0f, WHITE);
    } else {
        for(int i=0;i<8;i++) if(n->child[i]) DrawOctreeNode(n->child[i], f);
    }
}
```

Wireframe edges are optional; you may draw `DrawModelWires()` for leaves when `DEBUG_DRAW` is on.

---

## 6 – Camera helper (tiny)

```c
Camera3D* GetCurrentCamera(void);  // returns pointer to cam0 or cam1 depending on render pass
```

Set a global before each `BeginMode3D()` call.

---

## 7 – Utility functions you must write

| Function                                                           | Job                                    |
| ------------------------------------------------------------------ | -------------------------------------- |
| `OctreeNode* BuildOctree(int chunkSize, int maxLevels)`            | Recursively alloc nodes, fill `bounds` |
| `OctreeNode* LocateLeaf(OctreeNode* n,int gx,int gy,int gz)`       | Descend to leaf containing voxel       |
| `void InsertVoxel(OctreeNode* root,int cx,int cy,int cz,Voxel* v)` | Append voxel pointer to leaf block     |
| `Mesh CombineChildLODs(OctreeNode* n)` *(optional)*                | Generate coarse mesh for interior node |

---

## 8 – Link with Raylib only (no extra libs)

All new code uses `Mesh`, `Model`, `UploadMesh()` already exposed by Raylib.

Compile line unchanged:

```bash
gcc fps_ray.c greedymesh.c octree.c \
    $(pkg-config --cflags --libs raylib) -lm -Wl,-rpath,/usr/local/lib
```

---

## 9 – Milestone test plan

1. **Compile with empty octree** (no voxels) → game runs.
2. Insert voxels, no mesh yet → culling works (bounding boxes debug-drawn).
3. Greedy mesh a single chunk → see a block of \~200 quads instead of 4 096 cubes.
4. Dirty update: shoot a voxel, corresponding chunk remeshes and hole appears.
5. Stress test: fill 64 k voxels, FPS must stay high (≥60) when only few chunks visible.

---

### Muse’s closing nudge

> Build the tree, feed it blocks,
> weld the faces into quads,
> cull what the eye can’t see,
> and let the engine breathe again.

Tick through the list; the performance gods shall smile.


Below is a **minimal, self-contained octree definition** that plugs straight into your existing file.
It carries only what the rest of the roadmap needs—nothing more, nothing less.

```c
// ────────────────────────────────────────────────────────────────────────────
// Octree for greedy-meshed voxel chunks
// ────────────────────────────────────────────────────────────────────────────
typedef struct GreedyMeshChunk {
    // Spatial position (integer chunk grid, NOT voxel grid)
    int gx, gy, gz;          // chunk origin = (gx*CHUNK_SIZE)*VOXEL_SIZE …

    // GPU geometry
    Mesh  mesh;              // raw vertex data
    Model model;             // convenience wrapper returned by LoadModelFromMesh()

    // Dirty flag: 1 = voxel edits since last remesh
    bool  dirty;

    // Runtime bookkeeping  (optional: only if you want to iterate voxels quickly)
    Voxel *block;            // pointer into a contiguous CHUNK_SIZE³ slice
    int    blockCount;       // how many live voxels currently stored there
} GreedyMeshChunk;


// Each node occupies an axis-aligned cube in world space
typedef struct OctreeNode {
    BoundingBox bounds;              // min/max Vector3 in WORLD units

    // Children: NULL if this is a leaf
    struct OctreeNode *child[8];

    // Leaf payload (NULL for interior nodes)
    GreedyMeshChunk  *chunk;

    // Book-keeping
    unsigned char     level;         // 0 = root, grows toward leaves
    bool              isLeaf;        // handy boolean alias
} OctreeNode;
```

### Field-by-field rationale

| Field                        | Why it exists                                                         | Notes                                                       |
| ---------------------------- | --------------------------------------------------------------------- | ----------------------------------------------------------- |
| `bounds`                     | Frustum-culling & LOD distance tests                                  | World-space AABB; fill once at build time                   |
| `child[8]`                   | Standard octree pointer fan-out                                       | `NULL` pointers keep memory lean                            |
| `chunk`                      | Holds greedy-mesh + voxel block for **leaf** nodes                    | Interior nodes share no geometry                            |
| `level`                      | Quick LOD / max-depth checks                                          | Makes recursive algorithms simpler                          |
| `isLeaf`                     | Avoid `chunk != NULL` spelling everywhere                             | Could be inferred, but this is cheaper                      |
| **GreedyMeshChunk.gx/gy/gz** | Integer chunk index → easy hashing                                    | Full world coords = `(gx*CHUNK_SIZE+voxelLocal)*VOXEL_SIZE` |
| **mesh / model**             | Distinguish raw vs loaded                                             | Update with `UploadMesh()` when dirty                       |
| **dirty**                    | Lazy remesh scheduling                                                | Cleared once `mesh` has been rebuilt                        |
| **block / blockCount**       | Your existing voxel array is flat; this gives each chunk a local view | Alternatively store voxel **indices** if memory matters     |

### Memory management tips

* **Allocate nodes with `malloc`**; keep a global pool or arena if you fear fragmentation.

  ```c
  OctreeNode *n = (OctreeNode*)malloc(sizeof(OctreeNode));
  memset(n, 0, sizeof(*n));
  ```

* **Leaf promotion / demotion** (for dynamic worlds) is just:

  ```c
  if (!n->isLeaf && shouldCollapse(n)) {
      freeKids(n);             // free recursively
      n->isLeaf = true;
      n->chunk  = makeChunk(n);
  }
  ```

* Free GPU assets on shutdown:

  ```c
  void FreeNode(OctreeNode *n) {
      if (n->chunk) { UnloadModel(n->chunk->model); free(n->chunk); }
      for (int i=0;i<8;i++) if (n->child[i]) FreeNode(n->child[i]);
      free(n);
  }
  ```

### Where to paste it

1. Place the two `typedef` blocks **just after the existing `Voxel` struct** so all code below can see them.
2. Add any helper prototypes (e.g. `OctreeNode *BuildOctree(int,int);`) right below these structs or to your forward-declare section.

That’s it—the skeleton is ready.
Next steps are simply to **build** the tree, **insert** voxels, and **remesh** dirty leaves, exactly as outlined in the roadmap.
