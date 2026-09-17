# PBF particle-count scaling

Measured on 2026-09-13 with the `pbf-scaling` debug scenario at 1280x720.
Every timed frame includes a fixed PBF physics step and an offscreen convex-hull
rebuild/draw pass. Each run uses 120 steps, discards 30 warm-up steps, and averages
the remaining 90. These post-optimization values are one complete sweep. Fluid is
initialized as an `N x N x N` cube; each fluid cell has eight independent particles.
VSync and presentation are not included.

The GPU rows use the compute surface-net hull. The GPU report for each run also
separates timer-query compute/hull time from CPU packing, upload, and compact
readback time. The sweep now extends past 16K particles.

Backends:

- `gpu-gl43`: OpenGL 4.3 compute backend
- `cpu-mt`: persistent multithreaded CPU backend
- `cpu-st`: single-threaded CPU backend

Reproduce one point with:

```powershell
.\.build\bin\fps_ray.exe --debug-matrix --debug-scenario=pbf-scaling `
  --debug-pbf-side=8 --debug-steps=120 --debug-capture-steps=0 `
  --debug-output=benchmarks/pbf-scaling/side-8
```

The derived frame rate is `1000 / (averagePhysicsMs + averageRenderMs)`.
