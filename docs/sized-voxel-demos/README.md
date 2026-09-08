# Sized voxel CPU experiment

The legacy CPU solver now has an explicit-size voxel constructor and a bounded
one-level coarse group used by three debug fixtures. A double-edge cube is still
represented by eight ordinary unit voxels and their 27 shared particles. Its
eight outer particles are integrated independently; the other 19 are exact
trilinear samples of the coarse shape. Four additional unit voxels can attach
through the existing shared-pointer glue.

The parent is a constraint representation only. It adds no rendered voxel,
collision sample, or duplicate mass. Corrections on dependent child particles
are distributed back to the eight parent controls with the same interpolation
weights.

Run the validation with:

```sh
tools/build_sized_tests.sh
python3 tools/test_sized_harness.py .build/bin/fps_ray_sized_cpu --gif
```

The fixtures intentionally reject GPU backends. Automatic group construction,
partial group deletion, maps, networking, and editor spawning remain outside
this experiment.
