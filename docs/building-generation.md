# Building generation qualification

The five building grammars publish through `BuildingBlueprint`. A blueprint owns
voxel coordinates, material, source node, structural role, anchor flags, and a
structural group. Rasterization writes into this temporary object; voids are
applied there; validation and deterministic repair complete before any voxel is
added to the live world.

Generation uses independent integer streams derived from seed, style, family,
and purpose. Layout, dimensions, optional features, decoration, and Citadel
gameplay placement can therefore evolve without perturbing one another. Seed
compatibility with the old generators is intentionally not preserved.

The validator runs at most eight repair passes. It verifies six-connected paths
to foundations and checks vertical bearing for columns, walls, buttresses,
beams, and slabs. Repairs extend the originating member with visible piers; they
do not fill arbitrary enclosed volume. The complete sorted blueprint is
published transactionally. Structural group, source node, role, and anchor flags
survive activation and static restoration.

Run the deterministic generation and diversity checks with:

```sh
sh tools/test_building_generation.sh
```

Run a fully activated building against a physics backend with:

```sh
FPS_BUILDING_STYLE=greek FPS_BUILDING_SEED=1337 \
  .build/bin/fps_ray_cpu --physics=cpu-mt \
  --debug-scenario=building-stability --debug-steps=600
```

Valid styles are `greek`, `megalith`, `hyperborean`, `forerunner`, and
`citadel`. The scenario activates only the building's structural group, retains
explicit foundation and fantasy anchors, prevents test-time static recycling,
and checks fracture, anchor state, displacement, center-of-mass drift, and
settling energy.
