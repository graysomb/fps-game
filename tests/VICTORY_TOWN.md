# Reclaimed Foundry — standalone victory-town preview

This is a presentation prototype, deliberately not connected to the Firefight results screen. The game map, enemies, combat, and networking are unchanged.

Four eight-second scenes form a 32-second loop, with one-second crossfades included in each scene's duration:

1. Four player-colored mice sharing stew on the central platform.
2. Round wooden doors, warm windows, hanging lights and laundry on the existing Foundry bays; reading and carrying mice.
3. Garden beds, workshop, cosmetic tether lifting and a swimming mallard in an irrigation pool.
4. A wide view of the reclaimed arena beneath the original generated celestial texture, with dusk fog and distant industrial silhouettes.

## Build and run

From the repository root, after the normal build has produced the raylib libraries:

```powershell
./tests/build_victory_town.ps1 -GL 43
./tests/build_victory_town.ps1 -GL 33
./.build/victory-town-gl43/victory-town.exe
```

The separate packages contain shaders and required MinGW runtime DLLs. Launch from the repository root or from the package directory. Normal game executables are not replaced.

Controls: **1–4** select a scene; **0** returns to montage; **Space** pauses; **Left/Right** scrub; **B** toggles bloom; **Escape** exits.

## Reproduce the media

```powershell
New-Item -ItemType Directory -Force artifacts/victory-town/frames,artifacts/victory-town/scenes,artifacts/victory-town/qa
./.build/victory-town-gl43/victory-town.exe --stills
./.build/victory-town-gl43/victory-town.exe --capture
./.build/victory-town-gl43/victory-town.exe --scenes
python tests/encode_victory_town.py
```

`--capture` exports the montage, `--scenes` exports four uninterrupted scene loops. Both use explicit 20 fps timestamps. PNGs are 1280×720; the five GIFs are 768×432. The encoder also generates an overview and, when the committed concept artifact is available, a concept/render comparison sheet.

`--benchmark` measures render submission plus GPU completion over 600 frames after 40 warm-up frames. It includes crossfade frames and excludes image readback/encoding. `--no-bloom` and `--fallback` can be combined with these modes. Fallback tests the new town material; mouse fallback has its own existing regression.

## Implementation

`victory_town.inc` exposes init/sample/draw/shutdown functions to the standalone harness. It snapshots exposed surfaces from the real deterministic Foundry generator, caches the geometry and props, and releases the simulation world before rendering. The resulting cache contains 158,942 triangles across nine explicit material groups. All ongoing movement is analytic presentation animation: eleven mice total, one duck, cloth/foliage sway, water ripples and two cosmetic tether beams. No AI, physics, voxel creation, terrain rebuilding or gameplay state updates run during the montage. The harness asserts the world remains empty each rendered frame.

Mouse body, feet and hand poses are sampled from the same explicit scene time. Held food, book, basket and tools follow hand anchors. Visors are sorted back-to-front. The camp renderer now offers `camp_draw_at(camera, position, time, emission)`; its existing gameplay wrapper preserves visibility, position, scale and timing. Town shaders provide seven local warm lights, material-specific wood grain/stone wear/cloth weave/metal highlights, emissive fixtures, sway, ripples and distance fog. A 2048�2048 directional shadow map is baked once from the static scenery, with filtered sampling during playback. Mice have layered soft contact shadows. The optional `character_presentation_light(tint, strength)` interface changes only the preview lighting and is reset after both render passes. The original sky generator is reused with a presentation-specific horizon transition.

## Detail pass

Homes now have thick masonry surrounds, inset planked doors, hinges and ring handles, mullioned glass with curtains and flower boxes, braced awnings, plank porches and rails, hanging lantern housings, crates, log piles and branching vines. The feast has a lantern canopy, patterned rug, cushions, provisions and crockery. The workshop has a tool board and salvaged robot mechanism. The garden has leaf clusters, flowers, produce, reeds and an irrigation outlet.

The mice use staggered spoon/bowl gestures, an animated book page, a carried produce basket with body-relative orientation, watering and repair tools, and hand-anchored tether beams. The carrier follows a continuous closed path with analytical planted stance intervals and lifted swing phases. The reader stands on the porch surface. Looping cameras have been reframed to keep the activities visible.

`--qa` exports twelve frames (start, middle and end of each shot) to `artifacts/victory-town/qa`. The harness validates finite cached positions/normals and camera clearance at nine times per scene. The encoder creates `before-after.png` if the previous overview is available in `before-detail`.

## Validation (Windows, Quadro RTX 4000 / i7-2600)

At 1280×720, bloom enabled:

| Preview | Median | p95 |
|---|---:|---:|
| OpenGL 4.3 | 6.493 ms | 11.137 ms |
| OpenGL 3.3 | 6.553 ms | 11.214 ms |

The previous preview measured approximately 5.6 ms median / 10 ms p95. The detailed version adds about 1 ms median render time. These are standalone rendering costs, not match performance or a promise for other machines. Both are below the 16.7 ms median target.

Camp timing, placement, render and fallback regression passes. Its fixtures now finish the final round before expecting the loot/camp phase, matching the current multi-round waves. Character mesh and material-fallback regressions pass; the existing grounding regression was also verified for the initial preview. The new character presentation-lighting regression verifies pixel-identical restoration of default lighting and unchanged emission output. A repeated fixed-time PNG capture matches byte-for-byte after seeding both cosmetic random generators. Captures were inspected for camera obstructions, reversed faces, scale, water/duck placement, hand props, windows and transparent visor ordering. The concept's warm domestic motifs are represented in the actual blocky Foundry rather than changing its layout to match the illustration.

Linux and macOS were not tested. Audio, match integration and result-screen UI are intentionally deferred.
