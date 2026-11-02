# Reclaimed Foundry victory preview

The town montage is available from the results screen after completing all 12 Firefight waves. A three-second victory banner leads to results; Watch epilogue plays the trimmed song and credits full-screen. Hold B or Escape to skip back to results. Play again is controlled by the host in LAN games. Practice, endless and Enemy Lab sessions do not trigger the epilogue.

## Eight scenes

Each scene now lasts about 17.07 seconds, matching the supplied 2:16.5 soundtrack. The original eight-second animation timelines and one-second crossfades are uniformly stretched for playback.

1. **Restoring power** — mice turn the generator wheel and connect a cable outside the northwest bay library. Coils and lamps brighten; a reader has a book beside the shelves.
2. **Homes** — round doors, flower boxes, laundry and a mouse carrying produce.
3. **Market** — stalls below the west blast shelter, a balance scale, bread and a salvaged component exchanged between vendor and customer.
4. **Garden and workshop** — growing food, watering, repairs and cosmetic tether lifting.
5. **Playground** — a hollow robot head, an eye entrance, stairs, flowers and a rope swing suspended from a salvaged arm; three smaller mice play.
6. **Supper** — the original camp pot, table, cushions, food and four mice eating.
7. **Observatory** — a rooftop telescope, chart table, cushions and railing; two mice look up at a deterministic shooting star.
8. **Town under stars** — the wide arena ending, retaining the generated celestial sky and distant industrial silhouettes.

## Build and run

From the repository root with the existing raylib libraries:

```powershell
./tests/build_victory_town.ps1 -GL 43
./tests/build_victory_town.ps1 -GL 33
./.build/victory-town-gl43/victory-town.exe
```

The packages include shaders and MinGW runtime DLLs. Launch from the repository root or the package directory. Normal game executables are not replaced.

**1–8** select a scene; **0** selects the montage; **Space** pauses; **Left/Right** scrub; **B** toggles bloom; **Escape** exits.

## Reproduce media

```powershell
New-Item -ItemType Directory -Force artifacts/victory-town/frames,artifacts/victory-town/scenes,artifacts/victory-town/qa
./.build/victory-town-gl43/victory-town.exe --stills
./.build/victory-town-gl43/victory-town.exe --qa
./.build/victory-town-gl43/victory-town.exe --scenes
./.build/victory-town-gl43/victory-town.exe --capture
python tests/encode_victory_town.py
```

The harness writes `scenes.json` from the shared scene table. The encoder consumes this manifest rather than maintaining a separate name/order list.

- PNG captures: 1280×720, deterministic 20 fps.
- Eight individual GIFs: 768×432, about 17.07 seconds each.
- `reclaimed-foundry-eight-scenes.gif`: complete soundtrack-length montage.
- `eight-scene-overview.png`: all eight camera compositions.
- `--qa`: 32 frames, sampling each scene at 0, 3, 5.6 and 7.95 seconds.
- One-way power/market GIFs fade back to their opening over the final second. In the montage their completed actions hold while fading into the next scene.
- `--fallback` selects vertex-color town/prop materials; `--no-bloom` disables bloom. These combine with captures or benchmarks.

## Rendering and animation

`victory_town_scenes.inc` owns scene IDs, titles, cameras, durations and one-way flags. `victory_town_stories.inc` owns cached story geometry and explicit pose/prop sampling. `victory_town.inc` renders the shared town.

The actual Foundry's exposed surfaces are cached at initialization, then the simulation world is released. The cache contains 173,142 triangles, nine world material groups and five moving prop models. Static shadows are baked once, excluding moving props. Model and normal transforms are applied consistently to dynamic prop lighting. Materials, meshes and shadow targets are released on shutdown.

There are eleven background mice and at most three additional scene actors. The new actors use scaled local poses and explicit hand anchors. Swing seat and rope endpoints share an analytic pendulum sample; the climber resolves feet against stair heights. Market objects have explicit ownership and a common transfer position. Telescope rotation, generator power and meteor state are sampled once for both rendering passes. Generator power affects only its library area. Moving props use the same material lighting as the town; bloom-off and material fallback remain readable.

No runtime mesh generation, simulated constraints, enemy logic or terrain updates are used. The harness asserts the simulation world stays empty throughout rendering.

## Validation

The harness checks finite world/prop geometry, scene camera clearance, deterministic story samples, the 16-mouse bound, fixed swing rope length, market transfer attachment and scheduled meteor visibility. QA captures cover scene starts, interactions, meteor timing and endings. Benchmarks include crossfades and GPU completion but exclude PNG export.

Character mesh, presentation-lighting isolation and camp regressions have passed. Measurements for the final preview are recorded below after capture validation. Linux and macOS are not tested.

Measured at 1280×720 with bloom on Windows, Quadro RTX 4000 / i7-2600 (1,240 measured frames after 40 warm-up frames):

| Backend | Median | p95 |
|---|---:|---:|
| OpenGL 4.3 | 7.373 ms | 12.883 ms |
| OpenGL 3.3 | 7.603 ms | 13.257 ms |

Both meet the 16.7 ms median target. These are standalone preview render costs, not match frame times. Eight repeated still captures matched byte-for-byte. OpenGL 3.3 vertex-color fallback, OpenGL 4.3 bloom-disabled captures, character fallback, character lighting isolation and camp checks passed. Capture output contains 1,280 frames per full sequence. Non-Windows platforms remain untested.

## Soundtrack and stretched playback

The supplied **pretty good.mp3**, with its first 40 seconds removed, is packaged at `assets/audio/pretty-good.mp3`. Its duration is **136.5355 seconds (2:16.5)**. Interactive playback and the music video stretch the original 64 animation seconds by **2.13337×**, giving each scene **17.06694 seconds**, including its stretched transition. The music retains its original speed and pitch.

Live playback follows the streamed audio clock. Space pauses/resumes both; arrows seek both; selecting a scene seeks to its matching section of the song. The complete montage repeats with the song. If an audio device or the track is unavailable, the preview keeps the same stretched timing silently. Music resources are released on exit.

The original 20 fps PNG sequences remain compact authored-time samples for QA and GIF generation. The manifest now records `playback_seconds`; GIFs distribute their delay rounding to match the stretched duration (they remain silent). The MP4 freshly renders at 20 fps of actual playback time, rather than duplicating the old frames.

To export the music video, install FFmpeg and run from the repository root:

```powershell
$env:TOWN_FFMPEG='C:/path/to/ffmpeg.exe' # optional if ffmpeg is on PATH
./.build/victory-town-gl43/victory-town.exe --video
```

Output: `artifacts/victory-town/reclaimed-foundry-credits-trimmed.mp4`, 1280×720 H.264 with AAC audio. This export streams frames directly to FFmpeg without generating thousands of additional PNGs. `--audio-test` checks the loaded duration, pause, seeking, resume and cleanup.

## Credits overlay

The preview and newly rendered exports include five fading lower-left credit cards: Game by Michael; Music by Michael; Slop by chat; the Tim McGraw physics-engine attribution and work title; and all six play-testing names supplied by the creator. The overlay uses the same authored timeline as the soundtrack, so pause and scrubbing keep it synchronized. It is composited after bloom and scene crossfades for readable text with a subtle dark backdrop.

`--credits-stills` captures one frame per credit card to `artifacts/victory-town/credit-0.png` through `credit-4.png`. The soundtrack video with credits is exported by `--video` to `artifacts/victory-town/reclaimed-foundry-credits-trimmed.mp4`. Older exported GIFs/videos are historical captures; render them again to include the overlay.

The soundtrack trim removes exactly 1,764,000 stereo sample frames at 44,100 Hz (40 seconds). The remaining 6,021,216 frames play for 136.53551 seconds; the montage and credits scale to that duration. The original uploaded song is unchanged.

## In-game integration

The epilogue keeps the authoritative match in GAMEOVER and freezes gameplay input. Presentation is local: each LAN participant can watch or skip independently. A host restart cancels any client playback. Results are retained across replays; the completion timer counts active play on the local machine (clients label this Time in run). Revive totals are shown offline and on the host; they are omitted on clients because the existing protocol does not replicate them.

Scenery uses `assets/victory/foundry-shell.bin`, an immutable mesh snapshot, rather than rebuilding or replacing the live voxel world. Build scripts package this asset and the soundtrack. Missing scenery returns to results with an error; missing audio uses the same timed sequence silently. The banner uses a frozen final gameplay frame.

To regenerate the shell after intentionally changing the source Foundry design, build the standalone preview, set `FPS_EXPORT_TOWN_SHELL=1`, then run it with `--credits-stills` from the repository root. Clear the environment variable afterward. The binary contains magic/version `TOW1`, a vertex count, positions, normals and colors; it targets the supported little-endian platforms.

Regression harnesses: `tests/firefight_victory_test.c` tests frozen world geometry, results, keyboard/controller skip, replay, natural completion, exclusions and simulated disconnect/restart handling. `tests/firefight_victory_lan.c` runs as separate host/client processes (argument `host` for host, no argument for client), using port 27943, to test victory replication, independent playback and authoritative restart. Compile with the same renderer/network sources and libraries as the game, without NDEBUG so assertions remain enabled.
