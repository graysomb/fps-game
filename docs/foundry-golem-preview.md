# Foundry golem presentation prototype

The approved [concept](concepts/foundry-golem.png) is implemented as a cached,
articulated voxel model with blue-gray armor, orange shoulder/forearm plates,
dark mechanical joints, and a gold reactor behind a removable chest plate.
The preview shows the existing rat beside it for scale. It deliberately stays
outside normal matches: no Foundry AI, damage rules, physical debris, regeneration
simulation, or network changes are included in this phase.

## Run and capture

Run from the repository root on Windows (GCC and the existing raylib builds):

```powershell
./tools/preview_golem.ps1
./tools/preview_golem.ps1 -Mode checks -Backend 33 -Output artifacts/golem/gl33
./tools/preview_golem.ps1 -Mode checks -Backend 43
./tools/preview_golem.ps1 -Mode fallback
./tools/preview_golem.ps1 -Mode capture
python tools/golem_media.py
./tools/preview_golem.ps1 -Mode benchmark
```

The standalone executables are `.build/golem/preview-gl33.exe` and
`preview-gl43.exe`. Outputs go to `artifacts/golem/rendered` by default.
Captures comprise front/side/rear/three-quarter PNGs, a skeleton view,
bloom-disabled view, four-view montage, a concept comparison sheet, and nine
six-second, 20 fps GIFs: walking, slam, sweep, forearm throw, arm loss, leg loss,
chest loss, both-leg loss, and localized hit response. Each GIF starts from a
fresh deterministic timeline. Raw PNGs and GIFs are generated artifacts and are
not tracked; the source concept and capture tools are tracked.

| Control | Action |
| --- | --- |
| Tab | Next scripted scenario |
| Space | Pause/resume |
| Left / Right | Scrub time |
| Comma / Period | Step backward/forward by 0.05 s |
| Right mouse drag | Orbit |
| 1 / 2 | Toggle absent left/right arm |
| 3 / 4 | Toggle absent left/right lower leg |
| 5 | Toggle absent chest plate |
| R | Reset time and manual damage |
| S | Full model / skeleton |
| B / E | Toggle bloom / cosmetic effects |

## Presentation contract

Link `golem_renderer.c`, with raylib and the existing renderer/bloom interfaces.
`golem_renderer_init()` caches 18 rigid-group meshes, **283 cubes / 3,396
triangles**, below the 384-cube / 5,000-triangle budget. The 16-joint rig has root,
pelvis, torso, head, and three joints for each arm and leg. Mesh bounds are cached
alongside geometry; no meshes or entities are allocated during sampling/drawing.

Start with `golem_input()`. Positions and optional ankle contacts are world-space;
yaw and foot yaw are degrees. `travel` is accumulated distance in unscaled model
units, with a 1.2-unit gait cycle. The preview uses scale 1.5, approximately four
rat heights. Optional resolved ankle contacts hold the feet fixed through stance
on its curved walking path. Ground/contact selection belongs to the caller.

Call `golem_sample()` once per simulation/presentation frame, then share the
returned pose across all viewports and passes. Draw `WEAPON_COLOR` into the world
target, `WEAPON_EMISSION` using the existing depth-aware weapon bloom target, and
`WEAPON_GLASS` after opaque geometry. Glass outlines test world depth without
writing depth. `golem_draw_skeleton()` supplies the inspection overlay. The
renderer only consumes caller-owned state and never queries gameplay.

Capture a limb event using `golem_capture_damage()` **before** removing it. Keep
the returned matrices, yaw and floor height with the event. Advance its age and
provide normalized regeneration/recovery progress. Set `forearmOnly` for the
throw pose; ordinary arm loss detaches from the shoulder. Lower-leg loss removes
the shin and foot, leaving the thigh connected to the pelvis. Detached groups
use a bounded, scripted trajectory and floor clamp, not physics. Their lifetime
is four seconds, including a one-second fade.

Animation timing defaults are presentation tuning:

| Animation | Duration / behavior |
| --- | --- |
| Slam | 2.2 s; overhead windup, deep crouch, planted fist, impact ring |
| Sweep | 1.8 s; torso twist and side-to-side arm sweep |
| Throw | 2.4 s; windup and release of the forearm/hand group |
| Hit / arm stagger | 0.6 s; localized warm flash and recoil |
| Knee transition | 0.6 s; available hand braces, both knees fold if needed |
| Regeneration | 2 s; cyan target edges, staggered cube assembly into the pose |
| Standing recovery | 0.8 s after the lower leg is fully reconstructed |

The preview waits 1.25 s after damage before rebuilding. Gameplay will own that
decision later. Losing the attacking arm interrupts its attack; losing a leg
suppresses attacks and walking. Rebuilding lower legs are folded above ground
before standing recovery. Chest removal exposes the pulsing reactor; it does
not define future vulnerability or damage rules.

`shaders/golem.vert` carries per-cube assembly rank/center and flat normals;
`golem.frag` provides matte lighting, square face edges, hit flashes, gold emission,
and cyan ghost outlines. Solid/emission passes use identical transforms. Missing
shaders fall back to vertex colors (no ghost outlines or bloom), preserving the
animated model. Call shutdown before closing the GL context; repeated init and
shutdown are supported. Existing build scripts already copy the entire shader
directory, so the new assets follow Windows, shell, and macOS packaging paths.

## Validation

The preview checks deterministic replay over 6,000 samples, finite transforms,
unit normals, geometry budgets, arm/leg joint connectivity, stance contacts,
attack interruption, all 32 combinations of absent sections, recovery, idempotent
initialization and resource cleanup. A rendered wall comparison checks opaque,
emissive, and ghost occlusion together. Four-view captures use 640 x 360 tiles.
The checks run with both Windows OpenGL 3.3 and 4.3 and with deliberately missing
shader files. The slam hand reaches the floor; the kneeling brace is within
4 cm at the preview's enlarged scale.

The benchmark compares the same four 640 x 360 views, including bloom, with
cosmetic effects off/on during simultaneous lower-leg regeneration. It warms
60 frames and records 180 GPU-completed frames for each condition, repeated
three times. CSV contains individual frame times; capture/image encoding is
excluded. This measures the added ghost/debris effects, not the eventual cost
of physical voxel groups or boss gameplay.

Measured on Windows / NVIDIA Quadro RTX 4000 (driver 32.0.15.9671), OpenGL 4.3:

| Run | Effects off median | Effects on median | Added median | Off p95 | On p95 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.802 ms | 4.168 ms | 0.366 ms | 4.780 ms | 4.436 ms |
| 2 | 3.837 ms | 4.106 ms | 0.269 ms | 3.988 ms | 4.252 ms |
| 3 | 3.864 ms | 4.142 ms | 0.278 ms | 4.884 ms | 4.391 ms |

Pooled median: **3.834 ms off / 4.139 ms on**, an additional **0.305 ms**
across all four views. Pooled p95: 4.392 / 4.365 ms; the small negative delta
reflects timing noise, not a speedup from adding effects. The largest per-run
p95 increase was 0.264 ms. This meets the 0.5 ms median / 1 ms p95 target for
this preview workload. Uniform updates are cached, and fully opaque non-emissive
armor/debris skip redundant bloom-pass draws because world depth is already
copied into that pass.

Linux/macOS execution and gameplay integration are not validated in this phase.
Existing game executables and unrelated duck changes are preserved.
