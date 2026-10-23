# Weapon expansion gameplay validation

Validated on Windows, October 21, 2026. These results cover the gameplay integration,
not just the previously committed scripted visual demos.

## Behavior

`tests/weapon_expansion_gameplay.c` passes with both OpenGL 3.3 / CPU physics and
OpenGL 4.3 / GPU physics. Coverage includes:

- DMR four-hit lethality against 100 Matter; Glock exposure followed by execution;
  invulnerability, friendly-fire protection, terrain/armor occlusion and periodic rays.
- Glob projectile allocation failure without charging Matter, five-second refresh
  and expiry, blocked restoration, unconsumed Matter pickups, moving attachments
  and the three-splat-per-actor limit.
- Separate grenade inventories, cycling, automatic selection, in-flight gravity
  reservations, minimum/bounce fuse, death/revive preservation and resets.
- Gravity pulling owner and opponents, intact cover blocking actor pull, real
  dynamic-voxel motion on both backends, persistent AI pull, 32-unit per-frame and
  256-unit per-field activation limits.
- Atomic igloo placement, unobstructed doorway and chimney, failed-placement
  refunds, incarnation checks, capacity failure, and reinforced persistent cells.
  A fully activated 267-cell igloo remains standing after 30 seconds with physics,
  sleeping and recycling enabled. Static cells also survive an artificially expired
  ordinary construction lifetime. No shelter cells enter the regeneration queue.
- Two copies of every new pickup in Foundry and four in Crosscut, with valid ground
  clearance and symmetric placement.

The existing equipment suite, aiming suite and Firefight integration suite also
pass. Existing equipment LAN tests pass, including switching, charge/release,
rockets, sticky grenades, building, sword animation and free melee.

## Multiplayer

`tests/weapon_expansion_lan.c` uses two real ENet processes on loopback. Host and
client both pass. It observes Glob projectiles, DMR traces, gravity fields, igloo
construction, suppression, grenade selection, equipment reset and the final pickup
slot in the enlarged pickup table. Every observed packet fits the 1100-byte limit.

Protocol 24 rejects mismatched builds. Pickup snapshots are sent in 48-item batches;
new effects use independent bounded chunks with stale-tick rejection. New grenade
cycle commands are deduplicated on the authority. WAN latency/loss and multiple
physical computers were not exercised in this session.

## Rendering and builds

- OpenGL 3.3 and 4.3 visual checks pass: finite geometry, model budgets, connected
  shelter, open chimney, depleted emission, occlusion, effect expiry and renderer
  shutdown/reinitialization. The shared igloo surface mesh has 1080 triangles.
- First-person captures include Glob, Glock and DMR with paws. The aiming regression
  captures all guns at full and split-screen sizes. Four-view Foundry captures show
  equipped weapons, pickups, fields and grenade indicators using the real renderer
  with bloom. A separate capture shows the actual constructed igloo and rat for scale.
- Windows CPU and GPU release binaries are built with `-O2 -DNDEBUG`; shaders are
  copied to the existing package directory. Linux/macOS execution and release builds
  were unavailable here. The unchanged Windows launcher selects the updated binaries.

Reproducible commands:

```powershell
./tools/validate_weapon_expansion.ps1 -Backend 33
./tools/validate_weapon_expansion.ps1 -Backend 43 -Stress
./tools/preview_weapon_expansion.ps1 -Mode checks -Backend 33
./tools/preview_weapon_expansion.ps1 -Mode checks -Backend 43
```

Generated images are ignored by Git:
`artifacts/expansion-gameplay-{Glob,Glock,DMR,Igloo}.png`,
`artifacts/expansion-four-{0,1,2,3}.png`, and the existing aiming captures.

## Four-player cost and remaining performance limit

The fixed Foundry stress harness uses four 640×360 views, bloom, real equipped
weapons, GPU synchronization and GPU physics on the Quadro RTX 4000. Its executable
uses `-O1` with assertions; these are validation-harness timings, not a release-build
frame-rate promise. No other tests or compiles ran during the recorded measurement.

| Condition | Median | p95 | Worst measured after warm-up |
| --- | ---: | ---: | ---: |
| Before effects | 16.958 ms | 18.269 ms | 20.225 ms |
| Four fields plus repeated Glob/DMR fire | 30.293 ms | 94.718 ms | 130.891 ms |
| After expiry, cleanup and regeneration | 17.698 ms | 20.520 ms | 29.765 ms |

The stress test retains bounded pools, finite transforms and returns to its original
19,071-voxel count. After cleanup the median is 0.740 ms above baseline and p95 is 2.251 ms higher;
frame time recovers substantially, though this short run cannot establish a zero
residual performance cost. Four simultaneous destructive fields
are still expensive while active. Instrumented slow frames spent roughly 44–61 ms
in physics, with additional rendering and occasional DMR terrain-activation cost.
This remains a performance limitation rather than a claimed fix for all activation
hitching. The earlier sub-millisecond visual-only benchmark does not include these
simulation costs.

The harness prints phase medians/tails and per-stage timing for slow frames. Local
raw output for this run is `.build/expansion-clean-benchmark.log`; rerun with `-Stress`
to reproduce on another machine.
# Playtest tuning and fixes

The latest playtest pass makes Glob hits Matter-only, including repeated hits on
exposed actors. The Foundry core is left at a minimal positive value so Glob
cannot trigger its special death rule. Five-second restoration suppression is
unchanged. Attached globs are hidden from their victim's viewport in color and
emission passes; nearby glob geometry is also omitted to keep the camera clear.

Glock damage is now 8 (previously 15), switching is 0.12 seconds (previously
0.25), and Anchor healing is 16 Matter/second (previously 8). Weapon cooldowns,
Anchor exposure delay, and restoration suppression still apply.

Both 3x3x3 and upgraded 5x5x5 building allow actor overlap. Successful construction
lifts affected actors onto the roof; failed placement/allocation spends no Matter
and moves nobody. Existing terrain, loose debris and insufficient headroom still
block placement. Worn equipment does not reject an otherwise valid build.

Grenade supply pickups fill frag, gravity and igloo inventories to two each.
In Firefight only team-unlocked types refill, so early supplies cannot bypass
wave progression. Foundry retains two authored refill positions (near z +/-28),
and Crosscut retains four symmetric positions; each respawns after 30 seconds.
Carrier weapon drops and unlock celebrations are unchanged. Firefight's gold
tether pickup now renders the claw, matching the normal pickup path.

Gravity fields ignore dynamic debris/armor for actor-attraction visibility, while
intact static cover still blocks them. Actor acceleration is 45 (24 for Goliaths),
with the existing six-meter falloff and 12 m/s approach cap; bot pull momentum
decays at 0.5/s instead of 1.5/s. The Foundry remains immovable.

Reproduce focused gameplay and render checks with
`./tools/validate_playtest_fixes.ps1 -Backend 43` (or `33`). The test exercises
nonlethal projectile hits, exposed actors, the Foundry, switch cooldowns, Anchor
healing/suppression, overlapping builds, refill unlock gating, running bots,
armor and static occlusion, and unchanged victim-view pixels with an attached
glob. It writes `artifacts/playtest-tether-pickup.png`.

Validation for this pass: focused checks pass on OpenGL 3.3 and 4.3; expansion
gameplay and localhost host/client tests pass on both builds. GPU equipment,
sword, aiming, four-player equipment stress and equipment LAN regressions pass.
Firefight progression/carrier/checkpoint regressions pass. The tether pickup
capture was visually inspected, and both Windows game executables were rebuilt.
