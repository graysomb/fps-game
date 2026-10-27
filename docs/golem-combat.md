# Enemy Lab golem combat

The golem remains an Enemy Lab encounter. Firefight and LAN boss integration are unchanged.
The physical body retains its original 283 cells, with eight additional cells for two
independently destructible shoulder-gun mounts.

## Shields and counterplay

Every attached physical segment has independent Matter. Hands and feet have 50;
arm/leg segments, eye and gun mounts have 75; head and pelvis have 100; torso,
chest and reactor have 150. The reactor also has 200 permanent health.

Shields wait eight seconds after their last hit, then restore five Matter per second.
Glob suppression lasts five seconds. Shielded segments block cutting, mining and
physical tether grabs. Tether locks onto the contacted part and drains at the normal
two-bullets-per-second rate. Release and grab again after depletion to rip it off.
The empty lock prevents recharge until released.

The hit that breaks a body shield cannot also cut its voxels. Shotgun pellets aggregate
per part; sword damage is deduplicated per swing. Rail shots stop at shielded parts;
exposed terrain retains normal penetration. Rockets evaluate protection before terrain
damage. Reactor damage still overflows into permanent health, which does not regenerate.

Detached parts have no shields. Regrown parts start empty and wait before recharging.
Losing a parent releases dependent parts; grouped physical repairs and kneeling remain.
Breaking legs opens the reactor and suspends ranged attacks.

## Attacks

- Slams: 40 Matter, one-second windup, four-second cooldown.
- Sweeps: 25 Matter, 0.8-second windup, three-second cooldown.
- Shoulder guns: existing default projectile construction and impact behavior,
  three-shot bursts at three shots/second, two-second gaps, alternating shoulders,
  30 m range and deterministic spread within four degrees.
- Head rockets: existing rocket pool, speed, lifetime, blast and terrain behavior.
  One-second warning, final 0.3 seconds of locked aim, 15-second cooldown, 8–35 m range.
- Right-hand structure throw: replaces the autonomous fist launch. Acquires an eligible
  visible structure within the player's upgraded-tether range, lifts for 1.5 seconds,
  and releases through the shared upgraded-tether implementation. Ten-second cooldown;
  failed acquisitions retry after two seconds. No golem-specific voxel cap or Lobber
  extraction shortcut. Defensive Matter is not weapon ammunition.

Arm/mount loss disables the corresponding gun; head loss disables rockets; right-arm
loss drops a held payload without a throw. Ranged attacks pause during melee, throws,
recovery and kneeling. Holograms remain valid targets. Failed projectile allocation
does not register a shot. Existing projectiles finish naturally after boss death.

The body, turret recoil, head aperture, part shields, target-part HUD and upgraded
tether beam are sampled from the same simulation state. Effects use existing cached
weapon materials and depth-tested geometry; no new full-screen pass or runtime mesh
generation is required.

## Implementation

golem_arsenal.inc owns part shields, drain locks, ranged timers and attack presentation.
The reactor proxy retains its position: explicit hand-anchor and release-aim hooks
adapt the shared player tether without teleporting the proxy. Default bullets share
a projectile-construction helper; rockets share the existing bounded pool.

Combat timers and event creation run outside GPU retry attempts. Authored golem welds
are preserved during unrelated world glue rebuilds, so grabbing terrain does not
silently detach the body. Existing cuts are preserved too.

## Controls and reproduction

F3 / D-pad Right cycles combat and manual demonstrations. F4 resets; F5 pauses;
F6 single-steps; F7 toggles god mode; F8 toggles repair; F9 shows the skeleton;
F10 removes an arm; F11 removes a shin; F12 removes the chest. The original fist launch
remains a manual demonstration.

- ./tools/validate_golem_combat.ps1 -Backend 43
- ./tools/validate_golem_combat.ps1 -Backend 33
- ./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -Fallback
- ./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -Preview
- python tools/golem_combat_media.py (Pillow required)
- ./tools/validate_golem_combat.ps1 -Backend 43 -SkipBuild -StressSeconds 600
- ./tools/validate_golem_physics.ps1 -Backend 33 -Repair

Captures are written to artifacts/golem/combat: part-shields.gif,
shoulder-guns.gif, head-rocket.gif, and structure-throw.gif.
Stress timing measures simulation/update in the lab fixture, excluding rendering.

## Validation (Windows, 2026-10-26)

- CPU and GPU combat checks passed: independent part shields, delayed regeneration,
  shield drain/release, shotgun aggregation, sword deduplication, rocket protection,
  projectile pool exhaustion, timed gun bursts, head loss and shared upgraded tether.
- Existing equipment regression and CPU grouped-repair regression passed, including
  allocation rollback, held debris, actor obstruction and damage during recovery.
- Forced GPU validation failure during a structure throw produced one throw and one
  clock advance, retaining the hand; CPU fallback completed the combat suite.
- The 600-second GPU fixture completed 45 slams, 32 sweeps, 459 turret shots,
  29 rockets and 30 structure throws. Peak voxel count was 575, final count 500;
  all checked particle positions remained finite.
- Whole fixture simulation/update timing: median 8.591 ms, p95 12.871 ms;
  first/last two-minute means 9.308/9.093 ms. Rendering is excluded, and preview
  compilation overlapped the early portion, so these are endurance-run observations,
  not a controlled before/after performance benchmark.
- Optimized Windows OpenGL 3.3 CPU and OpenGL 4.3 GPU builds were rebuilt.
  Other operating systems and LAN boss gameplay are outside this local-lab phase.
