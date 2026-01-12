# CPU AI Refactor (Constructor-Brawler Update)

The current CPU AI logic in `fps_ray.c` (`UpdateBot`) is designed for a traditional Health/Ammo pickup loop. It needs to be overhauled to support the new Matter-based economy and "Exposed" death mechanics.

## 1. Resource Management (Harvesting)
- [ ] **Deprecate Pickup Seeking:** Remove the loop that searches for `pickups[]` array (Ammo/Health boxes).
- [ ] **Implement "Harvest" State:**
    - **Trigger:** When `bot->matter` is low (e.g., < 25% or < `MATTER_SHOT_COST * 5`).
    - **Search:** Scan for the nearest active **Static Voxel** (Environment).
    - **Pathing:** Move within `MELEE_RANGE` (2.0 units) of the target voxel.
    - **Action:** Synthesize Melee Input (`KEY_Z` equivalent) to destroy the voxel and gain +10 Matter.

## 2. Combat Logic Upgrade
The bot must distinguish between "Shielded" enemies (Matter > 0) and "Exposed" enemies (Matter == 0).

### Phase A: Shield Breaking (Target Matter > 0)
- [ ] **Behavior:** Maintain medium range.
- [ ] **Action:** Use `FireVoxel` (Shoot) to deplete enemy Matter.
- [ ] **Constraint:** Stop shooting if `bot->matter` is critical (save reserve for Panic Builds or Shield buffer).

### Phase B: Finisher (Target `isExposed` == true)
- [ ] **Behavior:** Aggressive closing of distance. Shooting deals 0 damage now, so behavior must change.
- [ ] **Melee Kill:** Move to close range (< 2.0 units) and trigger Melee to knockback/kill.
- [ ] **Physics Kill (Hard Bots):**
    - Identify loose debris or rip a chunk from the wall using Gravity Tether (`KEY_R`).
    - Aim and throw object at the exposed player.

## 3. Self-Preservation (Exposed State)
- [ ] **Trigger:** High priority override when `bot->isExposed` is true.
- [ ] **Flee:** Move directly away from nearest enemy.
- [ ] **Panic Build:** Occasionally trigger `perform_build` while fleeing to place blocks behind (breaking line of sight).
- [ ] **Emergency Harvest:** Aggressively seek any nearby voxel to Melee harvest, restoring Matter to > 0 to regain "Shields".

## 4. Architecture Refactor
- [ ] **State Machine / Utility System:** Refactor the monolithic `UpdateBot` function into weighted utility functions:
    - `CalculateUtility_Harvest()`
    - `CalculateUtility_Combat()`
    - `CalculateUtility_Flee()`
- [ ] **Input Synthesis:** Ensure the bot function can simulate the new key inputs:
    - Melee (`KEY_Z` / `KEY_M`)
    - Build (`KEY_E` / `KEY_O`)
    - Tether (`KEY_R` / `KEY_P`)
