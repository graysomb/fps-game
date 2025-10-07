# Project Tasks

This file tracks the current and future development tasks for the voxel FPS game.

## Gameplay Features

- [ ] Implement more weapon types (e.g., shotgun, grenade launcher).
- [ ] Add power-ups (e.g., speed boost, invincibility, quad damage).
- [ ] Implement more game modes (e.g., Capture the Flag, King of the Hill, Team Deathmatch).
- [ ] Add AI-controlled bots to play against.
- [ ] sliding, grappling hooks
- [ ] tune jump height with k/d
- [ ] what other stuff can be tuned?
- [ ] add sheilds
- [ ] health is add to the winner not the loser



## Graphics and Rendering

- [ ] Implement a skybox with a proper texture instead of a solid color.
- [ ] Improve voxel texturing and lighting.
- [ ] Add particle effects for explosions, bullet impacts, and destruction.
- [ ] Add post-processing effects (e.g., bloom, motion blur).

## Audio

- [ ] Add sound effects for shooting, jumping, walking, and voxel collisions.
- [ ] Add background music and a simple audio mixer.

## UI/UX

- [ ] Create a main menu with options to start the game and exit.
- [ ] Implement a pause menu.
- [ ] Add support for gamepads.
- [ ] Improve the in-game HUD to be more visually appealing.

## Code and Architecture

- [ ] Refactor the code to improve organization (e.g., separate files for player, voxel physics, rendering).
- [ ] Add more comments to the code to explain complex parts.
- [ ] Set up a build system (e.g., Makefile, CMake) to simplify compilation.

## Networking

- [ ] Implement basic networking for online multiplayer.
- [ ] ENet or make my own?
- [ ] Add a server browser and lobby system.

## Fractal Mode
- [ ] Implement modular fractal world where each node is essentially this code
- [ ] dirty chunking the greedymeshing with modular greedy meshing
- [ ] Implement basic networking for online multiplayer.

## Level building
- [ ] whats the vibe
- [ ] biomes?
- [ ] fractal generation
- [ ] super position

## physics
- [ ] add ray casting to player collisions
- [ ] Implement rigid body physics (floating bodies fall)
- [ ] Implement torque
- [ ] Implement intervoxel sticking

i treid octree it didn't really work, now im tryng the neighbor cascade check

# voxel types
- [ ] goop you move slowly throught
- [ ] bouncy
- [ ] slippy
- [ ] portal