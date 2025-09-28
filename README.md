# Split-Screen Voxel FPS Game

A simple split-screen first-person shooter game built with raylib. The game features a voxel-based world where players can destroy and create blocks. The game also features a Kill/Death ratio system that affects player stats.

## How to Compile and Run

To compile the game, you need to have raylib installed. Then run the following command:

```bash
gcc -g fps_ray.c -o fps_ray $(pkg-config --cflags --libs raylib) -lm -Wl,-rpath,/usr/local/lib
```

To run the game, execute the following command:

```bash
./fps_ray
```

## Controls

### Player 1

-   **Move:** WASD
-   **Turn:** F/H
-   **Look up/down:** T/G
-   **Shoot:** Left Control
-   **Toggle build/destroy mode:** Q

### Player 2

-   **Move:** I/K/J/L
-   **Turn:** Left/Right arrow keys
-   **Look up/down:** Up/Down arrow keys
-   **Shoot:** Right Control
-   **Toggle build/destroy mode:** U

## KD-Ratio System

The game features a Kill/Death (KD) ratio system that affects player stats.

-   The KD ratio is calculated as `(kills + 1) / (deaths + 1)`.
-   A higher KD ratio increases a player's movement speed, acceleration, and turning speed.
-   A lower KD ratio decreases these stats.
-   Player health upon respawning is also affected by the KD ratio. A higher KD ratio results in more health upon respawn.
