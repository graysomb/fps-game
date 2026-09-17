#!/bin/sh
# Isolated developer build: never rewrites the tracked application bundle.
set -eu
cd "$(dirname "$0")/.."
mkdir -p .build/adaptive
clang -std=c11 -O2 -Wall -Wextra -Werror -I. adaptive_mesh.c adaptive_metal.m \
    tests/adaptive_mesh_test.m -framework Foundation -framework Metal -o .build/adaptive/mesh_test
clang -std=c11 -O2 -Wall -Wextra -Werror -I. adaptive_mesh.c adaptive_metal.m adaptive_physics.m \
    tests/adaptive_physics_test.m -framework Foundation -framework Metal -o .build/adaptive/physics_test
clang -std=c11 -O2 -DNDEBUG -DFPS_GPU_METAL -I. -I.build/raylib-macos/src -Ithird_party/enet/include \
    fps_ray.c physics_gpu_metal.m adaptive_mesh.c adaptive_metal.m adaptive_physics.m \
    .build/raylib-macos/macos_gamepad.o net_protocol.c net_transport.c \
    third_party/enet/callbacks.c third_party/enet/compress.c third_party/enet/host.c \
    third_party/enet/list.c third_party/enet/packet.c third_party/enet/peer.c \
    third_party/enet/protocol.c third_party/enet/unix.c -L.build/raylib-macos/src -lraylib \
    -pthread -lm -framework Foundation -framework AppKit -framework IOKit \
    -framework CoreVideo -framework CoreGraphics -framework QuartzCore -framework OpenGL \
    -framework GameController -framework Metal -o .build/adaptive/fps_ray_gpu
