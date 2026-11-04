#!/bin/sh
# Build after ./build-macos.sh has prepared raylib and the gamepad object.
# Usage: sh tests/build_tether_profile_macos.sh [timing|diagnostics]
# Run from the repo root:
# FPS_TETHER_REALTIME=1 FPS_TETHER_ROUND=3 FPS_TETHER_TARGET=shield \
# FPS_TETHER_CYCLES=15 FPS_TETHER_QUIET=100 .build/wave7_tether_profile --physics=gpu
# Set FPS_TETHER_DIAGNOSTICS=1 with a diagnostics build for operation counters.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
test -f .build/raylib-macos/src/libraylib.a &&
test -f .build/raylib-macos/macos_gamepad.o || {
    echo "Run build-macos.sh first to prepare the macOS dependencies." >&2
    exit 2
}
case ${1:-timing} in
    timing) define=-DDEBRIS_BENCHMARK ;;
    diagnostics) define=-DDEBRIS_DIAGNOSTICS ;;
    *) echo "Expected timing or diagnostics" >&2; exit 2 ;;
esac
set -- tests/wave7_tether_profile.c weapon_renderer.c character_renderer.c duck_renderer.c \
    net_protocol.c net_transport.c physics_gpu_metal.m .build/raylib-macos/macos_gamepad.o
for source in callbacks compress host list packet peer protocol unix; do
    set -- "$@" "third_party/enet/$source.c"
done
clang "$@" "$define" -DFPS_GPU_METAL -std=c11 -O2 -g \
    -mmacosx-version-min="${MACOSX_DEPLOYMENT_TARGET:-12.0}" \
    -I. -Ithird_party/enet/include -I.build/raylib-macos/src \
    -L.build/raylib-macos/src -lraylib -pthread -lm \
    -framework Foundation -framework AppKit -framework IOKit -framework CoreVideo \
    -framework CoreGraphics -framework QuartzCore -framework OpenGL \
    -framework GameController -framework Metal -o .build/wave7_tether_profile
echo "Built .build/wave7_tether_profile (run from the repository root)."
