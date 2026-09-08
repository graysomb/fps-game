#!/bin/sh
# Separate development binaries: does not overwrite the packaged game.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
raylib_dir=${RAYLIB_DIR:-.build/raylib-macos/src}
mkdir -p .build/bin .build/sized-tests
binary_suffix=${BUILD_SUFFIX:-}
flags='-O2'
if [ "${SANITIZE:-0}" = 1 ]; then flags='-O1 -g -fsanitize=address,undefined'; fi
platform=''
if [ "$(uname -s)" = Darwin ]; then
    platform='-framework Foundation -framework AppKit -framework IOKit -framework CoreVideo -framework OpenGL'
fi
# Flags/source lists intentionally split; paths are quoted separately.
common='net_protocol.c net_transport.c third_party/enet/callbacks.c third_party/enet/compress.c third_party/enet/host.c third_party/enet/list.c third_party/enet/packet.c third_party/enet/peer.c third_party/enet/protocol.c third_party/enet/unix.c'
${CC:-clang} -std=c11 $flags fps_ray.c $common -I. -I"$raylib_dir" -Ithird_party/enet/include -L"$raylib_dir" -lraylib -pthread -lm $platform -o .build/bin/fps_ray_sized_cpu"$binary_suffix"
${CC:-clang} -std=c11 $flags tests/sized_physics_test.c $common -I. -I"$raylib_dir" -Ithird_party/enet/include -L"$raylib_dir" -lraylib -pthread -lm $platform -o .build/sized-tests/invariants"$binary_suffix"
.build/sized-tests/invariants"$binary_suffix"
