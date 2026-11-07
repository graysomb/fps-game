#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
raylib_source=${RAYLIB_SOURCE_DIR:-.build/raylib-macos}
clang -c -fobjc-arc -O1 -I. -I"$raylib_source/src" macos/macos_gamepad.m -o .build/firefight-save-gamepad.o
for test in firefight_cache_test firefight_sticky_drop_test firefight_save_test firefight_progression_test firefight_integration; do
    clang -std=c11 -O0 -g -DMACOS_GAMEPAD_IMPLEMENTATION -Wno-absolute-value -I. -I"$raylib_source/src" -Ithird_party/enet/include \
        "tests/$test.c" weapon_renderer.c character_renderer.c duck_renderer.c .build/firefight-save-gamepad.o \
        net_protocol.c net_transport.c third_party/enet/callbacks.c third_party/enet/compress.c \
        third_party/enet/host.c third_party/enet/list.c third_party/enet/packet.c third_party/enet/peer.c \
        third_party/enet/protocol.c third_party/enet/unix.c "$raylib_source/src/libraylib.a" -pthread -lm \
        -framework Foundation -framework AppKit -framework IOKit -framework CoreVideo -framework CoreGraphics \
        -framework QuartzCore -framework OpenGL -framework GameController -o ".build/$test"
done
.build/firefight_cache_test
.build/firefight_sticky_drop_test
.build/firefight_save_test
.build/firefight_save_test --read
.build/firefight_progression_test
.build/firefight_integration
