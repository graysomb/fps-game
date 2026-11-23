#!/bin/sh
# Pass the compiled tests/foundry_trap_preview binary and an optional destination.
set -eu
preview=${1:?usage: capture_trap_lab.sh PREVIEW_BINARY [OUTPUT_DIRECTORY]}
out=${2:-artifacts/traplab}
mkdir -p "$out"
frames=$(mktemp -d)
trap 'rm -rf "$frames"' EXIT
for scene in logs swing pit tower forest; do
    mkdir "$frames/$scene"
    FPS_TRAP_LAB="$scene" FPS_TRAP_FRAMES="$frames/$scene" "$preview" --physics=cpu-st
    ffmpeg -v error -y -framerate 20 -i "$frames/$scene/frame-%03d.png" \
        -filter_complex 'split[a][b];[a]palettegen=reserve_transparent=0[p];[b][p]paletteuse=dither=bayer:bayer_scale=4' \
        -gifflags 0 -loop 0 "$frames/$scene.gif"
    mv "$frames/$scene.gif" "$out/$scene.gif"
done
