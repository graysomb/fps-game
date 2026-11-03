#!/bin/sh
set -eu

project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
raylib_source=${RAYLIB_SOURCE_DIR:-"$project_root/.build/raylib-macos"}
deployment_target=${MACOSX_DEPLOYMENT_TARGET:-12.0}
app_dir="$project_root/.build/Reclaimed Foundry Credits.app"
bin_dir="$app_dir/Contents/MacOS"
resources="$app_dir/Contents/Resources"

if [ ! -f "$raylib_source/src/libraylib.a" ]; then
    echo "Build raylib first with ./build-macos.sh, or set RAYLIB_SOURCE_DIR to a built raylib tree." >&2
    exit 2
fi

mkdir -p "$bin_dir" "$resources/artifacts/victory-town/frames" \
    "$resources/artifacts/victory-town/scenes" "$resources/artifacts/victory-town/qa"
clang -c -fobjc-arc -O2 -mmacosx-version-min="$deployment_target" \
    -I"$project_root" -I"$raylib_source/src" \
    "$project_root/macos/macos_gamepad.m" -o "$bin_dir/macos_gamepad.o"

clang -std=c11 -O2 -mmacosx-version-min="$deployment_target" \
    -I"$project_root" -I"$raylib_source/src" -I"$project_root/third_party/enet/include" \
    "$project_root/tests/render_victory_town.c" \
    "$project_root/weapon_renderer.c" "$project_root/character_renderer.c" \
    "$project_root/duck_renderer.c" "$bin_dir/macos_gamepad.o" \
    "$project_root/net_protocol.c" "$project_root/net_transport.c" \
    "$project_root/third_party/enet/callbacks.c" "$project_root/third_party/enet/compress.c" \
    "$project_root/third_party/enet/host.c" "$project_root/third_party/enet/list.c" \
    "$project_root/third_party/enet/packet.c" "$project_root/third_party/enet/peer.c" \
    "$project_root/third_party/enet/protocol.c" "$project_root/third_party/enet/unix.c" \
    "$raylib_source/src/libraylib.a" -pthread -lm \
    -framework Foundation -framework AppKit -framework IOKit -framework CoreVideo \
    -framework CoreGraphics -framework QuartzCore -framework OpenGL -framework GameController \
    -o "$bin_dir/victory-town"
rm "$bin_dir/macos_gamepad.o"

cp -R "$project_root/shaders" "$resources/"
cp -R "$project_root/assets" "$resources/"
cat > "$bin_dir/credits" <<'LAUNCHER'
#!/bin/sh
set -eu
bin_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$bin_dir/../Resources"
exec "$bin_dir/victory-town" "$@"
LAUNCHER
chmod +x "$bin_dir/credits"
cat > "$app_dir/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
<key>CFBundleExecutable</key><string>credits</string>
<key>CFBundleIdentifier</key><string>local.fps-game.reclaimed-foundry-credits</string>
<key>CFBundleName</key><string>Reclaimed Foundry Credits</string>
<key>CFBundlePackageType</key><string>APPL</string>
<key>CFBundleVersion</key><string>1</string>
<key>NSHighResolutionCapable</key><true/>
</dict></plist>
PLIST
codesign --force --deep --sign - "$app_dir"
echo "Built $app_dir"
