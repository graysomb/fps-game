#!/bin/sh
set -eu

project_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
raylib_source=${RAYLIB_SOURCE_DIR:-}
configuration=${1:-release}
architecture=${2:-native}
deployment_target=${MACOSX_DEPLOYMENT_TARGET:-12.0}

if [ -z "$raylib_source" ] || [ ! -d "$raylib_source/src" ]; then
    echo "Set RAYLIB_SOURCE_DIR to a raylib source tree" >&2
    exit 2
fi
if [ "$configuration" != release ] && [ "$configuration" != debug ]; then
    echo "Usage: ./build-macos.sh [release|debug] [native|universal]" >&2
    exit 2
fi
if [ "$architecture" = universal ]; then
    arch_flags="-arch arm64 -arch x86_64"
elif [ "$architecture" = native ]; then
    arch_flags="-arch $(uname -m)"
else
    echo "Unknown architecture '$architecture' (expected native or universal)" >&2
    exit 2
fi

command -v xcrun >/dev/null 2>&1 || { echo "Install the Xcode command-line tools" >&2; exit 2; }
command -v clang >/dev/null 2>&1 || { echo "clang is required" >&2; exit 2; }
if command -v python3 >/dev/null 2>&1; then
    python3 "$project_root/tools/validate_gpu_layout.py"
fi

build_root="$project_root/.build"
bin_dir="$build_root/bin"
raylib_build="$build_root/raylib-macos"
app_dir="$build_root/FPS Game.app"
mkdir -p "$bin_dir" "$raylib_build" "$app_dir/Contents/MacOS" "$app_dir/Contents/Resources"

if [ ! -f "$raylib_build/src/Makefile" ]; then
    cp -R "$raylib_source/." "$raylib_build/"
fi

make -C "$raylib_build/src" clean PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC
make -C "$raylib_build/src" PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC \
    CUSTOM_CFLAGS="$arch_flags -mmacosx-version-min=$deployment_target"

if [ "$configuration" = debug ]; then
    opt_flags="-O0 -g"
else
    opt_flags="-O2 -DNDEBUG"
fi

frameworks="-framework Foundation -framework AppKit -framework IOKit -framework CoreVideo -framework OpenGL"
common_flags="-std=c11 $opt_flags $arch_flags -mmacosx-version-min=$deployment_target -I$project_root -I$raylib_build/src -L$raylib_build/src -lraylib -pthread -lm $frameworks"

# shellcheck disable=SC2086
clang "$project_root/fps_ray.c" "$project_root/physics_gpu_metal.m" -DFPS_GPU_METAL -o "$bin_dir/fps_ray_gpu" $common_flags -framework Metal
# shellcheck disable=SC2086
clang "$project_root/fps_ray.c" -o "$bin_dir/fps_ray_cpu" $common_flags
# shellcheck disable=SC2086
clang "$project_root/fps_launcher.c" -std=c11 $opt_flags $arch_flags \
    -mmacosx-version-min="$deployment_target" -o "$bin_dir/fps_ray"

mkdir -p "$bin_dir/shaders/pbd"
xcrun -sdk macosx metal -c "$project_root/shaders/pbd/pbd_pipeline.metal" \
    -o "$build_root/pbd_pipeline.air" -mmacosx-version-min="$deployment_target"
xcrun -sdk macosx metallib "$build_root/pbd_pipeline.air" \
    -o "$bin_dir/shaders/pbd/pbd_pipeline.metallib"
cp -R "$project_root/shaders/." "$bin_dir/shaders/"

cp "$bin_dir/fps_ray" "$bin_dir/fps_ray_gpu" "$bin_dir/fps_ray_cpu" "$app_dir/Contents/MacOS/"
cp "$project_root/macos/Info.plist" "$app_dir/Contents/Info.plist"
mkdir -p "$app_dir/Contents/Resources/shaders"
cp -R "$bin_dir/shaders/." "$app_dir/Contents/Resources/shaders/"
codesign --force --deep --sign - "$app_dir"

echo "Built $architecture macOS binaries in $bin_dir"
echo "Built app bundle at $app_dir"
