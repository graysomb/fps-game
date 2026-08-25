#!/usr/bin/env sh
set -eu

project_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
raylib_source=${RAYLIB_SOURCE_DIR:-}
configuration=${1:-release}

if [ -z "$raylib_source" ]; then
    echo "Set RAYLIB_SOURCE_DIR to a raylib 5.5 source tree" >&2
    exit 2
fi

build_root="$project_root/.build"
bin_dir="$build_root/bin"
raylib33="$build_root/raylib-gl33"
raylib43="$build_root/raylib-gl43"
mkdir -p "$bin_dir"

if [ ! -d "$raylib33/src" ]; then
    mkdir -p "$raylib33"
    cp -a "$raylib_source/." "$raylib33/"
fi
if [ ! -d "$raylib43/src" ]; then
    mkdir -p "$raylib43"
    cp -a "$raylib_source/." "$raylib43/"
fi

make -C "$raylib33/src" clean PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC
make -C "$raylib33/src" PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC
make -C "$raylib43/src" clean PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_43 RAYLIB_LIBTYPE=STATIC
make -C "$raylib43/src" PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_43 RAYLIB_LIBTYPE=STATIC

if [ "$configuration" = debug ]; then
    opt_flags="-O0 -g"
else
    opt_flags="-O2 -DNDEBUG"
fi

common_flags="-std=c11 $opt_flags -I$project_root -lGL -lm -lpthread -ldl -lrt -lX11"
# shellcheck disable=SC2086
cc "$project_root/fps_ray.c" -I"$raylib33/src" -L"$raylib33/src" -o "$bin_dir/fps_ray_cpu" -lraylib $common_flags
# shellcheck disable=SC2086
cc "$project_root/fps_ray.c" -DGRAPHICS_API_OPENGL_43 -I"$raylib43/src" -L"$raylib43/src" -o "$bin_dir/fps_ray_gpu" -lraylib $common_flags
cc "$project_root/fps_launcher.c" -std=c11 -O2 -o "$bin_dir/fps_ray"
cp -a "$project_root/shaders" "$bin_dir/"
echo "Built launcher and both physics backends in $bin_dir"
