param(
    [string]$RaylibSource = $env:RAYLIB_SOURCE_DIR,
    [string]$ToolchainBin = "D:\Program Files\msys2\ucrt64\bin",
    [ValidateSet("debug", "release")]
    [string]$Configuration = "release"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not $RaylibSource) {
    $RaylibSource = "D:\raylib-master\raylib-master"
}
$RaylibSource = (Resolve-Path -LiteralPath $RaylibSource).Path
$BuildRoot = Join-Path $ProjectRoot ".build"
$BinDir = Join-Path $BuildRoot "bin"
$Raylib33 = Join-Path $BuildRoot "raylib-gl33"
$Raylib43 = Join-Path $BuildRoot "raylib-gl43"
$Make = Join-Path $ToolchainBin "mingw32-make.exe"
$Gcc = Join-Path $ToolchainBin "gcc.exe"

New-Item -ItemType Directory -Force -Path $BuildRoot, $BinDir | Out-Null

function Copy-RaylibTree([string]$Destination) {
    if (Test-Path -LiteralPath $Destination) { return }
    New-Item -ItemType Directory -Force -Path $Destination | Out-Null
    Get-ChildItem -LiteralPath $RaylibSource -Force | Copy-Item -Destination $Destination -Recurse -Force
}

Copy-RaylibTree $Raylib33
Copy-RaylibTree $Raylib43

$PreviousPath = $env:PATH
$env:PATH = "$ToolchainBin;$PreviousPath"
try {
    & $Make -C (Join-Path $Raylib33 "src") clean PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC
    if ($LASTEXITCODE -ne 0) { throw "OpenGL 3.3 raylib clean failed" }
    & $Make -C (Join-Path $Raylib33 "src") PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_33 RAYLIB_LIBTYPE=STATIC
    if ($LASTEXITCODE -ne 0) { throw "OpenGL 3.3 raylib build failed" }
    & $Make -C (Join-Path $Raylib43 "src") clean PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_43 RAYLIB_LIBTYPE=STATIC
    if ($LASTEXITCODE -ne 0) { throw "OpenGL 4.3 raylib clean failed" }
    & $Make -C (Join-Path $Raylib43 "src") PLATFORM=PLATFORM_DESKTOP GRAPHICS=GRAPHICS_API_OPENGL_43 RAYLIB_LIBTYPE=STATIC
    if ($LASTEXITCODE -ne 0) { throw "OpenGL 4.3 raylib build failed" }

    $Optimization = if ($Configuration -eq "debug") { "-O0", "-g" } else { "-O2", "-DNDEBUG" }
    $Common = @("-std=c11") + $Optimization + @("-I$ProjectRoot", "-lopengl32", "-lgdi32", "-lwinmm", "-lpthread", "-lm")
    & $Gcc (Join-Path $ProjectRoot "fps_ray.c") "-I$(Join-Path $Raylib33 'src')" "-L$(Join-Path $Raylib33 'src')" "-o$(Join-Path $BinDir 'fps_ray_cpu.exe')" "-lraylib" @Common
    if ($LASTEXITCODE -ne 0) { throw "CPU game build failed" }
    & $Gcc (Join-Path $ProjectRoot "fps_ray.c") "-DGRAPHICS_API_OPENGL_43" "-I$(Join-Path $Raylib43 'src')" "-L$(Join-Path $Raylib43 'src')" "-o$(Join-Path $BinDir 'fps_ray_gpu.exe')" "-lraylib" @Common
    if ($LASTEXITCODE -ne 0) { throw "GPU game build failed" }
    & $Gcc (Join-Path $ProjectRoot "fps_launcher.c") "-std=c11" "-O2" "-o$(Join-Path $BinDir 'fps_ray.exe')"
    if ($LASTEXITCODE -ne 0) { throw "Launcher build failed" }
} finally {
    $env:PATH = $PreviousPath
}

Copy-Item -LiteralPath (Join-Path $ProjectRoot "shaders") -Destination $BinDir -Recurse -Force
Write-Host "Built launcher and both physics backends in $BinDir"
