param([switch]$Gpu)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    $raylib = if ($Gpu) { '.build/raylib-gl43/src' } else { '.build/raylib-gl33/src' }
    $variant = if ($Gpu) { 'gpu' } else { 'cpu' }
    $sources = @('weapon_renderer.c','net_protocol.c','net_transport.c')
    $sources += @('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object { "third_party/enet/$_.c" }
    $flags = @('-std=c11','-O1','-g',"-I$raylib",'-Ithird_party/enet/include',"-L$raylib",'-lraylib','-lopengl32','-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
    if ($Gpu) { $flags += '-DGRAPHICS_API_OPENGL_43' }
    & gcc tests/firefight_rules_test.c -std=c11 -o .build/firefight_rules_test.exe
    if ($LASTEXITCODE -ne 0) { throw 'Rules build failed' }
    & .build/firefight_rules_test.exe
    if ($LASTEXITCODE -ne 0) { throw 'Rules test failed' }
    & gcc tests/firefight_integration.c @sources @flags -o ".build/firefight_integration_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Integration build failed' }
    $backend = if ($Gpu) { 'gpu' } else { 'cpu-mt' }
    & ".build/firefight_integration_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Integration test failed' }
    $gameDir = '.build/couch'
    New-Item -ItemType Directory -Force -Path $gameDir | Out-Null
    Copy-Item -LiteralPath shaders -Destination $gameDir -Recurse -Force
    Copy-Item -LiteralPath '.build/bin/libwinpthread-1.dll' -Destination $gameDir -Force
    Copy-Item -LiteralPath 'game_song.mp3' -Destination $gameDir -Force
    & gcc fps_ray.c @sources @flags -o "$gameDir/fps_ray_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Game build failed' }
} finally { Pop-Location }
