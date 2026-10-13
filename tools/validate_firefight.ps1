param([switch]$Gpu)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    $raylib = if ($Gpu) { '.build/raylib-gl43/src' } else { '.build/raylib-gl33/src' }
    $variant = if ($Gpu) { 'gpu' } else { 'cpu' }
    $sources = @('weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
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
    & gcc tests/launcher_recycling_test.c @sources @flags -o ".build/launcher_recycling_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Launcher recycling build failed' }
    & ".build/launcher_recycling_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Launcher recycling regression failed' }
    & gcc tests/static_cache_test.c @sources @flags -o ".build/static_cache_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Static cache build failed' }
    & ".build/static_cache_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Static cache regression failed' }
    & gcc tests/activation_batch_test.c @sources @flags -o ".build/activation_batch_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Activation batching build failed' }
    & ".build/activation_batch_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Activation batching regression failed' }
    & gcc tests/chunk_mesh_test.c @sources @flags -o ".build/chunk_mesh_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Chunk mesh build failed' }
    & ".build/chunk_mesh_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Chunk mesh regression failed' }
    & gcc tests/chunk_regeneration_test.c @sources @flags -o ".build/chunk_regeneration_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Chunk regeneration build failed' }
    & ".build/chunk_regeneration_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Chunk regeneration regression failed' }
    & ".build/firefight_integration_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Integration test failed' }
    & gcc tests/firefight_playtest.c @sources @flags -o ".build/firefight_playtest_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Playtest regression build failed' }
    & ".build/firefight_playtest_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Playtest regression failed' }
    & gcc tests/combat_clarity_test.c @sources @flags -o ".build/combat_clarity_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Combat clarity build failed' }
    & ".build/combat_clarity_$variant.exe" "--physics=$backend"
    if ($LASTEXITCODE -ne 0) { throw 'Combat clarity regression failed' }
    & gcc tests/combat_clarity_lan.c @sources @flags -o ".build/combat_clarity_lan_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'LAN regression build failed' }
    $lanHost = Start-Process -FilePath (Join-Path $root ".build/combat_clarity_lan_$variant.exe") -ArgumentList 'host' -WindowStyle Hidden -PassThru -RedirectStandardOutput ".build/clarity-lan-host-$variant.log" -RedirectStandardError ".build/clarity-lan-host-$variant.err"
    $null = $lanHost.Handle # Retain the process handle so Windows PowerShell can read ExitCode after exit.
    & ".build/combat_clarity_lan_$variant.exe" client
    $clientExit = $LASTEXITCODE
    $lanHost.WaitForExit()
    if ($clientExit -ne 0 -or $lanHost.ExitCode -ne 0) { throw 'LAN regression failed' }
    $gameDir = '.build/couch'
    New-Item -ItemType Directory -Force -Path $gameDir | Out-Null
    Copy-Item -LiteralPath shaders -Destination $gameDir -Recurse -Force
    Copy-Item -LiteralPath '.build/bin/libwinpthread-1.dll' -Destination $gameDir -Force
    Copy-Item -LiteralPath 'game_song.mp3' -Destination $gameDir -Force
    & gcc fps_ray.c @sources @flags -o "$gameDir/fps_ray_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Game build failed' }
    & "$gameDir/fps_ray_$variant.exe" "--physics=$backend" --enemy-smoke=all
    if ($LASTEXITCODE -ne 0) { throw 'Enemy behavior tests failed' }
} finally { Pop-Location }
