param(
    [string]$ToolchainBin = "D:\Program Files\msys2\ucrt64\bin"
)

$ErrorActionPreference = 'Stop'
$weaponProjectRoot = Split-Path -Parent $PSScriptRoot
$weaponBuildDir = Join-Path $weaponProjectRoot '.build\weapon'
$weaponArtifactDir = Join-Path $weaponProjectRoot 'artifacts\weapon'
$weaponCompiler = Join-Path $ToolchainBin 'gcc.exe'
$weaponOldPath = $env:PATH
$weaponOldLocation = Get-Location
New-Item -ItemType Directory -Force -Path $weaponBuildDir, $weaponArtifactDir | Out-Null

function Invoke-WeaponCompiler([string[]]$CompilerArguments) {
    & $weaponCompiler @CompilerArguments
    if ($LASTEXITCODE -ne 0) { throw 'Weapon test compilation failed' }
}

try {
    $env:PATH = "$ToolchainBin;$weaponOldPath"
    Set-Location -LiteralPath $weaponProjectRoot
    Invoke-WeaponCompiler @('-std=c11', '-O2', '-Wall', '-Wextra', '-I.',
        'tests/net_protocol_test.c', 'net_protocol.c', '-o', "$weaponBuildDir\net_protocol_test.exe")
    & "$weaponBuildDir\net_protocol_test.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Network serialization test failed' }

    foreach ($api in @('33', '43')) {
        $raylibDir = ".build/raylib-gl$api/src"
        if (-not (Test-Path -LiteralPath "$raylibDir/libraylib.a")) {
            throw 'Build raylib first with build.ps1, then run this check again.'
        }
        $preview = "$weaponBuildDir\render_weapon_gl$api.exe"
        Invoke-WeaponCompiler @('-std=c11', '-O2', '-Wall', '-Wextra', "-DGRAPHICS_API_OPENGL_$api",
            '-I.', "-I$raylibDir", 'tests/render_weapon_preview.c', 'weapon_renderer.c',
            "-L$raylibDir", '-lraylib', '-lopengl32', '-lgdi32', '-lwinmm', '-lm', '-o', $preview)
        $destination = Join-Path $weaponArtifactDir "gl$api"
        New-Item -ItemType Directory -Force -Path $destination | Out-Null
        & $preview $destination *> "$weaponBuildDir\gl$api.log"
        if ($LASTEXITCODE -ne 0) { throw "OpenGL $api preview failed; see .build/weapon/gl$api.log" }
    }
    $noBloomDir = Join-Path $weaponArtifactDir 'no-bloom'
    New-Item -ItemType Directory -Force -Path $noBloomDir | Out-Null
    & "$weaponBuildDir\render_weapon_gl33.exe" $noBloomDir --no-bloom *> "$weaponBuildDir\no-bloom.log"
    if ($LASTEXITCODE -ne 0) { throw 'No-bloom preview failed' }

    $networkSources = @('net_protocol.c', 'net_transport.c',
        'third_party/enet/callbacks.c', 'third_party/enet/compress.c',
        'third_party/enet/host.c', 'third_party/enet/list.c', 'third_party/enet/packet.c',
        'third_party/enet/peer.c', 'third_party/enet/protocol.c', 'third_party/enet/win32.c')
    Invoke-WeaponCompiler (@('tests/render_weapon_gameplay.c', 'weapon_renderer.c') + $networkSources +
        @('-std=c11', '-O2', '-I.', '-Ithird_party/enet/include', '-I.build/raylib-gl33/src',
          '-L.build/raylib-gl33/src', '-lraylib', '-lopengl32', '-lgdi32', '-lwinmm',
          '-lws2_32', '-lpthread', '-lm', '-o', "$weaponBuildDir\render_weapon_gameplay.exe"))
    $gameplayDir = Join-Path $weaponArtifactDir 'gameplay'
    New-Item -ItemType Directory -Force -Path $gameplayDir | Out-Null
    & "$weaponBuildDir\render_weapon_gameplay.exe" $gameplayDir *> "$weaponBuildDir\gameplay.log"
    if ($LASTEXITCODE -ne 0) { throw 'Gameplay integration failed; see .build/weapon/gameplay.log' }
    Select-String -Path "$weaponBuildDir\gl33.log", "$weaponBuildDir\gl43.log",
        "$weaponBuildDir\no-bloom.log", "$weaponBuildDir\gameplay.log" -Pattern 'Four-view|passed'
    Write-Host "Weapon checks passed. Screenshots: $weaponArtifactDir"
} finally {
    $env:PATH = $weaponOldPath
    Set-Location -LiteralPath $weaponOldLocation.Path
}
