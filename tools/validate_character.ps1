param([string]$ToolchainBin = 'D:\Program Files\msys2\ucrt64\bin')
$ErrorActionPreference = 'Stop'
$oldPath = $env:PATH
Push-Location (Split-Path -Parent $PSScriptRoot)
try {
    $env:PATH = "$ToolchainBin;$oldPath"
    New-Item -ItemType Directory -Force '.build/rat','artifacts/rat' | Out-Null
    $network = @('duck_renderer.c','net_protocol.c','net_transport.c')
    $network += @('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object { "third_party/enet/$_.c" }
    foreach ($api in @('33','43')) {
        $raylib = ".build/raylib-gl$api/src"
        $flags = @('-std=c11','-O1','-g',"-DGRAPHICS_API_OPENGL_$api","-I$raylib",'-Ithird_party/enet/include',"-L$raylib",'-lraylib','-lopengl32','-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
        & gcc tests/character_mesh_test.c @flags -o ".build/rat/mesh-$api.exe"
        if ($LASTEXITCODE) { throw 'Mesh build failed' }
        & ".build/rat/mesh-$api.exe" *> ".build/rat/mesh-$api.log"
        if ($LASTEXITCODE) { throw 'Mesh validation failed' }
        & ".build/rat/mesh-$api.exe" fallback *> ".build/rat/fallback-$api.log"
        if ($LASTEXITCODE) { throw 'Shader fallback validation failed' }
        & gcc tests/render_character_preview.c character_renderer.c weapon_renderer.c @flags -o ".build/rat/preview-$api.exe"
        if ($LASTEXITCODE) { throw 'Preview build failed' }
        & ".build/rat/preview-$api.exe" artifacts/rat *> ".build/rat/preview-$api.log"
        if ($LASTEXITCODE) { throw 'Preview or emission occlusion check failed' }
        & gcc tests/character_integration.c character_renderer.c weapon_renderer.c @network @flags -o ".build/rat/integration-$api.exe"
        if ($LASTEXITCODE) { throw 'Integration build failed' }
        & ".build/rat/integration-$api.exe" *> ".build/rat/integration-$api.log"
        if ($LASTEXITCODE) { throw 'Character integration failed' }
        Get-Content ".build/rat/integration-$api.log" | Select-String 'Rat '
        & gcc tests/character_grounding.c character_renderer.c weapon_renderer.c @network @flags -o ".build/rat/grounding-$api.exe"
        if ($LASTEXITCODE) { throw 'Grounding build failed' }
        $backend = if ($api -eq '43') { 'gpu' } else { 'cpu-mt' }
        & ".build/rat/grounding-$api.exe" "--physics=$backend" *> ".build/rat/grounding-$api.log"
        if ($LASTEXITCODE) { throw 'Grounding validation failed' }
    }
    & gcc tests/net_protocol_test.c net_protocol.c -I. -std=c11 -o .build/rat/net-test.exe
    if ($LASTEXITCODE) { throw 'Protocol build failed' }
    & .build/rat/net-test.exe
    if ($LASTEXITCODE) { throw 'Protocol regression failed' }
    Write-Host 'Rat renderer checks passed. Previews: artifacts/rat'
} finally {
    $env:PATH = $oldPath
    Pop-Location
}
