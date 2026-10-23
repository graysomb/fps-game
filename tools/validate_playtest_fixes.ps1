param([ValidateSet('33','43')][string]$Backend='43')
$ErrorActionPreference='Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    New-Item -ItemType Directory -Force .build,artifacts | Out-Null
    $lib=".build/raylib-gl$Backend/src"
    $sources=@('weapon_renderer.c','character_renderer.c','duck_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"}
    & gcc tests/playtest_fixes.c @sources -std=c11 -O1 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$lib" -Ithird_party/enet/include "-L$lib" -o ".build/playtest-fixes-$Backend.exe" -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lpthread -lm
    if($LASTEXITCODE){throw 'Playtest build failed'}
    & ".build/playtest-fixes-$Backend.exe"
    if($LASTEXITCODE){throw 'Playtest checks failed'}
} finally {Pop-Location}
