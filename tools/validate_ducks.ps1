param([ValidateSet('33','43')][string]$Backend='43')
$ErrorActionPreference='Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    New-Item -ItemType Directory -Force -Path '.build/duck','artifacts/duck' | Out-Null
    $raylib=".build/raylib-gl$Backend/src"
    $sources=@('duck_renderer.c','weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object { "third_party/enet/$_.c" }
    & gcc tests/duck_gameplay_test.c @sources -std=c11 -O0 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$raylib" -Ithird_party/enet/include "-L$raylib" -o ".build/duck/gameplay-$Backend.exe" -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lpthread -lm
    if($LASTEXITCODE){throw 'Duck gameplay build failed'}
    $physics=if($Backend -eq '43'){'gpu'}else{'cpu-mt'}
    & ".build/duck/gameplay-$Backend.exe" "--physics=$physics"
    if($LASTEXITCODE){throw 'Duck gameplay checks failed'}
} finally { Pop-Location }
