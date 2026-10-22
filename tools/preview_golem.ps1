param(
    [ValidateSet('interactive','capture','checks','benchmark','fallback')][string]$Mode='interactive',
    [ValidateSet('33','43')][string]$Backend='43',
    [string]$Output='artifacts/golem/rendered'
)
$ErrorActionPreference='Stop'
$golemRoot=Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $golemRoot
try {
    New-Item -ItemType Directory -Force -Path '.build/golem', $Output | Out-Null
    $raylib=".build/raylib-gl$Backend/src"
    $exe=".build/golem/preview-gl$Backend.exe"
    & gcc tests/render_golem_preview.c character_renderer.c weapon_renderer.c -std=c11 -O2 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$raylib" "-L$raylib" -o $exe -lraylib -lopengl32 -lgdi32 -lwinmm -lm
    if($LASTEXITCODE -ne 0){throw 'Golem preview compilation failed'}
    if($Mode -eq 'interactive'){ & $exe } else { & $exe "--$Mode" $Output }
    if($LASTEXITCODE -ne 0){throw 'Golem preview validation failed'}
} finally { Pop-Location }
