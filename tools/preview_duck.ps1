param(
    [ValidateSet('interactive','capture','checks','fallback')][string]$Mode='interactive',
    [ValidateSet('33','43')][string]$Backend='43',
    [string]$Output='artifacts/duck'
)
$ErrorActionPreference='Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    New-Item -ItemType Directory -Force -Path '.build/duck', $Output | Out-Null
    $raylib=".build/raylib-gl$Backend/src"
    $exe=".build/duck/preview-gl$Backend.exe"
    & gcc tests/render_duck_preview.c -std=c11 -O2 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$raylib" "-L$raylib" -o $exe -lraylib -lopengl32 -lgdi32 -lwinmm -lm
    if($LASTEXITCODE -ne 0){throw 'Mallard preview compilation failed'}
    if($Mode -eq 'interactive'){ & $exe } else { & $exe "--$Mode" $Output }
    if($LASTEXITCODE -ne 0){throw 'Mallard preview validation failed'}
} finally { Pop-Location }
