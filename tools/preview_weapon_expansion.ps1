param(
    [ValidateSet('interactive','capture','checks','benchmark','fallback')][string]$Mode='interactive',
    [ValidateSet('33','43')][string]$Backend='43',
    [string]$Output='artifacts/weapon-expansion/rendered'
)
$ErrorActionPreference='Stop'
$weaponRoot=Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $weaponRoot
try {
    New-Item -ItemType Directory -Force -Path '.build/weapon', $Output | Out-Null
    $raylib=".build/raylib-gl$Backend/src"
    $exe=".build/weapon/expansion-gl$Backend.exe"
    # The harness includes the production renderer for private geometry assertions.
    & gcc tests/render_weapon_expansion.c character_renderer.c -std=c11 -O2 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$raylib" "-L$raylib" -o $exe -lraylib -lopengl32 -lgdi32 -lwinmm -lm
    if($LASTEXITCODE -ne 0){throw 'Preview compilation failed'}
    if($Mode -eq 'interactive'){ & $exe } else { & $exe "--$Mode" $Output }
    if($LASTEXITCODE -ne 0){throw 'Preview validation failed'}
} finally { Pop-Location }
