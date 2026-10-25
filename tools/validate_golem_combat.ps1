param([ValidateSet('33','43')][string]$Backend='43',[switch]$SkipBuild,[switch]$Capture,[switch]$Preview,[switch]$Fallback,[ValidateRange(0,3600)][int]$StressSeconds=0)
$ErrorActionPreference='Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$oldCapture=$env:GOLEM_COMBAT_CAPTURE
$oldPreview=$env:GOLEM_COMBAT_PREVIEW
$oldStress=$env:GOLEM_COMBAT_STRESS
$oldFallback=$env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE
try {
    New-Item -ItemType Directory -Force '.build/golem','artifacts/golem/combat' | Out-Null
    $lib=".build/raylib-gl$Backend/src"
    $sources=@('duck_renderer.c','weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object { "third_party/enet/$_.c" }
    $exe=".build/golem/combat-$Backend.exe"
    if(!$SkipBuild){
        & gcc tests/golem_combat_test.c @sources -std=c11 -O1 "-DGRAPHICS_API_OPENGL_$Backend" -I. "-I$lib" -Ithird_party/enet/include "-L$lib" -o $exe -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lpthread -lm
        if($LASTEXITCODE){throw 'Golem combat test build failed'}
    }
    if($Preview){$env:GOLEM_COMBAT_PREVIEW="1"}
    if($Capture){$env:GOLEM_COMBAT_CAPTURE='1'}
    if($StressSeconds){$env:GOLEM_COMBAT_STRESS="$StressSeconds"}
    if($Fallback){$env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE='1'}
    $physics=if($Backend -eq '43'){'gpu'}else{'cpu-st'}
    & $exe "--physics=$physics"
    if($LASTEXITCODE){throw 'Golem combat checks failed'}
} finally {
    $env:GOLEM_COMBAT_PREVIEW=$oldPreview;$env:GOLEM_COMBAT_CAPTURE=$oldCapture;$env:GOLEM_COMBAT_STRESS=$oldStress
    $env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE=$oldFallback
    Pop-Location
}
