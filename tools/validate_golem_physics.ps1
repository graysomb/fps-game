param([ValidateSet('33','43')][string]$Backend='43', [switch]$Capture, [switch]$SkipBuild, [switch]$Accuracy,
      [ValidateRange(0,3600)][int]$StressSeconds=0,
      [ValidateRange(0.9,1.0)][double]$Damping=0.99)
$ErrorActionPreference='Stop'
$oldCapture=$env:GOLEM_CAPTURE
$oldStress=$env:GOLEM_STRESS_SECONDS
$oldAccuracy=$env:GOLEM_ACCURACY
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    if($Accuracy){$env:GOLEM_ACCURACY='1'}
    New-Item -ItemType Directory -Force -Path '.build/golem','artifacts/golem/physics' | Out-Null
    if($Capture){
        $env:GOLEM_CAPTURE='1'
        @('idle','motion-1','motion-2','motion-3','motion-4','leg-loss','regeneration') | ForEach-Object {
            New-Item -ItemType Directory -Force -Path "artifacts/golem/physics/$_" | Out-Null
        }
    }
    if($StressSeconds){$env:GOLEM_STRESS_SECONDS="$StressSeconds"}
    $raylib=".build/raylib-gl$Backend/src"
    $dampingText=$Damping.ToString('0.######',[System.Globalization.CultureInfo]::InvariantCulture)
    $executable=".build/golem/test-$Backend-damping-$dampingText.exe"
    $sources=@('duck_renderer.c','weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object { "third_party/enet/$_.c" }
    if(!$SkipBuild){
    & gcc tests/golem_physics_test.c @sources -std=c11 -O2 "-DGRAPHICS_API_OPENGL_$Backend" "-DVELOCITY_DAMPING=$dampingText" -I. "-I$raylib" -Ithird_party/enet/include "-L$raylib" -o $executable -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lpthread -lm
    if($LASTEXITCODE){throw 'Golem physics build failed'}
    }
    $physics=if($Backend -eq '43'){'gpu'}else{'cpu-st'}
    & $executable "--physics=$physics"
    if($LASTEXITCODE){throw 'Golem physics checks failed'}
} finally { $env:GOLEM_CAPTURE=$oldCapture;$env:GOLEM_STRESS_SECONDS=$oldStress;$env:GOLEM_ACCURACY=$oldAccuracy;Pop-Location }
