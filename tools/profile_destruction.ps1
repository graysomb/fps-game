param(
    [ValidateSet('foundry','versus')][string]$Map='foundry',
    [ValidateSet('gpu','cpu-mt')][string]$Backend='gpu',
    [ValidateRange(0,24)][int]$Grenades=24,
    [ValidateRange(120,36000)][int]$Frames=3600,
    [switch]$RestorationOnly,
    [switch]$ActivateStructure
)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
$keys=@('FPS_STRESS_MAP','FPS_STRESS_GRENADES','FPS_STRESS_FRAMES','FPS_STRESS_RESTORE_ONLY','FPS_STRESS_ACTIVATE')
$previous=@{}
foreach($key in $keys){$previous[$key]=[Environment]::GetEnvironmentVariable($key,'Process')}
Push-Location $root
try {
    New-Item -ItemType Directory -Force .build,artifacts | Out-Null
    $sources=@('weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"}
    $flags=@('-std=c11','-O1','-g','-DGRAPHICS_API_OPENGL_43','-I.build/raylib-gl43/src',
        '-Ithird_party/enet/include','-L.build/raylib-gl43/src','-lraylib','-lopengl32',
        '-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
    & gcc tests/destruction_profile.c @sources @flags -o .build/destruction_profile.exe
    if($LASTEXITCODE -ne 0){throw 'Profiler build failed'}
    $env:FPS_STRESS_ACTIVATE=if($ActivateStructure){'1'}else{$null}
    $env:FPS_STRESS_MAP=$Map
    $env:FPS_STRESS_GRENADES="$Grenades"
    $env:FPS_STRESS_FRAMES="$Frames"
    $env:FPS_STRESS_RESTORE_ONLY=if($RestorationOnly){'1'}else{$null}
    if($RestorationOnly){$env:FPS_STRESS_GRENADES='0';$env:FPS_STRESS_FRAMES='120'}
    if($ActivateStructure){$env:FPS_STRESS_GRENADES='0'}
    $suffix=if($ActivateStructure){'activation'}elseif($RestorationOnly){'restoration'}else{'grenades'}
    $report="artifacts/profile-$Map-$Backend-$suffix"
    # The harness prints measurements to stderr; keep native output out of PS error handling.
    $process=Start-Process -FilePath (Join-Path $root '.build/destruction_profile.exe') `
        -ArgumentList "--physics=$Backend" -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput "$report-render.log" -RedirectStandardError "$report.log"
    $null=$process.Handle
    $process.WaitForExit()
    if($process.ExitCode -ne 0){throw "Profiler failed; see $report.log"}
    Get-Content "$report.log" | Where-Object {$_ -match '^(frame,|[0-9]+,)'} | Set-Content "$report.csv"
    Get-Content "$report.log" | Select-Object -Last 5
    Write-Output "Saved $report.csv and $report.log"
} finally {
    foreach($key in $keys){[Environment]::SetEnvironmentVariable($key,$previous[$key],'Process')}
    Pop-Location
}
