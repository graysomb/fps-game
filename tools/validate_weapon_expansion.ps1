param([ValidateSet('33','43')][string]$Backend='43',[switch]$Stress)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $root
try {
    New-Item -ItemType Directory -Force .build/weapon,artifacts | Out-Null
    $lib=".build/raylib-gl$Backend/src"
    $sources=@('weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"}
    $flags=@('-std=c11','-O1',"-DGRAPHICS_API_OPENGL_$Backend",'-I.',"-I$lib",'-Ithird_party/enet/include',"-L$lib",'-lraylib','-lopengl32','-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
    $game=".build/weapon/gameplay-$Backend.exe"
    & gcc tests/weapon_expansion_gameplay.c @sources @flags -o $game
    if($LASTEXITCODE){throw 'Gameplay build failed'}
    $physics=if($Backend -eq '43'){'gpu'}else{'cpu-mt'}
    if($Stress){$env:FPS_EXPANSION_STRESS='1'}
    try { & $game "--physics=$physics"; if($LASTEXITCODE){throw 'Gameplay checks failed'} }
    finally {if($Stress){Remove-Item Env:FPS_EXPANSION_STRESS}}
    $lan=".build/weapon/lan-$Backend.exe"
    & gcc tests/weapon_expansion_lan.c @sources @flags -o $lan
    if($LASTEXITCODE){throw 'LAN build failed'}
    $hostProcess=Start-Process -FilePath "$root/$lan" -ArgumentList host -WindowStyle Hidden -PassThru -RedirectStandardOutput "$root/.build/weapon/host-$Backend.log" -RedirectStandardError "$root/.build/weapon/host-$Backend.err"
    & $lan client
    $clientCode=$LASTEXITCODE
    $hostProcess.WaitForExit()
    if($clientCode -or $hostProcess.ExitCode){throw 'LAN checks failed'}
} finally {Pop-Location}
