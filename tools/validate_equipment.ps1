param([switch]$Gpu)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
 $api=if($Gpu){'43'}else{'33'}
 $backend=if($Gpu){'gpu-gl43'}else{'cpu-mt'}
 $lib=".build/raylib-gl$api/src"
 $sources=@('weapon_renderer.c','character_renderer.c','net_protocol.c','net_transport.c')
 $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"}
 $flags=@('-std=c11','-O1',"-DGRAPHICS_API_OPENGL_$api",'-I.',"-I$lib",'-Ithird_party/enet/include',"-L$lib",'-lraylib','-lopengl32','-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
 New-Item -ItemType Directory -Force .build,artifacts | Out-Null
 & gcc tests/equipment_test.c @sources @flags -o ".build/equipment-$api.exe"
 if($LASTEXITCODE){throw 'Equipment build failed'}
 & ".build/equipment-$api.exe" "--physics=$backend"
 if($LASTEXITCODE){throw 'Equipment checks failed'}
 & gcc tests/sword_wave_test.c @sources @flags -o ".build/sword-wave-$api.exe"
 if($LASTEXITCODE){throw "Sword wave build failed"}
 & ".build/sword-wave-$api.exe" "--physics=$backend"
 if($LASTEXITCODE){throw "Sword wave checks failed"}
 & gcc tests/aiming_test.c @sources @flags -o ".build/aiming-$api.exe"
 if($LASTEXITCODE){throw 'Aiming build failed'}
 & ".build/aiming-$api.exe" "--physics=$backend"
 if($LASTEXITCODE){throw 'Aiming checks failed'}
 & gcc tests/equipment_lan.c @sources @flags -o ".build/equipment-lan-$api.exe"
 if($LASTEXITCODE){throw 'LAN build failed'}
 $hostProcess=Start-Process -FilePath "$root/.build/equipment-lan-$api.exe" -ArgumentList host -WindowStyle Hidden -PassThru -RedirectStandardOutput "$root/.build/equipment-host-$api.log" -RedirectStandardError "$root/.build/equipment-host-$api.err"
 & ".build/equipment-lan-$api.exe" client
 $clientCode=$LASTEXITCODE
 $hostProcess.WaitForExit()
 if($clientCode -or $hostProcess.ExitCode){throw 'Equipment LAN checks failed'}
} finally {Pop-Location}
