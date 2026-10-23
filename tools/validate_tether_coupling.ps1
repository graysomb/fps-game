param([ValidateSet('33','43')][string]$Backend='43',[switch]$SkipBuild)
$ErrorActionPreference='Stop'
Push-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
try {
    New-Item -ItemType Directory -Force .build | Out-Null
    $lib=".build/raylib-gl$Backend/src"
    $sources=@('weapon_renderer.c','character_renderer.c','duck_renderer.c','net_protocol.c','net_transport.c')
    $sources+=@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"}
    $flags=@('-std=c11','-O1',"-DGRAPHICS_API_OPENGL_$Backend",'-I.',"-I$lib",'-Ithird_party/enet/include',"-L$lib",'-lraylib','-lopengl32','-lgdi32','-lwinmm','-lws2_32','-lpthread','-lm')
    $physics=if($Backend -eq '43'){'gpu'}else{'cpu-st'}
    foreach($test in @('tether_grapple_test','tether_coupling_physics','tether_grapple_lan')){
        $exe=".build/$test-coupled-$Backend.exe"
        if(!$SkipBuild){& gcc "tests/$test.c" @sources @flags -o $exe;if($LASTEXITCODE){throw "$test build failed"}}
        if($test -ne 'tether_grapple_lan'){
            & $exe "--physics=$physics";if($LASTEXITCODE){throw "$test failed"}
            if($Backend -eq '43'){
                & $exe --physics=cpu-st;if($LASTEXITCODE){throw "$test CPU comparison failed"}
                $old=$env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE
                try {$env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE='1';& $exe --physics=gpu;if($LASTEXITCODE){throw "$test GPU fallback failed"}}
                finally {$env:FPS_FORCE_GPU_VALIDATION_FAILURE_ONCE=$old}
            }
        }else{
            $root=(Get-Location).Path
            $server=Start-Process -FilePath (Join-Path $root $exe) -ArgumentList host -WorkingDirectory $root -WindowStyle Hidden -PassThru -RedirectStandardOutput ".build/tether-host-$Backend.log" -RedirectStandardError ".build/tether-host-$Backend.err"
            $null=$server.Handle
            & $exe client;$clientCode=$LASTEXITCODE;$server.WaitForExit()
            if($clientCode -ne 0 -or $server.ExitCode -ne 0){Get-Content ".build/tether-host-$Backend.err";throw 'Grapple LAN failed'}
        }
    }
} finally {Pop-Location}
