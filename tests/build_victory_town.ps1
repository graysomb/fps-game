param([ValidateSet('33','43')][string]$GL='43',[string]$ToolchainBin='D:\Program Files\msys2\ucrt64\bin')
$ErrorActionPreference='Stop'
Push-Location $PSScriptRoot/..
$oldPath=$env:PATH
try {
 $env:PATH="$ToolchainBin;$oldPath"
 $dest=".build/victory-town-gl$GL"
 New-Item -ItemType Directory -Force $dest | Out-Null
 $src=@('weapon_renderer.c','character_renderer.c','duck_renderer.c','net_protocol.c','net_transport.c')+(@('callbacks','compress','host','list','packet','peer','protocol','win32') | ForEach-Object {"third_party/enet/$_.c"})
 $defs=@();if($GL -eq '43'){$defs=@('-DGRAPHICS_API_OPENGL_43')}
 & gcc tests/render_victory_town.c @src -std=c11 -O2 @defs "-I.build/raylib-gl$GL/src" -Ithird_party/enet/include "-L.build/raylib-gl$GL/src" -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lpthread -lm -o "$dest/victory-town.exe"
 if($LASTEXITCODE -ne 0){throw 'Preview build failed'}
 Copy-Item shaders $dest -Recurse -Force
 foreach($dll in @('libwinpthread-1.dll','libgcc_s_seh-1.dll')){if(Test-Path "$ToolchainBin/$dll"){Copy-Item "$ToolchainBin/$dll" $dest -Force}}
 New-Item -ItemType Directory -Force "$dest/artifacts/victory-town/frames","$dest/artifacts/victory-town/scenes","$dest/artifacts/victory-town/qa" | Out-Null
 Write-Host "Standalone preview: $dest/victory-town.exe"
} finally {$env:PATH=$oldPath;Pop-Location}
