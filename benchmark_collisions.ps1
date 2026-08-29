param(
    [string]$BinDir = (Join-Path $PSScriptRoot ".build\bin"),
    [int[]]$VoxelCounts = @(256, 1024, 4096),
    [string[]]$Shapes = @("dense", "sparse", "colliding", "static"),
    [int]$Warmup = 60,
    [int]$Steps = 300,
    [int]$Repeats = 5,
    [string[]]$Backends = @("cpu-st", "cpu-mt", "gpu"),
    [string[]]$Profiles = @("wide-grid", "baseline", "hash", "cell-linked", "cell-span"),
    [int[]]$HashFactors = @(1, 2, 4, 8, 16),
    [switch]$DetailedCounters,
    [string]$OutputDir = (Join-Path $PSScriptRoot ("artifacts\collision-bench-" + (Get-Date -Format "yyyyMMdd-HHmmss")))
)

$ErrorActionPreference = "Continue"
if (Test-Path variable:PSNativeCommandUseErrorActionPreference) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$cpuExe = Join-Path $BinDir "fps_ray_cpu.exe"
$gpuExe = Join-Path $BinDir "fps_ray_gpu.exe"
if (-not (Test-Path -LiteralPath $cpuExe) -or -not (Test-Path -LiteralPath $gpuExe)) {
    throw "Build both executables first with .\build.ps1"
}

$shaderSource = Join-Path $PSScriptRoot "shaders\pbd\pbd_pipeline.comp"
$shaderTarget = Join-Path $BinDir "shaders\pbd\pbd_pipeline.comp"
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $shaderTarget) | Out-Null
Copy-Item -LiteralPath $shaderSource -Destination $shaderTarget -Force
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

$configs = foreach ($profile in $Profiles) {
    if ($profile -eq "hash") {
        foreach ($factor in $HashFactors) { [pscustomobject]@{ Profile = $profile; Factor = $factor } }
    } else {
        [pscustomobject]@{ Profile = $profile; Factor = if ($profile -in @("wide-grid", "baseline")) { 16 } else { 4 } }
    }
}

$rows = [System.Collections.Generic.List[object]]::new()
$oldLocation = Get-Location
try {
    Set-Location -LiteralPath $BinDir
    foreach ($backend in $Backends) {
        $exe = if ($backend -eq "gpu") { $gpuExe } else { $cpuExe }
        foreach ($shape in $Shapes) {
            foreach ($voxels in $VoxelCounts) {
                foreach ($config in $configs) {
                    for ($repeat = 1; $repeat -le $Repeats; ++$repeat) {
                        $arguments = @(
                            "--physics=$backend",
                            "--physics-smoke=$Steps",
                            "--physics-smoke-warmup=$Warmup",
                            "--physics-smoke-voxels=$voxels",
                            "--physics-smoke-shape=$shape",
                            "--collision-profile=$($config.Profile)",
                            "--collision-hash-factor=$($config.Factor)"
                        )
                        if ($DetailedCounters) { $arguments += "--collision-stats" }
                        $output = (& $exe @arguments 2>&1 | ForEach-Object { $_.ToString() }) -join "`n"
                        $exitCode = $LASTEXITCODE
                        $smoke = [regex]::Match($output, 'physics-smoke .*?simParticles=(?<particles>\d+).*?y=(?<y>-?[0-9.]+) msPerStep=(?<ms>[0-9.]+) fallback=(?<fallback>\d+)')
                        $stats = [regex]::Match($output, 'collision-stats occupied=(?<occupied>\d+) buildMs=(?<build>[0-9.]+) pairMs=(?<pair>[0-9.]+) staticMs=(?<static>[0-9.]+).*?aliases=(?<aliases>\d+).*?distance=(?<distance>\d+) contacts=(?<contacts>\d+)')
                        if ($exitCode -ne 0 -or -not $smoke.Success) {
                            throw "Benchmark failed: backend=$backend shape=$shape voxels=$voxels profile=$($config.Profile) factor=$($config.Factor)`n$output"
                        }
                        $rows.Add([pscustomobject]@{
                            Backend = $backend
                            Shape = $shape
                            Voxels = $voxels
                            Profile = $config.Profile
                            HashFactor = $config.Factor
                            Repeat = $repeat
                            Particles = [int]$smoke.Groups['particles'].Value
                            FinalY = [double]$smoke.Groups['y'].Value
                            MsPerStep = [double]$smoke.Groups['ms'].Value
                            OccupiedCells = if ($stats.Success) { [long]$stats.Groups['occupied'].Value } else { 0 }
                            BuildMsTotal = if ($stats.Success) { [double]$stats.Groups['build'].Value } else { 0 }
                            PairMsTotal = if ($stats.Success) { [double]$stats.Groups['pair'].Value } else { 0 }
                            StaticMsTotal = if ($stats.Success) { [double]$stats.Groups['static'].Value } else { 0 }
                            HashAliases = if ($stats.Success) { [long]$stats.Groups['aliases'].Value } else { 0 }
                            DistanceTests = if ($stats.Success) { [long]$stats.Groups['distance'].Value } else { 0 }
                            Contacts = if ($stats.Success) { [long]$stats.Groups['contacts'].Value } else { 0 }
                        })
                        Write-Host "$backend $shape $voxels $($config.Profile)/$($config.Factor) run $repeat`: $($smoke.Groups['ms'].Value) ms"
                    }
                }
            }
        }
    }
} finally {
    Set-Location $oldLocation
}

function Get-Median([double[]]$Values) {
    $sorted = $Values | Sort-Object
    $count = $sorted.Count
    if ($count -eq 0) { return 0.0 }
    $middle = [int][Math]::Floor($count / 2.0)
    if (($count % 2) -eq 1) { return [double]$sorted[$middle] }
    return 0.5 * ([double]$sorted[$middle - 1] + [double]$sorted[$middle])
}

$summary = $rows | Group-Object Backend, Shape, Voxels, Profile, HashFactor | ForEach-Object {
    $first = $_.Group[0]
    [pscustomobject]@{
        Backend = $first.Backend
        Shape = $first.Shape
        Voxels = $first.Voxels
        Profile = $first.Profile
        HashFactor = $first.HashFactor
        MedianMs = Get-Median ([double[]]$_.Group.MsPerStep)
        MinMs = ($_.Group.MsPerStep | Measure-Object -Minimum).Minimum
        MaxMs = ($_.Group.MsPerStep | Measure-Object -Maximum).Maximum
        Particles = $first.Particles
        FinalYSpread = (($($_.Group.FinalY | Measure-Object -Maximum).Maximum) - ($(($_.Group.FinalY | Measure-Object -Minimum).Minimum)))
    }
}

$comparisons = [System.Collections.Generic.List[object]]::new()
foreach ($case in ($summary | Group-Object Backend, Shape, Voxels)) {
    $items = $case.Group
    $baseline = $items | Where-Object Profile -eq "baseline" | Select-Object -First 1
    $wide = $items | Where-Object Profile -eq "wide-grid" | Select-Object -First 1
    $bestHash = $items | Where-Object Profile -eq "hash" | Sort-Object MedianMs | Select-Object -First 1
    $linked = $items | Where-Object Profile -eq "cell-linked" | Select-Object -First 1
    $span = $items | Where-Object Profile -eq "cell-span" | Select-Object -First 1
    $pairs = @(
        @("diameter-grid", $wide, $baseline),
        @("reduced-hash-collisions", $baseline, $bestHash),
        @("occupied-cell-iteration", $bestHash, $linked),
        @("contiguous-cell-spans", $linked, $span),
        @("cumulative", $baseline, $span)
    )
    foreach ($pair in $pairs) {
        $from = $pair[1]; $to = $pair[2]
        if ($null -eq $from -or $null -eq $to -or [double]$from.MedianMs -le 0.0) { continue }
        $comparisons.Add([pscustomobject]@{
            Backend = $items[0].Backend
            Shape = $items[0].Shape
            Voxels = $items[0].Voxels
            Upgrade = $pair[0]
            FromProfile = $from.Profile
            FromHashFactor = $from.HashFactor
            ToProfile = $to.Profile
            ToHashFactor = $to.HashFactor
            FromMedianMs = [double]$from.MedianMs
            ToMedianMs = [double]$to.MedianMs
            SpeedupPercent = 100.0 * ([double]$from.MedianMs - [double]$to.MedianMs) / [double]$from.MedianMs
        })
    }
}

$rows | Export-Csv -LiteralPath (Join-Path $OutputDir "runs.csv") -NoTypeInformation
$summary | Export-Csv -LiteralPath (Join-Path $OutputDir "summary.csv") -NoTypeInformation
$comparisons | Export-Csv -LiteralPath (Join-Path $OutputDir "comparisons.csv") -NoTypeInformation
[pscustomobject]@{
    GeneratedAt = (Get-Date).ToString("o")
    Warmup = $Warmup
    Steps = $Steps
    Repeats = $Repeats
    Runs = $rows
    Summary = $summary
    Comparisons = $comparisons
} | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $OutputDir "report.json") -Encoding utf8

Write-Host "Collision benchmark report: $OutputDir"
