$ErrorActionPreference = 'Stop'
$Url = 'https://raw.githubusercontent.com/manaakiwhenua/SleeperPests/6e4b9f032bec9001bf7241a9d59625894c776f9b/INApest%20INApestMeta%20core%20function%20code/INApestAnalytical.R'
$Out = Join-Path $PSScriptRoot 'R\INApestAnalytical.R'
Invoke-WebRequest -Uri $Url -OutFile $Out
$Expected = '49238e62a9f4c99445072d08356ca33081c60b37f1a4b268b328794990d163c5'
$Actual = (Get-FileHash -Algorithm SHA256 $Out).Hash.ToLowerInvariant()
if ($Actual -ne $Expected) { throw "INApestAnalytical.R SHA-256 mismatch: $Actual" }
Write-Host "PASS INApestAnalytical.R $Actual"
