param(
  [string]$Rscript = ""
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path

if ([string]::IsNullOrWhiteSpace($Rscript)) {
  $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
  if ($cmd) { $Rscript = $cmd.Source }
  elseif (Test-Path "C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe") { $Rscript = "C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe" }
  else { throw "Rscript.exe not found. Supply -Rscript with the full path." }
}

Write-Host "INApest PoA frozen validation suite"
Write-Host "Validation root: $Root"
Write-Host "Rscript: $Rscript"

& (Join-Path $Root "Layer1\RUN_VALIDATION.ps1") -Rscript $Rscript
if ($LASTEXITCODE -ne 0) { throw "Layer 1 runner failed." }

& (Join-Path $Root "Layer2\RUN_LAYER2_VALIDATION.ps1") -Rscript $Rscript
if ($LASTEXITCODE -ne 0) { throw "Layer 2 runner failed." }

& (Join-Path $Root "Layer3\RUN_LAYER3_VALIDATION.ps1") -Rscript $Rscript
if ($LASTEXITCODE -ne 0) { throw "Layer 3 runner failed." }

Write-Host ""
Write-Host "PASS: Layers 1-3 completed with no reproducible benchmark failure."
