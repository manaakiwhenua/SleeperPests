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
Write-Host "INApest PoA Layer 3 external published validation"
Write-Host "Root: $Root"
& $Rscript --version
$out = Join-Path $Root "validation_output"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$log = Join-Path $out "layer3_console.log"
$old = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& $Rscript (Join-Path $Root "validation\validate_PoA_layer3_external.R") $Root 2>&1 | Tee-Object -FilePath $log
$code = $LASTEXITCODE
$ErrorActionPreference = $old
if ($code -ne 0) { throw "Layer 3 validation failed with R exit code $code. See $log" }
Write-Host "PASS: no reproducible Layer 3 external benchmark failed."
Write-Host "Note: SUPPORTED and REVIEW are evidence classifications, not failed tests."
Write-Host "Outputs: $out"
