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
Write-Host "INApest PoA Layer 2 exact internal validation"
Write-Host "Root: $Root"
& $Rscript --version
$log = Join-Path $Root "validation_output\layer2_console.log"
New-Item -ItemType Directory -Force -Path (Join-Path $Root "validation_output") | Out-Null
$old = $ErrorActionPreference
$ErrorActionPreference = "Continue"
& $Rscript (Join-Path $Root "validation\validate_PoA_layer2_exact_internal.R") $Root 2>&1 | Tee-Object -FilePath $log
$code = $LASTEXITCODE
$ErrorActionPreference = $old
if ($code -ne 0) { throw "Layer 2 validation failed with R exit code $code. See $log" }
Write-Host "PASS: Layer 2 exact internal validation completed successfully."
Write-Host "Outputs: $(Join-Path $Root 'validation_output')"
