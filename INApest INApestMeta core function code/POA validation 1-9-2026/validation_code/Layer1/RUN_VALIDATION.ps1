param(
  [string]$Rscript = "Rscript.exe"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Out = Join-Path $Root "validation_output"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

function Run-RTest {
  param([string]$Name, [string]$Script)
  Write-Host ""
  Write-Host "=== $Name ==="
  $Log = Join-Path $Out ($Name + ".log")
  # Windows PowerShell 5.x can promote ordinary stderr from a native command
  # (including harmless R warnings) to a terminating NativeCommandError when
  # $ErrorActionPreference is Stop. Run R with native stderr non-terminating,
  # retain both stdout/stderr in the log, and use Rscript's process exit code
  # as the authoritative PASS/FAIL signal.
  $OldErrorActionPreference = $ErrorActionPreference
  $ExitCode = 1
  try {
    $ErrorActionPreference = "Continue"
    & $Rscript (Join-Path $Root $Script) $Root 2>&1 |
      ForEach-Object { $_.ToString() } |
      Tee-Object -FilePath $Log
    $ExitCode = $LASTEXITCODE
  } finally {
    $ErrorActionPreference = $OldErrorActionPreference
  }
  if ($ExitCode -ne 0) {
    throw "$Name failed with R exit code $ExitCode. See $Log"
  }
}

Write-Host "INApest PoA observation-engine patch validation"
Write-Host "Root: $Root"
& $Rscript --version
if ($LASTEXITCODE -ne 0) { throw "Rscript is not available. Pass -Rscript with the full Rscript.exe path." }

Run-RTest "01_observation_architecture" "validation\validate_PoA_observation_architecture_v2.R"
Run-RTest "02_engine_observation_patch" "validation\validate_PoA_engine_observation_patch.R"
Run-RTest "03_default_backward_compatibility" "validation\validate_default_backward_compatibility.R"

Write-Host ""
Write-Host "PASS: all supplied native-R validation scripts completed successfully."
Write-Host "Outputs: $Out"
