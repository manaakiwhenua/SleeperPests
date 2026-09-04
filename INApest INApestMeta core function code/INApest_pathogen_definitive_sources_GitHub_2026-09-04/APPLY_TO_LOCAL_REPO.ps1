param(
  [Parameter(Mandatory=$true)]
  [string]$RepoRoot
)

$ErrorActionPreference = 'Stop'
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Source = Join-Path $Here 'github_upload\INApest INApestMeta core function code'
$Target = Join-Path $RepoRoot 'INApest INApestMeta core function code'

if (-not (Test-Path $Source)) { throw "Upload source folder not found: $Source" }
if (-not (Test-Path $Target)) { throw "Repository core source folder not found: $Target" }

Get-ChildItem -Path $Source -File | ForEach-Object {
  Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $Target $_.Name) -Force
  Write-Host "Updated $($_.Name)"
}

Write-Host "Definitive commented pathogen simulation-engine source overlay complete."
