#Requires -Version 5.1
<#
.SYNOPSIS
  Build ridgerunner-*-windows-x64.zip from a MinGW install prefix + windows/bundle.

.PARAMETER Prefix
  ridge-prefix with bin\ridgerunner*.exe

.PARAMETER RepoRoot
  ridgerunner repo root (contains windows\bundle)

.PARAMETER MingwBin
  MSYS2 mingw64\bin (for DLL collection via ntldd). Optional if DLLs already beside exes.

.PARAMETER OutDir
  Where to write the zip (default: RepoRoot\dist)
#>
param(
  [Parameter(Mandatory = $true)][string]$Prefix,
  [Parameter(Mandatory = $true)][string]$RepoRoot,
  [string]$MingwBin = "C:\msys64\mingw64\bin",
  [string]$OutDir = "",
  [string]$Version = "2.3.1"
)

$ErrorActionPreference = "Stop"

$Prefix = [IO.Path]::GetFullPath($Prefix)
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$bundleSrc = Join-Path $RepoRoot "windows\bundle"
$binSrc = Join-Path $Prefix "bin"

if (-not (Test-Path $bundleSrc)) { throw "Missing bundle: $bundleSrc" }
if (-not (Test-Path (Join-Path $binSrc "ridgerunner.exe"))) {
  throw "Missing ridgerunner.exe under $binSrc - build first (run_all.cmd)"
}
if (-not (Test-Path (Join-Path $binSrc "ridgerunner_multithread.exe"))) {
  throw "Missing ridgerunner_multithread.exe under $binSrc"
}

if (-not $OutDir) { $OutDir = Join-Path $RepoRoot "dist" }
$OutDir = [IO.Path]::GetFullPath($OutDir)
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$stageName = "ridgerunner-$Version-windows-x64"
$stage = Join-Path $OutDir $stageName
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force -Path (Join-Path $stage "bin") | Out-Null

Copy-Item (Join-Path $bundleSrc "*") $stage -Force
Copy-Item (Join-Path $binSrc "ridgerunner.exe") (Join-Path $stage "bin") -Force
Copy-Item (Join-Path $binSrc "ridgerunner_multithread.exe") (Join-Path $stage "bin") -Force
if (Test-Path (Join-Path $binSrc "residual.exe")) {
  Copy-Item (Join-Path $binSrc "residual.exe") (Join-Path $stage "bin") -Force
}

function Get-MingwDeps {
  param([string]$Exe, [string]$Mingw)
  $ntldd = Join-Path $Mingw "ntldd.exe"
  $deps = New-Object "System.Collections.Generic.HashSet[string]"
  if (Test-Path $ntldd) {
    $old = $env:PATH
    $env:PATH = "$Mingw;$old"
    try {
      $out = & $ntldd -R $Exe 2>$null
      foreach ($line in $out) {
        if ($line -match '=>\s+(.+?)\s+\(') {
          $p = $Matches[1].Trim()
          if ($p -and (Test-Path $p) -and ($p -match 'mingw64')) {
            [void]$deps.Add([IO.Path]::GetFullPath($p))
          }
        }
      }
    } finally {
      $env:PATH = $old
    }
  }
  return $deps
}

$fallbackDlls = @(
  "libgomp-1.dll",
  "libgfortran-5.dll",
  "libquadmath-0.dll",
  "libgcc_s_seh-1.dll",
  "libwinpthread-1.dll",
  "libopenblas.dll",
  "libgsl-28.dll",
  "libgslcblas-0.dll"
)

$collected = New-Object "System.Collections.Generic.HashSet[string]"
foreach ($exeName in @("ridgerunner.exe", "ridgerunner_multithread.exe", "residual.exe")) {
  $exe = Join-Path $binSrc $exeName
  if (-not (Test-Path $exe)) { continue }
  foreach ($d in (Get-MingwDeps -Exe $exe -Mingw $MingwBin)) {
    [void]$collected.Add($d)
  }
}

if ($collected.Count -eq 0) {
  Write-Host "ntldd found no mingw deps (or missing); using fallback DLL list from $MingwBin"
  foreach ($name in $fallbackDlls) {
    $p = Join-Path $MingwBin $name
    if (Test-Path $p) { [void]$collected.Add([IO.Path]::GetFullPath($p)) }
  }
  foreach ($name in $fallbackDlls) {
    $p = Join-Path $binSrc $name
    if (Test-Path $p) { [void]$collected.Add([IO.Path]::GetFullPath($p)) }
  }
}

foreach ($dll in $collected) {
  Copy-Item $dll (Join-Path $stage "bin") -Force
  Write-Host ("DLL: " + (Split-Path $dll -Leaf))
}

$zipPath = Join-Path $OutDir ($stageName + ".zip")
if (Test-Path $zipPath) { Remove-Item -Force $zipPath }
Compress-Archive -Path $stage -DestinationPath $zipPath -Force
Write-Host "Wrote $zipPath"

$manifestFiles = @(Get-ChildItem $stage -Recurse -File | ForEach-Object {
  $_.FullName.Substring($stage.Length + 1).Replace("\", "/")
})
$manifest = [ordered]@{
  version = $Version
  zip = $zipPath
  files = $manifestFiles
}
$manifestPath = Join-Path $OutDir ($stageName + ".manifest.json")
($manifest | ConvertTo-Json -Depth 4) | Set-Content -Encoding utf8 $manifestPath
Write-Host "Wrote $manifestPath"
