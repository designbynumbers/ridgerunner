#Requires -Version 5.1
<#
.SYNOPSIS
  Add ridgerunner\scripts\windows to the current user's PATH (idempotent).

.DESCRIPTION
  After this, new terminals can run:
    ridgerunner -a -s 1000 C:\path\to\knot.txt

  The wrapper sets MinGW/ridge-prefix PATH internally for DLLs.
#>
$ErrorActionPreference = "Stop"

$wrapperDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$wrapperDir = [System.IO.Path]::GetFullPath($wrapperDir)

$userPath = [Environment]::GetEnvironmentVariable("Path", "User")
if ($null -eq $userPath) { $userPath = "" }

$parts = $userPath -split ";" | Where-Object { $_ -ne "" }
$already = $parts | Where-Object {
  try {
    [System.IO.Path]::GetFullPath($_) -ieq $wrapperDir
  } catch {
    $false
  }
}

if ($already) {
  Write-Host "Already on User PATH: $wrapperDir"
} else {
  $newPath = if ($userPath.TrimEnd(";")) { "$userPath;$wrapperDir" } else { $wrapperDir }
  [Environment]::SetEnvironmentVariable("Path", $newPath, "User")
  Write-Host "Added to User PATH: $wrapperDir"
}

# Also update this session so the user can try immediately.
$sessionParts = $env:Path -split ";" | Where-Object { $_ -ne "" }
$inSession = $sessionParts | Where-Object {
  try { [System.IO.Path]::GetFullPath($_) -ieq $wrapperDir } catch { $false }
}
if (-not $inSession) {
  $env:Path = "$wrapperDir;$env:Path"
  Write-Host "Prepended to this session PATH."
}

Write-Host ""
Write-Host "Try (new terminal, or this session):"
Write-Host "  ridgerunner -a -s 20 --NoOutputFiles path\to\knot.txt"
Write-Host "Output: path\to\knot_ridgerunned.txt"
