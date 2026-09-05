@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem One-click / copy-paste Windows build entrypoint.
rem Builds tsnnls → plcurve → ridgerunner via MSYS2 MinGW64.
rem Visual Studio Community / MSVC is NOT used and NOT required.

set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

where python >nul 2>&1
if errorlevel 1 (
  echo run_all.cmd: Python 3 not found on PATH.
  echo   Install: https://www.python.org/downloads/windows/
  echo   Needed for this build script and for .txt KnotPlot wrappers.
  exit /b 1
)

if not defined MSYS2_ROOT set "MSYS2_ROOT=C:\msys64"
if not exist "%MSYS2_ROOT%\usr\bin\bash.exe" (
  echo run_all.cmd: MSYS2 not found at "%MSYS2_ROOT%".
  echo.
  echo   To BUILD ridgerunner on Windows you need MSYS2 MinGW64:
  echo     https://www.msys2.org/
  echo   Default install path: C:\msys64
  echo   Or set MSYS2_ROOT / pass --msys to run_all.py
  echo.
  echo   You do NOT need Visual Studio Community for this project.
  echo   If you only want to run binaries, download the GitHub Release zip instead.
  exit /b 1
)

echo Using MSYS2 at %MSYS2_ROOT%  ^(MinGW64 — not Visual Studio^)
python "%SCRIPT_DIR%run_all.py" %*
exit /b %ERRORLEVEL%
