@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Entry point for ridgerunner on Windows.
rem - For .txt KnotPlot XYZ: convert -> run -> write {stem}_rr_{tag}.txt
rem - Otherwise: forward to ridgerunner.exe (VECT mode)
rem
rem This directory (scripts\windows) should be on User PATH (see install-user-path.ps1).
rem Do NOT put this .cmd next to ridgerunner.exe — .exe wins in PATHEXT.

set "SCRIPT_DIR=%~dp0"
rem scripts\windows -> scripts -> ridgerunner repo root
for %%I in ("%SCRIPT_DIR%..\..") do set "REPO=%%~fI"

if defined RIDGERUNNER_PREFIX (
  set "PREFIX=%RIDGERUNNER_PREFIX%"
) else (
  for %%I in ("%REPO%\..\ridge-prefix") do set "PREFIX=%%~fI"
)

if defined MINGW64_BIN (
  set "MINGW_BIN=%MINGW64_BIN%"
) else (
  set "MINGW_BIN=C:\msys64\mingw64\bin"
)

set "PATH=%PREFIX%\bin;%MINGW_BIN%;%PATH%"

rem Optional override: RIDGERUNNER_EXE (e.g. ridgerunner_multithread.exe)
if defined RIDGERUNNER_EXE (
  set "RR_EXE=%RIDGERUNNER_EXE%"
) else (
  set "RR_EXE=%PREFIX%\bin\ridgerunner.exe"
)

rem No args / wrapper help only -> Python wrapper help (not native exe).
if "%~1"=="" goto :wrapper_help
if /I "%~1"=="-h" if "%~2"=="" goto :wrapper_help
if /I "%~1"=="--help" if "%~2"=="" goto :wrapper_help
if "%~1"=="/?" if "%~2"=="" goto :wrapper_help

set "HAS_TXT="
for %%A in (%*) do (
  echo %%~A| findstr /I /E ".txt" >nul && set "HAS_TXT=1"
)

if defined HAS_TXT goto :run_txt

if not exist "%RR_EXE%" (
  echo ridgerunner.cmd: not found: "%RR_EXE%" 1>&2
  echo Build/install first, or set RIDGERUNNER_PREFIX / RIDGERUNNER_EXE. 1>&2
  exit /b 1
)

"%RR_EXE%" %*
exit /b !ERRORLEVEL!

:run_txt
where python >nul 2>&1
if errorlevel 1 (
  echo ridgerunner.cmd: python not found on PATH. Install Python 3 and retry. 1>&2
  exit /b 1
)
python "%REPO%\tools\run_knotplot_txt.py" %*
exit /b !ERRORLEVEL!

:wrapper_help
where python >nul 2>&1
if errorlevel 1 (
  echo ridgerunner.cmd: python not found on PATH. Install Python 3 and retry. 1>&2
  exit /b 1
)
python "%REPO%\tools\run_knotplot_txt.py" --help
exit /b !ERRORLEVEL!
