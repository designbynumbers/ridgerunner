@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Portable bundle entry point. bin\ must sit next to this file.
rem End users: no MSYS2, no Visual Studio ? only this folder + Python 3 for .txt.

set "BUNDLE=%~dp0"
set "PATH=%BUNDLE%bin;%PATH%"

if defined RIDGERUNNER_EXE (
  set "RR_EXE=%RIDGERUNNER_EXE%"
) else (
  set "RR_EXE=%BUNDLE%bin\ridgerunner.exe"
)

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
  exit /b 1
)
"%RR_EXE%" %*
exit /b !ERRORLEVEL!

:run_txt
where python >nul 2>&1
if errorlevel 1 (
  echo ridgerunner.cmd: python not found on PATH.
  echo   Install Python 3: https://www.python.org/downloads/windows/
  echo   MSYS2 and Visual Studio are NOT required to run this zip.
  exit /b 1
)
python "%BUNDLE%run_knotplot_txt.py" %*
exit /b !ERRORLEVEL!

:wrapper_help
where python >nul 2>&1
if errorlevel 1 (
  echo ridgerunner.cmd: python not found on PATH.
  echo   Install Python 3: https://www.python.org/downloads/windows/
  exit /b 1
)
python "%BUNDLE%run_knotplot_txt.py" --help
exit /b !ERRORLEVEL!
