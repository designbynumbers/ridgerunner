# Ridgerunner portable Windows bundle

## Do I need MSYS2 or Visual Studio?

| Goal | MSYS2 | Visual Studio Community | Python 3 |
|------|-------|-------------------------|----------|
| **Run this zip** | No | No | Yes (only for `.txt` KnotPlot files) |
| **Build from source** | **Yes** (MinGW64) | **No** (MSVC unsupported) | Yes (for `run_all`) |

Visual Studio may be on your PC for other projects — that is fine. This zip does not use it.

## One-time setup

```powershell
powershell -ExecutionPolicy Bypass -File .\install-user-path.ps1
```

## Run

```bat
ridgerunner -a -s 1000 path\to\knot.txt
ridgerunner --Threads=8 -a -s 1000 path\to\knot.vect
```

Multithread binary:

```bat
set RIDGERUNNER_EXE=%CD%\bin\ridgerunner_multithread.exe
ridgerunner --Threads=8 -a -s 1000 path\to\knot.txt
```

## License

ridgerunner is GPL-3.0. Bundled MinGW/OpenBLAS/GSL runtime DLLs have their own licenses; see upstream projects.
