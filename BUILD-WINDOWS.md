# Building and running ridgerunner on Windows

## Do I need MSYS2 or Visual Studio Community?

| Goal | Install MSYS2? | Install Visual Studio Community? | Install Python 3? |
|------|----------------|----------------------------------|-------------------|
| **Run a prebuilt zip** (GitHub Release) | **No** | **No** | **Yes** (only for KnotPlot `.txt` inputs; VECT-only can skip) |
| **Build from source** | **Yes** — MinGW64 toolchain | **No** — MSVC / `cl.exe` is **not** supported | **Yes** (for `run_all.cmd`) |

Visual Studio may already be installed for other projects. That is fine. This
tree still builds with **MSYS2 MinGW gcc**, not with a “x64 Native Tools”
Developer Prompt. If both are installed, `run_all` calls
`C:\msys64\mingw64\bin\cmake.exe` explicitly so generic / VS CMake on `PATH`
cannot hijack the build.

- MSYS2 installer: <https://www.msys2.org/> (default path `C:\msys64`)
- Python installer: <https://www.python.org/downloads/windows/>
- Do **not** install VS Community *for ridgerunner* — it will not help this build

## Quick start (build from source)

Sibling layout (same parent folder):

- `ridge_tsnnls` — [designbynumbers/tsnnls](https://github.com/designbynumbers/tsnnls) (or SST fork)
- `ridge_plcurve` — [designbynumbers/plcurve](https://github.com/designbynumbers/plcurve) (or SST fork)
- `ridgerunner` — this repo
- `ridge-prefix` — created by the build (`bin\ridgerunner.exe`, `ridgerunner_multithread.exe`)

From **cmd.exe** or PowerShell (not a VS developer shell):

```bat
cd C:\workspace\projects\ridgerunner
run_all.cmd
```

Useful flags:

```bat
run_all.cmd --fresh
run_all.cmd --skip-deps
run_all.cmd --package
run_all.cmd --msys C:\msys64 --prefix ..\ridge-prefix
```

`run_all` will:

1. Print the tooling banner (MSYS2 vs Visual Studio)
2. Fail with an install URL if MSYS2 is missing
3. `pacman -S --needed` the MinGW packages below
4. CMake/Ninja install **tsnnls → plcurve → ridgerunner**
5. Run ctest / OpenMP correctness when present
6. With `--package`, write `dist\ridgerunner-2.3.1-windows-x64.zip`

### One-time MSYS2 packages (also done by run_all)

In a **MinGW64** shell (`C:\msys64\mingw64.exe`), or via `run_all`:

```sh
pacman -S --needed \
  mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-cmake \
  mingw-w64-x86_64-ninja \
  mingw-w64-x86_64-pkgconf \
  mingw-w64-x86_64-openblas \
  mingw-w64-x86_64-gsl \
  mingw-w64-x86_64-argtable \
  mingw-w64-x86_64-ntldd
```

## Quick start (run portable zip — no MSYS2)

1. Download `ridgerunner-*-windows-x64.zip` from GitHub Releases
2. Unzip anywhere
3. Optional: `powershell -ExecutionPolicy Bypass -File .\install-user-path.ps1`
4. `ridgerunner -a -s 1000 path\to\knot.txt`

See [windows/bundle/README.md](windows/bundle/README.md). Runtime MinGW DLLs
ship inside `bin\`; you do **not** need MSYS2 or Visual Studio on the machine.

## PATH for a local (non-zip) build

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\windows\install-user-path.ps1
```

That adds `scripts\windows` to User PATH. The wrapper sets `ridge-prefix\bin`
and `C:\msys64\mingw64\bin` for DLLs.

```bat
ridgerunner -a -s 1000 C:\pad\naar\knotplot.txt
set RIDGERUNNER_EXE=C:\workspace\projects\ridge-prefix\bin\ridgerunner_multithread.exe
ridgerunner --Threads=8 -a -s 1000 C:\pad\naar\knotplot.txt
```

## Appendix: manual cmake (same as run_all)

From a MinGW64 shell, with `/mingw64/bin` on `PATH`:

```sh
export PREFIX=/c/workspace/projects/ridge-prefix
export PATH="/mingw64/bin:$PATH"

cmake -S ../ridge_tsnnls -B ../ridge_tsnnls/build-mingw -G Ninja \
  -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=$PREFIX
cmake --build ../ridge_tsnnls/build-mingw --target install

cmake -S ../ridge_plcurve -B ../ridge_plcurve/build-mingw -G Ninja \
  -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=$PREFIX
cmake --build ../ridge_plcurve/build-mingw --target install

cmake -S . -B build-mingw -G Ninja \
  -DCMAKE_PREFIX_PATH=$PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-mingw --target install
```

Or: `./scripts/build-mingw.sh` (bash-only; prefer `run_all.cmd` on Windows).

## CI / overlay

GitHub Actions (`.github/workflows/windows-release.yml`) checks out pinned
tsnnls/plcurve tags, applies `windows/deps-overlay/`, then runs the same
`run_all.py --package`. Tags matching `*-win*` skip the Ubuntu source-tarball
workflow.

## Notes

- Curses display stays off on Windows (stdout progress).
- Autotools (`./configure && make`) remains the Linux/macOS path; CMake is the
  Windows/MinGW path.
- KnotPlot campaign scripts live in SST-Workbench; they are **not** in the
  portable zip.
