#!/usr/bin/env python3
"""Build tsnnls → plcurve → ridgerunner on Windows via MSYS2 MinGW64.

Visual Studio / MSVC is NOT supported. This script refuses to use VS cmake
and tells you to install MSYS2 when it is missing.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

MSYS2_DEFAULT = Path(r"C:\msys64")
MSYS2_INSTALL_URL = "https://www.msys2.org/"
PYTHON_INSTALL_URL = "https://www.python.org/downloads/windows/"

MINGW_PACMAN_PACKAGES = (
    "mingw-w64-x86_64-toolchain",
    "mingw-w64-x86_64-cmake",
    "mingw-w64-x86_64-ninja",
    "mingw-w64-x86_64-pkgconf",
    "mingw-w64-x86_64-openblas",
    "mingw-w64-x86_64-gsl",
    "mingw-w64-x86_64-argtable",
    "mingw-w64-x86_64-ntldd",
)


@dataclass(frozen=True)
class Paths:
    rr_root: Path
    tsnnls: Path
    plcurve: Path
    prefix: Path
    msys: Path

    @property
    def mingw_bin(self) -> Path:
        return self.msys / "mingw64" / "bin"

    @property
    def bash(self) -> Path:
        return self.msys / "usr" / "bin" / "bash.exe"

    @property
    def mingw_cmake(self) -> Path:
        return self.mingw_bin / "cmake.exe"

    @property
    def mingw_ninja(self) -> Path:
        return self.mingw_bin / "ninja.exe"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build ridgerunner + deps on Windows with MSYS2 MinGW64. "
            "Do NOT use Visual Studio / MSVC for this project."
        )
    )
    p.add_argument(
        "--fresh",
        action="store_true",
        help="Delete build-mingw dirs in the three repos before building",
    )
    p.add_argument(
        "--skip-deps",
        action="store_true",
        help="Skip tsnnls/plcurve; only rebuild ridgerunner (prefix must exist)",
    )
    p.add_argument(
        "--package",
        action="store_true",
        help="After install, create the portable Windows zip",
    )
    p.add_argument(
        "--prefix",
        type=Path,
        default=None,
        help="Install prefix (default: ../ridge-prefix next to this repo)",
    )
    p.add_argument(
        "--msys",
        type=Path,
        default=None,
        help=r"MSYS2 root (default: C:\msys64 or %%MSYS2_ROOT%%)",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=0,
        help="Ninja parallelism (0 = ninja default)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned cmake commands without running them",
    )
    return p.parse_args(argv)


def resolve_paths(
    rr_root: Path | None = None,
    prefix: Path | None = None,
    msys: Path | None = None,
) -> Paths:
    root = (rr_root or Path(__file__).resolve().parent).resolve()
    parent = root.parent
    msys_root = (
        msys
        or (Path(os.environ["MSYS2_ROOT"]) if os.environ.get("MSYS2_ROOT") else None)
        or MSYS2_DEFAULT
    )
    return Paths(
        rr_root=root,
        tsnnls=(parent / "ridge_tsnnls").resolve(),
        plcurve=(parent / "ridge_plcurve").resolve(),
        prefix=(prefix or (parent / "ridge-prefix")).resolve(),
        msys=msys_root.resolve(),
    )


def _tooling_banner() -> str:
    return "\n".join(
        [
            "=" * 72,
            "TOOLING (read this once)",
            "=" * 72,
            "",
            "  Who are you?",
            "",
            "  A) I only want to RUN ridgerunner",
            "     -> Download the Windows zip from GitHub Releases.",
            "     -> You need: Python 3 (for .txt KnotPlot files).",
            "     -> You do NOT need MSYS2.",
            "     -> You do NOT need Visual Studio Community.",
            "",
            "  B) I want to BUILD from source on Windows",
            "     -> Install MSYS2 from " + MSYS2_INSTALL_URL,
            "     -> Use the MinGW64 toolchain (this script).",
            "     -> You do NOT need Visual Studio Community / MSVC.",
            "     -> Do NOT use 'x64 Native Tools' / VS cmake for this tree.",
            "",
            "  Visual Studio may already be installed for other projects - that is fine.",
            "  This build still uses MinGW gcc/g++ from MSYS2, never cl.exe.",
            "=" * 72,
        ]
    )


def check_prereqs(paths: Paths, *, skip_deps: bool) -> list[str]:
    """Return human-readable fatal errors (empty = ok). Also prints warnings."""
    print(_tooling_banner())
    errors: list[str] = []

    if not paths.msys.is_dir() or not paths.bash.is_file():
        errors.append(
            "MSYS2 not found at {root}.\n"
            "  BUILDING requires MSYS2 MinGW64 (not Visual Studio).\n"
            "  Install: {url}\n"
            "  Default path: C:\\msys64\n"
            "  Or pass --msys DIR / set MSYS2_ROOT.\n"
            "  If you only want to run binaries, use the GitHub Release zip instead.".format(
                root=paths.msys, url=MSYS2_INSTALL_URL
            )
        )
    else:
        if not paths.mingw_cmake.is_file():
            errors.append(
                f"MinGW cmake missing: {paths.mingw_cmake}\n"
                "  Open MSYS2 and run the pacman packages listed in BUILD-WINDOWS.md,\n"
                "  or re-run this script (it will try pacman -S --needed)."
            )
        if not paths.mingw_ninja.is_file():
            errors.append(f"MinGW ninja missing: {paths.mingw_ninja}")

    # Warn if a non-MinGW cmake would normally win on PATH.
    which_cmake = shutil.which("cmake")
    if which_cmake:
        cmake_path = Path(which_cmake)
        if "mingw64" not in str(cmake_path).lower():
            print(
                "NOTE: PATH currently resolves cmake to:\n"
                f"  {cmake_path}\n"
                "  That is typically CMake for Visual Studio / a generic install.\n"
                "  run_all will IGNORE it and call MinGW cmake explicitly:\n"
                f"  {paths.mingw_cmake}\n"
                "  You do not need to uninstall Visual Studio Community.\n"
            )

    vs_roots = [
        Path(r"C:\Program Files\Microsoft Visual Studio"),
        Path(r"C:\Program Files (x86)\Microsoft Visual Studio"),
    ]
    if any(p.is_dir() for p in vs_roots):
        print(
            "NOTE: Visual Studio appears installed on this machine.\n"
            "  That is OK for other work. For ridgerunner Windows builds, still use\n"
            "  MSYS2 MinGW64 only - do not open a 'Developer PowerShell for VS' and\n"
            "  expect cl.exe / MSBuild to compile this project.\n"
        )

    siblings = [] if skip_deps else [("ridge_tsnnls", paths.tsnnls), ("ridge_plcurve", paths.plcurve)]
    for name, path in siblings:
        cmake_lists = path / "CMakeLists.txt"
        if not path.is_dir():
            errors.append(
                f"Sibling repo missing: {path}\n"
                f"  Clone {name} next to ridgerunner (same parent folder)."
            )
        elif not cmake_lists.is_file():
            errors.append(
                f"Missing CMakeLists.txt in {path}\n"
                "  Checkout the Windows/CMake branch or apply windows/deps-overlay."
            )

    if not shutil.which("python") and not sys.executable:
        errors.append(f"Python 3 required. Install: {PYTHON_INSTALL_URL}")

    return errors


def mingw_env(paths: Paths) -> dict[str, str]:
    env = os.environ.copy()
    mingw_bin = str(paths.mingw_bin)
    usr_bin = str(paths.msys / "usr" / "bin")
    env["PATH"] = os.pathsep.join([mingw_bin, usr_bin, env.get("PATH", "")])
    pkg = str(paths.msys / "mingw64" / "lib" / "pkgconfig")
    prev = env.get("PKG_CONFIG_PATH", "")
    env["PKG_CONFIG_PATH"] = pkg if not prev else f"{pkg}{os.pathsep}{prev}"
    # Discourage picking up VS generators accidentally.
    env.pop("VCINSTALLDIR", None)
    env.pop("VSINSTALLDIR", None)
    return env


def planned_cmake_commands(
    paths: Paths, *, skip_deps: bool, jobs: int
) -> list[list[str]]:
    """Return the cmake configure/build/install command lists in order."""
    cmake = str(paths.mingw_cmake)
    prefix = str(paths.prefix)
    jobs_args = ["-j", str(jobs)] if jobs > 0 else []

    def one(src: Path, build: Path) -> list[list[str]]:
        return [
            [
                cmake,
                "-S",
                str(src),
                "-B",
                str(build),
                "-G",
                "Ninja",
                "-DCMAKE_BUILD_TYPE=Release",
                f"-DCMAKE_INSTALL_PREFIX={prefix}",
                f"-DCMAKE_PREFIX_PATH={prefix}",
            ],
            [cmake, "--build", str(build), *jobs_args],
            [cmake, "--install", str(build)],
        ]

    cmds: list[list[str]] = []
    if not skip_deps:
        cmds.extend(one(paths.tsnnls, paths.tsnnls / "build-mingw"))
        cmds.extend(one(paths.plcurve, paths.plcurve / "build-mingw"))
    cmds.extend(one(paths.rr_root, paths.rr_root / "build-mingw"))
    return cmds


def run(cmd: Sequence[str], *, env: dict[str, str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(list(cmd), check=True, env=env, cwd=cwd)


def pacman_ensure(paths: Paths, env: dict[str, str]) -> None:
    pkgs = " ".join(MINGW_PACMAN_PACKAGES)
    # pacman lives under usr/bin; run via bash -lc so MSYS path translation works.
    script = f"pacman -S --needed --noconfirm {pkgs}"
    run([str(paths.bash), "-lc", script], env=env)


def fresh_clean(paths: Paths, *, skip_deps: bool) -> None:
    targets = [paths.rr_root / "build-mingw"]
    if not skip_deps:
        targets = [
            paths.tsnnls / "build-mingw",
            paths.plcurve / "build-mingw",
            *targets,
        ]
    for build in targets:
        if build.is_dir():
            print(f"Removing {build}")
            shutil.rmtree(build)


def run_ctests(paths: Paths, env: dict[str], jobs: int) -> None:
    ctest = paths.mingw_bin / "ctest.exe"
    if not ctest.is_file():
        # cmake --build may still work; try cmake -E or skip softly
        print(f"WARN: ctest not found at {ctest}; skipping tests")
        return
    rr_build = paths.rr_root / "build-mingw"
    if rr_build.is_dir():
        args = [str(ctest), "--test-dir", str(rr_build), "--output-on-failure"]
        if jobs > 0:
            args.extend(["-j", str(jobs)])
        run(args, env=env)
    pl_build = paths.plcurve / "build-mingw"
    omp_test = pl_build / "octrope_omp_correctness.exe"
    if omp_test.is_file():
        run([str(omp_test)], env=env)
    elif pl_build.is_dir():
        args = [
            str(ctest),
            "--test-dir",
            str(pl_build),
            "-R",
            "octrope_omp_correctness",
            "--output-on-failure",
        ]
        try:
            run(args, env=env)
        except subprocess.CalledProcessError:
            print("WARN: plcurve octrope_omp_correctness ctest failed or missing")


def verify_exes(paths: Paths) -> None:
    bin_dir = paths.prefix / "bin"
    required = ["ridgerunner.exe", "ridgerunner_multithread.exe"]
    missing = [name for name in required if not (bin_dir / name).is_file()]
    if missing:
        raise SystemExit(
            "Build finished but missing:\n  "
            + "\n  ".join(str(bin_dir / m) for m in missing)
        )
    print(f"OK: installed into {bin_dir}")
    for name in required + ["residual.exe"]:
        p = bin_dir / name
        if p.is_file():
            print(f"  - {p}")


def print_run_hint(paths: Paths) -> None:
    wrapper = paths.rr_root / "scripts" / "windows" / "ridgerunner.cmd"
    mt = paths.prefix / "bin" / "ridgerunner_multithread.exe"
    print()
    print("Run (dev prefix):")
    print(f"  {wrapper} -a -s 20 --NoOutputFiles path\\to\\knot.txt")
    print("Or multithread:")
    print(f"  set RIDGERUNNER_EXE={mt}")
    print(f"  {wrapper} --Threads=8 -a -s 20 --NoOutputFiles path\\to\\knot.txt")
    print()
    print("End users without MSYS2: use the portable zip from GitHub Releases")
    print("(run_all.cmd --package), not this prefix.")


def run_package(paths: Paths) -> None:
    ps1 = paths.rr_root / "scripts" / "package-windows.ps1"
    if not ps1.is_file():
        raise SystemExit(f"Missing packager: {ps1}")
    cmd = [
        "powershell",
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        str(ps1),
        "-Prefix",
        str(paths.prefix),
        "-RepoRoot",
        str(paths.rr_root),
        "-MingwBin",
        str(paths.mingw_bin),
    ]
    run(cmd, env=os.environ.copy())


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    paths = resolve_paths(prefix=args.prefix, msys=args.msys)

    errors = check_prereqs(paths, skip_deps=args.skip_deps)
    if errors:
        print("\nFATAL — cannot build yet:\n", file=sys.stderr)
        for err in errors:
            print(err, file=sys.stderr)
            print(file=sys.stderr)
        return 2

    cmds = planned_cmake_commands(paths, skip_deps=args.skip_deps, jobs=args.jobs)
    if args.dry_run:
        print("Dry-run cmake plan:")
        for c in cmds:
            print(" ", " ".join(c))
        return 0

    env = mingw_env(paths)
    paths.prefix.mkdir(parents=True, exist_ok=True)

    print("Ensuring MinGW packages via pacman …")
    try:
        pacman_ensure(paths, env)
    except subprocess.CalledProcessError as exc:
        print(f"WARN: pacman failed ({exc}); continuing if packages already present")

    if args.fresh:
        fresh_clean(paths, skip_deps=args.skip_deps)

    for cmd in cmds:
        run(cmd, env=env)

    run_ctests(paths, env, args.jobs)
    verify_exes(paths)
    print_run_hint(paths)

    if args.package:
        run_package(paths)
    return 0


if __name__ == "__main__":
    sys.exit(main())
