#!/usr/bin/env python3
"""Unit tests for run_all.py path resolution and cmake command order."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import run_all


class ResolvePathsTests(unittest.TestCase):
    def test_defaults_use_siblings_and_msys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            paths = run_all.resolve_paths(rr_root=rr)
            self.assertEqual(paths.rr_root, rr.resolve())
            self.assertEqual(paths.tsnnls, (parent / "ridge_tsnnls").resolve())
            self.assertEqual(paths.plcurve, (parent / "ridge_plcurve").resolve())
            self.assertEqual(paths.prefix, (parent / "ridge-prefix").resolve())
            self.assertEqual(paths.msys, run_all.MSYS2_DEFAULT.resolve())

    def test_msys_and_prefix_override(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            msys = parent / "custom-msys"
            prefix = parent / "pfx"
            paths = run_all.resolve_paths(rr_root=rr, prefix=prefix, msys=msys)
            self.assertEqual(paths.msys, msys.resolve())
            self.assertEqual(paths.prefix, prefix.resolve())


class PrereqTests(unittest.TestCase):
    def test_missing_msys_is_fatal_and_mentions_vs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            (parent / "ridge_tsnnls").mkdir()
            (parent / "ridge_tsnnls" / "CMakeLists.txt").write_text("x", encoding="utf-8")
            (parent / "ridge_plcurve").mkdir()
            (parent / "ridge_plcurve" / "CMakeLists.txt").write_text("x", encoding="utf-8")
            paths = run_all.resolve_paths(
                rr_root=rr, msys=parent / "no-msys-here"
            )
            errs = run_all.check_prereqs(paths, skip_deps=False)
            self.assertTrue(errs)
            blob = "\n".join(errs)
            self.assertIn("MSYS2", blob)
            self.assertIn("Visual Studio", blob)
            self.assertIn(run_all.MSYS2_INSTALL_URL, blob)

    def test_missing_sibling_fail_fast(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            msys = parent / "msys"
            (msys / "usr" / "bin").mkdir(parents=True)
            (msys / "usr" / "bin" / "bash.exe").write_bytes(b"")
            (msys / "mingw64" / "bin").mkdir(parents=True)
            (msys / "mingw64" / "bin" / "cmake.exe").write_bytes(b"")
            (msys / "mingw64" / "bin" / "ninja.exe").write_bytes(b"")
            paths = run_all.resolve_paths(rr_root=rr, msys=msys)
            errs = run_all.check_prereqs(paths, skip_deps=False)
            self.assertTrue(any("ridge_tsnnls" in e for e in errs))


class CMakePlanTests(unittest.TestCase):
    def test_order_tsnnls_plcurve_ridgerunner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            paths = run_all.Paths(
                rr_root=rr.resolve(),
                tsnnls=(parent / "ridge_tsnnls").resolve(),
                plcurve=(parent / "ridge_plcurve").resolve(),
                prefix=(parent / "ridge-prefix").resolve(),
                msys=(parent / "msys").resolve(),
            )
            cmds = run_all.planned_cmake_commands(paths, skip_deps=False, jobs=0)
            # 3 repos × (configure, build, install)
            self.assertEqual(len(cmds), 9)
            joined = [" ".join(c) for c in cmds]
            self.assertIn("ridge_tsnnls", joined[0])
            self.assertIn("ridge_plcurve", joined[3])
            self.assertIn(str(rr.resolve()), joined[6])
            self.assertTrue(all("-G" in c and "Ninja" in c for c in cmds[0::3]))

    def test_skip_deps_only_ridgerunner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            paths = run_all.Paths(
                rr_root=rr.resolve(),
                tsnnls=(parent / "ridge_tsnnls").resolve(),
                plcurve=(parent / "ridge_plcurve").resolve(),
                prefix=(parent / "ridge-prefix").resolve(),
                msys=(parent / "msys").resolve(),
            )
            cmds = run_all.planned_cmake_commands(paths, skip_deps=True, jobs=4)
            self.assertEqual(len(cmds), 3)
            self.assertIn(str(rr.resolve()), " ".join(cmds[0]))
            self.assertNotIn("ridge_tsnnls", " ".join(cmds[0]))
            self.assertIn("-j", cmds[1])
            self.assertIn("4", cmds[1])


class MainDryRunTests(unittest.TestCase):
    def test_dry_run_returns_zero_when_prereqs_ok(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp)
            rr = parent / "ridgerunner"
            rr.mkdir()
            for name in ("ridge_tsnnls", "ridge_plcurve"):
                d = parent / name
                d.mkdir()
                (d / "CMakeLists.txt").write_text("project(x)\n", encoding="utf-8")
            msys = parent / "msys"
            (msys / "usr" / "bin").mkdir(parents=True)
            (msys / "usr" / "bin" / "bash.exe").write_bytes(b"")
            (msys / "mingw64" / "bin").mkdir(parents=True)
            (msys / "mingw64" / "bin" / "cmake.exe").write_bytes(b"")
            (msys / "mingw64" / "bin" / "ninja.exe").write_bytes(b"")

            # Point resolve_paths at our fake tree via chdir + argv overrides.
            argv = ["--dry-run", "--msys", str(msys), "--prefix", str(parent / "pfx")]
            # resolve_paths uses __file__ parent as rr_root — patch it.
            with mock.patch.object(run_all, "resolve_paths") as rp:
                rp.return_value = run_all.Paths(
                    rr_root=rr.resolve(),
                    tsnnls=(parent / "ridge_tsnnls").resolve(),
                    plcurve=(parent / "ridge_plcurve").resolve(),
                    prefix=(parent / "pfx").resolve(),
                    msys=msys.resolve(),
                )
                rc = run_all.main(argv)
            self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main()
