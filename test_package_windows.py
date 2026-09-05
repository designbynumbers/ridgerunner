#!/usr/bin/env python3
"""Tests for Windows zip packaging expectations."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path


REPO = Path(__file__).resolve().parent
PS1 = REPO / "scripts" / "package-windows.ps1"
BUNDLE = REPO / "windows" / "bundle"

REQUIRED_ZIP_MEMBERS = {
    "ridgerunner-2.3.1-windows-x64/ridgerunner.cmd",
    "ridgerunner-2.3.1-windows-x64/run_knotplot_txt.py",
    "ridgerunner-2.3.1-windows-x64/install-user-path.ps1",
    "ridgerunner-2.3.1-windows-x64/README.md",
    "ridgerunner-2.3.1-windows-x64/bin/ridgerunner.exe",
    "ridgerunner-2.3.1-windows-x64/bin/ridgerunner_multithread.exe",
}


class BundleLayoutTests(unittest.TestCase):
    def test_bundle_files_exist(self) -> None:
        for name in (
            "ridgerunner.cmd",
            "run_knotplot_txt.py",
            "install-user-path.ps1",
            "README.md",
        ):
            self.assertTrue((BUNDLE / name).is_file(), name)

    def test_readme_mentions_msys_and_vs(self) -> None:
        text = (BUNDLE / "README.md").read_text(encoding="utf-8")
        self.assertIn("MSYS2", text)
        self.assertIn("Visual Studio", text)
        self.assertIn("Run this zip", text)


class PackageScriptTests(unittest.TestCase):
    def test_package_creates_expected_manifest(self) -> None:
        prefix_bin = Path(r"C:\workspace\projects\ridge-prefix\bin")
        if not (prefix_bin / "ridgerunner.exe").is_file():
            self.skipTest("ridge-prefix binaries not present")
        if not PS1.is_file():
            self.fail(f"missing {PS1}")

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            cmd = [
                "powershell",
                "-NoProfile",
                "-ExecutionPolicy",
                "Bypass",
                "-File",
                str(PS1),
                "-Prefix",
                str(prefix_bin.parent),
                "-RepoRoot",
                str(REPO),
                "-OutDir",
                str(out),
                "-Version",
                "2.3.1",
            ]
            subprocess.run(cmd, check=True)
            zip_path = out / "ridgerunner-2.3.1-windows-x64.zip"
            manifest_path = out / "ridgerunner-2.3.1-windows-x64.manifest.json"
            self.assertTrue(zip_path.is_file())
            self.assertTrue(manifest_path.is_file())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
            files = set(manifest["files"])
            for member in (
                "ridgerunner.cmd",
                "run_knotplot_txt.py",
                "README.md",
                "bin/ridgerunner.exe",
                "bin/ridgerunner_multithread.exe",
            ):
                self.assertIn(member, files)
            with zipfile.ZipFile(zip_path) as zf:
                names = set(zf.namelist())
            for member in REQUIRED_ZIP_MEMBERS:
                self.assertIn(member, names)


if __name__ == "__main__":
    unittest.main()
