#!/usr/bin/env python3
"""
Run ridgerunner on a KnotPlot XYZ .txt file (portable bundle edition).

Converts txt → VECT, invokes bin/ridgerunner.exe, then writes checkpointed:
  {stem}_rr_{tag}.txt
  {stem}_rr_{tag}.metrics.json
next to the input (plus .vect and .rr/ from ridgerunner).

Example:
    python run_knotplot_txt.py -a -s 1000 path/to/knot.txt
    → knot_rr_001k.txt + knot_rr_001k.metrics.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

Point = tuple[float, float, float]

DEFAULT_COLORS: tuple[tuple[float, float, float, float], ...] = (
    (0.20, 0.65, 1.00, 1.00),
    (1.00, 0.45, 0.25, 1.00),
    (0.35, 0.85, 0.45, 1.00),
    (0.80, 0.45, 1.00, 1.00),
    (1.00, 0.80, 0.20, 1.00),
    (0.25, 0.90, 0.85, 1.00),
)

BUNDLE_ROOT = Path(__file__).resolve().parent


def parse_xyz_txt(path: Path) -> list[list[Point]]:
    components: list[list[Point]] = []
    current: list[Point] = []

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8-sig").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line:
            if current:
                components.append(current)
                current = []
            continue
        if line.startswith("#"):
            continue

        fields = line.replace(",", " ").split()
        if len(fields) < 3:
            raise ValueError(
                f"{path}:{line_number}: expected at least three numeric fields"
            )
        try:
            point = tuple(float(value) for value in fields[:3])
        except ValueError as exc:
            raise ValueError(
                f"{path}:{line_number}: invalid XYZ coordinate: {raw_line!r}"
            ) from exc
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path}:{line_number}: coordinates must be finite")
        current.append(point)  # type: ignore[arg-type]

    if current:
        components.append(current)
    if not components:
        raise ValueError(f"{path}: no XYZ coordinates found")
    return components


def split_components(
    components: list[list[Point]],
    *,
    component_size: int | None,
    component_count: int | None,
) -> list[list[Point]]:
    if component_size is not None and component_count is not None:
        raise ValueError("use either --component-size or --component-count, not both")

    if len(components) != 1:
        if component_size is not None or component_count is not None:
            raise ValueError(
                "component splitting options are only valid when the input "
                "contains no blank-line component separators"
            )
        return components

    points = components[0]

    if component_size is not None:
        if component_size < 3:
            raise ValueError("--component-size must be at least 3")
        if len(points) % component_size != 0:
            raise ValueError(
                f"{len(points)} vertices cannot be divided into components "
                f"of size {component_size}"
            )
        return [
            points[start : start + component_size]
            for start in range(0, len(points), component_size)
        ]

    if component_count is not None:
        if component_count < 1:
            raise ValueError("--component-count must be at least 1")
        if len(points) % component_count != 0:
            raise ValueError(
                f"{len(points)} vertices cannot be divided evenly over "
                f"{component_count} components"
            )
        size = len(points) // component_count
        if size < 3:
            raise ValueError("each component must contain at least 3 vertices")
        return [
            points[index * size : (index + 1) * size]
            for index in range(component_count)
        ]

    return components


def remove_duplicate_closure(
    components: list[list[Point]], tolerance: float = 1e-9
) -> list[list[Point]]:
    cleaned: list[list[Point]] = []
    for index, component in enumerate(components, start=1):
        points = list(component)
        if len(points) >= 2 and math.dist(points[0], points[-1]) <= tolerance:
            points.pop()
        if len(points) < 3:
            raise ValueError(
                f"component {index} has fewer than 3 distinct vertices"
            )
        cleaned.append(points)
    return cleaned


def txt_to_vect_text(components: list[list[Point]], source_name: str) -> str:
    component_count = len(components)
    total_vertices = sum(len(c) for c in components)
    vertex_counts = [-len(c) for c in components]
    color_counts = [1] * component_count

    lines: list[str] = [
        "VECT",
        f"# Converted from KnotPlot XYZ: {source_name}",
        f"# {component_count} closed component(s), {total_vertices} vertices",
        f"{component_count} {total_vertices} {component_count}",
        " ".join(str(v) for v in vertex_counts),
        " ".join(str(v) for v in color_counts),
    ]
    for component_index, component in enumerate(components):
        lines.append(f"# component {component_index}")
        for x, y, z in component:
            lines.append(f"{x:.17g} {y:.17g} {z:.17g}")
    lines.append("# one RGBA color per component")
    for index in range(component_count):
        r, g, b, a = DEFAULT_COLORS[index % len(DEFAULT_COLORS)]
        lines.append(f"{r:.6g} {g:.6g} {b:.6g} {a:.6g}")
    return "\n".join(lines) + "\n"


def parse_vect_components(path: Path) -> list[list[Point]]:
    tokens: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            tokens.extend(line.split())

    if not tokens or tokens[0] != "VECT":
        raise ValueError(f"{path}: not a VECT file")

    cursor = 1
    component_count = int(tokens[cursor])
    total_vertices = int(tokens[cursor + 1])
    total_colors = int(tokens[cursor + 2])
    cursor += 3

    vertex_counts = [abs(int(tokens[cursor + i])) for i in range(component_count)]
    cursor += component_count
    cursor += component_count

    if sum(vertex_counts) != total_vertices:
        raise ValueError(f"{path}: inconsistent vertex counts in VECT header")

    components: list[list[Point]] = []
    for count in vertex_counts:
        points: list[Point] = []
        for _ in range(count):
            x = float(tokens[cursor])
            y = float(tokens[cursor + 1])
            z = float(tokens[cursor + 2])
            cursor += 3
            points.append((x, y, z))
        components.append(points)

    cursor += total_colors * 4
    if cursor > len(tokens):
        raise ValueError(f"{path}: VECT file truncated")

    return components


def components_to_txt(components: list[list[Point]]) -> str:
    blocks: list[str] = []
    for component in components:
        lines = [f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in component]
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks) + "\n"


def format_checkpoint_tag(steps: int | None) -> str:
    if steps is None:
        return "default"
    if steps >= 1000 and steps % 1000 == 0:
        return f"{steps // 1000:03d}k"
    return f"s{steps}"


def parse_stop_steps(rr_args: list[str]) -> int | None:
    """Extract -s N / --StopSteps N / --StopSteps=N from forwarded args."""
    i = 0
    while i < len(rr_args):
        arg = rr_args[i]
        if arg in ("-s", "--StopSteps") and i + 1 < len(rr_args):
            return int(rr_args[i + 1])
        if arg.startswith("--StopSteps="):
            return int(arg.split("=", 1)[1])
        if arg.startswith("-s") and len(arg) > 2 and arg[2:].isdigit():
            return int(arg[2:])
        i += 1
    return None


def find_ridgerunner_exe() -> Path:
    env_exe = os.environ.get("RIDGERUNNER_EXE")
    if env_exe:
        env_path = Path(env_exe)
        if env_path.is_file():
            return env_path
        raise FileNotFoundError(
            f"RIDGERUNNER_EXE is set but not found: {env_path}"
        )
    candidates = [
        BUNDLE_ROOT / "bin" / "ridgerunner.exe",
        BUNDLE_ROOT / "bin" / "ridgerunner",
    ]
    # Repo-layout fallback (tools/run_knotplot_txt.py → ../ridge-prefix).
    repo_root = BUNDLE_ROOT.parent if BUNDLE_ROOT.name == "tools" else None
    if repo_root is not None:
        candidates.append(repo_root.parent / "ridge-prefix" / "bin" / "ridgerunner.exe")
        candidates.append(repo_root / "build-mingw" / "ridgerunner.exe")
    env = os.environ.get("RIDGERUNNER_PREFIX")
    if env:
        candidates.append(Path(env) / "bin" / "ridgerunner.exe")
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        "ridgerunner.exe not found.\n  " + "\n  ".join(str(p) for p in candidates)
    )


def split_args(argv: list[str]) -> tuple[Path, list[str]]:
    txt_args = [a for a in argv if a.lower().endswith(".txt")]
    if len(txt_args) != 1:
        raise SystemExit(
            "expected exactly one .txt input path among the arguments"
        )
    txt_path = Path(txt_args[0]).resolve()
    if not txt_path.is_file():
        raise SystemExit(f"input not found: {txt_path}")
    forwarded = [a for a in argv if a != txt_args[0]]
    return txt_path, forwarded


def last_logfile_value(path: Path) -> float | None:
    if not path.is_file():
        return None
    last: float | None = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = raw.strip().split()
        if len(parts) >= 2:
            try:
                last = float(parts[-1])
            except ValueError:
                continue
    return last


def parse_mr_struts_from_text(text: str) -> int | None:
    # Progress lines look like: "  20 Rop: ... Str:  17 MRstruts:   0 Thi: ..."
    pattern = re.compile(r"MRstruts:\s*(\d+)", re.IGNORECASE)
    last: int | None = None
    for line in text.splitlines():
        match = pattern.search(line)
        if match:
            last = int(match.group(1))
    return last


def edge_length_stats(
    components: list[list[Point]],
) -> dict[str, float | None]:
    """Closed-curve edge lengths across all components -> min/max/ratio/mean/cv."""
    lengths: list[float] = []
    for comp in components:
        n = len(comp)
        if n < 2:
            continue
        for i in range(n):
            x0, y0, z0 = comp[i]
            x1, y1, z1 = comp[(i + 1) % n]
            dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
            lengths.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    empty: dict[str, float | None] = {
        "edge_length_min": None,
        "edge_length_max": None,
        "edge_length_ratio": None,
        "edge_length_mean": None,
        "edge_length_cv": None,
    }
    if not lengths:
        return empty
    mn = min(lengths)
    mx = max(lengths)
    mean = sum(lengths) / len(lengths)
    var = sum((L - mean) ** 2 for L in lengths) / len(lengths)
    std = math.sqrt(var)
    return {
        "edge_length_min": mn,
        "edge_length_max": mx,
        "edge_length_ratio": (mx / mn) if mn > 0.0 else None,
        "edge_length_mean": mean,
        "edge_length_cv": (std / mean) if abs(mean) > 1e-30 else None,
    }


def flatness_rms(components: list[list[Point]]) -> float | None:
    """RMS distance of all vertices to the best-fit plane (via covariance PCA)."""
    points: list[Point] = [p for comp in components for p in comp]
    n = len(points)
    if n < 3:
        return None

    cx = sum(p[0] for p in points) / n
    cy = sum(p[1] for p in points) / n
    cz = sum(p[2] for p in points) / n

    # Covariance matrix (3x3).
    c = [[0.0] * 3 for _ in range(3)]
    for x, y, z in points:
        d = (x - cx, y - cy, z - cz)
        for i in range(3):
            for j in range(3):
                c[i][j] += d[i] * d[j]
    for i in range(3):
        for j in range(3):
            c[i][j] /= n

    # Power iteration for smallest eigenvector ≈ plane normal.
    # Start from (1,1,1); repeatedly apply (I - C) to favor small eigenvalues.
    # Better: inverse iteration on C. Use simple Jacobi-ish: try coordinate axes
    # and a few random starts for the eigenvector of the smallest eigenvalue.
    def mat_vec(m: list[list[float]], v: tuple[float, float, float]) -> tuple[float, float, float]:
        return (
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        )

    def norm(v: tuple[float, float, float]) -> float:
        return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    def normalize(v: tuple[float, float, float]) -> tuple[float, float, float]:
        nrm = norm(v)
        if nrm < 1e-30:
            return (0.0, 0.0, 1.0)
        return (v[0] / nrm, v[1] / nrm, v[2] / nrm)

    # Inverse iteration: repeatedly solve roughly by multiplying with
    # adjugate of (C + eps I) — for 3x3 use cross-product trick on columns.
    # Simpler robust approach: evaluate Rayleigh quotient on a dense set of
    # directions on the sphere (Fibonacci lattice).
    best_n = (0.0, 0.0, 1.0)
    best_lambda = float("inf")
    samples = 64
    golden = (1.0 + 5.0 ** 0.5) / 2.0
    for k in range(samples):
        z = 1.0 - 2.0 * (k + 0.5) / samples
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = 2.0 * math.pi * k / golden
        v = normalize((r * math.cos(theta), r * math.sin(theta), z))
        cv = mat_vec(c, v)
        lam = cv[0] * v[0] + cv[1] * v[1] + cv[2] * v[2]
        if lam < best_lambda:
            best_lambda = lam
            best_n = v

    # Refine with a few inverse-power steps using (C + eps I)^{-1} via Cramer's.
    eps = 1e-12
    a = [
        [c[0][0] + eps, c[0][1], c[0][2]],
        [c[1][0], c[1][1] + eps, c[1][2]],
        [c[2][0], c[2][1], c[2][2] + eps],
    ]
    det = (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )
    if abs(det) > 1e-30:
        inv = [
            [
                (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det,
                (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det,
                (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
            ],
            [
                (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det,
                (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det,
                (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det,
            ],
            [
                (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det,
                (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det,
                (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det,
            ],
        ]
        v = best_n
        for _ in range(8):
            v = normalize(mat_vec(inv, v))
        best_n = v

    dist2 = 0.0
    for x, y, z in points:
        d = (x - cx) * best_n[0] + (y - cy) * best_n[1] + (z - cz) * best_n[2]
        dist2 += d * d
    return math.sqrt(dist2 / n)


# Progress lines look like: "  20 Rop: ... Str:  17 MRstruts:   0 Thi: ..."
_STEP_RE = re.compile(r"^\s*(\d+)\s+.*\bRop:", re.IGNORECASE)
_FIELD_RE = {
    "rop": re.compile(r"\bRop:\s*([0-9.eE+-]+)", re.IGNORECASE),
    "str": re.compile(r"\bStr:\s*(\d+)", re.IGNORECASE),
    "mr": re.compile(r"\bMRstruts:\s*(\d+)", re.IGNORECASE),
    "thi": re.compile(r"\bThi:\s*([0-9.eE+-]+)", re.IGNORECASE),
}


def format_progress_bar(
    step: int,
    total: int | None,
    *,
    rop: str | None = None,
    strut: str | None = None,
    mr: str | None = None,
    thi: str | None = None,
    width: int = 40,
) -> str:
    if total and total > 0:
        frac = min(max(step / total, 0.0), 1.0)
        filled = int(round(frac * width))
        bar = "#" * filled + "-" * (width - filled)
        pct = int(round(frac * 100))
        head = f"[{bar}] {step}/{total}  {pct:3d}%"
    else:
        head = f"[....] step {step}"
    bits = [head]
    if rop is not None:
        try:
            bits.append(f"Rop:{float(rop):.2f}")
        except ValueError:
            bits.append(f"Rop:{rop}")
    if strut is not None:
        bits.append(f"Str:{strut}")
    if mr is not None:
        bits.append(f"MRstruts:{mr}")
    if thi is not None:
        try:
            bits.append(f"Thi:{float(thi):.4f}")
        except ValueError:
            bits.append(f"Thi:{thi}")
    return "  ".join(bits)


def parse_progress_line(line: str) -> dict[str, str | int] | None:
    if not _STEP_RE.match(line):
        return None
    step_m = re.match(r"^\s*(\d+)\s+", line)
    if not step_m:
        return None
    out: dict[str, str | int] = {"step": int(step_m.group(1))}
    for key, pat in _FIELD_RE.items():
        m = pat.search(line)
        if m:
            out[key] = m.group(1)
    return out


def run_ridgerunner_live(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    total_steps: int | None,
    verbose: bool,
) -> tuple[int, str]:
    """Run ridgerunner, stream output, return (returncode, full_stdout)."""
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    chunks: list[str] = []
    bar_active = False

    def finish_bar() -> None:
        nonlocal bar_active
        if bar_active:
            print(flush=True)
            bar_active = False

    try:
        for line in proc.stdout:
            chunks.append(line)
            if verbose:
                print(line, end="", flush=True)
                continue

            prog = parse_progress_line(line)
            if prog is not None:
                bar = format_progress_bar(
                    int(prog["step"]),
                    total_steps,
                    rop=str(prog["rop"]) if "rop" in prog else None,
                    strut=str(prog["str"]) if "str" in prog else None,
                    mr=str(prog["mr"]) if "mr" in prog else None,
                    thi=str(prog["thi"]) if "thi" in prog else None,
                )
                print(f"\r{bar}", end="", flush=True)
                bar_active = True
                continue

            # Non-progress line: end bar then print normally
            finish_bar()
            print(line, end="", flush=True)

        finish_bar()
    finally:
        returncode = proc.wait()

    return returncode, "".join(chunks)


def build_metrics(
    *,
    source_txt: Path,
    checkpoint_tag: str,
    steps: int | None,
    rr_args: list[str],
    out_components: list[list[Point]],
    rr_dir: Path,
    final_vect: Path,
    out_txt: Path,
    stdout_text: str = "",
) -> dict[str, object]:
    logfiles = rr_dir / "logfiles"
    checkpoint_stem = final_vect.name.replace(".final.vect", "")
    log_path = rr_dir / f"{checkpoint_stem}.log"

    mr_struts = parse_mr_struts_from_text(stdout_text)
    if mr_struts is None and log_path.is_file():
        mr_struts = parse_mr_struts_from_text(
            log_path.read_text(encoding="utf-8", errors="replace")
        )

    # Prefer last progress "Str:" from stdout when strutcount.dat looks empty/zero.
    strut_from_stdout: float | None = None
    strut_pat = re.compile(r"Str:\s*(\d+)", re.IGNORECASE)
    for line in stdout_text.splitlines():
        match = strut_pat.search(line)
        if match:
            strut_from_stdout = float(match.group(1))
    strutcount = last_logfile_value(logfiles / "strutcount.dat")
    if strut_from_stdout is not None and (strutcount is None or strutcount == 0.0):
        strutcount = strut_from_stdout

    edges = edge_length_stats(out_components)

    metrics: dict[str, object] = {
        "source_txt": str(source_txt),
        "checkpoint_tag": checkpoint_tag,
        "steps": steps,
        "component_count": len(out_components),
        "vertices_per_component": [len(c) for c in out_components],
        "length": last_logfile_value(logfiles / "length.dat"),
        "thickness": last_logfile_value(logfiles / "thickness.dat"),
        "ropelength": last_logfile_value(logfiles / "ropelength.dat"),
        "residual": last_logfile_value(logfiles / "residual.dat"),
        "strutcount": strutcount,
        "edge_length_variance": last_logfile_value(logfiles / "edgelenvariance.dat"),
        "edge_length_min": edges["edge_length_min"],
        "edge_length_max": edges["edge_length_max"],
        "edge_length_ratio": edges["edge_length_ratio"],
        "edge_length_mean": edges["edge_length_mean"],
        "edge_length_cv": edges["edge_length_cv"],
        "mr_struts": mr_struts,
        "flatness": flatness_rms(out_components),
        "ridgerunner_args": rr_args,
        "final_vect": str(final_vect),
        "output_txt": str(out_txt),
        "rr_dir": str(rr_dir),
    }
    return metrics


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(
        description="Convert KnotPlot .txt -> ridgerunner -> {stem}_rr_{tag}.txt",
        add_help=False,
    )
    parser.add_argument(
        "--full-output",
        action="store_true",
        help="keep full ridgerunner output (do not pass --NoOutputFiles)",
    )
    parser.add_argument(
        "--keep-vectfiles",
        action="store_true",
        help="deprecated alias for --full-output",
    )
    parser.add_argument(
        "--component-count",
        type=int,
        default=None,
        metavar="K",
        help="split a single unseparated XYZ stream into K equal components",
    )
    parser.add_argument(
        "--component-size",
        type=int,
        default=None,
        metavar="N",
        help="split a single unseparated XYZ stream into components of size N",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print every ridgerunner progress line (default: ASCII progress bar)",
    )
    parser.add_argument(
        "--label",
        default=None,
        metavar="NAME",
        help="optional suffix on checkpoint stem: {stem}_rr_{tag}_{NAME}",
    )
    parser.add_argument("-h", "--help", action="store_true")

    ours, rest = parser.parse_known_args(argv)
    if ours.help:
        parser.print_help()
        print(
            "\nAll other arguments are forwarded to ridgerunner.exe.\n"
            "Checkpoint tag comes from -s / --StopSteps "
            "(1000->001k, 20->s20).\n"
            "Default UI is a progress bar; use --verbose for full Rop lines.\n"
            "Native ridgerunner options:  bin\\ridgerunner.exe -h\n"
            "Example:\n"
            "  ridgerunner -a -s 1000 knot.txt\n"
            "  -> knot_rr_001k.txt + knot_rr_001k.metrics.json\n"
            "Continue + label (A/B without overwrite):\n"
            "  ridgerunner -c -s 10000 --label plain prev_rr_005k.txt\n"
            "  -> prev_rr_005k_rr_010k_plain.txt\n"
        )
        return 0

    txt_path, rr_args = split_args(rest)
    full_output = ours.full_output or ours.keep_vectfiles
    if not full_output and "--NoOutputFiles" not in rr_args:
        rr_args = list(rr_args) + ["--NoOutputFiles"]

    steps = parse_stop_steps(rr_args)
    tag = format_checkpoint_tag(steps)
    out_stem = f"{txt_path.stem}_rr_{tag}"
    if ours.label:
        label = ours.label.strip()
        if not label or any(c in label for c in r'\/:*?"<>|'):
            print(f"error: invalid --label {ours.label!r}", file=sys.stderr)
            return 1
        out_stem = f"{out_stem}_{label}"

    exe = find_ridgerunner_exe()
    bin_dir = exe.parent

    workdir = txt_path.parent
    vect_path = workdir / f"{out_stem}.vect"
    rr_dir = workdir / f"{out_stem}.rr"
    final_vect = rr_dir / f"{out_stem}.final.vect"
    out_txt = workdir / f"{out_stem}.txt"
    metrics_path = workdir / f"{out_stem}.metrics.json"

    try:
        components = split_components(
            parse_xyz_txt(txt_path),
            component_size=ours.component_size,
            component_count=ours.component_count,
        )
        components = remove_duplicate_closure(components)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    vect_path.write_text(
        txt_to_vect_text(components, txt_path.name), encoding="utf-8"
    )
    print(f"Wrote VECT: {vect_path}", flush=True)
    label_note = f", label={ours.label}" if ours.label else ""
    print(f"Checkpoint: {tag} (steps={steps}{label_note})", flush=True)

    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(bin_dir), env.get("PATH", "")])

    cmd = [str(exe), *rr_args, vect_path.name]
    print("Running:", " ".join(cmd), flush=True)
    print("  cwd:", workdir, flush=True)
    returncode, stdout_text = run_ridgerunner_live(
        cmd,
        cwd=workdir,
        env=env,
        total_steps=steps,
        verbose=ours.verbose,
    )
    if returncode != 0:
        print(
            f"ridgerunner failed with exit code {returncode}",
            file=sys.stderr,
        )
        return returncode

    if not final_vect.is_file():
        print(f"missing output: {final_vect}", file=sys.stderr)
        return 1

    out_components = parse_vect_components(final_vect)
    out_txt.write_text(components_to_txt(out_components), encoding="utf-8")
    print(f"Wrote: {out_txt}", flush=True)

    metrics = build_metrics(
        source_txt=txt_path,
        checkpoint_tag=tag,
        steps=steps,
        rr_args=rr_args,
        out_components=out_components,
        rr_dir=rr_dir,
        final_vect=final_vect,
        out_txt=out_txt,
        stdout_text=stdout_text,
    )
    metrics_path.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote: {metrics_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
