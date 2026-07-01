# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Ridgerunner computes minimum-ropelength ("tight") knots and links by constrained gradient
descent. It reads a knot/link as a Geomview VECT polygon, then repeatedly steps the vertices
in a direction that reduces length while enforcing the thickness (self-distance + curvature)
constraint, writing the tightened curve and diagnostic logs to a `<basename>.rr/` output
directory.

## Build, test, install

This is a GNU autotools (autoconf/automake, non-recursive) C project.

```sh
./configure
make
make check      # runs the test suite
make install
```

On Apple Silicon / macOS the dependency headers and libs live under Homebrew prefixes, so
configure needs explicit paths — use the provided helper instead of bare `./configure`:

```sh
./configsonoma.sh   # ./configure with LDFLAGS/CPPFLAGS for /usr/local + /opt/homebrew + openblas
```

If `configure` is missing or `configure.ac`/`Makefile.am` changed, regenerate with
`autoreconf -i` (or `./bootstrap` if present).

### Running a single test

`make check` runs the `TESTS` list (`csteptest mrgradtest strutgradtest rownumcvctest`).
Each is an ordinary executable built into the tree root, so run one directly:

```sh
make csteptest && ./csteptest
```

`lookahead` is a `check_PROGRAM` but not in `TESTS` — it builds under `make check` but is a
manual diagnostic, not an automated pass/fail test.

## Dependencies

Non-standard libraries that must be installed (headers are checked by `configure.ac`):
`plCurve` (provides both `plCurve.h`/`octrope.h` — the curve + octree-based thickness/strut
library), `tsnnls` (sparse constrained-least-squares solver), `argtable2` (CLI parsing),
OpenBLAS (cblas + lapacke), `gsl`, and `ncurses` (optional; enables the curses display).
The easiest install path is `brew tap designbynumbers/cantarellalab && brew install ridgerunner`.

## Architecture

All three `bin_PROGRAMS` (`ridgerunner`, `strutplot`, `residual`) share the same core object
set — the tools are thin front-ends that link the stepper/geometry code. `ridgerunner.h` is
the single shared header: it declares every prototype **and** the `extern` globals for the
whole project.

Core source (`src/`):

- **`ridgerunner_main.c`** — CLI entry. Defines all the global variables, parses ~40 argtable
  options, builds the initial `search_state`, and calls `bsearch_stepper`. Also where the
  autoscale / reresolution / eq preprocessing of the input curve happens.
- **`stepper.c`** — the heart of the program (by far the largest file). Owns the main loop
  `bsearch_stepper` → `bsearch_step`, the force pipeline, and the constraint solver:
  - Force assembly: `inputForce` combines sub-forces `dlenForce` (gradient of length),
    `eqForce` (equilateralization), `spinForce`, `specialForce` (stub), and `constraintForce`.
  - `buildRigidityMatrix` assembles the sparse strut/kink/constraint matrix; `resolveForce`
    uses tsnnls to project the raw force onto the feasible cone; `stepDirection` ties it
    together. `stepScore` is the scalar objective being minimized (length when strut-free,
    ropelength when struts are active).
  - `correct_thickness` (Newton) and `correct_constraints` pull the curve back onto the
    constraint surface after an overstep.
- **`ridgerunner.h`** — the `search_state` struct is the god-object threaded through
  everything: stopping criteria, filenames, current geometry (ropelength, thickness, minrad),
  the last step's struts/minrad-locs, logging config, and the `graphing[]`/`GraphTypes` enum
  that selects which per-step logs are written.
- **`linklib_additions.c`** — plCurve helpers (resolution/length fixing, tangents, torsion,
  drawing) that arguably belong upstream in plCurve.
- **`free_edge.c`** — "strut-free region" logic: identifies vertices far from any strut and
  (with `--Timewarp`) accelerates their motion, plus strut-free residual accounting.
- **`errors.c`** — `*_or_die` wrappers (fopen/malloc/rename/...), fatal/nonfatal error
  reporting, `logprintf`, and the `dumpAxb*`/`dump*` debugging dumps of the linear system.
- **`display.c`** — runtime terminal/curses display and geomview pipe output.
- **`dlen.c`**, **`settings.c`**, **`globals.c`**, **`mangle.c`** — length-gradient helper,
  runtime settings, remaining globals, and "MangleMode" (perturb configuration instead of
  minimizing) support.

Tools (`tools/`): **`strutplot.c`** and **`residual.c`** are companion binaries;
`lookahead.c` is a diagnostic. Other files here (`ratsrunner`, `povanimate`, `runvalgrind`,
`runmassif`, scripts) are run-management and visualization helpers.

Layout convention (from the Sept-2022 restructure): `src/` = core program, `test/` = test
programs, `tools/` = companion binaries + scripts, `attic/` = retired experiments, `data/` =
a bundled KnotPlot library of loose knots/links used for testing and installed as pkgdata.

### Flat-buffer convention

The code constantly converts between a `plCurve` (per-component vertex arrays) and a "flat"
representation — one contiguous array of `plc_vector`s (or of doubles) holding all vertices of
all components in dictionary `(component, vertex)` order. `search_state.compOffsets[i]` gives
the flat index of component *i*'s first vertex. When touching force/gradient/matrix code,
respect this ordering; helpers like `dlenPos` map a `(cmp,vt)` position to its flat address.

## Output & logging

A run creates `<basename>.rr/` containing the original and `.atstart` VECT files, the running
`.final.vect`/`.final.struts`, a log, and subdirectories `logfiles/`, `vectfiles/`, and
`snapshots/`. Which per-step quantities get logged is driven by the `GraphTypes` enum and the
`state.graphing[]` / `state.logfiles[]` arrays. `vectfiles/` and `logfiles/` are size-capped
(`--MaxVectDirSize`, `--MaxLogSize`); snapshots are taken every `--SnapshotInterval` steps.
See `doc/rr_file_behavior.txt` and `README` for the full description.

## Conventions

- `#define DEBUG 1` in `ridgerunner.h` keeps asserts enabled. `CURSES_DISPLAY` is off by
  default (plain stdout display); define it to enable the ncurses screen.
- Don't hand-edit `config.h`, `Makefile`, or `configure` — they're generated. Edit
  `configure.ac` / `Makefile.am` and re-run autotools. When adding a source file, add it to
  the relevant `*_SOURCES` list(s) in `Makefile.am` (note the shared core list is duplicated
  across each program/test target).
- The library version in `configure.ac` exists in two synced forms (`2.2.2` and `2:2:2`);
  keep them in sync when bumping.
