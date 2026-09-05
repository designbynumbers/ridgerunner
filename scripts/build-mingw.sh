#!/usr/bin/env bash
# Build tsnnls → plCurve → ridgerunner for MinGW64 and install to a shared prefix.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PREFIX="${PREFIX:-$ROOT/ridge-prefix}"
TSNNLS_SRC="${TSNNLS_SRC:-$ROOT/ridge_tsnnls}"
PLCURVE_SRC="${PLCURVE_SRC:-$ROOT/ridge_plcurve}"
RR_SRC="${RR_SRC:-$ROOT/ridgerunner}"

export PATH="/mingw64/bin:${PATH}"
export PKG_CONFIG_PATH="/mingw64/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

echo "Install prefix: $PREFIX"
mkdir -p "$PREFIX"

build_one() {
  local src="$1" build="$2"
  cmake -S "$src" -B "$build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_PREFIX_PATH="$PREFIX"
  cmake --build "$build"
  cmake --install "$build"
}

build_one "$TSNNLS_SRC" "$TSNNLS_SRC/build-mingw"
build_one "$PLCURVE_SRC" "$PLCURVE_SRC/build-mingw"

cmake -S "$RR_SRC" -B "$RR_SRC/build-mingw" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_PREFIX_PATH="$PREFIX"
cmake --build "$RR_SRC/build-mingw"
cmake --install "$RR_SRC/build-mingw"

echo
echo "Build complete. Add to PATH before running:"
echo "  export PATH=\"$PREFIX/bin:/mingw64/bin:\$PATH\""
echo "Then e.g.:"
echo "  ridgerunner -a -s 1000 knot.vect"
