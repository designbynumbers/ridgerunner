#!/bin/zsh
#
# Configure ridgerunner in a plain Homebrew dev shell (Apple Silicon or Intel).
#
# Two things a bare `./configure` needs on Homebrew that this script supplies:
#
#  1. OpenBLAS is keg-only, so it is not on the default pkg-config search path. We add its
#     keg to PKG_CONFIG_PATH; configure's PKG_CHECK_MODULES([openblas]) then supplies the
#     correct -I/-L/-fopenmp flags. No manual OpenBLAS flags are needed anymore.
#
#  2. plCurve and tsnnls do not ship pkg-config files and install onto the Homebrew prefix,
#     which the compiler does not search by default outside Homebrew's own build sandbox.
#     So we pass the prefix include/lib explicitly.
#
# NOTE: the Homebrew formula needs NONE of this -- brew's build environment adds dependency
# include/lib paths and keg-only pkgconfig dirs automatically, so the formula just runs
# `./configure`. This script is only for building by hand from a git checkout.
#
# Run ./autogen.sh first if ./configure does not yet exist.

PREFIX="$(brew --prefix)"
export PKG_CONFIG_PATH="$(brew --prefix openblas)/lib/pkgconfig:$PKG_CONFIG_PATH"

./configure CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" "$@"
