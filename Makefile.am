#
# Updated to non-recursive build and new directory structure Sept. 2022. -- JasonC
# Modernized 2026: static convenience library, generated data list, pkg-config OpenBLAS.
#


ACLOCAL_AMFLAGS = -I m4
dist_man1_MANS = doc/strutplot.1 doc/ridgerunner.1

# Developer helper scripts (not needed to build from a release tarball, which already
# ships ./configure, but handy to have travel with the sources).
EXTRA_DIST = autogen.sh configsonoma.sh

# The bundled KnotPlot library of loose example knots, installed as package data. The file
# list lives in data/knotdata.am, which autogen.sh regenerates from the contents of data/
# so it can never drift out of sync with the directory.
include data/knotdata.am

AM_CFLAGS = -g -finline-functions -Wall $(OPENBLAS_CFLAGS)

# The core src/ sources are shared by ridgerunner and every companion tool and test. Compile
# them once into a static convenience library rather than repeating the list per target.
# (strutplot is deliberately NOT built from this library: it uses only a small subset of the
# core objects, listed explicitly below.)
noinst_LIBRARIES = librrcore.a
librrcore_a_SOURCES = \
	src/errors.c src/errors.h \
	src/dlen.c src/dlen.h \
	src/linklib_additions.c \
	src/settings.c \
	src/stepper.c \
	src/display.c \
	src/free_edge.c \
	src/globals.c \
	src/mangle.c src/mangle.h \
	src/ridgerunner.h \
	src/portability.h

bin_PROGRAMS = ridgerunner strutplot residual

ridgerunner_SOURCES = src/ridgerunner_main.c
ridgerunner_LDADD = librrcore.a

residual_SOURCES = tools/residual.c
residual_LDADD = librrcore.a

# strutplot needs only a handful of the core objects, so it keeps its own explicit list.
strutplot_SOURCES = tools/strutplot.c \
	src/free_edge.c \
	src/globals.c \
	src/ordie.c

TESTS = csteptest mrgradtest strutgradtest rownumcvctest
check_PROGRAMS = csteptest mrgradtest lookahead strutgradtest rownumcvctest

lookahead_SOURCES = tools/lookahead.c
lookahead_LDADD = librrcore.a

csteptest_SOURCES = test/csteptest.c
csteptest_LDADD = librrcore.a

mrgradtest_SOURCES = test/mrgradtest.c
mrgradtest_LDADD = librrcore.a

strutgradtest_SOURCES = test/strutgradtest.c
strutgradtest_LDADD = librrcore.a

rownumcvctest_SOURCES = test/rownumcvctest.c
rownumcvctest_LDADD = librrcore.a
