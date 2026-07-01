#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl

AC_PREREQ([2.57])
AC_COPYRIGHT([
------------------------------- ridgerunner ---------------------------------

                  Copyright (c) Jason Cantarella, University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
])

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(ridgerunner,2.3.0,jason.cantarella@gmail.com)
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])

#
# The version number below has a very specific meaning.  It is in the form
# a:b:c.  The first number, a, is the interface version number.  Any time we
# add or remove a function, or change the order of parameters or change the way
# that the library reads or writes files, this number needs to be incremented.
# The second number, b, is the implementation number *of this interface*.  When
# a changes, b goes back to zero.  This is the number that needs to change
# every time we put out a new library distribution.  The third number, c, is
# the interface backward-compatibility count.  Suppose that we simply add a
# needed function to the library.  The new interface still provides everything
# that the old version provided and in the same way.  So the
# backward-compatibility count gets incremented.  Of course it gets reset to 0
# when a new interface is not binary-compatible with the old one.
#
# Due to various broken ideas, we have to have two versions of this number, one
# with colons below and one with periods above.  Please keep them in sync!
#
AC_SUBST(LIB_VERSION, 2:3:0)
#AC_CONFIG_SRCDIR([configure])
AC_CONFIG_SRCDIR([src/ridgerunner_main.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC([clang gcc])

# Tools for the static convenience library (the shared src/ objects are compiled once into
# librrcore.a and linked by every program). AM_PROG_AR must appear before it is used.
AM_PROG_AR
AC_PROG_RANLIB

# pkg-config: used to locate OpenBLAS, which is keg-only under Homebrew and so is not
# found on the default header/library search path. Homebrew adds the keg's pkgconfig
# directory to PKG_CONFIG_PATH automatically, so no manual flags are needed.
PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([OPENBLAS], [openblas])

# Checks for functions
AC_CHECK_FUNCS([strncpy strstr])

# Checks for header files.

AC_CHECK_HEADER([argtable2.h], ,AC_MSG_ERROR([ridgerunner requires headers from the argtable library]))
AC_CHECK_HEADER([octrope.h], ,AC_MSG_ERROR([ridgerunner requires headers from the plcurve library]))
AC_CHECK_HEADER([libtsnnls/tsnnls.h], ,AC_MSG_ERROR([ridgerunner requires headers from the tsnnls library]))
AC_CHECK_HEADER([libtsnnls/lsqr.h], ,AC_MSG_ERROR([ridgerunner requires headers from the tsnnls library]))
AC_CHECK_HEADER([plCurve.h], ,AC_MSG_ERROR([ridgerunner requires headers from the plcurve library]))

# The OpenBLAS headers (cblas.h, lapacke.h) are provided via OPENBLAS_CFLAGS from the
# pkg-config check above, so we do not probe for them on the default search path here.

AC_CHECK_HEADERS([malloc.h malloc/malloc.h])

# Checks for system functions.
AC_CHECK_FUNCS([asctime localtime time difftime mkstemp fdopen remove rename opendir getpid mallinfo]) 

# Checks for libraries.
AC_SEARCH_LIBS([sqrt],[m])
AC_SEARCH_LIBS([arg_parse],[argtable2],[],AC_MSG_ERROR([ridgerunner requires the argtable2 library]))
AC_SEARCH_LIBS([octrope],[plCurve],[],AC_MSG_ERROR([ridgerunner requires the plcurve library]))
AC_SEARCH_LIBS([t_snnls],[tsnnls],[],AC_MSG_ERROR([ridgerunner requires the tsnnls library]))
AC_SEARCH_LIBS([initscr],[ncurses],[have_curses=true],[have_curses=false])
AC_SEARCH_LIBS([gsl_fit_linear],[gsl], ,AC_MSG_ERROR([ridgerunner requires the gsl library]),[-lgslcblas])

# OpenBLAS is located via pkg-config (PKG_CHECK_MODULES above). Append its libraries after
# the plCurve/tsnnls libraries so that static linking on Linux resolves symbols correctly
# (a library must precede the libraries it depends on).
LIBS="$LIBS $OPENBLAS_LIBS"

# At this point, we used to check for various other packages, such as gnuplot and povray.
# This kind of thing is handled by the homebrew package manager now, so we no longer try
# to implement it through autoconf.

# Checks for header files.


# Checks for typedefs, structures, and compiler characteristics.


# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo \
" ----------------------------------------------------------------------

  ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

  Installing to prefix: ${prefix}

  Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
  
  Now type 'make @<:@<target>@:>@'

    where the optional <target> is:

    all					- build all binaries
    install				- install everything
    check 				- run self-tests


  -----------------------------------------------------------------------"
