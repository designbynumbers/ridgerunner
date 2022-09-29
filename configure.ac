#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ(2.57)
AC_COPYRIGHT(
------------------------------- ridgerunner ---------------------------------

                  Copyright (c) Jason Cantarella, University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(ridgerunner,2.1.1,jason.cantarella@gmail.com)
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
AC_SUBST(LIB_VERSION, 2:1:1)
AC_CONFIG_SRCDIR([src/ridgerunner_main.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

LT_INIT
# Checks for programs.
AC_PROG_CC([clang gcc])

# Checks for functions
AC_CHECK_FUNCS([strncpy strstr])

# Checks for header files.
AC_CHECK_HEADERS([cblas.h lapacke.h malloc.h malloc/malloc.h argtable2.h octrope.h libtsnnls/tsnnls.h plCurve.h])

# Checks for system functions.
AC_CHECK_FUNCS([asctime localtime time difftime mkstemp fdopen remove rename opendir getpid mallinfo]) 

# Checks for libraries.
AC_SEARCH_LIBS([sqrt],[m])
AC_SEARCH_LIBS([arg_parse],[argtable2],[],AC_MSG_ERROR([ridgerunner requires the argtable2 library]))
AC_SEARCH_LIBS([octrope],[plcurve],[],AC_MSG_ERROR([ridgerunner requires the plcurve library]))
AC_SEARCH_LIBS([t_snnls],[tsnnls],[],AC_MSG_ERROR([ridgerunner requires the tsnnls library]))
AC_SEARCH_LIBS([cblas_dgemv],[openblas], ,AC_MSG_ERROR([ridgerunner requires the OpenBlas library]))

AC_SEARCH_LIBS([initscr],[ncurses],[have_curses=true],[have_curses=false])
AC_SEARCH_LIBS([gsl_fit_linear],[gsl], ,AC_MSG_ERROR([ridgerunner requires the gsl library],[-lgslcblas])

#Checks for programs
AC_CHECK_PROG([HAVE_TUBE],[tube],[true],[false])
if test "$HAVE_TUBE" == true; then
   AC_DEFINE([HAVE_TUBE],[1],[we have the tube vecttools program])
else 
   AC_MSG_WARN([it is recommended that you install the Vecttools library with ridgerunner])
fi   

AC_CHECK_PROG([HAVE_GNUPLOT],[gnuplot],[true],[false])
if test "$HAVE_GNUPLOT" == true; then
   AC_DEFINE([HAVE_GNUPLOT],[1],[we have gnuplot on this system])
else 
   AC_MSG_WARN([it is recommended that you install gnuplot with ridgerunner])
fi   

AC_CHECK_PROG([HAVE_POVSNAP],[povsnap],[true],[false])
if test "$HAVE_POVSNAP" == true; then
   AC_DEFINE([HAVE_POVSNAP],[1],[we have povsnap and the povrayutils library on this system])
else 
   AC_MSG_WARN([ It is recommended that you install the povray utilities with ridgerunner.])
fi   

AC_CHECK_PROG([HAVE_POVRAY],[povray],[true],[false])
if test "$HAVE_POVRAY" == true; then
   AC_DEFINE([HAVE_POVRAY],[1],[we have povray on this system])
else 
   AC_MSG_WARN([it is recommended that you install povray with ridgerunner. Povray is NOT included in the metapackage because it has too many dependencies to be easily built from source on most systems. You should install it using the package manager native to your system (or from binaries).])
fi   

AC_DEFINE_UNQUOTED([SVNVERSION],["][`svnversion`]["],[the svn revision number (or "exported" if this was not built from svn)])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h math.h stdio.h assert.h stdarg.h \
                  string.h stdbool.h float.h \
                  time.h time.h \
		  sys/types.h sys/stat.h unistd.h dirent.h errno.h] )

# Checks for typedefs, structures, and compiler characteristics.


# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Add dmalloc support if desired
AM_WITH_DMALLOC

# Output
AC_CONFIG_FILES([Makefile src/Makefile])
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
