#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ([2.71])
AC_COPYRIGHT([
------------------------------- libplCurve.a ---------------------------------

                  Copyright (c) 2006 University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
])

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT([libplCurve],[11.2.3],[jason.cantarella@uga.edu])
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE()

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
AC_SUBST(LIB_VERSION, 11:2:3)
AC_CONFIG_SRCDIR([src/plCurve.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

#AC_USE_SYSTEM_EXTENSIONS
LT_INIT

# Checks for programs.
#AC_PROG_AWK
#AC_PROG_YACC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

# The bundled judy/planarity/planarmap sources use pre-C23 constructs (e.g. K&R function
# definitions) that modern compilers reject under their default -std=gnu23. Unless the user
# pinned a standard explicitly, build under -std=gnu17, which still accepts them.
case " $CFLAGS " in
  *" -std="*) ;;
  *)
    save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -std=gnu17"
    AC_MSG_CHECKING([whether $CC accepts -std=gnu17])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [])],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no]); CFLAGS="$save_CFLAGS"])
    ;;
esac


#AM_PROG_LEX
AC_PROG_LN_S

AC_PROG_MAKE_SET

# Checks for libraries. The syntax here is
#
# AC_SEARCH_LIBS([a function in the library],[the name of the library])
#
# if you find it, it will be automatically added to paths and the link commands
# for the project.
AC_SEARCH_LIBS([sqrt],[m])
#AC_SEARCH_LIBS([arg_parse],[argtable2],[have_argtable2=true],[have_argtable2=false])
AC_SEARCH_LIBS([cblas_dgemm],[gslcblas])
AC_SEARCH_LIBS([gsl_ran_gaussian],[gsl], ,AC_MSG_ERROR([plCurve 6.x.x and higher require gsl (for random number generation)]),[-lgslcblas])

#if test "x${have_argtable2}" = xtrue; then
#   AC_CHECK_HEADERS([argtable2.h], [], [have_argtable2=false])
#fi

#if test "x${have_argtable2}" = xfalse; then
#   AC_MSG_WARN([
#  -----------------------------------------------------------
#  Unable to find argtable2 on this system. Some demo programs
#  will not be installed in this build of plCurve.
#  -----------------------------------------------------------])
#fi

# Python extension library checks (Python, Cython, numpy)
#AC_ARG_ENABLE([python],
#        [AS_HELP_STRING([--enable-python],
#                [build Python bindings (currently deprecated until distutils dependency fixed)])],
#        [enable_python=$enableval],
#        [enable_python=no])

#AC_ARG_ENABLE([cython],
#        [AS_HELP_STRING([--enable-cython],
#        [build libpl .c Python extension source from .pyx (deprecated until distutils dependency fixed)])],
#        [enable_cython=$enableval],
#        [enable_cython=no])

#AC_ARG_ENABLE([python-otherdir],
#        [AS_HELP_STRING([--enable-python-otherdir],
#                [install Python to another directory (options: user, prefix)])],
#        [enable_pyotherdir=$enableval],
#        [enable_pyotherdir=global])

#if test "x${enable_python}" = xyes; then
#   AM_PATH_PYTHON(2.7)
#   AX_PYTHON_DEVEL(>= '2.7')
#   if test "x${enable_cython}" = xyes; then
#      CYTHON_CHECK_VERSION([0.21.0],
#           [have_cython=true],
#           [have_cython=false])
      # We might still be able to build if we're pre-cythonized
#      if test "x${have_cython}" = xfalse; then
         # Check for a Cython-generated C file to see if build should succeed
#         AC_CHECK_FILE([${srcdir}/pysrc/libpl/pdcode/diagram.c],
#            [AC_MSG_RESULT([Cython was disabled and pre-Cythonized C files were found.
#               Trying to build Python extensions without Cython.])],
#            [AC_MSG_FAILURE([Cython was disabled and no pre-Cythonized C files were found.
#               You can omit building the Python extensions by passing --disable-python.]
#         )])
#      fi
#   else
#      # Check for a Cython-generated C file to see if build should succeed
#      AC_CHECK_FILE([${srcdir}/pysrc/libpl/pdcode/diagram.c],
#        [AC_MSG_RESULT([Cython was not found, but pre-Cythonized C files were.
#            Trying to build Python extensions without Cython.])],
#        [AC_MSG_FAILURE([Neither Cython nor pre-Cythonized C files were found.
#            You can omit building the Python extensions by passing --disable-python.]
#        )])
#   fi
#   CIT_NUMPY_PYTHON_MODULE(1.8.0,
#        [have_numpy=true;
#         CIT_NUMPY_INCDIR],
#        [AC_MSG_WARN([Python bindings were not disabled but Numpy was not found.]);
#         have_numpy=false]])
#fi

#AM_CONDITIONAL([BUILD_PYTHON],[test "x${enable_python}" = xyes])
#AM_CONDITIONAL([BUILD_CYTHON],[test "x${enable_cython}" = xyes && test "x${have_cython}" = xtrue])
#AM_CONDITIONAL([INSTALL_PYTHON_USER],
#        [test "x${enable_pyotherdir}" = xyes || test "x${enable_pyotherdir}" = xuser])
#AM_CONDITIONAL([INSTALL_PYTHON_PREFIX],
#        [test "x${enable_pyotherdir}" = xprefix])
#AM_CONDITIONAL([INSTALL_PYTHON_GLOBAL],
#        [test "x${enable_pyotherdir}" = xglobal])

AC_DEFINE_UNQUOTED([SVNVERSION],["][`svnversion`]["],[the svn revision number (or "exported" if this was not built from svn)])

# Checks for header files. Just add the name of the header file that you
# require inside the [ ] in AC_CHECK_HEADERS. If you have it, you'll have
# a preprocessor symbol set in the form HAVE_MYHEADER_H, so you use it with
#
# #ifdef HAVE_MYHEADER_H
#    #include <myheader.h>
# #endif
#
# in your source files. Paths are all handled for you, so you shouldn't need
# to specify more about where it is.

#m4_warn([obsolete],
#[The preprocessor macro `STDC_HEADERS' is obsolete.
#  Except in unusual embedded environments, you can safely include all
#  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
#AC_CHECK_INCLUDES_DEFAULT
#AC_PROG_EGREP
# The C90 headers are:
#<assert.h>  <limits.h>  <signal.h>  <stdlib.h>
#<ctype.h>   <locale.h>  <stdarg.h>  <string.h>
#<errno.h>   <math.h>    <stddef.h>  <time.h>
#<float.h>   <setjmp.h>  <stdio.h>

AC_HEADER_ASSERT
AC_CHECK_HEADERS([stdbool.h wait.h sys/wait.h unistd.h \
                  sys/times.h sys/types.h fcntl.h signal.h complex.h \
		  gsl/gsl_rng.h gsl/gsl_randist.h gsl/gsl_combination.h \
		  gsl/gsl_permutation.h \
		  limits.h stdint.h strings.h])


# Signal to AUTOMAKE if we have argtable2

# AM_CONDITIONAL([HAVE_ARGTABLE2],[test "x${have_argtable2}" = xtrue])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_RESTRICT
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_STRFTIME

AC_CHECK_FUNCS([strcasestr mkstemp strcasecmp time floor sqrt ceil srandomdev sranddev random srandom])
AC_CHECK_FUNCS([bzero memset mkdir regcomp strchr strcspn strerror strrchr strstr strtol getline fstat64 memcpy])

AC_CHECK_FUNCS([bzero strchr strcspn strrchr strtol])
AC_CHECK_FUNC(getopt_long, SYS_GETOPTLONG=1, SYS_GETOPTLONG=0)
AC_CHECK_HEADERS([getopt.h])
AC_CHECK_FUNC(regcomp,     SYS_REGEX=1,      SYS_REGEX=0)
AC_CHECK_FUNC(strptime,    SYS_STRPTIME=1,   SYS_STRPTIME=0)

# Define automake conditionals
AM_CONDITIONAL(USE_SYS_GETOPTLONG, test "$SYS_GETOPTLONG" = "1")
AM_CONDITIONAL(USE_ARGREX,         test "$SYS_REGEX" = "1")
AM_CONDITIONAL(USE_ARGDATE,        test "$SYS_STRPTIME" = "1")

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile
		src/Makefile
                planarmap/Makefile
                tools/Makefile
		test/Makefile])
#		tools/pdstor_classify/Makefile])
#		tools/diagram_to_img/Makefile])

AC_CONFIG_SUBDIRS([judy])
# TODO: planarity is bundled but unused — either update to upstream v4.0.1.0
# (https://github.com/graph-algorithms/edge-addition-planarity-suite),
# integrate it into pd_embed, or remove it from the build entirely.
AC_CONFIG_SUBDIRS([planarity])
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
