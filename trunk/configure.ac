#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ(2.57)
AC_COPYRIGHT(
------------------------------- ridgerunner ---------------------------------

                  Copyright (c) 2006 University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(ridgerunner, 1.0.1, cantarel@math.uga.edu)
AM_INIT_AUTOMAKE

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
AC_SUBST(LIB_VERSION, 1:0:1)
AC_CONFIG_SRCDIR([configure])
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for functions
AC_CHECK_FUNCS([strncpy strstr])

# Checks for FORTRAN.
AC_PROG_FC([gfortran g77])		
AC_FC_LIBRARY_LDFLAGS

# Checks for Darwin.
AC_CANONICAL_HOST
AC_MSG_CHECKING(how to link with BLAS and LAPACK)

case "${host_os}" in
  darwin*)
    echo -framework vecLib
    AM_CONDITIONAL([DARWIN],[test 1 == 1])
    AM_CONDITIONAL([HAVE_LAPACK],[test 1 == 1])
    AC_DEFINE([HAVE_DARWIN],[1],[Defined if we are in the Apple environment.])
#    LIBS="$LIBS -framework vecLib"
    AC_SUBST([LAPACK_LIBS],["-framework vecLib"])
    AC_SUBST([BLAS_LIBS],[" "])
    # $FLIBS was set above by AC_F77_PROG  
    # How can I work around the damaged g77 on Darwin by not linking lcrt2.o?
    # By editing flibs to remove -lcrt2.o, of course!  
    # (I realize that this is a total hack.)

    FLIBS=`echo $FLIBS | sed s/-lcrt2.o//g` 
    ;;
  *)
    echo linux
    AM_CONDITIONAL([DARWIN],[test 0 == 1])

    ACX_BLAS
    if test "$acx_blas_ok" = no; then 
	  AC_MSG_ERROR([ridgerunner requires an F77BLAS. The ATLAS project
will provide these on any given system if no other is
available. Configure has searched the usual locations and can't find
one. You can use --with-blas=<lib> to point to a blas library in a
nonstandard spot. This test may fail even if you have BLAS if you
don't have a FORTRAN compiler (preferably g77) on your default path.])
    fi

    AC_FC_FUNC(dpotrf)
    ACX_LAPACK(AM_CONDITIONAL([HAVE_LAPACK],[test x1 = x1]),AC_CHECK_LIB([lapack],$dpotrf,[atlas_lapack=yes 
       LAPACK_LIBS="-llapack -lcblas -lf77blas -latlas"
       AC_DEFINE_UNQUOTED([HAVE_ATLAS_LAPACK],[1],[ATLAS (instead of full) LAPACK])
       AM_CONDITIONAL([HAVE_LAPACK],[test 1 == 0])
       ],[atlas_lapack=no],[-lcblas -lf77blas -latlas $FLIBS -lm]))
    AC_MSG_CHECKING([for minimal ATLAS lapack])
    echo $atlas_lapack

 #  Now look for a full (or minimal ATLAS) lapack
 #   AC_F77_FUNC(dpotrf)
 #   ACX_LAPACK([],AC_CHECK_LIB([lapack],$dpotrf,[atlas_lapack=yes 
 #      LAPACK_LIBS="-llapack -lcblas -lf77blas -latlas"
 #      AC_DEFINE_UNQUOTED([HAVE_ATLAS_LAPACK],[1],[ATLAS (instead of full) LAPACK])
 #      ],[atlas_lapack=no],[-lcblas -lf77blas -latlas $FLIBS -lm]))
 #   AC_MSG_CHECKING([for minimal ATLAS lapack])
 #   echo $atlas_lapack
    
  esac

# Checks for header files.
AC_CHECK_HEADERS([clapack.h atlas/clapack.h vecLib/clapack.h malloc.h malloc/malloc.h])

# Checks for system functions.
AC_CHECK_FUNCS([asctime localtime time difftime mkstemp fdopen remove rename opendir getpid mallinfo]) 

# Checks for libraries.
AC_CHECK_LIB(m,sqrt)
AC_CHECK_LIB(plCurve,plc_new,[],AC_MSG_ERROR([ridgerunner requires the plCurve library]))
AC_CHECK_LIB(octrope,octrope,[],AC_MSG_ERROR([ridgerunner requires the octrope library]))
AC_CHECK_LIB(argtable2,arg_parse,[have_argtable2=true],[have_argtable2=false])
AC_CHECK_LIB(ncurses,initscr,[have_curses=true],[have_curses=false])
AC_CHECK_LIB(tsnnls,t_snnls,[],AC_MSG_ERROR([ridgerunner requires the tsnnls library]),[$LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h math.h stdio.h assert.h stdarg.h \
                  string.h stdbool.h float.h \
                  time.h octrope.h plCurve.h libtsnnls/tsnnls.h time.h \
		  sys/types.h sys/stat.h unistd.h dirent.h errno.h] )

# Checks for typedefs, structures, and compiler characteristics.


# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Add dmalloc support if desired
AM_WITH_DMALLOC

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
