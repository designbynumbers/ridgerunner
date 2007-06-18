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
AC_CHECK_FUNCS([memset memcpy])

# Checks for FORTRAN.
AC_PROG_F77		
AC_F77_LIBRARY_LDFLAGS

# Checks for Darwin.
AC_CANONICAL_HOST
AC_MSG_CHECKING(how to link with BLAS and LAPACK)

echo ${host_os}

case "${host_os}" in
  darwin*)
    echo -framework vecLib
    AM_CONDITIONAL([DARWIN],[test 1 == 1])
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

    # Now look for a full (or minimal ATLAS) lapack
    AC_F77_FUNC(dpotrf)
    ACX_LAPACK([],AC_CHECK_LIB([lapack],$dpotrf,[atlas_lapack=yes 
       LAPACK_LIBS="-llapack -lcblas -lf77blas -latlas"
       AC_DEFINE_UNQUOTED([HAVE_ATLAS_LAPACK],[1],[ATLAS (instead of full) LAPACK])
       ],[atlas_lapack=no],[-lcblas -lf77blas -latlas $FLIBS -lm]))
    AC_MSG_CHECKING([for minimal ATLAS lapack])
    echo $atlas_lapack
    
  esac

# We'll need to call F77 functions from C, find out how to mangle their names
AC_F77_FUNC([daxpy],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="daxpy_"
fi
AC_DEFINE_UNQUOTED([DAXPY_F77],[$fname],[Name of daxpy after mangling])
AC_MSG_CHECKING([name of daxpy after mangling])
echo "$fname"

AC_F77_FUNC([dsyrk],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dsyrk_"
fi
AC_DEFINE_UNQUOTED([DSYRK_F77],[$fname],[Name of dsyrk after mangling])
AC_MSG_CHECKING([name of dsyrk after mangling])
echo "$fname"

AC_F77_FUNC([dsymv],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dsymv_"
fi
AC_DEFINE_UNQUOTED([DSYMV_F77],[$fname],[Name of dsymv after mangling])
AC_MSG_CHECKING([name of dsymv after mangling])
echo "$fname"

# dgemm

AC_F77_FUNC([dgemm],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dgemm_"
fi
AC_DEFINE_UNQUOTED([DGEMM_F77],[$fname],[Name of dgemm after mangling])
AC_MSG_CHECKING([name of dgemm after mangling])
echo "$fname"

# dtrsm
 
AC_F77_FUNC([dtrsm],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dtrsm_"
fi
AC_DEFINE_UNQUOTED([DTRSM_F77],[$fname],[Name of dtrsm after mangling])
AC_MSG_CHECKING([name of dtrsm after mangling])
echo "$fname"

# dpotrf

AC_F77_FUNC([dpotrf],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dpotrf_"
fi
AC_DEFINE_UNQUOTED([DPOTRF_F77],[$fname],[Name of dpotrf after mangling])
AC_MSG_CHECKING([name of dpotrf after mangling])
echo "$fname"

# dnrm2

AC_F77_FUNC([dnrm2],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dnrm2_"
fi
AC_DEFINE_UNQUOTED([DNRM2_F77],[$fname],[Name of dnrm2 after mangling])
AC_MSG_CHECKING([name of dnrm2 after mangling])
echo "$fname"

# ddot

AC_F77_FUNC([ddot],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="ddot_"
fi
AC_DEFINE_UNQUOTED([DDOT_F77],[$fname],[Name of ddot after mangling])
AC_MSG_CHECKING([name of ddot after mangling])
echo "$fname"

# dscal

AC_F77_FUNC([dscal],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dscal_"
fi
AC_DEFINE_UNQUOTED([DSCAL_F77],[$fname],[Name of dscal after mangling])
AC_MSG_CHECKING([name of dscal after mangling])
echo "$fname"

# Checks for header files.
AC_CHECK_HEADERS([clapack.h atlas/clapack.h vecLib/clapack.h])

# Checks for system functions.
AC_CHECK_FUNCS([asctime localtime time]) 

# Checks for libraries.
AC_CHECK_LIB(m,sqrt)
AC_CHECK_LIB(plCurve,plc_new)
AC_CHECK_LIB(octrope,octrope)
AC_CHECK_LIB(argtable2,arg_parse,[have_argtable2=true],[have_argtable2=false])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h math.h malloc.h stdarg.h ctype.h stdio.h assert.h \
                  string.h stdbool.h float.h wait.h sys/wait.h unistd.h \
                  time.h octrope.h plCurve.h unistd.h signal.h time.h \
		  sys/types.h sys/stat.h sys/dir.h])

# Checks for typedefs, structures, and compiler characteristics.


# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT