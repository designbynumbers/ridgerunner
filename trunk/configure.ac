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
AC_INIT(ridgerunner,2.0.0,cantarel@math.uga.edu)
AC_CONFIG_AUX_DIR([.])
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
AC_SUBST(LIB_VERSION, 2:0:0)
AC_CONFIG_SRCDIR([src/ridgerunner_main.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])


LT_INIT
# Checks for programs.
AC_PROG_CC
AC_CHECK_LIB([m],[cos])
AC_PROG_LIBTOOL

# Checks for functions
AC_CHECK_FUNCS([strncpy strstr])

# Checks for FORTRAN.
AC_PROG_FC([gfortran g77])		
AC_FC_LIBRARY_LDFLAGS
#LIBS="$LIBS $FCLIBS -lm"

# Checks for Darwin.
AC_CANONICAL_HOST
AC_MSG_CHECKING(how to link with BLAS and LAPACK)

case "${host_os}" in
  darwin*)
    echo -framework vecLib
    AM_CONDITIONAL([DARWIN],[test 1 = 1])
    AM_CONDITIONAL([HAVE_LAPACK],[test x1 = x1])
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

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# Checks for header files.
AC_CHECK_HEADERS([clapack.h atlas/clapack.h vecLib/clapack.h malloc.h malloc/malloc.h argtable2.h])

# Checks for system functions.
AC_CHECK_FUNCS([asctime localtime time difftime mkstemp fdopen remove rename opendir getpid mallinfo]) 

# Checks for libraries.
AC_SEARCH_LIBS([sqrt],[m])
AC_SEARCH_LIBS([arg_parse],[argtable2],[],AC_MSG_ERROR([ridgerunner requires the argtable2 library to compile]))
AC_SEARCH_LIBS([initscr],[ncurses],[have_curses=true],[have_curses=false])
AC_SEARCH_LIBS([gsl_fit_linear],[gsl],AM_CONDITIONAL([HAVEGSL],[test 1 == 1]),[AM_CONDITIONAL([HAVEGSL],[test 1 == 0]) AC_MSG_WARN([some self-tests require the GNU Scientific Library. Build will proceed, but fewer self-tests will be performed.])],[-lgslcblas])


# Now the plCurve library must work both as a prior install AND as a "sibling" library during a metapackage install.

AC_ARG_VAR([PLCINC],[The PATH where plCurve.h can be found])
AC_ARG_VAR([PLCLIB],[The PATH where libplCurve.la can be found])

# Must either specify both or none. 

if (test -n "$PLCINC" && test -z "$PLCLIB") || \
   (test -n "$PLCLIB" && test -z "$PLCINC"); then 

   AC_MSG_ERROR([Must specify both PLCINC and PLCLIB or neither.])

fi

# If no location is specified on the command line, we search for an installed copy.

if test -z "$PLCLIB"; then 
   AC_SEARCH_LIBS([plc_new],[plCurve], [have_installed_plcurve=true], [have_installed_plcurve=false])
   if test "x${have_installed_plcurve}" = xtrue; then
      AC_CHECK_HEADERS([plCurve.h], [ ], [have_installed_plcurve=false])
   fi

   if test "x${have_installed_plcurve}" = xfalse; then 
   # it may be the case that we're installing from a tarball. In this case, plCurve won't be 
   # installed yet, but it should live in a sister directory from octrope.
      if test -d "$srcdir/../plCurve/"; then 
         PLCINC='$(top_srcdir)/../plCurve'
         PLCLIB='$(top_builddir)/../plCurve'
      fi
   fi
fi

  
if (test -z "$PLCLIB" && test "x${have_installed_plCurve}" = xfalse); then 
   AC_MSG_ERROR([
-----------------------------------------------------------------------
ridgerunner requires the plCurve library to compile. You can set the
location of plCurve manually with the PLCINC and PLCLIB environment
variables, or install a copy in a standard place like /usr/local/lib
-----------------------------------------------------------------------
])
fi

# Now make these symbols available to automake.

if test -z "$PLCLIB"; then 
 AC_SUBST([PLC_LTLIB], [" "])
 AC_SUBST([PLC_INCLUDE], [" "])
else 
 AC_SUBST([PLC_LTLIB], ["$PLCLIB/libplCurve.la"])
 AC_SUBST([PLC_INCLUDE], ["-I$PLCINC"])
fi

#################### END OF PLCURVE CHECKING #####################

#AC_SEARCH_LIBS([plc_new],[plCurve], [have_installed_plcurve=true], [have_installed_plcurve=false])
#  if test "x${have_installed_plcurve}" = xtrue; then
#      AC_CHECK_HEADERS([plCurve.h], [ ], [have_installed_plcurve=false])
#  fi

# Now the octrope library must work both as a prior install AND as a "sibling" library during a metapackage install.

AC_ARG_VAR([OCTINC],[The PATH where octrope.h can be found])
AC_ARG_VAR([OCTLIB],[The PATH where liboctrope.la can be found])

# Must either specify both or none. 

if (test -n "$OCTINC" && test -z "$OCTLIB") || \
   (test -n "$OCTLIB" && test -z "$OCTINC"); then 

   AC_MSG_ERROR([Must specify both OCTINC and OCTLIB or neither.])

fi

# If no location is specified on the command line, we search for an installed copy.

if test -z "$OCTLIB"; then 
   AC_SEARCH_LIBS([octrope],[octrope], [have_installed_octrope=true], [have_installed_octrope=false],[-lplCurve])
   if test "x${have_installed_octrope}" = xtrue; then
      AC_CHECK_HEADERS([octrope.h], [ ], [have_installed_octrope=false])
   fi

   if test "x${have_installed_octrope}" = xfalse; then 
   # it may be the case that we're installing from a tarball. In this case, plCurve won't be 
   # installed yet, but it should live in a sister directory from octrope.
      if test -d "$srcdir/../octrope/src/octropesrc"; then 
         OCTINC='$(top_srcdir)/../octrope/src/octropesrc'
         OCTLIB='$(top_builddir)/../octrope/src/octropesrc'
      fi
   fi
fi

  
if (test -z "$OCTLIB" && test "x${have_installed_octrope}" = xfalse); then 
   AC_MSG_ERROR([
-----------------------------------------------------------------------
ridgerunner requires the octrope library. You can set the
location of octrope manually with the OCTINC and OCTLIB environment
variables, or install a copy in a standard place like /usr/local/lib
-----------------------------------------------------------------------
])
fi

# Now make these symbols available to automake.

if test -z "$OCTLIB"; then 
 AC_SUBST([OCT_LTLIB], [" "])
 AC_SUBST([OCT_INCLUDE], [" "])
else 
 AC_SUBST([OCT_LTLIB], ["$OCTLIB/liboctrope.la"])
 AC_SUBST([OCT_INCLUDE], ["-I$OCTINC"])
fi

#################### END OF OCTROPE CHECKING #####################

 
#AC_SEARCH_LIBS([octrope],[octrope], [have_installed_octrope=true], [have_installed_octrope=false],[-lplCurve])
#   if test "x${have_installed_octrope}" = xtrue; then
#      AC_CHECK_HEADERS([octrope.h], [ ], [have_installed_octrope=false])
#   fi

# Now the octrope library must work both as a prior install AND as a "sibling" library during a metapackage install.

AC_ARG_VAR([TSNNLSINC],[The PATH where tsnnls.h can be found])
AC_ARG_VAR([TSNNLSLIB],[The PATH where libtsnnls.la can be found])

# Must either specify both or none. 

if (test -n "$TSNNLSINC" && test -z "$TSNNLSLIB") || \
   (test -n "$TSNNLSLIB" && test -z "$TSNNLSINC"); then 

   AC_MSG_ERROR([Must specify both TSNNLSINC and TSNNLSLIB or neither.])

fi

# If no location is specified on the command line, we search for an installed copy.

if test -z "$TSNNLSLIB"; then 
   AC_SEARCH_LIBS([t_snnls],[tsnnls], [have_installed_tsnnls=true], [have_installed_tsnnls=false],[])
   if test "x${have_installed_tsnnls}" = xtrue; then
      AC_CHECK_HEADERS([libtsnnls/tsnnls.h], [ ], [have_installed_tsnnls=false])
   fi

   if test "x${have_installed_tsnnls}" = xfalse; then 
   # it may be the case that we're installing from a tarball. In this case, tsnnls won't be 
   # installed yet, but it should live in a sister directory from ridgerunner.
      if test -d "$srcdir/../tsnnls/tsnnls"; then 
         TSNNLSINC='$(top_srcdir)/../tsnnls/tsnnls'
         TSNNLSLIB='$(top_builddir)/../tsnnls/tsnnls'
      fi
   fi
fi

  
if (test -z "$TSNNLSLIB" && test "x${have_installed_tsnnls}" = xfalse); then 
   AC_MSG_ERROR([
-----------------------------------------------------------------------
ridgerunner requires the tsnnls library. You can set the
location of octrope manually with the TSNNLSINC and TSNNLSLIB environment
variables, or install a copy in a standard place like /usr/local/lib
-----------------------------------------------------------------------
])
fi

# Now make these symbols available to automake.

if test -z "$TSNNLSLIB"; then 
 AC_SUBST([TSNNLS_LTLIB], [" "])
 AC_SUBST([TSNNLS_INCLUDE], [" "])
else 
 AC_SUBST([TSNNLS_LTLIB], ["$TSNNLSLIB/libtsnnls.la"])
 AC_SUBST([TSNNLS_INCLUDE], ["-I$TSNNLSINC"])
fi

## END OF TSNNLS CHECKING
 
#AC_SEARCH_LIBS([t_snnls],[tsnnls], [have_installed_tsnnls=true], [have_installed_tsnnls=false],[])
#   if test "x${have_installed_tsnnls}" = xtrue; then
#      AC_CHECK_HEADERS([libtsnnls/tsnnls.h], [ ], [have_installed_tsnnls=false])
#   fi

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
AC_HEADER_STDC
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
  
  Package Features:
     
      plCurve location: ${PLCINC:-INSTALLED}
      octrope location: ${OCTINC:-INSTALLED}
      tsnnls location: ${TSNNLSINC:-INSTALLED}

      plCurve libraries: ${PLCLIB}
      octrope libraries: ${OCTLIB}
      tsnnls libraries: ${TSNNLSLIB}
      
      LAPACK libraries: ${LAPACK_LIBS}
      BLAS libraries: ${BLAS_LIBS}

      All libraries: ${LIBS}

  Now type 'make @<:@<target>@:>@'

    where the optional <target> is:

    all					- build all binaries
    install				- install everything
    check 				- run self-tests


  -----------------------------------------------------------------------"
