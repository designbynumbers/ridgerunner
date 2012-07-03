/*

Portability.h 

Part of the RidgeRunner project. The idea of this file is that all of
the machine-specific initializations and headers needed by the code in
general should appear in this file. This includes making sure that the
needed BLAS/Lapack calls appear in this file, as well as any functions
which need to be provided if system libraries are lacking.

*/

#include "config.h"

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#ifdef HAVE_TIME_H
  #include <time.h>
#endif

#ifdef HAVE_SYS_TYPES_H
  #include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif

#ifdef HAVE_ERRNO_H
  #include <errno.h>
#endif

#ifdef HAVE_DIRENT_H
  #include <dirent.h>
#endif 

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif

#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif

#ifdef HAVE_STDARG_H
  #include <stdarg.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include <stdbool.h>
#endif

#ifdef HAVE_UNISTD_H
  #include <unistd.h>
#endif

#ifdef HAVE_CLAPACK_H
  #include <clapack.h>
#else 
  #ifdef HAVE_ATLAS_CLAPACK_H
     #include <atlas/clapack.h>
  #else
     #ifdef HAVE_VECLIB_CLAPACK_H
       #include <vecLib/clapack.h>
     #endif
  #endif
#endif


#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif

#ifdef HAVE_LIBTSNNLS_TSNNLS_H    /* We are including a built tsnnls */
  #include "libtsnnls/tsnnls.h"
  #include "libtsnnls/lsqr.h"
#else   /* We are install tsnnls from a nonstandard location. Assume that there's an -I pointing directly to both headers. */
  #include "tsnnls.h"
  #include "lsqr.h"
#endif
