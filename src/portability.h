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

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#ifdef HAVE_FLOAT_H
  #include <float.h>
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

