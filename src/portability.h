/*

Portability.h 

Part of the RidgeRunner project. The idea of this file is that all of
the machine-specific initializations and headers needed by the code in
general should appear in this file. This includes making sure that the
needed BLAS/Lapack calls appear in this file, as well as any functions
which need to be provided if system libraries are lacking.

Copyright Jason Cantarella.

This file is part of ridgerunner. ridgerunner is free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

ridgerunner is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.  You should have received a copy of the GNU General
Public License along with ridgerunner. If not, see
<https://www.gnu.org/licenses/>.

*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdarg.h>
#include <stdbool.h>
#include <unistd.h>

#include "cblas.h"
#include "lapacke.h"
#include "libtsnnls/tsnnls.h"
#include "libtsnnls/lsqr.h"

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#else
  #ifdef HAVE_MALLOC_MALLOC_H
    #include<malloc/malloc.h>
  #endif
#endif

