/*

Portability.h

Part of the RidgeRunner project. Machine-specific headers and small
cross-platform helpers (directory create/remove, temp paths, setenv).

*/

#ifndef RIDGERUNNER_PORTABILITY_H
#define RIDGERUNNER_PORTABILITY_H

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdarg.h>
#include <stdbool.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#  include <direct.h>
#  include <io.h>
#  include <process.h>
#  ifndef getpid
#    define getpid _getpid
#  endif
#  ifndef S_ISDIR
#    define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#  endif
#else
#  include <unistd.h>
#endif

#include <dirent.h>

#include "cblas.h"
#include "lapacke.h"
#include "libtsnnls/tsnnls.h"
#include "libtsnnls/lsqr.h"

#ifndef TRUE
#  define TRUE 1
#endif
#ifndef FALSE
#  define FALSE 0
#endif

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#else
  #ifdef HAVE_MALLOC_MALLOC_H
    #include <malloc/malloc.h>
  #endif
#endif

/* Cross-platform helpers implemented in portability.c */
void rr_mkdir_or_die(const char *path, const char *file, int line);
void rr_rmtree_or_die(const char *path, const char *file, int line);
int  rr_setenv(const char *name, const char *value, int overwrite);
void rr_temp_path(char *buf, size_t buflen, const char *filename);

#endif /* RIDGERUNNER_PORTABILITY_H */
