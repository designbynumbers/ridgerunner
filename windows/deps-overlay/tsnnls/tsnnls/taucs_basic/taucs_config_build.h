/* Platform configuration for TAUCS.
 * Autotools historically shipped a Darwin-hardcoded copy; this portable
 * version selects OSTYPE from the compiler target. */

#ifndef TAUCS_CONFIG_BUILD_H
#define TAUCS_CONFIG_BUILD_H

#if defined(_WIN32) || defined(__MINGW32__) || defined(__CYGWIN__)
#  define TAUCS_OSTYPE win32
#  define OSTYPE_win32
#elif defined(__APPLE__)
#  define TAUCS_OSTYPE darwin
#  define OSTYPE_darwin
#else
#  define TAUCS_OSTYPE linux
#  define OSTYPE_linux
#endif

#define TAUCS_VARIANT none
#define OSTYPE_VARIANT_none
#define TAUCS_CONFIG_DREAL
#define TAUCS_CONFIG_BASE
#define TAUCS_CONFIG_MATRIX_IO
#define TAUCS_CONFIG_AMD
#define TAUCS_CONFIG_COLAMD
#define TAUCS_CONFIG_GENMMD
#define TAUCS_CONFIG_ORDERING
#define TAUCS_CONFIG_FACTOR
#define TAUCS_CONFIG_LLT
#define TAUCS_CONFIG_OOC_LLT
#define TAUCS_CONFIG_OOC_LU
#define TAUCS_CONFIG_ADVANCED_MEMORY_OPS
#define TAUCS_CONFIG_VAIDYA
#define TAUCS_CONFIG_REC_VAIDYA
#define TAUCS_CONFIG_GREMBAN
#define TAUCS_CONFIG_INCOMPLETE_CHOL
#define TAUCS_CONFIG_ITER
#define TAUCS_CONFIG_INVERSE_FACTOR
#define TAUCS_CONFIG_TESTING_PROGRAMS
#define TAUCS_CONFIG_TEST_DIRECT
#define TAUCS_CONFIG_TEST_RUN
#define TAUCS_CONFIG_TEST_ITER
#define TAUCS_CONFIG_MATRIX_GENERATORS
#define TAUCS_CONFIG_MALLOC_STUBS

#endif /* TAUCS_CONFIG_BUILD_H */
