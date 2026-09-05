cmake_minimum_required(VERSION 3.16)
project(libtsnnls VERSION 2.5.1 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckTypeSize)
include(GNUInstallDirs)

find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)

check_include_file(malloc.h HAVE_MALLOC_H)
check_function_exists(memset HAVE_MEMSET)
check_function_exists(memcpy HAVE_MEMCPY)
check_function_exists(mallinfo HAVE_MALLINFO)
check_function_exists(isnormal HAVE_ISNORMAL)
check_function_exists(isnan HAVE_ISNAN)

check_type_size("int" SIZEOF_INT BUILTIN_TYPES_ONLY)
check_type_size("long int" SIZEOF_LONG_INT BUILTIN_TYPES_ONLY)
if(NOT SIZEOF_INT)
  set(SIZEOF_INT 4)
endif()
if(NOT SIZEOF_LONG_INT)
  set(SIZEOF_LONG_INT 4)
endif()

set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${GEN_DIR})
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
  ${GEN_DIR}/config.h
  @ONLY
)

set(TSNNLS_SOURCES
  tsnnls/taucs_basic/taucs_logging.c
  tsnnls/taucs_basic/taucs_malloc.c
  tsnnls/taucs_basic/taucs_ccs_order.c
  tsnnls/taucs_basic/taucs_ccs_ops.c
  tsnnls/taucs_basic/taucs_vec_base.c
  tsnnls/taucs_basic/taucs_complex.c
  tsnnls/taucs_basic/colamd.c
  tsnnls/taucs_basic/amdbar.c
  tsnnls/taucs_basic/amdexa.c
  tsnnls/taucs_basic/amdtru.c
  tsnnls/taucs_basic/genmmd.c
  tsnnls/taucs_basic/taucs_timer.c
  tsnnls/taucs_basic/taucs_sn_llt.c
  tsnnls/taucs_basic/taucs_ccs_base.c
  tsnnls/taucs_basic/cblas_wrappers.c
  tsnnls/tlsqr.c
  tsnnls/tsnnls.c
  tsnnls/lsqr.c
)

add_library(tsnnls STATIC ${TSNNLS_SOURCES})
add_library(tsnnls::tsnnls ALIAS tsnnls)

target_include_directories(tsnnls
  PUBLIC
    $<BUILD_INTERFACE:${GEN_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tsnnls>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libtsnnls>
)

target_compile_options(tsnnls PRIVATE -std=gnu17 -Wall)

target_link_libraries(tsnnls PUBLIC PkgConfig::OPENBLAS m)

install(TARGETS tsnnls
  EXPORT tsnnlsTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES
  tsnnls/tsnnls.h
  tsnnls/lsqr.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libtsnnls
)
install(FILES
  tsnnls/taucs_basic/taucs.h
  tsnnls/taucs_basic/colamd.h
  tsnnls/taucs_basic/f2c.h
  tsnnls/taucs_basic/taucs_config_tests.h
  tsnnls/taucs_basic/taucs_private.h
  tsnnls/taucs_basic/taucs_config_build.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libtsnnls/taucs_basic
)
install(EXPORT tsnnlsTargets
  FILE tsnnlsTargets.cmake
  NAMESPACE tsnnls::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tsnnls
)
