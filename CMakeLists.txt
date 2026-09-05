cmake_minimum_required(VERSION 3.16)
project(ridgerunner VERSION 2.3.1 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

include(CheckIncludeFile)
include(CheckFunctionExists)
include(GNUInstallDirs)

find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)
find_package(GSL REQUIRED)
find_package(OpenMP)

find_path(PLCURVE_INCLUDE_DIR plCurve.h
  HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES include)
find_library(PLCURVE_LIBRARY NAMES plCurve
  HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES lib)
find_library(JUDY_LIBRARY NAMES Judy
  HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES lib)

find_path(TSNNLS_INCLUDE_DIR libtsnnls/tsnnls.h
  HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES include)
find_library(TSNNLS_LIBRARY NAMES tsnnls
  HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH
  PATH_SUFFIXES lib)

find_path(ARGTABLE2_INCLUDE_DIR argtable2.h)
find_library(ARGTABLE2_LIBRARY NAMES argtable2)

if(NOT PLCURVE_INCLUDE_DIR OR NOT PLCURVE_LIBRARY)
  message(FATAL_ERROR "plCurve not found. Set CMAKE_PREFIX_PATH to the ridge-prefix install.")
endif()
if(NOT TSNNLS_INCLUDE_DIR OR NOT TSNNLS_LIBRARY)
  message(FATAL_ERROR "tsnnls not found. Set CMAKE_PREFIX_PATH to the ridge-prefix install.")
endif()
if(NOT ARGTABLE2_INCLUDE_DIR OR NOT ARGTABLE2_LIBRARY)
  message(FATAL_ERROR "argtable2 not found (install mingw-w64-x86_64-argtable).")
endif()

check_function_exists(asctime HAVE_ASCTIME)
check_function_exists(localtime HAVE_LOCALTIME)
check_function_exists(time HAVE_TIME)
check_function_exists(difftime HAVE_DIFFTIME)
check_function_exists(mkstemp HAVE_MKSTEMP)
check_function_exists(fdopen HAVE_FDOPEN)
check_function_exists(remove HAVE_REMOVE)
check_function_exists(rename HAVE_RENAME)
check_function_exists(opendir HAVE_OPENDIR)
check_function_exists(getpid HAVE_GETPID)
check_function_exists(mallinfo HAVE_MALLINFO)
check_function_exists(setenv HAVE_SETENV)
check_function_exists(strncpy HAVE_STRNCPY)
check_function_exists(strstr HAVE_STRSTR)
check_include_file(malloc.h HAVE_MALLOC_H)
check_include_file(string.h HAVE_STRING_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(stdio.h HAVE_STDIO_H)
set(HAVE_CONFIG_H 1)

set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${GEN_DIR})
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
  ${GEN_DIR}/config.h
  @ONLY
)

set(RRCORE_SOURCES
  src/errors.c
  src/dlen.c
  src/linklib_additions.c
  src/settings.c
  src/stepper.c
  src/display.c
  src/free_edge.c
  src/globals.c
  src/mangle.c
  src/portability.c
)

add_library(rrcore STATIC ${RRCORE_SOURCES})
target_include_directories(rrcore
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${GEN_DIR}
    ${PLCURVE_INCLUDE_DIR}
    ${TSNNLS_INCLUDE_DIR}
    ${ARGTABLE2_INCLUDE_DIR}
)
target_compile_definitions(rrcore PUBLIC HAVE_CONFIG_H)
target_compile_options(rrcore PRIVATE -Wall)
target_link_libraries(rrcore
  PUBLIC
    ${PLCURVE_LIBRARY}
    ${TSNNLS_LIBRARY}
    ${ARGTABLE2_LIBRARY}
    PkgConfig::OPENBLAS
    GSL::gsl
    GSL::gslcblas
    m
)
if(JUDY_LIBRARY)
  target_link_libraries(rrcore PUBLIC ${JUDY_LIBRARY})
endif()
if(OpenMP_C_FOUND)
  target_link_libraries(rrcore PUBLIC OpenMP::OpenMP_C)
  message(STATUS "ridgerunner: OpenMP enabled")
else()
  message(STATUS "ridgerunner: OpenMP not found")
endif()

add_executable(ridgerunner src/ridgerunner_main.c)
target_link_libraries(ridgerunner PRIVATE rrcore)

# Same OpenMP-enabled binary under a distinct name for KnotPlot --threads runs.
add_executable(ridgerunner_multithread src/ridgerunner_main.c)
target_link_libraries(ridgerunner_multithread PRIVATE rrcore)

add_executable(residual tools/residual.c)
target_link_libraries(residual PRIVATE rrcore)

install(TARGETS ridgerunner ridgerunner_multithread residual RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS rrcore ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

enable_testing()
foreach(testname csteptest mrgradtest strutgradtest rownumcvctest)
  add_executable(${testname} test/${testname}.c)
  target_link_libraries(${testname} PRIVATE rrcore)
  add_test(NAME ${testname} COMMAND ${testname})
endforeach()
