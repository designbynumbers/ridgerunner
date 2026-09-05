cmake_minimum_required(VERSION 3.16)
project(libplCurve VERSION 11.2.3 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

include(CheckIncludeFile)
include(CheckFunctionExists)
include(GNUInstallDirs)

find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(OpenMP)

set(CMAKE_REQUIRED_INCLUDES ${GSL_INCLUDE_DIRS})
check_include_file(sys/stat.h HAVE_SYS_STAT_H)
check_include_file(sys/types.h HAVE_SYS_TYPES_H)
check_include_file(unistd.h HAVE_UNISTD_H)
check_include_file(stdint.h HAVE_STDINT_H)
check_include_file(stdbool.h HAVE_STDBOOL_H)
check_include_file(malloc.h HAVE_MALLOC_H)
check_include_file(limits.h HAVE_LIMITS_H)
check_include_file(complex.h HAVE_COMPLEX_H)
check_include_file(signal.h HAVE_SIGNAL_H)
check_include_file(gsl/gsl_rng.h HAVE_GSL_GSL_RNG_H)
check_include_file(gsl/gsl_randist.h HAVE_GSL_GSL_RANDIST_H)
check_include_file(gsl/gsl_combination.h HAVE_GSL_GSL_COMBINATION_H)
check_function_exists(time HAVE_TIME)
check_function_exists(strdup HAVE_STRDUP)
set(HAVE_CONFIG_H 1)

set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${GEN_DIR})
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
  ${GEN_DIR}/config.h
  @ONLY
)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Judy.cmake)

set(PLCURVE_SOURCES
  src/plCurve.c
  src/vector.c
  src/splines.c
  src/pllmpoly02.c
  src/plcTopology.c
  src/plwhitten.c
  src/nearest_neighbor.c
  src/symmetries.c
  src/matrix.c
  src/plcRandomPolygon.c
  src/plcRandomSap.c
  src/plc_dihedrals.c
  src/plc_csv_reader.c
  src/octrope.c
  src/octmem.c
  src/minrad.c
  src/ordie.c
  src/pdcode.c
  src/pd_isomorphisms.c
  src/pd_container.c
  src/pd_multidx.c
  src/pdcode_testcodes.c
  src/pd_dihedral.c
  src/pd_cyclic.c
  src/pd_orientation.c
  src/pd_perm.c
  src/plcToricSymplectic.c
  src/plcToricSymplecticTriangulation.c
  src/geyer_estimators.c
  src/polynomials.c
  src/pd_operations.c
  src/pd_crossingmoves.c
  src/pd_sortedbuf.c
  src/pd_tangles.c
  src/pd_splitdiagram.c
  src/pd_deletions.c
  src/kascentpermutation.c
  src/hypersimplex_sampling.c
  src/pd_pack.c
  src/pd_storage.c
  src/pd_invariants.c
  src/pd_simplify.c
)

add_library(plCurve STATIC ${PLCURVE_SOURCES})
add_library(plCurve::plCurve ALIAS plCurve)

target_include_directories(plCurve
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${GEN_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${JUDY_SRC_DIR}
)

target_compile_definitions(plCurve PRIVATE HAVE_CONFIG_H)
target_compile_options(plCurve PRIVATE -Wall)
target_link_libraries(plCurve PUBLIC Judy GSL::gsl GSL::gslcblas m)
if(OpenMP_C_FOUND)
  target_link_libraries(plCurve PUBLIC OpenMP::OpenMP_C)
  message(STATUS "plCurve: OpenMP enabled for octrope")
else()
  message(STATUS "plCurve: OpenMP not found; octrope stays single-threaded")
endif()

install(TARGETS plCurve Judy
  EXPORT plCurveTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES
  src/plCurve.h
  src/octrope.h
  src/nplCurve.h
  src/plcTopology.h
  src/plcRandomPolygon.h
  src/pd_invariants.h
  ${JUDY_SRC_DIR}/Judy.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT plCurveTargets
  FILE plCurveTargets.cmake
  NAMESPACE plCurve::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/plCurve
)

enable_testing()
add_executable(octrope_omp_correctness test/octrope_omp_correctness.c)
target_include_directories(octrope_omp_correctness PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src ${GEN_DIR})
target_compile_definitions(octrope_omp_correctness PRIVATE HAVE_CONFIG_H)
target_link_libraries(octrope_omp_correctness PRIVATE plCurve)
add_test(NAME octrope_omp_correctness COMMAND octrope_omp_correctness)
