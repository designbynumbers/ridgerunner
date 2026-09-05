# Judy static library for MinGW/CMake (mirrors judy/src/*/Makefile.am).
# Common sources are compiled twice: once with JUDYL and once with JUDY1.

set(JUDY_COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/judy/src/JudyCommon)
set(JUDY_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/judy/src)
set(JUDY_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/judy-gen)
file(MAKE_DIRECTORY ${JUDY_GEN_DIR}/JudyL ${JUDY_GEN_DIR}/Judy1)

set(JUDY_INCLUDES
  ${JUDY_SRC_DIR}
  ${JUDY_COMMON_DIR}
  ${JUDY_SRC_DIR}/JudyL
  ${JUDY_SRC_DIR}/Judy1
  ${JUDY_GEN_DIR}
)

# --- JudyL / Judy1 from JudyCommon (unique object names via OBJECT libs) ---
set(JUDY_COMMON_BASE
  JudyCascade
  JudyCount
  JudyCreateBranch
  JudyDecascade
  JudyDel
  JudyFirst
  JudyFreeArray
  JudyGet
  JudyInsArray
  JudyIns
  JudyInsertBranch
  JudyMallocIF
  JudyMemActive
  JudyMemUsed
  JudyByCount
)

function(judy_objects prefix define out_var)
  set(_objs "")
  foreach(base ${JUDY_COMMON_BASE})
    set(_name ${prefix}_${base})
    add_library(${_name} OBJECT ${JUDY_COMMON_DIR}/${base}.c)
    target_compile_definitions(${_name} PRIVATE ${define})
    target_include_directories(${_name} PRIVATE ${JUDY_INCLUDES})
    if("${base}" STREQUAL "JudyByCount")
      target_compile_definitions(${_name} PRIVATE NOSMARTJBB NOSMARTJBU NOSMARTJLB)
    endif()
    list(APPEND _objs $<TARGET_OBJECTS:${_name}>)
  endforeach()

  # Prev/Next variants
  add_library(${prefix}_Next OBJECT ${JUDY_COMMON_DIR}/JudyPrevNext.c)
  target_compile_definitions(${prefix}_Next PRIVATE ${define} JUDYNEXT)
  target_include_directories(${prefix}_Next PRIVATE ${JUDY_INCLUDES})
  list(APPEND _objs $<TARGET_OBJECTS:${prefix}_Next>)

  add_library(${prefix}_Prev OBJECT ${JUDY_COMMON_DIR}/JudyPrevNext.c)
  target_compile_definitions(${prefix}_Prev PRIVATE ${define} JUDYPREV)
  target_include_directories(${prefix}_Prev PRIVATE ${JUDY_INCLUDES})
  list(APPEND _objs $<TARGET_OBJECTS:${prefix}_Prev>)

  add_library(${prefix}_NextEmpty OBJECT ${JUDY_COMMON_DIR}/JudyPrevNextEmpty.c)
  target_compile_definitions(${prefix}_NextEmpty PRIVATE ${define} JUDYNEXT)
  target_include_directories(${prefix}_NextEmpty PRIVATE ${JUDY_INCLUDES})
  list(APPEND _objs $<TARGET_OBJECTS:${prefix}_NextEmpty>)

  add_library(${prefix}_PrevEmpty OBJECT ${JUDY_COMMON_DIR}/JudyPrevNextEmpty.c)
  target_compile_definitions(${prefix}_PrevEmpty PRIVATE ${define} JUDYPREV)
  target_include_directories(${prefix}_PrevEmpty PRIVATE ${JUDY_INCLUDES})
  list(APPEND _objs $<TARGET_OBJECTS:${prefix}_PrevEmpty>)

  # Inline Get variant
  add_library(${prefix}_GetInline OBJECT ${JUDY_COMMON_DIR}/JudyGet.c)
  target_compile_definitions(${prefix}_GetInline PRIVATE ${define} JUDYGETINLINE)
  target_include_directories(${prefix}_GetInline PRIVATE ${JUDY_INCLUDES})
  list(APPEND _objs $<TARGET_OBJECTS:${prefix}_GetInline>)

  set(${out_var} "${_objs}" PARENT_SCOPE)
endfunction()

judy_objects(JudyL JUDYL JUDYL_OBJS)
judy_objects(Judy1 JUDY1 JUDY1_OBJS)

# --- Generated JudyLTables.c / Judy1Tables.c ---
add_executable(JudyLTablesGen ${JUDY_COMMON_DIR}/JudyTables.c)
target_compile_definitions(JudyLTablesGen PRIVATE JUDYL)
target_include_directories(JudyLTablesGen PRIVATE ${JUDY_INCLUDES})

add_executable(Judy1TablesGen ${JUDY_COMMON_DIR}/JudyTables.c)
target_compile_definitions(Judy1TablesGen PRIVATE JUDY1)
target_include_directories(Judy1TablesGen PRIVATE ${JUDY_INCLUDES})

add_custom_command(
  OUTPUT ${JUDY_GEN_DIR}/JudyL/JudyLTables.c
  COMMAND JudyLTablesGen
  WORKING_DIRECTORY ${JUDY_GEN_DIR}/JudyL
  DEPENDS JudyLTablesGen
  COMMENT "Generating JudyLTables.c"
)
add_custom_command(
  OUTPUT ${JUDY_GEN_DIR}/Judy1/Judy1Tables.c
  COMMAND Judy1TablesGen
  WORKING_DIRECTORY ${JUDY_GEN_DIR}/Judy1
  DEPENDS Judy1TablesGen
  COMMENT "Generating Judy1Tables.c"
)

add_library(JudyLTables OBJECT ${JUDY_GEN_DIR}/JudyL/JudyLTables.c)
target_compile_definitions(JudyLTables PRIVATE JUDYL)
target_include_directories(JudyLTables PRIVATE ${JUDY_INCLUDES})

add_library(Judy1Tables OBJECT ${JUDY_GEN_DIR}/Judy1/Judy1Tables.c)
target_compile_definitions(Judy1Tables PRIVATE JUDY1)
target_include_directories(Judy1Tables PRIVATE ${JUDY_INCLUDES})

add_library(JudyMalloc OBJECT ${JUDY_COMMON_DIR}/JudyMalloc.c)
target_include_directories(JudyMalloc PRIVATE ${JUDY_INCLUDES})

add_library(JudySL OBJECT ${JUDY_SRC_DIR}/JudySL/JudySL.c)
target_include_directories(JudySL PRIVATE ${JUDY_INCLUDES})

add_library(JudyHS OBJECT ${JUDY_SRC_DIR}/JudyHS/JudyHS.c)
target_include_directories(JudyHS PRIVATE ${JUDY_INCLUDES})

add_library(Judy STATIC
  ${JUDYL_OBJS}
  ${JUDY1_OBJS}
  $<TARGET_OBJECTS:JudyLTables>
  $<TARGET_OBJECTS:Judy1Tables>
  $<TARGET_OBJECTS:JudyMalloc>
  $<TARGET_OBJECTS:JudySL>
  $<TARGET_OBJECTS:JudyHS>
)
target_include_directories(Judy PUBLIC
  $<BUILD_INTERFACE:${JUDY_SRC_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
