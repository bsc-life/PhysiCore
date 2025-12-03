add_library(mechanics.micromechanics STATIC)
add_library(mechanics::micromechanics ALIAS mechanics.micromechanics)
add_library(mechanics.micromechanics_internal_iface INTERFACE)

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  mechanics.micromechanics
  PRIVATE src/environment.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  mechanics.micromechanics_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/micromechanics"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(mechanics.micromechanics_internal_iface INTERFACE common)

target_link_libraries(mechanics.micromechanics
                      PRIVATE mechanics.micromechanics_internal_iface)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
