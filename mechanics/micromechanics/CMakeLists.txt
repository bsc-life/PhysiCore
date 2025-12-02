add_library(physicore.mechanics.micromechanics STATIC)
add_library(physicore::mechanics::micromechanics ALIAS
            physicore.mechanics.micromechanics)
add_library(physicore.mechanics.micromechanics_internal_iface INTERFACE)

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  physicore.mechanics.micromechanics
  PRIVATE src/environment.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  physicore.mechanics.micromechanics_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/micromechanics"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(physicore.mechanics.micromechanics_internal_iface
                      INTERFACE physicore::common)

target_link_libraries(physicore.mechanics.micromechanics
                      PRIVATE physicore.mechanics.micromechanics_internal_iface)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
