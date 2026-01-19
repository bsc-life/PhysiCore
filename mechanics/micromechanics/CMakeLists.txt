add_library(physicore.mechanics.micromechanics)
add_library(physicore::mechanics::micromechanics ALIAS
            physicore.mechanics.micromechanics)
add_library(physicore.mechanics.micromechanics_internal_iface INTERFACE)
add_library(physicore.mechanics.micromechanics_iface INTERFACE)

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  physicore.mechanics.micromechanics
  PRIVATE src/environment.cpp src/cell_data.cpp src/simulation_parameters.cpp
    src/uniform_grid_spatial_index.cpp src/solver_registry_sole.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  physicore.mechanics.micromechanics_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/micromechanics"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_sources(
  physicore.mechanics.micromechanics_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.mechanics.micromechanics PUBLIC physicore::common)

target_include_directories(
  physicore.mechanics.micromechanics
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
