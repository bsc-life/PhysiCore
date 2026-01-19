add_library(physicore.mechanics.micromechanics)
add_library(physicore::mechanics::micromechanics ALIAS
            physicore.mechanics.micromechanics)
add_library(physicore.mechanics.micromechanics_internal_iface INTERFACE)
add_library(physicore.mechanics.micromechanics_iface INTERFACE)

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  physicore.mechanics.micromechanics
  PRIVATE src/environment.cpp src/cell_data.cpp
          src/uniform_grid_spatial_index.cpp src/vtk_mechanics_serializer.cpp
          src/solver_registry_sole.cpp
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

target_link_libraries(physicore.mechanics.micromechanics_internal_iface
                      INTERFACE physicore::common)

find_package(VTK CONFIG REQUIRED)

target_link_libraries(physicore.mechanics.micromechanics_internal_iface
                      INTERFACE physicore::common::vtk_serializer VTK::IOXML)

target_link_libraries(
  physicore.mechanics.micromechanics
  PRIVATE physicore.mechanics.micromechanics_internal_iface
  PUBLIC physicore::common)

target_include_directories(
  physicore.mechanics.micromechanics
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Add OpenMP kernel
add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
