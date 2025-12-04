add_library(
  physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
  INTERFACE)
add_library(physicore.mechanics.micromechanics.kernels.openmp_solver)

find_package(OpenMP 4)

target_sources(
  physicore.mechanics.micromechanics.kernels.openmp_solver
  PRIVATE src/openmp_solver.cpp
          src/neighbor_solver.cpp
          src/force_solver.cpp
          src/position_solver.cpp
          src/motility_solver.cpp
          src/basement_membrane_solver.cpp
          src/spring_solver.cpp
          src/register_solver.cpp
          src/potentials/standard_potential.cpp
          src/potentials/morse_potential.cpp
          src/potentials/kelvin_voigt_potential.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_sources(
  physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(
  physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
  INTERFACE physicore.mechanics.micromechanics_internal_iface
            physicore.mechanics.micromechanics_iface physicore::common)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
    INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  physicore.mechanics.micromechanics.kernels.openmp_solver
  PRIVATE
    physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface)

target_link_libraries(
  physicore.mechanics.micromechanics
  PRIVATE
    $<BUILD_INTERFACE:physicore.mechanics.micromechanics.kernels.openmp_solver>)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
