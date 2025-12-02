add_library(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  INTERFACE)
add_library(physicore.reactions-diffusion.biofvm.kernels.openmp_solver OBJECT)

if(${BUILD_SHARED_LIBS})
  set_target_properties(
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver
    PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

find_package(hwy CONFIG REQUIRED)
find_package(OpenMP 4)

target_sources(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver
  PRIVATE src/bulk_solver.cpp src/cell_solver.cpp src/diffusion_solver.cpp
          src/dirichlet_solver.cpp src/openmp_solver.cpp src/register_solver.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:physicore.reactions-diffusion.biofvm_iface>
            physicore::common $<BUILD_INTERFACE:hwy::hwy>)

find_path(NOARR_STRUCTURES_INCLUDE_DIRS "noarr/introspection.hpp")
target_include_directories(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${NOARR_STRUCTURES_INCLUDE_DIRS}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
    INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver
  PRIVATE
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface)

target_link_libraries(
  physicore.reactions-diffusion.biofvm
  PRIVATE physicore.reactions-diffusion.biofvm.kernels.openmp_solver)

install(
  TARGETS
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  EXPORT BioFVMTargets
  FILE_SET HEADERS)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
