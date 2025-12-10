add_library(reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
            INTERFACE)
add_library(reactions-diffusion.biofvm.kernels.openmp_solver OBJECT)

if(${BUILD_SHARED_LIBS})
  set_target_properties(reactions-diffusion.biofvm.kernels.openmp_solver
                        PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

find_package(hwy CONFIG REQUIRED)
find_package(OpenMP 4)

target_sources(
  reactions-diffusion.biofvm.kernels.openmp_solver
  PRIVATE src/bulk_solver.cpp src/cell_solver.cpp src/diffusion_solver.cpp
          src/dirichlet_solver.cpp src/openmp_solver.cpp src/register_solver.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:reactions-diffusion.biofvm_iface> common hwy::hwy)

find_path(NOARR_STRUCTURES_INCLUDE_DIRS "noarr/introspection.hpp")
target_include_directories(
  reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${NOARR_STRUCTURES_INCLUDE_DIRS}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
    INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  reactions-diffusion.biofvm.kernels.openmp_solver
  PRIVATE reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface)

target_link_libraries(reactions-diffusion.biofvm
                      PRIVATE reactions-diffusion.biofvm.kernels.openmp_solver)

install(
  TARGETS reactions-diffusion.biofvm.kernels.openmp_solver
          reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
  EXPORT BioFVMTargets
  FILE_SET HEADERS)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
