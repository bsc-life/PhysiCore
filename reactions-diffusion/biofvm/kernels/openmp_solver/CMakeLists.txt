add_library(physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal
            OBJECT)
add_library(physicore.reactions-diffusion.biofvm.kernels.openmp_solver)

find_package(hwy CONFIG REQUIRED)
find_package(OpenMP 4)

target_sources(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal
  PRIVATE src/bulk_solver.cpp src/cell_solver.cpp src/diffusion_solver.cpp
          src/dirichlet_solver.cpp src/openmp_solver.cpp src/register_solver.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
         "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal
  PUBLIC physicore::common hwy::hwy)

find_path(NOARR_STRUCTURES_INCLUDE_DIRS "noarr/introspection.hpp")
target_include_directories(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal
  PUBLIC ${NOARR_STRUCTURES_INCLUDE_DIRS})

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal
    PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver
  PRIVATE physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal)

target_sources(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  physicore.reactions-diffusion.biofvm
  PRIVATE physicore.reactions-diffusion.biofvm.kernels.openmp_solver)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
