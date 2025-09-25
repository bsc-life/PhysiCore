add_library(physicore.reactions-diffusion.kernels.openmp_solver)
add_library(physicore::reactions-diffusion::kernels::openmp_solver ALIAS
            physicore.reactions-diffusion.kernels.openmp_solver)

find_package(hwy CONFIG REQUIRED)

find_package(OpenMP 4)

target_sources(
  physicore.reactions-diffusion.kernels.openmp_solver
  PRIVATE src/bulk_solver.cpp src/cell_solver.cpp src/diffusion_solver.cpp
          src/dirichlet_solver.cpp src/openmp_solver.cpp
          src/register_solver.cpp)

target_sources(
  physicore.reactions-diffusion.kernels.openmp_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.reactions-diffusion.kernels.openmp_solver
                      PUBLIC physicore::common)

target_link_libraries(physicore.reactions-diffusion.kernels.openmp_solver
                      PRIVATE hwy::hwy)

if(OpenMP_CXX_FOUND)
  target_link_libraries(physicore.reactions-diffusion.kernels.openmp_solver
                        PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(physicore.reactions-diffusion.kernels.openmp_solver
                           PRIVATE noarr-structures/include)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
