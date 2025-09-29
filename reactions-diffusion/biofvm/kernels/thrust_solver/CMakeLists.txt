add_library(physicore.reactions-diffusion.kernels.thrust_solver)
add_library(physicore::reactions-diffusion::kernels::thrust_solver ALIAS
            physicore.reactions-diffusion.kernels.thrust_solver)

set(Thrust_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cccl/lib/cmake/thrust")

# add_subdirectory(cccl)

find_package(Thrust REQUIRED CONFIG)
# find_package(libcudacxx REQUIRED CONFIG)

thrust_create_target(ThrustTBB HOST CPP DEVICE TBB)
# thrust_create_target(ThrustCUDA DEVICE CUDA)

# enable_language(CUDA)

# find_package(CUDAToolkit 12)

set_target_properties(physicore.reactions-diffusion.kernels.thrust_solver
                      PROPERTIES CUDA_ARCHITECTURES native)

target_include_directories(physicore.reactions-diffusion.kernels.thrust_solver
                           PRIVATE ../openmp_solver/noarr-structures/include)

target_sources(
  physicore.reactions-diffusion.kernels.thrust_solver
  PRIVATE src/bulk_solver.cpp src/dirichlet_solver.cpp src/cell_solver.cpp
          src/diffusion_solver.cpp src/thrust_solver.cpp)

target_sources(
  physicore.reactions-diffusion.kernels.thrust_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

find_package(hwy CONFIG REQUIRED)

target_link_libraries(physicore.reactions-diffusion.kernels.thrust_solver
                      PUBLIC physicore::common ThrustTBB hwy::hwy)

target_include_directories(physicore.reactions-diffusion.kernels.thrust_solver
                           PRIVATE noarr-structures/include)

# if(PHYSICORE_BUILD_EXAMPLES) add_subdirectory(examples) endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
