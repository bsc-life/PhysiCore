add_library(physicore.reactions-diffusion.kernels.thrust_solver)
add_library(physicore::reactions-diffusion::kernels::thrust_solver ALIAS
            physicore.reactions-diffusion.kernels.thrust_solver)

set(Thrust_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cccl/lib/cmake/thrust")

find_package(Thrust REQUIRED CONFIG)

thrust_create_target(ThrustTBB HOST TBB DEVICE TBB)
# thrust_create_target(ThrustCUDA DEVICE CUDA)

target_include_directories(physicore.reactions-diffusion.kernels.thrust_solver
                           PRIVATE ../openmp_solver/noarr-structures/include)

target_sources(
  physicore.reactions-diffusion.kernels.thrust_solver
  PRIVATE src/bulk_solver.cpp src/dirichlet_solver.cpp src/cell_solver.cpp
          src/diffusion_solver.cpp src/thrust_solver.cpp)

target_sources(
  physicore.reactions-diffusion.kernels.thrust_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.reactions-diffusion.kernels.thrust_solver
                      PUBLIC physicore::common ThrustTBB)

find_package(CUDAToolkit 12)

# if(${CUDAToolkit_FOUND}) enable_language(CUDA)

# add_library(physicore.reactions-diffusion.kernels.thrust_solver_cuda)
# add_library(physicore::reactions-diffusion::kernels::thrust_solver_cuda ALIAS
# physicore.reactions-diffusion.kernels.thrust_solver_cuda)

# target_sources( physicore.reactions-diffusion.kernels.thrust_solver_cuda
# PRIVATE src/bulk_solver.cpp src/dirichlet_solver.cpp src/cell_solver.cpp
# src/diffusion_solver.cpp src/thrust_solver.cpp)

# #
# set_target_properties(physicore.reactions-diffusion.kernels.thrust_solver_cuda
# PROPERTIES CUDA_ARCHITECTURES native)

# # set_target_properties(physicore.reactions-diffusion.kernels.thrust_solver #
# PROPERTIES CUDA_ARCHITECTURES native) Get all source files from a target
# get_target_property( source_files
# physicore.reactions-diffusion.kernels.thrust_solver_cuda SOURCES)

# thrust_create_target(ThrustCUDA HOST TBB DEVICE CUDA)

# target_include_directories(
# physicore.reactions-diffusion.kernels.thrust_solver_cuda PRIVATE
# ../openmp_solver/noarr-structures/include)

# # Mark them all as CUDA set_source_files_properties(${source_files} PROPERTIES
# LANGUAGE CUDA)

# target_sources(physicore.reactions-diffusion.kernels.thrust_solver_cuda
# PRIVATE ${source_files})

# target_compile_options(
# physicore.reactions-diffusion.kernels.thrust_solver_cuda PUBLIC
# --extended-lambda --expt-relaxed-constexpr)

# target_sources( physicore.reactions-diffusion.kernels.thrust_solver_cuda
# PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

# target_link_libraries(physicore.reactions-diffusion.kernels.thrust_solver_cuda
# PUBLIC physicore::common ThrustCUDA)

# endif()

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
