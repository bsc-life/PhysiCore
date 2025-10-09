# Core biofvm library without solver registering code
add_library(physicore.reactions-diffusion.biofvm_core OBJECT)
add_library(physicore::reactions-diffusion::biofvm_core ALIAS
            physicore.reactions-diffusion.biofvm_core)

target_sources(
  physicore.reactions-diffusion.biofvm_core
  PRIVATE src/microenvironment.cpp src/microenvironment_builder.cpp
          src/mesh.cpp)

target_sources(
  physicore.reactions-diffusion.biofvm_core
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

find_package(VTK CONFIG REQUIRED)

target_link_libraries(physicore.reactions-diffusion.biofvm_core
                      PUBLIC physicore::common VTK::IOXML)

# biofvm library with solver registering code but without the code that
# explicitly registers each kernels suitable for linking against kernel test
# executables
add_library(physicore.reactions-diffusion.biofvm_no_kernels)
add_library(physicore::reactions-diffusion::biofvm_no_kernels ALIAS
            physicore.reactions-diffusion.biofvm_no_kernels)

target_sources(physicore.reactions-diffusion.biofvm_no_kernels
               PRIVATE src/solver_registry_sole.cpp)

target_link_libraries(physicore.reactions-diffusion.biofvm_no_kernels
                      PUBLIC physicore::reactions-diffusion::biofvm_core)

# main biofvm library containing all supported kernels
add_library(physicore.reactions-diffusion.biofvm)
add_library(physicore::reactions-diffusion::biofvm ALIAS
            physicore.reactions-diffusion.biofvm)

target_sources(physicore.reactions-diffusion.biofvm
               PRIVATE src/solver_registry_attach.cpp)

target_link_libraries(physicore.reactions-diffusion.biofvm
                      PUBLIC physicore::reactions-diffusion::biofvm_core)

add_subdirectory(kernels/openmp_solver)

# Thrust does not support AppleClang
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  add_subdirectory(kernels/thrust_solver)
endif()

get_target_property(HAS_CUDA physicore.reactions-diffusion.biofvm
                    INTERFACE_HAS_CUDA)
if(HAS_CUDA)
  enable_language(CUDA)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
