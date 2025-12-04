add_library(physicore.reactions-diffusion.biofvm)
add_library(physicore::reactions-diffusion::biofvm ALIAS
            physicore.reactions-diffusion.biofvm)
add_library(physicore.reactions-diffusion.biofvm_internal_iface INTERFACE)
add_library(physicore.reactions-diffusion.biofvm_iface INTERFACE)

file(GLOB PUBLIC_INC "include/biofvm/*.h")
target_sources(
  physicore.reactions-diffusion.biofvm
  PRIVATE src/microenvironment.cpp
          src/microenvironment_builder.cpp
          src/mesh.cpp
          src/vtk_serializer.cpp
          src/config_reader.cpp
          src/vtk_agents_serializer.cpp
          src/solver_registry.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  physicore.reactions-diffusion.biofvm_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/biofvm"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_sources(
  physicore.reactions-diffusion.biofvm_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.reactions-diffusion.biofvm_internal_iface
                      INTERFACE physicore::common)

find_package(VTK CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)

target_link_libraries(
  physicore.reactions-diffusion.biofvm_internal_iface
  INTERFACE physicore::common::vtk_serializer VTK::IOXML pugixml::pugixml)

target_link_libraries(
  physicore.reactions-diffusion.biofvm
  PRIVATE $<BUILD_INTERFACE:physicore.reactions-diffusion.biofvm_internal_iface>
  PUBLIC physicore::common)

# Thrust does not support AppleClang
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  add_subdirectory(kernels/thrust_solver)
endif()

get_target_property(HAS_CUDA physicore.reactions-diffusion.biofvm
                    INTERFACE_HAS_CUDA)
if(HAS_CUDA)
  enable_language(CUDA)
endif()

add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

install(
  TARGETS physicore.reactions-diffusion.biofvm
  EXPORT BioFVMTargets
  FILE_SET HEADERS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BioFVMConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfigVersion.cmake"
        DESTINATION lib/cmake/BioFVM)

install(
  EXPORT BioFVMTargets
  NAMESPACE BioFVM::
  FILE BioFVMTargets.cmake
  DESTINATION lib/cmake/BioFVM)
