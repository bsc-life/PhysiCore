# Core biofvm library without solver registering code
add_library(physicore.reactions-diffusion.biofvm_core OBJECT)
# biofvm library with solver registering code but without the code that
# explicitly registers each kernels. Suitable for linking against kernel test
# executables
add_library(physicore.reactions-diffusion.biofvm_no_kernels)
# main biofvm library containing all supported kernels
add_library(physicore.reactions-diffusion.biofvm)
add_library(physicore::reactions-diffusion::biofvm ALIAS
            physicore.reactions-diffusion.biofvm)
add_library(physicore.reactions-diffusion.biofvm_internal_iface INTERFACE)

target_sources(
  physicore.reactions-diffusion.biofvm_core
  PRIVATE src/microenvironment.cpp
          src/microenvironment_builder.cpp
          src/mesh.cpp
          src/vtk_serializer.cpp
          src/vtk_serializer_base.cpp
          src/config_reader.cpp
          src/vtk_agents_serializer.cpp)

target_sources(
  physicore.reactions-diffusion.biofvm_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(physicore.reactions-diffusion.biofvm_internal_iface
                      INTERFACE physicore::common)

find_package(VTK CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)

target_link_libraries(physicore.reactions-diffusion.biofvm_internal_iface
                      INTERFACE VTK::IOXML pugixml::pugixml)

target_link_libraries(
  physicore.reactions-diffusion.biofvm_core
  PUBLIC physicore.reactions-diffusion.biofvm_internal_iface)

target_link_libraries(physicore.reactions-diffusion.biofvm_no_kernels
                      PUBLIC physicore.reactions-diffusion.biofvm_core)

target_sources(physicore.reactions-diffusion.biofvm_no_kernels
               PRIVATE src/solver_registry_sole.cpp)

target_link_libraries(
  physicore.reactions-diffusion.biofvm
  PRIVATE $<TARGET_OBJECTS:physicore.reactions-diffusion.biofvm_core>
          $<BUILD_INTERFACE:physicore.reactions-diffusion.biofvm_internal_iface>
  PUBLIC physicore::common)

file(GLOB PUBLIC_INC "include/*.h")
target_sources(
  physicore.reactions-diffusion.biofvm
  PRIVATE src/solver_registry_attach.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

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

install(
  TARGETS physicore.reactions-diffusion.biofvm
  EXPORT MyLibTargets
  FILE_SET HEADERS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MyLibConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MyLibConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/MyLibConfigVersion.cmake"
        DESTINATION lib/cmake/MyLib)
install(
  EXPORT MyLibTargets
  NAMESPACE MyLib::
  FILE MyLibTargets.cmake
  DESTINATION lib/cmake/MyLib)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
