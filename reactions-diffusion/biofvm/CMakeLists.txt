add_library(reactions-diffusion.biofvm)
add_library(reactions-diffusion::biofvm ALIAS reactions-diffusion.biofvm)
add_library(reactions-diffusion.biofvm_internal_iface INTERFACE)
add_library(reactions-diffusion.biofvm_iface INTERFACE)

include(GenerateExportHeader)
generate_export_header(
  reactions-diffusion.biofvm BASE_NAME biofvm EXPORT_FILE_NAME
  ${CMAKE_CURRENT_BINARY_DIR}/include/biofvm/biofvm_export.h)

if(${BUILD_SHARED_LIBS})
  set_target_properties(
    reactions-diffusion.biofvm
    PROPERTIES CXX_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON
               POSITION_INDEPENDENT_CODE ON
               INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

file(GLOB PUBLIC_INC "include/biofvm/*.h")
target_sources(
  reactions-diffusion.biofvm
  PRIVATE src/microenvironment.cpp
          src/microenvironment_builder.cpp
          src/mesh.cpp
          src/vtk_serializer.cpp
          src/vtk_serializer_base.cpp
          src/config_reader.cpp
          src/vtk_agents_serializer.cpp
          src/solver_registry.cpp
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         "${CMAKE_CURRENT_SOURCE_DIR}/include"
         "${CMAKE_CURRENT_BINARY_DIR}/include"
         FILES
         ${PUBLIC_INC}
         "${CMAKE_CURRENT_BINARY_DIR}/include/biofvm/biofvm_export.h")

target_include_directories(
  reactions-diffusion.biofvm_internal_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/biofvm>)

target_sources(
  reactions-diffusion.biofvm_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_BINARY_DIR}/include")

target_compile_definitions(reactions-diffusion.biofvm_iface INTERFACE BIOFVM_STATIC_DEFINE)

target_link_libraries(reactions-diffusion.biofvm_internal_iface
                      INTERFACE common)

find_package(VTK CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)

target_link_libraries(reactions-diffusion.biofvm_internal_iface
                      INTERFACE VTK::IOXML pugixml::pugixml)

target_link_libraries(
  reactions-diffusion.biofvm
  PRIVATE reactions-diffusion.biofvm_internal_iface
  PUBLIC common)

add_subdirectory(kernels/thrust_solver)

get_target_property(HAS_CUDA reactions-diffusion.biofvm INTERFACE_HAS_CUDA)
if(HAS_CUDA)
  set(BIOFVM_HAS_CUDA ON)
  enable_language(CUDA)
else()
  set(BIOFVM_HAS_CUDA OFF)
endif()

add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

install(
  TARGETS reactions-diffusion.biofvm reactions-diffusion.biofvm_internal_iface
  EXPORT BioFVMTargets
  FILE_SET HEADERS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  cmake/BioFVMConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfig.cmake"
  INSTALL_DESTINATION lib/cmake/BioFVM)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/BioFVMConfigVersion.cmake"
        DESTINATION lib/cmake/BioFVM)

install(
  EXPORT BioFVMTargets
  NAMESPACE PhysiCore::
  FILE BioFVMTargets.cmake
  DESTINATION lib/cmake/BioFVM)
