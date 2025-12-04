# Header-only common library (no VTK dependency)
add_library(physicore.common INTERFACE)
add_library(physicore::common ALIAS physicore.common)

file(GLOB PUBLIC_INC "include/common/*.h")
# Exclude vtk_serializer_base.h from header-only library
list(FILTER PUBLIC_INC EXCLUDE REGEX "vtk_serializer_base\\.h$")

target_sources(
  physicore.common
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

# VTK-dependent serializer base library (optional)
find_package(VTK CONFIG QUIET)
if(VTK_FOUND)
  add_library(physicore.common.vtk_serializer)
  add_library(physicore::common::vtk_serializer ALIAS
              physicore.common.vtk_serializer)

  target_sources(
    physicore.common.vtk_serializer
    PRIVATE src/vtk_serializer_base.cpp
    PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
           FILES include/common/vtk_serializer_base.h)

  target_link_libraries(physicore.common.vtk_serializer PUBLIC physicore::common
                                                               VTK::IOXML)

  install(
    TARGETS physicore.common.vtk_serializer
    EXPORT BioFVMTargets
    FILE_SET HEADERS)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(
  TARGETS physicore.common
  EXPORT BioFVMTargets
  FILE_SET HEADERS)
