# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.24)
project(
  PhysiCore
  LANGUAGES CXX
  VERSION 0.1.0)

add_library(libA)
add_library(libB OBJECT)
add_library(iface_obj3 INTERFACE)
target_link_libraries(iface_obj3 INTERFACE $<TARGET_OBJECTS:libB>)

target_sources(libA PRIVATE src/libA.cpp)
target_sources(libB PRIVATE src/libB.cpp)

target_link_libraries(libA PRIVATE iface_obj3)

install(TARGETS libA iface_obj3 EXPORT BioFVMTargets)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MyLibConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MyLibConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/MyLibConfigVersion.cmake"
        DESTINATION lib/cmake/MyLib)
install(
  EXPORT BioFVMTargets
  NAMESPACE MyLib::
  FILE BioFVMTargets.cmake
  DESTINATION lib/cmake/MyLib)
