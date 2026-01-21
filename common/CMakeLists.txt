add_library(common)

file(GLOB PUBLIC_INC "include/common/*.h")
file(GLOB SOURCES "src/*.cpp")

target_sources(
  common
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC}
  PRIVATE ${SOURCES})

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(
  TARGETS common
  EXPORT BioFVMTargets
  FILE_SET HEADERS)
