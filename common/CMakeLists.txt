add_library(physicore.common INTERFACE)
add_library(physicore::common ALIAS physicore.common)

file(GLOB PUBLIC_INC "include/*.h")
target_sources(
  physicore.common
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

install(
  TARGETS physicore.common
  EXPORT MyLibTargets
  FILE_SET HEADERS)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
