add_library(physicore.common INTERFACE)
add_library(physicore::common ALIAS physicore.common)

target_sources(physicore.common PUBLIC FILE_SET HEADERS BASE_DIRS
                                       "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
