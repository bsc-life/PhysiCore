add_library(physicore.common INTERFACE)
add_library(physicore::common ALIAS physicore.common)

target_sources(
  physicore.common
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         "${CMAKE_CURRENT_SOURCE_DIR}/include"
         FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/include/base_agent_container.h"
         "${CMAKE_CURRENT_SOURCE_DIR}/include/base_agent_data.h"
         "${CMAKE_CURRENT_SOURCE_DIR}/include/base_agent.h"
         "${CMAKE_CURRENT_SOURCE_DIR}/include/types.h")

if(PHYSICORE_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(physicore.common.tests)

  file(GLOB TEST_SOURCES tests/*.cpp)
  target_sources(physicore.common.tests PRIVATE ${TEST_SOURCES})

  target_link_libraries(
    physicore.common.tests PRIVATE physicore::common GTest::gtest
                                   GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(physicore.common.tests)
endif()
