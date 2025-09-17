add_library(physicore.mechanics.physicell)
add_library(physicore::mechanics::physicell ALIAS physicore.mechanics.physicell)

target_sources(physicore.mechanics.physicell PRIVATE src/environment.cpp)

target_sources(
  physicore.mechanics.physicell
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/include/environment.h")

target_link_libraries(physicore.mechanics.physicell PUBLIC physicore::common)

if(PHYSICORE_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(physicore.mechanics.physicell.tests)

  file(GLOB TEST_SOURCES tests/*.cpp)
  target_sources(physicore.mechanics.physicell.tests PRIVATE ${TEST_SOURCES})

  target_link_libraries(
    physicore.mechanics.physicell.tests PRIVATE physicore::mechanics::physicell
                                                GTest::gtest GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(physicore.mechanics.physicell.tests)
endif()
