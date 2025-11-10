add_library(physicore.mechanics.micromechanics)
add_library(physicore::mechanics::micromechanics ALIAS physicore.mechanics.micromechanics)

target_sources(physicore.mechanics.micromechanics PRIVATE src/environment.cpp)

target_sources(
  physicore.mechanics.micromechanics
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/include/environment.h")

target_link_libraries(physicore.mechanics.micromechanics PUBLIC physicore::common)

if(PHYSICORE_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(physicore.mechanics.micromechanics.tests)

  file(GLOB TEST_SOURCES tests/*.cpp)
  target_sources(physicore.mechanics.micromechanics.tests PRIVATE ${TEST_SOURCES})

  target_link_libraries(
    physicore.mechanics.micromechanics.tests PRIVATE physicore::mechanics::micromechanics
                                                GTest::gtest GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(physicore.mechanics.micromechanics.tests)
endif()
