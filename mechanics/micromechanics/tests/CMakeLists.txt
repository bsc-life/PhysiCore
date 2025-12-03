find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(mechanics.micromechanics.tests ${TEST_SOURCES})

target_link_libraries(
  mechanics.micromechanics.tests
  PRIVATE mechanics.micromechanics mechanics.micromechanics_internal_iface
          GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(mechanics.micromechanics.tests)
