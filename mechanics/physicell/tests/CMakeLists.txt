find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(mechanics.physicell.tests ${TEST_SOURCES})

target_link_libraries(
  mechanics.physicell.tests
  PRIVATE mechanics.physicell mechanics.physicell_internal_iface GTest::gtest
          GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(mechanics.physicell.tests)
