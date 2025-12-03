find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(common.tests ${TEST_SOURCES})

target_link_libraries(common.tests PRIVATE common GTest::gtest
                                           GTest::gtest_main)

target_include_directories(common.tests PRIVATE "../include/common")

include(GoogleTest)
gtest_discover_tests(common.tests)
