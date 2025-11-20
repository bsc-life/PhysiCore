find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.common.tests ${TEST_SOURCES})

target_link_libraries(physicore.common.tests
                      PRIVATE physicore::common GTest::gtest GTest::gtest_main)

target_include_directories(physicore.common.tests PRIVATE "../include/common")

include(GoogleTest)
gtest_discover_tests(physicore.common.tests)
