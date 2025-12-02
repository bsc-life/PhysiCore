find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.mechanics.physicell.tests ${TEST_SOURCES})

target_compile_options(physicore.mechanics.physicell.tests
                       PRIVATE -fno-sanitize=leak)
target_link_options(physicore.mechanics.physicell.tests PRIVATE
                    -fno-sanitize=leak)

target_link_libraries(
  physicore.mechanics.physicell.tests
  PRIVATE physicore.mechanics.physicell
          physicore.mechanics.physicell_internal_iface GTest::gtest
          GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(physicore.mechanics.physicell.tests TEST_EXECUTOR
                     ${CMAKE_COMMAND} -E env LSAN_OPTIONS=detect_leaks=0)
