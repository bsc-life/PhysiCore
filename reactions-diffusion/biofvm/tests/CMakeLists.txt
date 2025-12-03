find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(reactions-diffusion.biofvm.tests ${TEST_SOURCES})

target_link_libraries(
  reactions-diffusion.biofvm.tests
  PRIVATE reactions-diffusion.biofvm reactions-diffusion.biofvm_internal_iface
          GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(reactions-diffusion.biofvm.tests)
