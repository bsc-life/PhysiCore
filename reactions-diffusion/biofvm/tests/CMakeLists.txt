find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.reactions-diffusion.biofvm.tests ${TEST_SOURCES})

target_link_libraries(
  physicore.reactions-diffusion.biofvm.tests
  PRIVATE physicore.reactions-diffusion.biofvm_core
          physicore.reactions-diffusion.biofvm GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(physicore.reactions-diffusion.biofvm.tests)
