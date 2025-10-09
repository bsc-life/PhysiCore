find_package(GTest REQUIRED)
add_executable(physicore.reactions-diffusion.biofvm.tests)

file(GLOB TEST_SOURCES *.cpp)
target_sources(physicore.reactions-diffusion.biofvm.tests
               PRIVATE ${TEST_SOURCES})

target_link_libraries(
  physicore.reactions-diffusion.biofvm.tests
  PRIVATE physicore::reactions-diffusion::biofvm
          physicore::reactions-diffusion::kernels::openmp_solver GTest::gtest
          GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(physicore.reactions-diffusion.biofvm.tests)
