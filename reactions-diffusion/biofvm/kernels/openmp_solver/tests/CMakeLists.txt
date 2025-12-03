find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(reactions-diffusion.biofvm.kernels.openmp_solver.tests
               ${TEST_SOURCES})

target_link_libraries(
  reactions-diffusion.biofvm.kernels.openmp_solver.tests
  PRIVATE reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
          reactions-diffusion.biofvm GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(reactions-diffusion.biofvm.kernels.openmp_solver.tests)
