find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.reactions-diffusion.biofvm.kernels.openmp_solver.tests
               ${TEST_SOURCES})

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver.tests
  PRIVATE
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
    physicore.reactions-diffusion.biofvm_no_kernels
    GTest::gtest
    GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(
  physicore.reactions-diffusion.biofvm.kernels.openmp_solver.tests)
