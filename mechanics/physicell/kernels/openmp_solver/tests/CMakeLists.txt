find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(mechanics.physicell.kernels.openmp_solver.tests ${TEST_SOURCES})

target_link_libraries(
  mechanics.physicell.kernels.openmp_solver.tests
  PRIVATE mechanics.physicell.kernels.openmp_solver_internal_iface
          mechanics.physicell GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(mechanics.physicell.kernels.openmp_solver.tests)
