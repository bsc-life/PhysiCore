find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.mechanics.micromechanics.tests ${TEST_SOURCES})

target_link_libraries(
  physicore.mechanics.micromechanics.tests
  PRIVATE
    physicore.mechanics.micromechanics
    physicore.mechanics.micromechanics_internal_iface
    physicore.mechanics.micromechanics.kernels.openmp_solver
    physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
    GTest::gtest
    GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(physicore.mechanics.micromechanics.tests)
