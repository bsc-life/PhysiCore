find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(physicore.mechanics.micromechanics.kernels.openmp_solver.tests
               ${TEST_SOURCES})

target_link_libraries(
  physicore.mechanics.micromechanics.kernels.openmp_solver.tests
  PRIVATE
    physicore.mechanics.micromechanics.kernels.openmp_solver
    physicore.mechanics.micromechanics.kernels.openmp_solver_internal_iface
    physicore.mechanics.micromechanics
    physicore::common
    GTest::gtest
    GTest::gtest_main)

target_include_directories(
  physicore.mechanics.micromechanics.kernels.openmp_solver.tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src
          ${CMAKE_CURRENT_SOURCE_DIR}/../../../include)

include(GoogleTest)
gtest_discover_tests(
  physicore.mechanics.micromechanics.kernels.openmp_solver.tests)
