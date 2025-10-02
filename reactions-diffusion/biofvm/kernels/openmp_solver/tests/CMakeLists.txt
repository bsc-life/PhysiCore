find_package(GTest REQUIRED)

add_executable(physicore.reactions-diffusion.kernels.openmp_solver.tests)

file(GLOB TEST_SOURCES *.cpp)
target_sources(physicore.reactions-diffusion.kernels.openmp_solver.tests
               PRIVATE ${TEST_SOURCES})

target_link_libraries(
  physicore.reactions-diffusion.kernels.openmp_solver.tests
  PRIVATE physicore::reactions-diffusion::kernels::openmp_solver
          physicore::reactions-diffusion::biofvm GTest::gtest GTest::gtest_main)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    physicore.reactions-diffusion.kernels.openmp_solver.tests
    PRIVATE OpenMP::OpenMP_CXX)
endif()

include(GoogleTest)
gtest_discover_tests(physicore.reactions-diffusion.kernels.openmp_solver.tests)
