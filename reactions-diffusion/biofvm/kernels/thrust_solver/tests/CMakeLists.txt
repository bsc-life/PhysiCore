find_package(GTest REQUIRED)

add_executable(physicore.reactions-diffusion.kernels.thrust_solver.tests)
enable_language(CUDA)

file(GLOB TEST_SOURCES *.cpp)
target_sources(physicore.reactions-diffusion.kernels.thrust_solver.tests
               PRIVATE ${TEST_SOURCES})

set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)
set_target_properties(physicore.reactions-diffusion.kernels.thrust_solver.tests
                      PROPERTIES CUDA_ARCHITECTURES native)

target_link_libraries(
  physicore.reactions-diffusion.kernels.thrust_solver.tests
  PRIVATE physicore::reactions-diffusion::kernels::thrust_solver_cuda
          physicore::reactions-diffusion::biofvm GTest::gtest GTest::gtest_main)

target_include_directories(
  physicore.reactions-diffusion.kernels.thrust_solver.tests
  PRIVATE ../../openmp_solver/noarr-structures/include)

include(GoogleTest)
gtest_discover_tests(physicore.reactions-diffusion.kernels.thrust_solver.tests)
