find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(
  physicore.reactions-diffusion.biofvm.kernels.tbb_thrust_solver.tests
  ${TEST_SOURCES})

target_link_libraries(
  physicore.reactions-diffusion.biofvm.kernels.tbb_thrust_solver.tests
  PRIVATE
    physicore.reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal
    physicore.reactions-diffusion.biofvm_no_kernels GTest::gtest
    GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(
  physicore.reactions-diffusion.biofvm.kernels.tbb_thrust_solver.tests)

if(${CUDAToolkit_FOUND})
  set(CUDA_TEST_SOURCES "")

  # Copy source into build dir with .cu extension
  foreach(src ${TEST_SOURCES})
    file(RELATIVE_PATH relpath "${CMAKE_CURRENT_SOURCE_DIR}" "${src}")
    get_filename_component(fname ${relpath} NAME_WE)

    set(dst "${CMAKE_CURRENT_BINARY_DIR}/cuda_sources/${fname}.cu")

    add_custom_command(
      OUTPUT ${dst}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${CMAKE_CURRENT_SOURCE_DIR}/${relpath} ${dst}
      DEPENDS ${src})

    list(APPEND CUDA_TEST_SOURCES ${dst})
  endforeach()

  add_executable(
    physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver.tests)

  target_sources(
    physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver.tests
    PRIVATE ${CUDA_TEST_SOURCES})

  set_target_properties(
    physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver.tests
    PROPERTIES CUDA_ARCHITECTURES native CUDA_SEPARABLE_COMPILATION ON)

  target_link_libraries(
    physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver.tests
    PRIVATE
      physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal
      physicore.reactions-diffusion.biofvm_no_kernels GTest::gtest
      GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(
    physicore.reactions-diffusion.biofvm.kernels.cuda_thrust_solver.tests)
endif()
