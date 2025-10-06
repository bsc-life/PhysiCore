set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(
    physicore.reactions-diffusion.kernels.tbb_thrust_solver.${example})
  target_sources(
    physicore.reactions-diffusion.kernels.tbb_thrust_solver.${example}
    PRIVATE ${example}.cpp)
  target_link_libraries(
    physicore.reactions-diffusion.kernels.tbb_thrust_solver.${example}
    PRIVATE physicore::reactions-diffusion::kernels::tbb_thrust_solver
            ThrustTBB physicore::reactions-diffusion::biofvm_no_kernels)
endforeach()

if(${CUDAToolkit_FOUND})
  foreach(example ${ALL_EXAMPLES})
    set(cuda_example "${CMAKE_CURRENT_BINARY_DIR}/cuda_sources/${example}.cu")

    add_custom_command(
      OUTPUT ${cuda_example}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp ${cuda_example}
      DEPENDS ${example}.cpp)

    add_executable(
      physicore.reactions-diffusion.kernels.cuda_thrust_solver.${example})

    target_sources(
      physicore.reactions-diffusion.kernels.cuda_thrust_solver.${example}
      PRIVATE ${cuda_example})

    set_target_properties(
      physicore.reactions-diffusion.kernels.cuda_thrust_solver.${example}
      PROPERTIES CUDA_ARCHITECTURES native CUDA_SEPARABLE_COMPILATION ON)

    target_link_libraries(
      physicore.reactions-diffusion.kernels.cuda_thrust_solver.${example}
      PRIVATE physicore::reactions-diffusion::kernels::cuda_thrust_solver
              ThrustCUDA physicore::reactions-diffusion::biofvm_no_kernels)
  endforeach()
endif()
