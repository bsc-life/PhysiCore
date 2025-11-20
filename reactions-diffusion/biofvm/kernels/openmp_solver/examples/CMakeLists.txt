set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver.${example}
    ${example}.cpp)
  target_link_libraries(
    physicore.reactions-diffusion.biofvm.kernels.openmp_solver.${example}
    PRIVATE
      physicore.reactions-diffusion.biofvm.kernels.openmp_solver
      physicore.reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
      physicore.reactions-diffusion.biofvm_no_kernels)
endforeach()
