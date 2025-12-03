set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(reactions-diffusion.biofvm.kernels.openmp_solver.${example}
                 ${example}.cpp)
  target_link_libraries(
    reactions-diffusion.biofvm.kernels.openmp_solver.${example}
    PRIVATE reactions-diffusion.biofvm.kernels.openmp_solver_internal_iface
            reactions-diffusion.biofvm)
endforeach()
