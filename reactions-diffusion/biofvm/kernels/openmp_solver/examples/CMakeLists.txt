set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(physicore.reactions-diffusion.kernels.openmp_solver.${example})
  target_sources(physicore.reactions-diffusion.kernels.openmp_solver.${example}
                 PRIVATE ${example}.cpp)
  target_link_libraries(
    physicore.reactions-diffusion.kernels.openmp_solver.${example}
    PRIVATE physicore::reactions-diffusion::kernels::openmp_solver
            physicore::reactions-diffusion::biofvm_no_kernels hwy::hwy)

  if(OpenMP_CXX_FOUND)
    target_link_libraries(
      physicore.reactions-diffusion.kernels.openmp_solver.${example}
      PRIVATE OpenMP::OpenMP_CXX)
  endif()
endforeach()
