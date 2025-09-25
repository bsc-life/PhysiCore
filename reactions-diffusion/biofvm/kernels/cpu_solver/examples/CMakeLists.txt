set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})
  add_executable(physicore.reactions-diffusion.kernels.cpu_solver.${example})
  target_sources(physicore.reactions-diffusion.kernels.cpu_solver.${example}
                 PRIVATE ${example}.cpp)
  target_link_libraries(
    physicore.reactions-diffusion.kernels.cpu_solver.${example}
    PRIVATE physicore::reactions-diffusion::kernels::cpu_solver
            physicore::reactions-diffusion::biofvm hwy::hwy)

  if(OpenMP_CXX_FOUND)
    target_link_libraries(
      physicore.reactions-diffusion.kernels.cpu_solver.${example}
      PRIVATE OpenMP::OpenMP_CXX)
  endif()

  target_include_directories(
    physicore.reactions-diffusion.kernels.cpu_solver.${example}
    PRIVATE ../noarr-structures/include)
endforeach()
