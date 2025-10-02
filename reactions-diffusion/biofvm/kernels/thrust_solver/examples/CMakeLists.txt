set(ALL_EXAMPLES benchmark)

foreach(example ${ALL_EXAMPLES})

  # set_source_files_properties(${example} PROPERTIES LANGUAGE CUDA)
  # set_target_properties(
  # physicore.reactions-diffusion.kernels.thrust_solver.${example} PROPERTIES
  # CUDA_ARCHITECTURES native)

  add_executable(physicore.reactions-diffusion.kernels.thrust_solver.${example})
  target_sources(physicore.reactions-diffusion.kernels.thrust_solver.${example}
                 PRIVATE ${example}.cpp)
  target_link_libraries(
    physicore.reactions-diffusion.kernels.thrust_solver.${example}
    PRIVATE physicore::reactions-diffusion::kernels::thrust_solver
            physicore::reactions-diffusion::biofvm)

  target_include_directories(
    physicore.reactions-diffusion.kernels.thrust_solver.${example}
    PRIVATE ../../openmp_solver/noarr-structures/include)
endforeach()
