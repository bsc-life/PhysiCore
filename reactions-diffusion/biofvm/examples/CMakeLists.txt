set(ALL_EXAMPLES diffuse)

foreach(example ${ALL_EXAMPLES})
  add_executable(physicore.reactions-diffusion.biofvm.${example})
  target_sources(physicore.reactions-diffusion.biofvm.${example}
                 PRIVATE ${example}.cpp)

  target_link_libraries(physicore.reactions-diffusion.biofvm.${example}
                        PRIVATE physicore::reactions-diffusion::biofvm)
endforeach()
