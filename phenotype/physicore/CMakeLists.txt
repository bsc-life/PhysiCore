add_executable(physicore.phenotype.physicore)

target_sources(physicore.phenotype.physicore PRIVATE src/main.cpp)

target_link_libraries(
  physicore.phenotype.physicore PUBLIC physicore::mechanics::physicell
                                       physicore::reactions-diffusion::biofvm)

get_target_property(HAS_CUDA physicore::reactions-diffusion::biofvm
                    INTERFACE_HAS_CUDA)
if(HAS_CUDA)
  enable_language(CUDA)
endif()
