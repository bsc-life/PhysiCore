add_executable(physicore.phenotype.physicore)

target_sources(physicore.phenotype.physicore PRIVATE src/main.cpp)

target_link_libraries(
  physicore.phenotype.physicore PRIVATE physicore::mechanics::physicell
                                        physicore::reactions-diffusion::biofvm)
