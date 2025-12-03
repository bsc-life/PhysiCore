add_executable(phenotype.physicore)

target_sources(phenotype.physicore PRIVATE src/main.cpp)

target_link_libraries(phenotype.physicore PRIVATE mechanics::micromechanics
                                                  reactions-diffusion::biofvm)
