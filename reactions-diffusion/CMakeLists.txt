add_library(reactions-diffusion.interface INTERFACE)
add_library(reactions-diffusion::interface ALIAS reactions-diffusion.interface)

target_sources(
  reactions-diffusion.interface
  PUBLIC
    FILE_SET
    HEADERS
    BASE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/reactions_diffusion/reaction_diffusion_interface.h"
)

target_link_libraries(reactions-diffusion.interface INTERFACE common)

add_subdirectory(biofvm)

install(
  TARGETS reactions-diffusion.interface
  EXPORT BioFVMTargets
  FILE_SET HEADERS)
