set(ALL_EXAMPLES diffuse)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/settings.xml
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/settings.xml
    ${CMAKE_CURRENT_BINARY_DIR}/settings.xml
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/settings.xml)

add_custom_target(copy_xml ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/settings.xml)

foreach(example ${ALL_EXAMPLES})
  add_executable(physicore.reactions-diffusion.biofvm.${example} ${example}.cpp)

  target_link_libraries(physicore.reactions-diffusion.biofvm.${example}
                        PRIVATE physicore.reactions-diffusion.biofvm)

  add_dependencies(physicore.reactions-diffusion.biofvm.${example} copy_xml)
endforeach()
