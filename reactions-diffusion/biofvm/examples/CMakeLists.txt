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
  add_executable(reactions-diffusion.biofvm.${example} ${example}.cpp)

  target_link_libraries(reactions-diffusion.biofvm.${example}
                        PRIVATE reactions-diffusion.biofvm)

  add_dependencies(reactions-diffusion.biofvm.${example} copy_xml)
endforeach()
