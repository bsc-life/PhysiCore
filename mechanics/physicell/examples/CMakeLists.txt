set(ALL_EXAMPLES mechanize)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/settings.xml
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/settings.xml
    ${CMAKE_CURRENT_BINARY_DIR}/settings.xml
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/settings.xml)

add_custom_target(mechanics.physicell.copy_xml ALL
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/settings.xml)

foreach(example ${ALL_EXAMPLES})
  add_executable(mechanics.physicell.${example} ${example}.cpp)

  target_link_libraries(mechanics.physicell.${example}
                        PRIVATE mechanics.physicell)

  add_dependencies(mechanics.physicell.${example} mechanics.physicell.copy_xml)
endforeach()
