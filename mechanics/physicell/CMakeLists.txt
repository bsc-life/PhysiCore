add_library(physicore.mechanics.physicell)
add_library(physicore::mechanics::physicell ALIAS physicore.mechanics.physicell)
add_library(physicore.mechanics.physicell_internal_iface INTERFACE)

file(GLOB PUBLIC_INC "include/physicell/*.h")
target_sources(
  physicore.mechanics.physicell
  PRIVATE src/environment.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  physicore.mechanics.physicell_internal_iface
  INTERFACE FILE_SET HEADERS
            BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include/physicell"
                      "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(physicore.mechanics.physicell_internal_iface
                      INTERFACE physicore::common)

target_link_libraries(physicore.mechanics.physicell
                      PRIVATE physicore.mechanics.physicell_internal_iface)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
