add_library(mechanics.physicell)
add_library(physicore::mechanics::physicell ALIAS mechanics.physicell)
add_library(mechanics.physicell_internal_iface INTERFACE)

file(GLOB PUBLIC_INC "include/physicell/*.h")
target_sources(
  mechanics.physicell
  PRIVATE src/environment.cpp src/config_reader.cpp src/agent_data.cpp
          src/mechanical_agent.cpp src/vtk_serializer_base.cpp
          src/vtk_agents_serializer.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  mechanics.physicell_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/physicell"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(pugixml CONFIG REQUIRED)
find_package(VTK CONFIG REQUIRED)

target_link_libraries(mechanics.physicell_internal_iface INTERFACE common
                                                                   VTK::IOXML)

target_link_libraries(
  mechanics.physicell
  PUBLIC common
  PRIVATE mechanics.physicell_internal_iface pugixml::pugixml VTK::IOXML)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
