add_library(mechanics.physicell)
add_library(mechanics::physicell ALIAS mechanics.physicell)
add_library(mechanics.physicell_internal_iface INTERFACE)

if(${BUILD_SHARED_LIBS})
  set_target_properties(
    mechanics.physicell
    PROPERTIES CXX_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON
               POSITION_INDEPENDENT_CODE ON
               INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

file(GLOB PUBLIC_INC "include/physicell/*.h")
target_sources(
  mechanics.physicell
  PRIVATE src/environment.cpp
          src/environment_builder.cpp
          src/config_reader.cpp
          src/mechanical_agent_data.cpp
          src/mechanical_agent.cpp
          src/vtk_serializer_base.cpp
          src/vtk_agents_serializer.cpp
          src/solver_registry.cpp
  PUBLIC
    FILE_SET
    HEADERS
    BASE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES
    ${PUBLIC_INC}
    BASE_DIRS
    "${CMAKE_SOURCE_DIR}/reactions-diffusion/include"
    FILES
    ${CMAKE_SOURCE_DIR}/reactions-diffusion/include/reactions_diffusion/reactions_diffusion_interface.h
)

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

add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_EXAMPLES AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples")
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
