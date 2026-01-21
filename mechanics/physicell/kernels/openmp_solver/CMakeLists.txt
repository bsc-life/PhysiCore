add_library(mechanics.physicell.kernels.openmp_solver_internal_iface INTERFACE)
add_library(mechanics.physicell.kernels.openmp_solver OBJECT)

if(${BUILD_SHARED_LIBS})
  set_target_properties(mechanics.physicell.kernels.openmp_solver
                        PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

find_package(OpenMP 4)

target_sources(
  mechanics.physicell.kernels.openmp_solver
  PRIVATE src/openmp_solver.cpp src/register_solver.cpp src/position_solver.cpp
          src/common_solver.cpp
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         "${CMAKE_CURRENT_SOURCE_DIR}/include"
         FILES
         include/physicell/openmp_solver/register_solver.h
         include/physicell/openmp_solver/position_solver.h)

target_link_libraries(mechanics.physicell.kernels.openmp_solver_internal_iface
                      INTERFACE mechanics.physicell_internal_iface common)

target_include_directories(
  mechanics.physicell.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

if(OpenMP_CXX_FOUND)
  target_link_libraries(mechanics.physicell.kernels.openmp_solver_internal_iface
                        INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  mechanics.physicell.kernels.openmp_solver
  PRIVATE mechanics.physicell.kernels.openmp_solver_internal_iface)

target_link_libraries(mechanics.physicell
                      PRIVATE mechanics.physicell.kernels.openmp_solver)

if(PHYSICORE_BUILD_EXAMPLES AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples")
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
