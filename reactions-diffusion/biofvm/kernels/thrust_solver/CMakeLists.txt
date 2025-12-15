add_library(reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface
            INTERFACE)
add_library(reactions-diffusion.biofvm.kernels.tbb_thrust_solver OBJECT)

if(${BUILD_SHARED_LIBS})
  set_target_properties(reactions-diffusion.biofvm.kernels.tbb_thrust_solver
                        PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

find_package(Thrust REQUIRED CONFIG)

thrust_create_target(ThrustTBB HOST TBB DEVICE TBB)

# Configure namespace header for TBB backend
set(THRUST_SOLVER_NAMESPACE "tbb_thrust_solver")
set(THRUST_SOLVER_REGISTRY_NAME "tbb_thrust_solver")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/namespace_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/tbb/namespace_config.h" @ONLY)

target_sources(
  reactions-diffusion.biofvm.kernels.tbb_thrust_solver
  PRIVATE src/bulk_solver.cpp
          src/dirichlet_solver.cpp
          src/cell_solver.cpp
          src/diffusion_solver.cpp
          src/thrust_solver.cpp
          src/register_solver.cpp
          src/data_manager.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:reactions-diffusion.biofvm_iface> common
            ThrustTBB)

target_include_directories(
  reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/tbb>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

target_link_libraries(
  reactions-diffusion.biofvm.kernels.tbb_thrust_solver
  PRIVATE reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface)

find_path(NOARR_STRUCTURES_INCLUDE_DIRS "noarr/introspection.hpp")
target_include_directories(
  reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${NOARR_STRUCTURES_INCLUDE_DIRS}>)

find_package(CUDAToolkit 12)

if(${CUDAToolkit_FOUND})
  enable_language(CUDA)
  add_library(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
    INTERFACE)
  add_library(reactions-diffusion.biofvm.kernels.cuda_thrust_solver OBJECT)

  if(${BUILD_SHARED_LIBS})
    set_target_properties(reactions-diffusion.biofvm.kernels.cuda_thrust_solver
                          PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif()

  get_target_property(
    TBB_SOURCES reactions-diffusion.biofvm.kernels.tbb_thrust_solver SOURCES)

  set(CUDA_SOURCES "")

  # Copy source into build dir with .cu extension
  foreach(src ${TBB_SOURCES})
    get_filename_component(fname ${src} NAME_WE)
    set(dst "${CMAKE_CURRENT_BINARY_DIR}/cuda_sources/${fname}.cu")

    add_custom_command(
      OUTPUT ${dst}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${CMAKE_CURRENT_SOURCE_DIR}/${src} ${dst}
      DEPENDS ${src})

    list(APPEND CUDA_SOURCES ${dst})
  endforeach()

  # Configure namespace header for CUDA backend
  set(THRUST_SOLVER_NAMESPACE "cuda_thrust_solver")
  set(THRUST_SOLVER_REGISTRY_NAME "cuda_thrust_solver")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/namespace_config.h.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/cuda/namespace_config.h" @ONLY)

  target_sources(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver
    PRIVATE ${CUDA_SOURCES}
    PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

  thrust_create_target(ThrustCUDA HOST TBB DEVICE CUDA)

  target_compile_options(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver
    PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda
           --expt-relaxed-constexpr>)

  target_link_libraries(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
    INTERFACE $<BUILD_INTERFACE:reactions-diffusion.biofvm_iface> common
              ThrustCUDA)

  target_include_directories(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/cuda>
              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

  target_link_libraries(
    reactions-diffusion.biofvm.kernels.cuda_thrust_solver
    PRIVATE reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
  )
  
find_path(NOARR_STRUCTURES_INCLUDE_DIRS "noarr/introspection.hpp")
target_include_directories(
  reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${NOARR_STRUCTURES_INCLUDE_DIRS}>)
endif()

target_link_libraries(
  reactions-diffusion.biofvm
  PRIVATE reactions-diffusion.biofvm.kernels.tbb_thrust_solver)

target_compile_definitions(reactions-diffusion.biofvm
                           PRIVATE PHYSICORE_HAS_TBB_THRUST)

install(
  TARGETS reactions-diffusion.biofvm.kernels.tbb_thrust_solver
          reactions-diffusion.biofvm.kernels.tbb_thrust_solver_internal_iface
  EXPORT BioFVMTargets
  FILE_SET HEADERS)

if(${CUDAToolkit_FOUND})
  set_target_properties(
    reactions-diffusion.biofvm
    PROPERTIES LINKER_LANGUAGE CUDA
               CUDA_SEPARABLE_COMPILATION ON
               INTERFACE_HAS_CUDA TRUE)

  target_link_libraries(
    reactions-diffusion.biofvm
    PRIVATE reactions-diffusion.biofvm.kernels.cuda_thrust_solver)

  target_compile_definitions(reactions-diffusion.biofvm
                             PRIVATE PHYSICORE_HAS_CUDA_THRUST)

  install(
    TARGETS reactions-diffusion.biofvm.kernels.cuda_thrust_solver
            reactions-diffusion.biofvm.kernels.cuda_thrust_solver_internal_iface
    EXPORT BioFVMTargets
    FILE_SET HEADERS)
endif()

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
