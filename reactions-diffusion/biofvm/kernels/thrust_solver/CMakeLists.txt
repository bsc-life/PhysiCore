add_library(physicore.reactions-diffusion.kernels.tbb_thrust_solver)
add_library(physicore::reactions-diffusion::kernels::tbb_thrust_solver ALIAS
            physicore.reactions-diffusion.kernels.tbb_thrust_solver)

find_package(Thrust REQUIRED CONFIG)

thrust_create_target(ThrustTBB HOST TBB DEVICE TBB)

target_sources(
  physicore.reactions-diffusion.kernels.tbb_thrust_solver
  PRIVATE src/bulk_solver.cpp
          src/dirichlet_solver.cpp
          src/cell_solver.cpp
          src/diffusion_solver.cpp
          src/thrust_solver.cpp
          src/register_solver.cpp
          src/data_manager.cpp)

target_sources(
  physicore.reactions-diffusion.kernels.tbb_thrust_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.reactions-diffusion.kernels.tbb_thrust_solver
                      PRIVATE physicore::common ThrustTBB)

find_package(CUDAToolkit 12)

if(${CUDAToolkit_FOUND})
  enable_language(CUDA)

  add_library(physicore.reactions-diffusion.kernels.cuda_thrust_solver)
  add_library(physicore::reactions-diffusion::kernels::cuda_thrust_solver ALIAS
              physicore.reactions-diffusion.kernels.cuda_thrust_solver)

  get_target_property(
    TBB_SOURCES physicore.reactions-diffusion.kernels.tbb_thrust_solver SOURCES)

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

  target_sources(physicore.reactions-diffusion.kernels.cuda_thrust_solver
                 PRIVATE ${CUDA_SOURCES})

  set_target_properties(
    physicore.reactions-diffusion.kernels.cuda_thrust_solver
    PROPERTIES CUDA_ARCHITECTURES native CUDA_SEPARABLE_COMPILATION ON)

  thrust_create_target(ThrustCUDA HOST TBB DEVICE CUDA)

  target_compile_options(
    physicore.reactions-diffusion.kernels.cuda_thrust_solver
    PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda
           --expt-relaxed-constexpr>)

  target_sources(
    physicore.reactions-diffusion.kernels.cuda_thrust_solver
    PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

  target_link_libraries(physicore.reactions-diffusion.kernels.cuda_thrust_solver
                        PRIVATE physicore::common ThrustCUDA)
endif()

if(${CUDAToolkit_FOUND})
  set_target_properties(
    physicore.reactions-diffusion.biofvm
    PROPERTIES LINKER_LANGUAGE CUDA
               CUDA_ARCHITECTURES native
               CUDA_SEPARABLE_COMPILATION ON
               INTERFACE_HAS_CUDA TRUE)

  target_link_libraries(
    physicore.reactions-diffusion.biofvm
    PRIVATE physicore::reactions-diffusion::kernels::cuda_thrust_solver)
else()
  target_link_libraries(
    physicore.reactions-diffusion.biofvm
    PRIVATE physicore::reactions-diffusion::kernels::tbb_thrust_solver)
endif()

target_compile_definitions(physicore.reactions-diffusion.biofvm
                           PUBLIC PHYSICORE_HAS_THRUST)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
