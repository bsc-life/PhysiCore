add_library(physicore.reactions-diffusion.kernels.cpu_solver)
add_library(physicore::reactions-diffusion::kernels::cpu_solver ALIAS
            physicore.reactions-diffusion.kernels.cpu_solver)

find_package(hwy CONFIG REQUIRED)

find_package(OpenMP 4)

target_sources(
  physicore.reactions-diffusion.kernels.cpu_solver
  PRIVATE src/bulk_solver.cpp src/cell_solver.cpp src/diffusion_solver.cpp
          src/dirichlet_solver.cpp src/cpu_solver.cpp src/register_solver.cpp)

target_sources(
  physicore.reactions-diffusion.kernels.cpu_solver
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(physicore.reactions-diffusion.kernels.cpu_solver
                      PUBLIC physicore::common)

target_link_libraries(physicore.reactions-diffusion.kernels.cpu_solver
                      PRIVATE hwy::hwy)

if(OpenMP_CXX_FOUND)
  target_link_libraries(physicore.reactions-diffusion.kernels.cpu_solver
                        PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(physicore.reactions-diffusion.kernels.cpu_solver
                           PRIVATE noarr-structures/include)

if(PHYSICORE_BUILD_EXAMPLES)

  set(ALL_EXAMPLES benchmark)

  foreach(example ${ALL_EXAMPLES})
    add_executable(physicore.reactions-diffusion.kernels.cpu_solver.${example})
    target_sources(physicore.reactions-diffusion.kernels.cpu_solver.${example}
                   PRIVATE examples/${example}.cpp)
    target_link_libraries(
      physicore.reactions-diffusion.kernels.cpu_solver.${example}
      PRIVATE physicore::reactions-diffusion::kernels::cpu_solver
              physicore::reactions-diffusion::biofvm hwy::hwy)

    if(OpenMP_CXX_FOUND)
      target_link_libraries(
        physicore.reactions-diffusion.kernels.cpu_solver.${example}
        PRIVATE OpenMP::OpenMP_CXX)
    endif()

    target_include_directories(
      physicore.reactions-diffusion.kernels.cpu_solver.${example}
      PRIVATE noarr-structures/include)
  endforeach()

endif()

if(PHYSICORE_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(physicore.reactions-diffusion.kernels.cpu_solver.tests)

  file(GLOB TEST_SOURCES tests/*.cpp)
  target_sources(physicore.reactions-diffusion.kernels.cpu_solver.tests
                 PRIVATE ${TEST_SOURCES})

  target_link_libraries(
    physicore.reactions-diffusion.kernels.cpu_solver.tests
    PRIVATE physicore::reactions-diffusion::kernels::cpu_solver
            physicore::reactions-diffusion::biofvm GTest::gtest
            GTest::gtest_main)

  target_include_directories(
    physicore.reactions-diffusion.kernels.cpu_solver.tests
    PRIVATE noarr-structures/include)

  if(OpenMP_CXX_FOUND)
    target_link_libraries(physicore.reactions-diffusion.kernels.cpu_solver.tests
                          PRIVATE OpenMP::OpenMP_CXX)
  endif()

  include(GoogleTest)
  gtest_discover_tests(physicore.reactions-diffusion.kernels.cpu_solver.tests)
endif()
