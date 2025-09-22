add_library(physicore.reactions-diffusion.biofvm)
add_library(physicore::reactions-diffusion::biofvm ALIAS
            physicore.reactions-diffusion.biofvm)

target_sources(
  physicore.reactions-diffusion.biofvm
  PRIVATE src/microenvironment.cpp src/microenvironment_builder.cpp
          src/agent_data.cpp src/agent.cpp src/mesh.cpp)

target_sources(
  physicore.reactions-diffusion.biofvm
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/include/microenvironment.h")

target_link_libraries(physicore.reactions-diffusion.biofvm
                      PUBLIC physicore::common)

add_subdirectory(src/kernels/cpu_solver)

if(PHYSICORE_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(physicore.reactions-diffusion.biofvm.tests)

  file(GLOB TEST_SOURCES tests/*.cpp)
  target_sources(physicore.reactions-diffusion.biofvm.tests
                 PRIVATE ${TEST_SOURCES})

  target_link_libraries(
    physicore.reactions-diffusion.biofvm.tests
    PRIVATE physicore::reactions-diffusion::biofvm GTest::gtest
            GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(physicore.reactions-diffusion.biofvm.tests)
endif()
