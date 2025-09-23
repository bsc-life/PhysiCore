# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(
  PhysiCore
  LANGUAGES CXX
  VERSION 0.1.0)

option(
  PHYSICORE_BUILD_TESTS
  "Enable building tests and test infrastructure. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL})

option(
  PHYSICORE_BUILD_EXAMPLES
  "Enable building examples. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL})

include(CTest)

add_subdirectory(common)
add_subdirectory(reactions-diffusion/biofvm)
add_subdirectory(mechanics/physicell)
add_subdirectory(phenotype/physicore)
