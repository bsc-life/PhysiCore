# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.24)

# Get version from git tag
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  if(GIT_TAG)
    # Remove 'v' prefix if present
    string(REGEX REPLACE "^v" "" PROJECT_VERSION_FROM_GIT "${GIT_TAG}")
  endif()

  # Fallback version if git tag not found
  if(NOT PROJECT_VERSION_FROM_GIT)
    set(PROJECT_VERSION_FROM_GIT "9.9.9")
  endif()
endif()

project(
  PhysiCore
  LANGUAGES CXX
  VERSION ${PROJECT_VERSION_FROM_GIT})

option(
  PHYSICORE_BUILD_TESTS
  "Enable building tests and test infrastructure. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL})

option(
  PHYSICORE_BUILD_EXAMPLES
  "Enable building examples. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL})

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

if(BUILD_SHARED_LIBS)
  set(PHYSICORE_BUILD_TESTS
      OFF
      CACHE
        BOOL
        "Shared libraries are incompatible with tests due to symbol visibility."
        FORCE)
endif()

include(CTest)

add_subdirectory(common)
add_subdirectory(reactions-diffusion/biofvm)
add_subdirectory(mechanics/micromechanics)
add_subdirectory(mechanics/physicell)
add_subdirectory(phenotype/physicore)
