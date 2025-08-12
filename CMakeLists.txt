# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
project(PhysiCore LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(COVERAGE_BUILD OFF "Enable coverage build")

if(${COVERAGE_BUILD})
  set(CMAKE_CXX_FLAGS
      " ${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping -mllvm -runtime-counter-relocation"
  )
endif()

# FetchContent for GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/heads/main.zip
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
enable_testing()

add_subdirectory(common)
add_subdirectory(reactions-diffusion/biofvm)
add_subdirectory(mechanics/physicell)

# Add main executable
add_executable(physicore src/main.cpp)
target_link_libraries(physicore PRIVATE biofvm micromechanics)
