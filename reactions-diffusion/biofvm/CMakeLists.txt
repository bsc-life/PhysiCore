cmake_minimum_required(VERSION 3.14)
project(BioFVM LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable all warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
  add_compile_options(/W4)
endif()

# Add source files (update as needed)
file(GLOB_RECURSE BIOFVM_SOURCES src/*.cpp)
add_library(biofvm ${BIOFVM_SOURCES})

target_include_directories(biofvm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(biofvm PUBLIC physicore-common)

file(GLOB_RECURSE TEST_SOURCES tests/*.cpp)

add_executable(biofvm_tests ${TEST_SOURCES})
target_link_libraries(biofvm_tests PRIVATE biofvm gtest_main)
gtest_discover_tests(biofvm_tests)
