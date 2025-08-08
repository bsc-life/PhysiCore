cmake_minimum_required(VERSION 3.14)

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
file(GLOB_RECURSE COMMON_SOURCES src/*.cpp)
add_library(physicore-common ${COMMON_SOURCES})

target_include_directories(physicore-common
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
