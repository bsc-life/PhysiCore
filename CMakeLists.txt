# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
project(PhysiCore LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add subdirectory for biofvm (which also adds common)
add_subdirectory(reactions-diffusion/biofvm)

# Add main executable
add_executable(physicore src/main.cpp)
target_link_libraries(physicore PRIVATE biofvm)

# Make ctest see all tests
enable_testing()
