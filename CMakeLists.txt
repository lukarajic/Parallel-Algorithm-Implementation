cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)
project(ParallelNBody VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP
if(APPLE)
    set(OpenMP_ROOT /opt/homebrew/opt/libomp)
endif()
find_package(OpenMP REQUIRED)

add_executable(nbody src/main.cpp src/init.cpp src/utils.cpp)

target_include_directories(nbody PUBLIC include)

# Optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(nbody PRIVATE -O3 -march=native)
endif()

target_link_libraries(nbody PUBLIC OpenMP::OpenMP_CXX)

# Testing
enable_testing()
add_executable(test_nbody tests/test_nbody.cpp src/utils.cpp)
target_include_directories(test_nbody PUBLIC include)
target_link_libraries(test_nbody PUBLIC OpenMP::OpenMP_CXX)
add_test(NAME test_nbody COMMAND test_nbody)
