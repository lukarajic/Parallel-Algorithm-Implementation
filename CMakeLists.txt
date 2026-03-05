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

add_executable(nbody src/main.cpp)

target_include_directories(nbody PUBLIC include)

target_link_libraries(nbody PUBLIC OpenMP::OpenMP_CXX)
