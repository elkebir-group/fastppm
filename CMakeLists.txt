cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native")

# turn off JSON library tests
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# compile libraries
add_subdirectory(extern/spdlog)
add_subdirectory(extern/pprint)
add_subdirectory(extern/csv-parser)
add_subdirectory(extern/argparse)
add_subdirectory(extern/json)
add_subdirectory(extern/pybind11)

# compile the project
project(fastppm VERSION 1.0)
add_subdirectory(src/)
