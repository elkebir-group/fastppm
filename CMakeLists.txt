cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# turn off JSON library tests
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# setup PyBind11
set(PYBIND11_FINDPYTHON ON)

# compile libraries
add_subdirectory(extern/spdlog)
add_subdirectory(extern/pprint)
add_subdirectory(extern/csv-parser)
add_subdirectory(extern/argparse)
add_subdirectory(extern/yyjson)
add_subdirectory(extern/pybind11)

# compile the project
project(fastppm VERSION 1.0)
add_subdirectory(include/)
add_subdirectory(src/)
