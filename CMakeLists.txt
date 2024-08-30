cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native")

project(fastppm VERSION 1.0)

add_subdirectory(src/)

target_include_directories(fastppm-cli PUBLIC "${fastppm_SOURCE_DIR}/third-party/csv-parser/include")

# turn off JSON library tests
set(JSON_BuildTests OFF CACHE INTERNAL "")

# compile libraries
add_subdirectory(extern/spdlog)
add_subdirectory(extern/pprint)
add_subdirectory(extern/csv-parser)
add_subdirectory(extern/argparse)
add_subdirectory(extern/json)
add_subdirectory(extern/pybind11)

# TODO: fix at end
# pybind11_add_module(
# fastppm MODULE 
# src/pymain.cpp 
# src/Solver.cpp src/Solver.h 
# src/PWL.cpp src/PWL.h 
# src/Func.h src/Func.cpp 
# src/Primal.cpp src/Primal.h 
# src/Dual.cpp src/Dual.h
# )
