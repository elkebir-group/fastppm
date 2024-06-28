cmake_minimum_required(VERSION 3.21)
project(EfficientLLHEstimator)

set(CMAKE_CXX_STANDARD 14)

add_executable(EfficientLLHEstimator src/main.cpp src/Solver.cpp src/Solver.h src/PWL.cpp src/PWL.h src/Func.h src/Func.cpp)
