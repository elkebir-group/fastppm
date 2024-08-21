cmake_minimum_required(VERSION 3.21)
project(EfficientLLHEstimator)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(extern/pybind11)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(test MODULE src/pymain.cpp src/Solver.cpp src/Solver.h src/PWL.cpp src/PWL.h src/Func.h src/Func.cpp src/Primal.cpp src/Primal.h src/Dual.cpp src/Dual.h)

add_executable(EfficientLLH src/main.cpp src/Solver.cpp src/Solver.h src/PWL.cpp src/PWL.h src/Func.h src/Func.cpp src/Primal.cpp src/Primal.h src/Dual.cpp src/Dual.h)
