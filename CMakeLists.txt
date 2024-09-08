cmake_minimum_required(VERSION 3.21)
project(fastppm)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(extern/pybind11)

pybind11_add_module(fastppm MODULE src/pymain.cpp src/Solver.cpp src/Solver.h src/State.cpp src/State.h src/Func.h src/Func.cpp src/Primal.cpp src/Primal.h src/Dual.cpp src/Dual.h)

add_executable(fastppm-test src/main.cpp src/Solver.cpp src/Solver.h src/State.cpp src/State.h src/Func.h src/Func.cpp src/Primal.cpp src/Primal.h src/Dual.cpp src/Dual.h src/Tree.cpp src/Tree.h src/Node.cpp src/Node.h src/Element.h)
