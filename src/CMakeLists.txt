add_executable(fastppm-cli 
    Main.cpp
    Lib.cpp 
    Solvers/LogBinomialPiecewise/Dual.cpp 
    Solvers/LogBinomialPiecewise/Primal.cpp
    Solvers/LogBinomialPiecewise/Solver.cpp
    Solvers/LogBinomialPiecewise/Node.cpp
    Solvers/LogBinomialPiecewise/State.cpp
    Solvers/LogBinomialPiecewise/Tree.cpp
    Solvers/LogBinomialPiecewise/Func.cpp
    Solvers/L2/Solver.cpp
)

# add libraries to the executable
target_link_libraries(fastppm-cli PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(fastppm-cli PRIVATE pprint)
target_link_libraries(fastppm-cli PRIVATE spdlog)
target_link_libraries(fastppm-cli PRIVATE csv)
target_link_libraries(fastppm-cli PRIVATE argparse)

# compile python bindings
pybind11_add_module(
    fastppm MODULE 
    PyMain.cpp 
    Lib.cpp
    Solvers/LogBinomialPiecewise/Dual.cpp 
    Solvers/LogBinomialPiecewise/Primal.cpp
    Solvers/LogBinomialPiecewise/Solver.cpp
    Solvers/LogBinomialPiecewise/Node.cpp
    Solvers/LogBinomialPiecewise/State.cpp
    Solvers/LogBinomialPiecewise/Tree.cpp
    Solvers/LogBinomialPiecewise/Func.cpp
    Solvers/L2/Solver.cpp
)

# add libraries to the python bindings
target_link_libraries(fastppm PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(fastppm PRIVATE pprint)
target_link_libraries(fastppm PRIVATE spdlog)
target_link_libraries(fastppm PRIVATE csv)
target_link_libraries(fastppm PRIVATE argparse)
