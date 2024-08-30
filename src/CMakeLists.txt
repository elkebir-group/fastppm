add_executable(fastppm-cli 
    Main.cpp 
    Solvers/LogBinomialPiecewise/Dual.cpp 
    Solvers/LogBinomialPiecewise/Primal.cpp
    Solvers/LogBinomialPiecewise/Solver.cpp
    Solvers/LogBinomialPiecewise/PWL.cpp
    Solvers/LogBinomialPiecewise/Func.cpp
)


# add libraries
target_link_libraries(fastppm-cli PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(fastppm-cli PRIVATE pprint)
target_link_libraries(fastppm-cli PRIVATE spdlog)
target_link_libraries(fastppm-cli PRIVATE csv)
target_link_libraries(fastppm-cli PRIVATE argparse)
