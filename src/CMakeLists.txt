add_executable(fastppm-cli Main.cpp Solver.cpp Solver.h PWL.cpp PWL.h Func.h Func.cpp Primal.cpp Primal.h Dual.cpp Dual.h)

# add libraries
target_link_libraries(fastppm-cli PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(fastppm-cli PRIVATE pprint)
target_link_libraries(fastppm-cli PRIVATE spdlog)
target_link_libraries(fastppm-cli PRIVATE csv)
target_link_libraries(fastppm-cli PRIVATE argparse)
