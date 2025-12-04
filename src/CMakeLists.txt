# Build the fastppm library

add_library(fastppm-lib STATIC 
    Lib.cpp
    Solvers/PiecewiseLinear/Dual.cpp
    Solvers/PiecewiseLinear/Primal.cpp
    Solvers/PiecewiseLinear/Solver.cpp
    Solvers/PiecewiseLinear/Node.cpp
    Solvers/PiecewiseLinear/State.cpp
    Solvers/PiecewiseLinear/Tree.cpp
    Solvers/PiecewiseLinear/Func.cpp
    Solvers/SeparableADMM/Solver.cpp
    Solvers/L2/Solver.cpp
)

# add libraries to the executable
target_link_libraries(fastppm-lib PRIVATE
    yyjson
    pprint
    spdlog
    csv
)

# Pass CMake version variables to the C++ preprocessor
target_compile_definitions(fastppm-lib PRIVATE
        FASTPPM_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        FASTPPM_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        FASTPPM_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)


target_include_directories(fastppm-lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Build the fastppm CLI
add_executable(fastppm-cli
    Main.cpp
)

target_link_libraries(fastppm-cli PRIVATE
    fastppm-lib
    argparse
    yyjson
    pprint
    spdlog
    csv
)

# Pass CMake version variables to the C++ preprocessor
target_compile_definitions(fastppm-cli PRIVATE
        FASTPPM_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        FASTPPM_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        FASTPPM_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)


# Build the fastppm Python bindings
pybind11_add_module(
    fastppm MODULE 
    PyMain.cpp 
)

# Pass CMake version variables to the C++ preprocessor
target_compile_definitions(fastppm PRIVATE
        FASTPPM_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        FASTPPM_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        FASTPPM_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)

# add libraries to the python bindings
target_link_libraries(fastppm PRIVATE
    fastppm-lib
    yyjson
    pprint
    spdlog
    csv
    argparse
)
