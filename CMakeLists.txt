cmake_minimum_required(VERSION 3.10)
# --- VERSION DETECTION LOGIC ---
# Priority 1: Check if passed via command line (e.g. -DFASTPPM_VERSION=1.1.0)
if(NOT DEFINED FASTPPM_VERSION)
  # Priority 2: Try to retrieve version from Git tags (matches meta.yaml behavior)
  find_package(Git)
  if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
    )

    # If we found a tag (e.g. "v1.1.0"), strip the leading 'v'
    if(GIT_TAG_VERSION MATCHES "^v")
      string(SUBSTRING "${GIT_TAG_VERSION}" 1 -1 GIT_TAG_VERSION)
    endif()
  endif()

  # Priority 3: Fallback to default
  if(GIT_TAG_VERSION)
    set(FASTPPM_VERSION "${GIT_TAG_VERSION}")
  else()
    set(FASTPPM_VERSION "0.0.0") # Fallback for local builds with no git
  endif()
endif()
project(fastppm VERSION ${FASTPPM_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
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

add_subdirectory(include/)
add_subdirectory(src/)
