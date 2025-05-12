cmake_minimum_required(VERSION 3.24)
project(ObliviousRouting)

set(CMAKE_CXX_STANDARD 20)


add_executable(oblivious_routing main.cpp
        src/graph.h
        src/graph.cpp)

include(FetchContent)

# Try to find system-installed OR-Tools first
find_package(ortools QUIET)

if(NOT ortools_FOUND)
    message(STATUS "OR-Tools not found. Fetching from GitHub...")

    FetchContent_Declare(
            ortools
            GIT_REPOSITORY https://github.com/google/or-tools.git
            GIT_TAG        v9.9 # you can update this to the latest stable release
    )

    set(BUILD_DEPS ON)
    set(BUILD_SAMPLES OFF)
    set(BUILD_TESTING OFF)
    set(BUILD_PYTHON OFF)
    set(BUILD_JAVA OFF)
    set(BUILD_CSHARP OFF)


    FetchContent_MakeAvailable(ortools)
endif()

target_link_libraries(oblivious_routing PRIVATE ortools::ortools)

add_subdirectory(src)
add_subdirectory(extern/googletest)
enable_testing()
add_subdirectory(test)
