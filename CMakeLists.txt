cmake_minimum_required(VERSION 3.24)
project(OblivousRouting)

set(CMAKE_CXX_STANDARD 20)

add_executable(OblivousRouting main.cpp
        src/graph.h
        src/graph.cpp)


add_subdirectory(src)
add_subdirectory(extern/googletest)
enable_testing()
add_subdirectory(test)
