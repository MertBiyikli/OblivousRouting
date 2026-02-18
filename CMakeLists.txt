cmake_minimum_required(VERSION 3.20)
project(ObliviousRouting LANGUAGES CXX)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Options ---
option(BUILD_TESTS "Build test targets" OFF)
option(USE_OPENMP "Enable OpenMP support" ON)

option(OR_ENABLE_WARNINGS "Enable strict compiler warnings" ON)
option(OR_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

option(OR_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(OR_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(OR_ENABLE_TSAN "Enable ThreadSanitizer" OFF)

option(OR_ENABLE_LTO "Enable Link Time Optimization" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prefer consistent output dirs (handy for scripts)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Dependencies ---
find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(ortools REQUIRED CONFIG)
find_package(Eigen3 REQUIRED)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# --- Target ---
add_executable(oblivious_routing)
target_sources(oblivious_routing PRIVATE
        main.cpp
        src/utils/my_math.cpp
        src/datastructures/GraphADJ.cpp
        src/datastructures/GraphCSR.cpp
        experiments/performance/demands/GravityModel.cpp
        experiments/performance/demands/BimodalModel.cpp
        experiments/performance/demands/GaussianModel.cpp
        experiments/performance/demands/UniformModel.cpp
        src/tree_based/utils/ultrametric_tree.cpp
        src/electrical/electrical.cpp
        src/electrical/parallel/electrical_flow_par_batches.cpp
        src/electrical/laplacian_solvers/AMGSolver.cpp
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
)

# Includes
target_include_directories(oblivious_routing
        PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/experiments
        ${CMAKE_SOURCE_DIR}/extern/amgcl
)

# Link libs
target_link_libraries(oblivious_routing PRIVATE
        ortools::ortools
        Eigen3::Eigen
        Boost::program_options
        Boost::serialization
)

if(USE_OPENMP)
    target_link_libraries(oblivious_routing PRIVATE OpenMP::OpenMP_CXX)
endif()


if(OR_ENABLE_WARNINGS)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(oblivious_routing PRIVATE
                -Wall -Wextra -Wpedantic
                -Wshadow -Wconversion -Wsign-conversion
                -Wnon-virtual-dtor -Wold-style-cast
                -Woverloaded-virtual -Wnull-dereference
                -Wdouble-promotion -Wformat=2
        )
        if(OR_WARNINGS_AS_ERRORS)
            target_compile_options(oblivious_routing PRIVATE -Werror)
        endif()
    elseif(MSVC)
        target_compile_options(oblivious_routing PRIVATE /W4)
        if(OR_WARNINGS_AS_ERRORS)
            target_compile_options(oblivious_routing PRIVATE /WX)
        endif()
    endif()
endif()


if(OR_ENABLE_TSAN AND (OR_ENABLE_ASAN OR OR_ENABLE_UBSAN))
    message(FATAL_ERROR "TSAN cannot be combined with ASAN/UBSAN")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    if(OR_ENABLE_ASAN)
        target_compile_options(oblivious_routing PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_options(oblivious_routing PRIVATE -fsanitize=address)
    endif()

    if(OR_ENABLE_UBSAN)
        target_compile_options(oblivious_routing PRIVATE -fsanitize=undefined -fno-omit-frame-pointer)
        target_link_options(oblivious_routing PRIVATE -fsanitize=undefined)
    endif()

    if(OR_ENABLE_TSAN)
        target_compile_options(oblivious_routing PRIVATE -fsanitize=thread -fno-omit-frame-pointer)
        target_link_options(oblivious_routing PRIVATE -fsanitize=thread)
    endif()
endif()

# --- LTO ---
if(OR_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(ipo_supported)
        set_property(TARGET oblivious_routing PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "IPO/LTO not supported: ${ipo_error}")
    endif()
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP: ${USE_OPENMP}")
message(STATUS "Warnings: ${OR_ENABLE_WARNINGS} (Werror=${OR_WARNINGS_AS_ERRORS})")
message(STATUS "Sanitizers: ASAN=${OR_ENABLE_ASAN} UBSAN=${OR_ENABLE_UBSAN} TSAN=${OR_ENABLE_TSAN}")
message(STATUS "LTO: ${OR_ENABLE_LTO}")

# ---- Tests ----
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
