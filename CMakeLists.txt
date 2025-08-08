cmake_minimum_required(VERSION 3.24)
project(ObliviousRouting)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")




# Tell CMake where Homebrew put Armadillo
set(ARMADILLO_INCLUDE_DIR "/usr/local/opt/armadillo/include")
set(ARMADILLO_LIBRARY     "/usr/local/opt/armadillo/lib/libarmadillo.dylib")

set(JULIA_DIR "/usr/local/opt/julia")
set(JULIA_INCLUDE_DIR "${JULIA_DIR}/include/julia"
        CACHE PATH "Path to the directory containing julia.h and julia_version.h"
)

set(JULIA_LIBRARY "${JULIA_DIR}/lib/libjulia.dylib")
message(STATUS "Using Julia in ${JULIA_DIR}")
# make sure our exe can find libjulia at runtime
link_directories(${JULIA_DIR}/lib)


# Add this after your other includes
# AMG solver

include_directories(/Users/halilibrahim/amgcl)


add_executable(oblivious_routing main.cpp
        src/graph.h
        src/graph.cpp
        src/utils/hash.h
        src/lp_solver/LPSolver.h
        src/lp_solver/LPSolver.cpp
        src/lp_solver/MCCF_lp_solver.cpp
        src/lp_solver/MCCF_lp_solver.h
        src/tree_based/mcct/Node.h
        src/tree_based/mcct/Tree.h
        src/tree_based/mcct/mcct_derandomized_weighted_solver.h
        src/tree_based/raecke_tree_decomp.cpp
        src/tree_based/raecke_tree_decomp.h
        src/tree_based/mcct/Node.cpp
        src/tree_based/mcct/mcct_derandomized_weighted_solver.cpp
        src/tree_based/raecke_transform.cpp
        src/tree_based/raecke_transform.h
        src/lp_solver/LP_Base.h
        src/electrical/electrical_flow_naive.cpp
        src/electrical/electrical_flow_naive.h
        src/electrical/laplcian_solver.cpp
        src/electrical/laplcian_solver.h
        src/electrical/AMGSolver.cpp
        src/electrical/AMGSolver.h
)

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

target_link_libraries(oblivious_routing PRIVATE ortools::ortools ${JULIA_LIBRARY})

# Include Armadillo headers and link against the library
target_include_directories(oblivious_routing PRIVATE ${ARMADILLO_INCLUDE_DIR})
target_link_libraries(oblivious_routing PRIVATE ${ARMADILLO_LIBRARY})
target_link_libraries(oblivious_routing PRIVATE ${JULIA_LIBRARY})
# target_link_libraries(oblivious_routing PRIVATE amgcl)
# On macOS you also want to embed the RPATH so it looks in JULIA_DIR/lib
# now set the RPATH to include both Julia *and* OR-Tools
set_target_properties(oblivious_routing PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH
        "${JULIA_ROOT}/lib;${JULIA_ROOT}/lib/julia;$<TARGET_FILE_DIR:ortools::ortools>"
)

add_link_options(-ljulia)

# Add the include directory for Julia
target_include_directories(oblivious_routing PRIVATE ${JULIA_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(extern/googletest)
enable_testing()
add_subdirectory(test)
