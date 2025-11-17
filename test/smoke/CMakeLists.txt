# --------------------------------------------------
# SMoke, dependency-free correctness tests
# --------------------------------------------------

cmake_minimum_required(VERSION 3.16)

# Find GoogleTest
find_package(GTest REQUIRED)
find_package(Boost REQUIRED)

# --------------------------------------------------
# Fix RPATH issues on macOS for GTest dynamic libs
# --------------------------------------------------
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH
        "${CMAKE_BUILD_RPATH};/Users/halilibrahim/opt/anaconda3/lib"
)

# --------------------------------------------------
# Smoke test-only library (no OR-Tools, no AMGCL)
# --------------------------------------------------

file(GLOB SMOKE_SRC
        ${CMAKE_SOURCE_DIR}/src/utils/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/mcct/*.cpp
        ${CMAKE_SOURCE_DIR}/src/electrical/*.cpp
        ${CMAKE_SOURCE_DIR}/src/electrical/laplacian_solvers/*.cpp
        ${CMAKE_SOURCE_DIR}/src/graph/*.cpp
)

add_library(smoke_core STATIC ${SMOKE_SRC})

target_include_directories(smoke_core PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# --------------------------------------------------
# Automatically build smoke test files
# --------------------------------------------------

file(GLOB SMOKE_TEST_FILES "*.cpp")

foreach(TEST_FILE ${SMOKE_TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)

    add_executable(${TEST_NAME} ${TEST_FILE})

    target_link_libraries(${TEST_NAME} PRIVATE
            smoke_core
            GTest::gtest
            GTest::gtest_main
            Boost::boost
    )



    target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
    )

    target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/extern/amgcl)

    gtest_discover_tests(${TEST_NAME})
endforeach()
