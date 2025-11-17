# --------------------------------------------------
# Integration, dependency-free correctness tests
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
# INTEGRATION test-only library (no OR-Tools, no AMGCL)
# --------------------------------------------------

file(GLOB INTEGRATION_SRC
        ${CMAKE_SOURCE_DIR}/src/utils/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/mcct/*.cpp
        ${CMAKE_SOURCE_DIR}/src/graph/*.cpp
)

add_library(integration_core STATIC ${INTEGRATION_SRC})

target_include_directories(integration_core PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# --------------------------------------------------
# Automatically build integration test files
# --------------------------------------------------

file(GLOB INTEGRATION_TEST_FILES "*.cpp")

foreach(TEST_FILE ${INTEGRATION_TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)

    add_executable(${TEST_NAME} ${TEST_FILE})

    target_link_libraries(${TEST_NAME} PRIVATE
            integration_core
            GTest::gtest
            GTest::gtest_main
            Boost::boost
    )



    target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
    )

    gtest_discover_tests(${TEST_NAME})
endforeach()
