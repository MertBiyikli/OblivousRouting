# --------------------------------------------------
# Tiny, dependency-free correctness tests
# --------------------------------------------------

cmake_minimum_required(VERSION 3.16)

# Find GoogleTest
find_package(GTest REQUIRED)

# --------------------------------------------------
# Fix RPATH issues on macOS for GTest dynamic libs
# --------------------------------------------------
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH
        "${CMAKE_BUILD_RPATH};/Users/halilibrahim/opt/anaconda3/lib"
)

# --------------------------------------------------
# Tiny test-only library (no OR-Tools, no AMGCL)
# --------------------------------------------------

file(GLOB TINY_SRC
        ${CMAKE_SOURCE_DIR}/src/utils/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/*.cpp
        ${CMAKE_SOURCE_DIR}/src/tree_based/frt/mcct/*.cpp
        ${CMAKE_SOURCE_DIR}/src/graph/*.cpp
)

add_library(tiny_core STATIC ${TINY_SRC})

target_include_directories(tiny_core PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

# --------------------------------------------------
# Automatically build tiny test files
# --------------------------------------------------

file(GLOB TINY_TEST_FILES "*.cpp")

foreach(TEST_FILE ${TINY_TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)

    add_executable(${TEST_NAME} ${TEST_FILE})

    target_link_libraries(${TEST_NAME} PRIVATE
            tiny_core
            GTest::gtest
            GTest::gtest_main
    )

    target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
    )

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
