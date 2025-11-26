# --------------------------------------------------
# Integration, dependency-free correctness tests
# --------------------------------------------------

cmake_minimum_required(VERSION 3.16)

# Find GoogleTest
find_package(GTest REQUIRED)


 # Find Boost
# only require what the main binary needs
find_package(Boost REQUIRED COMPONENTS program_options serialization)

# if you have tests that need unit_test_framework, guard it:
option(BUILD_TESTS "Build test targets" ON)

if(BUILD_TESTS)
    find_package(Boost 1.70 REQUIRED COMPONENTS unit_test_framework)
endif()

# --------------------------------------------------
# Fix RPATH issues on macOS for GTest dynamic libs
# --------------------------------------------------
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


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
            Boost::boost
            Boost::program_options
            Boost::serialization
            Boost::unit_test_framework
            integration_core
            GTest::gtest
            GTest::gtest_main
    )



    target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
            # add boost source as well
    )

    target_include_directories(integration_core PUBLIC
            ${CMAKE_SOURCE_DIR}/src
            ${Boost_INCLUDE_DIRS}     # <-- ADD THIS
    )

    target_link_libraries(integration_core PUBLIC
            Boost::boost              # <-- ADD THIS
    )


    gtest_discover_tests(${TEST_NAME})
endforeach()
