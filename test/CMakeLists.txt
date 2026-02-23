include(FetchContent)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

# Windows users: avoid overriding parent project runtime settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(oblivious_tests
        main.cpp   # or test_*.cpp files
)

target_link_libraries(oblivious_tests
        PRIVATE
        oblivious_core
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oblivious_tests)