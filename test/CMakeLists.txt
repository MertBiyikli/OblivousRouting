file(GLOB_RECURSE TEST_SRC_FILES smoke/*.cpp)

add_executable(smoke_test ${TEST_SRC_FILES})
target_link_libraries(smoke_test PRIVATE core gtest_main)
target_include_directories(smoke_test PRIVATE ${CMAKE_SOURCE_DIR}/src)

include(GoogleTest)
gtest_discover_tests(smoke_test)
