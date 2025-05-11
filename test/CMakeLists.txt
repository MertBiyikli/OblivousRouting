add_executable(smoke_test smoke/smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE core gtest_main)

# Add correct include paths
target_include_directories(smoke_test PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${gtest_SOURCE_DIR}/include
)
