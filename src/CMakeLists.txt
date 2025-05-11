file(GLOB SRC_FILES *.cpp)
add_library(core ${SRC_FILES}
        solver/solver.h)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
