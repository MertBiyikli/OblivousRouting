include_guard(GLOBAL)

# Keep dependency discovery in one place. CI pins the actual Linux package
# versions through the toolchain container; local developer machines may use
# compatible versions discovered by CMake.

find_package(fmt CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Cbc       REQUIRED IMPORTED_TARGET cbc)
pkg_check_modules(OsiCbc    REQUIRED IMPORTED_TARGET osi-cbc)
pkg_check_modules(Clp       REQUIRED IMPORTED_TARGET clp)
pkg_check_modules(OsiClp    REQUIRED IMPORTED_TARGET osi-clp)
pkg_check_modules(Cgl       REQUIRED IMPORTED_TARGET cgl)
pkg_check_modules(Osi       REQUIRED IMPORTED_TARGET osi)
pkg_check_modules(CoinUtils REQUIRED IMPORTED_TARGET coinutils)

find_package(ortools CONFIG REQUIRED)

# OR-Tools 9.12 installed-package workaround.
# Some 9.12 exports contain INTERFACE_SOURCES references to internal OBJECT
# libraries that are not exported. The project only consumes the installed
# OR-Tools libraries, so these interface object sources are unnecessary.
foreach(target
        ortools::ortools_math_opt
        ortools::ortools_math_opt_constraints)
    if(TARGET ${target})
        set_property(TARGET ${target} PROPERTY INTERFACE_SOURCES "")
    endif()
endforeach()

get_filename_component(ORTOOLS_PREFIX "${ortools_DIR}/../../.." REALPATH)
set(ORTOOLS_INCLUDE_DIR "${ORTOOLS_PREFIX}/include")

if(NOT EXISTS "${ORTOOLS_INCLUDE_DIR}/ortools/linear_solver/linear_solver.h")
    message(FATAL_ERROR
            "OR-Tools headers were not found under ${ORTOOLS_INCLUDE_DIR}. "
            "Check ortools_DIR/CMAKE_PREFIX_PATH.")
endif()

find_package(Boost REQUIRED COMPONENTS program_options serialization)
find_package(Eigen3 REQUIRED)

if(OR_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Prefer an already-installed Catch2 (CI installs the pinned version). For
# developer machines without Catch2, fetch the exact v3.5.0 commit.
if(BUILD_TESTING)
    find_package(Catch2 3 CONFIG QUIET)
    if(NOT Catch2_FOUND)
        include(FetchContent)
        set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
        message(STATUS "Catch2 not installed; fetching pinned revision")
        FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG 53d0d913a422d356b23dd927547febdf69ee9081
        )
        FetchContent_MakeAvailable(Catch2)
    endif()
endif()

# xcut expects spdlog to use the same fmt installation as this project.
set(SPDLOG_FMT_EXTERNAL ON CACHE BOOL "Use external fmt for spdlog" FORCE)

# Do not allow a developer-only, untracked xcut checkout to silently make the
# local build succeed while a clean CI checkout fails.
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/xcut/CMakeLists.txt")
    message(FATAL_ERROR
            "extern/xcut is missing. It must be tracked by the repository either "
            "as a pinned Git submodule or as vendored source.")
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/xcut" EXCLUDE_FROM_ALL)

if(NOT TARGET xcut_lib)
    message(FATAL_ERROR "extern/xcut did not define the expected target xcut_lib")
endif()

target_link_libraries(xcut_lib PUBLIC fmt::fmt)
