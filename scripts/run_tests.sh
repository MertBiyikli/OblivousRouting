#!/bin/bash

# ObliviousRouting Test Runner Script
# This script helps build and run the test suite

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Defaults
BUILD_DIR="cmake-build-tests"
CMAKE_FLAGS="-DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release"
RUN_TESTS=true
VERBOSE=false
PARALLEL=""
TEST_FILTER=""

# Usage information
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help              Show this help message"
    echo "  -b, --build-dir DIR     Custom build directory (default: $BUILD_DIR)"
    echo "  -j, --parallel N        Number of parallel jobs (default: auto-detect)"
    echo "  -f, --filter PATTERN    Run only tests matching pattern"
    echo "  -v, --verbose           Verbose test output"
    echo "  -c, --clean             Clean build directory before building"
    echo "  --no-test               Build only, don't run tests"
    echo "  --debug                 Build in Debug mode instead of Release"
    echo ""
    echo "Examples:"
    echo "  # Build and run all tests"
    echo "  $0"
    echo ""
    echo "  # Run only graph tests with verbose output"
    echo "  $0 -f GraphCSR -v"
    echo ""
    echo "  # Build with 4 parallel jobs, don't run tests"
    echo "  $0 -j 4 --no-test"
    echo ""
    echo "  # Clean build with debug info"
    echo "  $0 -c --debug"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -j|--parallel)
            PARALLEL="-j$2"
            shift 2
            ;;
        -f|--filter)
            TEST_FILTER="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -c|--clean)
            echo "Cleaning $BUILD_DIR..."
            rm -rf "$BUILD_DIR"
            shift
            ;;
        --no-test)
            RUN_TESTS=false
            shift
            ;;
        --debug)
            CMAKE_FLAGS="${CMAKE_FLAGS/Release/Debug}"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Detect number of parallel jobs if not specified
if [ -z "$PARALLEL" ]; then
    if command -v nproc &> /dev/null; then
        PARALLEL="-j$(nproc)"
    elif [ "$(uname)" = "Darwin" ]; then
        PARALLEL="-j$(sysctl -n hw.ncpu)"
    else
        PARALLEL="-j4"  # Fallback
    fi
fi

echo -e "${GREEN}=== ObliviousRouting Test Suite ===${NC}"
echo ""
echo -e "Configuration:"
echo "  Build Directory: $BUILD_DIR"
echo "  CMake Flags: $CMAKE_FLAGS"
echo "  Parallel Jobs: $PARALLEL"
if [ -n "$TEST_FILTER" ]; then
    echo "  Test Filter: $TEST_FILTER"
fi
echo ""

# Step 1: Configure with CMake
echo -e "${YELLOW}[1/3] Configuring CMake...${NC}"
cmake -S . -B "$BUILD_DIR" -G Ninja $CMAKE_FLAGS
if [ $? -ne 0 ]; then
    echo -e "${RED}CMake configuration failed!${NC}"
    exit 1
fi
echo -e "${GREEN}CMake configuration succeeded${NC}"
echo ""

# Step 2: Build
echo -e "${YELLOW}[2/3] Building tests...${NC}"
cmake --build "$BUILD_DIR" $PARALLEL
if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi
echo -e "${GREEN}Build succeeded${NC}"
echo ""

# Step 3: Run tests (if enabled)
if [ "$RUN_TESTS" = true ]; then
    echo -e "${YELLOW}[3/3] Running tests...${NC}"
    echo ""

    cd "$BUILD_DIR"

    # Build ctest command
    CTEST_CMD="ctest --output-on-failure"

    if [ "$VERBOSE" = true ]; then
        CTEST_CMD="$CTEST_CMD -VV"
    fi

    if [ -n "$TEST_FILTER" ]; then
        CTEST_CMD="$CTEST_CMD -R '$TEST_FILTER'"
    fi

    # Add parallel jobs
    CTEST_CMD="$CTEST_CMD $PARALLEL"

    echo "Running: $CTEST_CMD"
    echo ""

    eval "$CTEST_CMD"
    TEST_RESULT=$?

    cd - > /dev/null

    if [ $TEST_RESULT -eq 0 ]; then
        echo ""
        echo -e "${GREEN}=== All tests passed! ===${NC}"
    else
        echo ""
        echo -e "${RED}=== Some tests failed ===${NC}"
        echo ""
        echo "To debug specific test:"
        echo "  $BUILD_DIR/unit_tests '[TestPattern]'"
        exit 1
    fi
else
    echo -e "${GREEN}=== Build completed successfully ===${NC}"
    echo ""
    echo "To run tests:"
    echo "  cd $BUILD_DIR && ctest --output-on-failure"
    echo "Or:"
    echo "  ./$BUILD_DIR/unit_tests"
fi

