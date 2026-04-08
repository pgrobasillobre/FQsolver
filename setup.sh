#!/bin/sh
# Setup script for building FQSolver

README=README.md
ompvar=0
build="build"

# Parse CLI options
while test $# -gt 0; do
    case "$1" in
        -h|--help)
            echo "Usage: ./setup.sh [options]"
            echo "Options:"
            echo "  -b, --build <dir>    Build directory name (default: build)"
            echo "  -omp, --omp          Enable OpenMP"
            exit 0
            ;;
        -b|--build)
            shift
            build=$1
            shift
            ;;
        -omp|--omp)
            shift
            ompvar=1
            ;;
        *)
            break
            ;;
    esac
done

# Set final build directory
buildir=$build

# Create and enter build directory
mkdir -p "$buildir"
cd "$buildir"

# Compose cmake arguments
cmake_args="-DENABLE_AUTO_BLAS=ON -DENABLE_AUTO_LAPACK=ON"
if [ $ompvar -eq 1 ]; then
    cmake_args="$cmake_args -DENABLE_OMP=ON"
fi

# Run cmake
cmake .. $cmake_args

# Check if cmake succeeded
if [ $? -ne 0 ]; then
    echo "❌ CMake configuration failed. Check the errors above."
    exit 1
fi

# Final user instructions
echo ""
echo "✅ CMake configuration complete."
echo ""
echo "To compile:"
echo "$ cd $buildir"
echo "$ make -j"
echo ""
echo "To test:"
echo "$ cd $buildir"
echo "$ ctest"
echo ""

