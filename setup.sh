#!/bin/sh
# Setup script for building FQSolver

README=README.md
omp_mode="auto"
build="build"

# Parse CLI options
while test $# -gt 0; do
    case "$1" in
        -h|--help)
            echo "Usage: ./setup.sh [options]"
            echo "Options:"
            echo "  -b, --build <dir>    Build directory name (default: build)"
            echo "  -omp, --omp          Require OpenMP; fail if it is not available"
            echo "  --serial             Disable OpenMP"
            echo ""
            echo "By default, setup tries OpenMP first and falls back to serial if needed."
            exit 0
            ;;
        -b|--build)
            shift
            build=$1
            shift
            ;;
        -omp|--omp)
            shift
            omp_mode="required"
            ;;
        --serial)
            shift
            omp_mode="serial"
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

# Run cmake
if [ "$omp_mode" = "serial" ]; then
    cmake .. $cmake_args -DENABLE_OMP=OFF
elif [ "$omp_mode" = "required" ]; then
    cmake .. $cmake_args -DENABLE_OMP=ON
else
    echo "Trying OpenMP build..."
    cmake .. $cmake_args -DENABLE_OMP=ON
    if [ $? -ne 0 ]; then
        echo ""
        echo "OpenMP was not found. Falling back to serial build..."
        cmake .. $cmake_args -DENABLE_OMP=OFF
    fi
fi

cmake_status=$?
if [ $cmake_status -ne 0 ]; then
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
