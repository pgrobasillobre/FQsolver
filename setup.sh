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

# Detect available processors for parallel test execution.
if command -v getconf >/dev/null 2>&1; then
    nproc=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
elif command -v sysctl >/dev/null 2>&1; then
    nproc=$(sysctl -n hw.ncpu 2>/dev/null)
else
    nproc=1
fi

case "$nproc" in
    ''|*[!0-9]*)
        nproc=1
        ;;
esac

if [ "$nproc" -lt 1 ]; then
    nproc=1
fi

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

# Configure convenient parallel test execution.
cat > run_tests.sh <<EOF
#!/bin/sh
CTEST_PARALLEL_LEVEL=${nproc} ctest -j ${nproc} --output-on-failure "\$@"
EOF
chmod +x run_tests.sh

cat > CTestCustom.cmake <<EOF
set(CTEST_PARALLEL_LEVEL ${nproc})
EOF

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
echo "$ ./run_tests.sh"
echo ""
echo "Parallel test level: $nproc"
echo ""
