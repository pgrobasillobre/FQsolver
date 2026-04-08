# ========================
# Top-level CMakeLists.txt
# ========================

cmake_minimum_required(VERSION 3.10)

project(FQSolver)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Extra module paths (FindMKL.cmake, etc.)
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
)

# ------------------------
# Options
# ------------------------
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_THREADED_MKL   "Enable OpenMP parallelization in MKL"             ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    OFF)

# ------------------------
# Compiler flags by vendor
# ------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops -Wall")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w1 -xHost")
  if(ENABLE_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -parallel")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Minform=inform -Mstandard -Wall")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Mbounds -Mchkptr")
endif()

# ------------------------
# Testing & Math config
# ------------------------
include(ConfigTesting)

set(BLAS_LANG "CXX")
set(LAPACK_LANG "CXX")
set(MKL_COMPILER_BINDINGS ${CMAKE_CXX_COMPILER_ID})
include(ConfigMath)

if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

# ------------------------
# Target
# ------------------------
add_executable(FQSolver "")
add_subdirectory(src)

# ------------------------
# OpenMP (portable)
# ------------------------
if(ENABLE_OMP)
  if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    find_path(HOMEBREW_LIBOMP_INCLUDE_DIR
      NAMES omp.h
      PATHS
        /opt/homebrew/opt/libomp/include
        /usr/local/opt/libomp/include
    )
    find_library(HOMEBREW_LIBOMP_LIBRARY
      NAMES omp
      PATHS
        /opt/homebrew/opt/libomp/lib
        /usr/local/opt/libomp/lib
    )

    if(HOMEBREW_LIBOMP_INCLUDE_DIR AND HOMEBREW_LIBOMP_LIBRARY)
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
      set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
      set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP_LIBRARY}" CACHE FILEPATH "" FORCE)
    endif()
  endif()

  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(FQSolver PRIVATE OpenMP::OpenMP_CXX)
    if(HOMEBREW_LIBOMP_INCLUDE_DIR)
      target_include_directories(FQSolver PRIVATE ${HOMEBREW_LIBOMP_INCLUDE_DIR})
    endif()
  else()
    message(FATAL_ERROR "ENABLE_OMP=ON but OpenMP not found")
  endif()
else()
  # Create a private shim so #include <omp.h> works when OpenMP is OFF.
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/omp_stub)
  file(WRITE ${CMAKE_BINARY_DIR}/omp_stub/omp.h
"// Auto-generated shim used when ENABLE_OMP=OFF
#ifndef FQSolver_OMP_STUB_H
#define FQSolver_OMP_STUB_H
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num()  { return 0; }
#endif // FQSolver_OMP_STUB_H
")
  # Prepend shim include path so it wins only when OpenMP is OFF.
  target_include_directories(FQSolver PRIVATE ${CMAKE_BINARY_DIR}/omp_stub)
endif()

# ------------------------
# Link BLAS/LAPACK (and OpenMP if ON above)
# ------------------------
target_link_libraries(FQSolver
  PRIVATE
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

# ------------------------
# Tests assets
# ------------------------
file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})

if(EXISTS ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest)
  file(COPY ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest DESTINATION ${CMAKE_BINARY_DIR}/tests)
endif()
