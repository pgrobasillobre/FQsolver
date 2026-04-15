# src/CMakeLists.txt

# Add all relevant C++ source files to the FQSolver target
target_sources(FQSolver
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/algorithm/algorithm.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/density_solvent/integrals.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/density_solvent/solvent.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/density_solvent/density.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/input_output/output.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/input_output/input.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/string_manipulation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/timer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/parameters.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pol_embed/fq.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# Make headers in src/ and subfolders accessible
target_include_directories(FQSolver
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/input_output
        ${CMAKE_CURRENT_SOURCE_DIR}/tools
        ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
        ${CMAKE_CURRENT_SOURCE_DIR}/density_solvent
        ${CMAKE_CURRENT_SOURCE_DIR}/pol_embed
)

