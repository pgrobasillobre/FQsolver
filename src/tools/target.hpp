#ifndef TARGET_HPP
#define TARGET_HPP

#include "enum.hpp"

#include <string>
#include <array>

//----------------------------------------------------------------------
/// @struct Target
/// @brief Represents all user-defined input options and flags controlling the calculation.
///
/// This struct aggregates all input parameters needed by the FQSolver application.
struct Target
{
    // Density integration
    std::string input_filename;                 ///< Name of the input file used
    std::string density_file_integration;       ///< Full path to the density cube for integration
    std::string density_file_integration_input; ///< Filename as given in the input file

    // What to calculate
    std::string what; ///< What to calculate ("pot", "pot_fld")

    // Solvent -  geometry
    bool is_solvent_present = false; ///< Whether a solvent file is provided.

    std::string solvent_file;       ///< Full path to the solvent geometry file
    std::string solvent_input_file; ///< Solvent file name from input file

    // Solute - density

    bool is_cube_density_present = false;

    std::string solute_density_file;       ///< Full path to the density cube file
    std::string solute_density_input_file; ///< Density cube name from input file

    // Target + other options
    TargetMode mode = TargetMode::None; ///< Selected calculation target (main mode of operation)

    bool integrate_density = false; ///< Whether to integrate a single cube density

    bool is_what_present = false; ///< Whether the "what" keyword is present in the input file

    bool is_cutoff_present = false; ///< Flag to apply distance cutoff for grid reduction

    double cutoff = 0.0; ///< Cutoff energy (Hartree) for density grid reduction

    bool is_debug_present; ///< Debug mode enabled
    int debug = 0;

    int n_threads_OMP = 1; ///< Number of threads for OpenMP parallelism
};

#endif // TARGET_HPP
//----------------------------------------------------------------------
