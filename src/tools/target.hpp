#ifndef TARGET_HPP
#define TARGET_HPP

#include "enum.hpp"

#include <string>
#include <array>

//----------------------------------------------------------------------
/// @struct Target
/// @brief Represents all user-defined input options and flags controlling the calculation.
///
/// This struct aggregates all input parameters needed by the FQSolver application,
/// including file paths, calculation modes, dipole alignment settings, rotation axes,
/// and integration/interaction flags.
struct Target
{
    TargetMode mode = TargetMode::None; ///< Selected calculation target mode

    // Density integration
    std::string input_filename; ///< Name of the input file used
    std::string density_file_integration;       ///< Full path to the density cube for integration
    std::string density_file_integration_input; ///< Filename as given in the input file

    // Acceptor
    bool is_acceptor_density_present = false;
    bool is_acceptor_transition_dipole_present = false;
    bool is_acceptor_transition_dipole_align_present = false;
    bool rotate_acceptor = false;

    std::string acceptor_density_file;       ///< File for acceptor density (full path)
    std::string acceptor_density_input_file; ///< Acceptor cube name from input file

    std::array<double, 3> acceptor_transdip = {0.0, 0.0, 0.0};      ///< Acceptor transition dipole (input)
    std::array<double, 3> acceptor_transdip_rot = {0.0, 0.0, 0.0};   ///< Rotated acceptor dipole
    std::array<double, 3> acceptor_ref_vector = {0.0, 0.0, 0.0};     ///< Alignment reference vector

    double acceptor_density_rotation_angle = 0.0;              ///< Requested alignment rotation angle (rad)
    double acceptor_density_rotation_angle_check = 0.0;        ///< Actual applied rotation angle (rad)

    // Donor
    bool is_donor_density_present = false;
    bool is_donor_transition_dipole_present = false;
    bool is_donor_transition_dipole_align_present = false;
    bool rotate_donor = false;

    std::string donor_density_file;       ///< Full path to the donor density cub
    std::string donor_density_input_file; ///< Donor cube name from input file

    std::array<double, 3> donor_transdip = {0.0, 0.0, 0.0};     ///< Donor transition dipole (input)
    std::array<double, 3> donor_transdip_rot = {0.0, 0.0, 0.0}; ///< Rotated donor dipole
    std::array<double, 3> donor_ref_vector = {0.0, 0.0, 0.0};   ///< Alignment reference vector

    double donor_density_rotation_angle = 0.0;           ///< Requested alignment rotation angle (rad)
    double donor_density_rotation_angle_check = 0.0;     ///< Actual applied rotation angle (rad)

    // Nanoparticle
    bool is_nanoparticle_present = false;
    bool rotate_nanoparticle = false;

    std::string nanoparticle_file;       ///< Full path to the nanoparticle definition file
    std::string nanoparticle_input_file; ///< File name as given in the input

    // Target + other options
    TargetMode target_mode = TargetMode::None; ///< Selected calculation target (main mode of operation)

    std::string rotation_axys;          ///< Axis around which rotations are performed
    bool is_rotation_axys_present = false;

    bool integrate_density = false;     ///< Whether to integrate a single cube density

    bool calc_overlap_int = false;      ///< Whether to calculate overlap integral

    bool is_cutoff_present = false;     ///< Flag to apply distance cutoff for grid reduction

    double cutoff = 0.0;;               ///< Cutoff radius (Bohr)

    bool is_omega_0_present = false;
    double omega_0 = 0.0;               ///< Energy difference between donor and acceptor (Hartree)

    bool is_spectral_overlap_present = false;
    double spectral_overlap = 0.0;      ///< Spectral overlap integral (unitless)

    bool is_debug_present;              ///< Debug mode enabled
    int debug = 0;

    int n_threads_OMP = 1;              ///< Number of threads for OpenMP parallelism
};

#endif // TARGET_HPP
//----------------------------------------------------------------------
