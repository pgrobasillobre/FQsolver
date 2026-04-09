#ifndef DENSITY_HPP
#define DENSITY_HPP

#include "target.hpp"

#include <string>
#include <vector>
#include <array>

//----------------------------------------------------------------------
class Output; // Forward declaration
//----------------------------------------------------------------------
/// @class Density
/// @brief Represents and processes electron density data from a cube file.
///
/// @details
/// Stores cube metadata (headers, grid, voxel vectors), atomic data,
/// and the density grid. Provides routines to read densities from inputs
/// described by @ref Target, integrate the density, and perform basic
/// geometric operations (centers, rotations, etc.).
class Density {
public:

    // ---- Raw metadata from the cube file (descriptive header lines) ----
    std::string str1, str2;

    // ---- Grid and system metadata ----
    int natoms = 0, nx = 0, ny = 0, nz = 0;
    int nelectrons = 0;
    int n_points_reduced = 0;

    // ---- Atomic data (size natoms) ----
    std::vector<int> atomic_number;              ///< Atomic numbers of the atoms
    std::vector<std::string> atomic_label;       ///< Element symbols
    std::vector<double> atomic_charge;           ///< Atomic charges
    std::vector<double> x, y, z;                 ///< Atomic positions

    // ---- Grid origin and voxel vectors ----
    double xmin = 0.0, ymin = 0.0, zmin = 0.0;
    std::array<double, 3> dx{}, dy{}, dz{};      ///< Voxel vectors in each direction

    // ---- Density data ----
    /// 3D density grid: rho[ix][iy][iz] in e/Å^3 (units may depend on source)
    std::vector<std::vector<std::vector<double>>> rho;

    /// Reduced density values and their coordinates (optional, if reduction used)
    std::vector<double> rho_reduced;
    std::vector<std::array<double, 3>> xyz;

    // ---- Derived quantities ----
    double maxdens = 0.0;  ///< Maximum density value in the grid
    double volume = 0.0;   ///< Total grid volume (a.u.^3)
    
    std::array<double, 3> geom_center = {0.0, 0.0, 0.0}; ///< Geometric center of the density
    std::array<double, 3> geom_center_mol = {0.0, 0.0, 0.0}; ///< Geometric center of the molecule

    double integral = 0.0;  ///< Integral of the density over the full grid


    // ---- API ----

    /// @brief Reads cube-like density and metadata from input described by Target.
    /// @param target  Input descriptor providing file names and calculation context.
    /// @param out     Output/logging sink.
    /// @param what_dens Optional role tag (e.g., "Acceptor", "Donor", "Cube").
    void read_density(Target& target, const Output& out, const std::string& what_dens = "");

    /// @brief Integrates the full 3D density grid to obtain total electrons.
    void int_density();

private:

    /// @brief Maps atomic number to the element symbol.
    /// @param Z Atomic number.
    /// @return Element symbol (e.g., "H", "He", "C").
    std::string map_atomic_number_to_label(int Z) const;

};

#endif // DENSITY_HPP
//----------------------------------------------------------------------
