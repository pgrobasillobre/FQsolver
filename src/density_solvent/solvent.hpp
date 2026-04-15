#ifndef SOLVENT_HPP
#define SOLVENT_HPP

#include "target.hpp"

#include <string>
#include <vector>
#include <array>

class Output; // Forward declaration

//----------------------------------------------------------------------
/// @class Solvent
/// @brief Represents a solvent with atomic positions, potential/field at atomic sites and charges.
///
/// @details
/// This class stores data (atomic positions, potential/field at atomic sites, and charges)
/// for a solvent system and provides methods for reading that data.

class Solvent
{
public:
  int natoms = 0; ///< Number of atoms in the solvent.
  int nmol = 0;   ///< Number of solvent molecules (for FQ charge calculation).
  int ntypes = 0; ///< Number of unique atom types (for FQ charge calculation).

  bool potential = false;           ///< True if only potential calculation is requested.
  bool field = false;               ///< True if only field calculation is requested.
  bool potential_and_field = false; ///< True if potential and field calculation is requested.

  // Arrays
  std::vector<std::string> typeName;        // FQ atom type names
  std::vector<std::array<int, 1>> MolIndex; // Index of the molecule each atom belongs to (for FQ charge calculation)

  std::vector<int> typeIndex; // Atom index to FQ atom type index

  std::vector<double> typeChi;                 // FQ electronegativity by atom type
  std::vector<double> typeEta;                 // FQ hardness by atom type
  std::vector<double> typeRq;                  // Gaussian width for charges probability distribution by atom type (for FQ charge calculation)
  
  std::vector<std::array<double, 1>> solv_pot; // Scalar potential at each atomic site
  std::vector<std::array<double, 3>> solv_fld; // Electric field vector at each atomic site
  std::vector<std::array<double, 3>> xyz;      // XYZ coordinates

  std::vector<std::vector<double>> tempTqq; // Temporary storage for Tqq tensor components (for FQ charge calculation) natoms natoms

  std::vector<std::string> atomic_label; //< Atomic labels (e.g., "C", "O").

  std::string solvent_file_extension; ///< Extension of the solvent geometry file (e.g., ".xyz").

  /// @brief Loads solvent coordinates from the specified file.
  /// @param target Target object providing file paths.
  /// @param out    Output object used for logging/debugging.
  /// @details
  /// This method reads the solvent geometry from the file specified in the Target object.
  void read_solvent(const Target &target, const Output &out);

private:
  /// @brief Checks if the specified file exists.
  /// @param filepath Path to the file to check.
  void check_file_exists(const std::string &filepath) const;

  /// @brief Verifies that the solvent file has the expected extension.
  /// @param filepath Path to the solvent geometry file.
  /// @param target Target object providing expected file extension information.
  void check_solvent_file_extension(const std::string &str, const Target &target);

  /// @brief Selects the appropriate method to read solvent geometry based on file extension.
  /// @param filepath Path to the solvent geometry file.
  void read_solvent_geometry(const std::string &filepath, const Target &target);

  /// @brief Reads solvent geometry from an XYZ file.
  /// @param filepath Path to the XYZ file containing solvent geometry.
  void read_solvent_geometry_xyz(const std::string &filepath);

  /// @brief Reads solvent geometry from a PDB file.
  /// @param filepath Path to the PDB file containing solvent geometry.
  /// @param target Target object providing optional PDB filters.
  void read_solvent_geometry_pdb(const std::string &filepath, const Target &target);

  /// @brief Assigns FQ parameters to solvent atoms based on the specified parametrization.
  /// @param parametrization Name of the FQ parametrization to use (e.g., "giovannini").
  void assign_solvent_parameters(const std::string &parametrization);
};

#endif // SOLVENT_HPP
//----------------------------------------------------------------------
