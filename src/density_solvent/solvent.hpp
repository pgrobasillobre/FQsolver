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

  int natoms = 0;  ///< Number of atoms in the nanoparticle.

  bool charges = false;              ///< True if nanoparticle contains only charges.
  bool charges_and_dipoles = false;  ///< True if nanoparticle includes dipoles.

  std::string nanoparticle_model;    ///< Input model name (used to determine input structure).

  std::array<double, 3> geom_center; ///< Geometric center of the nanoparticle.

  std::vector<std::array<double, 2>> q;    // Charges with real + imaginary part
  std::vector<std::array<double, 6>> mu;   // Dipoles with 3 components each for real + imaginary part

  std::vector<std::array<double, 3>> xyz;  // XYZ coordinates

  /// @brief Loads nanoparticle coordinates and multipole data from input.
  /// @param target Target object providing file paths and rotation options.
  /// @param out    Output object used for logging/debugging.
  /// @details
  /// Automatically detects the model type and loads charge and/or dipole data accordingly.
  /// Applies rotation if specified in the target configuration.
  void read_nanoparticle(const Target &target, const Output &out);

};

#endif // NANOPARTICCLE_HPP
//----------------------------------------------------------------------
