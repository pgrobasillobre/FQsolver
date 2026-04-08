#ifndef NANOPARTICLE_HPP
#define NANOPARTICLE_HPP

#include "target.hpp"

#include <string>
#include <vector>
#include <array>


class Output; // Forward declaration

//----------------------------------------------------------------------
/// @class Nanoparticle
/// @brief Represents a nanoparticle with atomic positions, charges, and dipoles.
///
/// @details
/// This class stores coordinate and multipole data for a nanoparticle system
/// (charges, dipoles) and provides methods for reading and rotating that data.
/// It supports multiple input formats based on the model defined in Target.
/// Charges and dipoles are stored as real and imaginary components when needed.
class Nanoparticle
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

private:

  /// @brief Applies rotation to nanoparticle coordinates and dipoles (if present).
  /// @param target Target object with rotation settings.
  /// @param out    Output object used for debug output.
  void rotate_np_coords_and_dipoles(const Target &target, const Output &out);

  /// @brief Rotates all dipoles in the nanoparticle.
  /// @param angle Rotation angle (radians).
  /// @param axis  Axis of rotation ("x", "y", or "z").
  void rotate_np_dipoles(const double angle, const std::string& axis);

  /// @brief Rotates all atomic coordinates.
  /// @param angle Rotation angle (radians).
  /// @param axis  Axis of rotation ("x", "y", or "z").
  /// @return Rotated coordinates.
  std::vector<std::array<double, 3>> rotate_np_coords(const double angle, const std::string& axis) const;
};

#endif // NANOPARTICCLE_HPP
//----------------------------------------------------------------------
