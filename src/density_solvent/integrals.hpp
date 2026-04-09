#ifndef INTEGRALS_HPP
#define INTEGRALS_HPP

#include "target.hpp"
#include "density.hpp"
#include "solvent.hpp"

//----------------------------------------------------------------------
/// @class Integrals
/// @brief Provides routines to compute potential and field at solvent sites from solute density,
/// and to compute FQ solvent charges from solute density.
///
/// @details
/// The Integrals class stores computed values of electronic couplings
/// (Coulomb and overlap) between donor/acceptor densities and
/// acceptor/nanoparticle systems. It offers functions to evaluate
/// these integrals based on density cubes or nanoparticle data,
/// and makes the results available through public data members.
class Integrals
{
public:
  // ---- Computed values ----
  double coulomb_acceptor_donor = 0.0;

  // ---- API ----

  //  /// @brief Computes solute-solvent potential.
  //  /// @param target        Target object containing calculation options.
  //  /// @param cube_acceptor Density of the acceptor (cube grid).
  // /// @param np            Nanoparticle object (atoms, charges, dipoles).
  // /// @details
  // /// Populates @ref overlap_acceptor_nanoparticle with the computed
  // /// real and imaginary parts of the integral.
  // void acceptor_np(const Target &target, const Density &cube_acceptor, const Nanoparticle &np);

  void solute_solvent_pot_field(const Target &target, const Density &solute, const Solvent &solv);
};

#endif // INTEGRALS_HPP
//----------------------------------------------------------------------
