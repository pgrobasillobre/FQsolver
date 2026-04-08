#ifndef INTEGRALS_HPP
#define INTEGRALS_HPP

#include "target.hpp"
#include "density.hpp"
#include "nanoparticle.hpp"

//----------------------------------------------------------------------
/// @class Integrals
/// @brief Provides routines to compute donor–acceptor and acceptor–nanoparticle integrals.
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
  double overlap_acceptor_donor = 0.0;

  /// Overlap between acceptor density and nanoparticle:
  /// [0] = real part, [1] = imaginary part
  std::array<double, 2> overlap_acceptor_nanoparticle = {0.0, 0.0};



  // ---- API ----

  /// @brief Computes donor–acceptor Coulomb and overlap integrals.
  /// @param target         Target object containing calculation options.
  /// @param cube_acceptor  Density of the acceptor (cube grid).
  /// @param cube_donor     Density of the donor (cube grid).
  /// @details
  /// Populates @ref coulomb_acceptor_donor and @ref overlap_acceptor_donor.
  void acceptor_donor(const Target &target, const Density &cube_acceptor, const Density &cube_donor);

  /// @brief Computes acceptor–nanoparticle coupling integrals.
  /// @param target        Target object containing calculation options.
  /// @param cube_acceptor Density of the acceptor (cube grid).
  /// @param np            Nanoparticle object (atoms, charges, dipoles).
  /// @details
  /// Populates @ref overlap_acceptor_nanoparticle with the computed
  /// real and imaginary parts of the integral.
  void acceptor_np(const Target &target, const Density &cube_acceptor, const Nanoparticle &np);
};

#endif // INTEGRALS_HPP
//----------------------------------------------------------------------
