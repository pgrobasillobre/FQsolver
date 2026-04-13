#include "integrals.hpp"
#include "target.hpp"
#include "density.hpp"
#include "parameters.hpp"

#include <cmath>
#include <omp.h>
#include <ostream>
#include <stdexcept>
#include <iostream>

//----------------------------------------------------------------------
// Solute-solvent potential and field:
//
void Integrals::solute_solvent_pot_fld(const Target &target, const Density &solute, const Solvent &solv)
{

  // This function computes the potential/field at the solvent coordinates due to the solute density.
  // This involves evaluating integrals of the form:
  //
  // V(r) = ∫ ρ(r') / |r - r'| dr'
  // E(r) = -∇V(r)
  //
  // where ρ(r') is the solute electron density, and r are the coordinates of the solvent atoms.
  //

  // Retrieve the number of points in the solute density and the number of solvent atoms
  const int n_solute = solute.n_points_reduced;
  const int n_solvent = solv.natoms;

  // Retrieve the reduced solute density values and their coordinates, and the solvent coordinates
  const auto &rho_solute = solute.rho_reduced;
  const auto &xyz_solute = solute.xyz;
  const auto &xyz_solvent = solv.xyz;

  // Precompute constants for the screened Coulomb potential
  const double inv_QMscrnFact = 1.0 / Parameters::QMscrnFact;
  const double sqrt_pi = Parameters::sqrt_pi;

  if (target.what == "potential")
  {
    solv_pot.assign(n_solvent, {0.0});
    solv_fld.clear();

    // Compute only the potential at the solvent coordinates
     // #pragma omp parallel for schedule(static) default(none)                                       \
    shared(n_solute, n_solvent, xyz_solute, xyz_solvent, rho_solute, inv_QMscrnFact, sqrt_pi) \
    reduction(+ : solv_pot)

    for (int i = 0; i < n_solvent; ++i)
    {
      double pot_i = 0.0;
      for (int j = 0; j < n_solute; ++j)
      {
        const double dx = xyz_solvent[i][0] - xyz_solute[j][0];
        const double dy = xyz_solvent[i][1] - xyz_solute[j][1];
        const double dz = xyz_solvent[i][2] - xyz_solute[j][2];

        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist = std::sqrt(dist2);

        if (dist <= 1.0e-14)
          continue;

        const double invdist = 1.0 / dist;

        const double sf = dist * inv_QMscrnFact;
        const double screen_pot = std::erf(sf);

        // Change sign: ADF prints densities with opposite sign
        pot_i += -rho_solute[j] * invdist * screen_pot;
      }
      solv_pot[i][0] = pot_i;
    }
  }
  else if (target.what == "field")
  {
    solv_fld.assign(n_solvent, {0.0, 0.0, 0.0});
    solv_pot.clear();
    // Compute only the field at the solvent coordinates
     // #pragma omp parallel for schedule(static) default(none)                                       \
    shared(n_solute, n_solvent, xyz_solute, xyz_solvent, rho_solute, inv_QMscrnFact, sqrt_pi) \
    reduction(+ : solv_fld)
    for (int i = 0; i < n_solvent; ++i)
    {
      std::array<double, 3> fld_i = {0.0, 0.0, 0.0};
      for (int j = 0; j < n_solute; ++j)
      {
        const double dx = xyz_solvent[i][0] - xyz_solute[j][0];
        const double dy = xyz_solvent[i][1] - xyz_solute[j][1];
        const double dz = xyz_solvent[i][2] - xyz_solute[j][2];
        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist = std::sqrt(dist2);
        if (dist <= 1.0e-14)
          continue;
        const double invdist = 1.0 / dist;
        const double sf = dist * inv_QMscrnFact;
        const double screen_pot = std::erf(sf);
        const double sf1 = (2.0 * sf / sqrt_pi) * std::exp(-sf * sf);
        const double screen_fld = screen_pot - sf1;
        // Change sign: ADF prints densities with opposite sign
        fld_i[0] += -rho_solute[j] * dx * (invdist * invdist * invdist) * screen_fld;
        fld_i[1] += -rho_solute[j] * dy * (invdist * invdist * invdist) * screen_fld;
        fld_i[2] += -rho_solute[j] * dz * (invdist * invdist * invdist) * screen_fld;
      }
      solv_fld[i][0] = fld_i[0];
      solv_fld[i][1] = fld_i[1];
      solv_fld[i][2] = fld_i[2];
    }
  }
  else if (target.what == "potential+field")
  {

    solv_pot.assign(n_solvent, {0.0});
    solv_fld.assign(n_solvent, {0.0, 0.0, 0.0});

    // Compute both potential and field at the solvent coordinates
     // #pragma omp parallel for schedule(static) default(none)                                       \
    shared(n_solute, n_solvent, xyz_solute, xyz_solvent, rho_solute, inv_QMscrnFact, sqrt_pi) \
    reduction(+ : solv_pot, solv_fld)
    for (int i = 0; i < n_solvent; ++i)
    {
      double pot_i = 0.0;
      std::array<double, 3> fld_i = {0.0, 0.0, 0.0};
      for (int j = 0; j < n_solute; ++j)
      {
        // Check how to perform this substraction !!!!! first solvent or first solute ???????????????
        const double dx = xyz_solvent[i][0] - xyz_solute[j][0];
        const double dy = xyz_solvent[i][1] - xyz_solute[j][1];
        const double dz = xyz_solvent[i][2] - xyz_solute[j][2];

        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist = std::sqrt(dist2);

        if (dist <= 1.0e-14)
          continue;

        const double invdist = 1.0 / dist;
        const double sf = dist * inv_QMscrnFact;
        const double screen_pot = std::erf(sf);
        const double sf1 = (2.0 * sf / sqrt_pi) * std::exp(-sf * sf);
        const double screen_fld = screen_pot - sf1;

        // Change sign: ADF prints densities with opposite sign
        pot_i += -rho_solute[j] * invdist * screen_pot;
        fld_i[0] += -rho_solute[j] * dx * (invdist * invdist * invdist) * screen_fld;
        fld_i[1] += -rho_solute[j] * dy * (invdist * invdist * invdist) * screen_fld;
        fld_i[2] += -rho_solute[j] * dz * (invdist * invdist * invdist) * screen_fld;
      }
      solv_pot[i][0] = pot_i;
      solv_fld[i][0] = fld_i[0];
      solv_fld[i][1] = fld_i[1];
      solv_fld[i][2] = fld_i[2];
    }
    // debugpgi: Print results for debugging
    for (int i = 0; i < n_solvent; ++i)
    {
      std::cout << "For atom #" << i << ": Potential = " << solv_pot[i][0] << ", Field = (" << solv_fld[i][0] << ", " << solv_fld[i][1] << ", " << solv_fld[i][2] << std::endl;
    }
  }
  else
  {
    throw std::runtime_error("Invalid 'what' option for solute-solvent calculation.");
  }
}
//----------------------------------------------------------------------
