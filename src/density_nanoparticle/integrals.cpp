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
// Donor–acceptor integrals:
// - int_coulomb: screened Coulomb interaction between reduced densities
// - int_overlap: optional overlap contribution (scaled by -omega_0)
//
// Uses OpenMP for the double loop when available.
void Integrals::acceptor_donor(const Target &target, const Density &acceptor, const Density &donor)
{
  int n_acc = acceptor.n_points_reduced;
  int n_don = donor.n_points_reduced;

  // Direct references (no deep copies). Needed for multithreaded performance.
  const auto &rho_acc = acceptor.rho_reduced;
  const auto &rho_don = donor.rho_reduced;
  const auto &xyz_acc = acceptor.xyz;
  const auto &xyz_don = donor.xyz;

  double int_coulomb = 0.0;
  double int_overlap = 0.0;

  bool calc_overlap = target.calc_overlap_int;

  // QMscrnFact is a constant defined in Parameters
  const double inv_QMscrnFact = 1.0 / Parameters::QMscrnFact;

  // For parallel computation if OMP is ON
  #pragma omp parallel for reduction(+ : int_coulomb, int_overlap) schedule(static)
  for (int i = 0; i < n_acc; ++i)
  {
    for (int j = 0; j < n_don; ++j)
    {
      if (calc_overlap && i == j)
      {
        int_overlap += rho_acc[i] * rho_don[j];
      }

      double dx = xyz_acc[i][0] - xyz_don[j][0];
      double dy = xyz_acc[i][1] - xyz_don[j][1];
      double dz = xyz_acc[i][2] - xyz_don[j][2];

      double dist2 = dx * dx + dy * dy + dz * dz;
      double dist = std::sqrt(dist2);

      if (dist <= 1.0e-14)
        continue;

      // Screening factor to avoid singularities
      double invdist = 1.0 / dist;
      double sf = dist * inv_QMscrnFact;
      double screen_pot = std::erf(sf);

      int_coulomb += rho_acc[i] * rho_don[j] * invdist * screen_pot;
    }
  }

  coulomb_acceptor_donor = int_coulomb;
  if (calc_overlap)
    overlap_acceptor_donor = -target.omega_0 * int_overlap;
}
//----------------------------------------------------------------------
// Acceptor–nanoparticle coupling:
// - charges-only model: uses mm_q (real, imag)
// - charges+dipoles model: adds dipolar field contribution via mm_mu
//
// Both branches compute real/imag parts and store them in
// overlap_acceptor_nanoparticle[0] / [1].
void Integrals::acceptor_np(const Target &target, const Density &acceptor, const Nanoparticle &np)
{
  const int n_acc = acceptor.n_points_reduced;
  const int n_np = np.natoms;

  const auto &rho_acc = acceptor.rho_reduced;
  const auto &xyz_acc = acceptor.xyz;
  const auto &xyz_np = np.xyz;

  const double inv_QMscrnFact = 1.0 / Parameters::QMscrnFact;
  const double sqrt_pi = Parameters::sqrt_pi;

  if (np.charges)
  {
    const auto &mm_q = np.q; // vector<array<double,2>>

    double acceptor_np_int_re_q = 0.0;
    double acceptor_np_int_im_q = 0.0;

    // Parallel computation of the integrals
    #pragma omp parallel for collapse(2) schedule(static) default(none) \
            shared(n_acc, n_np, xyz_acc, xyz_np, rho_acc, mm_q, inv_QMscrnFact) \
            reduction(+ : acceptor_np_int_re_q, acceptor_np_int_im_q)
    for (int i = 0; i < n_acc; ++i)
    {
      for (int j = 0; j < n_np; ++j)
      {

        const double dx = xyz_acc[i][0] - xyz_np[j][0];
        const double dy = xyz_acc[i][1] - xyz_np[j][1];
        const double dz = xyz_acc[i][2] - xyz_np[j][2];

        const double dist2 = dx * dx + dy * dy + dz * dz;
        const double dist = std::sqrt(dist2);

        if (dist <= 1.0e-14)
          continue;

        const double invdist = 1.0 / dist;
        const double sf = dist * inv_QMscrnFact;
        const double screen_pot = std::erf(sf);

        // Change sign: ADF prints densities with opposite sign
        acceptor_np_int_re_q += -rho_acc[i] * mm_q[j][0] * invdist * screen_pot;
        acceptor_np_int_im_q += -rho_acc[i] * mm_q[j][1] * invdist * screen_pot;
      }
    }

    overlap_acceptor_nanoparticle[0] = acceptor_np_int_re_q;
    overlap_acceptor_nanoparticle[1] = acceptor_np_int_im_q;
  }
  else if (np.charges_and_dipoles)
  {
    const auto &mm_q = np.q;   // vector<array<double,2>>
    const auto &mm_mu = np.mu; // vector<array<double,6>>

    double acceptor_np_int_re_q = 0.0;
    double acceptor_np_int_im_q = 0.0;
    double acceptor_np_int_re_mu = 0.0;
    double acceptor_np_int_im_mu = 0.0;

    // Parallel computation of the integrals
    #pragma omp parallel for collapse(2) schedule(static) default(none) \
            shared(n_acc, n_np, xyz_acc, xyz_np, rho_acc, mm_q, inv_QMscrnFact, sqrt_pi, mm_mu) \
            reduction(+:acceptor_np_int_re_q, acceptor_np_int_im_q, \
                        acceptor_np_int_re_mu, acceptor_np_int_im_mu)
    for (int i = 0; i < n_acc; ++i)
    {
      for (int j = 0; j < n_np; ++j)
      {
        const double dx = xyz_acc[i][0] - xyz_np[j][0];
        const double dy = xyz_acc[i][1] - xyz_np[j][1];
        const double dz = xyz_acc[i][2] - xyz_np[j][2];

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
        acceptor_np_int_re_q += -rho_acc[i] * mm_q[j][0] * invdist * screen_pot;
        acceptor_np_int_im_q += -rho_acc[i] * mm_q[j][1] * invdist * screen_pot;

        acceptor_np_int_re_mu += (-rho_acc[i] * mm_mu[j][0] * dx * (invdist*invdist*invdist) * screen_fld
                                  -rho_acc[i] * mm_mu[j][1] * dy * (invdist*invdist*invdist) * screen_fld
                                  -rho_acc[i] * mm_mu[j][2] * dz * (invdist*invdist*invdist) * screen_fld);
                                  
        acceptor_np_int_im_mu += (-rho_acc[i] * mm_mu[j][3] * dx * (invdist*invdist*invdist) * screen_fld
                                  -rho_acc[i] * mm_mu[j][4] * dy * (invdist*invdist*invdist) * screen_fld
                                  -rho_acc[i] * mm_mu[j][5] * dz * (invdist*invdist*invdist) * screen_fld);

      }
    }

    overlap_acceptor_nanoparticle[0] = acceptor_np_int_re_q + acceptor_np_int_re_mu;
    overlap_acceptor_nanoparticle[1] = acceptor_np_int_im_q + acceptor_np_int_im_mu;
  }
  else
  {
    throw std::runtime_error(
        "Nanoparticle model not recognized. Check input file: " + target.nanoparticle_input_file);
  }
}
//----------------------------------------------------------------------
