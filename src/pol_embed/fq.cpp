#include "fq.hpp"
#include "solvent.hpp"

#include <cmath>
#include <stdexcept>

//----------------------------------------------------------------------
// Compute FQ charges from the potential and solvent parameters.
void FQ::calc_charges(Solvent &solv)
{
    // Compute Tqq matrix (charge-charge interaction tensor).
    calc_tqq(solv);
}
//  ----------------------------------------------------------------------
// Compute the Tqq matrix (charge-charge interaction tensor).
void FQ::calc_tqq(Solvent &solv)
{

    // Initial checks
    if (static_cast<int>(solv.typeIndex.size()) != solv.natoms ||
        static_cast<int>(solv.xyz.size()) != solv.natoms ||
        static_cast<int>(solv.tempTqq.size()) != solv.natoms)
    {
        throw std::runtime_error("calc_tqq: Solvent data are not ready for FQ Tqq calculation.");
    }

    // Select the appropriate kernel for Tqq calculation based on the input file.
    solv.tempTqq.assign(solv.natoms, std::vector<double>(solv.natoms, 0.0));

    // Diagonal elements: Tqq_ii = eta_i
    for (int i = 0; i < solv.natoms; ++i)
    {
        const int type_i = solv.typeIndex[i];
        if (type_i < 0 ||
            type_i >= static_cast<int>(solv.typeEta.size()) ||
            type_i >= static_cast<int>(solv.typeRq.size()))
        {
            throw std::runtime_error("Invalid FQ atom type index while building Tqq.");
        }

        solv.tempTqq[i][i] = solv.typeEta[type_i];
    }

    // Off-diagonal elements: Tqq_ij dependent on the selected interaction kernel

    if (solv.fq_kernel == "gaus")
    {
        for (int i = 0; i < solv.natoms; ++i)
        {
            const int type_i = solv.typeIndex[i];
            const double rq_i2 = solv.typeRq[type_i] * solv.typeRq[type_i];

            for (int j = 0; j < i; ++j)
            {
                const int type_j = solv.typeIndex[j];
                if (type_j < 0 ||
                    type_j >= static_cast<int>(solv.typeRq.size()))
                {
                    throw std::runtime_error("Invalid FQ atom type index while building Tqq.");
                }

                const double dx = solv.xyz[i][0] - solv.xyz[j][0];
                const double dy = solv.xyz[i][1] - solv.xyz[j][1];
                const double dz = solv.xyz[i][2] - solv.xyz[j][2];
                const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (dist <= 1.0e-14)
                {
                    throw std::runtime_error("Cannot build FQ Tqq matrix: two solvent atoms have the same coordinates.");
                }

                const double rq_j2 = solv.typeRq[type_j] * solv.typeRq[type_j];
                solv.tempTqq[i][j] = std::erf(dist / std::sqrt(rq_i2 + rq_j2)) / dist;
                solv.tempTqq[j][i] = solv.tempTqq[i][j];
            }
        }
    }
    else if (solv.fq_kernel == "coul")
    {
        throw std::runtime_error("FQ kernel \"coul\" is not implemented yet. Only \"gaus\" is currently available.");
    }
    else if (solv.fq_kernel == "ohno")
    {
        throw std::runtime_error("FQ kernel \"ohno\" is not implemented yet. Only \"gaus\" is currently available.");
    }
    else
    {
        throw std::runtime_error("Unsupported FQ kernel \"" + solv.fq_kernel + "\". Only \"gaus\" is currently available.");
    }
}
