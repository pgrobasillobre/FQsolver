#include "fq.hpp"
#include "output.hpp"
#include "solvent.hpp"
#include "target.hpp"
#include "integrals.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>

//----------------------------------------------------------------------
// Compute FQ charges from the potential and solvent parameters.
void FQ::calc_charges(Solvent &solv, const Target &target, const Output &out)
{
    // Compute Tqq matrix (charge-charge interaction tensor).
    calc_tqq(solv, target, out);

    // Calculate FQ LHS DMatrix
    calc_dmat(solv, target, out);

    // Calculate FQ RHS vector
    calc_rhs(solv, target, out);
}
//  ----------------------------------------------------------------------
// Compute the Tqq matrix (charge-charge interaction tensor).
void FQ::calc_tqq(Solvent &solv, const Target &target, const Output &out)
{

    // Initial checks
    if (static_cast<int>(solv.typeIndex.size()) != solv.natoms ||
        static_cast<int>(solv.xyz.size()) != solv.natoms)
    {
        throw std::runtime_error("calc_tqq: Solvent data are not ready for FQ Tqq calculation.");
    }

    // Select the appropriate kernel for Tqq calculation based on the input file.
    tempTqq.assign(solv.natoms, std::vector<double>(solv.natoms, 0.0));

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

        tempTqq[i][i] = solv.typeEta[type_i];
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
                tempTqq[i][j] = std::erf(dist / std::sqrt(rq_i2 + rq_j2)) / dist;
                tempTqq[j][i] = tempTqq[i][j];
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

    if (target.debug >= 2)
    {
        out.print_matrix("Tqq Matrix", tempTqq);
    }
}
// ----------------------------------------------------------------------
// Calculation of FQ LHS matrix: DMatrix
void FQ::calc_dmat(Solvent &solv, const Target &target, const Output &out)
{

    // Initial checks.
    // 1. Check molecule indices are ready for DMatrix construction.
    if (static_cast<int>(solv.atomsToIndex.size()) != solv.natoms)
    {
        throw std::runtime_error("calc_dmat: Solvent molecule indices are not ready for FQ DMatrix calculation.");
    }

    // Check order of DMatrix
    if (target.what == "fq")
    {
        norder = solv.natoms + solv.nmol; // number of solvent atoms + number of solvent molecules (for charge conservation constraint)
    }
    else if (target.what == "fqfmu")
    {
        throw std::runtime_error("FQFmu charge calculation is not implemented yet. Only \"fq\" is currently available.");
        // norder = 4*solv.natoms + solv.nmol; // 4 times the number of solvent atoms + number of solvent molecules (for charge conservation constraint)
    }
    else
    {
        throw std::runtime_error("Invalid 'what' keyword for FQ charge calculation. Expected 'fq' or 'fq+'.");
    }

    // 2. Check Tqq matrix is ready for DMatrix construction.
    if (static_cast<int>(tempTqq.size()) != solv.natoms ||
        static_cast<int>(tempTqq[0].size()) != solv.natoms)
    {
        throw std::runtime_error("calc_dmat: Tqq matrix is not ready for FQ DMatrix calculation.");
    }

    //
    // DMatrix is the packed lower-triangular form of the full FQ linear system:
    //
    //        | Tqq   C^T |
    //    D = |           |
    //        | C      0  |
    //
    // where:
    //   - Tqq is natoms x natoms
    //   - C is nmol x natoms
    //   - C(iMol, iAtom) = 1 if atom iAtom belongs to molecule iMol
    //   - the 0 block is nmol x nmol
    //
    // Full dense shape:
    //
    //              q0   q1   q2   ...  qN-1            l0   l1   ...  lM-1
    //        -----------------------------------------------------------------
    //   q0  |      Tqq block                         |   C^T block
    //   q1  |                                        |
    //   q2  |                                        |
    //  ...  |                                        |
    // qN-1  |                                        |
    //        -----------------------------------------------------------------
    //   l0  |      C block                           |   0 block
    //   l1  |                                        |
    //  ...  |                                        |
    // lM-1  |                                        |
    //        -----------------------------------------------------------------
    //
    // Packed lower-triangular indexing uses zero-based C++ indices:
    //
    //   packed_index(row, col) = row * (row + 1) / 2 + col,  with col <= row
    //
    // Therefore:
    //   - Tqq(i,j) goes to DMatrix[packed_index(i,j)]
    //   - constraint 1s go to DMatrix[packed_index(natoms + iMol, iAtom)]
    //
    //
    // Practical example for 2 water molecules (natoms = 6 and nmol = 2, so norder = 8).
    //
    // Each molecule has 3 atoms, and each atom is represented by a charge variable (q0 to q5), and each molecule has a Lambda multiplier for the charge conservation constraint (l0 and l1):
    //
    //   molecule l0: q0, q1, q2
    //   molecule l1: q3, q4, q5
    //
    // Full dense DMatrix:
    //
    //              q0        q1        q2        q3        q4        q5        l0        l1
    //        ---------------------------------------------------------------------------------
    //   q0  |   Tqq(0,0)  Tqq(0,1)  Tqq(0,2)  Tqq(0,3)  Tqq(0,4)  Tqq(0,5)  1         0       |
    //   q1  |   Tqq(1,0)  Tqq(1,1)  Tqq(1,2)  Tqq(1,3)  Tqq(1,4)  Tqq(1,5)  1         0       |
    //   q2  |   Tqq(2,0)  Tqq(2,1)  Tqq(2,2)  Tqq(2,3)  Tqq(2,4)  Tqq(2,5)  1         0       |
    //   q3  |   Tqq(3,0)  Tqq(3,1)  Tqq(3,2)  Tqq(3,3)  Tqq(3,4)  Tqq(3,5)  0         1       |
    //   q4  |   Tqq(4,0)  Tqq(4,1)  Tqq(4,2)  Tqq(4,3)  Tqq(4,4)  Tqq(4,5)  0         1       |
    //   q5  |   Tqq(5,0)  Tqq(5,1)  Tqq(5,2)  Tqq(5,3)  Tqq(5,4)  Tqq(5,5)  0         1       |
    //   l0  |   1         1         1         0         0         0         0         0       |
    //   l1  |   0         0         0         1         1         1         0         0       |
    //        ---------------------------------------------------------------------------------
    //
    // Dmat stores only the lower triangle of this symmetric matrix.
    // The packed index for zero-based C++ indices is:
    //
    //   packed_index(row, col) = row * (row + 1) / 2 + col, with col <= row
    //
    // Therefore the constraint entries stored explicitly are:
    //
    //   Dmat[packed_index(6,0)] = D(l0,q0) = 1
    //   Dmat[packed_index(6,1)] = D(l0,q1) = 1
    //   Dmat[packed_index(6,2)] = D(l0,q2) = 1
    //
    //   Dmat[packed_index(7,3)] = D(l1,q3) = 1
    //   Dmat[packed_index(7,4)] = D(l1,q4) = 1
    //   Dmat[packed_index(7,5)] = D(l1,q5) = 1

    // 1. Compute nlength to store the D matrix as lower triangular
    nlength = (norder * (norder + 1)) / 2;
    Dmat.assign(nlength, 0.0);

    // 2. Introduce Tqq
    for (int i = 0; i < solv.natoms; ++i)
    {
        const int ii = (i * (i + 1)) / 2;
        for (int j = 0; j <= i; ++j)
        {
            Dmat[ii + j] = tempTqq[i][j];
        }
    }

    // 3. Eliminate Tqq temporary storage to free memory, as it is no longer needed after being copied to Dmat.
    tempTqq.clear();         // Clear the contents of tempTqq to free memory.
    tempTqq.shrink_to_fit(); // Reduce capacity to fit the new size (which is zero after clear).

    // 4. Introduce Lambda multipliers for charge conservation constraint at each solvent molecule.
    for (int i = 0; i < solv.natoms; ++i)
    {
        const int molecule_index = solv.atomsToIndex[i];
        if (molecule_index < 0 || molecule_index >= solv.nmol)
        {
            throw std::runtime_error("calc_dmat: Invalid compact molecule index.");
        }

        const int lambda_index = solv.natoms + molecule_index;
        const int packed_index = (lambda_index * (lambda_index + 1)) / 2 + i;
        Dmat[packed_index] = 1.0;
    }

    // Print Dmatrix reconstructing the full dense form
    if (target.debug >= 2)
    {
        out.print_matrix("FQ DMatrix", Dmat, norder, "L");
    }
}
// ----------------------------------------------------------------------
// Calculation of FQ RHS vector
void FQ::calc_rhs(Solvent &solv, const Target &target, const Output &out)
{
    // Initial checks
    if (static_cast<int>(solv.solv_pot.size()) != solv.natoms)
    {
        throw std::runtime_error("calc_rhs: Solvent potential data are not ready for FQ RHS calculation.");
    }

    // debugpgi
    //  Introduce potential from ADF to debug
    // solv.solv_pot[0][0] = -0.0000267235957726;
    // solv.solv_pot[1][0] = -0.0000281716086517;
    // solv.solv_pot[2][0] = -0.0000291111489726;
    // solv.solv_pot[3][0] = -0.0000067889794474;
    // solv.solv_pot[4][0] = -0.0000049062211308;
    // solv.solv_pot[5][0] = -0.0000057162322666;
    // enddebugpgi

    rhs.assign(norder, 0.0);
    // RHS atomic block: -potential at each atomic site - FQ electronegativity by atom type
    for (int i = 0; i < solv.natoms; ++i)
    {
        rhs[i] = -solv.solv_pot[i][0] - solv.typeChi[solv.typeIndex[i]];
    }

    // RHS constraint block: molecular charge for each solvent molecule
    for (int i = solv.natoms; i < norder; ++i)
    {
        rhs[i] = solv.MolCharge;
    }

    // print rhs for debuggin
    for (int i = 0; i < norder; ++i)
    {
        std::cout << rhs[i] << std::endl;
    }

    // Print RHS vector for debugging.
    if (target.debug >= 2)
    {
        out.print_matrix_rhs("FQ RHS Vector", solv, rhs);
    }
}
