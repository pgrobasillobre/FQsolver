#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "target.hpp"
#include "output.hpp"
#include "density.hpp"
#include "integrals.hpp"
#include "solvent.hpp"
#include "fq.hpp"

//----------------------------------------------------------------------
/// @class Algorithm
/// @brief High-level driver for potential, field, and FQ calculations.
///
/// @details
/// The Algorithm class orchestrates the workflow of FQSolver,
/// including reading input data, performing density integrations,
/// computing potentials, fields, charges, and outputting results.
class Algorithm
{
public:
    /// @brief Construct a new Algorithm object.
    /// @param out    Output handler for logging and reporting.
    /// @param target Target system to be used in calculations.
    Algorithm(Output &out, Target &target);

    /// @brief Integrates the electron density from a cube file.
    /// @param target Target system providing the cube file name and calculation type.
    void integrate_density(Target &target);

    ///
    /// @brief Compute electronic energy transfer rate between donor and acceptor.
    ///

    /// @brief Computes the potential/field at the solvent coordinates from a solute density.
    /// @param target Target system containing solute and solvent file names and calculation type.
    void solute_dens_solvent_pot_fld(Target &target);

    /// @brief Computes FQ charges from the potential and solvent parameters.
    /// @param target Target system containing necessary input for FQ charge calculation.
    void compute_fq_charges(Target &target);

private:
    Output &out;         ///< Reference to output handler.
    Target &target;      ///< Reference to target system.
    Density cube;        ///< General-purpose density object.
    Density cube_solute; ///< Density associated with the solute
    Integrals integrals; ///< Integral evaluator for couplings.
    Solvent solv;        ///< Solvent representation.
    FQ fq;               ///< FQ charge calculator.
};

#endif
//----------------------------------------------------------------------
