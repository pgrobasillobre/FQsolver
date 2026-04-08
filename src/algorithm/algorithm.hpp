#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "target.hpp"
#include "output.hpp"
#include "density.hpp"
#include "integrals.hpp"
#include "nanoparticle.hpp"

//----------------------------------------------------------------------
/// @class Algorithm
/// @brief High-level driver for electronic energy transfer calculations
///
/// @details
/// The Algorithm class orchestrates the workflow of density integration,
/// donor–acceptor coupling, nanoparticle–acceptor interactions, and 
/// nanoparticle-mediated transfer rates. 
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

    /// @brief Computes the electronic energy transfer rate 
    ///        between donor and acceptor.
    /// @param target Target system containing donor and acceptor file names 
    ///        and calculation type.
    void acceptor_donor(Target &target);

    /// @brief Computes the coupling between an acceptor and a nanoparticle.
    /// @param target Target system containing acceptor and nanoparticle file names
    ///        and calculation type.
    void acceptor_np(Target &target);

 
    /// @brief Computes the donor–acceptor electronic energy transfer rate
    ///        mediated by a nanoparticle.
    /// @param target Target system containing donor, acceptor, and nanoparticle file names
    ///        and calculation type.
    void acceptor_np_donor(Target &target);

private:
    Output &out;        ///< Reference to output handler.
    Target &target;     ///< Reference to target system.
    Density cube;       ///< General-purpose density object.
    Density cube_acceptor; ///< Density associated with the acceptor.
    Density cube_donor;    ///< Density associated with the donor.
    Integrals integrals;   ///< Integral evaluator for couplings.
    Nanoparticle np;       ///< Nanoparticle representation.
};

#endif
//----------------------------------------------------------------------

