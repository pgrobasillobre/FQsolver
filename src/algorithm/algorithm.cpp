#include "algorithm.hpp"
#include "output.hpp"
#include "target.hpp"
#include "density.hpp"
#include "parameters.hpp"
#include "integrals.hpp"
#include "solvent.hpp"

#include <iostream>

//----------------------------------------------------------------------
// Constructor
Algorithm::Algorithm(Output &out, Target &target) : out(out), target(target) {}
//----------------------------------------------------------------------
// Integrates the density from the input cube file:
// 1) read -> 2) integrate -> 3) print.
void Algorithm::integrate_density(Target &target)
{

    cube.read_density(target, out, "Cube");

    cube.int_density();

    out.print_density(target, cube);
}
//----------------------------------------------------------------------
// Solute-solvent coupling:
// 1) read solute density + solvent geometry
// 2) print characteristics
// 3) compute integrals
// 4) print results
void Algorithm::solute_dens_solvent_pot_fld(Target &target)
{
    //
    //  Read input files
    //
    cube_solute.read_density(target, out, "Solute");

    solv.read_solvent(target, out);
    //
    //  Print acceptor / donor density characteristics
    //
    out.print_solvent(target, solv);

    out.print_density(target, cube_solute, Parameters::solute_header);
    //
    //  Compute potential/field at solvent coordinates
    //
    integrals.solute_solvent_pot_fld(target, cube_solute, solv);
    //
    //  Print results
    //
    //    out.print_results_integrals(target, integrals);
}
//----------------------------------------------------------------------
