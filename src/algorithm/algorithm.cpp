#include "algorithm.hpp"
#include "output.hpp"
#include "target.hpp"
#include "density.hpp"
#include "parameters.hpp"
#include "integrals.hpp"
#include "nanoparticle.hpp"

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
// Donor–acceptor electronic energy transfer:
// 1) read both densities
// 2) print characteristics
// 3) compute integrals
// 4) print results
void Algorithm::acceptor_donor(Target &target)
{
    //
    //  Read input files
    //
    cube_acceptor.read_density(target, out, "Acceptor");

    cube_donor.read_density(target, out, "Donor");
    //
    //   Print acceptor / donor density characteristics
    //
    out.print_density(target, cube_acceptor, Parameters::acceptor_header);
                        
    out.print_density(target, cube_donor, Parameters::donor_header);

    //
    //  Compute integrals
    //
    integrals.acceptor_donor(target, cube_acceptor, cube_donor);
    //
    //  Print results
    //
    out.print_results_integrals(target, integrals);
}
//----------------------------------------------------------------------
// Acceptor–nanoparticle coupling:
// 1) read acceptor density + nanoparticle
// 2) print characteristics
// 3) compute integrals
// 4) print results
void Algorithm::acceptor_np(Target &target)
{
    //
    //  Read input files
    //
    cube_acceptor.read_density(target, out, "Acceptor");

    np.read_nanoparticle(target, out);
    //
    //  Print acceptor / donor density characteristics
    //
    out.print_nanoparticle(np);

    out.print_density(target, cube_acceptor, Parameters::acceptor_header);
    //
    //  Compute integrals
    //
    integrals.acceptor_np(target, cube_acceptor, np);
    //
    //  Print results
    //
    out.print_results_integrals(target, integrals);
}
//----------------------------------------------------------------------
// Nanoparticle-mediated donor–acceptor electronic energy transfer:
// 1) read acceptor + donor + nanoparticle
// 2) print characteristics
// 3) compute integrals (DA and A–NP)
// 4) print results
void Algorithm::acceptor_np_donor(Target &target)
{
    //
    //  Read input files
    //
    cube_acceptor.read_density(target, out, "Acceptor");

    cube_donor.read_density(target, out, "Donor");

    np.read_nanoparticle(target, out);
    //
    //  Print nanoparticle / acceptor / donor density characteristics
    //
    out.print_nanoparticle(np);

    out.print_density(target, cube_acceptor, Parameters::acceptor_header);

    out.print_density(target,cube_donor, Parameters::donor_header);
    //
    //  Compute integrals
    //
    integrals.acceptor_donor(target, cube_acceptor, cube_donor);

    integrals.acceptor_np(target, cube_acceptor, np);
    //
    //  Print results
    //
    out.print_results_integrals(target, integrals);
}
//----------------------------------------------------------------------