#include "fq.hpp"
#include "output.hpp"
#include "target.hpp"
#include "solvent.hpp"

#include <iostream>
#include <stdexcept>
//----------------------------------------------------------------------
// Compute FQ charges from the potential and solvent parameters.
void FQ::calc_charges(Target &target, Solvent &solv)
{
    // Check solvent has FQ parameters.
    //if (!solv.has_fq_parameters)
    //{
    //    throw std::runtime_error("Solvent does not have FQ parameters for charge calculation.");
    //} 
}