#include "fq.hpp"
#include "solvent.hpp"

#include <stdexcept>
//----------------------------------------------------------------------
// Compute FQ charges from the potential and solvent parameters.
void FQ::calc_charges(Target &, Solvent &solv)
{

    
    for (int i = 0; i < solv.natoms; ++i)
    {
        solv.tempTqq[i][i] = solv.typeEta[solv.typeIndex[i]];
    }
}
