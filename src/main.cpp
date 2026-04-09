#include <exception>
#include <iostream>

#include "target.hpp"
#include "input.hpp"
#include "output.hpp"
#include "timer.hpp"
#include "algorithm.hpp"

//---------------------------------------------------------------------------
//                 __________      _____       __
//                / ____/ __ \    / ___/____  / /   _____  _____
//               / /_  / / / /    \__ \/ __ \/ / | / / _ \/ ___/
//              / __/ / /_/ /    ___/ / /_/ / /| |/ /  __/ /
//             /_/    \___\_\   /____/\____/_/ |___/\___/_/
//
//---------------------------------------------------------------------------
//
//                       Program by Pablo Grobas Illobre
//
//                         For any problem write to:
//                         pgrobasillobre@gmail.com
//
//---------------------------------------------------------------------------
/// @file main.cpp
/// @brief Main entry point of the FQSolver program.
///
/// This function sets up and coordinates the key components:
/// - Parses input arguments
/// - Reads and checks input
/// - Dispatches the appropriate computational algorithm based on the target
/// - Manages performance timing
/// - Handles error reporting and output
int main(int argc, char *argv[])
{

    Output out;

    try
    {
        // Instantiate components
        Timer timer;
        Input inp;
        Target target;

        // Parse input arguments
        inp.get_arguments(argc, argv, out, target);
        timer.initialize();
        timer.start("total");

        // Open output file. Check input file existence and extension
        out.open();
        inp.check_input_file(out);

        // Print banner, read input and print input info.
        out.print_banner();

        inp.read(target);
        inp.print_input_info(out, target);

        // Initialize algorithm instance with output and target references.
        Algorithm algorithm(out, target);

        switch (target.mode)
        {
        case TargetMode::IntegrateCube:
            algorithm.integrate_density(target);
            break;
        
        case TargetMode::Solute_Solvent_Pot_Field:
            algorithm.solute_dens_solvent_pot(target);
            break;

        default:
            throw std::runtime_error("No valid calculation target specified in input.");
        }

        // Finalize timing and output
        timer.finish("total");
        timer.conclude(out);

        out.close();
    }
    catch (const std::exception &e)
    {
        out.stream() << " Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
