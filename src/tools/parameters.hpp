#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>

//----------------------------------------------------------------------
/// @brief Defines global physical constants and header strings used throughout the application.
/// 
/// Contains the Parameters namespace, which centralizes physical constants,
/// unit conversions, numerical thresholds, and string headers relevant to
/// density parsing and output formatting in the FQSolver codebase.
namespace Parameters
{

    // Physical constants
    constexpr double pi = 3.14159265358979323846;
    constexpr double sqrt_pi = 1.77245385090551602730; // Square root of pi

    constexpr double to_radians = pi / 180.0; 
    constexpr double to_degrees = 1.0 / to_radians; 

    constexpr double angle_thres = 0.1 * to_radians; // Threshold for angle comparison in radians = 6 degrees

    constexpr double ToBohr = 1.8897261254578281;      // Conversion from Angstrom to Bohr
    constexpr double ToAng = 1.0 / ToBohr;             // Conversion from Bohr to Angstrom

    constexpr double QMscrnFact = 0.2; // Screening factor for Coulomb integrals

    // Constraint for reduce density file
    constexpr int ncellmax = 10000000;

    // Header strings (declared here, defined in parameters.cpp)
    extern const std::string acceptor_header;
    extern const std::string donor_header;

    extern const std::string fret_start;
    extern const std::string fret_end;

    extern const std::string charges_header;
    extern const std::string charges_and_dipoles_header;

} // namespace Parameters

#endif // PARAMETERS_HPP
//----------------------------------------------------------------------
