#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <array>
#include <string>
#include <string_view>

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

    constexpr double ToBohr = 1.8897261254578281; // Conversion from Angstrom to Bohr
    constexpr double ToAng = 1.0 / ToBohr;        // Conversion from Bohr to Angstrom

    constexpr double QMscrnFact = 0.2; // Screening factor for Coulomb integrals

    // Constraint for reduce density file
    constexpr int ncellmax = 10000000;

    // Accepted entries for the "what" keyword
    inline constexpr std::array<std::string_view, 4> accepted_what_entries = {
        "potential",
        "field",
        "potential+field",
        "fq"};

    // Accepted entries for the "parametrization" keyword
    inline constexpr std::array<std::string_view, 1> accepted_parametrization_entries = {
        "giovannini"};

    // Accepted solvent file extensions
    inline constexpr std::array<std::string_view, 1> accepted_solvent_file_extensions = {
        ".xyz"};

    // Header strings (declared here, defined in parameters.cpp)
    extern const std::string solute_header;

} // namespace Parameters

#endif // PARAMETERS_HPP
//----------------------------------------------------------------------
