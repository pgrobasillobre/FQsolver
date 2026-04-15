#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <array>
#include <string>
#include <string_view>
#include <unordered_map>

//----------------------------------------------------------------------
/// @brief Defines global physical constants and header strings used throughout the application.
///
/// Contains the Parameters namespace, which centralizes physical constants,
/// unit conversions, numerical thresholds, and string headers relevant to
/// density parsing and output formatting in the FQSolver codebase.
namespace Parameters
{
    struct FQAtomParams
    {
        double chi;
        double eta;
    };

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

    // Accepted entries for the "kernel" keyword
    inline constexpr std::array<std::string_view, 3> accepted_kernel_entries = {
        "gaus",
        "ohno",
        "coul"};

    // Accepted solvent file extensions
    inline constexpr std::array<std::string_view, 2> accepted_solvent_file_extensions = {
        ".xyz",
        ".pdb"};

    inline const std::unordered_map<std::string, std::unordered_map<std::string, FQAtomParams>> fq_params = {
        {"giovannini",
         {
             {"water oxygen", {0.2908429850, 0.5625181140}},
             {"water hydrogen", {0.1675711970, 0.6093265770}},
         }}};

    inline const std::unordered_map<std::string, std::unordered_map<std::string, std::string>> fq_label_to_type = {
        {"giovannini",
         {
             {"OW", "water oxygen"},
             {"O", "water oxygen"},
             {"HW1", "water hydrogen"},
             {"HW2", "water hydrogen"},
             {"H1", "water hydrogen"},
             {"H2", "water hydrogen"},
             {"H", "water hydrogen"},
         }}};

    // Header strings (declared here, defined in parameters.cpp)
    extern const std::string solute_header;

} // namespace Parameters

#endif // PARAMETERS_HPP
//----------------------------------------------------------------------
