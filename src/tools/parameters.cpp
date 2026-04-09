#include "parameters.hpp"

/// @brief Defines global string constants declared in Parameters namespace.
///
/// This file contains the definitions for header labels and formatting strings
/// used throughout the FQSolver output system and density processing.
namespace Parameters
{
    const std::string solute_header = std::string(25, ' ') + "Solute Density Information";
    const std::string acceptor_header = std::string(25, ' ') + "Acceptor Density Information";
    const std::string donor_header = std::string(25, ' ') + "Donor Density Information";

    const std::string fret_start = " # fret quantities ------------------------#";
    const std::string fret_end = " # end fret quantities ------------------------";

    const std::string charges_header = " #        q_re                     q_im                coords_x                 coords_y                 coords_z";
    const std::string charges_and_dipoles_header = " # q_re,   q_im,   mu_re_x,   mu_re_y,   mu_re_z,  mu_im_x,   mu_im_y,   mu_im_z, coords_x, coords_y, coords_z";
}
