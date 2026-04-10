#include "output.hpp"
#include "parameters.hpp"
#include "integrals.hpp"
#include "solvent.hpp"

#include <iostream>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <string>
#include <ostream>
#include <numeric>
#include <cmath>

// Constructor: Initializes internal fields and file stream.
Output::Output() {}

//----------------------------------------------------------------------
// Creates the output filename by replacing `.inp` with `.log`.
// For example, "test.inp" → "test.log".
void Output::out_file_fill(const std::string &in_file)
{
    output_filename = in_file.substr(0, in_file.size() - 4) + ".log";
}
//----------------------------------------------------------------------
// Opens the output file stream for writing.
void Output::open()
{
    log_stream.open(output_filename, std::ios::out);
    if (!log_stream.is_open())
    {
        throw std::runtime_error("Failed to open output file: " + output_filename);
    }
}
//----------------------------------------------------------------------
// Returns the internal output stream object used for writing.
std::ofstream &Output::stream() const
{
    return const_cast<std::ofstream &>(log_stream);
}
//----------------------------------------------------------------------
// Closes the output stream if it is open.
void Output::close()
{
    if (log_stream.is_open())
    {
        log_stream.close();
    }
}
//----------------------------------------------------------------------
// Prints a banner for FQSolver at the beginning of the output.
void Output::print_banner()
{
    const std::string indent = std::string(18, ' ');

    log_stream << " " << sticks << "\n \n";

    log_stream << indent << "    __________      _____       __                \n";
    log_stream << indent << "   / ____/ __ \\    / ___/____  / /   _____  _____ \n";
    log_stream << indent << "  / /_  / / / /    \\__ \\/ __ \\/ / | / / _ \\/ ___/ \n";
    log_stream << indent << " / __/ / /_/ /    ___/ / /_/ / /| |/ /  __/ /     \n";
    log_stream << indent << "/_/    \\___\\_\\   /____/\\____/_/ |___/\\___/_/      \n";
    log_stream << indent << "                                                  \n";

    log_stream << " \n " << sticks << "\n \n";
    log_stream << std::string(25, ' ') << "Program by Pablo Grobas Illobre\n \n";
    log_stream << " " << sticks << "\n \n";
    log_stream.flush();
}
//----------------------------------------------------------------------
// Logs density file information: grid size, atoms, and dipole rotation (if requested).
// Also prints integral of the density, if present.
void Output::print_density(const Target &target, const Density &cube, std::optional<std::string> header)
{

    std::string solute_header = Parameters::solute_header;
    std::string filepath;

    if (header == solute_header)
    {
        filepath = target.solute_density_file;
    }
    else
    {
        filepath = target.density_file_integration;
    }

    // Extract the filename from the full path
    std::string filename = std::filesystem::path(filepath).filename().string();

    if (header.has_value())
    {
        log_stream << header.value() << "\n \n";
    }
    else
    {
        log_stream << std::string(29, ' ') << "Density Information                    \n \n";
    }
    log_stream << " " << sticks << "\n \n";
    log_stream << std::string(3, ' ') << "Density File: " << filename << "\n \n";
    log_stream << std::string(3, ' ') << "Density Grid (CUBE format): \n \n";

    print_formatted_line1(log_stream, cube.natoms, cube.xmin, cube.ymin, cube.zmin);
    print_formatted_line1(log_stream, cube.nx, cube.dx[0], cube.dx[1], cube.dx[2]);
    print_formatted_line1(log_stream, cube.ny, cube.dy[0], cube.dy[1], cube.dy[2]);
    print_formatted_line1(log_stream, cube.nz, cube.dz[0], cube.dz[1], cube.dz[2]);
    log_stream << " \n";
    log_stream << "     Total number of grid points: " << cube.nx * cube.ny * cube.nz << "\n";
    if (!header.has_value())
    {
        log_stream << " \n"; // Integrate cube case
    }
    else if (header.has_value())
    {
        log_stream << "     ---> Reduced density points: " << cube.n_points_reduced << "\n \n"; // Integrate cube case
    }

    log_stream << std::string(3, ' ') << "Associated molecular coordinates (Å): \n \n";
    for (int i = 0; i < cube.natoms; ++i)
    {
        print_formatted_line2(log_stream, std::string(cube.atomic_label[i]),
                              cube.x[i] * Parameters::ToAng,
                              cube.y[i] * Parameters::ToAng,
                              cube.z[i] * Parameters::ToAng);
    }

    if (cube.integral > 0.0)
    {
        log_stream << " \n";
        log_stream << "    ============================================================\n";
        log_stream << "     Integrated electron density -->    " << std::fixed << std::setprecision(14) << cube.integral << "\n";
        log_stream << "    ============================================================\n";
    }

    log_stream << " \n " << sticks << "\n\n";
}
//----------------------------------------------------------------------
/// Prints solvent geometry information: number of solvent atoms and their coordinates.
void Output::print_solvent(const Target &target, const Solvent &solv)
{
    log_stream << std::string(29, ' ') << "Solvent Information                    \n \n";
    log_stream << " " << sticks << "\n \n";
    log_stream << std::string(3, ' ') << "Solvent Geometry File: " << std::filesystem::path(target.solvent_file).filename().string() << "\n \n";
    log_stream << std::string(3, ' ') << "Number of solvent atoms: " << solv.natoms << "\n \n";
    log_stream << std::string(3, ' ') << "Solvent Atomic Coordinates (Å): \n \n";
    for (int i = 0; i < solv.natoms; ++i)
    {
        print_formatted_line2(log_stream, std::string(solv.atomic_label[i]),
                              solv.xyz[i][0]*Parameters::ToAng,
                              solv.xyz[i][1]*Parameters::ToAng,
                              solv.xyz[i][2]*Parameters::ToAng);
    }
    log_stream << " \n " << sticks << "\n\n";
}
//----------------------------------------------------------------------
// Prints a single formatted line for grid or voxel information (used in density headers).
void Output::print_formatted_line1(std::ostream &out, int i, double a, double b, double c)
{
    char line[100];
    std::snprintf(line, sizeof(line), "   %5d %15.7E %15.7E %15.7E\n", i, a, b, c);
    out << line;
}
// ----------------------------------------------------------------------
// Prints a formatted atomic line for XYZ positions (used in CUBE or NMD formats).
void Output::print_formatted_line2(std::ostream &out, const std::string atom, double x, double y, double z) const
{
    char line[100];
    std::snprintf(line, sizeof(line), "       %-2s  %12.6f  %12.6f  %12.6f\n", atom.c_str(), x, y, z);
    out << line;
}
// ----------------------------------------------------------------------
// Prints a formatted atomic line with larger spacing (used for nanoparticle logging).
void Output::print_formatted_line3(std::ostream &out, const std::string atom, double x, double y, double z)
{
    char line[100];
    std::snprintf(line, sizeof(line), "              %-2s  %19.6f %19.6f %19.6f\n", atom.c_str(), x, y, z);
    out << line;
}
//----------------------------------------------------------------------
// Logs the computed interaction integrals depending on the selected target mode.
// Includes: Coulomb, overlap, modulus, and Keet rates.
void Output::print_results_integrals(const Target &target, const Integrals &integrals)
{
    std::array<double, 2> v_tot = {0.0, 0.0};
    double v_mod = 0.0;

    // Print header
    log_stream << std::string(36, ' ') << "RESULTS\n\n";
    log_stream << " " << sticks << " \n\n";

    switch (target.mode)
    {

        // case TargetMode::Acceptor_Donor:

        //    log_stream << std::string(5, ' ') << "Acceptor-Donor Coulomb  : " << std::fixed << std::setw(25) << std::setprecision(16) << integrals.coulomb_acceptor_donor << "  a.u.\n";

        //    v_tot[0] = integrals.coulomb_acceptor_donor + integrals.overlap_acceptor_donor;

        //    v_mod = std::sqrt(std::inner_product(v_tot.begin(), v_tot.end(), v_tot.begin(), 0.0));

        //    log_stream << std::string(37, ' ') << std::string(26, '-') << "\n";
        //    log_stream
        //        << std::string(5, ' ') << "Total Potential         : " << std::fixed << std::setw(25) << std::setprecision(16) << v_tot[0] << "  a.u.\n\n";
        //    log_stream << std::string(5, ' ') << "Total Potential Modulus : " << std::fixed << std::setw(25) << std::setprecision(16) << v_mod << "  a.u.\n\n";

        //    log_stream << " " << sticks << "\n\n";
        //    log_stream.flush();

        //    break;

    case TargetMode::Solute_Solvent_Pot_Fld:

        // log_stream << std::string(5, ' ') << "Acceptor-NP Interaction : " << std::fixed << std::setw(25) << std::setprecision(16) << integrals.overlap_acceptor_nanoparticle[0] << " + " << integrals.overlap_acceptor_nanoparticle[1] << " i  a.u.\n\n";
        // log_stream << " " << sticks << "\n\n";
        // log_stream.flush();

        // break;

        // case TargetMode::Acceptor_NP_Donor:

        //    log_stream << std::string(5, ' ') << "Acceptor-Donor Coulomb  : " << std::fixed << std::setw(25) << std::setprecision(16) << integrals.coulomb_acceptor_donor << "  a.u.\n";

        //    log_stream << std::string(5, ' ') << "Acceptor-NP Interaction : " << std::fixed << std::setw(25) << std::setprecision(16) << integrals.overlap_acceptor_nanoparticle[0] << " + " << integrals.overlap_acceptor_nanoparticle[1] << " i  a.u.\n";
        //    log_stream.flush();

        //    v_tot[0] = integrals.coulomb_acceptor_donor + integrals.overlap_acceptor_donor + integrals.overlap_acceptor_nanoparticle[0];
        //    v_tot[1] = integrals.overlap_acceptor_nanoparticle[1];

        //    v_mod = std::sqrt(std::inner_product(v_tot.begin(), v_tot.end(), v_tot.begin(), 0.0));

        //    log_stream << std::string(37, ' ') << std::string(26, '-') << "\n";
        //    log_stream << std::string(5, ' ') << "Total Potential         : " << std::fixed << std::setw(25) << std::setprecision(16) << v_tot[0] << " + " << v_tot[1] << " i  a.u.\n\n";
        //    log_stream << std::string(5, ' ') << "Total Potential Modulus : " << std::fixed << std::setw(25) << std::setprecision(16) << v_mod << "  a.u.\n\n";

        //    log_stream << " " << sticks << "\n\n";
        //    log_stream.flush();

        //    break;

    case TargetMode::None:
    default:
        throw std::runtime_error("No valid calculation target specified in input.");
    }
}
//----------------------------------------------------------------------
// Writes cube point XYZ coordinates to a debug .xyz file.
void Output::print_cube_coordinates(const std::string what_dens,
                                    const int n_points,
                                    const std::vector<std::array<double, 3>> &xyz) const
{
    double ToAng = Parameters::ToAng;
    std::string infile = "debug/" + what_dens + "_cube_points.xyz";

    std::ofstream cubefile(infile, std::ios::out);
    if (!cubefile)
    {
        throw std::runtime_error("Cannot open file: " + infile + ".xyz");
    }

    cubefile << n_points;
    cubefile << "\ncube coordinates\n";

    for (int i = 0; i < n_points; ++i)
    {
        print_formatted_line2(cubefile, "H",
                              xyz[i][0] * ToAng,
                              xyz[i][1] * ToAng,
                              xyz[i][2] * ToAng);
    }
}
//----------------------------------------------------------------------
// Writes nanoparticle atom positions to .xyz file and dipole data to .nmd files.
// Includes real and imaginary components if available.
// void Output::print_np_coords_dipoles(const std::string infile, const Nanoparticle &np) const
//{
//    double ToAng = Parameters::ToAng;
//
//    std::ofstream npfile(infile + ".xyz", std::ios::out);
//    if (!npfile)
//    {
//        throw std::runtime_error("Cannot open file: " + infile + ".xyz");
//    }
//
//    npfile << np.natoms << "\n";
//    npfile << "NP coordinates\n";
//
//    for (int i = 0; i < np.natoms; ++i)
//    {
//        print_formatted_line2(npfile, "Xx",
//                              np.xyz[i][0] * ToAng,
//                              np.xyz[i][1] * ToAng,
//                              np.xyz[i][2] * ToAng);
//    }
//
//    npfile.close();
//
//    //
//    // Print dipoles, if present
//    //
//    if (np.charges_and_dipoles)
//    {
//
//        //
//        // Real part
//        //
//        std::ofstream dip_re(infile + "_re.nmd", std::ios::out);
//        if (!dip_re)
//        {
//            throw std::runtime_error("Cannot open file: " + infile + ".nmd");
//        }
//
//        dip_re << "coordinates";
//        for (int i = 0; i < np.natoms; ++i)
//        {
//            char line[100];
//            std::snprintf(line, sizeof(line), "%10.5f  %10.5f  %10.5f  ", np.xyz[i][0] * ToAng, np.xyz[i][1] * ToAng, np.xyz[i][2] * ToAng);
//            dip_re << line;
//        }
//        dip_re << "\nmode 1";
//        for (int i = 0; i < np.natoms; ++i)
//        {
//            char line[100];
//            std::snprintf(line, sizeof(line), "%10.5f  %10.5f  %10.5f  ", np.mu[i][0], np.mu[i][1], np.mu[i][2]);
//            dip_re << line;
//        }
//        dip_re << "\n";
//        dip_re.close();
//
//        //
//        // Imaginary part
//        //
//        std::ofstream dip_im(infile + "_im.nmd", std::ios::out);
//        if (!dip_im)
//        {
//            throw std::runtime_error("Cannot open file: " + infile + ".nmd");
//        }
//
//        dip_im << "coordinates";
//        for (int i = 0; i < np.natoms; ++i)
//        {
//            char line[100];
//            std::snprintf(line, sizeof(line), "%10.5f  %10.5f  %10.5f  ", np.xyz[i][0] * ToAng, np.xyz[i][1] * ToAng, np.xyz[i][2] * ToAng);
//            dip_im << line;
//        }
//        dip_im << "\nmode 1";
//        for (int i = 0; i < np.natoms; ++i)
//        {
//            char line[100];
//            std::snprintf(line, sizeof(line), "%10.5f  %10.5f  %10.5f  ", np.mu[i][3], np.mu[i][4], np.mu[i][5]);
//            dip_im << line;
//        }
//        dip_im << "\n";
//        dip_im.close();
//    }
//}
//----------------------------------------------------------------------
