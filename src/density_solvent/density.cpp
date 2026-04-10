#include "density.hpp"
#include "parameters.hpp"
#include "target.hpp"
#include "output.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

//----------------------------------------------------------------------
// Map atomic number Z to element label.
// Returns "X" for unknown/unsupported Z.
std::string Density::map_atomic_number_to_label(int Z) const
{
    static const std::vector<std::string> periodic_table = {
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"};
    return (Z > 0 && Z < (int)periodic_table.size()) ? periodic_table[Z] : "X";
}
//----------------------------------------------------------------------
// Read cube file, load metadata and density grid, build reduced density
// representation, and compute geometric centers.
void Density::read_density(Target &target, const Output &out, const std::string &what_dens)
{

    // Check density file final purpose: cube integration, acceptor, or donor density.
    std::string filepath;

    if (what_dens == "Cube")
    {
        filepath = target.density_file_integration;
    }
    else if (what_dens == "Solute")
    {
        filepath = target.solute_density_file;
    }
    else
    {
        throw std::runtime_error("Unknown density file mode to read.");
    }

    // Check file existance.
    std::ifstream infile(filepath);
    if (!infile)
    {
        throw std::runtime_error("File: " + filepath + "not found.");
    }

    // Skip header lines
    std::getline(infile, str1);
    std::getline(infile, str2);

    // Read grid and origin info
    infile >> natoms >> xmin >> ymin >> zmin;
    infile >> nx >> dx[0] >> dx[1] >> dx[2];
    infile >> ny >> dy[0] >> dy[1] >> dy[2];
    infile >> nz >> dz[0] >> dz[1] >> dz[2];

    // Ensure voxel matrix is diagonal. Compute voxel volume.
    if (dx[1] != 0.0 || dx[2] != 0.0 ||
        dy[0] != 0.0 || dy[2] != 0.0 ||
        dz[0] != 0.0 || dz[1] != 0.0)
    {
        throw std::runtime_error("Cube file conflict: dx, dy, dz matrix is not diagonal.");
    }
    volume = dx[0] * dy[1] * dz[2];

    // Initialize atom vectors
    atomic_number.resize(natoms);
    atomic_label.resize(natoms);
    atomic_charge.resize(natoms);
    x.resize(natoms);
    y.resize(natoms);
    z.resize(natoms);

    // Read atom block
    nelectrons = 0;
    for (int i = 0; i < natoms; ++i)
    {
        infile >> atomic_number[i] >> atomic_charge[i] >> x[i] >> y[i] >> z[i];
        nelectrons += atomic_number[i];
        atomic_label[i] = map_atomic_number_to_label(atomic_number[i]);
    }

    // Allocate 3D density grid
    rho.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz)));

    // Read density values, weight by voxel volume, and find maximum density value
    maxdens = 0.0;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                infile >> rho[i][j][k];
                rho[i][j][k] *= volume;
                maxdens = std::max(maxdens, std::abs(rho[i][j][k]));
            }
        }
    }

    //
    // Save reduced density of the cube, and calculate associated coordinates
    //
    if (!target.integrate_density)
    {
        rho_reduced.resize(Parameters::ncellmax);
        xyz.resize(Parameters::ncellmax);

        double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;

        for (int i = 0; i < nx; ++i)
        {
            x_tmp = xmin + dx[0] * i;

            for (int j = 0; j < ny; ++j)
            {
                y_tmp = ymin + dy[1] * j;

                for (int k = 0; k < nz; ++k)
                {
                    z_tmp = zmin + dz[2] * k;

                    if (std::abs(rho[i][j][k]) > maxdens * target.cutoff)
                    {

                        if (n_points_reduced >= Parameters::ncellmax)
                        {
                            throw std::runtime_error("Too many points (" + std::to_string(n_points_reduced) + ") in " + what_dens + " density file. Increase cutoff or ncellmax.");
                        }
                        rho_reduced[n_points_reduced] = rho[i][j][k];
                        xyz[n_points_reduced] = {x_tmp, y_tmp, z_tmp};
                        ++n_points_reduced;
                    }
                }
            }
        }
    }

    // Resize vectors to actual number of points
    rho_reduced.resize(n_points_reduced);
    xyz.resize(n_points_reduced);

    infile.close();

    // Calculate geometric center of the density and the molecule
    for (const auto &coord : xyz)
    {
        geom_center[0] += coord[0];
        geom_center[1] += coord[1];
        geom_center[2] += coord[2];
    }
    geom_center[0] /= n_points_reduced;
    geom_center[1] /= n_points_reduced;
    geom_center[2] /= n_points_reduced;

    for (int i = 0; i < natoms; ++i)
    {
        geom_center_mol[0] += x[i];
        geom_center_mol[1] += y[i];
        geom_center_mol[2] += z[i];
    }
    geom_center_mol[0] /= natoms;
    geom_center_mol[1] /= natoms;
    geom_center_mol[2] /= natoms;

    // Check both centers are the same up to a threshold of 0.3)
    // debugpgi: uncomment this after debug
    // if (std::abs(geom_center[0] - geom_center_mol[0]) > 0.3 ||
    //    std::abs(geom_center[1] - geom_center_mol[1]) > 0.3 ||
    //    std::abs(geom_center[2] - geom_center_mol[2]) > 0.3)
    //{
    // throw std::runtime_error("Geometric center of the density and molecule do not match.");
    //}
}
//----------------------------------------------------------------------
// Integrate the density grid (sum over all voxels).
void Density::int_density()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                integral += rho[i][j][k];
            }
        }
    }
}
