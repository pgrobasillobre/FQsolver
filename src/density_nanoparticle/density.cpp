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
std::string Density::map_atomic_number_to_label(int Z) const {
    static const std::vector<std::string> periodic_table = {
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
    };
    return (Z > 0 && Z < (int)periodic_table.size()) ? periodic_table[Z] : "X";
}
//----------------------------------------------------------------------
// Compute angle between two vectors using dot product formula.
// Returns angle in radians. Throws if a zero-length vector is given.
double Density::compute_angle_between_vectors(const std::array<double, 3>& vec1, const std::array<double, 3>& vec2) const {
    double dot_product = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    double norm1 = std::sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
    double norm2 = std::sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);

    if (norm1 == 0.0 || norm2 == 0.0) {
        throw std::runtime_error("Cannot compute angle with zero-length vector.");
    }

    // Handle decimal points to avoid floating-point precision issues
    double division = dot_product / (norm1 * norm2);

    if (division < -1.0) {
        division = -1.0; // Clamp to avoid NaN
    } else if (division > 1.0) {
        division = 1.0; // Clamp to avoid NaN
    }

    return std::acos(division); // Returns angle in radians;
}
//----------------------------------------------------------------------
// Rotate a 3D vector by an angle around a given axis ("x", "y", "z").
std::array<double, 3> Density::rotate_vector(const std::array<double, 3>& vec, const double angle, const std::string& axis) const {

    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    std::array<double, 3> vec_rot = {0.0, 0.0, 0.0}; 

    if (axis == "x") {
        vec_rot[0] = vec[0];
        vec_rot[1] = vec[1] * cos_angle - vec[2] * sin_angle;
        vec_rot[2] = vec[1] * sin_angle + vec[2] * cos_angle;
    }
    else if (axis == "y") {
        vec_rot[0] = vec[0] * cos_angle + vec[2] * sin_angle;
        vec_rot[1] = vec[1];
        vec_rot[2] = -vec[0] * sin_angle + vec[2] * cos_angle;
    }
    else if (axis == "z") {
        vec_rot[0] = vec[0] * cos_angle - vec[1] * sin_angle;
        vec_rot[1] = vec[0] * sin_angle + vec[1] * cos_angle;
        vec_rot[2] = vec[2];
    }
    else {
        throw std::runtime_error("Unknown rotation axys: " + axis);
    }
    
    return vec_rot;
}
//----------------------------------------------------------------------
// Rotate all reduced density coordinates by an angle around an axis.
std::vector<std::array<double, 3>> Density::rotate_density(const double angle, const std::string& axis) const {

    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    std::vector<std::array<double, 3>> xyz_rot(n_points_reduced, {0.0, 0.0, 0.0});
    
    if (axis == "x") {
        for (int i = 0; i < n_points_reduced; ++i)
        {
            xyz_rot[i][0] = xyz[i][0];
            xyz_rot[i][1] = xyz[i][1] * cos_angle - xyz[i][2] * sin_angle;
            xyz_rot[i][2] = xyz[i][1] * sin_angle + xyz[i][2] * cos_angle;
        }
    }
    else if (axis == "y") {
        for (int i = 0; i < n_points_reduced; ++i)
        {
            xyz_rot[i][0] = xyz[i][0] * cos_angle + xyz[i][2] * sin_angle;
            xyz_rot[i][1] = xyz[i][1];
            xyz_rot[i][2] = -xyz[i][0] * sin_angle + xyz[i][2] * cos_angle;
        }
    }
    else if (axis == "z") {
        for (int i = 0; i < n_points_reduced; ++i)
        {
            xyz_rot[i][0] = xyz[i][0] * cos_angle - xyz[i][1] * sin_angle;
            xyz_rot[i][1] = xyz[i][0] * sin_angle + xyz[i][1] * cos_angle;
            xyz_rot[i][2] = xyz[i][2];
        }
    }
    else {
        throw std::runtime_error("Unknown rotation axys: " + axis);
    }
    
    return xyz_rot;
}
//----------------------------------------------------------------------
// Read cube file, load metadata and density grid, build reduced density
// representation, and compute geometric centers.
void Density::read_density(Target& target, const Output &out, const std::string& what_dens) {

    // Check density file final purpose: cube integration, acceptor, or donor density.
    std::string filepath;

    if (what_dens=="Cube"){
       filepath = target.density_file_integration;
    } else if (what_dens=="Acceptor"){
       filepath = target.acceptor_density_file;
    } else if (what_dens=="Donor"){
       filepath = target.donor_density_file;
    } else {
        throw std::runtime_error("Unknown density file mode to read.");
    }

    // Check file existance.
    std::ifstream infile(filepath);
    if (!infile) {
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
        dz[0] != 0.0 || dz[1] != 0.0) {
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
    for (int i = 0; i < natoms; ++i) {
        infile >> atomic_number[i] >> atomic_charge[i] >> x[i] >> y[i] >> z[i];
        nelectrons += atomic_number[i];
        atomic_label[i] = map_atomic_number_to_label(atomic_number[i]);
    }

    // Allocate 3D density grid
    rho.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz)));

    // Read density values, weight by voxel volume, and find maximum density value
    maxdens = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                infile >> rho[i][j][k];
                rho[i][j][k] *= volume;
                maxdens = std::max(maxdens, std::abs(rho[i][j][k]));
            }
        }
    }

    //
    // Save reduced density of the cube, and calculate associated coordinates
    //
    if (!target.integrate_density) {
        rho_reduced.resize(Parameters::ncellmax);
        xyz.resize(Parameters::ncellmax);
    
        double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
    
        for (int i = 0; i < nx; ++i) {
            x_tmp = xmin + dx[0] * i;
    
            for (int j = 0; j < ny; ++j) {
                y_tmp = ymin + dy[1] * j;
    
                for (int k = 0; k < nz; ++k) {
                    z_tmp = zmin + dz[2] * k;
    
                    if (std::abs(rho[i][j][k]) > maxdens * target.cutoff || target.calc_overlap_int) {

                        if (n_points_reduced >= Parameters::ncellmax) {
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
    if (std::abs(geom_center[0] - geom_center_mol[0]) > 0.3 ||
        std::abs(geom_center[1] - geom_center_mol[1]) > 0.3 ||
        std::abs(geom_center[2] - geom_center_mol[2]) > 0.3) {
        throw std::runtime_error("Geometric center of the density and molecule do not match.");
    }

    // Rotate acceptor and/or donor density based on the transition dipole alignment vector, if requested
    //   1. Compute angle between the transition dipole and the reference vector
    //   2. Rotate cube coordinates
    if (target.rotate_acceptor || target.rotate_donor)
    {
        rotate_transition_dipole(target, out, what_dens);    
        rotate_cube_coordinates(target, out, what_dens);
    } 
}
//----------------------------------------------------------------------
// Integrate the density grid (sum over all voxels).
void Density::int_density() {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                integral += rho[i][j][k];
            }
        }
    }
}
//---------------------------------------------------------------------- 
// Rotate the transition dipole moment of donor/acceptor to align with
// a reference vector. Handles both acceptor and donor cases, performs
// consistency checks, and prints debug output if enabled.
void Density::rotate_transition_dipole(Target &target, const Output & out, const std::string& what_dens) {
    //
    // Move reference vector to the molecular origin of coordinates
    //
    double ToAng = Parameters::ToAng; // Conversion factor from Bohr to Angstrom
    std::array<double, 3> ref_vector_trans = {0.0, 0.0, 0.0}; ///< Geometric center of the molecule
    
    if (what_dens == "Acceptor") {
        ref_vector_trans[0] = target.acceptor_ref_vector[0] - geom_center_mol[0] * ToAng;
        ref_vector_trans[1] = target.acceptor_ref_vector[1] - geom_center_mol[1] * ToAng;
        ref_vector_trans[2] = target.acceptor_ref_vector[2] - geom_center_mol[2] * ToAng;
    }
    else if (what_dens == "Donor") {
        ref_vector_trans[0] = target.donor_ref_vector[0] - geom_center_mol[0] * ToAng;
        ref_vector_trans[1] = target.donor_ref_vector[1] - geom_center_mol[1] * ToAng;
        ref_vector_trans[2] = target.donor_ref_vector[2] - geom_center_mol[2] * ToAng;
    }
    else {
        throw std::runtime_error("Unknown density type for transition dipole rotation: " + what_dens);
    }

    //
    // Calculate the angle between the transition dipole and the reference vector
    //
    double theta = 0.0;
    if (what_dens == "Acceptor") theta = compute_angle_between_vectors(target.acceptor_transdip, ref_vector_trans);
    if (what_dens == "Donor")    theta = compute_angle_between_vectors(target.donor_transdip, ref_vector_trans);

    //
    // Rotate the transition dipole moment
    //
    std::array<double, 3> transdip_rot = {0.0, 0.0, 0.0}; 

    if (what_dens == "Acceptor") transdip_rot = rotate_vector(target.acceptor_transdip, theta, target.rotation_axys);
    if (what_dens == "Donor")    transdip_rot = rotate_vector(target.donor_transdip, theta, target.rotation_axys);

    // Check angle between rotated transdip and reference vector
    // If they are not coincident, change rotaton direction and apply the rotation again
    double theta_check = compute_angle_between_vectors(transdip_rot, ref_vector_trans);

    if (std::abs(theta_check) > Parameters::angle_thres) {
        theta = -theta; // Change rotation direction

        transdip_rot = {0.0, 0.0, 0.0}; 

        if (what_dens == "Acceptor") transdip_rot = rotate_vector(target.acceptor_transdip, theta, target.rotation_axys);
        if (what_dens == "Donor")    transdip_rot = rotate_vector(target.donor_transdip, theta, target.rotation_axys);

        // Check angle again and throw error in case still not aligned
        theta_check = compute_angle_between_vectors(transdip_rot, ref_vector_trans);
        if (std::abs(theta_check) > Parameters::angle_thres) {
           out.stream() << "\n " << "Calculated angle: " << std::fixed << std::setprecision(3) << std::setw(8)
            << (theta_check * Parameters::to_degrees) << "\n \n";

            out.stream() << " Geom-centered reference vector: "
            << std::setw(8) << ref_vector_trans[0] << "  "
            << std::setw(8) << ref_vector_trans[1] << "  "
            << std::setw(8) << ref_vector_trans[2] << '\n' 
            << std::flush;

            throw std::runtime_error("Alignment with reference vector was not possible");
        }
    }

    if (what_dens == "Acceptor") {
        target.acceptor_density_rotation_angle = theta;
        target.acceptor_density_rotation_angle_check = theta_check;

        target.acceptor_transdip_rot[0] = transdip_rot[0];
        target.acceptor_transdip_rot[1] = transdip_rot[1];
        target.acceptor_transdip_rot[2] = transdip_rot[2];
    }
    else if (what_dens == "Donor") {
        target.donor_density_rotation_angle = theta;
        target.donor_density_rotation_angle_check = theta_check;

        target.donor_transdip_rot[0] = transdip_rot[0];
        target.donor_transdip_rot[1] = transdip_rot[1];
        target.donor_transdip_rot[2] = transdip_rot[2];
    }

    //
    // Print rotated transition dipole
    //
    if (target.debug >= 1) {
        std::string nmd_file = "debug/transition_dipole_" + what_dens;
        std::string nmd_file_rotated = "debug/transition_dipole_" + what_dens + "_ROTATED";

        if (what_dens == "Acceptor") 
        {
            out.print_transdip_nmd(nmd_file,         target.acceptor_transdip,     geom_center_mol);
            out.print_transdip_nmd(nmd_file_rotated, target.acceptor_transdip_rot, geom_center_mol);
        }
        else if (what_dens == "Donor") 
        {
            out.print_transdip_nmd(nmd_file,         target.donor_transdip,     geom_center_mol);
            out.print_transdip_nmd(nmd_file_rotated, target.donor_transdip_rot, geom_center_mol);
        }
    }
}
//---------------------------------------------------------------------- 
// Rotate cube coordinates by saved rotation angle, translate to origin,
// rotate, then translate back. Used when aligning densities.
void Density::rotate_cube_coordinates(const Target &target, const Output & out, const std::string& what_dens)  {

    if (target.debug >= 1) out.print_cube_coordinates(what_dens, n_points_reduced, xyz);

    // Save angle
    double angle = 0.0;

    if (what_dens == "Acceptor") {
        angle = target.acceptor_density_rotation_angle;
    }
    else if (what_dens == "Donor") {
        angle = target.donor_density_rotation_angle;
    }
    else {
        throw std::runtime_error("Unknown density type for transition dipole rotation: " + what_dens);
    }

    //
    // Translate density center to the origin of coordinates
    //
    for (int i = 0; i < n_points_reduced; ++i)
    {
        xyz[i][0] = xyz[i][0] - geom_center[0];
        xyz[i][1] = xyz[i][1] - geom_center[1];
        xyz[i][2] = xyz[i][2] - geom_center[2];       
    }

    //
    // Rotate translated density
    //
    std::vector<std::array<double, 3>> xyz_rot(n_points_reduced, {0.0, 0.0, 0.0});

    xyz_rot = rotate_density(angle, target.rotation_axys);

    //
    // Translate rotated density to initial position and save
    //
    for (int i = 0; i < n_points_reduced; ++i)
    {
        xyz[i][0] = xyz_rot[i][0] + geom_center[0];
        xyz[i][1] = xyz_rot[i][1] + geom_center[1];
        xyz[i][2] = xyz_rot[i][2] + geom_center[2];       
    }

    if (target.debug >= 1) out.print_cube_coordinates(what_dens + "_ROTATED_", n_points_reduced, xyz);

}
//----------------------------------------------------------------------
