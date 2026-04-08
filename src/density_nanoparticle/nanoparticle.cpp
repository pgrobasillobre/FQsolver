#include "nanoparticle.hpp"
#include "output.hpp"
#include "target.hpp"
#include "parameters.hpp"
#include "string_manipulation.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

//----------------------------------------------------------------------
// Load nanoparticle data: coordinates, charges, and optionally dipoles.
// Performs automatic format detection and rotation if requested.
void Nanoparticle::read_nanoparticle(const Target &target, const Output &out)
{
  // Check file existance.
  std::string filepath = target.nanoparticle_file;

  std::ifstream infile(filepath);
  if (!infile)
  {
    throw std::runtime_error("File: " + filepath + "not found.");
  }

  // Check if FRET quantities are present
  bool is_fret_quantities_present = false;
  is_fret_quantities_present = go_to_string(infile, Parameters::fret_start);

  // Check if FRET end marker is present
  bool is_fret_end_marker_present = false;
  is_fret_end_marker_present = go_to_string(infile, Parameters::fret_end);

  // Rewind and come back to fret_start
  infile.clear();
  infile.seekg(0, std::ios::beg);
  is_fret_quantities_present = go_to_string(infile, Parameters::fret_start);

  // Check charges / charges + dipoles
  std::string line;
  std::getline(infile, line);
  if (line == Parameters::charges_header)
  {
    charges = true;
    nanoparticle_model = "charges";
  }
  else if (line == Parameters::charges_and_dipoles_header)
  {
    charges_and_dipoles = true;
    nanoparticle_model = "charges + dipoles";
  }
  else
  {
    throw std::runtime_error("Expected header line with charges or charges and dipoles, got: " + line);
  }

  // Read charges / charges and dipoles.
  while (std::getline(infile, line))
  {
    if (line == Parameters::fret_end)
      break;

    std::istringstream iss(line);
    std::array<double, 3> coords{};

    if (charges)
    {
      std::array<double, 2> charge{};
      iss >> charge[0] >> charge[1] >> coords[0] >> coords[1] >> coords[2];
      q.push_back(charge);
    }
    else if (charges_and_dipoles)
    {
      std::array<double, 2> charge{};
      std::array<double, 6> dipole{};
      iss >> charge[0] >> charge[1] >> dipole[0] >> dipole[1] >> dipole[2] >> dipole[3] >> dipole[4] >> dipole[5] >> coords[0] >> coords[1] >> coords[2];
      q.push_back(charge);
      mu.push_back(dipole);
    }

    xyz.push_back(coords);
  }
  natoms = xyz.size();

  // Compute geometrical center
  for (const auto &coord : xyz)
  {
    geom_center[0] += coord[0];
    geom_center[1] += coord[1];
    geom_center[2] += coord[2];
  }
  geom_center[0] /= natoms;
  geom_center[1] /= natoms;
  geom_center[2] /= natoms;

  // Rotate nanoparticle, if required.
  if (target.rotate_nanoparticle) rotate_np_coords_and_dipoles(target, out);
}
//----------------------------------------------------------------------
// Rotate nanoparticle coordinates and dipoles (if present)
// using donor's rotation parameters.
void Nanoparticle::rotate_np_coords_and_dipoles(const Target &target, const Output &out)
{
  //
  // Print rotated nanoparticle coordinates and dipoles for debug
  //
  if (target.debug >= 1) out.print_np_coords_dipoles("debug/np",*this);

  // Rotate with the same angle of the donor, which has induced the charges and dipoles
  double angle = target.donor_density_rotation_angle;
  
  //
  // Translate density center to the origin of coordinates
  //
  for (int i = 0; i < natoms; ++i)
  {
    xyz[i][0] = xyz[i][0] - geom_center[0];
    xyz[i][1] = xyz[i][1] - geom_center[1];
    xyz[i][2] = xyz[i][2] - geom_center[2];       
  }
  
  //
  // Rotate translated density
  //
  std::vector<std::array<double, 3>> xyz_rot(natoms, {0.0, 0.0, 0.0});

  xyz_rot = rotate_np_coords(angle, target.rotation_axys);

  //
  // Translate rotated nanoparticle to initial position and save
  //
  for (int i = 0; i < natoms; ++i)
  {
    xyz[i][0] = xyz_rot[i][0] + geom_center[0];
    xyz[i][1] = xyz_rot[i][1] + geom_center[1];
    xyz[i][2] = xyz_rot[i][2] + geom_center[2];       
  }

  //
  // Rotate nanoparticle dipoles, if present
  //
  if (charges_and_dipoles) rotate_np_dipoles(angle, target.rotation_axys);

  //
  // Print rotated nanoparticle coordinates and dipoles for debug
  //
  if (target.debug >= 1) out.print_np_coords_dipoles("debug/np_rot",*this);
}
//----------------------------------------------------------------------
// Rotate nanoparticle coordinates around specified axis.
std::vector<std::array<double, 3>> Nanoparticle::rotate_np_coords(const double angle, const std::string& axis) const {

    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    std::vector<std::array<double, 3>> xyz_rot(natoms, {0.0, 0.0, 0.0});
    
    if (axis == "x") {
        for (int i = 0; i < natoms; ++i)
        {
            xyz_rot[i][0] = xyz[i][0];
            xyz_rot[i][1] = xyz[i][1] * cos_angle - xyz[i][2] * sin_angle;
            xyz_rot[i][2] = xyz[i][1] * sin_angle + xyz[i][2] * cos_angle;
        }
    }
    else if (axis == "y") {
        for (int i = 0; i < natoms; ++i)
        {
            xyz_rot[i][0] = xyz[i][0] * cos_angle + xyz[i][2] * sin_angle;
            xyz_rot[i][1] = xyz[i][1];
            xyz_rot[i][2] = -xyz[i][0] * sin_angle + xyz[i][2] * cos_angle;
        }
    }
    else if (axis == "z") {
        for (int i = 0; i < natoms; ++i)
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
// Rotate real and imaginary components of each dipole
// using the same angle/axis as the coordinates.
void Nanoparticle::rotate_np_dipoles(const double angle, const std::string& axis) 
{
    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    double x_tmp = 0.0;
    double y_tmp = 0.0;
    double z_tmp = 0.0;
        
    if (axis == "x") {
        for (int i = 0; i < natoms; ++i)
        {
            // Real part
            x_tmp = mu[i][0];
            y_tmp = mu[i][1] * cos_angle - mu[i][2] * sin_angle;
            z_tmp = mu[i][1] * sin_angle + mu[i][2] * cos_angle;

            mu[i][0] = x_tmp;
            mu[i][1] = y_tmp;
            mu[i][2] = z_tmp;

            // Imaginary part
            x_tmp = mu[i][3];
            y_tmp = mu[i][4] * cos_angle - mu[i][5] * sin_angle;
            z_tmp = mu[i][4] * sin_angle + mu[i][5] * cos_angle;

            mu[i][3] = x_tmp;
            mu[i][4] = y_tmp;
            mu[i][5] = z_tmp;
        }
    }
    else if (axis == "y") {
        for (int i = 0; i < natoms; ++i)
        {
            // Real part
            x_tmp = mu[i][0] * cos_angle + mu[i][2] * sin_angle;
            y_tmp = mu[i][1];
            z_tmp = -mu[i][0] * sin_angle + mu[i][2] * cos_angle;

            mu[i][0] = x_tmp;
            mu[i][1] = y_tmp;
            mu[i][2] = z_tmp;

            // Imaginary part
            x_tmp = mu[i][3] * cos_angle + mu[i][5] * sin_angle;
            y_tmp = mu[i][4];
            z_tmp = -mu[i][3] * sin_angle + mu[i][5] * cos_angle;

            mu[i][3] = x_tmp;
            mu[i][4] = y_tmp;
            mu[i][5] = z_tmp;
        }
    }
    else if (axis == "z") {
        for (int i = 0; i < natoms; ++i)
        {
            // Real part
            x_tmp = mu[i][0] * cos_angle - mu[i][1] * sin_angle;
            y_tmp = mu[i][0] * sin_angle + mu[i][1] * cos_angle;
            z_tmp = mu[i][2];

            mu[i][0] = x_tmp;
            mu[i][1] = y_tmp;
            mu[i][2] = z_tmp;

            // Imaginary part
            x_tmp = mu[i][3] * cos_angle - mu[i][4] * sin_angle;
            y_tmp = mu[i][3] * sin_angle + mu[i][4] * cos_angle;
            z_tmp = mu[i][5];

            mu[i][3] = x_tmp;
            mu[i][4] = y_tmp;
            mu[i][5] = z_tmp;
        }
    }
    else {
        throw std::runtime_error("Unknown rotation axys: " + axis);
    }    
}
//----------------------------------------------------------------------