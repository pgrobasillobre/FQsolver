#include "solvent.hpp"
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
//void Nanoparticle::read_nanoparticle(const Target &target, const Output &out)
//{
//  // Check file existance.
//  std::string filepath = target.nanoparticle_file;
//
//  std::ifstream infile(filepath);
//  if (!infile)
//  {
//    throw std::runtime_error("File: " + filepath + "not found.");
//  }
//
//  // Check if FRET quantities are present
//  bool is_fret_quantities_present = false;
//  is_fret_quantities_present = go_to_string(infile, Parameters::fret_start);
//
//  // Check if FRET end marker is present
//  bool is_fret_end_marker_present = false;
//  is_fret_end_marker_present = go_to_string(infile, Parameters::fret_end);
//
//  // Rewind and come back to fret_start
//  infile.clear();
//  infile.seekg(0, std::ios::beg);
//  is_fret_quantities_present = go_to_string(infile, Parameters::fret_start);
//
//  // Check charges / charges + dipoles
//  std::string line;
//  std::getline(infile, line);
//  if (line == Parameters::charges_header)
//  {
//    charges = true;
//    nanoparticle_model = "charges";
//  }
//  else if (line == Parameters::charges_and_dipoles_header)
//  {
//    charges_and_dipoles = true;
//    nanoparticle_model = "charges + dipoles";
//  }
//  else
//  {
//    throw std::runtime_error("Expected header line with charges or charges and dipoles, got: " + line);
//  }
//
//  // Read charges / charges and dipoles.
//  while (std::getline(infile, line))
//  {
//    if (line == Parameters::fret_end)
//      break;
//
//    std::istringstream iss(line);
//    std::array<double, 3> coords{};
//
//    if (charges)
//    {
//      std::array<double, 2> charge{};
//      iss >> charge[0] >> charge[1] >> coords[0] >> coords[1] >> coords[2];
//      q.push_back(charge);
//    }
//    else if (charges_and_dipoles)
//    {
//      std::array<double, 2> charge{};
//      std::array<double, 6> dipole{};
//      iss >> charge[0] >> charge[1] >> dipole[0] >> dipole[1] >> dipole[2] >> dipole[3] >> dipole[4] >> dipole[5] >> coords[0] >> coords[1] >> coords[2];
//      q.push_back(charge);
//      mu.push_back(dipole);
//    }
//
//    xyz.push_back(coords);
//  }
//  natoms = xyz.size();
//
//  // Compute geometrical center
//  for (const auto &coord : xyz)
//  {
//    geom_center[0] += coord[0];
//    geom_center[1] += coord[1];
//    geom_center[2] += coord[2];
//  }
//  geom_center[0] /= natoms;
//  geom_center[1] /= natoms;
//  geom_center[2] /= natoms;
//}
//----------------------------------------------------------------------
