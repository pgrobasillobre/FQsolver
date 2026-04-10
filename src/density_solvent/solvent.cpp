#include "solvent.hpp"
#include "output.hpp"
#include "target.hpp"
#include "parameters.hpp"
#include "string_manipulation.hpp"

#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>

//----------------------------------------------------------------------
// Check file exists
void Solvent::check_file_exists(const std::string &filepath) const
{
  std::ifstream file(filepath);
  if (!file)
  {
    throw std::runtime_error("File: " + filepath + " not found.");
  }
}

//----------------------------------------------------------------------
// Check file extension is .xyz
void Solvent::check_solvent_file_extension(const std::string &str)
{
  std::string extension = str.substr(str.find_last_of('.'));

  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char c)
                 { return std::tolower(c); });

  const auto it = std::find(Parameters::accepted_solvent_file_extensions.begin(),
                            Parameters::accepted_solvent_file_extensions.end(),
                            extension);

  if (it == Parameters::accepted_solvent_file_extensions.end())
  {
    std::ostringstream accepted;
    for (std::size_t i = 0; i < Parameters::accepted_solvent_file_extensions.size(); ++i)
    {
      if (i > 0)
        accepted << ", ";
      accepted << Parameters::accepted_solvent_file_extensions[i];
    }

    throw std::runtime_error(
        "Error: '" + str + "' does not have a valid solvent file extension. Accepted values are: " +
        accepted.str() + ".\n");
  }

  // Store the extension for later use
  solvent_file_extension = extension;
}
//----------------------------------------------------------------------
// Read solvent geometry file
void Solvent::read_solvent_geometry(const std::string &filepath)
{

  // Read file according to extension (currently only .xyz is supported)
  if (solvent_file_extension == ".xyz")
  {
    read_solvent_geometry_xyz(filepath);
  }

}
//----------------------------------------------------------------------
// Read solvent geometry from an XYZ file
void Solvent::read_solvent_geometry_xyz(const std::string &filepath)
{
  std::string line;

  // Open the file
  std::ifstream infile(filepath);

  // First line: number of atoms.
  // Check the lines is present, and contains a single integer, otherwise raise an error.
  if (!std::getline(infile, line))
  {
    throw std::runtime_error("File \"" + filepath + "\" is empty or corrupted.");
  }

  {
    std::istringstream iss(line);
    int parsed_natoms = 0;
    std::string trailing;

    if (!(iss >> parsed_natoms) || (iss >> trailing) || parsed_natoms < 0)
    {
      throw std::runtime_error(
          "File \"" + filepath + "\" corrupted: first line must contain a single integer with the number of atoms.");
    }

    natoms = parsed_natoms;
  }

  // Second line: title/comment.
  // Check it's present, otherwise raise an error. 
  // We will not use it, but we want to ensure the file is not corrupted (e.g., missing lines).
  if (!std::getline(infile, line))
  {
    throw std::runtime_error(
        "File \"" + filepath + "\" is corrupted: missing title line after the number of atoms.");
  }

  atomic_label.clear();
  xyz.clear();

  // Reserve memory for atomic labels and coordinates based on the number of atoms, to improve performance and avoid unnecessary reallocations.
  atomic_label.reserve(natoms);
  xyz.reserve(natoms);

  // Read atomic coordinates
  for (int i = 0; i < natoms; ++i)
  {
    // Raise error if line is missing. 
    // Check the line contains an atomic label followed by three floating-point numbers, and nothing else.
    if (!std::getline(infile, line))
    {
      throw std::runtime_error(
          "File \"" + filepath + "\" is corrupted: expected " + std::to_string(natoms) +
          " coordinate lines, but found only " + std::to_string(i) + ".");
    }

    std::istringstream iss(line); // Create a string stream to parse the line.
    std::string label; // Variable to hold the atomic label.
    std::string trailing; // Variable to check for any extra content after the expected fields.

    std::array<double, 3> coords = {0.0, 0.0, 0.0};

    if (!(iss >> label >> coords[0] >> coords[1] >> coords[2]) || (iss >> trailing))
    {
      throw std::runtime_error(
          "File \"" + filepath + "\" is corrupted: invalid XYZ coordinate line \"" + line + "\".");
    }

    atomic_label.push_back(label);
    xyz.push_back(coords);
  }

  // Final check: ensure we read exactly the number of atoms specified in the first line.
  if (static_cast<int>(xyz.size()) != natoms)
  {
    throw std::runtime_error(
        "File \"" + filepath + "\" is corrupted: number of coordinates read does not match the number of atoms.");
  }
}
//----------------------------------------------------------------------
// Load solvent data: coordinates.
void Solvent::read_solvent(const Target &target, const Output &out)
{
  // Check file existance.
  std::string filepath = target.solvent_file;
  check_file_exists(filepath);

  // Check file extension.
  check_solvent_file_extension(filepath);

  // Check file with given extension is not corrupted (e.g., missing lines, wrong format).
  read_solvent_geometry(filepath);

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
}
//----------------------------------------------------------------------
