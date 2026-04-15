#include "solvent.hpp"
#include "output.hpp"
#include "target.hpp"
#include "parameters.hpp"
#include "string_manipulation.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

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
// Check solvent file extension
void Solvent::check_solvent_file_extension(const std::string &str, const Target &target)
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
        "'" + str + "' does not have a valid solvent file extension. Accepted values are: " +
        accepted.str() + ".\n");
  }

  // Store the extension for later use
  solvent_file_extension = extension;

  // If file extension is pdb, "group" and "read atoms" filters must be present in the input file, otherwise raise an error.
  if (extension == ".pdb")
  {
    if (!target.is_group_present)
    {
      throw std::runtime_error(
          "Solvent file has .pdb extension, but no 'group' keyword found in the input file. Please specify a residue/group name to read from the PDB solvent file using the 'group' keyword.\n");
    }

    if (!target.is_read_atoms_present)
    {
      throw std::runtime_error(
          "Solvent file has .pdb extension, but no 'read atoms' keyword found in the input file. Please specify atom names to read from the PDB solvent file using the 'read atoms' keyword.\n");
    }
  }
}
//----------------------------------------------------------------------
std::vector<std::string> split_whitespace(const std::string &line)
{
  std::istringstream stream(line);
  std::vector<std::string> tokens;
  std::string token;
  while (stream >> token)
  {
    tokens.push_back(token);
  }

  return tokens;
}

bool parse_double(const std::string &token, double &value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(token, &parsed);
    return parsed == token.size();
  }
  catch (const std::exception &)
  {
    return false;
  }
}

bool parse_int(const std::string &token, int &value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stoi(token, &parsed);
    return parsed == token.size();
  }
  catch (const std::exception &)
  {
    return false;
  }
}

std::string uppercase_copy(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c)
                 { return std::toupper(c); });
  return value;
}

//----------------------------------------------------------------------
void Solvent::read_solvent_geometry(const std::string &filepath, const Target &target)
{

  // Read file according to extension.
  if (solvent_file_extension == ".xyz")
  {
    read_solvent_geometry_xyz(filepath);
  }
  else if (solvent_file_extension == ".pdb")
  {
    read_solvent_geometry_pdb(filepath, target);
  }
  else
  {
    throw std::runtime_error("Unsupported solvent file extension: " + solvent_file_extension);
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
    std::string label;            // Variable to hold the atomic label.
    std::string trailing;         // Variable to check for any extra content after the expected fields.

    std::array<double, 3> coords = {0.0, 0.0, 0.0};

    if (!(iss >> label >> coords[0] >> coords[1] >> coords[2]) || (iss >> trailing))
    {
      throw std::runtime_error(
          "File \"" + filepath + "\" is corrupted: invalid XYZ coordinate line \"" + line + "\".");
    }

    // Convert coordinates from Angstroms to Bohr
    coords[0] *= Parameters::ToBohr;
    coords[1] *= Parameters::ToBohr;
    coords[2] *= Parameters::ToBohr;

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
// Read solvent geometry from a PDB file.
void Solvent::read_solvent_geometry_pdb(const std::string &filepath, const Target &target)
{
  std::ifstream infile(filepath);
  std::string line;

  atomic_label.clear();
  xyz.clear();
  MolIndex.clear();

  std::vector<std::string> read_atoms_upper;
  read_atoms_upper.reserve(target.read_atoms.size());
  MolIndex.reserve(target.read_atoms.size());

  for (const auto &atom : target.read_atoms)
  {
    read_atoms_upper.push_back(uppercase_copy(atom));
  }

  const std::string read_group_upper = uppercase_copy(target.read_group);
  int previous_molecule_id = 0;
  bool has_previous_molecule_id = false;
  nmol = 0;

  while (std::getline(infile, line))
  {
    if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0)
    {
      continue;
    }

    const std::vector<std::string> tokens = split_whitespace(line);
    if (tokens.size() < 4)
    {
      throw std::runtime_error("File \"" + filepath + "\" is corrupted: invalid PDB coordinate line \"" + line + "\".");
    }

    // PDB-like whitespace columns: atom name = 3. The configured group may span columns.
    const std::string atom_name = tokens[2];

    if (!read_atoms_upper.empty() &&
        std::find(read_atoms_upper.begin(), read_atoms_upper.end(), uppercase_copy(atom_name)) == read_atoms_upper.end())
    {
      continue;
    }

    const std::string line_upper = uppercase_copy(line);
    const std::size_t group_pos = line_upper.find(read_group_upper);
    if (group_pos == std::string::npos)
    {
      continue;
    }

    const std::string after_group = line.substr(group_pos + target.read_group.size());
    const std::vector<std::string> after_group_tokens = split_whitespace(after_group);

    int molecule_id = 0;
    std::size_t molecule_id_index = after_group_tokens.size();
    for (std::size_t i = 0; i < after_group_tokens.size(); ++i)
    {
      if (parse_int(after_group_tokens[i], molecule_id))
      {
        molecule_id_index = i;
        break;
      }
    }

    if (molecule_id_index == after_group_tokens.size() ||
        molecule_id_index + 3 >= after_group_tokens.size())
    {
      throw std::runtime_error("File \"" + filepath + "\" is corrupted: invalid PDB coordinate line \"" + line + "\".");
    }

    std::array<double, 3> coords = {0.0, 0.0, 0.0};
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!parse_double(after_group_tokens[molecule_id_index + 1], x) ||
        !parse_double(after_group_tokens[molecule_id_index + 2], y) ||
        !parse_double(after_group_tokens[molecule_id_index + 3], z))
    {
      throw std::runtime_error("File \"" + filepath + "\" is corrupted: invalid PDB coordinate line \"" + line + "\".");
    }

    if (!has_previous_molecule_id || molecule_id != previous_molecule_id)
    {
      ++nmol;
      previous_molecule_id = molecule_id;
      has_previous_molecule_id = true;
    }

    coords[0] = x * Parameters::ToBohr;
    coords[1] = y * Parameters::ToBohr;
    coords[2] = z * Parameters::ToBohr;

    MolIndex.push_back({nmol});
    atomic_label.push_back(atom_name);
    xyz.push_back(coords);
  }

  natoms = static_cast<int>(xyz.size());
  if (natoms == 0)
  {
    throw std::runtime_error("No PDB atoms matched the requested solvent filters in file \"" + filepath + "\".");
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
  check_solvent_file_extension(filepath, target);

  // Check file with given extension is not corrupted (e.g., missing lines, wrong format).
  read_solvent_geometry(filepath, target);
}
//----------------------------------------------------------------------
