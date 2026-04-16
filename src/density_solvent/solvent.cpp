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
          "Solvent file has .pdb extension, but no 'group' keyword found in the input file. \n\n Please specify a residue/group name to read from the PDB solvent file using the 'group' keyword.\n");
    }

    if (!target.is_read_atoms_present)
    {
      throw std::runtime_error(
          "Solvent file has .pdb extension, but no 'read atoms' keyword found in the input file. \n\n Please specify atom names to read from the PDB solvent file using the 'read atoms' keyword.\n");
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

  // Assign parameters for FQ charge calculation based on atom types.
  if (target.is_parametrization_present)
  {
    assign_solvent_parameters(target.fq_parametrization);
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
  atomsToMol.clear();
  indexToMol.clear();
  atomsToIndex.clear();

  std::vector<std::string> read_atoms_lower;
  read_atoms_lower.reserve(target.read_atoms.size());
  MolIndex.reserve(target.read_atoms.size());
  atomsToMol.reserve(target.read_atoms.size());

  for (const auto &atom : target.read_atoms)
  {
    std::string atom_lower = atom;
    std::transform(atom_lower.begin(), atom_lower.end(), atom_lower.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    read_atoms_lower.push_back(atom_lower);
  }

  std::string read_group_lower = target.read_group;
  std::transform(read_group_lower.begin(), read_group_lower.end(), read_group_lower.begin(),
                 [](unsigned char c)
                 { return std::tolower(c); });
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

    std::string atom_name_lower = atom_name;
    std::transform(atom_name_lower.begin(), atom_name_lower.end(), atom_name_lower.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    if (!read_atoms_lower.empty() &&
        std::find(read_atoms_lower.begin(), read_atoms_lower.end(), atom_name_lower) == read_atoms_lower.end())
    {
      continue;
    }

    std::string line_lower = line;
    std::transform(line_lower.begin(), line_lower.end(), line_lower.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    const std::size_t group_pos = line_lower.find(read_group_lower);
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

    coords[0] = x * Parameters::ToBohr;
    coords[1] = y * Parameters::ToBohr;
    coords[2] = z * Parameters::ToBohr;

    atomsToMol.push_back(molecule_id);
    atomic_label.push_back(atom_name);
    xyz.push_back(coords);
  }

  natoms = static_cast<int>(xyz.size());
  if (natoms == 0)
  {
    throw std::runtime_error("No PDB atoms matched the requested solvent filters in file \"" + filepath + "\".");
  }

  assign_molecule_indices();

  // Store MolCharge from target in solvent object
  MolCharge = target.MolCharge;
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

  // Store the FQ kernel selected in the input file.
  fq_kernel = target.fq_kernel;
}
// ----------------------------------------------------------------------
void Solvent::assign_solvent_parameters(const std::string &parametrization)
{
  const auto param_set = Parameters::fq_params.find(parametrization);
  if (param_set == Parameters::fq_params.end())
  {
    throw std::runtime_error("Unsupported FQ parametrization: " + parametrization);
  }

  const auto label_map = Parameters::fq_label_to_type.find(parametrization);
  if (label_map == Parameters::fq_label_to_type.end())
  {
    throw std::runtime_error("Missing atom label map for FQ parametrization: " + parametrization);
  }

  typeIndex.assign(natoms, -1);
  typeName.clear();
  typeChi.clear();
  typeEta.clear();

  std::unordered_map<std::string, int> type_name_to_index;

  for (int i = 0; i < natoms; ++i)
  {
    const std::string &label = atomic_label[i];
    const auto type_name = label_map->second.find(label);
    if (type_name == label_map->second.end())
    {
      throw std::runtime_error(
          "Unsupported atom label \"" + label + "\" for FQ parametrization \"" + parametrization + "\".\n"
                                                                                                   "\n"
                                                                                                   "Supported atom labels are:\n"
                                                                                                   "  - Oxygen:  OW, O\n"
                                                                                                   "  - Hydrogen: HW1, HW2, H1, H2, H\n"
                                                                                                   "\n"
                                                                                                   "Atom labels are case-sensitive.");
    }

    const auto params = param_set->second.find(type_name->second);
    if (params == param_set->second.end())
    {
      throw std::runtime_error(
          "Missing FQ parameters for atom type \"" + type_name->second +
          "\" in parametrization \"" + parametrization + "\".");
    }

    auto type_index = type_name_to_index.find(type_name->second);
    if (type_index == type_name_to_index.end())
    {
      const int new_type_index = static_cast<int>(typeName.size());
      type_name_to_index[type_name->second] = new_type_index;
      typeName.push_back(type_name->second);
      typeChi.push_back(params->second.chi);
      typeEta.push_back(params->second.eta);
      typeIndex[i] = new_type_index;
    }
    else
    {
      typeIndex[i] = type_index->second;
    }
  }

  // Final check: ensure all atoms were assigned a valid type index.
  for (int i = 0; i < natoms; ++i)
  {
    if (typeIndex[i] == -1)
    {
      throw std::runtime_error(
          "Failed to assign FQ parameters for atom \"" + atomic_label[i] + "\" at index " + std::to_string(i) +
          ". This should not happen if all atom labels are supported and parameters are provided, but something went wrong during parameter assignment.");
    }
  }

  // Store the number of unique atom types.
  ntypes = static_cast<int>(typeName.size());

  // Compute Gaussian widths for charge distribution based on the FQ parameters, using the formula Rq = sqrt(2/pi)/eta.
  typeRq.assign(ntypes, 0.0);
  for (int i = 0; i < ntypes; ++i)
  {
    typeRq[i] = std::sqrt(2.0 / Parameters::pi) / typeEta[i];
  }
}
//----------------------------------------------------------------------
// Builds compact molecule indices from original molecule labels.
void Solvent::assign_molecule_indices()
{
  indexToMol.clear();
  atomsToIndex.assign(atomsToMol.size(), -1);

  for (const int molecule_label : atomsToMol)
  {
    if (std::find(indexToMol.begin(), indexToMol.end(), molecule_label) == indexToMol.end())
    {
      indexToMol.push_back(molecule_label);
    }
  }

  for (std::size_t i = 0; i < atomsToMol.size(); ++i)
  {
    const auto molecule_index = std::find(indexToMol.begin(), indexToMol.end(), atomsToMol[i]);
    if (molecule_index == indexToMol.end())
    {
      throw std::runtime_error("Failed to assign compact molecule index for solvent atom.");
    }

    atomsToIndex[i] = static_cast<int>(std::distance(indexToMol.begin(), molecule_index));
  }

  nmol = static_cast<int>(indexToMol.size());

  MolIndex.clear();
  MolIndex.reserve(atomsToIndex.size());
  for (const int compact_index : atomsToIndex)
  {
    MolIndex.push_back({compact_index + 1});
  }
}