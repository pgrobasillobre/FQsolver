#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "density.hpp"
#include "integrals.hpp"
#include "solvent.hpp"

#include <optional>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <string>
#include <ostream>
#include <vector>

class Input; // Forward declaration of Input to avoid circular dependency

//----------------------------------------------------------------------
/// @class Output
/// @brief Manages output filenames, streams, and formatted printing.
///
/// Handles writing of results, densities, and diagnostics. Also provides
/// access to the main output stream.
class Output
{
public:
    /// @brief Constructor.
    Output();

    /// @brief Generates the output filename based on the input filename.
    /// @param in_file The input filename provided by the user.
    void out_file_fill(const std::string &in_file);

    /// @brief Opens the output file for writing.
    /// @throws std::runtime_error if the file cannot be opened.
    void open();

    /// @brief Closes the output file stream.
    void close();

    /// @brief Prints a standard FQSolver banner to the output stream.
    void print_banner();

    /// @brief Prints cube density data to the output stream.
    /// @param target The current target object containing input parameters.
    /// @param cube The Density object containing data to print.
    /// @param header Optional custom header for the density section.
    void print_density(const Target &target,
                       const Density &cube,
                       std::optional<std::string> header = std::nullopt);

    /// @brief Prints solvent data to the output stream.
    /// @param target The current target object containing input parameters.
    /// @param solv The Solvent object containing data to print.
    void print_solvent(const Target &target, const Solvent &solv);

    /// @brief Prints results of computed integrals.
    /// @param target Target configuration.
    /// @param integrals Integrals object containing computed values.
    void print_results_integrals(const Target &target, const Integrals &integrals);

    /// @brief Prints results of computed potentials and fields.
    /// @param target Target configuration.
    /// @param solv Solvent object containing solvent data.
    void print_results_pot_fld(const Target &target, const Solvent &solv);

    /// @brief Outputs the coordinates of the cube density.
    /// @param what_dens Label (e.g., "donor", "acceptor").
    /// @param n_points Number of grid points.
    /// @param xyz Vector of XYZ coordinates.
    void print_cube_coordinates(const std::string what_dens,
                                const int n_points,
                                const std::vector<std::array<double, 3>> &xyz) const;

    /// @brief Prints a matrix in 5-column blocks.
    /// @param title Section title printed above the matrix.
    /// @param matrix Matrix to print.
    void print_matrix(const std::string &title, const std::vector<std::vector<double>> &matrix) const;

    /// @brief Prints a packed triangular matrix by temporarily reconstructing its full dense form.
    /// @param title Section title printed above the matrix.
    /// @param packed_matrix Packed triangular matrix.
    /// @param order Order of the reconstructed square matrix.
    /// @param triangle Stored triangle: "L" for lower or "U" for upper.
    void print_matrix(const std::string &title,
                      const std::vector<double> &packed_matrix,
                      int order,
                      const std::string &triangle) const;

    /// @brief Prints the FQ RHS vector together with the quantities used to build it.
    /// @param title Section title printed above the matrix.
    /// @param solv Solvent object containing atom data, potentials, and FQ parameters.
    /// @param rhs RHS vector to print.
    void print_matrix_rhs(const std::string &title, const Solvent &solv, const std::vector<double> &rhs) const;

    /// @brief Prints FQ charges by atom.
    /// @param title Section title printed above the charge table.
    /// @param solv Solvent object containing atom labels and molecule indices.
    /// @param results Solution vector containing charges followed by Lagrange multipliers.
    void print_results(const std::string &title, const Solvent &solv, const std::vector<double> &results) const;

    /// @brief Writes FQ charges to the FQSolver_results directory.
    /// @param target Target configuration.
    /// @param solv Solvent object containing solvent data, potentials, and molecule indices.
    /// @param results Solution vector containing charges followed by Lagrange multipliers.
    void print_results(const Target &target, const Solvent &solv, const std::vector<double> &results) const;

    /// @brief Prints nanoparticle coordinates and dipoles, if present.
    /// @param infile Source file name.
    /// @param np Nanoparticle object.
    // void print_np_coords_dipoles(const std::string infile, const Nanoparticle& np) const;

    /// @brief Horizontal separator (80 dashes) used in reports.
    const std::string sticks = std::string(80, '-');

    /// @brief Returns a reference to the output stream.
    /// @return The file stream used for output.
    std::ofstream &stream() const;

    /// @brief Full path or name of the output file.
    std::string output_filename;

private:
    /// @brief Prints a density point: index and XYZ.
    /// @param out The output stream to write to.
    /// @param i Index of the point.
    /// @param a X-coordinate.
    /// @param b Y-coordinate.
    /// @param c Z-coordinate.
    void print_formatted_line1(std::ostream &out, int i, double a, double b, double c);

    /// @brief Prints atom name and coordinates (for cube files).
    /// @param out Output stream.
    /// @param atom Atom label (e.g., "C").
    /// @param x X-coordinate.
    /// @param y Y-coordinate.
    /// @param z Z-coordinate.
    void print_formatted_line2(std::ostream &out, const std::string atom, double x, double y, double z) const;

    /// @brief Prints atom name and coordinates (for nanoparticles).
    /// @param out Output stream.
    /// @param atom Atom label (e.g., "C").
    /// @param x X-coordinate.
    /// @param y Y-coordinate.
    /// @param z Z-coordinate.
    void print_formatted_line3(std::ostream &out, const std::string atom, double x, double y, double z);


    /// @brief Prints atom name, coordinates, and molecule index (for solvents loaded from PDB).
    /// @param out Output stream.
    /// @param atom Atom label (e.g., "C").
    /// @param x X-coordinate.
    /// @param y Y-coordinate.
    /// @param z Z-coordinate.
    /// @param mol_index Index of the molecule the atom belongs to.
    void print_formatted_line4(std::ostream &out, const std::string atom, double x, double y, double z, int mol_index);

    /// @brief Format string for printing index and coordinates.
    std::string format1 = "   {:5d} {:15.7E} {:15.7E} {:15.7E}\n";

    /// @brief File stream used for output (mutable to allow usage in const methods).
    mutable std::ofstream log_stream;
};

#endif // OUTPUT_HPP
//----------------------------------------------------------------------
