#ifndef INPUT_HPP
#define INPUT_HPP

#include "output.hpp"
#include "string_manipulation.hpp"
#include "target.hpp"

#include <string>

//----------------------------------------------------------------------
/// @class Input
/// @brief Input module for handling user input and input file parsing.
///
/// @details
/// The Input class manages user arguments, input file location,
/// and extraction of user-defined options. It is responsible for:
/// - parsing command-line arguments,
/// - validating and reading the input file,
/// - setting paths and values in the @ref Target object.
/// It also stores internal utilities (e.g., string parser).
class Input {
public:
    /// @brief Default constructor (initializes default filename to "input.inp").
    Input();

    /// @brief Parses command-line arguments to override input file name.
    /// @param argc Argument count from main().
    /// @param argv Argument vector from main().
    /// @param out Output handler used for logging.
    /// @param target Target object to store parsed values.
    void get_arguments(int argc, char* argv[], Output& out, Target& target);

    /// @brief Checks input file existence and extension.
    /// @param out Output stream for logging (default is std::cout).
    void check_input_file(const Output& out);

    /// @brief Prints input file header information and configuration summary.
    /// @param out Output object for printing.
    /// @param target Parsed Target object to display.
    void print_input_info(const Output& out, const Target& target);

    /// @brief Reads and parses the input file, populating the Target.
    /// @param target Target object to populate with parsed options.
    void read(Target& target); 

    /// @brief Determines the target calculation based on parsed input.
    /// @param target Target object to configure.
    void get_target(Target& target);

    std::string input_filename;  ///< Filename of the input file (default: "input.inp")

    Output out;                  ///< Output handler for logging
    Target target;               ///< Target object for simulation configuration
    String_manipulation str_manipulation; ///< Utility for parsing input strings

private:

    /// @brief Parses command-line arguments and sets the input filename.
    /// @param argc Argument count from main().
    /// @param argv Argument vector from main().
    /// @param out Output handler for error reporting.
    void parse_arguments(int argc, char* argv[], Output& out);

    /// @brief Checks if a file exists and is accessible.
    /// @param path Absolute or relative path to the file.
    /// @throws std::runtime_error if file is not found or not readable.
    void file_exists(const std::string& path) const;

    /// @brief Resolves relative path to an absolute one based on the input file's location.
    /// @param relative_path Path specified in the input file.
    /// @return Absolute normalized path string.
    std::string resolve_relative_to_input(const std::string& relative_path) const;

    /// @brief Validates and stores file input values from the input file.
    /// @param raw_input Filename or path as specified by user.
    /// @param input_field Reference to variable for raw (user-given) input.
    /// @param resolved_field Reference to store resolved absolute path.
    void check_and_store_file(const std::string& raw_input,
                              std::string& input_field,
                              std::string& resolved_field) const;

    /// @brief Parses and stores a transition dipole from the input file.
    /// @param raw_input String line containing dipole vector components.
    /// @param transdip Reference to array that will store the parsed vector.
    void check_and_store_transition_dipole(const std::string& raw_input,
                                           std::array<double, 
                                           3>& transdip) const;
};

#endif // INPUT_HPP

//----------------------------------------------------------------------
