#ifndef STRING_MANIPULATION_HPP
#define STRING_MANIPULATION_HPP

#include <string>
#include <iostream>
#include <istream>
#include <stdexcept>
#include <string_view>

//----------------------------------------------------------------------
/// @class String_manipulation
/// @brief Provides utility functions for converting strings to numerical values
///        and scanning input streams.
class String_manipulation
{
public:
    /// @brief Initializes internal state.
    void initialize();

    /// @brief Converts a string to a double-precision floating point number.
    /// @param str Input string to convert.
    /// @param out Output variable to store the parsed double value.
    void string_to_float(const std::string &str, double &out);

    /// @brief Converts a string to an integer.
    /// @param str Input string to convert.
    /// @param out Output variable to store the parsed integer value.
    void string_to_int(const std::string &str, int &out);

    /// @brief Checks if the provided string is one of the accepted entries for the "what" keyword.
    /// @param str Input string to check.
    /// @param out Output variable to store the validated "what" value.
    void string_what_accepted_entries(const std::string &str, std::string &out);

    /// @brief Checks if the provided string is one of the accepted entries for the "parametrization" keyword.
    /// @param str Input string to check.
    /// @param out Output variable to store the validated "parametrization" value.
    void string_parametrization_accepted_entries(const std::string &str, std::string &out);
};
//----------------------------------------------------------------------
/// @brief Moves the stream cursor to the line immediately after a specified marker string.
/// 
/// This function scans each line of the input stream `in`, looking for a line
/// that exactly matches `marker`. If found, the stream will be positioned on the next line.
/// If not found, the function returns false and leaves the stream at EOF.
///
/// @param in Input stream to scan.
/// @param marker String to match against a line in the stream.
/// @return true if marker is found and the stream is positioned after it, false otherwise.
inline bool go_to_string(std::istream &in, std::string_view marker)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r')
            line.pop_back(); // CRLF-safe
        if (line == marker)
            return true; // positioned after marker line
    }
    return false;
}

#endif // STRING_MANIPULATION_HPP
//----------------------------------------------------------------------
