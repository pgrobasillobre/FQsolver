#include "string_manipulation.hpp"
#include "parameters.hpp"

#include <string>
#include <iostream>
#include <stdexcept>
#include <istream>
#include <string_view>
#include <algorithm>
#include <sstream>


//----------------------------------------------------------------------
// Converts a string to a double-precision floating point number.
void String_manipulation::string_to_float(const std::string &str, double &out)
{
    try
    {
        out = std::stod(str);
    }
    catch (const std::invalid_argument &)
    {
        throw std::runtime_error("Error: '" + str + "' is not a valid float.\n");
    }
    catch (const std::out_of_range &)
    {
        throw std::runtime_error("Error: '" + str + "' is out of range for a float.\n");
    }
}
//----------------------------------------------------------------------
// Moves the stream cursor to the line immediately after a specified marker string.
void String_manipulation::string_to_int(const std::string &str, int &out)
{
    try
    {
        out = std::stoi(str);
    }
    catch (const std::invalid_argument &)
    {
        throw std::runtime_error("Error: '" + str + "' is not a valid integer.\n");
    }
    catch (const std::out_of_range &)
    {
        throw std::runtime_error("Error: '" + str + "' is out of range for an integer.\n");
    }
}
//----------------------------------------------------------------------
// Checks if the provided string is one of the accepted entries for the "what" keyword.
void String_manipulation::string_what_accepted_entries(const std::string &str, std::string &out)
{
    out = str;

    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    const auto it = std::find(Parameters::accepted_what_entries.begin(),
                              Parameters::accepted_what_entries.end(),
                              out);

    if (it == Parameters::accepted_what_entries.end())
    {
        std::ostringstream accepted;
        for (std::size_t i = 0; i < Parameters::accepted_what_entries.size(); ++i)
        {
            if (i > 0)
                accepted << ", ";
            accepted << Parameters::accepted_what_entries[i];
        }

        throw std::runtime_error(
            "Error: '" + str + "' is not a valid entry for 'what'. Accepted values are: " +
            accepted.str() + ".\n");
    }
}
//----------------------------------------------------------------------
