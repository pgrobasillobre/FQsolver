#include "string_manipulation.hpp"

#include <string>
#include <iostream>
#include <stdexcept>
#include <istream>
#include <string_view>

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
