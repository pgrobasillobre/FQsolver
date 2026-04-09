#ifndef ENUMS_HPP
#define ENUMS_HPP

/// @brief Defines calculation modes for the program
enum class TargetMode {
    None,             ///< No calculation selected
    IntegrateCube,    ///< Electron density integration
    Solute_Solvent_Pot_Field, ///< Solute density potential/field at solvent positions
};

#endif // ENUMS_HPP

