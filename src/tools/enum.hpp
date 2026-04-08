#ifndef ENUMS_HPP
#define ENUMS_HPP

/// @brief Defines calculation modes for the program
enum class TargetMode {
    None,             ///< No calculation selected
    IntegrateCube,    ///< Electron density integration
    Acceptor_Donor,   ///< Acceptor-Donor coulomb (optional: overlap integral)
    Acceptor_NP,      ///< Acceptor + Nanoparticle coulomb
    Acceptor_NP_Donor ///< Acceptor-Donor + Nanoparticle coulomb (optional: overlap integral)
};

#endif // ENUMS_HPP

