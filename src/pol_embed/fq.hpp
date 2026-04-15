#ifndef FQ_HPP
#define FQ_HPP

//----------------------------------------------------------------------
/// @class FQ
/// @brief Implements FQ matrices and charges calculation routines.
///

class Solvent;

class FQ
{
public:
  // ---- API ----
  /// @brief Computes FQ charges from the potential and solvent parameters.
  /// @param solv   Solvent system containing necessary parameters for FQ charge calculation.
  void calc_charges(Solvent &solv);

private:
  /// @brief Computes the Tqq matrix (charge-charge interaction tensor) for FQ charge calculation.
  /// @param solv   Solvent system containing necessary parameters for Tqq matrix calculation.
  void calc_tqq(Solvent &solv);
};

#endif // FQ_HPP
//----------------------------------------------------------------------
