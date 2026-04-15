#ifndef FQ_HPP
#define FQ_HPP

//----------------------------------------------------------------------
/// @class FQ
/// @brief Implements FQ matrices and charges calculation routines.
///

struct Target;
class Solvent;

class FQ
{
public:
  // ---- API ----
  /// @brief Computes FQ charges from the potential and solvent parameters.
  /// @param target Target system containing necessary input for FQ charge calculation.
  /// @param solv   Solvent system containing necessary parameters for FQ charge calculation.
  void calc_charges(Target &target, Solvent &solv);
};

#endif // FQ_HPP
//----------------------------------------------------------------------
