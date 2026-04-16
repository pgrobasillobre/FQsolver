#ifndef FQ_HPP
#define FQ_HPP

#include <vector>

//----------------------------------------------------------------------
/// @class FQ
/// @brief Implements FQ matrices and charges calculation routines.
///

class Solvent;
class Output;
struct Target;

class FQ
{
public:
  // Arrays
  std::vector<std::vector<double>> tempTqq; // Temporary storage for Tqq tensor components.
  std::vector<double> Dmat;   // FQ LHS matrix lower triangular part stored in a 1D array

  int norder = 0;  ///< Order of FQ Dmatrix
  int nlength = 0; ///< Length of to store FQ matrix lower triangular part in a 1D array

  // ---- API ----
  /// @brief Computes FQ charges from the potential and solvent parameters.
  /// @param solv   Solvent system containing necessary parameters for FQ charge calculation.
  /// @param target Target system containing debug settings.
  /// @param out    Output object used for debug printing.
  void calc_charges(Solvent &solv, const Target &target, const Output &out);

private:
  /// @brief Computes the Tqq matrix (charge-charge interaction tensor) for FQ charge calculation.
  /// @param solv   Solvent system containing necessary parameters for Tqq matrix calculation.
  /// @param target Target system containing debug settings.
  /// @param out    Output object used for debug printing.
  void calc_tqq(Solvent &solv, const Target &target, const Output &out);

  /// @brief  Calculates the FQ LHS DMatrix
  /// @param solv  Solvent system containing necessary parameters for DMatrix calculation.
  /// @param target Target system containing debug settings.
  /// @param out  Output object used for debug printing.
  void calc_dmat(Solvent &solv, const Target &target, const Output &out);
};

#endif // FQ_HPP
//----------------------------------------------------------------------
