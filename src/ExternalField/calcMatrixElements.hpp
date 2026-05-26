#pragma once
#include "CorePolarisation.hpp"
#include "Coulomb/meTable.hpp"
#include "MBPT/StructureRad.hpp"
#include "qip/String.hpp"
#include <iostream>
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class HartreeFock;
}

namespace ExternalField {

/*!
  @brief Small struct to store a single calculated reduced matrix element with RPA correction.
  @details
  Holds the result of a single matrix element calculation \f$ \redmatel{a}{h}{b} \f$,
  including the lowest-order value, the RPA correction, and the transition
  frequency.

  The printed output includes columns: state labels, frequency, lowest-order
  matrix element, and (if non-zero) the RPA-corrected value \f$ t_{ab} + \delta V_{ab} \f$.
*/
struct MEdata {

  std::string a, b;
  double w_ab;
  double hab, dv;

  static std::string title(bool rpaQ = true) {
    if (rpaQ)
      return "    a    b   w_ab      t_ab           RPA_ab";
    else
      return "    a    b   w_ab      t_ab";
  }
  static std::string title_noRPA() { return "    a    b   w_ab      t_ab"; }

  friend std::ostream &operator<<(std::ostream &os, const MEdata &m) {
    os << qip::fstring(" %4s %4s  %8.5f  %13.6e", m.a.c_str(), m.b.c_str(),
                       m.w_ab, m.hab);
    if (m.dv != 0.0) {
      os << qip::fstring("  %13.6e", m.hab + m.dv);
    }
    return os;
  }
};

/*!
  @brief Calculates reduced matrix elements <a||h||b> for lists of orbitals.
  @details
  Computes reduced matrix elements \f$ \redmatel{a}{h}{b} \f$ for all
  non-zero combinations of states from @p a_orbs and @p b_orbs. Optionally
  includes RPA corrections via @p dV.

  The transition frequency is set to \f$ \omega_{ab} = \varepsilon_a - \varepsilon_b \f$
  for off-diagonal elements, and 0 for diagonal elements.

  If @p each_freq is true, the operator and dV are updated at each transition
  frequency individually. Otherwise, they are updated once using @p omega.

  Selection rules are applied via the operator's `isZero()` method. For
  odd-parity operators, only elements with the even-parity state on the right
  are included (unless @p calculate_both is true). For even-parity operators,
  each off-diagonal pair is included once (upper triangle, unless
  @p calculate_both is true). When @p b_orbs != @p a_orbs,
  @p calculate_both is set to true automatically.

  @param b_orbs         Bra states (final states, index b).
  @param a_orbs         Ket states (initial states, index a).
  @param h              Pointer to the tensor operator.
  @param dV             Optional RPA/core-polarisation correction. If nullptr,
                        no correction is applied.
  @param omega          Frequency used to update frequency-dependent operators
                        and dV when @p each_freq is false.
  @param each_freq      If true, update the operator and dV at the natural
                        transition frequency for each element.
  @param diagonal       If true, include diagonal elements \f$ \redmatel{a}{h}{a} \f$.
  @param off_diagonal   If true, include off-diagonal elements.
  @param calculate_both If true, include both \f$ \redmatel{a}{h}{b} \f$ and
                        \f$ \redmatel{b}{h}{a} \f$ for off-diagonal pairs.
                        Forced true when @p b_orbs != @p a_orbs.

  @return Vector of MEdata structs, one per non-zero matrix element.
*/
std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &b_orbs,
                   const std::vector<DiracSpinor> &a_orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal = true,
                   bool off_diagonal = true, bool calculate_both = false);

/*!
  @brief Calculates reduced matrix elements for a single list of orbitals.
  @details
  Convenience overload; calls calcMatrixElements(orbs, orbs, ...) with both
  bra and ket taken from the same set @p orbs.
*/
inline std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal = true,
                   bool off_diagonal = true, bool calculate_both = false) {
  return calcMatrixElements(orbs, orbs, h, dV, omega, each_freq, diagonal,
                            off_diagonal, calculate_both);
}

/*!
  @brief Builds a lookup table of reduced matrix elements <a||h||b>.
  @details
  Fills and returns a `Coulomb::meTable<double>` with reduced matrix elements
  \f[ t_{ab} = \redmatel{a}{h}{b} + \delta V_{ab} + \delta_{\rm SR}^{ab} \f]
  for all non-zero pairs from @p a_orbs and @p b_orbs.

  Optionally includes RPA correction @p dV and structure radiation @p srn.
  The symmetry-conjugate \f$ \redmatel{b}{h}{a} \f$ is also stored via
  `symm_sign()`. Table is filled with OpenMP parallelisation.

  @param a_orbs  Bra states.
  @param b_orbs  Ket states.
  @param h       Pointer to the (const) tensor operator.
  @param dV      Optional RPA correction. If nullptr, not applied.
  @param srn     Optional structure radiation/normalisation correction.
                 If nullptr, not applied.
  @param omega   Only used for Structure Radiation. If set, use this fixed
                 frequency for all elements; otherwise uses |eb - ea|.

  @return meTable containing t_ab for all non-zero pairs (and conjugates).

  @note For frequency-dependent operators, @p omega should be set before
        calling this function.
*/
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h,
                                  const CorePolarisation *dV = nullptr,
                                  const MBPT::StructureRad *srn = nullptr,
                                  std::optional<double> omega = {});

/*!
  @brief Builds a matrix element table for a single set of orbitals.
  @details
  Convenience overload; calls me_table(a_orbs, a_orbs, ...) with both
  bra and ket taken from @p a_orbs.
*/
inline Coulomb::meTable<double> me_table(
  const std::vector<DiracSpinor> &a_orbs,
  const DiracOperator::TensorOperator *h, const CorePolarisation *dV = nullptr,
  const MBPT::StructureRad *srn = nullptr, std::optional<double> omega = {}) {
  return me_table(a_orbs, a_orbs, h, dV, srn, omega);
}

/*!
  @brief Factory function to construct a core-polarisation (RPA) object.
  @details
  Parses @p method and returns a `std::unique_ptr<CorePolarisation>` of the
  appropriate type. Returns nullptr if @p method is "none" or "false".

  Supported methods (case-insensitive):
  - `"TDHF"`: time-dependent Hartree-Fock (TDHF).
  - `"basis"`: TDHF solved in a basis set.
  - `"diagram"`: diagram RPA.
  - `"none"`, `"false"`, `""`: no RPA; returns nullptr.

  If the method string is not recognised, prints an error and defaults to none.

  @param method    String specifying the RPA method (see above).
  @param h         Pointer to the operator for which RPA is to be solved.
  @param vhf       Pointer to the Hartree-Fock object (provides core potential).
  @param print     If true, print a brief description of the chosen method.
  @param basis     Basis set for basis/diagram methods (ignored for TDHF).
  @param identity  Identifier string passed to DiagramRPA (e.g. for caching).

  @return Unique pointer to the constructed CorePolarisation object, or
          nullptr if RPA is disabled.

  @warning An unrecognised @p method string triggers an error message and
           falls through to no RPA rather than throwing.
*/
std::unique_ptr<CorePolarisation>
make_rpa(const std::string &method, const DiracOperator::TensorOperator *h,
         const HF::HartreeFock *vhf, bool print = false,
         const std::vector<DiracSpinor> &basis = {},
         const std::string &identity = "");

} // namespace ExternalField
