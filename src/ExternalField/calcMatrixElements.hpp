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

//! Small struct to store calculated matrix elements
/*!
@details 
a and b: shortSymbols for finail/initial states <a||h||b>

hab: Lowest-order matrix element

dv: RPA correction

w_wb: frequency: Ea-Eb
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

std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &b_orbs,
                   const std::vector<DiracSpinor> &a_orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal = true,
                   bool off_diagonal = true, bool calculate_both = false);

inline std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal = true,
                   bool off_diagonal = true, bool calculate_both = false) {
  return calcMatrixElements(orbs, orbs, h, dV, omega, each_freq, diagonal,
                            off_diagonal, calculate_both);
}

//! Fills me_table with MEs, <a||h||b> and <b||h||a>.
//! Required to set omega for freq. dependent operators initially
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h,
                                  const CorePolarisation *dV = nullptr,
                                  const MBPT::StructureRad *srn = nullptr,
                                  std::optional<double> omega = {});

//! Fills me_table with MEs.
//! Required to set omega for freq. dependent operators initially
inline Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const DiracOperator::TensorOperator *h,
         const CorePolarisation *dV = nullptr,
         const MBPT::StructureRad *srn = nullptr,
         std::optional<double> omega = {}) {
  return me_table(a_orbs, a_orbs, h, dV, srn, omega);
}

std::unique_ptr<CorePolarisation>
make_rpa(const std::string &method, const DiracOperator::TensorOperator *h,
         const HF::HartreeFock *vhf, bool print = false,
         const std::vector<DiracSpinor> &basis = {},
         const std::string &identity = "");

} // namespace ExternalField
