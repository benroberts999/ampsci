#pragma once
#include "CorePolarisation.hpp"
#include "qip/String.hpp"
#include <iostream>
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
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

  static std::string title() {
    return "    a    b   w_ab      t_ab           RPA_ab";
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
                   bool each_freq = false, bool diagonal_only = false,
                   bool print_both = false, bool radial_int = false);

inline std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal_only = false,
                   bool print_both = false, bool radial_int = false) {
  return calcMatrixElements(orbs, orbs, h, dV, omega, each_freq, diagonal_only,
                            print_both, radial_int);
}

} // namespace ExternalField
