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

struct MEdata {
  std::string a, b;
  double hab, dv1, dv;

  friend std::ostream &operator<<(std::ostream &os, const MEdata &m) {
    os << qip::fstring("<%4s,%4s>: %13.6e", m.a.c_str(), m.b.c_str(), m.hab);
    if (m.dv1 != 0.0) {
      os << qip::fstring("  %13.6e  %13.6e", m.hab + m.dv1, m.hab + m.dv);
    }
    return os;
  }
};

std::vector<MEdata>
calcMatrixElements(const std::vector<DiracSpinor> &orbs,
                   DiracOperator::TensorOperator *const h,
                   CorePolarisation *const dV = nullptr, double omega = 0.0,
                   bool each_freq = false, bool diagonal_only = false,
                   bool print_both = false, bool radial_int = false);

} // namespace ExternalField
