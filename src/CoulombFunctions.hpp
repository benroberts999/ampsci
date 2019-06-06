#pragma once
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

namespace Coulomb {
const int xyz = 3;
void calculate_v_abk(const DiracSpinor &phi_a, const DiracSpinor &phi_b,
                     const int k, std::vector<double> &vabk);
} // namespace Coulomb
