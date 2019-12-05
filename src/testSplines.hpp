#include "Dirac/DiracSpinor.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <algorithm>
#include <iostream>
#include <string>

inline std::vector<DiracSpinor> test_splines(const BSplines &bspl, int kappa) {

  bspl.write_splines("Bspl.txt");

  auto n_spl = bspl.get_n();
  std::vector<DiracSpinor> basis;
  const auto &rgrid = bspl.get_grid();

  for (auto i = 0ul; i < 2 * n_spl; i++) {

    basis.emplace_back(i + 1, kappa, rgrid);
    auto &phi = basis[i];
    if (i < n_spl) {
      phi.f = bspl.get_spline(i);
      auto [p0, pinf] = bspl.get_ends(i);
      phi.pinf = pinf;
      phi.p0 = p0;
    } else {
      phi.g = bspl.get_spline(i - n_spl);
      auto [p0, pinf] = bspl.get_ends(i - n_spl);
      phi.p0 = p0;
      phi.pinf = pinf;
    }
  }
  return basis;
}
