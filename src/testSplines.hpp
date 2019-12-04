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
      // std::cout << i << ": ";
      // std::cout << phi.f[10] << " " << phi.f[11] << " " << phi.f[12] << "\n";
    } else {
      phi.g = bspl.get_spline(i - n_spl);
      // std::cout << i - n_spl << "\n";
      // std::cout << phi.g[10] << " " << phi.g[11] << " " << phi.g[12] << "\n";
    }
  }
  return basis;
}
