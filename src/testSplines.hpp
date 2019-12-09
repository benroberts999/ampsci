#include "Dirac/DiracSpinor.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"

#include <algorithm>
#include <iostream>
#include <string>

inline std::vector<DiracSpinor> test_splines(int kappa, std::size_t n_spl,
                                             std::size_t k_spl, double r0_spl,
                                             double rmax_spl,
                                             const Grid &rgrid) {
  ChronoTimer sw("splines");

  BSplines bspl(n_spl, k_spl, rgrid, r0_spl, rmax_spl);

  bspl.derivitate();
  bspl.write_splines("Bspl.txt");
  bspl.write_splines("Bspl_deriv2.txt", true);

  auto alpha = 1.0 / 137.035999084; // XXX

  std::vector<DiracSpinor> basis;

  auto imin = static_cast<std::size_t>(std::abs(kappa));
  auto imax = n_spl - 1;
  auto half_n_orbs = (imax - imin);

  auto n_count = 1;

  for (auto i = imin; i < imax; i++) {

    basis.emplace_back(n_count++, kappa, rgrid);

    auto &phi = basis.back();
    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    phi.f = Bi;

    auto gtmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    NumCalc::scaleVec(gtmp, double(-kappa));
    phi.g = NumCalc::add_vectors(dBi, gtmp);
    NumCalc::scaleVec(phi.g, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    // std::cout << i + 1 << "/" << n_spl << ": " << phi.symbol() << "\n";
  }

  for (auto i = imin; i < imax; i++) {

    basis.emplace_back(n_count++, kappa, rgrid);
    auto &phi = basis.back();

    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    phi.g = Bi;
    auto ftmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    NumCalc::scaleVec(ftmp, double(kappa));
    phi.f = NumCalc::add_vectors(dBi, ftmp);
    NumCalc::scaleVec(phi.f, 0.5 * alpha);
    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    // std::cout << i + 1 << "/" << n_spl << ": " << phi.symbol() << "\n";
  }

  // std::cin.get();

  // for (auto i = 0ul; i < 2 * half_n_orbs; i++) {
  //
  //   basis.emplace_back(int(i + 1), kappa, rgrid);
  //
  //   auto &phi = basis.back();
  //   if (i < half_n_orbs) {
  //     auto Bi = bspl.get_spline(i + imin);
  //     auto dBi = bspl.get_spline_deriv(i + imin);
  //     phi.f = Bi;
  //
  //     auto kbor = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
  //     NumCalc::scaleVec(kbor, double(kappa));
  //     auto tmpg = NumCalc::add_vectors(dBi, kbor);
  //     NumCalc::scaleVec(tmpg, 0.5 * alpha);
  //     phi.g = NumCalc::add_vectors(tmpg);
  //
  //     // auto gtmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
  //     // NumCalc::scaleVec(gtmp, double(kappa));
  //     // phi.g = NumCalc::add_vectors(dBi, gtmp);
  //     // NumCalc::scaleVec(phi.g, 0.5 * alpha);
  //
  //     auto [p0, pinf] = bspl.get_ends(i + imin);
  //     phi.pinf = pinf;
  //     phi.p0 = p0;
  //   } else {
  //     auto Bi = bspl.get_spline(i - half_n_orbs + imin);
  //     auto dBi = bspl.get_spline_deriv(i - half_n_orbs + imin);
  //     phi.g = Bi;
  //     auto ftmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
  //     NumCalc::scaleVec(ftmp, double(-kappa));
  //     phi.f = NumCalc::add_vectors(dBi, ftmp);
  //     NumCalc::scaleVec(phi.f, 0.5 * alpha);
  //     auto [p0, pinf] = bspl.get_ends(i - half_n_orbs + imin);
  //     phi.pinf = pinf;
  //     phi.p0 = p0;
  //   }
  // }

  return basis;
}
