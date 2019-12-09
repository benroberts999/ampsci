#include "Dirac/DiracSpinor.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>

inline auto test_splines(int kappa, std::size_t n_spl, std::size_t k_spl,
                         double r0_spl, double rmax_spl, const Grid &rgrid,
                         const double alpha) {
  ChronoTimer sw("splines");

  BSplines bspl(n_spl, k_spl, rgrid, r0_spl, rmax_spl);

  bspl.derivitate();
  bspl.write_splines("Bspl.txt");
  bspl.write_splines("Bspl_deriv2.txt", true);

  // auto alpha = 1.0 / 137.035999084; // XXX

  std::vector<DiracSpinor> basis;

  auto imin = static_cast<std::size_t>(std::abs(kappa));
  auto imax = n_spl - 1;
  // auto half_n_orbs = (imax - imin);

  auto n_count = 1;
  for (auto i = imin; i < imax; i++) {

    basis.emplace_back(n_count++, kappa, rgrid);

    auto &phi = basis.back();
    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    phi.f = Bi;
    // NumCalc::scaleVec(phi.f, -1.0);

    auto gtmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    NumCalc::scaleVec(gtmp, double(kappa));
    phi.g = NumCalc::add_vectors(dBi, gtmp);
    NumCalc::scaleVec(phi.g, 0.5 * alpha);
    // phi *= -1;

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
    NumCalc::scaleVec(ftmp, double(-kappa));
    phi.f = NumCalc::add_vectors(dBi, ftmp);
    NumCalc::scaleVec(phi.f, 0.5 * alpha);
    // phi *= -1;

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    // std::cout << i + 1 << "/" << n_spl << ": " << phi.symbol() << "\n";
  }

  std::vector<DiracSpinor> d_basis_g;
  n_count = 1;
  for (auto i = imin; i < imax; i++) {

    d_basis_g.emplace_back(n_count++, kappa, rgrid);
    auto &phi = d_basis_g.back();

    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    auto d2Bi = bspl.get_spline_deriv2(i);

    auto dBior = NumCalc::mult_vectors(rgrid.inverse_r(), dBi);

    auto tmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    auto Bor2 = NumCalc::mult_vectors(rgrid.inverse_r(), tmp);

    NumCalc::scaleVec(dBior, double(kappa));
    NumCalc::scaleVec(Bor2, double(-kappa));

    phi.f = NumCalc::add_vectors(d2Bi, dBior, Bor2);
    NumCalc::scaleVec(phi.f, 0.5 * alpha);
    // NumCalc::scaleVec(phi.f, 0.5);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    // std::cout << i + 1 << "/" << n_spl << ": " << phi.symbol() << "\n";
  }
  for (auto i = imin; i < imax; i++) {

    d_basis_g.emplace_back(n_count++, kappa, rgrid);
    auto &phi = d_basis_g.back();

    auto dBi = bspl.get_spline_deriv(i);
    phi.f = dBi;
    // NumCalc::scaleVec(phi.f, 1.0 / alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    // std::cout << i + 1 << "/" << n_spl << ": " << phi.symbol() << "\n";
  }

  return std::make_pair(basis, d_basis_g);
}
