#include "DiracODE/AsymptoticSpinor.hpp"
#include "Angular/include.hpp"
#include "DiracODE/include.hpp"
#include "DiracOperator/include.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracContinuum.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: AsymptoticSpinor expansion", "[DiracODE][asymp][unit]") {

  std::cout << "AsymptoticSpinor expansion (large r)\n";

  fmt::print("{:<3s} {:>2s} {:1s} {:1s} {:>5s}  {:16s}  {:17s}   {}\n", "z",
             "k", "l", "n", "r", "f/g (asym)", "f/g (exact)", "eps");

  for (auto z : {0.1, 0.5, 1.0, 10.0, 100.0, 150.0}) {
    for (auto kappa : {-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6}) {
      const auto l = Angular::l_k(kappa);
      double w_eps = -1.0;
      int w_n{0};
      double w_asym{0.0}, w_exact{0.0}, w_r{0.0};

      if (z * PhysConst::alpha > 1.0 && std::abs(kappa) == 1)
        continue;

      const auto r0{10.0 / z};
      const auto rmax{150.0 / z};
      const auto num_grid_points{15ul};
      const auto b{(rmax + r0) / 4};
      const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                     GridType::loglinear, b);
      for (auto n : {1, 2, 3, 4, 5, 6, 7}) {

        if (l >= n)
          continue;

        const auto e = AtomData::diracen(z, n, kappa, PhysConst::alpha);
        const auto F1s = DiracSpinor::exactHlike(n, kappa, grid, z);
        REQUIRE(F1s.en() == Approx(e));
        DiracODE::AsymptoticSpinor x{kappa, z, e};
        for (std::size_t i = 0; i < num_grid_points; ++i) {
          const auto r = grid->r(i);
          const auto [f, g] = x.fg(r);
          const auto f0 = F1s.f(i);
          const auto g0 = F1s.g(i);
          if (g == 0.0 || f0 == 0.0 || g0 == 0.0) {
            continue;
          }

          const auto ratio_asym = f / g;
          const auto ratio_exact = f0 / g0;
          const auto eps = std::abs(ratio_asym / ratio_exact - 1.0);

          if (eps > w_eps) {
            w_eps = eps;
            w_n = n;
            w_asym = ratio_asym;
            w_exact = ratio_exact;
            w_r = r;
          }

          auto eps_targ = z > 99.0 ? 1.0e-10 :
                          z > 9.0  ? 1.0e-8 :
                          z > 0.9  ? 1.0e-7 :
                                     1.0e-5;

          if (l >= 5)
            eps_targ *= 100;
          if (l >= 5 && z < 1.0)
            eps_targ *= 100;

          REQUIRE(eps < eps_targ);
        }
      }

      // if (l >= 5)
      fmt::print("{:3g} {:2} {:1} {:1} {:5.1f} {:+.10e} [{:+.10e}]  {:.1e}\n",
                 z, kappa, l, w_n, w_r, w_asym, w_exact, w_eps);
    }
  }
}

//==============================================================================
// Relative residual of the radial Dirac-Coulomb equations, v(r) = -zeff/r,
// for the radial spinor {f(r), g(r)} returned by the callable fg_of_r.
// Derivatives are 5-point central finite differences of step h. Each
// equation's residual is scaled by its largest term; returns the larger of
// the two.
template <typename FGofR>
double asymptotic_dirac_residual(const FGofR &fg_of_r, int kappa, double zeff,
                                 double en, double alpha, double r, double h) {
  const auto [fp2, gp2] = fg_of_r(r + 2.0 * h);
  const auto [fp1, gp1] = fg_of_r(r + h);
  const auto [fm1, gm1] = fg_of_r(r - h);
  const auto [fm2, gm2] = fg_of_r(r - 2.0 * h);
  const auto df = (-fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2) / (12.0 * h);
  const auto dg = (-gp2 + 8.0 * gp1 - 8.0 * gm1 + gm2) / (12.0 * h);

  const auto [f, g] = fg_of_r(r);
  const auto v = -zeff / r;
  const auto c = 1.0 / alpha;
  // f' = -(kappa/r) f + (alpha (en - v) + 2c) g
  // g' = alpha (v - en) f + (kappa/r) g
  const std::array f_terms{df, (kappa / r) * f,
                           -(alpha * (en - v) + 2.0 * c) * g};
  const std::array g_terms{dg, -alpha * (v - en) * f, -(kappa / r) * g};

  double worst = 0.0;
  for (const auto &terms : {f_terms, g_terms}) {
    double sum = 0.0;
    double scale = 0.0;
    for (const auto term : terms) {
      sum += term;
      scale = std::max(scale, std::abs(term));
    }
    worst = std::max(worst, std::abs(sum) / scale);
  }
  return worst;
}

//==============================================================================
// Continuum (en > 0) tail: F^C and G^C must each satisfy the Dirac-Coulomb
// equation at large r, with Wronskian W[F^C, G^C] = -alpha/pi (energy
// normalisation). With FLINT, also projects the exact Dirac-Coulomb
// continuum function onto {F^C, G^C}: a^2 + b^2 = 1, with (a, b) constant
// in r.
TEST_CASE("DiracODE: AsymptoticSpinorContinuum tail",
          "[DiracODE][asymp][cntm][unit]") {

  const auto alpha = PhysConst::alpha;
  const auto W_expected = -alpha / M_PI;

  fmt::print("\nAsymptoticSpinorContinuum: large-r tail\n");
  fmt::print("{:>3s} {:>5s} {:>9s} {:>9s} {:>9s}\n", "z", "en", "W eps",
             "residual", "exact");

  for (const double z : {0.0, 1.0, 5.0}) {
    for (const double en : {0.05, 0.5, 5.0}) {
      double worst_w = 0.0;
      double worst_residual = 0.0;
      double worst_exact = 0.0;
      for (const int kappa : {-1, 1, -2, 2, -3, 3}) {
        const DiracODE::AsymptoticSpinorContinuum<15> asy{kappa, z, en, alpha};
        const auto p = asy.momentum();
        REQUIRE(p == Approx(DiracContinuum::pe(en, alpha)).epsilon(1.0e-14));
        REQUIRE(asy.beta() ==
                Approx(std::sqrt(en / (en + 2.0 / (alpha * alpha))))
                  .epsilon(1.0e-14));

        const auto FC = [&asy](double r) {
          const auto s = asy.fg(r);
          return std::pair{s.fC, s.gC};
        };
        const auto GC = [&asy](double r) {
          const auto s = asy.fg(r);
          return std::pair{s.fG, s.gG};
        };

        // The 1/r series converges for p r >> nu^2, nu = Z/p (Coulomb phase)
        const auto nu = z / p;
        const auto r_min = std::max(100.0, 20.0 * nu * nu) / p;
        const auto h = 1.0e-3 / p;

        // Exact solution projected onto {F^C, G^C}: F = a F^C + b G^C
        double a_prev = 0.0;
        double b_prev = 0.0;

        for (const double r_scale : {1.0, 3.0, 10.0}) {
          const auto r = r_scale * r_min;
          const auto [fC, gC, fG, gG] = asy.fg(r);

          const auto W = fC * gG - fG * gC;
          const auto w_eps = std::abs(W / W_expected - 1.0);
          const auto residual =
            std::max(asymptotic_dirac_residual(FC, kappa, z, en, alpha, r, h),
                     asymptotic_dirac_residual(GC, kappa, z, en, alpha, r, h));
          worst_w = std::max(worst_w, w_eps);
          worst_residual = std::max(worst_residual, residual);
          REQUIRE(w_eps < 1.0e-9);
          REQUIRE(residual < 1.0e-8);

          if (DiracContinuum::available) {
            const auto [f, g] = DiracContinuum::fg(r, en, kappa, z, alpha);
            const auto a = (f * gG - g * fG) / W;
            const auto b = (fC * g - gC * f) / W;
            const auto q_eps = std::abs(a * a + b * b - 1.0);
            worst_exact = std::max(worst_exact, q_eps);
            REQUIRE(q_eps < 1.0e-6);
            if (r_scale > 1.0) {
              worst_exact = std::max(
                {worst_exact, std::abs(a - a_prev), std::abs(b - b_prev)});
              REQUIRE(a == Approx(a_prev).margin(1.0e-6));
              REQUIRE(b == Approx(b_prev).margin(1.0e-6));
            }
            a_prev = a;
            b_prev = b;
          }
        }
      }
      if (DiracContinuum::available) {
        fmt::print("{:3g} {:5g} {:9.1e} {:9.1e} {:9.1e}\n", z, en, worst_w,
                   worst_residual, worst_exact);
      } else {
        fmt::print("{:3g} {:5g} {:9.1e} {:9.1e} {:>9s}\n", z, en, worst_w,
                   worst_residual, "no FLINT");
      }
    }
  }
}

//==============================================================================
//! The complex-energy instantiation, at real bound-state energy, must
//! reproduce the real (double) expansion exactly. Choosing the other branch
//! of lambda explicitly must give the growing solution of the same equation:
//! checks that the small-component amplitude follows the chosen branch.
TEST_CASE("DiracODE: AsymptoticSpinor complex energy reduces to bound",
          "[DiracODE][asymp][cntm][unit]") {
  using Complex = std::complex<double>;
  const auto alpha = PhysConst::alpha;

  fmt::print(
    "\nAsymptoticSpinor<complex> at real energy vs AsymptoticSpinor\n");
  fmt::print("{:>3s} {:>9s} {:>9s} {:>9s}\n", "z", "re eps", "im eps",
             "residual");

  for (const double z : {1.0, 10.0, 80.0}) {
    double worst_re = 0.0;
    double worst_im = 0.0;
    double worst_residual = 0.0;
    for (const int kappa : {-1, 1, -2, 2, -3}) {
      const auto l = Angular::l_k(kappa);
      for (const int n : {l + 1, l + 2}) {
        const auto en = AtomData::diracen(z, n, kappa, alpha);
        REQUIRE(en < 0.0);

        const DiracODE::AsymptoticSpinor<double> real_expansion{kappa, z, en,
                                                                alpha};
        const DiracODE::AsymptoticSpinor<Complex> complex_expansion{
          kappa, z, Complex{en}, alpha};

        // Same (principal) branch given explicitly: identical
        const auto lambda = std::sqrt(-en * (2.0 + en * alpha * alpha));
        const auto explicit_expansion =
          DiracODE::AsymptoticSpinor<double>::with_lambda(kappa, z, en, lambda,
                                                          alpha);
        // Other branch: the exponentially growing solution
        const auto growing_expansion =
          DiracODE::AsymptoticSpinor<double>::with_lambda(kappa, z, en, -lambda,
                                                          alpha);
        const auto growing_fg = [&growing_expansion](double r) {
          return growing_expansion.fg(r);
        };

        for (const double zr : {10.0, 50.0, 150.0}) {
          const auto r = zr / z;
          const auto [f, g] = real_expansion.fg(r);
          const auto [fc, gc] = complex_expansion.fg(r);
          const auto [fe, ge] = explicit_expansion.fg(r);

          const auto re_eps = std::max(std::abs(fc.real() / f - 1.0),
                                       std::abs(gc.real() / g - 1.0));
          const auto im_eps =
            std::max(std::abs(fc.imag() / f), std::abs(gc.imag() / g));
          worst_re = std::max(worst_re, re_eps);
          worst_im = std::max(worst_im, im_eps);
          REQUIRE(re_eps < 1.0e-13);
          REQUIRE(im_eps < 1.0e-14);
          REQUIRE(fe == f);
          REQUIRE(ge == g);
        }

        // Growing solution must satisfy the same Dirac equation. Its 1/r
        // coefficients go as (i + n)^2 rather than (i - n)^2, so the series
        // needs larger lambda r than the decaying one to converge
        for (const double lambda_r : {60.0, 120.0}) {
          const auto r = lambda_r / lambda;
          const auto h = 1.0e-3 / lambda;
          const auto residual =
            asymptotic_dirac_residual(growing_fg, kappa, z, en, alpha, r, h);
          worst_residual = std::max(worst_residual, residual);
          REQUIRE(residual < 1.0e-8);
          // ...and not be the decaying one
          const auto [fgrow, ggrow] = growing_expansion.fg(r);
          const auto [fdecay, gdecay] = real_expansion.fg(r);
          REQUIRE(std::abs(fgrow) > 10.0 * std::abs(fdecay));
        }
      }
    }
    fmt::print("{:3g} {:9.1e} {:9.1e} {:9.1e}\n", z, worst_re, worst_im,
               worst_residual);
  }
}
