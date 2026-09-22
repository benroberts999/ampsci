#include "DiracODE/FreeDirac.hpp"
#include "Angular/Wigner369j.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/DiracContinuum.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

//==============================================================================
TEST_CASE("DiracODE: freeDirac", "[DiracODE][cntm][unit]") {
  std::cout << "Free (V = 0) Dirac spherical waves\n";

  const auto grid = std::make_shared<const Grid>(1.0e-6, 40.0, 8000ul,
                                                 GridType::loglinear, 5.0);
  const auto alpha = PhysConst::alpha;

  // Free Dirac equation on the grid (V = 0):
  //   g = alpha (f' + kappa f/r) / (2 + alpha^2 en)
  //   f = -(g' - kappa g/r) / (alpha en)
  // Checks the small component, its sign, and the l' = l(-kappa) choice
  for (const auto en : {0.5, 5.0, 60.0}) {
    const auto waves = DiracODE::freeDirac(en, 0, 6, grid, alpha);
    // l = 0 has one kappa, l = 1..6 two each
    REQUIRE(waves.size() == 13);
    double worst = 0.0;
    for (const auto &Fk : waves) {
      REQUIRE(Fk.en() == en);
      REQUIRE(Fk.max_pt() > 100);
      const auto df = NumCalc::derivative(Fk.f(), grid->drdu(), grid->du());
      const auto dg = NumCalc::derivative(Fk.g(), grid->drdu(), grid->du());
      double err_g = 0.0, err_f = 0.0, max_f = 0.0, max_g = 0.0;
      // Interior of the resolved region (one-sided differences at the ends)
      for (std::size_t i = 10; i + 10 < Fk.max_pt(); ++i) {
        const auto r = grid->r(i);
        const auto kappa = Fk.kappa();
        const auto g_from_f =
          alpha * (df[i] + kappa * Fk.f(i) / r) / (2.0 + alpha * alpha * en);
        const auto f_from_g = -(dg[i] - kappa * Fk.g(i) / r) / (alpha * en);
        err_g = std::max(err_g, std::abs(g_from_f - Fk.g(i)));
        err_f = std::max(err_f, std::abs(f_from_g - Fk.f(i)));
        max_f = std::max(max_f, std::abs(Fk.f(i)));
        max_g = std::max(max_g, std::abs(Fk.g(i)));
      }
      worst = std::max({worst, err_g / max_g, err_f / max_f});
      // Zero beyond the resolved region
      for (auto i = Fk.max_pt(); i < grid->num_points(); ++i) {
        REQUIRE(Fk.f(i) == 0.0);
        REQUIRE(Fk.g(i) == 0.0);
      }
    }
    fmt::print("  en = {:5.1f}: Dirac equation residual (rel.) = {:.1e}\n", en,
               worst);
    REQUIRE(worst < 1.0e-6);
  }

  // Non-relativistic limit: f -> energy-normalised free Coulomb function with
  // zero charge, P_el = sqrt(2/(pi p)) p r j_l(pr); checks the normalisation
  {
    const auto alpha_nr = 1.0e-4 * alpha;
    double worst = 0.0;
    for (const auto en : {0.3, 8.0}) {
      for (const auto kappa : {-1, 1, -2, 3, -5}) {
        const auto Fk = DiracODE::freeDirac(en, kappa, grid, alpha_nr);
        REQUIRE(Fk.kappa() == kappa);
        double err = 0.0, max_f = 0.0;
        for (std::size_t i = 0; i < Fk.max_pt(); ++i) {
          const auto expected =
            DiracContinuum::P_el(grid->r(i), en, Fk.l(), 0.0);
          err = std::max(err, std::abs(Fk.f(i) - expected));
          max_f = std::max(max_f, std::abs(expected));
        }
        worst = std::max(worst, err / max_f);
      }
    }
    fmt::print("  NR limit vs P_el(Z = 0): rel. error = {:.1e}\n", worst);
    REQUIRE(worst < 1.0e-8);
  }

  // Relativistic: exact Dirac-Coulomb continuum function as Z -> 0 (needs
  // FLINT). Checks the relativistic normalisation and the small component
  if (DiracContinuum::available) {
    const auto zeff = 1.0e-7;
    double worst = 0.0;
    for (const auto en : {1.0, 20.0}) {
      for (const auto kappa : {-1, 1, -3, 3}) {
        const auto Fk = DiracODE::freeDirac(en, kappa, grid, alpha);
        double err_f = 0.0, err_g = 0.0, max_f = 0.0, max_g = 0.0;
        // every 20th point is plenty (FLINT is slow)
        for (std::size_t i = 0; i < Fk.max_pt(); i += 20) {
          const auto [f, g] =
            DiracContinuum::fg(grid->r(i), en, kappa, zeff, alpha);
          err_f = std::max(err_f, std::abs(Fk.f(i) - f));
          err_g = std::max(err_g, std::abs(Fk.g(i) - g));
          max_f = std::max(max_f, std::abs(f));
          max_g = std::max(max_g, std::abs(g));
        }
        worst = std::max({worst, err_f / max_f, err_g / max_g});
      }
    }
    fmt::print("  vs Dirac-Coulomb (Z = 1e-7): rel. error = {:.1e}\n", worst);
    REQUIRE(worst < 1.0e-5);
  } else {
    std::cout << "  (no FLINT: skipping Dirac-Coulomb comparison)\n";
  }

  // Single-kappa and all-l overloads must agree
  {
    const auto waves = DiracODE::freeDirac(3.0, 0, 5, grid, alpha);
    double worst = 0.0;
    for (const auto &Fk : waves) {
      const auto single = DiracODE::freeDirac(3.0, Fk.kappa(), grid, alpha);
      REQUIRE(single.max_pt() == Fk.max_pt());
      double err = 0.0, max_f = 0.0;
      for (std::size_t i = 0; i < Fk.max_pt(); ++i) {
        err = std::max({err, std::abs(single.f(i) - Fk.f(i)),
                        std::abs(single.g(i) - Fk.g(i))});
        max_f = std::max(max_f, std::abs(Fk.f(i)));
      }
      worst = std::max(worst, err / max_f);
    }
    fmt::print("  single kappa vs all l: rel. difference = {:.1e}\n", worst);
    REQUIRE(worst < 1.0e-12);
  }

  // Truncation where the grid does not resolve the oscillations: as
  // solveContinuum, max_pt marks the resolved region
  {
    const auto coarse = std::make_shared<const Grid>(1.0e-6, 40.0, 400ul,
                                                     GridType::loglinear, 5.0);
    const auto low = DiracODE::freeDirac(0.1, -1, coarse, alpha);
    REQUIRE(low.max_pt() == coarse->num_points());
    const auto high = DiracODE::freeDirac(200.0, -1, coarse, alpha);
    REQUIRE(high.max_pt() < coarse->num_points());
    REQUIRE(high.max_pt() > 0);
    // resolved: at least 10 points per wavelength at the last stored point
    const auto k = std::sqrt(200.0 * (2.0 + alpha * alpha * 200.0));
    const auto i = high.max_pt() - 1;
    REQUIRE(coarse->drdu(i) * coarse->du() <= 2.0 * M_PI / (10.0 * k));
    REQUIRE(coarse->drdu(i + 1) * coarse->du() > 2.0 * M_PI / (10.0 * k));
  }
}
