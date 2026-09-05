#include "Angular/include.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
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
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
TEST_CASE("DiracODE: continuum relativistic", "[DiracODE][cntm][integration]") {
  if (!DiracContinuum::available) {
    std::cout << "FLINT not available: No analytic continuum states\n";
    SUCCEED("Compiled without FLINT: skipping numerical tests");
    return;
  }

  // Set true to print detailed (expected vs found) table:
  constexpr bool print_table = true;
  // Optionally fill the unresolved (zeroed) tail with local average:
  // nb: Grid always big enough!
  constexpr bool average_tail = false;

  // Bound-continuum radial integrals <en kappa| r^k |n kappa>:
  // DiracODE::solveContinuum solution vs exact (analytic) Dirac Coulomb
  // continuum function (DiracContinuum). The same exact bound state
  // (DiracHydrogen) and quadrature are used on both sides, so this isolates
  // the ODE continuum solution (including its energy normalisation).

  fmt::print("\nDiracODE continuum vs exact Dirac Coulomb: "
             "<en kappa|r^p|n kappa>\n");

  // For each energy, build a grid dense enough to store the continuum
  // state (RequiredContinuumGrid), with at least min_points. If more than
  // max_points required, clamp, and instead reduce b (more linear grid):
  const double r0 = 1.0e-5;
  const double rmax = 50.0;
  const double b0 = 0.1;
  const std::size_t min_points = 2000;
  // const std::size_t max_points = 20000;

  if (print_table) {
    fmt::print("{:>5} {:>3} {:>2} {:>3} {:>5} {:>5} {:>6} {:>2} {:>16} "
               "{:>16} {:>8}\n",
               "α/α0", "Z", "n", "kap", "en", "b", "npts", "p", "expected",
               "found", "eps");
  }

  const double alpha0 = PhysConst::alpha;
  for (const double en : {0.1, 1e2, 1.0e4, 1.0e5}) {

    // Grid for this energy:
    const Grid grid_0(r0, rmax, min_points, GridType::loglinear, b0);
    const auto req_N = DiracODE::RequiredContinuumGrid(en, grid_0).num_points;
    auto num_points = std::max(req_N, min_points);
    double b = b0;
    const auto grid = std::make_shared<const Grid>(r0, rmax, num_points,
                                                   GridType::loglinear, b);

    // Reference grid for the 'expected' side: must fully resolve the
    // oscillations, else the quadrature of the (exact) integrand is
    // aliasing junk. Same as working grid when no clamping occurred:
    const auto ref_N = std::max(req_N, min_points);
    const auto grid_ref =
      (ref_N == num_points && b == b0) ?
        grid :
        std::make_shared<const Grid>(r0, rmax, ref_N, GridType::loglinear, b0);

    for (const double alpha : {PhysConst::alpha, 0.001 * PhysConst::alpha}) {
      for (const double z : {1.0, 10.0, 100.0}) {

        const auto v0 = Nuclear::sphericalNuclearPotential(z, 0.0, grid->r());

        int count = 0;
        for (const int n : {1, 3}) {
          for (const int kappa : {-1, 1, 2, -3, 3, 6}) { // not all
            const int l = kappa > 0 ? kappa : -kappa - 1;
            if (l > n - 1) {
              continue;
            }

            // Get reasonable coverage, skip half combos (just for speed)
            if (++count % 2 == 0)
              continue;

            // Numerical (ODE) continuum solution, on the working grid:
            DiracSpinor Fe{0, kappa, grid};
            DiracODE::solveContinuum(Fe, en, v0, alpha, nullptr, nullptr,
                                     average_tail);
            const auto Fb = DiracSpinor::exactHlike(n, kappa, grid, z, alpha);

            // Exact (analytic) Dirac-Coulomb bound + continuum, on the
            // reference grid:
            const auto Fb_ref =
              DiracSpinor::exactHlike(n, kappa, grid_ref, z, alpha);
            const auto Fc_ref =
              DiracSpinor::exactHlike_cntm(en, kappa, grid_ref, z, alpha);

            for (const int p : {-1, 1}) {
              // <Fb| r^p |Fe(ODE)> (working grid) vs
              // <Fb| r^p |exact continuum> (reference grid):
              const auto rk = grid->rpow(p);
              const auto rk_ref = grid_ref->rpow(p);
              const auto found = Fb * (rk * Fe);
              const auto expected = Fb_ref * (rk_ref * Fc_ref);

              if (print_table) {
                fmt::print("{:>5.0e} {:>3.0f} {:>2} {:>3} {:>5.0e} {:>5.2f} "
                           "{:>6} {:>2} {:>16.8e} {:>16.8e} {:>8.0e}\n",
                           alpha / alpha0, z, n, kappa, en, b, num_points, p,
                           expected, found, (found - expected) / expected);
              }

              // Tolerances (calibrated to observed worst cases, ~3-10x
              // headroom). High Z is hardest: interior (Coulomb-boosted)
              // oscillations are under-resolved on the smallest grids
              // (worst for the en ~ 1e2 grid, which has fewest points).
              // Tiny (cancellation-dominated) MEs are covered by the
              // absolute margin: they sit at the comparison noise floor.
              const double eps_t = en <= 1.0   ? (z > 50.0 ? 1e-2 : 1e-6) :
                                   en <= 1.0e3 ? (z > 50.0 ? 3e-1 : 3e-2) :
                                                 1e-2;
              REQUIRE(found == Approx(expected).epsilon(eps_t).margin(1e-9));
            }
          }
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("DiracODE: continuum relativistic - unit", "[DiracODE][cntm][unit]") {
  if (!DiracContinuum::available) {
    std::cout << "FLINT not available: No analytic continuum states\n";
    SUCCEED("Compiled without FLINT: skipping numerical tests");
    return;
  }

  // Fast version of "DiracODE: continuum relativistic": low energies,
  // single Zeff, all kappa up to l = 6 (lowest n for each kappa).
  // Bound-continuum radial integrals <en kappa| r^p |n kappa>, ODE solution
  // vs exact (analytic) Dirac Coulomb, same grid and quadrature both sides.

  constexpr bool print_table = false;

  const double r0 = 1.0e-5;
  const double rmax = 50.0;
  const double b0 = 0.1;
  const std::size_t min_points = 2000;
  const double zeff = 5.0;
  const double alpha = PhysConst::alpha;

  fmt::print("\nDiracODE continuum vs exact (unit): Zeff = {:.0f}\n", zeff);
  if (print_table) {
    fmt::print("{:>5} {:>3} {:>2} {:>6} {:>2} {:>16} {:>16} {:>8}\n", "en",
               "kap", "n", "npts", "p", "expected", "found", "eps");
  }

  for (const double en : {0.1, 1.0, 10.0}) {

    // Grid for this energy (min_points always sufficient at these en):
    const Grid grid_0(r0, rmax, min_points, GridType::loglinear, b0);
    const auto req_N = DiracODE::RequiredContinuumGrid(en, grid_0).num_points;
    const auto num_points = std::max(req_N, min_points);
    const auto grid = std::make_shared<const Grid>(r0, rmax, num_points,
                                                   GridType::loglinear, b0);

    const auto v0 = Nuclear::sphericalNuclearPotential(zeff, 0.0, grid->r());

    for (std::size_t k_i = 0;; ++k_i) {
      const auto kappa = Angular::kindex_to_kappa(k_i);
      const auto l = Angular::l_k(kappa);
      if (l > 6) {
        break;
      }
      const auto n = l + 1;

      DiracSpinor Fe{0, kappa, grid};
      DiracODE::solveContinuum(Fe, en, v0, alpha);

      const auto Fb = DiracSpinor::exactHlike(n, kappa, grid, zeff, alpha);
      const auto Fc =
        DiracSpinor::exactHlike_cntm(en, kappa, grid, zeff, alpha);

      for (const int p : {-1, 1}) {
        const auto rk = grid->rpow(p);
        const auto found = Fb * (rk * Fe);
        const auto expected = Fb * (rk * Fc);

        if (print_table) {
          fmt::print("{:>5.0e} {:>3} {:>2} {:>6} {:>2} {:>16.8e} {:>16.8e} "
                     "{:>8.0e}\n",
                     en, kappa, n, num_points, p, expected, found,
                     (found - expected) / expected);
        }

        REQUIRE(found == Approx(expected).epsilon(1.0e-6));
      }
    }
  }
}

//==============================================================================
TEST_CASE("DiracODE: continuum averageTail", "[DiracODE][cntm][unit]") {
  if (!DiracContinuum::available) {
    std::cout << "FLINT not available: No analytic continuum states\n";
    SUCCEED("Compiled without FLINT: skipping numerical tests");
    return;
  }

  // averageTail: on a grid too coarse to store the continuum oscillations
  // at large r (below the ~10 ppw zeroing threshold), the locally-averaged
  // tail must do as well or better than the zeroed tail, for matrix
  // elements against (smooth) bound states. Use a grid ~3x too coarse
  // (truncation radius well inside the bound-state extent), and compare
  // both against the exact result from a fully-resolved reference grid.

  constexpr bool print_table = false;

  const double r0 = 1.0e-5;
  const double rmax = 40.0;
  const double b0 = 0.1;
  const double zeff = 5.0;
  const double en = 1.0e3;
  const double alpha = PhysConst::alpha;

  // Fully-resolved reference grid, and ~3x too coarse working grid:
  const Grid grid_0(r0, rmax, 2000, GridType::loglinear, b0);
  const auto req_N = DiracODE::RequiredContinuumGrid(en, grid_0).num_points;
  const auto grid_ref =
    std::make_shared<const Grid>(r0, rmax, req_N, GridType::loglinear, b0);
  const auto grid_c =
    std::make_shared<const Grid>(r0, rmax, req_N / 3, GridType::loglinear, b0);

  const auto v0_ref =
    Nuclear::sphericalNuclearPotential(zeff, 0.0, grid_ref->r());
  const auto v0_c = Nuclear::sphericalNuclearPotential(zeff, 0.0, grid_c->r());

  fmt::print("\naverageTail: coarse grid ({} points; require {}), "
             "Zeff = {:.0f}, en = {:.0f}\n",
             grid_c->num_points(), req_N, zeff, en);
  if (print_table) {
    fmt::print("{:>3} {:>2} {:>2} {:>16} {:>16} {:>16} {:>8} {:>8}\n", "kap",
               "n", "p", "expected", "zeroed", "averaged", "eps_zero",
               "eps_avg");
  }

  double worst_zero = 0.0;
  double worst_avg = 0.0;
  for (std::size_t k_i = 0;; ++k_i) {
    const auto kappa = Angular::kindex_to_kappa(k_i);
    const auto l = Angular::l_k(kappa);
    if (l > 2) {
      break;
    }
    const auto n = l + 1;

    // Coarse-grid solutions, with zeroed and with averaged tail:
    DiracSpinor Fe_zero{0, kappa, grid_c};
    DiracODE::solveContinuum(Fe_zero, en, v0_c, alpha);
    DiracSpinor Fe_avg{0, kappa, grid_c};
    DiracODE::solveContinuum(Fe_avg, en, v0_c, alpha, nullptr, nullptr, true);

    const auto Fb = DiracSpinor::exactHlike(n, kappa, grid_c, zeff, alpha);

    // Exact, on the reference grid:
    const auto Fb_ref =
      DiracSpinor::exactHlike(n, kappa, grid_ref, zeff, alpha);
    const auto Fc_ref =
      DiracSpinor::exactHlike_cntm(en, kappa, grid_ref, zeff, alpha);

    for (const int p : {-1, 1}) {
      const auto rk = grid_c->rpow(p);
      const auto rk_ref = grid_ref->rpow(p);
      const auto expected = Fb_ref * (rk_ref * Fc_ref);
      const auto found_zero = Fb * (rk * Fe_zero);
      const auto found_avg = Fb * (rk * Fe_avg);

      const auto err_zero = std::abs(found_zero - expected);
      const auto err_avg = std::abs(found_avg - expected);

      if (print_table) {
        fmt::print("{:>3} {:>2} {:>2} {:>16.8e} {:>16.8e} {:>16.8e} "
                   "{:>8.0e} {:>8.0e}\n",
                   kappa, n, p, expected, found_zero, found_avg,
                   err_zero / std::abs(expected), err_avg / std::abs(expected));
      }

      // Averaging must do as well or better, for every matrix element
      // (small slack for near-tie flukes):
      REQUIRE(err_avg <= err_zero + 1.0e-9);

      worst_zero = std::max(worst_zero, err_zero);
      worst_avg = std::max(worst_avg, err_avg);
    }
  }

  // ... and substantially better overall:
  fmt::print("worst |err|: zeroed = {:.1e}, averaged = {:.1e}\n", worst_zero,
             worst_avg);
  REQUIRE(worst_avg < 0.1 * worst_zero);
}

//==============================================================================
// The irregular partner F_irr of the regular continuum solution F_reg, and
// the standing-wave inhomogeneous solve built on the pair. Checks:
//  - the Wronskian f_reg g_irr - f_irr g_reg is constant over the grid
//    (F_irr solves the same homogeneous equation), and equals alpha/pi for
//    the energy-normalised F_reg (the seed is the exact quarter-wave partner
//    of the asymptotic Coulomb pair; barely-open energies use the outward
//    extension);
//  - the forward solve of (h - en) phi = S gives phi -> K F_irr beyond the
//    source, with K = -pi <F_reg|S> (the Green's-function identity).
TEST_CASE("DiracODE: continuum irregular partner and forward solve",
          "[DiracODE][cntm][TDHFcntm][unit]") {

  const double alpha = PhysConst::alpha;
  const auto grid =
    std::make_shared<const Grid>(1.0e-6, 60.0, 3000, GridType::loglinear, 1.0);
  const auto num_points = grid->num_points();

  // Coulomb potential of a singly-charged ion (the tail an ejected electron
  // sees)
  std::vector<double> v(num_points);
  for (std::size_t i = 0; i < num_points; ++i) {
    v[i] = -1.0 / grid->r(i);
  }

  fmt::print("\nContinuum pair (F_reg, F_irr) and forward solve, Z_ion = 1:\n");
  fmt::print("{:>5s} {:>6s} {:>12s} {:>12s} {:>10s} {:>12s} {:>12s} {:>10s}\n",
             "kappa", "en", "W (found)", "alpha/pi", "W var.", "K (found)",
             "-pi<Freg|S>", "phi/K.Firr");

  for (const int kappa : {-1, 1, -2, 2}) {
    // A short-ranged source: bound-like radial functions
    DiracSpinor S(0, kappa, grid);
    for (std::size_t i = 0; i < num_points; ++i) {
      const auto r = grid->r(i);
      S.f(i) = r * std::exp(-r);
      S.g(i) = 0.5 * alpha * r * std::exp(-r);
    }
    S.min_pt() = 0;
    S.max_pt() = num_points;

    for (const double en : {0.05, 0.5, 3.0}) {
      DiracSpinor Freg(0, kappa, grid);
      DiracSpinor Firr(0, kappa, grid);
      DiracODE::solveContinuum(Freg, en, v, alpha);
      REQUIRE(Freg.norm2() != 0.0);
      DiracODE::solveContinuumIrregular(Firr, Freg, en, v, alpha);
      const auto pinf = Freg.max_pt();
      REQUIRE(Firr.max_pt() == pinf);

      // Wronskian: constant over the grid (away from the origin, where
      // F_irr grows like r^{-l-1}), and equal to alpha/pi
      const auto W_expected = alpha / M_PI;
      double W_min = 1.0e300;
      double W_max = -1.0e300;
      for (std::size_t i = pinf / 20; i < pinf; ++i) {
        const auto W = Freg.f(i) * Firr.g(i) - Firr.f(i) * Freg.g(i);
        W_min = std::min(W_min, W);
        W_max = std::max(W_max, W);
      }
      const auto W_mid = 0.5 * (W_min + W_max);
      const auto W_variation = (W_max - W_min) / std::abs(W_mid);

      // Forward solve: K from the boundary-condition extraction vs the
      // overlap identity, and phi = K F_irr in the outer region
      DiracSpinor phi(0, kappa, grid);
      const auto K =
        DiracODE::solveContinuumForward(phi, Freg, Firr, en, v, alpha, S);
      const auto K_overlap = -M_PI * (Freg * S);
      double phi_ratio = 0.0;
      {
        const auto i = std::size_t(0.85 * double(pinf));
        phi_ratio = (phi.f(i) * Firr.f(i) + phi.g(i) * Firr.g(i)) /
                    (K * (Firr.f(i) * Firr.f(i) + Firr.g(i) * Firr.g(i)));
      }

      fmt::print("{:>5d} {:6.2f} {:12.6e} {:12.6e} {:10.1e} {:12.5e} "
                 "{:12.5e} {:10.6f}\n",
                 kappa, en, W_mid, W_expected, W_variation, K, K_overlap,
                 phi_ratio);

      REQUIRE(W_variation < 1.0e-5);
      REQUIRE(W_mid == Approx(W_expected).epsilon(0.05));
      REQUIRE(K == Approx(K_overlap).epsilon(1.0e-3));
      REQUIRE(phi_ratio == Approx(1.0).epsilon(1.0e-3));
    }
  }
}
