#include "Kionisation/Kion_ridge.hpp"
#include "DiracODE/FreeDirac.hpp"
#include "Kionisation/Kion_functions.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

namespace {
// Minimal 4x4 complex matrices, for the reference Dirac traces
using Cx = std::complex<double>;
using M2 = std::array<Cx, 4>;
using M4 = std::array<Cx, 16>;

M2 pauli(int i) {
  const Cx I{0.0, 1.0};
  if (i == 0)
    return {0.0, 1.0, 1.0, 0.0};
  if (i == 1)
    return {0.0, -I, I, 0.0};
  return {1.0, 0.0, 0.0, -1.0};
}
M2 sigma_dot(const std::array<double, 3> &v) {
  M2 out{};
  for (int i = 0; i < 3; ++i) {
    const auto s = pauli(i);
    for (std::size_t j = 0; j < 4; ++j) {
      out[j] += v[std::size_t(i)] * s[j];
    }
  }
  return out;
}
M2 scale(Cx a, const M2 &m) { return {a * m[0], a * m[1], a * m[2], a * m[3]}; }
const M2 I2{1.0, 0.0, 0.0, 1.0};
const M2 Z2{0.0, 0.0, 0.0, 0.0};
M4 block(const M2 &a, const M2 &b, const M2 &c, const M2 &d) {
  M4 m{};
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      m[i * 4 + j] = a[i * 2 + j];
      m[i * 4 + j + 2] = b[i * 2 + j];
      m[(i + 2) * 4 + j] = c[i * 2 + j];
      m[(i + 2) * 4 + j + 2] = d[i * 2 + j];
    }
  }
  return m;
}
M4 mul(const M4 &a, const M4 &b) {
  M4 m{};
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      for (std::size_t k = 0; k < 4; ++k) {
        m[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
      }
    }
  }
  return m;
}
M4 dagger(const M4 &a) {
  M4 m{};
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      m[i * 4 + j] = std::conj(a[j * 4 + i]);
    }
  }
  return m;
}
Cx trace(const M4 &a) { return a[0] + a[5] + a[10] + a[15]; }
// Tr[Lam Ga rho Gb^dag]
Cx T(const M4 &Lam, const M4 &Ga, const M4 &rho, const M4 &Gb) {
  return trace(mul(mul(mul(Lam, Ga), rho), dagger(Gb)));
}

// The 11 traces by brute force, for a general direction of p (azimuth phi)
std::array<double, 11> numeric_traces(double p, double F, double G, int kappa,
                                      double pf, double q, double ef,
                                      double alpha, double phi) {
  const double c = 1.0 / alpha;
  const double mc2 = c * c;
  const double Ef = mc2 + ef;
  const double cos_t = (pf * pf - p * p - q * q) / (2.0 * p * q);
  const double sin_t = std::sqrt(1.0 - cos_t * cos_t);
  const std::array<double, 3> n{sin_t * std::cos(phi), sin_t * std::sin(phi),
                                cos_t};
  const std::array<double, 3> u{p * n[0], p * n[1], p * n[2] + q};
  const double s = kappa > 0 ? 1.0 : -1.0;

  const auto Lam = block(scale(Ef + mc2, I2), scale(c, sigma_dot(u)),
                         scale(c, sigma_dot(u)), scale(Ef - mc2, I2));
  const auto rho = block(scale(F * F, I2), scale(s * F * G, sigma_dot(n)),
                         scale(s * F * G, sigma_dot(n)), scale(G * G, I2));
  const auto one = block(I2, Z2, Z2, I2);
  const auto beta = block(I2, Z2, Z2, scale(-1.0, I2));
  const auto g5 = block(Z2, I2, I2, Z2);
  std::array<M4, 3> al, Sig;
  for (int i = 0; i < 3; ++i) {
    al[std::size_t(i)] = block(Z2, pauli(i), pauli(i), Z2);
    Sig[std::size_t(i)] = block(pauli(i), Z2, Z2, pauli(i));
  }
  const auto P = mul(
    mul(block(scale(Cx{0.0, 1.0}, I2), Z2, Z2, scale(Cx{0.0, 1.0}, I2)), beta),
    g5);

  std::array<double, 11> out;
  out[0] = std::real(T(Lam, one, rho, one));
  out[1] = std::real(T(Lam, al[2], rho, al[2]));
  out[2] = std::real(T(Lam, al[0], rho, al[0]) + T(Lam, al[1], rho, al[1]));
  out[3] = std::real(T(Lam, one, rho, al[2]));
  out[4] = std::real(T(Lam, g5, rho, g5));
  out[5] = std::real(T(Lam, Sig[2], rho, Sig[2]));
  out[6] = std::real(T(Lam, Sig[0], rho, Sig[0]) + T(Lam, Sig[1], rho, Sig[1]));
  out[7] = std::real(T(Lam, g5, rho, Sig[2]));
  out[8] = std::imag(T(Lam, al[0], rho, Sig[1]));
  out[9] = std::real(T(Lam, beta, rho, beta));
  out[10] = std::real(T(Lam, P, rho, P));
  return out;
}
} // namespace

//==============================================================================
TEST_CASE("Kion: plane-wave traces and momentum orbitals",
          "[Kion][ridge][unit]") {
  std::cout << "Plane-wave traces and momentum-space orbitals\n";
  const auto alpha = PhysConst::alpha;

  // Closed-form traces vs brute-force 4x4 traces, at random kinematics and
  // random azimuth of p (the returned combinations must not depend on it)
  {
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    double worst = 0.0;
    for (int trial = 0; trial < 200; ++trial) {
      const double ef = 1.0 + 3000.0 * uni(gen);
      const double pf = std::sqrt(ef * (2.0 + alpha * alpha * ef));
      const double q = pf * (0.2 + 3.0 * uni(gen));
      const double p = std::abs(pf - q) + 2.0 * std::min(pf, q) * uni(gen);
      const double F = uni(gen) - 0.5;
      const double G = 0.3 * (uni(gen) - 0.5);
      const int kappa = uni(gen) < 0.5 ? -2 : 3;
      const double phi = 2.0 * M_PI * uni(gen);
      const auto closed =
        Kion::planewave_traces(p, F, G, kappa, pf, q, ef, alpha);
      const auto numeric =
        numeric_traces(p, F, G, kappa, pf, q, ef, alpha, phi);
      const auto scale_ = std::abs(numeric[0]);
      for (std::size_t i = 0; i < 11; ++i) {
        worst = std::max(worst, std::abs(closed[i] - numeric[i]) / scale_);
      }
    }
    fmt::print("  closed-form vs numeric traces: {:.1e}\n", worst);
    REQUIRE(worst < 1.0e-12);
  }

  // Electron at rest (p -> 0, F = 1, G = 0), on shell (q = k): the
  // free-electron ratios V_L = (E/qc)^2 V_T, X = -(E/qc) V_T, A_L = V_T,
  // A_T = V_L, A_E + A_M = 2 V_T, V_E + V_M = 2 (E_f - m)/(E_f + m) V_T
  {
    const double c = 1.0 / alpha;
    double worst = 0.0;
    for (const auto ef : {10.0, 1000.0, 20000.0}) {
      const double pf = std::sqrt(ef * (2.0 + alpha * alpha * ef));
      const double Ef = c * c + ef;
      const auto T =
        Kion::planewave_traces(1.0e-9, 1.0, 0.0, -1, pf, pf, ef, alpha);
      const auto w_q = ef / (pf * c);
      const auto ratio = [&](std::size_t i) { return T[i] / T[0]; };
      worst = std::max(worst, std::abs(ratio(1) - w_q * w_q));
      worst = std::max(worst, std::abs(ratio(2) - 2.0 * ef / (Ef + c * c)));
      worst = std::max(worst, std::abs(ratio(3) - w_q));
      worst = std::max(worst, std::abs(ratio(4) - w_q * w_q));
      worst = std::max(worst, std::abs(ratio(5) - 1.0));
      worst = std::max(worst, std::abs(ratio(6) - 2.0));
      worst = std::max(worst, std::abs(ratio(7) - w_q));
      worst = std::max(worst, std::abs(ratio(9) - 1.0));
      worst = std::max(worst, std::abs(ratio(10) - w_q * w_q));
    }
    fmt::print("  electron-at-rest ratios: {:.1e}\n", worst);
    REQUIRE(worst < 1.0e-8);
  }

  // Momentum-space orbitals: normalisation, and hydrogen 1s analytic
  {
    const auto grid = std::make_shared<const Grid>(1.0e-6, 60.0, 5000ul,
                                                   GridType::loglinear, 4.0);
    double worst = 0.0;
    for (const auto &[n, kappa] :
         {std::pair{1, -1}, {2, -1}, {2, 1}, {2, -2}, {3, 2}, {4, -3}}) {
      auto Fa = DiracSpinor::exactHlike(n, kappa, grid, 2.0, alpha);
      Fa.occ_frac() = 1.0;
      const auto orb = Kion::momentum_orbital(Fa, Fa.en());
      REQUIRE(orb.p.size() > 200);
      // int (F^2 + G^2) p^2 dp = (2 pi)^3, trapezoid on the log grid
      double norm = 0.0;
      for (std::size_t i = 1; i < orb.p.size(); ++i) {
        const auto d0 =
          orb.F[i - 1] * orb.F[i - 1] + orb.G[i - 1] * orb.G[i - 1];
        const auto d1 = orb.F[i] * orb.F[i] + orb.G[i] * orb.G[i];
        norm += 0.5 * (orb.p[i] - orb.p[i - 1]) *
                (d0 * orb.p[i - 1] * orb.p[i - 1] + d1 * orb.p[i] * orb.p[i]);
      }
      norm /= std::pow(2.0 * M_PI, 3);
      worst = std::max(worst, std::abs(norm - 1.0));
    }
    fmt::print("  momentum density normalisation: 1 {:+.1e}\n", worst);
    REQUIRE(worst < 1.0e-4);

    // Hydrogen 1s: F(p) = 16 pi/(1+p^2)^2 (non-relativistic; O(alpha^2)
    // here)
    auto F1s = DiracSpinor::exactHlike(1, -1, grid, 1.0, alpha);
    F1s.occ_frac() = 1.0;
    const auto orb = Kion::momentum_orbital(F1s, F1s.en());
    double worst_1s = 0.0;
    for (std::size_t i = 0; i < orb.p.size(); i += 50) {
      const auto p = orb.p[i];
      if (p > 5.0)
        break;
      const auto expected = 16.0 * M_PI / std::pow(1.0 + p * p, 2);
      worst_1s = std::max(worst_1s, std::abs(orb.F[i] - expected) / expected);
    }
    fmt::print("  hydrogen 1s F(p) vs 16 pi/(1+p^2)^2: {:.1e}\n", worst_1s);
    REQUIRE(worst_1s < 2.0e-4);
  }
}

//==============================================================================
TEST_CASE("Kion: all-K plane-wave response vs multipole sum",
          "[Kion][ridge][unit]") {
  std::cout << "All-K plane-wave response vs plane-wave multipoles "
               "(K <= 40)\n";
  const auto alpha = PhysConst::alpha;
  // Kmax ~ q r_99 + 5 for convergence: up to ~35 here (Z = 2, 2p at q = 4.5)
  const int K_top = 40;
  const auto titles = std::vector<std::string>{
    "V_T", "V_E+M", "V_L", "X", "A_T", "A_E+M", "A_L", "Y", "Z", "S", "P"};

  // Exact H-like bound states; light (Z = 2) and heavy (Z = 30, where the
  // small component and the relativistic kinematics matter)
  struct Case {
    double Z;
    int n, kappa;
    std::vector<double> Es, qs;
    std::shared_ptr<const Grid> grid;
  };
  const std::vector<Case> cases{
    {2.0,
     1,
     -1,
     {5.0, 10.0, 20.0},
     {2.0, 4.5, 7.0},
     std::make_shared<const Grid>(1.0e-6, 40.0, 6000ul, GridType::loglinear,
                                  4.0)},
    {2.0,
     2,
     1,
     {5.0, 10.0},
     {2.0, 4.5},
     std::make_shared<const Grid>(1.0e-6, 40.0, 6000ul, GridType::loglinear,
                                  4.0)},
    {30.0,
     1,
     -1,
     {1200.0, 2250.0},
     {30.0, 60.0},
     std::make_shared<const Grid>(1.0e-7, 5.0, 8000ul, GridType::loglinear,
                                  2.0)}};

  double worst = 0.0;
  for (const auto &cs : cases) {
    auto Fa = DiracSpinor::exactHlike(cs.n, cs.kappa, cs.grid, cs.Z, alpha);
    Fa.occ_frac() = 1.0;
    const auto orb = Kion::momentum_orbital(Fa, Fa.en());

    const SphericalBessel::JL_table jK_tab(K_top + 1, cs.qs, cs.grid->r());
    auto multipoles = Kion::multipole_operators(*cs.grid, false, &jK_tab, true,
                                                true, true, true, true);
    std::vector<std::pair<std::size_t, double>> q_columns;
    for (std::size_t iq = 0; iq < cs.qs.size(); ++iq) {
      q_columns.emplace_back(iq, cs.qs[iq] * PhysConst::c);
    }
    auto PW = Kion::allocate_formFactors(cs.Es.size(), cs.qs.size(), true, true,
                                         true, true, true);
    const auto [lc_min, lc_max] = Kion::continuum_l_range(Fa, K_top, {});
    for (std::size_t iE = 0; iE < cs.Es.size(); ++iE) {
      const auto ec = cs.Es[iE] + Fa.en();
      const auto waves =
        DiracODE::freeDirac(ec, lc_min, lc_max, cs.grid, alpha);
      Kion::accumulate_multipole_sum(&PW, iE, Fa, 1.0, waves, multipoles, 0,
                                     K_top, q_columns);
    }

    fmt::print("  Z = {:g} n = {} kappa = {:+}:\n", cs.Z, cs.n, cs.kappa);
    fmt::print("  {:>7} {:>6} {:>6} {:>13} {:>13} {:>8}\n", "factor", "E", "q",
               "PW(K<=40)", "IA(all K)", "rel");
    for (std::size_t iE = 0; iE < cs.Es.size(); ++iE) {
      for (std::size_t iq = 0; iq < cs.qs.size(); ++iq) {
        const auto IA =
          Kion::planewave_formFactors(orb, cs.Es[iE], cs.qs[iq], alpha);
        // E and M compared as their sum; FormFactorSet order
        const std::array<double, 11> pw{
          PW[0](iE, iq),  PW[1](iE, iq) + PW[2](iE, iq),
          PW[3](iE, iq),  PW[4](iE, iq),
          PW[5](iE, iq),  PW[6](iE, iq) + PW[7](iE, iq),
          PW[8](iE, iq),  PW[9](iE, iq),
          PW[10](iE, iq), PW[11](iE, iq),
          PW[12](iE, iq)};
        const std::array<double, 11> ia{IA[0],  IA[1] + IA[2], IA[3], IA[4],
                                        IA[5],  IA[6] + IA[7], IA[8], IA[9],
                                        IA[10], IA[11],        IA[12]};
        for (std::size_t i = 0; i < 11; ++i) {
          // relative to the factor, with a floor relative to V_T
          const auto denom = std::max(std::abs(pw[i]), 1.0e-6 * pw[0]);
          const auto rel = std::abs(ia[i] - pw[i]) / denom;
          if (std::abs(pw[i]) > 1.0e-6 * pw[0]) {
            worst = std::max(worst, rel);
          }
          fmt::print("  {:>7} {:6.0f} {:6.1f} {:13.6e} {:13.6e} {:8.1e}\n",
                     titles[i], cs.Es[iE], cs.qs[iq], pw[i], ia[i], rel);
        }
      }
    }
  }
  fmt::print("  worst relative difference: {:.1e}\n", worst);
  REQUIRE(worst < 5.0e-4);
}

//==============================================================================
TEST_CASE("Kion: ridge correction, hydrogen Kmax independence",
          "[Kion][ridge][unit]") {
  std::cout << "Ridge correction: hydrogen 1s on the Bethe ridge\n";
  // Exact all-K distorted-wave values (non-relativistic, partial waves to
  // L = 70) on the ridge q = sqrt(2E): E = 10 au: 0.1884; E = 50 au: 0.0848.
  // With K <= 6 alone: 0.1300 and 0.0220.
  // Local (nuclear) potential for the core: the form factors use the Zeff
  // bound and continuum states (Zeff = 1 from the 1s energy), so only that
  // energy is taken from the core
  Wavefunction wf({4000, 1.0e-6, 80.0, 10.0, GridType::loglinear},
                  {"H", 1, "pointlike"});
  wf.solve_core("Local", "1s1", {}, 1.0e-13, false);
  REQUIRE(wf.core().size() == 1);
  REQUIRE(std::abs(wf.core().front().en() + 0.5) < 1.0e-4);

  const std::vector<double> Egrid{10.0, 50.0};
  const std::vector<double> qgrid{1.5, std::sqrt(20.0), 10.0, 15.0};
  const auto exact = std::array{0.1884, 0.0848};
  const auto ridge_iq = std::array<std::size_t, 2>{1, 2};

  fmt::print("  {:>4} {:>5} {:>8} {:>8} {:>8} {:>8} {:>8}\n", "Kmax", "E", "DW",
             "dK", "total", "exact", "rel");
  for (const auto &[Kmax, tol] : {std::pair{6, 0.03}, {12, 0.01}}) {
    const SphericalBessel::JL_table jK_tab(Kmax + 1, qgrid, wf.grid().r());
    const auto method = Kion::AtomicMethod::Zeff;
    const auto bound = Kion::model_bound_states(*wf.vHF(), method);
    const auto DW = Kion::calculate_formFactors(
      wf.vHF(), bound, {}, 0.0, 1.0e99, false, false, true, Egrid, qgrid, false,
      false, jK_tab, 0, Kmax, true, false, false, false, false, method);
    const auto dK = Kion::calculate_ridge_correction(
      wf.vHF(), bound, 0.0, 1.0e99, Egrid, qgrid, jK_tab, Kmax, true, false,
      false, false, false, 1.0e-4);
    for (std::size_t iE = 0; iE < 2; ++iE) {
      const auto iq = ridge_iq[iE];
      const auto dw = DW[0][0](iE, iq);
      const auto dk = dK[0][0](iE, iq);
      const auto total = dw + dk;
      const auto rel = (total - exact[iE]) / exact[iE];
      fmt::print("  {:4} {:5.0f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:+8.1e}\n",
                 Kmax, Egrid[iE], dw, dk, total, exact[iE], rel);
      REQUIRE(std::abs(rel) < tol);
      REQUIRE(dk > 0.0);
      // Off the ridge (q at 1/7 - 1/3, and 1.5 - 3.4, times the ridge
      // value) the completion is negligible on the scale of the response:
      // below 1e-3 of the on-ridge total. Relative to the (tiny) local
      // response it can reach the percent level, which is the plane-wave
      // K > Kmax content there; printed for information
      for (const auto iq_off : {std::size_t(0), std::size_t(3)}) {
        const auto dk_off = dK[0][0](iE, iq_off);
        const auto dw_off = DW[0][0](iE, iq_off);
        fmt::print("       off-ridge q = {:5.2f}: dK = {:.1e} ({:.1e} of local "
                   "DW, {:.1e} of ridge total)\n",
                   qgrid[iq_off], dk_off, dk_off / dw_off, dk_off / total);
        REQUIRE(std::abs(dk_off) < 1.0e-3 * total);
      }
    }
  }
}
