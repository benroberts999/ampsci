#include "MBPT/Sigma2.hpp"
#include "CI/ConfigurationInteraction.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <string>
#include <vector>

//==============================================================================
// Independent (m-scheme) check of the reduced two-body correlation integrals.
// Builds Sigma2_vwxy by brute-force summation over magnetic substates of the
// five Goldstone diagrams (a, b, c1, c2, d), using m-scheme Coulomb integrals
// g (from the Qk table), and compares against the angular-reduced form
//   Sigma2_vwxy = sum_k A^k_vwxy S^k_vwxy ,
// with S^k from Sk_vwxy(), returned to the m-scheme by LkTable::g().
// RS denominators: each is averaged with that of the bra-ket partner diagram,
// 1/D -> (1/D + 1/D')/2, as in Sk_vwxy().
// The core includes a p shell, so core and excited j differ by both even and
// odd amounts: the phase of the c1 and c2 box diagrams depends on j_a + j_n.
TEST_CASE("Sigma2: m-scheme check", "[Sigma2][MBPT][unit]") {

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;

  // H-like orbitals: pure angular/algebra test - no HF required
  const auto orbs = DiracSpinor::HlikeBasis("4sp3d", radial_grid, zeff);

  // "core" = n<=2 (1s, 2s, 2p); "excited" = all others
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fn : orbs) {
    (Fn.n() <= 2 ? core : excited).push_back(Fn);
  }
  REQUIRE(core.size() == 4);

  const Coulomb::YkTable yk(orbs);
  const auto &sjt = yk.SixJ();
  Coulomb::QkTable qk;
  qk.fill(orbs, yk, -1, false);

  const auto find = [&orbs](int n, int kappa) {
    const auto Fn = DiracSpinor::find(n, kappa, orbs);
    REQUIRE(Fn != nullptr);
    return *Fn;
  };

  // Brute-force m-scheme Sigma2_vwxy: diagrams a, b, c1, c2, d, with
  // internal magnetic substates summed explicitly
  const auto sigma2_brute = [&](const DiracSpinor &v, const DiracSpinor &w,
                                const DiracSpinor &x, const DiracSpinor &y,
                                int tmv, int tmw, int tmx, int tmy) {
    double sigma = 0.0;
    for (const auto &a : core) {
      for (const auto &n : excited) {
        const auto inv_de_a = 0.5 * (1.0 / (x.en() + a.en() - v.en() - n.en()) +
                                     1.0 / (w.en() + a.en() - y.en() - n.en()));
        const auto inv_de_b = 0.5 * (1.0 / (y.en() + a.en() - w.en() - n.en()) +
                                     1.0 / (v.en() + a.en() - x.en() - n.en()));
        const auto inv_de_c1 =
          0.5 * (1.0 / (y.en() + a.en() - v.en() - n.en()) +
                 1.0 / (w.en() + a.en() - x.en() - n.en()));
        const auto inv_de_c2 =
          0.5 * (1.0 / (x.en() + a.en() - w.en() - n.en()) +
                 1.0 / (v.en() + a.en() - y.en() - n.en()));
        for (int tma = -a.twoj(); tma <= a.twoj(); tma += 2) {
          for (int tmn = -n.twoj(); tmn <= n.twoj(); tmn += 2) {
            // a: g_vnxa (g_awny - g_awyn) - g_vnax g_awny
            const auto g_awny = qk.g(a, w, n, y, tma, tmw, tmn, tmy);
            const auto g_awyn = qk.g(a, w, y, n, tma, tmw, tmy, tmn);
            sigma += (qk.g(v, n, x, a, tmv, tmn, tmx, tma) * (g_awny - g_awyn) -
                      qk.g(v, n, a, x, tmv, tmn, tma, tmx) * g_awny) *
                     inv_de_a;
            // b: g_vaxn (g_nway - g_nwya) - g_vanx g_nway
            const auto g_nway = qk.g(n, w, a, y, tmn, tmw, tma, tmy);
            const auto g_nwya = qk.g(n, w, y, a, tmn, tmw, tmy, tma);
            sigma += (qk.g(v, a, x, n, tmv, tma, tmx, tmn) * (g_nway - g_nwya) -
                      qk.g(v, a, n, x, tmv, tma, tmn, tmx) * g_nway) *
                     inv_de_b;
            // c1: -g_vnay g_awxn
            sigma -= qk.g(v, n, a, y, tmv, tmn, tma, tmy) *
                     qk.g(a, w, x, n, tma, tmw, tmx, tmn) * inv_de_c1;
            // c2: -g_vany g_nwxa
            sigma -= qk.g(v, a, n, y, tmv, tma, tmn, tmy) *
                     qk.g(n, w, x, a, tmn, tmw, tmx, tma) * inv_de_c2;
          }
        }
      }
    }
    // d: g_vwab g_abxy
    for (const auto &a : core) {
      for (const auto &b : core) {
        const auto inv_de_d = 0.5 * (1.0 / (a.en() + b.en() - v.en() - w.en()) +
                                     1.0 / (a.en() + b.en() - x.en() - y.en()));
        for (int tma = -a.twoj(); tma <= a.twoj(); tma += 2) {
          for (int tmb = -b.twoj(); tmb <= b.twoj(); tmb += 2) {
            sigma += qk.g(v, w, a, b, tmv, tmw, tma, tmb) *
                     qk.g(a, b, x, y, tma, tmb, tmx, tmy) * inv_de_d;
          }
        }
      }
    }
    return sigma;
  };

  // Test tuples {v, w, x, y}: mixed-j pairs, both parities of l_v + l_x
  // (Coulomb-forbidden-parity k), v = x and w = y, and external states of
  // different energy (tests which external energy enters each denominator)
  const auto s3 = find(3, -1), s4 = find(4, -1);
  const auto p3 = find(3, 1), P3 = find(3, -2), p4 = find(4, 1);
  const auto d3 = find(3, 2), D3 = find(3, -3);
  using Four = std::array<DiracSpinor, 4>;
  const std::vector<Four> tuples{{s3, s3, P3, P3},
                                 {s3, p3, P3, s3},
                                 {p3, D3, P3, d3},
                                 {s4, P3, d3, p4},
                                 {P3, s4, P3, s4}};

  fmt::print("{:16s} {:>13s} {:>13s} {:>9s}\n", "Sigma2_vwxy", "expected",
             "found", "eps");
  for (const auto &[v, w, x, y] : tuples) {

    // Reduced S^k, returned to the m-scheme via LkTable::g()
    Coulomb::LkTable Sk;
    const auto [k0, kI] = MBPT::k_minmax_S(v, w, x, y);
    REQUIRE(kI >= k0);
    for (int k = k0; k <= kI; ++k) {
      Sk.add(k, v, w, x, y,
             MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, sjt,
                           MBPT::Denominators::RS));
    }

    // Compare brute-force and reduced, for every magnetic substate combo
    double worst_del = 0.0, scale = 0.0;
    double x_expect = 0.0, x_found = 0.0;
    for (int tmv = -v.twoj(); tmv <= v.twoj(); tmv += 2) {
      for (int tmw = -w.twoj(); tmw <= w.twoj(); tmw += 2) {
        for (int tmx = -x.twoj(); tmx <= x.twoj(); tmx += 2) {
          for (int tmy = -y.twoj(); tmy <= y.twoj(); tmy += 2) {
            const auto sb = sigma2_brute(v, w, x, y, tmv, tmw, tmx, tmy);
            const auto sr = Sk.g(v, w, x, y, tmv, tmw, tmx, tmy);
            worst_del = std::max(worst_del, std::abs(sb - sr));
            if (std::abs(sb) > scale) {
              scale = std::abs(sb);
              x_expect = sb;
              x_found = sr;
            }
          }
        }
      }
    }
    REQUIRE(scale > 0.0);
    const auto eps = worst_del / scale;
    const auto label = "S_" + v.shortSymbol() + w.shortSymbol() +
                       x.shortSymbol() + y.shortSymbol();
    fmt::print("{:16s} {:13.6e} {:13.6e} {:9.1e}\n", label, x_expect, x_found,
               eps);
    REQUIRE(eps < 1.0e-10);
  }
}

//==============================================================================
// Sigma2 is the screening of the valence-valence Coulomb interaction by the
// core. Non-relativistically it is independent of spin, so it cannot split an
// LS multiplet: it changes the Mg 3s3p 3P_J fine structure only through small
// relativistic corrections, and leaves the Lande interval ratio
// (E_2 - E_1)/(E_1 - E_0) essentially unchanged. An error in the angular
// reduction of Sigma2 generally violates this.
TEST_CASE("Sigma2: Mg fine structure", "[Sigma2][CI][integration]") {

  Wavefunction wf({2000, 1.0e-6, 120.0, 40.0, "loglinear"},
                  {"Mg", -1, "Fermi"});
  wf.solve_core("HartreeFock", "[Ne]");
  wf.formBasis({"15spdf", 30, 7, {1.0e-4}, 40.0});

  // 3s3p 3P_J energies, J = 0, 1, 2, from CI+MBPT (with or without Sigma2)
  const auto triplet_P = [&wf](bool include_sigma2) {
    const IO::InputBlock input{
      "CI", "ci_basis = 6spd; cis2_basis = 4spd; J- = 0, 1, 2;"
            "num_solutions = 1; sigma1 = true; sigma2 = " +
              std::string(include_sigma2 ? "true" : "false") +
              "; qk_file = false; sk_file = false; ci_file = false;"
              "bk_file = false;"};
    const auto levels = CI::configuration_interaction(input, wf).levels;
    REQUIRE(levels.size() == 3);
    std::array<double, 3> energies;
    for (const auto &level : levels) {
      REQUIRE(level.parity() == -1);
      REQUIRE(level.info(0).config == "3s3p");
      REQUIRE(level.info(0).S2 == Approx(2.0).margin(0.05));
      energies.at(std::size_t(level.twoJ() / 2)) = level.energy(0);
    }
    return energies;
  };

  const auto E_sigma1 = triplet_P(false);
  const auto E_sigma2 = triplet_P(true);

  const auto lande_ratio = [](const std::array<double, 3> &E) {
    return (E[2] - E[1]) / (E[1] - E[0]);
  };

  fmt::print("\nMg 3s3p 3P_J fine structure, CI+MBPT (cm^-1)\n");
  fmt::print("{:12s} {:>12s} {:>12s}\n", "", "Sigma1", "Sigma1+2");
  fmt::print("{:12s} {:12.4f} {:12.4f}\n", "E(3P1-3P0)",
             (E_sigma1[1] - E_sigma1[0]) * PhysConst::Hartree_invcm,
             (E_sigma2[1] - E_sigma2[0]) * PhysConst::Hartree_invcm);
  fmt::print("{:12s} {:12.4f} {:12.4f}\n", "E(3P2-3P1)",
             (E_sigma1[2] - E_sigma1[1]) * PhysConst::Hartree_invcm,
             (E_sigma2[2] - E_sigma2[1]) * PhysConst::Hartree_invcm);
  fmt::print("{:12s} {:12.6f} {:12.6f}\n", "ratio", lande_ratio(E_sigma1),
             lande_ratio(E_sigma2));

  REQUIRE(lande_ratio(E_sigma2) ==
          Approx(lande_ratio(E_sigma1)).epsilon(1.0e-3));
}
