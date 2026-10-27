#include "MBPT/Sigma2.hpp"
#include "CI/ConfigurationInteraction.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Feynman.hpp"
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
// with S^k from Sk_vwxy(), returned to the m-scheme by WkTable::g().
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

    // Reduced S^k, returned to the m-scheme via WkTable::g()
    Coulomb::WkTable Sk;
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

//==============================================================================
// calculate_Sk stores S^k_vwxy when both states of the pair (v, w), or of
// (x, y), are in the pair basis, or the integral is diagonal (xy in the same
// shells as vw, either j), and all four are in the external basis: exactly those (and only those)
// allowed by the selection rule, with the values of Sk_vwxy
TEST_CASE("Sigma2: S^k pair selection", "[Sigma2][MBPT][unit]") {

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto orbs = DiracSpinor::HlikeBasis("5sp4d", radial_grid, 1.0);

  // core: n <= 2; external: n = 3, 4; pair basis: n = 3
  std::vector<DiracSpinor> core, excited, external, pair_basis;
  for (const auto &Fn : orbs) {
    (Fn.n() <= 2 ? core : excited).push_back(Fn);
    if (Fn.n() == 3 || Fn.n() == 4) {
      external.push_back(Fn);
    }
    if (Fn.n() == 3) {
      pair_basis.push_back(Fn);
    }
  }

  const Coulomb::YkTable yk(orbs);
  Coulomb::QkTable qk;
  qk.fill(orbs, yk, -1, false);

  const auto Sk =
    MBPT::calculate_Sk("false", external, pair_basis, core, excited, qk, -1,
                       false, MBPT::Denominators::RS);

  const auto in = [](const std::vector<DiracSpinor> &set,
                     const DiracSpinor &a) {
    return std::find(set.cbegin(), set.cend(), a) != set.cend();
  };
  const auto same_shell = [](const DiracSpinor &a, const DiracSpinor &b) {
    return a.n() == b.n() && a.l() == b.l();
  };

  // All excited states (n = 3, 4, 5): includes legs outside external
  std::size_t num_selected = 0;
  std::size_t num_wrong = 0;
  double worst_value = 0.0;
  for (const auto &v : excited) {
    for (const auto &w : excited) {
      for (const auto &x : excited) {
        for (const auto &y : excited) {
          const auto all_in = in(external, v) && in(external, w) &&
                              in(external, x) && in(external, y);
          const auto pair_in = (in(pair_basis, v) && in(pair_basis, w)) ||
                               (in(pair_basis, x) && in(pair_basis, y));
          const auto diagonal = (same_shell(x, v) && same_shell(y, w)) ||
                                (same_shell(x, w) && same_shell(y, v));
          const auto [k0, k1] = MBPT::k_minmax_S(v, w, x, y);
          for (int k = k0; k <= k1; ++k) {
            const auto selected = all_in && (pair_in || diagonal) &&
                                  MBPT::Sk_vwxy_SR(k, v, w, x, y);
            if (Sk.contains(k, v, w, x, y) != selected) {
              ++num_wrong;
            }
            if (selected) {
              ++num_selected;
              const auto expected =
                MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, yk.SixJ(),
                              MBPT::Denominators::RS);
              worst_value =
                std::max(worst_value, std::abs(Sk.Q(k, v, w, x, y) - expected));
            }
          }
        }
      }
    }
  }
  REQUIRE(num_selected > 0);
  REQUIRE(num_wrong == 0);
  REQUIRE(worst_value < 1.0e-14);
}

//==============================================================================
// The radial two-body matrix element on the Feynman sub-grid, with the bare
// Coulomb matrix q^k as kernel, is the Coulomb radial integral R^k_abcd:
// checks the kernel convention (both integration measures) and the sub-grid
// quadrature, which the screening correction relies on
TEST_CASE("Sigma2: screening two-body ME", "[Sigma2][MBPT][Feynman]") {

  Wavefunction wf({900, 1.0e-6, 100.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  wf.formBasis({"4spd6fg8h8i",
                30,
                7,
                {1.0e-3},
                40.0,
                "",
                SplineBasis::SplineType::Derevianko,
                false,
                false});

  const auto &orbs = wf.basis();

  // Pairs of orbitals with |n_a - n_c| <= 1 (the others have no new
  // information); only integrals large enough for a relative comparison
  const int max_del_n = 1;
  const double R_min = 1.0e-3;

  fmt::print("two_body_ME(q^k) vs R^k_abcd, |R| > {:.0e}:\n", R_min);
  for (const std::size_t stride : {1ul, 2ul}) {
    // Sub-grid: whole radial grid. No Q*Pi*Q
    const auto num_subgrid_points = (wf.grid().num_points() - 1) / stride + 1;
    const MBPT::Feynman feynman(wf.vHF(), 0, stride, num_subgrid_points, {}, 1,
                                true, false, "", false);

    double max_rel_diff = 0.0;
    long num_compared = 0;
    std::string worst;
    for (int k = 0; k <= feynman.max_k(); ++k) {
      const auto qk_both_measures = feynman.get_qk(k).dri().real();
#pragma omp parallel for collapse(2) reduction(+ : num_compared)
      for (std::size_t ia = 0; ia < orbs.size(); ++ia) {
        for (std::size_t ib = 0; ib < orbs.size(); ++ib) {
          const auto &Fa = orbs[ia];
          const auto &Fb = orbs[ib];
          for (std::size_t ic = ia; ic < orbs.size(); ++ic) {
            const auto &Fc = orbs[ic];
            if (std::abs(Fa.n() - Fc.n()) > max_del_n)
              continue;
            if (!Angular::Ck_kk_SR(k, Fa.kappa(), Fc.kappa()))
              continue;
            for (std::size_t id = ib; id < orbs.size(); ++id) {
              const auto &Fd = orbs[id];
              if (std::abs(Fb.n() - Fd.n()) > max_del_n)
                continue;
              if (!Angular::Ck_kk_SR(k, Fb.kappa(), Fd.kappa()))
                continue;
              const auto R_exact = Coulomb::Rk_abcd(k, Fa, Fb, Fc, Fd);
              if (std::abs(R_exact) < R_min)
                continue;
              const auto R_subgrid =
                MBPT::two_body_ME(qk_both_measures, Fa, Fb, Fc, Fd);
              const auto rel_diff = std::abs(R_subgrid / R_exact - 1.0);
              ++num_compared;
#pragma omp critical(worst_R)
              if (rel_diff > max_rel_diff) {
                max_rel_diff = rel_diff;
                worst = fmt::format("k={} {} {} {} {}: {:.6e} vs {:.6e}", k,
                                    Fa.shortSymbol(), Fb.shortSymbol(),
                                    Fc.shortSymbol(), Fd.shortSymbol(),
                                    R_subgrid, R_exact);
              }
            }
          }
        }
      }
    }
    fmt::print("  stride {}: {} integrals, max rel diff = {:.2e} ({})\n",
               stride, num_compared, max_rel_diff, worst);
    REQUIRE(num_compared > 1000);
    // Stride 1 is the full grid: only the quadrature rule differs. Stride 2
    // tests the sub-grid measure (dr*stride); the quadrature error is then a
    // few parts in 10^3 (high-l integrals, on this coarse grid)
    const auto tolerance = stride == 1 ? 1.0e-12 : 2.0e-2;
    REQUIRE(max_rel_diff < tolerance);
  }
}

//==============================================================================
// Static kernels from the stored Q*Pi*Q on the Feynman contour (as used in
// the CI). (1) Without screening, the kernel is the second-order loop: it
// matches the Goldstone loop diagrams a1 + b1 with Fermi0 denominators,
//   (-1)^k/[k] sum_an Q^k_vnxa Q^k_awny * 2/(e_a - e_n),
// from a complete basis, up to the small (g) components (not in the Green's
// functions here) and the sub-grid quadrature. Checks the Cauchy formula for
// the static limit, the angular factor and the integration measures. Also
// the S^k symmetries (needed for the WkTable). (2) The all-orders chain
// screens less than the loop (the higher orders alternate in sign). (3) The
// hole-particle interaction increases the polarisability, so the screening
// (k = 1: k = 0 is numerically delicate with hole-particle)
TEST_CASE("Sigma2: Feynman static loop", "[Sigma2][MBPT][Feynman]") {

  Wavefunction wf({800, 1.0e-6, 80.0, 20.0, "loglinear"}, {"Mg", -1, "Fermi"});
  wf.solve_core("HartreeFock", "[Ne]");
  wf.solve_valence("3s");
  // Full spline set: complete basis (within the box) for the Goldstone sum
  wf.formBasis({"30spdf",
                30,
                7,
                {1.0e-4},
                40.0,
                "",
                SplineBasis::SplineType::Derevianko,
                false,
                false});

  // Named references (not structured bindings), so the lambdas below can
  // capture them in C++17
  const auto core_excited =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
  const auto &core = core_excited.first;
  const auto &excited = core_excited.second;

  // As Correlations: sub-grid over r = 1e-4 to 30
  const std::size_t stride = 2;
  const auto i0 = wf.grid().getIndex(1.0e-4);
  const auto num_subgrid_points = (wf.grid().getIndex(30.0) - i0) / stride + 1;
  const auto feynman = [&](MBPT::Screening screening, MBPT::HoleParticle hp) {
    MBPT::FeynmanOptions options;
    options.screening = screening;
    options.hole_particle = hp;
    options.max_l_internal = 3;
    options.omre = MBPT::best_omre(wf.core(), wf.valence());
    return MBPT::Feynman(wf.vHF(), i0, stride, num_subgrid_points, options, 1,
                         false, false, "", false);
  };
  auto second_order =
    feynman(MBPT::Screening::exclude, MBPT::HoleParticle::exclude);
  auto all_orders =
    feynman(MBPT::Screening::include, MBPT::HoleParticle::exclude);
  auto all_orders_hp =
    feynman(MBPT::Screening::include, MBPT::HoleParticle::include);
  // The polarisation loop does not depend on screening: formed once
  const auto pi_wk = second_order.polarisation_wk();
  second_order.form_qpiq(pi_wk);
  all_orders.form_qpiq(pi_wk);
  all_orders_hp.form_qpiq("");

  const auto dQ2 = second_order.dQ_0_k();
  const auto dQ_all = all_orders.dQ_0_k();
  const auto dQ_all_hp = all_orders_hp.dQ_0_k();

  // Goldstone loop diagrams a1 + b1, Fermi0 denominators. Q^k_awny is
  // looked up as Q^k_wayn (same integral), so only core-basis y^k are needed
  const Coulomb::YkTable yk(core, wf.basis());
  const auto loop = [&](int k, const DiracSpinor &v, const DiracSpinor &w,
                        const DiracSpinor &x, const DiracSpinor &y) {
    double sum = 0.0;
    for (const auto &a : core) {
      for (const auto &n : excited) {
        if (!Angular::Ck_kk_SR(k, a.kappa(), n.kappa()))
          continue;
        sum +=
          yk.Q(k, v, n, x, a) * yk.Q(k, w, a, y, n) * 2.0 / (a.en() - n.en());
      }
    }
    return Angular::neg1pow(k) / (2.0 * k + 1.0) * sum;
  };

  // Valence-like external legs: n = 3, and 4s
  std::vector<DiracSpinor> valence;
  for (const auto &Fn : excited) {
    if (Fn.n() == 3 || (Fn.n() == 4 && Fn.l() == 0)) {
      valence.push_back(Fn);
    }
  }
  REQUIRE(valence.size() == 6);

  double max_diff = 0.0;
  double max_value = 0.0;
  double max_asymmetry = 0.0;
  for (const auto &v : valence) {
    for (const auto &w : valence) {
      for (const auto &x : valence) {
        for (const auto &y : valence) {
          const auto [k0, k1] = Coulomb::k_minmax_Q(v, w, x, y);
          for (int k = k0; k <= k1; k += 2) {
            const auto &kernel = dQ2.at(std::size_t(k));
            const auto S_feynman =
              MBPT::Sigma2::S_Sigma2_screen(k, v, w, x, y, kernel);
            const auto S_goldstone = loop(k, v, w, x, y);
            max_diff = std::max(max_diff, std::abs(S_feynman - S_goldstone));
            max_value = std::max(max_value, std::abs(S_goldstone));
            // S^k_vwxy = S^k_wvyx = S^k_xyvw
            const auto S_wvyx =
              MBPT::Sigma2::S_Sigma2_screen(k, w, v, y, x, kernel);
            const auto S_xyvw =
              MBPT::Sigma2::S_Sigma2_screen(k, x, y, v, w, kernel);
            max_asymmetry =
              std::max({max_asymmetry, std::abs(S_wvyx - S_feynman),
                        std::abs(S_xyvw - S_feynman)});
          }
        }
      }
    }
  }
  fmt::print("Second-order static kernel vs Goldstone a1+b1 (Fermi0): max "
             "diff = {:.2e} (of max |S| = {:.2e})\n",
             max_diff / max_value, max_value);
  fmt::print("S^k symmetry (vwxy = wvyx = xyvw): max diff = {:.2e}\n",
             max_asymmetry);
  // Small components: at most ~(alpha Z)^2/4 ~ 2e-3 of the density for Mg
  REQUIRE(max_diff / max_value < 2.0e-3);
  REQUIRE(max_asymmetry < 1.0e-12);

  // (2) All orders vs second order, k = 0, 3s3s
  const auto &v3s = *DiracSpinor::find(3, -1, valence);
  const auto S2_3s =
    MBPT::Sigma2::S_Sigma2_screen(0, v3s, v3s, v3s, v3s, dQ2.at(0));
  const auto Sall_3s =
    MBPT::Sigma2::S_Sigma2_screen(0, v3s, v3s, v3s, v3s, dQ_all.at(0));
  fmt::print("<3s3s|dQ^0|3s3s>: second order {:.4e}, all orders {:.4e}\n",
             S2_3s, Sall_3s);
  REQUIRE(S2_3s < 0.0);
  REQUIRE(Sall_3s / S2_3s > 0.0);
  REQUIRE(Sall_3s / S2_3s < 1.0);

  // (3) Hole-particle, k = 1, <3s3p|dQ^1|3p3s>
  const auto &v3p = *DiracSpinor::find(3, 1, valence);
  const auto S1_all =
    MBPT::Sigma2::S_Sigma2_screen(1, v3s, v3p, v3p, v3s, dQ_all.at(1));
  const auto S1_all_hp =
    MBPT::Sigma2::S_Sigma2_screen(1, v3s, v3p, v3p, v3s, dQ_all_hp.at(1));
  fmt::print("<3s3p|dQ^1|3p3s>: all orders {:.4e}, all orders + hp {:.4e}\n",
             S1_all, S1_all_hp);
  REQUIRE(std::abs(S1_all_hp) > std::abs(S1_all));
}

//==============================================================================
// With screening kernels, Sk_vwxy takes the static loop diagrams (a1, b1)
// from the kernels and keeps only the energy dependence of their
// denominators at second order. So, without fk, for any denominators:
//   S^k(kernels) - S^k = <kernel loop> - <Goldstone a1 + b1, Fermi0>
// with the Fermi0 loop (-1)^k/[k] sum_an Q^k_vnxa Q^k_awny 2/(e_a - e_n).
// Any kernel tests this: the Coulomb kernel q^k is used
TEST_CASE("Sigma2: Feynman loop in S^k", "[Sigma2][MBPT][unit]") {

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{300, 1.0e-4, 200.0, 40.0, GridType::loglinear});
  const auto orbs = DiracSpinor::HlikeBasis("4sp3d", radial_grid, 1.0);

  // "core" = n<=2 (1s, 2s, 2p); "excited" = all others
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fn : orbs) {
    (Fn.n() <= 2 ? core : excited).push_back(Fn);
  }

  const Coulomb::YkTable yk(orbs);
  const auto &sjt = yk.SixJ();
  Coulomb::QkTable qk;
  qk.fill(orbs, yk, -1, false);

  // Kernels: q^k_ij dr_i dr_j on the full grid
  const int max_k = 4;
  const auto size = radial_grid->num_points();
  std::vector<MBPT::RMatrix> kernels(std::size_t(max_k + 1),
                                     MBPT::RMatrix{0, 1, size, radial_grid});
  for (std::size_t i = 0; i < size; ++i) {
    for (std::size_t j = 0; j < size; ++j) {
      const auto r_less = std::min(radial_grid->r(i), radial_grid->r(j));
      const auto r_greater = std::max(radial_grid->r(i), radial_grid->r(j));
      const auto dr_dr = radial_grid->drdu(i) * radial_grid->drdu(j) *
                         radial_grid->du() * radial_grid->du();
      for (int k = 0; k <= max_k; ++k) {
        kernels[std::size_t(k)](i, j) =
          std::pow(r_less / r_greater, k) / r_greater * dr_dr;
      }
    }
  }

  // Goldstone loop diagrams a1 + b1, Fermi0 denominators
  const auto loop_Fermi0 = [&](int k, const DiracSpinor &v,
                               const DiracSpinor &w, const DiracSpinor &x,
                               const DiracSpinor &y) {
    double sum = 0.0;
    for (const auto &a : core) {
      for (const auto &n : excited) {
        if (!Angular::Ck_kk_SR(k, a.kappa(), n.kappa()))
          continue;
        sum +=
          qk.Q(k, v, n, x, a) * qk.Q(k, a, w, n, y) * 2.0 / (a.en() - n.en());
      }
    }
    return Angular::neg1pow(k) / (2.0 * k + 1.0) * sum;
  };

  // External legs: n = 3, 4 (s and p_1/2), and 3d_3/2
  std::vector<DiracSpinor> external;
  for (const auto &Fn : excited) {
    if (Fn.kappa() == -1 || Fn.kappa() == 1 || Fn.kappa() == 2) {
      external.push_back(Fn);
    }
  }

  const double E0 = -0.15;
  for (const auto denominators :
       {MBPT::Denominators::Fermi0, MBPT::Denominators::DFK,
        MBPT::Denominators::BW, MBPT::Denominators::RS}) {
    double max_diff = 0.0;
    double max_value = 0.0;
    for (const auto &v : external) {
      for (const auto &w : external) {
        for (const auto &x : external) {
          for (const auto &y : external) {
            const auto [k0, k1] = MBPT::k_minmax_S(v, w, x, y);
            for (int k = k0; k <= std::min(k1, max_k); ++k) {
              if (!MBPT::Sk_vwxy_SR(k, v, w, x, y))
                continue;
              const auto Sk_Goldstone = MBPT::Sk_vwxy(
                k, v, w, x, y, qk, core, excited, sjt, denominators, {}, E0);
              const auto Sk_Feynman =
                MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, sjt,
                              denominators, {}, E0, &kernels);
              const auto expected = MBPT::Sigma2::S_Sigma2_screen(
                                      k, v, w, x, y, kernels[std::size_t(k)]) -
                                    loop_Fermi0(k, v, w, x, y);
              max_diff = std::max(
                max_diff, std::abs((Sk_Feynman - Sk_Goldstone) - expected));
              max_value = std::max(max_value, std::abs(expected));
            }
          }
        }
      }
    }
    REQUIRE(max_value > 1.0e-3);
    REQUIRE(max_diff < 1.0e-12 * max_value);
  }
}

//==============================================================================
// With all-orders screening and hole-particle off, the Feynman loop is the
// second-order loop. With a complete basis (internal l up to 3 in both),
// and no fk (it screens the whole Goldstone loop, but only the
// energy-dependence part of the Feynman one), CI with Feynman
// loops must reproduce regular Goldstone Sigma_2, up to the numerical
// differences: mainly the sub-grid quadrature (Mg, this grid: max diff 2.4e-3
// of the Sigma_2 shift at stride 2, 5e-4 at stride 1)
TEST_CASE("Sigma2: Feynman second order vs Goldstone",
          "[Sigma2][CI][integration]") {

  Wavefunction wf({800, 1.0e-6, 80.0, 20.0, "loglinear"}, {"Mg", -1, "Fermi"});
  wf.solve_core("HartreeFock", "[Ne]");
  wf.solve_valence("3s");
  // Full spline set: complete basis (within the box)
  wf.formBasis({"30spdf",
                30,
                7,
                {1.0e-4},
                40.0,
                "",
                SplineBasis::SplineType::Derevianko,
                false,
                false});
  // Feynman correlation potential: provides the sub-grid and parameters
  wf.formSigma(1, 1.0e-4, 30.0, 2, false, false, false, 0, {}, {}, {}, false,
               "", true, false, false, 3,
               MBPT::best_omre(wf.core(), wf.valence()));

  // Lowest level of each of J = 0 (even), J = 0, 1, 2 (odd), for 3sp
  const auto levels = [&wf](const std::string &options) {
    const IO::InputBlock input{
      "CI", "ci_basis = 3sp; cis2_basis = 3sp; J+ = 0; J- = 0, 1, 2; "
            "num_solutions = 1; fk = false; qk_file = false; sk_file = "
            "false; ci_file = false; bk_file = false; " +
              options};
    std::vector<double> energies;
    for (const auto &level : CI::configuration_interaction(input, wf).levels) {
      energies.push_back(level.energy(0));
    }
    return energies;
  };
  const auto E_no_sigma2 = levels("sigma2 = false;");
  const auto E_Goldstone = levels("sigma2 = true; Feynman_screen = false;");
  const auto E_Feynman =
    levels("sigma2 = true; Feynman_screen = true; all_order_screen = false; "
           "hp_screen = false;");
  REQUIRE(E_Goldstone.size() == 4);
  REQUIRE(E_Feynman.size() == 4);

  double max_diff = 0.0;
  double max_sigma2 = 0.0;
  for (std::size_t i = 0; i < E_Goldstone.size(); ++i) {
    max_diff = std::max(max_diff, std::abs(E_Feynman[i] - E_Goldstone[i]));
    max_sigma2 =
      std::max(max_sigma2, std::abs(E_Goldstone[i] - E_no_sigma2[i]));
  }
  fmt::print("CI levels, Feynman second-order loop vs Goldstone: max diff "
             "= {:.3f} /cm (Sigma_2 shifts up to {:.1f} /cm)\n",
             max_diff * PhysConst::Hartree_invcm,
             max_sigma2 * PhysConst::Hartree_invcm);
  REQUIRE(max_sigma2 > 0.0);
  REQUIRE(max_diff < 1.0e-2 * max_sigma2);
}

//==============================================================================
// The Feynman loop through the CI. (1) With a single CSF (3s^2), without
// fk and hole-particle, the energy shift is the change in <3s3s|S^0|3s3s>:
// the all-orders chain screens less than the second-order loop, so the
// valence repulsion is less reduced and the level moves up. (2) All-orders
// screening with hole-particle must not split LS multiplets: the Lande
// interval ratio of the 3s3p 3P_J levels is unchanged (as for the "Mg fine
// structure" test). (3) Without a Feynman correlation potential, the option
// does nothing (with a warning)
TEST_CASE("Sigma2: Mg Feynman screening", "[Sigma2][CI][integration]") {

  Wavefunction wf({2000, 1.0e-6, 120.0, 40.0, "loglinear"},
                  {"Mg", -1, "Fermi"});
  wf.solve_core("HartreeFock", "[Ne]");
  wf.solve_valence("3s");
  wf.formBasis({"15spdf", 30, 7, {1.0e-4}, 40.0});

  const auto common = std::string("qk_file = false; sk_file = false; "
                                  "ci_file = false; bk_file = false;");
  const auto E_3s2 = [&wf, &common](bool screening) {
    const IO::InputBlock input{
      "CI", "ci_basis = 3s; J+ = 0; num_solutions = 1; sigma2 = true; "
            "fk = false; hp_screen = false; Feynman_screen = " +
              std::string(screening ? "true; " : "false; ") + common};
    const auto levels = CI::configuration_interaction(input, wf).levels;
    REQUIRE(levels.size() == 1);
    REQUIRE(levels.front().info(0).config == "3s^2");
    return levels.front().energy(0);
  };

  // (3) No correlation potential yet
  const auto E_3s2_no_Sigma = E_3s2(false);
  std::cout << "\nExpect one warning below (deliberate): screening without a "
               "Feynman correlation potential\n";
  const auto E_3s2_no_Sigma_screening = E_3s2(true);
  REQUIRE(E_3s2_no_Sigma_screening == E_3s2_no_Sigma);

  // Feynman correlation potential: provides the sub-grid and parameters
  wf.formSigma(1, 1.0e-4, 30.0, 4, false, false, false, 0, {}, {}, {}, false,
               "", true, true, false, 3,
               MBPT::best_omre(wf.core(), wf.valence()));
  REQUIRE(wf.Sigma() != nullptr);
  REQUIRE(wf.Sigma()->feynman() != nullptr);

  // (1) Single CSF
  const auto E_3s2_sigma2 = E_3s2(false);
  const auto E_3s2_screened = E_3s2(true);
  fmt::print("\nMg 3s^2 (single CSF), Feynman screening vs second order: "
             "{:+.3f} /cm\n",
             (E_3s2_screened - E_3s2_sigma2) * PhysConst::Hartree_invcm);
  REQUIRE(E_3s2_sigma2 == E_3s2_no_Sigma);
  REQUIRE(E_3s2_screened > E_3s2_sigma2);

  // (2) Full CI in 6spd. Returns {3s2 1S0, 3s3p 3P0, 3P1, 3P2}
  const auto lowest_levels = [&wf, &common](bool screening) {
    const IO::InputBlock input{
      "CI", "ci_basis = 6spd; cis2_basis = 4spd; J+ = 0; J- = 0, 1, 2;"
            "num_solutions = 1; sigma1 = true; sigma2 = true; "
            "Feynman_screen = " +
              std::string(screening ? "true; " : "false; ") + common};
    const auto levels = CI::configuration_interaction(input, wf).levels;
    REQUIRE(levels.size() == 4);
    std::array<double, 4> energies;
    for (const auto &level : levels) {
      if (level.parity() == 1) {
        REQUIRE(level.info(0).config == "3s^2");
        energies.at(0) = level.energy(0);
      } else {
        REQUIRE(level.info(0).config == "3s3p");
        energies.at(std::size_t(level.twoJ() / 2) + 1) = level.energy(0);
      }
    }
    return energies;
  };

  const auto E_sigma2 = lowest_levels(false);
  const auto E_screened = lowest_levels(true);

  const auto lande_ratio = [](const std::array<double, 4> &E) {
    return (E[3] - E[2]) / (E[2] - E[1]);
  };

  fmt::print("Mg lowest levels: shift from second-order Sigma_2 (/cm)\n");
  const std::array<std::string, 4> names{"3s2 1S0", "3s3p 3P0", "3s3p 3P1",
                                         "3s3p 3P2"};
  for (std::size_t i = 0; i < 4; ++i) {
    fmt::print("{:8s} {:14.8f} {:+10.3f}\n", names[i], E_sigma2[i],
               (E_screened[i] - E_sigma2[i]) * PhysConst::Hartree_invcm);
  }
  fmt::print("3P Lande ratio: {:.6f} vs {:.6f}\n", lande_ratio(E_screened),
             lande_ratio(E_sigma2));
  REQUIRE(lande_ratio(E_screened) ==
          Approx(lande_ratio(E_sigma2)).epsilon(1.0e-3));
}
