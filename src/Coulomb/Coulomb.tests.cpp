#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace UnitTest {

//==============================================================================
// This is a Naiive (slow, but simple + correct) implementation of yk
// This forms the "baseline" to compare against
inline std::vector<double> yk_naive(const DiracSpinor &Fa,
                                    const DiracSpinor &Fb, int k);

// This takes two sets of spinors, and calculates each y^k_ab two ways (using
// YkTable, and using Coulomb:: functions). Also checks if y_ab = y_ba.
// Compares these, returns the worst y1(r)-y2(r). This is a test of the YkTable
// and symmetry, but not the Coulomb:: formula
inline double check_ykab_Tab(const std::vector<DiracSpinor> &a,
                             const std::vector<DiracSpinor> &b,
                             const Coulomb::YkTable &Yab);

// This takes a sets of spinors, and calculates each y^k_ab two ways (using
// Coulomb:: functions, and the above 'Naiive' formula). Compares
// these, returns the worst y1(r)-y2(r).
// Calculates only those with n values that differ by less than max_del_n.
// This is a test of the Coulomb:: formulas
inline std::vector<double> check_ykab(const std::vector<DiracSpinor> &orbs,
                                      int max_del_n = 99);

// This takes in a set of orbitals, and calculates a subset of the R^k radial
// integrals a number of ways. Subset only, otherwise takes too long. Calculates
// only those with n values that differ by less than max_del_n.
// Calculates R^k several ways: Using YkTable, and a few different Coulomb::
// functions. Compares them all, and resturns worst comparison.
inline double check_Rkabcd(const std::vector<DiracSpinor> &orbs,
                           int max_del_n = 99);

} // namespace UnitTest

//==============================================================================

//==============================================================================
//! First, test quadrature integration method:
TEST_CASE("Coulomb: quad integrate", "[Coulomb][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: quad integrate\n";

  // Define a wavefunction-like function, func,
  auto func = [](double x) {
    return x * std::exp(-0.2 * x) * (1.0 + 2.0 * std::sin(x));
  };
  // that has an exact integral, Intfunc = \int_0^x(func):
  auto Intfunc = [](double x) {
    return -(5.0 / 169.0) * std::exp(-0.2 * x) *
           (169.0 * (5.0 + x) + 5.0 * (5.0 + 13.0 * x) * std::cos(x) +
            (13.0 * x - 60.0) * std::sin(x));
  };

  // Perform the numerical integrals using a different grids, compare to exact
  const auto pts_lst = std::vector<std::size_t>{1000, 2000};
  for (const auto pts : pts_lst) {
    const Grid grll(1.0e-6, 100.0, pts, GridType::loglinear, 10);
    const Grid grlog(1.0e-6, 100.0, pts, GridType::logarithmic, 0);

    std::vector<double> vll, vlog;
    vll.reserve(grll.num_points());
    vlog.reserve(grlog.num_points());
    for (const auto &r : grll.r())
      vll.push_back(func(r));
    for (const auto &r : grlog.r())
      vlog.push_back(func(r));

    // numerical integration (on grid):
    const auto intll = NumCalc::integrate(grll.du(), 0, 0, vll, grll.drdu());
    const auto intlog =
        NumCalc::integrate(grlog.du(), 0, 0, vlog, grlog.drdu());

    // Account for possibility r0, rmax slightly different:
    const auto exactll = Intfunc(grll.r().back()) - Intfunc(grll.r().front());
    const auto exactlog =
        Intfunc(grlog.r().back()) - Intfunc(grlog.r().front());

    const auto error1 = std::abs((intll - exactll) / exactll);
    const auto eps1 = 1.0e-13;
    REQUIRE(error1 < eps1);

    const auto error2 = std::abs((intlog - exactlog) / exactlog);
    const auto eps2 = 1.0e-4 * (500.0 / double(pts));
    REQUIRE(error2 < eps2);
  }
}

//==============================================================================
TEST_CASE("Coulomb: yk,Rk formulas", "[Coulomb][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: yk,Rk formulas\n";

  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  // Compare yk formula to "Naive" formula
  using namespace qip::overloads;
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      if (Fb < Fa)
        continue;
      auto [k0, ki] = Angular::kminmax_Ck(Fa.kappa(), Fb.kappa());
      for (int k = k0; k <= ki; k += 2) {
        auto yk1 = UnitTest::yk_naive(Fa, Fb, k);
        auto yk2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto del = std::abs(qip::compare(yk1, yk2).first);
        REQUIRE(del < 1.0e-14);
      }
    }
  }

  // Check Rk_abcd formula
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      for (const auto &Fc : orbs) {
        if (Fc < Fa)
          continue;
        for (const auto &Fd : orbs) {
          if (Fd < Fb)
            continue;
          auto [k0, ki] = Coulomb::k_minmax_Q(Fa, Fb, Fc, Fd);
          for (int k = k0; k <= ki; k += 2) {
            const auto Rabcd = Coulomb::Rk_abcd(Fa, Fb, Fc, Fd, k);
            const auto ykbd = Coulomb::yk_ab(Fb, Fd, k);
            const auto Rabcd2 = Fa * (ykbd * Fc);
            const auto Rabcd3 = Coulomb::Rk_abcd(Fa, Fc, ykbd);
            const auto eps = std::abs((Rabcd - Rabcd2) / Rabcd);
            const auto eps2 = std::abs((Rabcd - Rabcd3) / Rabcd);
            REQUIRE(eps < 1.0e-14);
            REQUIRE(eps2 < 1.0e-14);

            // check symmetries:
            const auto Rbadc = Coulomb::Rk_abcd(Fb, Fa, Fd, Fc, k);
            const auto Rcbad = Coulomb::Rk_abcd(Fc, Fa, ykbd);
            const auto eps3 = std::abs((Rbadc - Rabcd) / Rabcd);
            const auto eps4 = std::abs((Rcbad - Rabcd) / Rabcd);
            REQUIRE(eps3 < 1.0e-14);
            REQUIRE(eps4 < 1.0e-14);

            // Check R^K(v)_bcd formula
            auto Rka_bcd = Coulomb::Rkv_bcd(Fa.kappa(), Fc, ykbd);
            const auto Rka_bcd2 = Coulomb::Rkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
            const auto Rabcd4 = Fa * Rka_bcd;
            const auto Rabcd5 = Fa * Rka_bcd2;
            const auto eps5 = std::abs((Rabcd4 - Rabcd) / Rabcd);
            const auto eps6 = std::abs((Rabcd5 - Rabcd) / Rabcd);
            REQUIRE(eps5 < 1.0e-14);
            REQUIRE(eps6 < 1.0e-14);
            Coulomb::Rkv_bcd(&Rka_bcd, Fc, ykbd);
            const auto Rabcd6 = Fa * Rka_bcd;
            const auto eps7 = std::abs((Rabcd6 - Rabcd) / Rabcd);
            REQUIRE(eps7 < 1.0e-14);

            // check Qk
            const auto qk1 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
            const auto qk2 = Rbadc * Angular::neg1pow(k) *
                             Angular::tildeCk_kk(k, Fa.kappa(), Fc.kappa()) *
                             Angular::tildeCk_kk(k, Fb.kappa(), Fd.kappa());
            const auto eps_q = std::abs((qk1 - qk2) / qk2);
            REQUIRE(eps_q < 1.0e-14);
            const auto qk3 = Fa * Coulomb::Qkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
            const auto eps_q3 = std::abs((qk3 - qk2) / qk2);
            REQUIRE(eps_q3 < 1.0e-14);
          }
        }
      }
    }
  }

  // Check other formula
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      for (const auto &Fc : orbs) {
        if (Fc < Fa)
          continue;
        for (const auto &Fd : orbs) {
          if (Fd < Fb)
            continue;
          auto [k01, ki1] = Coulomb::k_minmax_P(Fa, Fb, Fc, Fd);
          auto [k02, ki2] = Coulomb::k_minmax_P(Fa.kappa(), Fb, Fc, Fd);
          REQUIRE(k01 == k02);
          REQUIRE(ki1 == ki2);
          auto [k0, ki] = Coulomb::k_minmax_W(Fa, Fb, Fc, Fd);
          for (int k = k0; k <= ki; ++k) {
            // Pk

            const auto pk1 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k);
            const auto pk2 = Fa * Coulomb::Pkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
            const auto eps_p =
                pk2 == 0.0 ? std::abs(pk1 - pk2) : std::abs((pk1 - pk2) / pk2);
            REQUIRE(eps_p < 1.0e-10);

            // P^k_abcd = \sum_l [k] 6j * Q^l_abdc
            auto [l0, li] = Coulomb::k_minmax_Q(Fa, Fb, Fd, Fc);
            // use larger l range than required, more thorough test
            double pk3 = 0.0;
            for (int l = l0 - 1; l <= li + 1; l++) {
              pk3 += (2 * k + 1) * Coulomb::sixj(Fa, Fc, k, Fb, Fd, l) *
                     Coulomb::Qk_abcd(Fa, Fb, Fd, Fc, l);
            }
            const auto eps_p3 =
                pk2 == 0.0 ? std::abs(pk1 - pk3) : std::abs((pk1 - pk3) / pk3);
            REQUIRE(eps_p3 < 1.0e-10);

            // test W:
            const auto qk = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k); // may be 0
            const auto wk = Coulomb::Wk_abcd(Fa, Fb, Fc, Fd, k);
            const auto wk2 = Fa * Coulomb::Wkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
            const auto wk3 = qk + pk1;
            const auto eps_w = std::abs(wk) < 1.0e-8 ?
                                   std::abs(wk2 - wk) :
                                   std::abs((wk2 - wk) / wk);
            const auto eps_w3 = std::abs(wk) < 1.0e-8 ?
                                    std::abs(wk3 - wk) :
                                    std::abs((wk3 - wk) / wk);
            REQUIRE(eps_w < 1.0e-10);
            REQUIRE(eps_w3 < 1.0e-10);
          }
        }
      }
    }
  }

  // test"magic" 6j formulas
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      for (const auto &Fc : orbs) {
        if (Fc < Fa)
          continue;
        for (const auto &Fd : orbs) {
          auto [k0, ki] = Coulomb::k_minmax_W(Fa, Fb, Fd, Fc);
          auto [l0, li] = Coulomb::k_minmax_W(Fa, Fb, Fc, Fd);
          for (int k = k0; k <= ki; ++k) {
            for (int l = l0; l <= li; ++l) {
              const auto sj1 = Coulomb::sixj(Fa, Fb, k, Fc, Fd, l);
              const auto sj2 = Angular::sixj_2(Fa.twoj(), Fb.twoj(), 2 * k,
                                               Fc.twoj(), Fd.twoj(), 2 * l);
              REQUIRE(std::abs(sj1 - sj2) < 1.0e-10);

              const auto triad = Coulomb::sixjTriads(Fa, Fb, k, Fc, Fd, l);
              const auto triangle = Coulomb::triangle(Fa, Fb, k);
              if (!triad || !triangle) {
                REQUIRE(sj1 == 0.0);
              }
            }
          }
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("Coulomb: yk tables", "[Coulomb][yktable][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: yk tables\n";

  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> aorbs, borbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      aorbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
      borbs.push_back(DiracSpinor::exactHlike(n_min + 1, l, radial_grid, zeff));
    }
    aorbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
    borbs.push_back(
        DiracSpinor::exactHlike(n_min + 1, -l - 1, radial_grid, zeff));
  }

  // test same orbs case:
  Coulomb::YkTable yaa(aorbs);
  for (const auto &Fa : aorbs) {
    for (const auto &Faa : aorbs) {
      if (Faa < Fa)
        continue;
      auto [k0, ki] = Angular::kminmax_Ck(Fa.kappa(), Faa.kappa());
      for (int k = k0; k <= ki; k += 2) {
        auto yk = yaa.get(k, Fa, Faa);
        REQUIRE(yk != nullptr);
        auto yk2 = Coulomb::yk_ab(Fa, Faa, k);
        const auto del = std::abs(qip::compare(*yk, yk2).first);
        REQUIRE(del < 1.0e-14);
        // symmetric: should return same pointer:
        auto yk3 = Coulomb::yk_ab(Faa, Fa, k);
        REQUIRE(yk3 == yk2);
      }
    }
  }

  // test different orbs case:
  Coulomb::YkTable yab(aorbs, borbs);
  for (const auto &Fa : aorbs) {
    for (const auto &Fb : borbs) {
      auto [k0, ki] = Angular::kminmax_Ck(Fa.kappa(), Fb.kappa());
      for (int k = k0; k <= ki; k += 2) {
        auto yk = yab.get(k, Fa, Fb);
        REQUIRE(yk != nullptr);
        auto yk2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto del = std::abs(qip::compare(*yk, yk2).first);
        REQUIRE(del < 1.0e-14);
        // symmetric: should return same pointer:
        auto yk3 = Coulomb::yk_ab(Fb, Fa, k);
        REQUIRE(yk3 == yk2);
      }
    }
  }

  // Test the Q,P,W formulas
  for (const auto &Fa : aorbs) {
    for (const auto &Fb : aorbs) {
      for (const auto &Fc : aorbs) {
        if (Fc < Fa)
          continue;
        for (const auto &Fd : aorbs) {
          if (Fd < Fb)
            continue;

          // test Q formulas from YkTable
          {
            auto [k0, ki] = Coulomb::k_minmax_Q(Fa, Fb, Fc, Fd);
            for (int k = k0; k <= ki; k += 2) {
              // test q
              auto q1 = yaa.Q(k, Fa, Fb, Fc, Fd);
              auto q2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
              auto eps = std::abs((q1 - q2) / q2);
              REQUIRE(eps < 1.0e-12);

              auto q3 = Fa * Coulomb::Qkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
              auto eps3 = std::abs((q3 - q2) / q2);
              REQUIRE(eps3 < 1.0e-10);
            }
            // test of k_minmax_Q
            auto q3 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k0 - 1);
            auto q4 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, ki + 1);
            REQUIRE(q3 == 0.0);
            REQUIRE(q4 == 0.0);
          }

          // test W,P formulas from YkTable
          {
            auto [k0, ki] = Coulomb::k_minmax_W(Fa, Fb, Fc, Fd);
            for (int k = k0; k <= ki; k++) {
              // test P
              auto p1 = yaa.P(k, Fa, Fb, Fc, Fd);
              auto p2 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k);
              auto eps =
                  p2 == 0.0 ? std::abs(p1 - p2) : std::abs((p1 - p2) / p2);
              REQUIRE(eps < 1.0e-10);

              auto p3 = Fa * Coulomb::Pkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
              auto eps3 =
                  p2 == 0.0 ? std::abs(p3 - p2) : std::abs((p3 - p2) / p2);
              REQUIRE(eps3 < 1.0e-10);
              // test W
              // nb: sometimes W "non-zero" but extremely small. OK?
              auto q1 = yaa.Q(k, Fa, Fb, Fc, Fd);
              auto w1 = yaa.W(k, Fa, Fb, Fc, Fd);
              auto w2 = q1 + p2;
              auto epsw = std::abs(w2) < 1.0e-8 ? std::abs(w1 - w2) :
                                                  std::abs((w1 - w2) / w2);
              REQUIRE(epsw < 1.0e-10);
            }
            // test of k_minmax_P,W
            auto [k0p, kip] = Coulomb::k_minmax_P(Fa, Fb, Fc, Fd);
            auto p3 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k0p - 1);
            auto p4 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, kip + 1);
            REQUIRE(p3 == 0.0);
            REQUIRE(p4 == 0.0);
            auto w3 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k0 - 1);
            auto w4 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, ki + 1);

            REQUIRE(w3 == 0.0);
            REQUIRE(w4 == 0.0);
          }
          //
        }
      }
    }
  }

  yab.extend_angular(15);
  REQUIRE(yab.Ck().max_tj() == 15);
  REQUIRE(yab.SixJ().max_2jk() == 15);

  //
}

//==============================================================================
//==============================================================================
//! Unit tests for Coulomb integrals (y^k_ab, R^k_abcd, lookup tables etc).
//! Also: tests some 6J table things
TEST_CASE("Coulomb: formulas", "[Coulomb][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: formulas\n";

  // Test the Coulomb formulas
  // Don't need dense grid, and use a local potential (Hartree)
  Wavefunction wf({900, 1.0e-6, 100.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Hartree", 0.0, "[Ne]");
  wf.formBasis({"4spd6fg8h8i", 30, 7, 1.0e-3, 1.0e-3, 40.0, false});

  // Split basis into core/excited
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fb : wf.basis()) {
    if (wf.isInCore(Fb.n(), Fb.kappa())) {
      core.push_back(Fb);
    } else {
      excited.push_back(Fb);
    }
  }
  // Form the Coulomb lookup tables:
  const Coulomb::YkTable Yce(core, excited);
  const Coulomb::YkTable Yij(wf.basis());

  { // Check the Ykab lookup-tables
    // check the 'different' orbitals case:
    double del1 = UnitTest::check_ykab_Tab(core, excited, Yce);
    // check the 'same' orbitals case:
    double del2 = UnitTest::check_ykab_Tab(wf.basis(), wf.basis(), Yij);
    REQUIRE(std::abs(del1) < 1.0e-17);
    REQUIRE(std::abs(del2) < 1.0e-17);
  }

  { // Testing the Hartree Y functions formula against naiive:
    const auto delk_core = UnitTest::check_ykab(wf.core(), 1);
    const auto delk_basis = UnitTest::check_ykab(wf.basis(), 1);
    for (const auto &dk : delk_core) {
      REQUIRE(std::abs(dk) < 1.0e-14);
    }
    for (const auto &dk : delk_basis) {
      REQUIRE(std::abs(dk) < 1.0e-14);
    }
  }

  // test R^k_abcd:
  const double eps_R = UnitTest::check_Rkabcd(wf.core(), 1);
  const double eps_R2 = UnitTest::check_Rkabcd(wf.basis(), 1);
  REQUIRE(std::abs(eps_R) < 1.0e-13);
  REQUIRE(std::abs(eps_R2) < 1.0e-13);

  //============================================================================
  // Test "other" Coulomb integrals: P, Q, W (these are defined in terms of
  // sums over R and angular coeficients):
  {
    // Contruct a vector of DiracSpinors, with just single spinor of each
    // kappa:
    std::vector<DiracSpinor> torbs;
    for (int kappa_index = 0;; ++kappa_index) {
      auto k = Angular::kappaFromIndex(kappa_index);
      auto phi = std::find_if(cbegin(wf.basis()), cend(wf.basis()),
                              [k](auto x) { return x.kappa() == k; });
      if (phi == cend(wf.basis()))
        break;
      torbs.emplace_back(*phi);
    }

    // test Q
    const auto maxtj = std::max_element(wf.basis().cbegin(), wf.basis().cend(),
                                        DiracSpinor::comp_j)
                           ->twoj();
    const auto &Ck = Yij.Ck();
    const Angular::SixJTable sj(2 * maxtj);

    double worstQ = 0.0;
    double worstP = 0.0;
    double worstW = 0.0;
    for (const auto &Fa : torbs) {
      for (const auto &Fb : torbs) {
        for (const auto &Fc : torbs) {
          for (const auto &Fd : torbs) {
            for (int k = 0; k <= Ck.max_k(); ++k) {

              double Q1 = 0.0;
              if (Angular::Ck_kk_SR(k, Fa.kappa(), Fc.kappa()) &&
                  Angular::Ck_kk_SR(k, Fb.kappa(), Fd.kappa())) {

                Q1 = Yij.Q(k, Fa, Fb, Fc, Fd);
                const auto Q2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
                const auto Q3 = Fa * Yij.Qkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);

                const auto Q4 =
                    Angular::neg1pow_2(2 * k + Fa.twoj() + Fb.twoj() + 2) *
                    Ck(k, Fa.kappa(), Fc.kappa()) *
                    Ck(k, Fb.kappa(), Fd.kappa()) *
                    Coulomb::Rk_abcd(Fa, Fb, Fc, Fd, k);

                // test the 'Qk' version, including k_minmax_Q
                const auto [kmin, kmax] = Coulomb::k_minmax_Q(Fa, Fb, Fc, Fd);
                // This k_min should have correct parity rule too
                const auto Q5 =
                    (k >= kmin && k <= kmax && (kmin % 2 == k % 2)) ?
                        Yij.Q(k, Fa, Fb, Fc, Fd) :
                        0.0;

                const auto delQ =
                    std::abs(qip::max_difference(Q1, Q2, Q3, Q4, Q5));
                if (delQ > worstQ)
                  worstQ = delQ;
              }

              // Calc P
              const auto P1 = Yij.P(k, Fa, Fb, Fc, Fd);
              const auto P2 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k);
              double P4 = 0.0;
              for (int l = 0; l <= Ck.max_k(); ++l) {
                if (Angular::Ck_kk_SR(l, Fa.kappa(), Fd.kappa()) &&
                    Angular::Ck_kk_SR(l, Fb.kappa(), Fc.kappa())) {
                  P4 += (2 * k + 1) * sj.get(Fa, Fc, k, Fb, Fd, l) *
                        Yij.Q(l, Fa, Fb, Fd, Fc);
                }
              }

              // test the 'Pk' version, including k_minmax_P
              const auto [kminP, kmaxP] = Coulomb::k_minmax_P(Fa, Fb, Fc, Fd);
              // This k_min CANNOT contain correct parity rule
              const auto P5 =
                  (k >= kminP && k <= kmaxP) ? Yij.P(k, Fa, Fb, Fc, Fd) : 0.0;

              const auto delP = std::abs(qip::max_difference(P1, P2, P4, P5));
              if (delP > worstP)
                worstP = delP;

              // calc W
              const auto W1 = Q1 + P1;
              const auto W2 = Coulomb::Wk_abcd(Fa, Fb, Fc, Fd, k);

              // test the 'Wk' version, including k_minmax_P
              const auto [kmin, kmax] = Coulomb::k_minmax_W(Fa, Fb, Fc, Fd);
              // This k_min CANNOT contain correct parity rule
              const auto W5 =
                  (k >= kmin && k <= kmax) ? Yij.W(k, Fa, Fb, Fc, Fd) : 0.0;

              const auto delW = std::abs(qip::max_difference(W1, W2, W5));
              if (delW > worstW)
                worstW = delW;

            } // k
          }   // d
        }     // c
      }       // b
    }         // a
    REQUIRE(std::abs(worstQ) < 5.0e-13);
    REQUIRE(std::abs(worstP) < 5.0e-14);
    REQUIRE(std::abs(worstW) < 5.0e-14);

    // Test the "Magic" 6J functions
    {
      double worst = 0.0;
      for (const auto &Fa : torbs) {
        for (const auto &Fb : torbs) {
          for (const auto &Fc : torbs) {
            for (const auto &Fd : torbs) {
              for (int k = 0; k < 15; ++k) {
                for (int l = 0; l < 15; ++l) {
                  auto sj1 = Angular::sixj_2(Fa.twoj(), Fb.twoj(), Fc.twoj(),
                                             Fd.twoj(), 2 * k, 2 * l);
                  auto sj2 = Coulomb::sixjTriads(Fa, Fb, Fc, Fd, k, l) ?
                                 Coulomb::sixj(Fa, Fb, Fc, Fd, k, l) :
                                 0.0;
                  if (std::abs(sj1 - sj2) > worst)
                    worst = std::abs(sj1 - sj2);
                }
              }
            }
          }
        }
      }
      REQUIRE(std::abs(worst) < 1.0e-13);
    }
  }
}

//============================================================================
//============================================================================

//============================================================================
inline std::vector<double> UnitTest::yk_naive(const DiracSpinor &Fa,
                                              const DiracSpinor &Fb, int k) {
  const auto &gr = Fa.grid();
  std::vector<double> yk(gr.r().size());
#pragma omp parallel for
  for (auto i = 0ul; i < yk.size(); ++i) {
    auto r = gr.r(i);

    auto rtkr = [](double x, double y, int kk) {
      return x < y ? std::pow(x / y, kk) / y : std::pow(y / x, kk) / x;
    };

    std::vector<double> f;
    f.reserve(gr.r().size());
    for (auto j = 0ul; j < yk.size(); ++j) {
      f.push_back(rtkr(r, gr.r(j), k) *
                  (Fa.f(j) * Fb.f(j) + Fa.g(j) * Fb.g(j)));
    }

    // const auto p0 = std::max(Fa.min_pt(), Fb.min_pt());
    // const auto pi = std::min(Fa.max_pt(), Fb.max_pt());
    yk[i] = NumCalc::integrate(gr.du(), 0, 0, f, gr.drdu());
  }

  return yk;
}

//============================================================================
inline double UnitTest::check_ykab_Tab(const std::vector<DiracSpinor> &a,
                                       const std::vector<DiracSpinor> &b,
                                       const Coulomb::YkTable &Yab) {
  //
  double worst = 0.0;
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        // Only check if Angular factor is non-zero (since Ykab only calc'd
        // in this case)
        if (!Angular::Ck_kk_SR(k, Fa.kappa(), Fb.kappa()))
          continue;
        const auto y1 = Yab.get(k, Fa, Fb);
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y3 = Coulomb::yk_ab(Fb, Fa, k); // symmetric
        if (y1 == nullptr) {
          std::cout << k << " " << Fa.symbol() << " " << Fb.symbol() << "\n";
        }
        const auto del = std::abs(qip::compare(*y1, y2).first) +
                         std::abs(qip::compare(y2, y3).first);
        if (del > worst)
          worst = del;
      }
    }
  }

  return worst;
}

//============================================================================
inline std::vector<double>
UnitTest::check_ykab(const std::vector<DiracSpinor> &orbs, int max_del_n) {

  // Compared Yk as calculated by fast Coulomb::yk_ab routine (used in the
  // code) to UnitTest::yk_naive, a very slow, but simple version. In
  // theory, should be exactly the same

  std::vector<double> worst; // = 0.0;
  // nb: slow, so only check sub-set
  for (auto ia = 0ul; ia < orbs.size(); ++ia) {
    const auto &Fa = orbs[ia];
    for (auto ib = ia; ib < orbs.size(); ++ib) {
      const auto &Fb = orbs[ib];
      const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        if (!Angular::Ck_kk_SR(k, Fa.kappa(), Fb.kappa()))
          continue;
        if (std::abs(Fa.n() - Fb.n()) > max_del_n)
          continue;
        if (std::size_t(k + 1) > worst.size())
          worst.resize(std::size_t(k + 1));
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y4 = UnitTest::yk_naive(Fa, Fb, k); // slow

        const auto max = *std::max_element(cbegin(y4), cend(y4), qip::comp_abs);
        const auto del = std::abs(qip::compare(y2, y4).first / max);

        if (del > worst[std::size_t(k)])
          worst[std::size_t(k)] = del;
      }
    }
  }

  return worst;
}

//============================================================================
double UnitTest::check_Rkabcd(const std::vector<DiracSpinor> &orbs,
                              int max_del_n) {
  double eps_R = 0.0;
  const Coulomb::YkTable Yab(orbs);
#pragma omp parallel for
  for (auto ia = 0ul; ia < orbs.size(); ia++) {
    const auto &Fa = orbs[ia];
    for (auto ib = 0ul; ib < orbs.size(); ib += 2) {
      const auto &Fb = orbs[ib];
      for (auto ic = ia; ic < orbs.size(); ic++) {
        const auto &Fc = orbs[ic];
        if (std::abs(Fa.n() - Fc.n()) > max_del_n)
          continue;
        if (std::abs(Fb.n() - Fc.n()) > max_del_n)
          continue;
        for (auto id = ib; id < orbs.size(); id += 2) {
          const auto &Fd = orbs[id];
          if (std::abs(Fb.n() - Fd.n()) > max_del_n)
            continue;
          const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fc);
          for (int k = kmin; k <= kmax; ++k) {
            if (!Angular::Ck_kk_SR(k, Fa.kappa(), Fc.kappa()))
              continue;
            if (!Angular::Ck_kk_SR(k, Fb.kappa(), Fd.kappa()))
              continue;

            //---------
            const auto yac = Yab.get(k, Fa, Fc);
            const auto ybd = Yab.get(k, Fb, Fd);

            const auto r1a = Coulomb::Rk_abcd(Fa, Fb, Fc, Fd, k);
            const auto r1b = Coulomb::Rk_abcd(Fb, Fa, Fd, Fc, k);
            const auto r1c = Coulomb::Rk_abcd(Fc, Fd, Fa, Fb, k);
            const auto r2a = Coulomb::Rk_abcd(Fa, Fc, *ybd);
            const auto r2b = Coulomb::Rk_abcd(Fb, Fd, *yac);
            const auto r2c = Coulomb::Rk_abcd(Fc, Fa, *ybd);
            const auto r3 = Fa * Coulomb::Rkv_bcd(Fa.kappa(), Fb, Fc, Fd, k);
            const auto r4 = Fa * Coulomb::Rkv_bcd(Fa.kappa(), Fc, *ybd);
            const auto eps = std::max({r1a, r1b, r1c, r2a, r2b, r2c, r3, r4}) -
                             std::min({r1a, r1b, r1c, r2a, r2b, r2c, r3, r4});
#pragma omp critical(compare_epsR)
            if (eps > eps_R) {
              eps_R = eps;
            }
          }
        }
      }
    }
  }
  return eps_R;
}
