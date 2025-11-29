#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_version.h>
#include <string>
#include <utility>

// This should make code work with old versions of GSL
// This is for a work-around for what appears to be a bug in GSL v1:
#ifdef GSL_MAJOR_VERSION
#if GSL_MAJOR_VERSION == 1
#define GSL_VERSION_1
#endif
#endif

namespace UnitTest {
// Helper functions for unit tests:

// This loops through every possibly C^k and three-j symbol (up to given maximum
// j), and calculates {C^k, 3j} symbol multiple ways (using lookup table, +
// direct calculation from basic formulas etc.). Checks all these against each
// other, and returns the eps of the worst comparison. Ck must be pre-filled,
// ck_m may be empty initially (_mutable versions are called)
double ck_compare_direct(const Angular::CkTable &Ck, Angular::CkTable &Ck_m);

// Loops through all possible 6j tables (including non-physical) and compares to
// direct calculation. Returns maximum difference (absolute value)
double sj_compare_direct(const Angular::SixJTable &sjt);

// As above, but uses DiracSpinor interface to sjt
double sj_compare_DiracSpinor(const Angular::SixJTable &sjt,
                              const std::vector<DiracSpinor> &basis);

// Creates a dummy set of DiracSpinors
std::vector<DiracSpinor> dummy_basis(int l_max);

// Calculates typical speedup of table over direct calculaion
double speedup_6jt(const Angular::SixJTable &sjt,
                   const std::vector<DiracSpinor> &basis);

} // namespace UnitTest

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
TEST_CASE("Angular: Winger369j functions", "[Angular][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Angular: Winger369j functions\n";

  // list in form: {kappa_index, kappa, l, 2j}
  const std::vector<std::tuple<int, int, int, int>> test_data{
      {0, -1, 0, 1},   {1, 1, 1, 1},   {2, -2, 1, 3},   {3, 2, 2, 3},
      {4, -3, 2, 5},   {5, 3, 3, 5},   {6, -4, 3, 7},   {7, 4, 4, 7},
      {8, -5, 4, 9},   {9, 5, 5, 9},   {10, -6, 5, 11}, {11, 6, 6, 11},
      {12, -7, 6, 13}, {13, 7, 7, 13}, {14, -8, 7, 15}, {15, 8, 8, 15},
      {16, -9, 8, 17}, {17, 9, 9, 17}};

  for (const auto &[ki, k, l, tj] : test_data) {
    REQUIRE(Angular::l_k(k) == l);
    REQUIRE(Angular::twoj_k(k) == tj);
    const auto parity = int(std::round(std::pow(-1, l)));
    REQUIRE(Angular::parity_k(k) == parity);
    REQUIRE(Angular::parity_l(l) == Angular::parity_k(k));
    REQUIRE(Angular::l_tilde_k(k) == ((k > 0) ? l - 1 : l + 1));
    REQUIRE(Angular::kappa_twojl(tj, l) == k);
    REQUIRE(Angular::kappa_twojpi(tj, parity) == k);
    REQUIRE(Angular::zeroQ(1.0e-11));
    REQUIRE(Angular::zeroQ(-1.0e-11));
    REQUIRE(!Angular::zeroQ(1.0e-9));
    REQUIRE(!Angular::zeroQ(-1.0e-9));

    REQUIRE(Angular::indexFromKappa(k) == ki);
    REQUIRE(Angular::kappaFromIndex(ki) == k);
    REQUIRE(Angular::twojFromIndex(ki) == tj);
    REQUIRE(Angular::lFromIndex(ki) == l);
  }

  //===========================================
  REQUIRE(Angular::states_below_n(1) == 0);
  REQUIRE(Angular::states_below_n(2) == 1);
  REQUIRE(Angular::states_below_n(3) == 4);
  REQUIRE(Angular::states_below_n(4) == 9);
  REQUIRE(Angular::states_below_n(5) == 16);

  int index = 0;
  for (int n = 1; n < 15; ++n) {
    for (int l = 0; l < n; ++l) {
      const int k1 = l;
      const int k2 = -l - 1;
      if (l > 0) {
        REQUIRE(Angular::nk_to_index(n, k1) == index);
        const auto [nx, kx] = Angular::index_to_nk(index);
        REQUIRE(nx == n);
        REQUIRE(kx == k1);
        REQUIRE(Angular::nkindex_to_kappa(index) == kx);
        REQUIRE(Angular::nkindex_to_twoj(index) == Angular::twoj_k(kx));
        REQUIRE(Angular::nkindex_to_l(index) == Angular::l_k(kx));
        ++index;
      }
      REQUIRE(Angular::nk_to_index(n, k2) == index);
      const auto [ny, ky] = Angular::index_to_nk(index);
      REQUIRE(ny == n);
      REQUIRE(ky == k2);
      REQUIRE(Angular::nkindex_to_kappa(index) == ky);
      REQUIRE(Angular::nkindex_to_twoj(index) == Angular::twoj_k(ky));
      REQUIRE(Angular::nkindex_to_l(index) == Angular::l_k(ky));
      ++index;
    }
  }
  //===========================================

  REQUIRE(Angular::evenQ(0));
  REQUIRE(Angular::evenQ(-2));
  REQUIRE(Angular::evenQ(2));
  REQUIRE(!Angular::evenQ(-1));
  REQUIRE(!Angular::evenQ(1));

  REQUIRE(Angular::evenQ_2(2 * 0));
  REQUIRE(Angular::evenQ_2(-2 * 2));
  REQUIRE(Angular::evenQ_2(2 * 2));
  REQUIRE(!Angular::evenQ_2(-2 * 1));
  REQUIRE(!Angular::evenQ_2(2 * 1));

  REQUIRE(Angular::neg1pow(0) == 1);
  REQUIRE(Angular::neg1pow(1) == -1);
  REQUIRE(Angular::neg1pow(2) == 1);
  REQUIRE(Angular::neg1pow(-1) == -1);
  REQUIRE(Angular::neg1pow(-2) == 1);

  REQUIRE(Angular::neg1pow_2(2 * 0) == 1);
  REQUIRE(Angular::neg1pow_2(2 * 1) == -1);
  REQUIRE(Angular::neg1pow_2(2 * 2) == 1);
  REQUIRE(Angular::neg1pow_2(-2 * 1) == -1);
  REQUIRE(Angular::neg1pow_2(-2 * 2) == 1);

  REQUIRE(Angular::parity(0, 0, 0) == 1);
  REQUIRE(Angular::parity(0, 0, 1) == 0);
  REQUIRE(Angular::parity(1, 1, 0) == 1);
  REQUIRE(Angular::parity(1, 1, 1) == 0);
  REQUIRE(Angular::parity(1, 0, 1) == 1);
  REQUIRE(Angular::parity(0, 1, 1) == 1);

  int max_l = 6;
  for (int la = 0; la <= max_l; ++la) {
    for (int lb = 0; lb <= max_l; ++lb) {
      auto l_min = std::abs(la - lb);
      auto l_max = (la + lb);
      REQUIRE(Angular::triangle(la, lb, l_min) == 1);
      REQUIRE(Angular::triangle(la, lb, l_max) == 1);
      REQUIRE(Angular::triangle(la, lb, l_min - 1) == 0);
      REQUIRE(Angular::triangle(la, lb, l_max + 1) == 0);
      const auto m_q = -la - lb;
      REQUIRE(Angular::sumsToZero(la, lb, m_q) == 1);
      REQUIRE(Angular::sumsToZero(la, lb, m_q + 1) == 0);
      REQUIRE(Angular::sumsToZero(la, lb, m_q - 1) == 0);
    }
  }

  // ThreeJSymbol[{5/2, 1/2}, {1, 0}, {3/2, -1/2}] = -1/Sqrt[10]
  const auto eps = 1.0e-12;
  REQUIRE(std::abs(Angular::threej_2(5, 2, 3, 1, 0, -1) -
                   (-1 / std::sqrt(10.0))) < eps);
  REQUIRE(std::abs(Angular::threej_2(5, 3, 2, 1, -1, 0) -
                   (+1 / std::sqrt(10.0))) < eps);
  REQUIRE(std::abs(Angular::threej_2(5, 2, 9, 1, 0, -1) - (0.0)) < eps);
  REQUIRE(std::abs(Angular::special_threej_2(5, 3, 2) -
                   (-1 / std::sqrt(10.0))) < eps);
  REQUIRE(std::abs(Angular::special_threej_2(5, 9, 2) - (0.0)) < eps);

  // ClebschGordan[{5, 0}, {4, 0}, {1, 0}] = Sqrt[5/33]
  REQUIRE(std::abs(Angular::cg_2(10, 0, 8, 0, 2, 0) - (std::sqrt(5 / 33.0))) <
          eps);

  // <ka||C^k||kb>
  REQUIRE(std::abs(Angular::Ck_kk(0, -1, -1) - (std::sqrt(2))) < eps);
  REQUIRE(std::abs(Angular::Ck_kk(1, -1, -1) - (0.0)) < eps);
  REQUIRE(std::abs(Angular::Ck_kk(1, -1, 1) - (-std::sqrt(2.0 / 3))) < eps);
  REQUIRE(std::abs(Angular::Ck_kk(2, -2, 1) - (-2 / std::sqrt(5.0))) < eps);
  REQUIRE(std::abs(Angular::Ck_kk(2, 1, -2) - (+2 / std::sqrt(5.0))) < eps);

  // Test C^k, including selection rules, tilde version etc.
  const auto kap_list = {-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6};
  for (auto ka : kap_list) {
    for (auto kb : kap_list) {
      auto [kmin, kmax] = Angular::kminmax_Ck(ka, kb);
      REQUIRE(Angular::Ck_kk(kmin, ka, kb) != 0.0);
      REQUIRE(Angular::Ck_kk(kmax, ka, kb) != 0.0);
      REQUIRE(Angular::Ck_kk(kmin - 1, ka, kb) == 0.0);
      REQUIRE(Angular::Ck_kk(kmax + 1, ka, kb) == 0.0);
      // only every second should be non-zero
      for (auto k = kmin; k <= kmax; k += 2) {

        const auto ck1 = Angular::Ck_kk(k, ka, kb);
        const auto ck_next = Angular::Ck_kk(k + 1, ka, kb);
        REQUIRE(ck1 != 0.0);
        REQUIRE(ck_next == 0.0);
        REQUIRE(Angular::Ck_kk_SR(k, ka, kb) == true);
        REQUIRE(Angular::Ck_kk_SR(k + 1, ka, kb) == false);

        // test tilde version (which should be symmetric)
        const auto tja = Angular::twoj_k(ka);
        const auto ck_tilde1 = Angular::tildeCk_kk(k, ka, kb);
        const auto ck_tilde2 = Angular::tildeCk_kk(k, kb, ka);
        const auto ck2 = Angular::neg1pow_2(tja + 1) * ck_tilde1;
        REQUIRE(std::abs(ck_tilde1 - ck_tilde2) < eps);
        REQUIRE(std::abs(ck1 - ck2) < eps);

        // <ka||C^k||kb> = (-1)^(ja+1/2) * srt([ja][jb]) * 3js(ja jb k, -1/2 1/2
        // 0) * Pi
        const auto tjb = Angular::twoj_k(kb);
        const auto la = Angular::l_k(ka);
        const auto lb = Angular::l_k(kb);
        const auto ck3 = Angular::neg1pow_2(tja + 1) *
                         std::sqrt((tja + 1) * (tjb + 1)) *
                         Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0) *
                         Angular::parity(la, lb, k);
        REQUIRE(std::abs(ck1 - ck3) < eps);
      }
    }
  }
}

//------------------------------------------------------------------------------
//! Unit tests for angular functions/classes (threeJ symbols, lookup
//! tables etc)
TEST_CASE("Angular: Ck tables", "[Angular][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Angular: Ck tables\n";

  // Maximum value of 2*j (for initial run)
  const int max2j_1 = 5;
  // Form C^k lookup tables. One used statically, one dynamically re-sized
  // ("dynamic" one uses the _mutable lookup functions, that calclate the
  // angular factor if it doesn't exist already)
  Angular::CkTable Ck(max2j_1);
  Angular::CkTable Ck_dynamic(0);
  REQUIRE(Ck.max_tj() >= max2j_1);
  REQUIRE(Ck_dynamic.max_tj() >= 0);

  // SECTION("Compare Ck tables agaist direct calculation")
  {
    const auto max_eps = UnitTest::ck_compare_direct(Ck, Ck_dynamic);
    REQUIRE(max_eps < 1.0e-13);
  }

  // SECTION("Compare Ck tables agaist direct calculation - use fill() to "
  // "extend")
  {
    // Test the "extend" capability (extend tables to larger j values:)
    const int max2j_2 = 15;
    // "Extend" the tables:
    Ck.fill(max2j_2);
    INFO("testing info")
    REQUIRE(Ck.max_tj() >= max2j_2);
    const auto max_eps = UnitTest::ck_compare_direct(Ck, Ck_dynamic);
    REQUIRE(max_eps < 1.0e-13);

    // extend _again_
    const int max2j_3 = 17;
    Ck.fill(max2j_3);
    REQUIRE(Ck.max_tj() >= max2j_3);
    const auto max_eps3 = UnitTest::ck_compare_direct(Ck, Ck_dynamic);
    REQUIRE(max_eps3 < 1.0e-13);
  }
}

//------------------------------------------------------------------------------
//! Unit tests for 6J symbol lookup tables
TEST_CASE("Angular: 6j tables", "[Angular][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Angular: 6j tables\n";

  // 6j tables - all symbols. This is the slow paer
  {
    Angular::SixJTable sjt;
    for (const auto max_k : {3, 9}) {
      sjt.fill(max_k); // nb: each loop, "extends" the table
      REQUIRE(sjt.max_2jk() == max_k);
      const auto delta = UnitTest::sj_compare_direct(sjt);
      REQUIRE(delta < 1.0e-14);
    }
  }

  const auto l_max = 6;
  std::vector<DiracSpinor> basis = UnitTest::dummy_basis(l_max);
  // Fill 6J table as usually done:
  const auto max_k = DiracSpinor::max_tj(basis); // max_k = 2*max_j
  Angular::SixJTable sjt(max_k);

  // SECTION("6j tables - with DiracSpinor basis")
  {
    const auto delta = UnitTest::sj_compare_DiracSpinor(sjt, basis);
    REQUIRE(delta < 1.0e-14);
  }
}

//------------------------------------------------------------------------------
TEST_CASE("Angular: 6j tables - performance",
          "[Angular][!mayfail][performance]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Angular: 6j tables - performance\n";

  const auto l_max = 6;
  std::vector<DiracSpinor> basis = UnitTest::dummy_basis(l_max);
  // Fill 6J table as usually done:
  const auto max_k = DiracSpinor::max_tj(basis); // max_k = 2*max_j
  Angular::SixJTable sjt(max_k);
  const auto speedup = UnitTest::speedup_6jt(sjt, basis);
  REQUIRE(speedup > 2.0);
}

//==============================================================================
//==============================================================================
//==============================================================================
double UnitTest::ck_compare_direct(const Angular::CkTable &Ck,
                                   Angular::CkTable &Ck_m) {
  const auto max2j = Ck.max_tj();
  double max = 0.0;

  for (int kia = 0;; ++kia) {
    const auto tja = Angular::twojFromIndex(kia);
    if (tja > max2j)
      break;
    const auto ka = Angular::kappaFromIndex(kia);
    for (int kib = 0;; ++kib) {
      const auto tjb = Angular::twojFromIndex(kib);
      if (tjb > max2j)
        break;
      const auto kb = Angular::kappaFromIndex(kib);

      // loop through all k multipolarities:
      for (int k = 0; k <= max2j; ++k) {

        const auto c1 = Ck.get_Ckab(k, ka, kb);
        const auto c1_m = Ck_m.get_Ckab_mutable(k, ka, kb);
        const auto c2 = Angular::Ck_kk(k, ka, kb);
        const auto sign = ((tja + 1) / 2) % 2 == 0 ? 1 : -1;
        const auto c3 = Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k) *
                        Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0) * sign *
                        std::sqrt((tja + 1) * (tjb + 1));

        const auto tc1 = Ck.get_tildeCkab(k, ka, kb);
        const auto tc2 =
            Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k) *
            Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0) *
            std::sqrt((tja + 1) * (tjb + 1));

        const auto tj1 = Ck.get_3jkab(k, ka, kb);
        const auto tj2 = Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        const auto tj3 = gsl_sf_coupling_3j(tja, tjb, 2 * k, -1, 1, 0);

        const auto l1 = Ck.get_Lambdakab(k, ka, kb);
        const auto l2 =
            tj1 * tj1 * Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k);

        max = qip::max_abs(max, c1 - c2, c1_m - c1, c2 - c3, tc1 - tc2,
                           tj1 - tj2, tj2 - tj3, l1 - l2);
      }
    }
  }

  return max;
}

//==============================================================================
double UnitTest::sj_compare_direct(const Angular::SixJTable &sjt) {
  auto max_del = 0.0;
  const auto max_2k = 2 * sjt.max_2jk();
  // compare values from table to direct calculation, no short-circuit
  for (int a = 0; a <= max_2k; ++a) {
    for (int b = 0; b <= max_2k; ++b) {
      for (int c = 0; c <= max_2k; ++c) {
        for (int d = 0; d <= max_2k; ++d) {
          for (int e = 0; e <= max_2k; ++e) {
            for (int f = 0; f <= max_2k; ++f) {

#ifndef GSL_VERSION_1
              const auto sj = gsl_sf_coupling_6j(a, b, c, d, e, f);
#else
              // nb: there seems to be a bug in GSL:V1 where
              // gsl_sf_coupling_6j returns non-zero result when
              // symbol is non-triangular
              const auto sj = Angular::sixj_2(a, b, c, d, e, f);
#endif
              const auto sj2 = sjt.get_2(a, b, c, d, e, f);
              const auto sj3 = Angular::sixj_2(a, b, c, d, e, f);
              const auto del = std::max(std::abs(sj - sj2), std::abs(sj - sj3));
              if (del > max_del) {
                max_del = del;
              }
              if (Angular::sixj_zeroQ(a, b, c, d, e, f) &&
                  std::abs(sj) > 1.0e-12) {
                max_del = 1.0;
              }
              if (!Angular::sixjTriads(a, b, c, d, e, f) &&
                  std::abs(sj) > 1.0e-12) {
                max_del = 1.0;
              }
            }
          }
        }
      }
    }
  }
  return max_del;
}

//==============================================================================
std::vector<DiracSpinor> UnitTest::dummy_basis(int l_max) {
  std::vector<DiracSpinor> basis;
  for (auto l = 0; l <= l_max; ++l) {
    // 2j = 2l+1, 2l-1
    if (l != 0) {
      const auto kappa = Angular::kappa_twojl(2 * l - 1, l);
      basis.emplace_back(0, kappa, nullptr);
    }
    const auto kappa = Angular::kappa_twojl(2 * l + 1, l);
    basis.emplace_back(0, kappa, nullptr);
  }
  return basis;
}

//==============================================================================
double UnitTest::sj_compare_DiracSpinor(const Angular::SixJTable &sjt,
                                        const std::vector<DiracSpinor> &basis) {
  // Test DiracSpinor (same test as above, but with DiracSpinor)
  // Tests {a,b,k \ c, d, l} and {a,b,k, \ l,u,c} seperately
  // loop through all kinds of 6J symbols
  auto max_del = 0.0;
  for (const auto &a : basis) {
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        for (const auto &d : basis) {
          // add 5, ensure our table is not missing any non-zero
          const auto tjmax =
              std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()}) + 5;
          for (int k = 0; k <= tjmax; ++k) {
            for (int l = 0; l <= tjmax; ++l) {
              const auto sj1 = gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k,
                                                  c.twoj(), d.twoj(), 2 * l);
              const auto sj2 = sjt.get(a, b, k, c, d, l);

              const auto sj3 = Angular::SixJ(a, b, k, c, d, l);

              REQUIRE(sj3 == Approx(sj1));

              // Compare, find worst offender:
              const auto del = std::abs(sj1 - sj2);
              if (del > max_del) {
                max_del = del;
              }
            }
          }
        }
      }
    }
  }
  for (const auto &a : basis) {
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        const auto tjmax = 2 * std::max({a.twoj(), b.twoj(), c.twoj()});
        for (int k = 0; k <= tjmax; ++k) {
          for (int l = 0; l <= tjmax; ++l) {
            for (int u = 0; u <= tjmax; ++u) {
              const auto sj1 = gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k,
                                                  2 * l, 2 * u, c.twoj());
              const auto sj2 = sjt.get(a, b, k, l, u, c);

              // Compare, find worst offender:
              const auto del = std::abs(sj1 - sj2);
              if (del > max_del) {
                max_del = del;
              }
            }
          }
        }
      }
    }
  }
  return max_del;
}

//==============================================================================
double UnitTest::speedup_6jt(const Angular::SixJTable &sjt,
                             const std::vector<DiracSpinor> &basis) {

  // performance test
  std::cout << "6JTable performance test:\n";
  double t1 = 0.0, t2 = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  {
    IO::ChronoTimer t("New");
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            const auto tjmax =
                std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()});
            for (int k = 0; k <= tjmax; ++k) {
              for (int l = 0; l <= tjmax; ++l) {
                sum1 += sjt.get(a, b, k, c, d, l);
              }
            }
          }
          const auto tjmax = std::max({a.twoj(), b.twoj(), c.twoj()});
          for (int k = 0; k <= tjmax; ++k) {
            for (int l = 0; l <= tjmax; ++l) {
              for (int u = 0; u <= tjmax; ++u) {
                sum1 += sjt.get(a, b, k, l, u, c);
              }
            }
          }
        }
      }
    }
    t1 = t.lap_reading_ms();
  }
  {
    IO::ChronoTimer t("Old");
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            const auto tjmax =
                std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()});
            for (int k = 0; k <= tjmax; ++k) {
              for (int l = 0; l <= tjmax; ++l) {
                sum2 += gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k, c.twoj(),
                                           d.twoj(), 2 * l);
              }
            }
          }
          const auto tjmax = std::max({a.twoj(), b.twoj(), c.twoj()});
          for (int k = 0; k <= tjmax; ++k) {
            for (int l = 0; l <= tjmax; ++l) {
              for (int u = 0; u <= tjmax; ++u) {
                sum2 += gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k, 2 * l,
                                           2 * u, c.twoj());
              }
            }
          }
        }
      }
    }
    t2 = t.lap_reading_ms();
  }
  std::cout << sum1 << " " << sum2 << " " << 2.0 * (sum1 - sum2) / (sum1 + sum2)
            << "\n";
  std::cout << t2 / t1 << "x speedup\n";
  return t2 / t1;
}