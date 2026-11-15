#include "MBPT/LadderExternalLine.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/Ladder.hpp"
#include "MBPT/Ladder.testdata.hpp"
#include "MBPT/LadderIteration.hpp"
#include "MBPT/LadderReference.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace UnitTest::LadderTest {

//! One ladder integral L^k_abcd as a function of (k, a, b, c, d)
using LadderIntegral =
  std::function<double(int, const DiracSpinor &, const DiracSpinor &,
                       const DiracSpinor &, const DiracSpinor &)>;

//! Ladder table with i orbitals i_orbs: lowest order, then num_sweeps sweeps
inline Coulomb::LkTable dressed_table(const LadderTestData &data,
                                      const std::vector<DiracSpinor> &i_orbs,
                                      int num_sweeps) {
  const auto max_k = stored_max_k(data);
  Coulomb::LkTable lk;
  MBPT::add_Lk_mnib_entries(&lk, data.excited, data.core, i_orbs, max_k);
  MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, i_orbs, {},
                          nullptr, max_k, false, 0.0);
  for (int sweep = 0; sweep < num_sweeps; ++sweep) {
    const auto lk_prev = lk;
    MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, i_orbs, {},
                            &lk_prev, max_k, false, 0.0);
  }
  return lk;
}

//! The ladder integrals with the external line of v, every L_x (not only
//! x = v) taken from lk, which must hold them all
inline MBPT::ExternalLineLadder ladder_from_table(const DiracSpinor &v,
                                                  const LadderTestData &data,
                                                  const Coulomb::LkTable &lk) {
  auto ladder = MBPT::external_line_start(v, data.core, data.excited, lk);
  for (std::size_t k = 0; k < ladder.L_mnxb->size(); ++k) {
    for (auto &[key, value] : ladder.L_mnxb->at(k)) {
      value = lk.Q(int(k), key);
    }
  }
  for (std::size_t k = 0; k < ladder.L_xnab->size(); ++k) {
    for (auto &[key, value] : ladder.L_xnab->at(k)) {
      const auto [x, n, a, b] = ladder.L_xnab.UnFormIndex(key);
      value = lk.Q(int(k), x, n, a, b);
    }
  }
  return ladder;
}

//! One sweep of both tables, from start, with the external line at en_v
inline MBPT::ExternalLineLadder
external_line_sweep(const MBPT::ExternalLineLadder &start, double en_v,
                    const LadderTestData &data, const Coulomb::LkTable &lk,
                    const MBPT::LadderDiagrams &diagrams) {
  auto ladder = start;
  MBPT::calculate_Lk_mnxb(&ladder.L_mnxb, data.qk, data.excited, data.core,
                          ladder.x_orbitals, en_v, diagrams, lk, start.L_mnxb,
                          0.0);
  MBPT::calculate_Lk_xnab(&ladder.L_xnab, data.qk, data.excited, data.core,
                          ladder.x_orbitals, en_v, diagrams, lk, start.L_xnab,
                          0.0);
  return ladder;
}

//! Largest difference between L^k_mnxb of ladder (x in x_list) and
//! expected(k,m,n,x,b), relative to the largest expected; 1 if all are zero
inline double mnxb_difference(const MBPT::ExternalLineLadder &ladder,
                              const std::vector<DiracSpinor> &x_list,
                              const LadderTestData &data,
                              const LadderIntegral &expected) {
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &x : x_list) {
    for (const auto &m : data.excited) {
      for (const auto &n : data.excited) {
        for (const auto &b : data.core) {
          const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, x, b);
          for (int k = k_min; k <= std::min(k_max, ladder.L_mnxb.max_k());
               ++k) {
            const auto value = expected(k, m, n, x, b);
            max_difference = std::max(
              max_difference, std::abs(value - ladder.L_mnxb.Q(k, m, n, x, b)));
            largest = std::max(largest, std::abs(value));
          }
        }
      }
    }
  }
  return largest == 0.0 ? 1.0 : max_difference / largest;
}

//! Largest difference between L^k_xnab of ladder (x in x_list, n in n_list)
//! and expected(k,x,n,a,b), relative to the largest expected; 1 if all zero
inline double xnab_difference(const MBPT::ExternalLineLadder &ladder,
                              const std::vector<DiracSpinor> &x_list,
                              const std::vector<DiracSpinor> &n_list,
                              const LadderTestData &data,
                              const LadderIntegral &expected) {
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &x : x_list) {
    for (const auto &n : n_list) {
      for (const auto &a : data.core) {
        for (const auto &b : data.core) {
          const auto [k_min, k_max] = MBPT::k_minmax_L(x, n, a, b);
          for (int k = k_min; k <= std::min(k_max, ladder.L_xnab.max_k());
               ++k) {
            const auto value = expected(k, x, n, a, b);
            max_difference = std::max(
              max_difference, std::abs(value - ladder.L_xnab.Q(k, x, n, a, b)));
            largest = std::max(largest, std::abs(value));
          }
        }
      }
    }
  }
  return largest == 0.0 ? 1.0 : max_difference / largest;
}

} // namespace UnitTest::LadderTest

//==============================================================================
TEST_CASE("Ladder: external line at its own energy", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  // Every excited state as an i orbital: the table holds every L_x
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto lk = dressed_table(data, i_orbs, 1);
  auto lk_next = lk;
  MBPT::calculate_Lk_mnib(&lk_next, data.qk, data.excited, data.core, i_orbs,
                          {}, &lk, max_k, false, 0.0);

  // With en_v = e_x0, the x0 integrals are those of a normal sweep
  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  std::size_t num_lines = 0;
  for (const auto &v : data.valence) {
    const auto start = ladder_from_table(v, data, lk);
    num_lines += start.x_orbitals.size();
    for (const auto &x0 : start.x_orbitals) {
      const auto ladder = external_line_sweep(start, x0.en(), data, lk, {});
      const auto from_table = [&](int k, const DiracSpinor &a,
                                  const DiracSpinor &b, const DiracSpinor &c,
                                  const DiracSpinor &d) {
        return lk_next.Q(k, a, b, c, d);
      };
      worst_mnxb =
        std::max(worst_mnxb, mnxb_difference(ladder, {x0}, data, from_table));
      worst_xnab =
        std::max(worst_xnab,
                 xnab_difference(ladder, {x0}, data.excited, data, from_table));
    }
  }
  // Some valence state has more than one projection state
  REQUIRE(num_lines > data.valence.size());
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: external line matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto lk = dressed_table(data, i_orbs, 1);
  // The integrals with the external line start from a different table, so a
  // dressed integral read from the wrong table is seen
  const auto lk_line = dressed_table(data, i_orbs, 2);

  for (const auto &v : data.valence) {
    // An energy away from every e_x, and from every zero of the denominators
    const auto en_v = v.en() - 0.05;
    const auto start = ladder_from_table(v, data, lk_line);
    const auto &x_list = start.x_orbitals;

    // Each diagram alone: L^k_mnxb [x as i] and L^k_xnab [x as m]. The
    // dressed integrals that contain x come from lk_line: L1 and L3 of
    // L^k_mnxb, L2 and L4 of L^k_xnab
    const auto L1 =
      external_line_sweep(start, en_v, data, lk, {true, false, false});
    const auto L2_L3 =
      external_line_sweep(start, en_v, data, lk, {false, true, false});
    const auto L4 =
      external_line_sweep(start, en_v, data, lk, {false, false, true});

    const auto eps_L1_mnxb =
      mnxb_difference(L1, x_list, data,
                      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &x, const DiracSpinor &b) {
                        return MBPT::L1(k, m, n, x, b, data.qk, data.excited,
                                        sjt, &lk_line, en_v);
                      });
    const auto eps_L2_L3_mnxb =
      mnxb_difference(L2_L3, x_list, data,
                      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &x, const DiracSpinor &b) {
                        return MBPT::L2(k, m, n, x, b, data.qk, data.core,
                                        data.excited, sjt, &lk) +
                               MBPT::L3(k, m, n, x, b, data.qk, data.core,
                                        data.excited, sjt, &lk_line, en_v);
                      });
    const auto eps_L4_mnxb = mnxb_difference(
      L4, x_list, data,
      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &x, const DiracSpinor &b) {
        return MBPT::L4(k, m, n, x, b, data.qk, data.core, sjt, &lk);
      });

    const auto eps_L1_xnab = xnab_difference(
      L1, x_list, data.excited, data,
      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b) {
        return MBPT::L1(k, x, n, a, b, data.qk, data.excited, sjt, &lk);
      });
    const auto L2_L3_xnab = [&](int k, const DiracSpinor &x,
                                const DiracSpinor &n, const DiracSpinor &a,
                                const DiracSpinor &b) {
      return MBPT::L2(k, x, n, a, b, data.qk, data.core, data.excited, sjt,
                      &lk_line, {}, en_v) +
             MBPT::L3(k, x, n, a, b, data.qk, data.core, data.excited, sjt,
                      &lk);
    };
    const auto eps_L2_L3_xnab =
      xnab_difference(L2_L3, x_list, data.excited, data, L2_L3_xnab);
    // n also a projection state: n is an internal line in L3 and the line in
    // L2 of the same integral
    const auto eps_L2_L3_xnab_n_projection =
      xnab_difference(L2_L3, x_list, x_list, data, L2_L3_xnab);
    const auto eps_L4_xnab = xnab_difference(
      L4, x_list, data.excited, data,
      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b) {
        return MBPT::L4(k, x, n, a, b, data.qk, data.core, sjt, &lk_line, en_v);
      });

    REQUIRE(eps_L1_mnxb < 1.0e-12);
    REQUIRE(eps_L2_L3_mnxb < 1.0e-12);
    REQUIRE(eps_L4_mnxb < 1.0e-12);
    REQUIRE(eps_L1_xnab < 1.0e-12);
    REQUIRE(eps_L2_L3_xnab < 1.0e-12);
    REQUIRE(eps_L2_L3_xnab_n_projection < 1.0e-12);
    REQUIRE(eps_L4_xnab < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: external line potential", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data(-1, 8.0);
  const auto &sjt = data.yk.SixJ();
  // Converged table (the iteration converges geometrically)
  const auto lk = dressed_table(data, data.i_orbs, 40);

  fmt::print("\nLadder: <v|Sigma_L|v> vs de_valence_w (au); eps(<y|Sigma_L|x>) "
             "vs the sum over the tables\n");
  fmt::print("{:4s}  {:>13s}  {:>13s}  {:>8s}  {:>8s}  {:>8s}\n", "v",
             "expected", "found", "eps", "eps(L_v)", "eps(yx)");
  double worst_eps = 0.0;
  double worst_eps_L = 0.0;
  double worst_eps_yx = 0.0;
  for (const auto &v : data.valence) {
    const auto ladder = MBPT::external_line_ladder(
      v, data.core, data.excited, data.qk, lk, {}, 60, 1.0e-13, 0.0, false);

    // x = v: the converged integrals are those of the ladder table
    const auto eps_L = std::max(
      mnxb_difference(ladder, {v}, data,
                      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &x, const DiracSpinor &b) {
                        return lk.Q(k, m, n, x, b);
                      }),
      xnab_difference(ladder, {v}, data.excited, data,
                      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                          const DiracSpinor &a, const DiracSpinor &b) {
                        return lk.Q(k, x, n, a, b);
                      }));

    // The column Sigma_L|v> is exact: <v|Sigma_L|v> is the ladder energy
    // (lower components included in Sigma_L, as they are in the energy)
    const auto Sigma = MBPT::Sigma_ladder(v, ladder, data.core, data.excited,
                                          1.0e-4, 200.0, 1, true);
    const auto expected =
      MBPT::de_valence_w(v, data.qk, lk, data.core, data.excited, &sjt);
    const auto found = v * (Sigma * v);
    const auto eps = std::abs(found / expected - 1.0);

    // Every <y|Sigma_L|x>, x != v included:
    // sum_{amn,k} W^k_yamn L^k_mnxa / ([k][j_v](e_v+e_a-e_m-e_n))
    //   + sum_{nab,k} W^k_ynab L^k_xnab / ([k][j_v](e_v+e_n-e_a-e_b))
    double max_difference = 0.0;
    double largest = 0.0;
    for (const auto &x : ladder.x_orbitals) {
      const auto Sigma_x = Sigma * x;
      for (const auto &y : ladder.x_orbitals) {
        double sum = 0.0;
        for (const auto &n : data.excited) {
          for (const auto &a : data.core) {
            for (const auto &m : data.excited) {
              const auto de = v.en() + a.en() - m.en() - n.en();
              const auto [k_min, k_max] = Coulomb::k_minmax_W(y, a, m, n);
              for (int k = k_min; k <= k_max; ++k) {
                sum += data.qk.W(k, y, a, m, n, &sjt) *
                       ladder.L_mnxb.Q(k, m, n, x, a) / ((2 * k + 1) * de);
              }
            }
            for (const auto &b : data.core) {
              const auto de = v.en() + n.en() - a.en() - b.en();
              const auto [k_min, k_max] = Coulomb::k_minmax_W(y, n, a, b);
              for (int k = k_min; k <= k_max; ++k) {
                sum += data.qk.W(k, y, n, a, b, &sjt) *
                       ladder.L_xnab.Q(k, x, n, a, b) / ((2 * k + 1) * de);
              }
            }
          }
        }
        sum /= v.twojp1();
        max_difference = std::max(max_difference, std::abs(y * Sigma_x - sum));
        largest = std::max(largest, std::abs(sum));
      }
    }
    const auto eps_yx = max_difference / largest;

    fmt::print("{:4s}  {:+13.6e}  {:+13.6e}  {:8.1e}  {:8.1e}  {:8.1e}\n",
               v.shortSymbol(), expected, found, eps, eps_L, eps_yx);
    worst_eps = std::max(worst_eps, eps);
    worst_eps_L = std::max(worst_eps_L, eps_L);
    worst_eps_yx = std::max(worst_eps_yx, eps_yx);
  }
  REQUIRE(worst_eps_L < 1.0e-9);
  REQUIRE(worst_eps < 1.0e-8);
  REQUIRE(worst_eps_yx < 1.0e-8);
}
