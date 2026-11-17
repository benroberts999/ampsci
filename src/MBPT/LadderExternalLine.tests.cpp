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
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace UnitTest::LadderTest {

//! One ladder integral L^k_abcd as a function of (k, a, b, c, d)
using LadderIntegral =
  std::function<double(int, const DiracSpinor &, const DiracSpinor &,
                       const DiracSpinor &, const DiracSpinor &)>;

//! Ladder table with i orbitals i_orbs: lowest order, then num_sweeps sweeps
inline Coulomb::LkTable
dressed_table(const LadderTestData &data,
              const std::vector<DiracSpinor> &i_orbs, int num_sweeps,
              const MBPT::LadderDiagrams &diagrams = {}) {
  const auto max_k = stored_max_k(data);
  Coulomb::LkTable lk;
  MBPT::add_Lk_mnib_entries(&lk, data.excited, data.core, i_orbs, max_k);
  MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, i_orbs,
                          diagrams, nullptr, max_k, false, 0.0);
  for (int sweep = 0; sweep < num_sweeps; ++sweep) {
    const auto lk_prev = lk;
    MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, i_orbs,
                            diagrams, &lk_prev, max_k, false, 0.0);
  }
  return lk;
}

//! The ladder diagrams with the three pieces of P
inline MBPT::LadderDiagrams diagrams_with_P() {
  return {true, true, true, true, true, true, true};
}

//! The diagrams of x_F: piece (1) of P and the hole-particle term
inline MBPT::LadderDiagrams diagrams_of_x_F() {
  MBPT::LadderDiagrams diagrams{false, false, false, false};
  diagrams.P1 = true;
  diagrams.hole_particle = true;
  return diagrams;
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
                    const MBPT::LadderDiagrams &diagrams,
                    const LinAlg::Matrix<double> *rho_ma = nullptr) {
  auto ladder = start;
  MBPT::calculate_Lk_mnxb(&ladder.L_mnxb, data.qk, data.excited, data.core,
                          ladder.x_orbitals, en_v, diagrams, lk, start.L_mnxb,
                          0.0, rho_ma);
  MBPT::calculate_Lk_xnab(&ladder.L_xnab, data.qk, data.excited, data.core,
                          ladder.x_orbitals, en_v, diagrams, lk, start.L_xnab,
                          0.0, rho_ma);
  return ladder;
}

//! One sweep of L^k_xnab alone, from start, with the external line at en_v;
//! L^k_mnxb is that of start
inline MBPT::ExternalLineLadder
xnab_sweep(const MBPT::ExternalLineLadder &start, double en_v,
           const LadderTestData &data, const Coulomb::LkTable &lk,
           const MBPT::LadderDiagrams &diagrams) {
  auto ladder = start;
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

//! Largest difference between L^k_abcd of table and expected(k,a,b,c,d), for
//! a in a_list etc. and k <= max_k, relative to the largest expected; 1 if
//! all are zero
inline double table_vs_reference(const Coulomb::NkTable &table,
                                 const std::vector<DiracSpinor> &a_list,
                                 const std::vector<DiracSpinor> &b_list,
                                 const std::vector<DiracSpinor> &c_list,
                                 const std::vector<DiracSpinor> &d_list,
                                 int max_k, const LadderIntegral &expected) {
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &a : a_list) {
    for (const auto &b : b_list) {
      for (const auto &c : c_list) {
        for (const auto &d : d_list) {
          const auto [k_min, k_max] = MBPT::k_minmax_L(a, b, c, d);
          for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
            const auto value = expected(k, a, b, c, d);
            max_difference = std::max(max_difference,
                                      std::abs(value - table.Q(k, a, b, c, d)));
            largest = std::max(largest, std::abs(value));
          }
        }
      }
    }
  }
  return largest == 0.0 ? 1.0 : max_difference / largest;
}

//! Every orbital of the test data (core and excited) with the kappa of v
inline std::vector<DiracSpinor> orbitals_of_kappa(const DiracSpinor &v,
                                                  const LadderTestData &data) {
  return qip::select_if(
    qip::merge(data.core, data.excited),
    [&v](const DiracSpinor &x) { return x.kappa() == v.kappa(); });
}

//! The excited orbitals with n <= 4: a smaller set of internal lines, so that
//! the d and f orbitals with n = 5, 6 are outside it
inline std::vector<DiracSpinor> excited_up_to_n4(const LadderTestData &data) {
  return qip::select_if(data.excited,
                        [](const DiracSpinor &x) { return x.n() <= 4; });
}

//! <bra|Sigma_mirror|ket> at energy en, one integral at a time: the six
//! lowest-order terms with ket on the Coulomb vertex and bra in the ladder
//! integral (closed with Q for L2 and L3, with W for L1 and L4); internal
//! lines from excited and data.core
inline double mirror_terms_element(const DiracSpinor &bra,
                                   const DiracSpinor &ket, double en,
                                   const std::vector<DiracSpinor> &excited,
                                   const LadderTestData &data) {
  const auto &sjt = data.yk.SixJ();
  double sum = 0.0;
  for (const auto &n : excited) {
    for (const auto &a : data.core) {
      for (const auto &m : excited) {
        const auto de = en + a.en() - m.en() - n.en();
        const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, bra, a);
        for (int k = k_min; k <= k_max; ++k) {
          const auto f = (2 * k + 1) * ket.twojp1() * de;
          sum += data.qk.Q(k, ket, a, m, n) *
                 MBPT::L2(k, m, n, bra, a, data.qk, data.core, excited, sjt) /
                 f;
          sum += data.qk.W(k, ket, a, m, n, &sjt) *
                 MBPT::L4(k, m, n, bra, a, data.qk, data.core, sjt) / f;
        }
      }
      for (const auto &b : data.core) {
        const auto de = en + n.en() - a.en() - b.en();
        const auto [k_min, k_max] = MBPT::k_minmax_L(bra, n, a, b);
        for (int k = k_min; k <= k_max; ++k) {
          const auto f = (2 * k + 1) * ket.twojp1() * de;
          sum += data.qk.W(k, ket, n, a, b, &sjt) *
                 MBPT::L1(k, bra, n, a, b, data.qk, excited, sjt) / f;
          sum += data.qk.Q(k, ket, n, a, b) *
                 MBPT::L3(k, bra, n, a, b, data.qk, data.core, excited, sjt) /
                 f;
        }
      }
    }
  }
  return sum;
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

    // Each diagram alone (L2 and L3 together for L^k_mnxb): L^k_mnxb [x as
    // i] and L^k_xnab [x as m]. The dressed integrals that contain x come
    // from lk_line: L1 and L3 of L^k_mnxb, L2 and L4 of L^k_xnab
    const auto L1 =
      external_line_sweep(start, en_v, data, lk, {true, false, false, false});
    const auto L2_L3 =
      external_line_sweep(start, en_v, data, lk, {false, true, true, false});
    const auto L4 =
      external_line_sweep(start, en_v, data, lk, {false, false, false, true});
    // L2 and L3 alone: only L^k_xnab, which has no symmetry
    const auto L2 =
      xnab_sweep(start, en_v, data, lk, {false, true, false, false});
    const auto L3 =
      xnab_sweep(start, en_v, data, lk, {false, false, true, false});

    // L^k_mnxb: L2 and L3 together, since L_mnxb and L_nmbx are one entry
    // and L3^k_mnxb = L2^k_nmbx
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

    // L^k_xnab. For L2 and L3, n also a projection state: n is an internal
    // line in L3 and the line in L2 of the same integral
    const auto eps_L1_xnab = xnab_difference(
      L1, x_list, data.excited, data,
      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b) {
        return MBPT::L1(k, x, n, a, b, data.qk, data.excited, sjt, &lk);
      });
    const auto L2_xnab = [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                             const DiracSpinor &a, const DiracSpinor &b) {
      return MBPT::L2(k, x, n, a, b, data.qk, data.core, data.excited, sjt,
                      &lk_line, {}, en_v);
    };
    const auto L3_xnab = [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                             const DiracSpinor &a, const DiracSpinor &b) {
      return MBPT::L3(k, x, n, a, b, data.qk, data.core, data.excited, sjt,
                      &lk);
    };
    const auto eps_L2_xnab =
      std::max(xnab_difference(L2, x_list, data.excited, data, L2_xnab),
               xnab_difference(L2, x_list, x_list, data, L2_xnab));
    const auto eps_L3_xnab =
      std::max(xnab_difference(L3, x_list, data.excited, data, L3_xnab),
               xnab_difference(L3, x_list, x_list, data, L3_xnab));
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
    REQUIRE(eps_L2_xnab < 1.0e-12);
    REQUIRE(eps_L3_xnab < 1.0e-12);
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

//==============================================================================
TEST_CASE("Ladder: ratio2 potential", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data(-1, 8.0);
  const auto &sjt = data.yk.SixJ();
  const auto lk = dressed_table(data, data.i_orbs, 40);

  // Full sub-grid (stride 1), so the quadrature of Sigma_L is that of the
  // orbitals
  double worst_column = 0.0;
  double worst_diagonal = 0.0;
  double worst_asymmetry = 0.0;
  double smallest_change = 1.0;
  for (const auto &v : data.valence) {
    const auto Sigma = MBPT::Sigma_ladder_ratio2(
      v, data.core, data.excited, data.qk, lk, 1.0e-4, 200.0, 1, true);

    // Acting on v: the exact column, as from the projection onto the basis
    // before it is made symmetric (sum_x L_x|..><x|v> = L_v|..>)
    const auto ladder = MBPT::external_line_ladder(
      v, data.core, data.excited, data.qk, lk, {}, 60, 1.0e-13, 0.0, false);
    const auto column = MBPT::Sigma_ladder(v, ladder, data.core, data.excited,
                                           1.0e-4, 200.0, 1, true) *
                        v;
    const auto Sigma_v = Sigma * v;
    const auto difference = Sigma_v - column;
    const auto column_norm = std::sqrt(column * column);
    worst_column =
      std::max(worst_column, std::sqrt(difference * difference) / column_norm);

    const auto expected =
      MBPT::de_valence_w(v, data.qk, lk, data.core, data.excited, &sjt);
    const auto found = v * Sigma_v;
    worst_diagonal = std::max(worst_diagonal, std::abs(found / expected - 1.0));

    // Symmetric kernel: Sigma equals its transpose (measure on the column)
    const auto Sigma_T = Sigma.transpose_drj();
    double max_element = 0.0;
    double max_asymmetry = 0.0;
    for (std::size_t i = 0; i < Sigma.size(); ++i) {
      for (std::size_t j = 0; j < Sigma.size(); ++j) {
        max_element = std::max(max_element, std::abs(Sigma.ff(i, j)));
        max_asymmetry =
          std::max(max_asymmetry, std::abs(Sigma.ff(i, j) - Sigma_T.ff(i, j)));
      }
    }
    worst_asymmetry = std::max(worst_asymmetry, max_asymmetry / max_element);

    // The ratio method alone does not act exactly on v, so the correction
    // matters
    const auto Sigma_ratio = MBPT::Sigma_ladder(
      v, data.core, data.excited, data.qk, &lk, 1.0e-4, 200.0, 1, true);
    const auto ratio_difference = Sigma_ratio * v - column;
    smallest_change =
      std::min(smallest_change,
               std::sqrt(ratio_difference * ratio_difference) / column_norm);
  }
  REQUIRE(smallest_change > 1.0e-3);
  REQUIRE(worst_column < 1.0e-8);
  REQUIRE(worst_diagonal < 1.0e-8);
  REQUIRE(worst_asymmetry < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: lowest order with any x matches reference",
          "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  const auto max_k = stored_max_k(data);
  const auto excited = excited_up_to_n4(data);

  std::size_t num_core_x = 0;
  std::size_t num_outside_x = 0;
  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  for (const auto &v : data.valence) {
    // x: every orbital of kappa_v; a core x has the key of a core line
    const auto x_list = orbitals_of_kappa(v, data);
    for (const auto &x : x_list) {
      if (std::find(data.core.cbegin(), data.core.cend(), x) !=
          data.core.cend()) {
        ++num_core_x;
      } else if (std::find(excited.cbegin(), excited.cend(), x) ==
                 excited.cend()) {
        ++num_outside_x;
      }
    }

    // Each diagram alone, every entry, both orders of (m,n)
    const auto L2_mnxb = MBPT::lowest_order_Lk_mnxb(
      data.qk, excited, data.core, x_list, {false, true, false, false}, max_k);
    const auto L4_mnxb = MBPT::lowest_order_Lk_mnxb(
      data.qk, excited, data.core, x_list, {false, false, false, true}, max_k);
    const auto L1_xnab = MBPT::lowest_order_Lk_xnab(
      data.qk, excited, data.core, x_list, {true, false, false, false}, max_k);
    const auto L3_xnab = MBPT::lowest_order_Lk_xnab(
      data.qk, excited, data.core, x_list, {false, false, true, false}, max_k);

    const auto eps_L2 = table_vs_reference(
      L2_mnxb, excited, excited, x_list, data.core, max_k,
      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &x, const DiracSpinor &b) {
        return MBPT::L2(k, m, n, x, b, data.qk, data.core, excited, sjt);
      });
    const auto eps_L4 = table_vs_reference(
      L4_mnxb, excited, excited, x_list, data.core, max_k,
      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &x, const DiracSpinor &b) {
        return MBPT::L4(k, m, n, x, b, data.qk, data.core, sjt);
      });
    const auto eps_L1 = table_vs_reference(
      L1_xnab, x_list, excited, data.core, data.core, max_k,
      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b) {
        return MBPT::L1(k, x, n, a, b, data.qk, excited, sjt);
      });
    const auto eps_L3 = table_vs_reference(
      L3_xnab, x_list, excited, data.core, data.core, max_k,
      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b) {
        return MBPT::L3(k, x, n, a, b, data.qk, data.core, excited, sjt);
      });
    worst_mnxb = std::max({worst_mnxb, eps_L2, eps_L4});
    worst_xnab = std::max({worst_xnab, eps_L1, eps_L3});
  }
  REQUIRE(num_core_x > 0);
  REQUIRE(num_outside_x > 0);
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: mirror terms", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data(-1, 8.0);
  const auto max_k = stored_max_k(data);
  // Internal lines n <= 4: the projection states are every orbital of
  // kappa_v, so some are internal lines, some are not, and some are core
  const auto excited = excited_up_to_n4(data);

  fmt::print("\nLadder: <y|Sigma_mirror|v> and <v|Sigma_mirror|y> vs the six "
             "terms (au); x = every orbital of kappa_v\n");
  fmt::print("{:4s}  {:4s}  {:8s}  {:>13s}  {:>13s}  {:>8s}  {:>8s}\n", "v",
             "y", "y is", "expected", "found", "eps", "eps(vy)");
  // Worst eps for y = v, and for y an internal line, outside the internal
  // lines, and a core orbital
  double worst_diagonal = 0.0;
  double worst_internal = 0.0;
  double worst_outside = 0.0;
  double worst_core = 0.0;
  std::size_t num_internal = 0;
  std::size_t num_outside = 0;
  std::size_t num_core = 0;
  for (const auto &v : data.valence) {
    const auto x_list = orbitals_of_kappa(v, data);
    const auto Sigma =
      MBPT::Sigma_ladder_mirror(v, x_list, data.core, excited, data.qk, max_k,
                                true, 1.0e-4, 200.0, 1, true);
    const auto Sigma_v = Sigma * v;

    for (const auto &y : x_list) {
      // <y|Sigma_mirror|v>: the transpose of the terms with y in the ladder
      // integral (y one of the x). For y = v, the energy of the six terms,
      // whatever the other x
      const auto expected = mirror_terms_element(y, v, v.en(), excited, data);
      const auto found = y * Sigma_v;
      const auto eps = std::abs(found / expected - 1.0);
      // <v|Sigma_mirror|y>: the terms with y on the Coulomb vertex
      const auto expected_vy =
        mirror_terms_element(v, y, v.en(), excited, data);
      const auto eps_vy = std::abs((v * (Sigma * y)) / expected_vy - 1.0);
      const auto worst = std::max(eps, eps_vy);

      const bool is_core =
        std::find(data.core.cbegin(), data.core.cend(), y) != data.core.cend();
      const bool is_internal =
        std::find(excited.cbegin(), excited.cend(), y) != excited.cend();
      std::string kind = "v";
      if (y == v) {
        worst_diagonal = std::max(worst_diagonal, worst);
      } else if (is_core) {
        kind = "core";
        worst_core = std::max(worst_core, worst);
        ++num_core;
      } else if (is_internal) {
        kind = "internal";
        worst_internal = std::max(worst_internal, worst);
        ++num_internal;
      } else {
        kind = "outside";
        worst_outside = std::max(worst_outside, worst);
        ++num_outside;
      }
      fmt::print(
        "{:4s}  {:4s}  {:8s}  {:+13.6e}  {:+13.6e}  {:8.1e}  {:8.1e}\n",
        v.shortSymbol(), y.shortSymbol(), kind, expected, found, eps, eps_vy);
    }
  }
  REQUIRE(num_internal > 0);
  REQUIRE(num_outside > 0);
  REQUIRE(num_core > 0);
  REQUIRE(worst_diagonal < 1.0e-8);
  REQUIRE(worst_internal < 1.0e-8);
  REQUIRE(worst_outside < 1.0e-8);
  REQUIRE(worst_core < 1.0e-8);
}

//==============================================================================
TEST_CASE("Ladder: Coulomb integrals of the projection states",
          "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  // orbitals (the internal lines): core n = 2 and excited n <= 4. Of the
  // projection states, 1s (a core state) and 5d+ (above the excited orbitals)
  // are not among them
  const auto core_n2 =
    qip::select_if(data.core, [](const DiracSpinor &a) { return a.n() == 2; });
  const auto orbitals = qip::merge(core_n2, excited_up_to_n4(data));
  const auto basis = qip::merge(data.core, data.excited);

  Coulomb::QkTable qk;
  qk.fill(orbitals, data.yk, -1, false);
  const auto num_added =
    MBPT::fill_Qk_projection_states(&qk, orbitals, basis, data.valence, -1);

  // Every integral with one projection state x that is not in orbitals, and
  // three of orbitals: by the symmetry of Q, x may go first. Each
  // must be in the table, equal to that of the table of all orbitals. No
  // other integral is added: expected holds each distinct integral once
  Coulomb::QkTable expected;
  std::size_t num_core_x = 0;
  std::size_t num_outside_x = 0;
  std::size_t num_missing = 0;
  double worst = 0.0;
  double largest = 0.0;
  for (const auto &v : data.valence) {
    for (const auto &x : MBPT::projection_states(v, basis)) {
      if (std::find(orbitals.cbegin(), orbitals.cend(), x) != orbitals.cend())
        continue;
      if (std::find(data.core.cbegin(), data.core.cend(), x) !=
          data.core.cend()) {
        ++num_core_x;
      } else {
        ++num_outside_x;
      }
      for (const auto &a : orbitals) {
        for (const auto &b : orbitals) {
          for (const auto &c : orbitals) {
            const auto [k_min, k_max] = Coulomb::k_minmax_Q(x, a, b, c);
            for (int k = k_min; k <= k_max; k += 2) {
              const auto value = data.qk.Q(k, x, a, b, c);
              expected.add(k, x, a, b, c, value);
              if (!qk.contains(k, x, a, b, c)) {
                ++num_missing;
              }
              worst = std::max(worst, std::abs(qk.Q(k, x, a, b, c) - value));
              largest = std::max(largest, std::abs(value));
            }
          }
        }
      }
    }
  }
  const auto num_expected = expected.count();
  const auto eps = worst / largest;
  REQUIRE(num_core_x > 0);
  REQUIRE(num_outside_x > 0);
  REQUIRE(num_missing == 0);
  REQUIRE(num_added == num_expected);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
// Each piece of P with the external line, one sweep, every entry: the side
// with the line (dressed integrals from the integrals with the line, at
// en_v) plus the partner (the line only in a bare integral; dressed
// integrals from the ladder table)
TEST_CASE("Ladder: external line with P matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto lk = dressed_table(data, i_orbs, 1);
  // The integrals with the line start from a different table
  const auto lk_line = dressed_table(data, i_orbs, 2);

  // One piece of P, one side, with the reference function of that piece
  using Piece = std::function<double(
    int, const DiracSpinor &, const DiracSpinor &, const DiracSpinor &,
    const DiracSpinor &, const Coulomb::LkTable *, std::optional<double>,
    std::optional<double>)>;
  const Piece piece_1 =
    [&](int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
        const DiracSpinor &j, const Coulomb::LkTable *table,
        std::optional<double> e_i, std::optional<double> e_m) {
      return MBPT::P1(k, m, n, i, j, data.qk, data.core, data.excited, table,
                      e_i, e_m);
    };
  const Piece piece_2 =
    [&](int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
        const DiracSpinor &j, const Coulomb::LkTable *table,
        std::optional<double> e_i, std::optional<double> e_m) {
      return MBPT::P2(k, m, n, i, j, data.qk, data.core, data.excited, sjt,
                      table, e_i, e_m);
    };
  const Piece piece_3 =
    [&](int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
        const DiracSpinor &j, const Coulomb::LkTable *table,
        std::optional<double> e_i, std::optional<double> e_m) {
      return MBPT::P3(k, m, n, i, j, data.qk, data.core, data.excited, sjt,
                      table, e_i, e_m);
    };
  const std::vector<Piece> pieces{piece_1, piece_2, piece_3};

  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  for (const auto &v : data.valence) {
    const auto en_v = v.en() - 0.05;
    const auto start = ladder_from_table(v, data, lk_line);
    const auto &x_list = start.x_orbitals;
    for (std::size_t index = 0; index < pieces.size(); ++index) {
      MBPT::LadderDiagrams diagrams{false, false, false, false};
      diagrams.P1 = index == 0;
      diagrams.P2 = index == 1;
      diagrams.P3 = index == 2;
      const auto &piece = pieces[index];
      const auto ladder = external_line_sweep(start, en_v, data, lk, diagrams);

      // x as i: the row (m,x) holds the line
      const auto eps_mnxb =
        mnxb_difference(ladder, x_list, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &x, const DiracSpinor &b) {
                          return piece(k, m, n, x, b, &lk_line, en_v, {}) +
                                 piece(k, n, m, b, x, &lk, {}, {});
                        });
      // x as m: the row (x,a) holds the line. n also a projection state: in
      // the partner it is an internal line
      const auto eps_xnab =
        xnab_difference(ladder, x_list, data.excited, data,
                        [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                            const DiracSpinor &a, const DiracSpinor &b) {
                          return piece(k, x, n, a, b, &lk_line, {}, en_v) +
                                 piece(k, n, x, b, a, &lk, {}, {});
                        });
      worst_mnxb = std::max(worst_mnxb, eps_mnxb);
      worst_xnab = std::max(worst_xnab, eps_xnab);
    }
  }
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
// x_F with the external line, one sweep, every entry. x as i: the
// polarisation loop on either pair, and the hole-particle term on (n,b)
// only. x as m: the closing amplitude at en_v
TEST_CASE("Ladder: external line of x_F matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto diagrams = diagrams_of_x_F();
  const auto lk_F = dressed_table(data, i_orbs, 0, diagrams);
  const auto lk_F_line = dressed_table(data, i_orbs, 1, diagrams);

  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  double smallest_hole_particle = 1.0;
  for (const auto &v : data.valence) {
    const auto en_v = v.en() - 0.05;
    auto start = ladder_from_table(v, data, lk_F_line);
    // The closing amplitude has no symmetry: every entry its own value
    std::mt19937 generator(5);
    std::uniform_real_distribution<double> random_value(-0.05, 0.05);
    for (std::size_t k = 0; k < start.L_xnab->size(); ++k) {
      for (auto &entry : start.L_xnab->at(k)) {
        entry.second = random_value(generator);
      }
    }
    const auto &x_list = start.x_orbitals;

    auto ladder = start;
    MBPT::calculate_Lk_mnxb(&ladder.L_mnxb, data.qk, data.excited, data.core,
                            x_list, en_v, diagrams, lk_F, start.L_mnxb, 0.0);
    MBPT::calculate_Lk_vnab_F(&ladder.L_xnab, data.qk, data.excited, data.core,
                              x_list, diagrams, &lk_F, &start.L_xnab, 0.0,
                              en_v);

    // No hole-particle term on (m,x): there x is the valence line
    const auto x_F_mnxb = [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                              const DiracSpinor &x, const DiracSpinor &b) {
      return MBPT::P1(k, m, n, x, b, data.qk, data.core, data.excited,
                      &lk_F_line, en_v) +
             MBPT::Vhp(k, m, n, x, b, data.qk, data.core, data.excited,
                       &lk_F_line, en_v) +
             MBPT::P1(k, n, m, b, x, data.qk, data.core, data.excited, &lk_F);
    };
    worst_mnxb =
      std::max(worst_mnxb, mnxb_difference(ladder, x_list, data, x_F_mnxb));
    // The hole-particle term matters: without it the reference is far off
    smallest_hole_particle = std::min(
      smallest_hole_particle,
      mnxb_difference(ladder, x_list, data,
                      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &x, const DiracSpinor &b) {
                        return x_F_mnxb(k, m, n, x, b) -
                               MBPT::Vhp(k, m, n, x, b, data.qk, data.core,
                                         data.excited, &lk_F_line, en_v);
                      }));

    worst_xnab = std::max(
      worst_xnab,
      xnab_difference(ladder, x_list, data.excited, data,
                      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                          const DiracSpinor &a, const DiracSpinor &b) {
                        return MBPT::Lk_vnab_F(
                          k, x, n, a, b, data.qk, data.core, data.excited,
                          diagrams, &lk_F, &start.L_xnab, en_v);
                      }));
  }
  REQUIRE(smallest_hole_particle > 1.0e-3);
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
// The core relaxation with the external line, alone, every entry. x as i:
// the hole x is not dressed. x as m: the particle x is not dressed (no
// single excitations of the core into x), the particle n is
TEST_CASE("Ladder: external line with relaxation matches reference",
          "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto lk = dressed_table(data, i_orbs, 1);
  const auto rho_ma = random_core_singles(data);
  const MBPT::LadderDiagrams no_diagrams{false, false, false, false};

  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  double smallest_line_term = 1.0;
  for (const auto &v : data.valence) {
    const auto en_v = v.en() - 0.05;
    const auto start = ladder_from_table(v, data, lk);
    const auto &x_list = start.x_orbitals;
    const auto ladder =
      external_line_sweep(start, en_v, data, lk, no_diagrams, &rho_ma);

    worst_mnxb = std::max(
      worst_mnxb,
      mnxb_difference(ladder, x_list, data,
                      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &x, const DiracSpinor &b) {
                        return MBPT::Rkmnij(k, m, n, x, b, data.qk, rho_ma,
                                            data.core, data.excited);
                      }));
    worst_xnab = std::max(
      worst_xnab,
      xnab_difference(ladder, x_list, data.excited, data,
                      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                          const DiracSpinor &a, const DiracSpinor &b) {
                        return MBPT::Rkmnij(k, x, n, a, b, data.qk, rho_ma,
                                            data.core, data.excited, true);
                      }));
    // With the term of the line itself the integrals are different, for a
    // line with a core orbital of its kappa
    const auto with_line_term =
      xnab_difference(ladder, x_list, data.excited, data,
                      [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                          const DiracSpinor &a, const DiracSpinor &b) {
                        return MBPT::Rkmnij(k, x, n, a, b, data.qk, rho_ma,
                                            data.core, data.excited);
                      });
    const bool core_of_kappa_v = std::any_of(
      data.core.cbegin(), data.core.cend(),
      [&v](const DiracSpinor &c) { return c.kappa() == v.kappa(); });
    if (core_of_kappa_v) {
      smallest_line_term = std::min(smallest_line_term, with_line_term);
    }
  }
  REQUIRE(smallest_line_term > 1.0e-3);
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
// All diagrams, P and the core relaxation together: with en_v = e_x0 the x0
// integrals are those of a normal sweep, except for the single excitations
// of the core into x0, which the integrals with x0 as m do not have
TEST_CASE("Ladder: external line with P and relaxation at its own energy",
          "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto i_orbs = qip::merge(data.core, data.excited);
  const auto diagrams = diagrams_with_P();
  const auto lk = dressed_table(data, i_orbs, 1, diagrams);
  const auto rho_ma = random_core_singles(data);
  auto lk_next = lk;
  MBPT::calculate_Lk_mnib(&lk_next, data.qk, data.excited, data.core, i_orbs,
                          diagrams, &lk, max_k, false, 0.0, nullptr, &rho_ma);

  double worst_mnxb = 0.0;
  double worst_xnab = 0.0;
  for (const auto &v : data.valence) {
    const auto start = ladder_from_table(v, data, lk);
    for (const auto &x0 : start.x_orbitals) {
      const auto index_x0 =
        std::size_t(std::find(data.excited.cbegin(), data.excited.cend(), x0) -
                    data.excited.cbegin());
      const auto ladder =
        external_line_sweep(start, x0.en(), data, lk, diagrams, &rho_ma);
      worst_mnxb = std::max(
        worst_mnxb,
        mnxb_difference(ladder, {x0}, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &x, const DiracSpinor &b) {
                          return lk_next.Q(k, m, n, x, b);
                        }));
      // The normal sweep has -sum_c rho_xc Q^k_cnab: add it back
      worst_xnab = std::max(
        worst_xnab,
        xnab_difference(ladder, {x0}, data.excited, data,
                        [&](int k, const DiracSpinor &x, const DiracSpinor &n,
                            const DiracSpinor &a, const DiracSpinor &b) {
                          double line_term = 0.0;
                          for (std::size_t index_c = 0;
                               index_c < data.core.size(); ++index_c) {
                            const auto &c = data.core[index_c];
                            if (c.kappa() == x.kappa()) {
                              line_term += rho_ma(index_x0, index_c) *
                                           data.qk.Q(k, c, n, a, b);
                            }
                          }
                          return lk_next.Q(k, x, n, a, b) + line_term;
                        }));
    }
  }
  REQUIRE(worst_mnxb < 1.0e-12);
  REQUIRE(worst_xnab < 1.0e-12);
}

//==============================================================================
// The correction to the correlation potential from the full pair equation,
// converged: x = L + P + R closed with W, less x_F closed with Q, plus the
// core relaxation potential
TEST_CASE("Ladder: potential with P and relaxation", "[Ladder]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data(-1, 8.0);
  const auto &sjt = data.yk.SixJ();
  const auto max_k = stored_max_k(data);
  const auto &core = data.core;
  const auto &excited = data.excited;
  const auto diagrams = diagrams_with_P();
  const auto diagrams_F = diagrams_of_x_F();
  const int num_sweeps = 50;

  // x = L + P, and x = L + P + R with the core single excitations: both
  // iterated to convergence (the iteration converges geometrically)
  const auto lk_P = dressed_table(data, data.i_orbs, num_sweeps, diagrams);
  auto lk = dressed_table(data, data.i_orbs, 0, diagrams);
  auto rho_ma =
    MBPT::core_singles(MBPT::core_singles_source(data.qk, excited, core, &lk),
                       data.qk, excited, core, nullptr);
  for (int sweep = 0; sweep < num_sweeps; ++sweep) {
    const auto lk_prev = lk;
    MBPT::calculate_Lk_mnib(&lk, data.qk, excited, core, data.i_orbs, diagrams,
                            &lk_prev, max_k, false, 0.0, nullptr, &rho_ma);
    rho_ma =
      MBPT::core_singles(MBPT::core_singles_source(data.qk, excited, core, &lk),
                         data.qk, excited, core, &rho_ma);
  }
  const auto density = MBPT::relaxation_density(data.qk, core, excited);

  // Closing amplitude of x for each valence state, alone: no single
  // excitations of the core into v
  Coulomb::NkTable L_vnab;
  MBPT::add_Lk_vnab_F_entries(&L_vnab, excited, core, data.valence, max_k);
  for (int sweep = 0; sweep < num_sweeps; ++sweep) {
    const auto L_vnab_prev = L_vnab;
    for (const auto &v : data.valence) {
      MBPT::calculate_Lk_xnab(&L_vnab, data.qk, excited, core, {v}, v.en(),
                              diagrams, lk, L_vnab_prev, 0.0, &rho_ma);
    }
  }

  // x_F and its closing amplitude
  const auto lk_F = dressed_table(data, data.i_orbs, num_sweeps, diagrams_F);
  Coulomb::NkTable L_vnab_F;
  MBPT::add_Lk_vnab_F_entries(&L_vnab_F, excited, core, data.valence, max_k);
  for (int sweep = 0; sweep < num_sweeps; ++sweep) {
    const auto L_vnab_F_prev = L_vnab_F;
    MBPT::calculate_Lk_vnab_F(&L_vnab_F, data.qk, excited, core, data.valence,
                              diagrams_F, &lk_F, &L_vnab_F_prev, 0.0);
  }

  fmt::print("\nLadder: x = v integrals with the external line vs the valence "
             "integrals (eps);\n"
             "<v|dSigma|v> vs E_W[x] - E_Q[x_F] + <v|dV^R|v> (au);\n"
             "eps(<y|dSigma|x>) vs the sum over the tables\n");
  fmt::print("{:4s}  {:>7s}  {:>7s}  {:>7s}  {:>7s}  {:>13s}  {:>13s}  {:>7s}  "
             "{:>7s}\n",
             "v", "x(L+P)", "x(i)", "x(m)", "x_F", "expected", "found", "eps",
             "eps(yx)");
  double worst_column_P = 0.0;
  double worst_column = 0.0;
  double worst_column_F = 0.0;
  double worst_diagonal = 0.0;
  double worst_yx = 0.0;
  double smallest_line_term = 1.0;
  for (const auto &v : data.valence) {
    const auto from_table = [](const Coulomb::LkTable &table) {
      return [&table](int k, const DiracSpinor &a, const DiracSpinor &b,
                      const DiracSpinor &c,
                      const DiracSpinor &d) { return table.Q(k, a, b, c, d); };
    };
    const auto from_closing = [](const Coulomb::NkTable &table) {
      return [&table](int k, const DiracSpinor &a, const DiracSpinor &b,
                      const DiracSpinor &c,
                      const DiracSpinor &d) { return table.Q(k, a, b, c, d); };
    };

    // L + P: the x = v integrals are those of the ladder table, both slots
    const auto ladder_P = MBPT::external_line_ladder(
      v, core, excited, data.qk, lk_P, diagrams, 80, 1.0e-13, 0.0, false);
    const auto eps_column_P =
      std::max(mnxb_difference(ladder_P, {v}, data, from_table(lk_P)),
               xnab_difference(ladder_P, {v}, excited, data, from_table(lk_P)));

    // L + P + R: x as i from the ladder table; x as m from the closing table
    // (found for v alone), not from the ladder table
    const auto ladder =
      MBPT::external_line_ladder(v, core, excited, data.qk, lk, diagrams, 80,
                                 1.0e-13, 0.0, false, &rho_ma, &L_vnab);
    const auto eps_column_i =
      mnxb_difference(ladder, {v}, data, from_table(lk));
    const auto eps_column_m =
      xnab_difference(ladder, {v}, excited, data, from_closing(L_vnab));
    const auto line_term =
      xnab_difference(ladder, {v}, excited, data, from_table(lk));
    const bool core_of_kappa_v =
      std::any_of(core.cbegin(), core.cend(), [&v](const DiracSpinor &c) {
        return c.kappa() == v.kappa();
      });
    if (core_of_kappa_v) {
      smallest_line_term = std::min(smallest_line_term, line_term);
    }

    // x_F: x as i from the x_F table; x as m from its closing table
    const auto ladder_F =
      MBPT::external_line_ladder_F(v, core, excited, data.qk, lk_F, L_vnab_F,
                                   diagrams_F, 80, 1.0e-13, 0.0, false);
    const auto eps_column_F = std::max(
      mnxb_difference(ladder_F, {v}, data, from_table(lk_F)),
      xnab_difference(ladder_F, {v}, excited, data, from_closing(L_vnab_F)));

    // The potential, lower components included, on the full grid
    auto dSigma =
      MBPT::Sigma_ladder(v, ladder, core, excited, 1.0e-4, 200.0, 1, true);
    dSigma -= MBPT::Sigma_ladder(v, ladder_F.x_orbitals, ladder_F.L_mnxb,
                                 ladder_F.L_xnab, core, excited, 1.0e-4, 200.0,
                                 1, true, false);
    dSigma +=
      MBPT::relaxation_potential(ladder.x_orbitals, data.qk, rho_ma, density,
                                 core, excited, 1.0e-4, 200.0, 1, true);

    const auto [rho_part, occupation_part] =
      MBPT::de_valence_relaxation(v, data.qk, rho_ma, density, core, excited);
    const auto expected =
      MBPT::de_valence_w(v, data.qk, lk, L_vnab, core, excited, &sjt) -
      MBPT::de_valence_F(v, data.qk, lk_F, L_vnab_F, core, excited) + rho_part +
      occupation_part;
    const auto found = v * (dSigma * v);
    const auto eps = std::abs(found / expected - 1.0);

    // Every <y|dSigma|x>: W closing of x, less the Q closing of x_F, plus
    // <y|dV^R|x>
    const auto dV = MBPT::relaxation_potential_matrix(
      ladder.x_orbitals, data.qk, rho_ma, density, core, excited);
    double max_difference = 0.0;
    double largest = 0.0;
    for (std::size_t index_x = 0; index_x < ladder.x_orbitals.size();
         ++index_x) {
      const auto &x = ladder.x_orbitals[index_x];
      const auto dSigma_x = dSigma * x;
      for (std::size_t index_y = 0; index_y < ladder.x_orbitals.size();
           ++index_y) {
        const auto &y = ladder.x_orbitals[index_y];
        double sum = 0.0;
        for (const auto &n : excited) {
          for (const auto &a : core) {
            for (const auto &m : excited) {
              const auto de = v.en() + a.en() - m.en() - n.en();
              const auto [k_min, k_max] = Coulomb::k_minmax_W(y, a, m, n);
              for (int k = k_min; k <= k_max; ++k) {
                sum += (data.qk.W(k, y, a, m, n, &sjt) *
                          ladder.L_mnxb.Q(k, m, n, x, a) -
                        data.qk.Q(k, y, a, m, n) *
                          ladder_F.L_mnxb.Q(k, m, n, x, a)) /
                       ((2 * k + 1) * de);
              }
            }
            for (const auto &b : core) {
              const auto de = v.en() + n.en() - a.en() - b.en();
              const auto [k_min, k_max] = Coulomb::k_minmax_W(y, n, a, b);
              for (int k = k_min; k <= k_max; ++k) {
                sum += (data.qk.W(k, y, n, a, b, &sjt) *
                          ladder.L_xnab.Q(k, x, n, a, b) -
                        data.qk.Q(k, y, n, a, b) *
                          ladder_F.L_xnab.Q(k, x, n, a, b)) /
                       ((2 * k + 1) * de);
              }
            }
          }
        }
        sum = sum / v.twojp1() + dV(index_y, index_x);
        max_difference = std::max(max_difference, std::abs(y * dSigma_x - sum));
        largest = std::max(largest, std::abs(sum));
      }
    }
    const auto eps_yx = max_difference / largest;

    fmt::print("{:4s}  {:7.1e}  {:7.1e}  {:7.1e}  {:7.1e}  {:+13.6e}  "
               "{:+13.6e}  {:7.1e}  {:7.1e}\n",
               v.shortSymbol(), eps_column_P, eps_column_i, eps_column_m,
               eps_column_F, expected, found, eps, eps_yx);
    worst_column_P = std::max(worst_column_P, eps_column_P);
    worst_column = std::max({worst_column, eps_column_i, eps_column_m});
    worst_column_F = std::max(worst_column_F, eps_column_F);
    worst_diagonal = std::max(worst_diagonal, eps);
    worst_yx = std::max(worst_yx, eps_yx);
  }
  // With R, the integrals with x = v as m are not those of the ladder table
  REQUIRE(smallest_line_term > 1.0e-6);
  REQUIRE(worst_column_P < 1.0e-9);
  REQUIRE(worst_column < 1.0e-9);
  REQUIRE(worst_column_F < 1.0e-9);
  REQUIRE(worst_diagonal < 1.0e-8);
  REQUIRE(worst_yx < 1.0e-8);
}
