#include "MBPT/LadderIteration.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/Ladder.testdata.hpp"
#include "MBPT/LadderReference.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <vector>

namespace UnitTest::LadderTest {

//! Number of entries of table a that are not in table b
inline std::size_t missing_entries(Coulomb::LkTable &a, Coulomb::LkTable &b) {
  std::size_t num_missing = 0;
  for (std::size_t k = 0; k < a->size(); ++k) {
    for (const auto &[key, value] : a->at(k)) {
      if (b.get(int(k), key) == nullptr) {
        ++num_missing;
      }
    }
  }
  return num_missing;
}

//! Largest difference between the entries of two ladder tables, relative to
//! the largest entry of the first. Infinite if the sets of entries differ or
//! the first table is all zero
inline double ladder_table_difference(Coulomb::LkTable &expected,
                                      Coulomb::LkTable &found) {
  constexpr auto infinite = std::numeric_limits<double>::infinity();
  if (expected.count() != found.count() ||
      missing_entries(expected, found) != 0)
    return infinite;
  double max_difference = 0.0;
  double largest = 0.0;
  for (std::size_t k = 0; k < expected->size(); ++k) {
    for (const auto &[key, value] : expected->at(k)) {
      max_difference =
        std::max(max_difference, std::abs(value - *found.get(int(k), key)));
      largest = std::max(largest, std::abs(value));
    }
  }
  return largest == 0.0 ? infinite : max_difference / largest;
}

//! Adds a zero entry, one at a time, for every L^k_mnib (m,n excited; i in
//! i_orbs; b core; k <= max_k) not already in the table: the expected entry
//! set of add_Lk_mnib_entries()
inline void add_required_entries(Coulomb::LkTable *lk,
                                 const LadderTestData &data,
                                 const std::vector<DiracSpinor> &i_orbs,
                                 int max_k) {
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &i : i_orbs) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, i, b);
          for (int k = k0; k <= std::min(kI, max_k); ++k) {
            lk->add(k, m, n, i, b, 0.0);
          }
        }
      }
    }
  }
}

//! One sweep of the ladder equation on a copy of the start table, one
//! integral at a time with the reference MBPT::Lkmnij: the expected result of
//! calculate_Lk_mnib(). Entries from add_required_entries(); only_new and
//! damping rules as documented for calculate_Lk_mnib()
inline Coulomb::LkTable
reference_sweep(const Coulomb::LkTable &start, const LadderTestData &data,
                const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                const Coulomb::LkTable *lk_prev, int max_k, bool only_new,
                double a_damp) {
  auto lk_old = start;
  add_required_entries(&lk_old, data, i_orbs, max_k);
  auto lk = lk_old;
  const auto &sjt = data.yk.SixJ();
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &i : i_orbs) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, i, b);
          for (int k = k0; k <= std::min(kI, max_k); ++k) {
            const auto L_old = lk_old.Q(k, m, n, i, b);
            if (only_new && L_old != 0.0)
              continue;
            const auto L_new =
              MBPT::Lkmnij(k, m, n, i, b, data.qk, data.core, data.excited,
                           include_L4, sjt, lk_prev);
            // L_mnib and L_nmbi (both visited if i is core) are one entry,
            // and both visits write the same value
            lk.update(k, m, n, i, b,
                      (L_old == 0.0 || a_damp == 0.0) ?
                        L_new :
                        a_damp * L_old + (1.0 - a_damp) * L_new);
          }
        }
      }
    }
  }
  return lk;
}

//! One sweep of MBPT::calculate_Lk_mnib() on a copy of the start table
inline Coulomb::LkTable ladder_sweep(const Coulomb::LkTable &start,
                                     const LadderTestData &data,
                                     const std::vector<DiracSpinor> &i_orbs,
                                     const MBPT::LadderDiagrams &diagrams,
                                     const Coulomb::LkTable *lk_prev, int max_k,
                                     bool only_new, double a_damp) {
  auto lk = start;
  MBPT::add_Lk_mnib_entries(&lk, data.excited, data.core, i_orbs, max_k);
  MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, i_orbs,
                          diagrams, lk_prev, max_k, only_new, a_damp);
  return lk;
}

//! Lowest-order table of all diagrams, from the reference (one at a time)
inline Coulomb::LkTable reference_fill(const LadderTestData &data) {
  return reference_sweep({}, data, data.i_orbs, true, nullptr,
                         stored_max_k(data), true, 0.0);
}

//! Largest difference between the stored L^k_mnib (m,n excited; i in i_orbs;
//! b core) and the given one-at-a-time diagram, relative to the largest
inline double diagram_vs_scalar(
  const Coulomb::LkTable &lk, const LadderTestData &data,
  const std::function<double(int, const DiracSpinor &, const DiracSpinor &,
                             const DiracSpinor &, const DiracSpinor &)>
    &diagram) {
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &i : data.i_orbs) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, i, b);
          for (int k = k0; k <= std::min(kI, lk.max_k()); ++k) {
            const auto expected = diagram(k, m, n, i, b);
            const auto found = lk.Q(k, m, n, i, b);
            max_difference =
              std::max(max_difference, std::abs(expected - found));
            largest = std::max(largest, std::abs(expected));
          }
        }
      }
    }
  }
  return largest == 0.0 ? 1.0 : max_difference / largest;
}

} // namespace UnitTest::LadderTest

//==============================================================================
TEST_CASE("Ladder: 6j array", "[Ladder][unit]") {
  const int max_tj = 7;
  const int max_rank = 7;
  const MBPT::SixJArray sixj(max_tj, max_rank);
  double max_difference = 0.0;
  double largest = 0.0;
  for (int tja = 1; tja <= max_tj; tja += 2) {
    for (int tjb = 1; tjb <= max_tj; tjb += 2) {
      for (int c = 0; c <= max_rank; ++c) {
        for (int tjd = 1; tjd <= max_tj; tjd += 2) {
          for (int tje = 1; tje <= max_tj; tje += 2) {
            for (int f = 0; f <= max_rank; ++f) {
              const auto expected =
                Angular::sixj_2(tja, tjb, 2 * c, tjd, tje, 2 * f);
              const auto found = sixj.get(tja, tjb, c, tjd, tje, f);
              max_difference =
                std::max(max_difference, std::abs(expected - found));
              largest = std::max(largest, std::abs(expected));
            }
          }
        }
      }
    }
  }
  REQUIRE(largest > 0.1);
  REQUIRE(max_difference == 0.0);
}

//==============================================================================
TEST_CASE("Ladder: orbital pair index", "[Ladder][unit]") {
  using MBPT::OrbitalPairIndex;
  const auto data = UnitTest::LadderTest::ladder_test_data();
  const auto &excited = data.excited;
  const auto &core = data.core;
  const int max_rank = 3;
  const OrbitalPairIndex pairs(excited, core, max_rank);

  // Each pair is in exactly the blocks allowed by the triangle rule, and the
  // pairs of a block are numbered 0, 1, ..., size-1 (each number once)
  std::size_t wrong_membership = 0;
  std::size_t wrong_numbering = 0;
  for (int rank = 0; rank <= max_rank; ++rank) {
    for (int parity = 0; parity < 2; ++parity) {
      std::vector<int> times_used(pairs.size(rank, parity), 0);
      for (std::size_t ix = 0; ix < excited.size(); ++ix) {
        for (std::size_t iy = 0; iy < core.size(); ++iy) {
          const auto &x = excited[ix];
          const auto &y = core[iy];
          if ((x.l() + y.l()) % 2 != parity)
            continue;
          const bool triangle =
            Angular::triangle(x.twoj(), y.twoj(), 2 * rank) != 0;
          const auto index = pairs.index(rank, ix, iy);
          if (triangle != (index >= 0)) {
            ++wrong_membership;
          }
          if (index >= 0 && std::size_t(index) < times_used.size()) {
            ++times_used[std::size_t(index)];
          } else if (index >= 0) {
            ++wrong_numbering;
          }
        }
      }
      for (const auto times : times_used) {
        if (times != 1) {
          ++wrong_numbering;
        }
      }
    }
  }
  const auto num_blocks = pairs.num_blocks();
  REQUIRE(num_blocks == 2 * (std::size_t(max_rank) + 1));
  REQUIRE(wrong_membership == 0);
  REQUIRE(wrong_numbering == 0);

  // With a mask: only the allowed pairs are included
  std::vector<char> allowed(excited.size() * core.size(), 0);
  allowed[0] = 1;
  const OrbitalPairIndex one_pair(excited, core, max_rank, allowed);
  std::size_t total = 0;
  for (int rank = 0; rank <= max_rank; ++rank) {
    total += one_pair.size(rank, 0) + one_pair.size(rank, 1);
  }
  const auto [rank_min, rank_max] =
    Coulomb::k_minmax_tj(excited[0].twoj(), core[0].twoj());
  const auto expected_total =
    std::size_t(std::min(rank_max, max_rank) - rank_min + 1);
  REQUIRE(total == expected_total);
}

//==============================================================================
TEST_CASE("Ladder: multipoles", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  auto lk = reference_fill(data);

  // An integral with more k in the ladder range (both parities) than in the
  // Coulomb range: 3d+ 3d+ 3p+ 2p+
  const auto *m = DiracSpinor::find(3, -3, data.excited);
  const auto *i = DiracSpinor::find(3, -2, data.valence);
  const auto *b = DiracSpinor::find(2, -2, data.core);
  REQUIRE(m != nullptr);
  REQUIRE(i != nullptr);
  REQUIRE(b != nullptr);
  const auto Q = MBPT::coulomb_multipoles(*m, *m, *i, *b, data.qk);
  const auto QL = MBPT::QL_multipoles(*m, *m, *i, *b, data.qk, &lk);
  const auto [q_min, q_max] = Coulomb::k_minmax_Q(*m, *m, *i, *b);
  const auto [l_min, l_max] = MBPT::k_minmax_L(*m, *m, *i, *b);
  REQUIRE(q_max >= q_min);
  REQUIRE(l_max - l_min > (q_max - q_min) / 2);
  REQUIRE(Q.k_min == q_min);
  REQUIRE(Q.k_max == q_max);
  REQUIRE(Q.k_step == 2);
  REQUIRE(QL.k_min == l_min);
  REQUIRE(QL.k_max == l_max);
  REQUIRE(QL.k_step == 1);

  double max_difference_Q = 0.0;
  for (int k = q_min; k <= q_max; k += 2) {
    const auto expected = data.qk.Q(k, *m, *m, *i, *b);
    max_difference_Q =
      std::max(max_difference_Q, std::abs(expected - Q.value[std::size_t(k)]));
  }
  double max_difference_QL = 0.0;
  double largest_L = 0.0;
  for (int k = l_min; k <= l_max; ++k) {
    const auto expected =
      data.qk.Q(k, *m, *m, *i, *b) + lk.Q(k, *m, *m, *i, *b);
    max_difference_QL = std::max(max_difference_QL,
                                 std::abs(expected - QL.value[std::size_t(k)]));
    largest_L = std::max(largest_L, std::abs(lk.Q(k, *m, *m, *i, *b)));
  }
  REQUIRE(largest_L > 0.0);
  REQUIRE(max_difference_Q == 0.0);
  REQUIRE(max_difference_QL == 0.0);
}

//==============================================================================
TEST_CASE("Ladder: table entries", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const int max_k = 3;

  // An existing entry must be left unchanged
  const auto &m = data.excited.at(0);
  const auto &b = data.core.at(0);
  const auto [k_first, k_last] = MBPT::k_minmax_L(m, m, b, b);
  REQUIRE(k_last >= k_first);
  for (const auto *i_orbs : {&data.i_orbs, &data.core, &data.valence}) {
    Coulomb::LkTable lk_calculated;
    lk_calculated.add(k_first, m, m, b, b, 1.5);
    Coulomb::LkTable lk_expected;
    lk_expected.add(k_first, m, m, b, b, 1.5);
    MBPT::add_Lk_mnib_entries(&lk_calculated, data.excited, data.core, *i_orbs,
                              max_k);
    add_required_entries(&lk_expected, data, *i_orbs, max_k);
    const auto num_calculated = lk_calculated.count();
    const auto num_expected = lk_expected.count();
    const auto num_missing = missing_entries(lk_expected, lk_calculated);
    const auto kept_value = lk_calculated.Q(k_first, m, m, b, b);
    REQUIRE(num_calculated > 1);
    REQUIRE(num_calculated == num_expected);
    REQUIRE(num_missing == 0);
    REQUIRE(kept_value == 1.5);
  }
}

//==============================================================================
TEST_CASE("Ladder: fill matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  for (const bool include_L4 : {false, true}) {
    auto reference = reference_sweep({}, data, data.i_orbs, include_L4, nullptr,
                                     max_k, true, 0.0);
    auto calculated =
      ladder_sweep({}, data, data.i_orbs, {true, true, include_L4}, nullptr,
                   max_k, true, 0.0);
    const auto eps = ladder_table_difference(reference, calculated);
    REQUIRE(eps < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: core update matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  auto lk = reference_fill(data);

  // Core integrals only, from the filled table: valence integrals untouched
  auto reference =
    reference_sweep(lk, data, data.core, true, &lk, max_k, false, 0.0);
  auto calculated = ladder_sweep(lk, data, data.core, {true, true, true}, &lk,
                                 max_k, false, 0.0);
  const auto eps = ladder_table_difference(reference, calculated);
  const auto eps_changed = ladder_table_difference(lk, calculated);
  double max_change_valence = 0.0;
  for (const auto &m1 : data.excited) {
    for (const auto &n1 : data.excited) {
      for (const auto &v : data.valence) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m1, n1, v, b);
          for (int k = k0; k <= std::min(kI, max_k); ++k) {
            const auto change =
              std::abs(calculated.Q(k, m1, n1, v, b) - lk.Q(k, m1, n1, v, b));
            max_change_valence = std::max(max_change_valence, change);
          }
        }
      }
    }
  }
  REQUIRE(eps_changed > 1.0e-3);
  REQUIRE(eps < 1.0e-12);
  REQUIRE(max_change_valence == 0.0);
}

//==============================================================================
TEST_CASE("Ladder: valence update matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  auto lk = reference_fill(data);
  const auto lk_core =
    reference_sweep(lk, data, data.core, true, &lk, max_k, false, 0.0);

  // Valence orbitals in the i index, with the dressed core integrals
  for (const bool include_L4 : {false, true}) {
    auto reference = reference_sweep(lk_core, data, data.valence, include_L4,
                                     &lk_core, max_k, false, 0.0);
    auto calculated =
      ladder_sweep(lk_core, data, data.valence, {true, true, include_L4},
                   &lk_core, max_k, false, 0.0);
    const auto eps = ladder_table_difference(reference, calculated);
    REQUIRE(eps < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: damping matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const double a_damp = 0.4;
  auto lk = reference_fill(data);
  auto reference =
    reference_sweep(lk, data, data.i_orbs, true, &lk, max_k, false, a_damp);
  auto calculated = ladder_sweep(lk, data, data.i_orbs, {true, true, true}, &lk,
                                 max_k, false, a_damp);
  auto undamped = ladder_sweep(lk, data, data.i_orbs, {true, true, true}, &lk,
                               max_k, false, 0.0);
  const auto eps = ladder_table_difference(reference, calculated);
  const auto eps_undamped = ladder_table_difference(reference, undamped);
  REQUIRE(eps_undamped > 1.0e-3);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: only_new matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  auto lk = reference_fill(data);

  // Start: the filled table, with every third entry (by key) set to zero.
  // Only those are re-calculated, with the dressed integrals
  auto start = lk;
  std::size_t num_zeroed = 0;
  for (std::size_t k = 0; k < start->size(); ++k) {
    for (auto &[key, value] : start->at(k)) {
      if (key % 3 == 0) {
        value = 0.0;
        ++num_zeroed;
      }
    }
  }
  auto reference =
    reference_sweep(start, data, data.i_orbs, true, &lk, max_k, true, 0.0);
  auto calculated = ladder_sweep(start, data, data.i_orbs, {true, true, true},
                                 &lk, max_k, true, 0.0);
  auto all_new = ladder_sweep(start, data, data.i_orbs, {true, true, true}, &lk,
                              max_k, false, 0.0);
  const auto eps = ladder_table_difference(reference, calculated);
  const auto eps_all_new = ladder_table_difference(reference, all_new);
  REQUIRE(num_zeroed > 0);
  REQUIRE(num_zeroed < lk.count());
  REQUIRE(eps_all_new > 1.0e-3);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: multipole cut-off matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const int k_cut = 2;
  const auto data = ladder_test_data(k_cut);
  const auto max_k = stored_max_k(data);
  auto reference_lo =
    reference_sweep({}, data, data.i_orbs, true, nullptr, max_k, true, 0.0);
  auto calculated_lo = ladder_sweep({}, data, data.i_orbs, {true, true, true},
                                    nullptr, max_k, true, 0.0);
  auto reference = reference_sweep(reference_lo, data, data.core, true,
                                   &reference_lo, max_k, false, 0.0);
  auto calculated =
    ladder_sweep(reference_lo, data, data.core, {true, true, true},
                 &reference_lo, max_k, false, 0.0);
  const auto eps_fill = ladder_table_difference(reference_lo, calculated_lo);
  const auto eps_update = ladder_table_difference(reference, calculated);
  const auto lk_max_k = calculated.max_k();
  REQUIRE(max_k == k_cut);
  REQUIRE(lk_max_k == k_cut);
  REQUIRE(eps_fill < 1.0e-12);
  REQUIRE(eps_update < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: diagrams match reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto &sjt = data.yk.SixJ();
  const auto lk_filled = reference_fill(data);

  // Each diagram alone, at lowest order and with the dressed integrals
  const std::vector<const Coulomb::LkTable *> dressing{nullptr, &lk_filled};
  for (const auto *lk_prev : dressing) {
    const auto L1_only = ladder_sweep(
      {}, data, data.i_orbs, {true, false, false}, lk_prev, max_k, false, 0.0);
    const auto L2_L3_only = ladder_sweep(
      {}, data, data.i_orbs, {false, true, false}, lk_prev, max_k, false, 0.0);
    const auto L4_only = ladder_sweep(
      {}, data, data.i_orbs, {false, false, true}, lk_prev, max_k, false, 0.0);
    const auto eps_L1 = diagram_vs_scalar(
      L1_only, data,
      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &b) {
        return MBPT::L1(k, m, n, i, b, data.qk, data.excited, sjt, lk_prev);
      });
    const auto eps_L2_L3 =
      diagram_vs_scalar(L2_L3_only, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::L2(k, m, n, i, b, data.qk, data.core,
                                          data.excited, sjt, lk_prev) +
                                 MBPT::L3(k, m, n, i, b, data.qk, data.core,
                                          data.excited, sjt, lk_prev);
                        });
    const auto eps_L4 = diagram_vs_scalar(
      L4_only, data,
      [&](int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &b) {
        return MBPT::L4(k, m, n, i, b, data.qk, data.core, sjt, lk_prev);
      });
    REQUIRE(eps_L1 < 1.0e-12);
    REQUIRE(eps_L2_L3 < 1.0e-12);
    REQUIRE(eps_L4 < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: coupling round trips", "[Ladder][unit]") {
  using namespace MBPT;
  const auto data = UnitTest::LadderTest::ladder_test_data();

  // One orbital of each kappa: s, p-, p+, d-, d+, f-, f+
  std::vector<DiracSpinor> orbitals;
  for (const auto &group : MBPT::group_by_kappa(data.excited)) {
    orbitals.push_back(data.excited[group.front()]);
  }
  const auto max_tj = DiracSpinor::max_tj(orbitals);
  const SixJArray sixj(max_tj, max_tj);
  std::mt19937 generator(1);
  std::uniform_real_distribution<double> random_value(-1.0, 1.0);

  // Largest error of the round trip L^k -> coupled form -> L^k; and largest
  // error when the sum over the coupled rank is truncated (must be large)
  double error_pair = 0.0, error_pair_truncated = 0.0;
  double error_ph = 0.0, error_ph_truncated = 0.0;
  double error_coulomb = 0.0;
  std::size_t num_quadruples = 0;
  for (const auto &a : orbitals) {
    for (const auto &b : orbitals) {
      for (const auto &c : orbitals) {
        for (const auto &d : orbitals) {
          // all k of the triangles {a,c,k} and {b,d,k}, both parities
          Multipoles L_k;
          L_k.k_min = std::max(std::abs(a.twoj() - c.twoj()),
                               std::abs(b.twoj() - d.twoj())) /
                      2;
          L_k.k_max = std::min(a.twoj() + c.twoj(), b.twoj() + d.twoj()) / 2;
          if (L_k.k_max < L_k.k_min)
            continue;
          ++num_quadruples;
          for (int k = L_k.k_min; k <= L_k.k_max; ++k) {
            L_k.value[std::size_t(k)] = random_value(generator);
          }
          // Coulomb-like: every second k; the others hold values that must
          // not be used
          auto Q_k = L_k;
          Q_k.k_step = 2;
          for (int k = Q_k.k_min + 1; k <= Q_k.k_max; k += 2) {
            Q_k.value[std::size_t(k)] = 1.0e6;
          }

          // Pair coupling: J of the triangles {a,b,J} and {c,d,J}
          const auto J_min = std::max(std::abs(a.twoj() - b.twoj()),
                                      std::abs(c.twoj() - d.twoj())) /
                             2;
          const auto J_max =
            std::min(a.twoj() + b.twoj(), c.twoj() + d.twoj()) / 2;
          std::array<double, ladder_max_multipoles> L_J{}, Q_J{};
          for (int J = J_min; J <= J_max; ++J) {
            L_J[std::size_t(J)] = pair_coupled(J, a, b, c, d, L_k, sixj);
            Q_J[std::size_t(J)] = pair_coupled(J, a, b, c, d, Q_k, sixj);
          }
          // Particle-hole coupling: K of the triangles {a,d,K} and {b,c,K}
          const auto K_min = std::max(std::abs(a.twoj() - d.twoj()),
                                      std::abs(b.twoj() - c.twoj())) /
                             2;
          const auto K_max =
            std::min(a.twoj() + d.twoj(), b.twoj() + c.twoj()) / 2;
          std::array<double, ladder_max_multipoles> L_K{}, Q_K{};
          for (int K = K_min; K <= K_max; ++K) {
            L_K[std::size_t(K)] =
              particle_hole_coupled(K, a, b, c, d, L_k, sixj);
            Q_K[std::size_t(K)] =
              particle_hole_coupled(K, a, b, c, d, Q_k, sixj);
          }

          for (int k = L_k.k_min; k <= L_k.k_max; ++k) {
            const auto expected = L_k.value[std::size_t(k)];
            const auto expected_Q =
              (k - Q_k.k_min) % 2 == 0 ? Q_k.value[std::size_t(k)] : 0.0;
            error_pair = std::max(
              error_pair, std::abs(multipole_from_pair_coupled(
                                     k, a, b, c, d, L_J, J_min, J_max, sixj) -
                                   expected));
            error_ph = std::max(
              error_ph, std::abs(multipole_from_particle_hole_coupled(
                                   k, a, b, c, d, L_K, K_min, K_max, sixj) -
                                 expected));
            error_coulomb =
              std::max({error_coulomb,
                        std::abs(multipole_from_pair_coupled(
                                   k, a, b, c, d, Q_J, J_min, J_max, sixj) -
                                 expected_Q),
                        std::abs(multipole_from_particle_hole_coupled(
                                   k, a, b, c, d, Q_K, K_min, K_max, sixj) -
                                 expected_Q)});
            if (J_max > J_min) {
              error_pair_truncated = std::max(
                error_pair_truncated,
                std::abs(multipole_from_pair_coupled(k, a, b, c, d, L_J, J_min,
                                                     J_max - 1, sixj) -
                         expected));
            }
            if (K_max > K_min) {
              error_ph_truncated = std::max(
                error_ph_truncated,
                std::abs(multipole_from_particle_hole_coupled(
                           k, a, b, c, d, L_K, K_min, K_max - 1, sixj) -
                         expected));
            }
          }
        }
      }
    }
  }
  REQUIRE(orbitals.size() == 7);
  REQUIRE(num_quadruples > 100);
  REQUIRE(error_pair < 1.0e-13);
  REQUIRE(error_ph < 1.0e-13);
  REQUIRE(error_coulomb < 1.0e-13);
  REQUIRE(error_pair_truncated > 1.0e-2);
  REQUIRE(error_ph_truncated > 1.0e-2);
}
