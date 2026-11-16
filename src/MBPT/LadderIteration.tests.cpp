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
//! integral at a time with the reference MBPT::Lkmnij (and MBPT::Pkmnij if
//! include_P): the expected result of calculate_Lk_mnib(). Entries from
//! add_required_entries(); only_new and damping rules as documented for
//! calculate_Lk_mnib()
inline Coulomb::LkTable
reference_sweep(const Coulomb::LkTable &start, const LadderTestData &data,
                const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                const Coulomb::LkTable *lk_prev, int max_k, bool only_new,
                double a_damp, bool include_P = false) {
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
                           include_L4, sjt, lk_prev) +
              (include_P ? MBPT::Pkmnij(k, m, n, i, b, data.qk, data.core,
                                        data.excited, sjt, lk_prev) :
                           0.0);
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
      ladder_sweep({}, data, data.i_orbs, {true, true, true, include_L4},
                   nullptr, max_k, true, 0.0);
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
  auto calculated = ladder_sweep(lk, data, data.core, {true, true, true, true},
                                 &lk, max_k, false, 0.0);
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
      ladder_sweep(lk_core, data, data.valence, {true, true, true, include_L4},
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
  auto calculated = ladder_sweep(
    lk, data, data.i_orbs, {true, true, true, true}, &lk, max_k, false, a_damp);
  auto undamped = ladder_sweep(lk, data, data.i_orbs, {true, true, true, true},
                               &lk, max_k, false, 0.0);
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
  auto calculated = ladder_sweep(
    start, data, data.i_orbs, {true, true, true, true}, &lk, max_k, true, 0.0);
  auto all_new = ladder_sweep(start, data, data.i_orbs,
                              {true, true, true, true}, &lk, max_k, false, 0.0);
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
  auto calculated_lo = ladder_sweep(
    {}, data, data.i_orbs, {true, true, true, true}, nullptr, max_k, true, 0.0);
  auto reference = reference_sweep(reference_lo, data, data.core, true,
                                   &reference_lo, max_k, false, 0.0);
  auto calculated =
    ladder_sweep(reference_lo, data, data.core, {true, true, true, true},
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
    const auto L1_only =
      ladder_sweep({}, data, data.i_orbs, {true, false, false, false}, lk_prev,
                   max_k, false, 0.0);
    const auto L2_L3_only =
      ladder_sweep({}, data, data.i_orbs, {false, true, true, false}, lk_prev,
                   max_k, false, 0.0);
    const auto L4_only =
      ladder_sweep({}, data, data.i_orbs, {false, false, false, true}, lk_prev,
                   max_k, false, 0.0);
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
TEST_CASE("Ladder: P diagrams match reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto &sjt = data.yk.SixJ();
  const auto x_random = random_ladder_table(data);

  // Each piece of P alone (with its partner), i core and valence, at lowest
  // order and with a dressed integral that has both parities of k
  MBPT::LadderDiagrams P1_only{false, false, false, false};
  P1_only.P1 = true;
  MBPT::LadderDiagrams P2_only{false, false, false, false};
  P2_only.P2 = true;
  MBPT::LadderDiagrams P3_only{false, false, false, false};
  P3_only.P3 = true;
  const std::vector<const Coulomb::LkTable *> dressing{nullptr, &x_random};
  for (const auto *lk_prev : dressing) {
    const auto P1_table =
      ladder_sweep({}, data, data.i_orbs, P1_only, lk_prev, max_k, false, 0.0);
    const auto P2_table =
      ladder_sweep({}, data, data.i_orbs, P2_only, lk_prev, max_k, false, 0.0);
    const auto P3_table =
      ladder_sweep({}, data, data.i_orbs, P3_only, lk_prev, max_k, false, 0.0);
    const auto eps_P1 =
      diagram_vs_scalar(P1_table, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::P1(k, m, n, i, b, data.qk, data.core,
                                          data.excited, lk_prev) +
                                 MBPT::P1(k, n, m, b, i, data.qk, data.core,
                                          data.excited, lk_prev);
                        });
    const auto eps_P2 =
      diagram_vs_scalar(P2_table, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::P2(k, m, n, i, b, data.qk, data.core,
                                          data.excited, sjt, lk_prev) +
                                 MBPT::P2(k, n, m, b, i, data.qk, data.core,
                                          data.excited, sjt, lk_prev);
                        });
    const auto eps_P3 =
      diagram_vs_scalar(P3_table, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::P3(k, m, n, i, b, data.qk, data.core,
                                          data.excited, sjt, lk_prev) +
                                 MBPT::P3(k, n, m, b, i, data.qk, data.core,
                                          data.excited, sjt, lk_prev);
                        });
    REQUIRE(eps_P1 < 1.0e-12);
    REQUIRE(eps_P2 < 1.0e-12);
    REQUIRE(eps_P3 < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: P pieces match reference one side", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  using MBPT::OrbitalPairIndex;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  const auto x_random = random_ladder_table(data);

  // Factors and pieces as in calculate_Lk_mnib(), i core and valence, with a
  // dressed integral that has both parities of k
  const auto ij_orbitals =
    MBPT::ladder_ij_orbitals(data.core, data.excited, data.i_orbs);
  const auto max_tj = DiracSpinor::max_tj(qip::merge(data.core, data.excited));
  const MBPT::SixJArray sixj(max_tj, max_tj);
  const OrbitalPairIndex ph_pairs(data.excited, ij_orbitals, max_tj);
  const OrbitalPairIndex rc_pairs(data.excited, data.core, max_tj);
  const auto QLk_mi_rc_over_de =
    MBPT::QL_multipole_over_de(data.excited, data.core, ij_orbitals, ph_pairs,
                               rc_pairs, data.qk, &x_random);
  const auto Qk_rc_nj = MBPT::Q_multipole(data.excited, data.core, ij_orbitals,
                                          rc_pairs, ph_pairs, data.qk);
  const auto QK_cr_nj =
    MBPT::Q_particle_hole(data.excited, data.core, data.excited, ij_orbitals,
                          rc_pairs, ph_pairs, data.qk, sixj);
  const auto QLK_mi_rc_over_de = MBPT::QL_particle_hole_over_de(
    data.excited, data.core, ij_orbitals, ph_pairs, rc_pairs, data.qk,
    &x_random, sixj);
  const std::array<std::vector<LinAlg::Matrix<double>>, 3> pieces{
    MBPT::P1_multipole(QLk_mi_rc_over_de, Qk_rc_nj),
    MBPT::P2_multipole(QLk_mi_rc_over_de, QK_cr_nj),
    MBPT::P3_multipole(QLK_mi_rc_over_de, Qk_rc_nj)};
  MBPT::LadderDiagrams all_P{false, false, false, false};
  all_P.P1 = all_P.P2 = all_P.P3 = true;
  const auto P_sum =
    MBPT::P_multipole(data.excited, data.core, ij_orbitals, ph_pairs, data.qk,
                      &x_random, sixj, all_P);

  // The factors element by element: QL^k_mric/de at k of both parities, and
  // Q^k_cnrj with (n,j) = (m,i) at the Coulomb parity only (the blocks of the
  // other parity are empty)
  double max_difference_QL = 0.0, largest_QL = 0.0;
  double max_difference_Q = 0.0, largest_Q = 0.0;
  for (std::size_t index_m = 0; index_m < data.excited.size(); ++index_m) {
    const auto &m = data.excited[index_m];
    for (std::size_t index_i = 0; index_i < ij_orbitals.size(); ++index_i) {
      const auto &i = ij_orbitals[index_i];
      const auto parity = MBPT::pair_parity(m, i);
      for (std::size_t index_r = 0; index_r < data.excited.size(); ++index_r) {
        const auto &r = data.excited[index_r];
        for (std::size_t index_c = 0; index_c < data.core.size(); ++index_c) {
          const auto &c = data.core[index_c];
          if (MBPT::pair_parity(r, c) != parity)
            continue;
          const auto de = i.en() + c.en() - m.en() - r.en();
          const auto [k_min, k_max] = MBPT::k_minmax_L(m, r, i, c);
          for (int k = k_min; k <= k_max; ++k) {
            const auto block = OrbitalPairIndex::block(k, parity);
            const auto index_mi =
              std::size_t(ph_pairs.index(k, index_m, index_i));
            const auto index_rc =
              std::size_t(rc_pairs.index(k, index_r, index_c));
            const auto QL_over_de =
              (data.qk.Q(k, m, r, i, c) + x_random.Q(k, m, r, i, c)) / de;
            max_difference_QL =
              std::max(max_difference_QL,
                       std::abs(QLk_mi_rc_over_de[block](index_mi, index_rc) -
                                QL_over_de));
            largest_QL = std::max(largest_QL, std::abs(QL_over_de));
            if ((k + parity) % 2 != 0)
              continue;
            const auto Q = data.qk.Q(k, c, m, r, i);
            max_difference_Q =
              std::max(max_difference_Q,
                       std::abs(Qk_rc_nj[block](index_rc, index_mi) - Q));
            largest_Q = std::max(largest_Q, std::abs(Q));
          }
        }
      }
    }
  }
  std::size_t num_filled_other_parity = 0;
  for (int k = 0; k <= rc_pairs.max_rank(); ++k) {
    if (Qk_rc_nj[OrbitalPairIndex::block(k, (k + 1) % 2)].rows() != 0) {
      ++num_filled_other_parity;
    }
  }

  // Element (m,i),(n,j) at rank k: the reference piece of (m,n,i,j), one side
  // (the partner is the transposed element). j core: the first orbitals of
  // ij_orbitals. read_Pk(): the table entry P^k_mnij, one side plus partner
  std::array<double, 3> max_difference{}, largest{};
  double max_difference_sum = 0.0, largest_sum = 0.0;
  double max_difference_entry = 0.0, largest_entry = 0.0;
  for (std::size_t index_m = 0; index_m < data.excited.size(); ++index_m) {
    const auto &m = data.excited[index_m];
    for (std::size_t index_n = 0; index_n < data.excited.size(); ++index_n) {
      const auto &n = data.excited[index_n];
      for (std::size_t index_i = 0; index_i < ij_orbitals.size(); ++index_i) {
        const auto &i = ij_orbitals[index_i];
        for (std::size_t index_j = 0; index_j < data.core.size(); ++index_j) {
          const auto &j = ij_orbitals[index_j];
          const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, i, j);
          if (k_max < k_min)
            continue;
          const auto parity = MBPT::pair_parity(m, i);
          std::array<std::array<double, MBPT::ladder_max_multipoles>, 3> found;
          for (std::size_t p = 0; p < 3; ++p) {
            found[p] = MBPT::read_L2K_element(pieces[p], ph_pairs, ph_pairs,
                                              parity, index_m, index_i, index_n,
                                              index_j, k_min, k_max);
          }
          const auto found_sum =
            MBPT::read_L2K_element(P_sum, ph_pairs, ph_pairs, parity, index_m,
                                   index_i, index_n, index_j, k_min, k_max);
          const auto found_entry =
            MBPT::read_Pk(P_sum, ph_pairs, parity, index_m, index_i, index_n,
                          index_j, k_min, k_max);
          for (int k = k_min; k <= k_max; ++k) {
            const auto sk = std::size_t(k);
            const std::array<double, 3> expected{
              MBPT::P1(k, m, n, i, j, data.qk, data.core, data.excited,
                       &x_random),
              MBPT::P2(k, m, n, i, j, data.qk, data.core, data.excited, sjt,
                       &x_random),
              MBPT::P3(k, m, n, i, j, data.qk, data.core, data.excited, sjt,
                       &x_random)};
            for (std::size_t p = 0; p < 3; ++p) {
              max_difference[p] = std::max(
                max_difference[p], std::abs(found[p][sk] - expected[p]));
              largest[p] = std::max(largest[p], std::abs(expected[p]));
            }
            const auto sum = found[0][sk] + found[1][sk] + found[2][sk];
            max_difference_sum =
              std::max(max_difference_sum, std::abs(found_sum[sk] - sum));
            largest_sum = std::max(largest_sum, std::abs(sum));
            const auto entry = MBPT::Pkmnij(k, m, n, i, j, data.qk, data.core,
                                            data.excited, sjt, &x_random);
            max_difference_entry =
              std::max(max_difference_entry, std::abs(found_entry[sk] - entry));
            largest_entry = std::max(largest_entry, std::abs(entry));
          }
        }
      }
    }
  }
  REQUIRE(num_filled_other_parity == 0);
  REQUIRE(max_difference_QL / largest_QL < 1.0e-14);
  REQUIRE(max_difference_Q / largest_Q < 1.0e-14);
  REQUIRE(largest[0] > 0.0);
  REQUIRE(largest[1] > 0.0);
  REQUIRE(largest[2] > 0.0);
  REQUIRE(max_difference[0] / largest[0] < 1.0e-12);
  REQUIRE(max_difference[1] / largest[1] < 1.0e-12);
  REQUIRE(max_difference[2] / largest[2] < 1.0e-12);
  REQUIRE(max_difference_sum / largest_sum < 1.0e-14);
  REQUIRE(largest_entry > 0.0);
  REQUIRE(max_difference_entry / largest_entry < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: signed block products", "[Ladder][unit]") {
  using MBPT::OrbitalPairIndex;
  using Matrix = LinAlg::Matrix<double>;
  // Blocks of ranks 0 and 1: (1 2) times (3 4)^T = 11, times (-1)^rank. Rank 0
  // in the block of parity 1 and rank 1 in that of parity 0: the sign follows
  // the rank. A block with an empty factor stays empty
  std::vector<Matrix> left(4), right(4);
  left[OrbitalPairIndex::block(0, 1)] = Matrix{{1.0, 2.0}};
  right[OrbitalPairIndex::block(0, 1)] = Matrix{{3.0}, {4.0}};
  left[OrbitalPairIndex::block(1, 0)] = Matrix{{1.0, 2.0}};
  right[OrbitalPairIndex::block(1, 0)] = Matrix{{3.0}, {4.0}};
  left[OrbitalPairIndex::block(1, 1)] = Matrix{{1.0}};
  const auto products = MBPT::signed_block_products(left, right);
  const auto &rank_0 = products[OrbitalPairIndex::block(0, 1)];
  const auto &rank_1 = products[OrbitalPairIndex::block(1, 0)];
  REQUIRE(rank_0.rows() == 1);
  REQUIRE(rank_0.cols() == 1);
  REQUIRE(rank_1.rows() == 1);
  REQUIRE(rank_1.cols() == 1);
  REQUIRE(rank_0(0, 0) == 11.0);
  REQUIRE(rank_1(0, 0) == -11.0);
  REQUIRE(products[OrbitalPairIndex::block(0, 0)].rows() == 0);
  REQUIRE(products[OrbitalPairIndex::block(1, 1)].rows() == 0);
}

//==============================================================================
TEST_CASE("Ladder: P fill matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  // All diagrams and all pieces of P, lowest order, i core and valence
  const MBPT::LadderDiagrams all{true, true, true, true, true, true, true};
  auto reference = reference_sweep({}, data, data.i_orbs, true, nullptr, max_k,
                                   true, 0.0, true);
  auto calculated =
    ladder_sweep({}, data, data.i_orbs, all, nullptr, max_k, true, 0.0);
  auto without_P = reference_fill(data);
  const auto eps = ladder_table_difference(reference, calculated);
  const auto eps_without_P = ladder_table_difference(reference, without_P);
  REQUIRE(eps_without_P > 1.0e-3);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: P update matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto x_random = random_ladder_table(data);
  const double a_damp = 0.3;
  // All diagrams and all pieces of P, damped, with a dressed integral that has
  // both parities of k: the core sweep, then the valence sweep (which reads
  // the core entries too)
  const MBPT::LadderDiagrams all{true, true, true, true, true, true, true};
  for (const auto *i_orbs : {&data.core, &data.valence}) {
    auto reference = reference_sweep(x_random, data, *i_orbs, true, &x_random,
                                     max_k, false, a_damp, true);
    auto calculated = ladder_sweep(x_random, data, *i_orbs, all, &x_random,
                                   max_k, false, a_damp);
    const auto eps = ladder_table_difference(reference, calculated);
    REQUIRE(eps < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: Vhp diagram matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto x_random = random_ladder_table(data);

  // The hole-particle interaction alone, and with the polarisation loop (they
  // share a factor), i core and valence, at lowest order and with a dressed
  // integral that has both parities of k. Table entry: one side plus partner
  MBPT::LadderDiagrams Vhp_only{false, false, false, false};
  Vhp_only.hole_particle = true;
  auto P1_Vhp = Vhp_only;
  P1_Vhp.P1 = true;
  const std::vector<const Coulomb::LkTable *> dressing{nullptr, &x_random};
  for (const auto *lk_prev : dressing) {
    const auto Vhp_table =
      ladder_sweep({}, data, data.i_orbs, Vhp_only, lk_prev, max_k, false, 0.0);
    const auto P1_Vhp_table =
      ladder_sweep({}, data, data.i_orbs, P1_Vhp, lk_prev, max_k, false, 0.0);
    const auto eps_Vhp =
      diagram_vs_scalar(Vhp_table, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::Vhp(k, m, n, i, b, data.qk, data.core,
                                           data.excited, lk_prev) +
                                 MBPT::Vhp(k, n, m, b, i, data.qk, data.core,
                                           data.excited, lk_prev);
                        });
    const auto eps_P1_Vhp =
      diagram_vs_scalar(P1_Vhp_table, data,
                        [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &b) {
                          return MBPT::P1(k, m, n, i, b, data.qk, data.core,
                                          data.excited, lk_prev) +
                                 MBPT::P1(k, n, m, b, i, data.qk, data.core,
                                          data.excited, lk_prev) +
                                 MBPT::Vhp(k, m, n, i, b, data.qk, data.core,
                                           data.excited, lk_prev) +
                                 MBPT::Vhp(k, n, m, b, i, data.qk, data.core,
                                           data.excited, lk_prev);
                        });
    REQUIRE(eps_Vhp < 1.0e-12);
    REQUIRE(eps_P1_Vhp < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: Vhp update matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto x_random = random_ladder_table(data);
  const double a_damp = 0.3;

  // x_F (polarisation loop and hole-particle interaction), damped, with a
  // dressed integral that has both parities of k: the core sweep, then the
  // valence sweep (which reads the core entries too). Table entry: one side
  // plus partner, damped as a whole
  MBPT::LadderDiagrams P1_Vhp{false, false, false, false};
  P1_Vhp.P1 = true;
  P1_Vhp.hole_particle = true;
  double worst_eps = 0.0;
  for (const auto *i_orbs : {&data.core, &data.valence}) {
    const auto calculated = ladder_sweep(x_random, data, *i_orbs, P1_Vhp,
                                         &x_random, max_k, false, a_damp);
    double max_difference = 0.0, largest = 0.0;
    for (const auto &m : data.excited) {
      for (const auto &n : data.excited) {
        for (const auto &i : *i_orbs) {
          for (const auto &b : data.core) {
            const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, i, b);
            for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
              const auto L_new = MBPT::P1(k, m, n, i, b, data.qk, data.core,
                                          data.excited, &x_random) +
                                 MBPT::P1(k, n, m, b, i, data.qk, data.core,
                                          data.excited, &x_random) +
                                 MBPT::Vhp(k, m, n, i, b, data.qk, data.core,
                                           data.excited, &x_random) +
                                 MBPT::Vhp(k, n, m, b, i, data.qk, data.core,
                                           data.excited, &x_random);
              const auto L_old = x_random.Q(k, m, n, i, b);
              const auto expected =
                L_old == 0.0 ? L_new : a_damp * L_old + (1.0 - a_damp) * L_new;
              const auto found = calculated.Q(k, m, n, i, b);
              max_difference =
                std::max(max_difference, std::abs(found - expected));
              largest = std::max(largest, std::abs(expected));
            }
          }
        }
      }
    }
    worst_eps =
      std::max(worst_eps, largest == 0.0 ? 1.0 : max_difference / largest);
  }
  REQUIRE(worst_eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: Vhp factors", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  using MBPT::OrbitalPairIndex;
  const auto data = ladder_test_data();
  const auto x_random = random_ladder_table(data);

  // Factors as in calculate_Lk_mnib(), i core and valence, with a dressed
  // integral that has both parities of k
  const auto ij_orbitals =
    MBPT::ladder_ij_orbitals(data.core, data.excited, data.i_orbs);
  const auto max_tj = DiracSpinor::max_tj(qip::merge(data.core, data.excited));
  const OrbitalPairIndex ph_pairs(data.excited, ij_orbitals, max_tj);
  const OrbitalPairIndex rc_pairs(data.excited, data.core, max_tj);
  const auto Y_rc_nj = MBPT::Vhp_matrix(data.excited, data.core, ij_orbitals,
                                        rc_pairs, ph_pairs, data.qk);
  const auto QLk_mi_rc_over_de =
    MBPT::QL_multipole_over_de(data.excited, data.core, ij_orbitals, ph_pairs,
                               rc_pairs, data.qk, &x_random);
  const auto Hk_mi_nj = MBPT::Vhp_multipole(QLk_mi_rc_over_de, Y_rc_nj);

  // Y at every rank of both parities: -y^(c)_nr = -R^0_ncrc for the hole
  // unchanged (j = c) and kappa_r = kappa_n; zero otherwise, so zero in the
  // columns with a valence hole
  double max_difference_Y = 0.0, largest_Y = 0.0;
  std::size_t num_nonzero_valence_hole = 0;
  for (std::size_t index_r = 0; index_r < data.excited.size(); ++index_r) {
    const auto &r = data.excited[index_r];
    for (std::size_t index_c = 0; index_c < data.core.size(); ++index_c) {
      const auto &c = data.core[index_c];
      const auto parity = MBPT::pair_parity(r, c);
      for (std::size_t index_n = 0; index_n < data.excited.size(); ++index_n) {
        const auto &n = data.excited[index_n];
        for (std::size_t index_j = 0; index_j < ij_orbitals.size(); ++index_j) {
          const auto &j = ij_orbitals[index_j];
          if (MBPT::pair_parity(n, j) != parity)
            continue;
          const bool hole_unchanged = j == c && r.kappa() == n.kappa();
          const auto expected =
            hole_unchanged ? -data.qk.R(0, n, c, r, c) : 0.0;
          for (int K = 0; K <= ph_pairs.max_rank(); ++K) {
            const auto row = rc_pairs.index(K, index_r, index_c);
            const auto col = ph_pairs.index(K, index_n, index_j);
            if (row < 0 || col < 0)
              continue;
            const auto found = Y_rc_nj[OrbitalPairIndex::block(K, parity)](
              std::size_t(row), std::size_t(col));
            max_difference_Y =
              std::max(max_difference_Y, std::abs(found - expected));
            largest_Y = std::max(largest_Y, std::abs(expected));
            if (index_j >= data.core.size() && found != 0.0) {
              ++num_nonzero_valence_hole;
            }
          }
        }
      }
    }
  }

  // H^k_{mi,nj} = (A Y)^k: element (m,i),(n,j) at rank k is the reference
  // one side, j core
  double max_difference_H = 0.0, largest_H = 0.0;
  for (std::size_t index_m = 0; index_m < data.excited.size(); ++index_m) {
    const auto &m = data.excited[index_m];
    for (std::size_t index_n = 0; index_n < data.excited.size(); ++index_n) {
      const auto &n = data.excited[index_n];
      for (std::size_t index_i = 0; index_i < ij_orbitals.size(); ++index_i) {
        const auto &i = ij_orbitals[index_i];
        for (std::size_t index_j = 0; index_j < data.core.size(); ++index_j) {
          const auto &j = ij_orbitals[index_j];
          const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, i, j);
          if (k_max < k_min)
            continue;
          const auto found = MBPT::read_L2K_element(
            Hk_mi_nj, ph_pairs, ph_pairs, MBPT::pair_parity(m, i), index_m,
            index_i, index_n, index_j, k_min, k_max);
          for (int k = k_min; k <= k_max; ++k) {
            const auto expected = MBPT::Vhp(k, m, n, i, j, data.qk, data.core,
                                            data.excited, &x_random);
            max_difference_H = std::max(
              max_difference_H, std::abs(found[std::size_t(k)] - expected));
            largest_H = std::max(largest_H, std::abs(expected));
          }
        }
      }
    }
  }

  REQUIRE(largest_Y > 0.0);
  REQUIRE(max_difference_Y == 0.0);
  REQUIRE(num_nonzero_valence_hole == 0);
  REQUIRE(largest_H > 0.0);
  REQUIRE(max_difference_H / largest_H < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: closing amplitude factor", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  using MBPT::OrbitalPairIndex;
  const auto data = ladder_test_data();
  const auto L_vnab_random = random_closing_table(data);

  // Rows (v,a) of the closing amplitude, as in calculate_Lk_vnab_F():
  // (Q + L')^k_vrac/de, with a dressed integral that has k of both parities
  const auto max_tj = DiracSpinor::max_tj(qip::merge(data.core, data.excited));
  const OrbitalPairIndex rc_pairs(data.excited, data.core, max_tj);
  const OrbitalPairIndex va_pairs(data.valence, data.core, max_tj);
  const auto QLk_va_rc_over_de = MBPT::QL_multipole_vrac_over_de(
    data.valence, data.excited, data.core, va_pairs, rc_pairs, data.qk,
    &L_vnab_random);

  double max_difference_QL = 0.0, largest_QL = 0.0;
  for (std::size_t index_v = 0; index_v < data.valence.size(); ++index_v) {
    const auto &v = data.valence[index_v];
    for (std::size_t index_a = 0; index_a < data.core.size(); ++index_a) {
      const auto &a = data.core[index_a];
      const auto parity = MBPT::pair_parity(v, a);
      for (std::size_t index_r = 0; index_r < data.excited.size(); ++index_r) {
        const auto &r = data.excited[index_r];
        for (std::size_t index_c = 0; index_c < data.core.size(); ++index_c) {
          const auto &c = data.core[index_c];
          if (MBPT::pair_parity(r, c) != parity)
            continue;
          const auto de = a.en() + c.en() - v.en() - r.en();
          const auto [k_min, k_max] = MBPT::k_minmax_L(v, r, a, c);
          for (int k = k_min; k <= k_max; ++k) {
            const auto va = std::size_t(va_pairs.index(k, index_v, index_a));
            const auto rc = std::size_t(rc_pairs.index(k, index_r, index_c));
            const auto expected =
              (data.qk.Q(k, v, r, a, c) + L_vnab_random.Q(k, v, r, a, c)) / de;
            const auto found =
              QLk_va_rc_over_de[OrbitalPairIndex::block(k, parity)](va, rc);
            max_difference_QL =
              std::max(max_difference_QL, std::abs(found - expected));
            largest_QL = std::max(largest_QL, std::abs(expected));
          }
        }
      }
    }
  }
  REQUIRE(largest_QL > 0.0);
  REQUIRE(max_difference_QL / largest_QL < 1.0e-14);
}

//==============================================================================
TEST_CASE("Ladder: x_F closing matches reference", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);
  const auto x_random = random_ladder_table(data);
  const auto L_vnab_random = random_closing_table(data);

  // Every L'^k_vnab (v valence; n excited; a,b core) is an entry
  std::size_t num_expected = 0;
  for (const auto &v : data.valence) {
    for (const auto &n : data.excited) {
      for (const auto &a : data.core) {
        for (const auto &b : data.core) {
          const auto [k_min, k_max] = MBPT::k_minmax_L(v, n, a, b);
          num_expected += std::size_t(std::max(0, k_max - k_min + 1));
        }
      }
    }
  }

  // The polarisation loop, the hole-particle interaction, and both; at lowest
  // order, and from random tables (both parities of k) with damping
  MBPT::LadderDiagrams P1_only{false, false, false, false};
  P1_only.P1 = true;
  MBPT::LadderDiagrams Vhp_only{false, false, false, false};
  Vhp_only.hole_particle = true;
  auto P1_Vhp = P1_only;
  P1_Vhp.hole_particle = true;
  std::size_t num_wrong_count = 0;
  double worst_eps = 0.0;
  for (const auto *diagrams : {&P1_only, &Vhp_only, &P1_Vhp}) {
    for (const bool dressed : {false, true}) {
      const auto *lk_F = dressed ? &x_random : nullptr;
      const auto *L_vnab_prev = dressed ? &L_vnab_random : nullptr;
      const double a_damp = dressed ? 0.3 : 0.0;
      auto L_vnab = dressed ? L_vnab_random : Coulomb::NkTable{};
      MBPT::add_Lk_vnab_F_entries(&L_vnab, data.excited, data.core,
                                  data.valence, max_k);
      if (L_vnab.count() != num_expected) {
        ++num_wrong_count;
      }
      MBPT::calculate_Lk_vnab_F(&L_vnab, data.qk, data.excited, data.core,
                                data.valence, *diagrams, lk_F, L_vnab_prev,
                                a_damp);
      double max_difference = 0.0, largest = 0.0;
      for (const auto &v : data.valence) {
        for (const auto &n : data.excited) {
          for (const auto &a : data.core) {
            for (const auto &b : data.core) {
              const auto [k_min, k_max] = MBPT::k_minmax_L(v, n, a, b);
              for (int k = k_min; k <= k_max; ++k) {
                const auto L_new =
                  MBPT::Lk_vnab_F(k, v, n, a, b, data.qk, data.core,
                                  data.excited, *diagrams, lk_F, L_vnab_prev);
                const auto L_old =
                  dressed ? L_vnab_random.Q(k, v, n, a, b) : 0.0;
                const auto expected = L_old == 0.0 ?
                                        L_new :
                                        a_damp * L_old + (1.0 - a_damp) * L_new;
                const auto found = L_vnab.Q(k, v, n, a, b);
                max_difference =
                  std::max(max_difference, std::abs(found - expected));
                largest = std::max(largest, std::abs(expected));
              }
            }
          }
        }
      }
      worst_eps =
        std::max(worst_eps, largest == 0.0 ? 1.0 : max_difference / largest);
    }
  }
  REQUIRE(num_expected > 0);
  REQUIRE(num_wrong_count == 0);
  REQUIRE(worst_eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: x_F closing vs core amplitude", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);

  // Without the hole-particle interaction, the closing amplitude is the core
  // amplitude x_F_vnab at every order (two sweeps here). With it, at lowest
  // order, it is the core amplitude without the term on (v,a), H^k_nvba
  MBPT::LadderDiagrams P1_only{false, false, false, false};
  P1_only.P1 = true;
  auto P1_Vhp = P1_only;
  P1_Vhp.hole_particle = true;

  const auto lk_1 =
    ladder_sweep({}, data, data.i_orbs, P1_only, nullptr, max_k, false, 0.0);
  const auto lk_2 =
    ladder_sweep({}, data, data.i_orbs, P1_only, &lk_1, max_k, false, 0.0);
  Coulomb::NkTable L_1, L_2;
  MBPT::add_Lk_vnab_F_entries(&L_1, data.excited, data.core, data.valence,
                              max_k);
  MBPT::add_Lk_vnab_F_entries(&L_2, data.excited, data.core, data.valence,
                              max_k);
  MBPT::calculate_Lk_vnab_F(&L_1, data.qk, data.excited, data.core,
                            data.valence, P1_only, nullptr, nullptr, 0.0);
  MBPT::calculate_Lk_vnab_F(&L_2, data.qk, data.excited, data.core,
                            data.valence, P1_only, &lk_1, &L_1, 0.0);

  const auto lk_Vhp =
    ladder_sweep({}, data, data.i_orbs, P1_Vhp, nullptr, max_k, false, 0.0);
  Coulomb::NkTable L_Vhp;
  MBPT::add_Lk_vnab_F_entries(&L_Vhp, data.excited, data.core, data.valence,
                              max_k);
  MBPT::calculate_Lk_vnab_F(&L_Vhp, data.qk, data.excited, data.core,
                            data.valence, P1_Vhp, nullptr, nullptr, 0.0);

  double max_difference_loop = 0.0, largest_loop = 0.0;
  double max_difference_Vhp = 0.0, largest_Vhp_va = 0.0;
  for (const auto &v : data.valence) {
    for (const auto &n : data.excited) {
      for (const auto &a : data.core) {
        for (const auto &b : data.core) {
          const auto [k_min, k_max] = MBPT::k_minmax_L(v, n, a, b);
          for (int k = k_min; k <= k_max; ++k) {
            const auto core_loop = lk_2.Q(k, v, n, a, b);
            max_difference_loop = std::max(
              max_difference_loop, std::abs(L_2.Q(k, v, n, a, b) - core_loop));
            largest_loop = std::max(largest_loop, std::abs(core_loop));
            const auto Vhp_va = MBPT::Vhp(k, n, v, b, a, data.qk, data.core,
                                          data.excited, nullptr);
            const auto expected = lk_Vhp.Q(k, v, n, a, b) - Vhp_va;
            max_difference_Vhp = std::max(
              max_difference_Vhp, std::abs(L_Vhp.Q(k, v, n, a, b) - expected));
            largest_Vhp_va = std::max(largest_Vhp_va, std::abs(Vhp_va));
          }
        }
      }
    }
  }
  REQUIRE(largest_loop > 0.0);
  REQUIRE(max_difference_loop / largest_loop < 1.0e-12);
  REQUIRE(largest_Vhp_va > 0.0);
  REQUIRE(max_difference_Vhp / largest_Vhp_va < 1.0e-12);
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
