#include "MBPT/LadderIteration.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/Ladder.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace UnitTest {

//! Small set of hydrogen-like orbitals, with Coulomb integrals: a pure
//! angular/algebra test bed for the ladder iteration (no Hartree-Fock needed).
struct LadderTestData {
  //! core: n <= 2 (1s, 2s, 2p-, 2p+), so hole pairs have several J
  std::vector<DiracSpinor> core;
  //! excited: n = 3, 4 (up to f states)
  std::vector<DiracSpinor> excited;
  //! "valence": 3s, 3p+, 3d+ (members of excited, as in the ladder driver)
  std::vector<DiracSpinor> valence;
  //! core + valence
  std::vector<DiracSpinor> i_orbs;
  Coulomb::YkTable yk;
  Coulomb::QkTable qk;
};

//! Forms the test data; Qk is filled up to k_cut (<= 0: no cut-off).
//! Hydrogen-like orbitals of nuclear charge zeff: the ladder iteration
//! converges for large zeff (Coulomb integrals ~ zeff, energy gaps ~ zeff^2)
inline LadderTestData ladder_test_data(int k_cut = -1, double zeff = 1.0,
                                       std::size_t num_points = 500) {
  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{num_points, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto orbitals = DiracSpinor::HlikeBasis(3, 3, radial_grid, zeff);
  LadderTestData data{{}, {}, {}, {}, Coulomb::YkTable(orbitals), {}};
  for (const auto &Fn : orbitals) {
    (Fn.n() <= 2 ? data.core : data.excited).push_back(Fn);
    const bool valence_kappa =
      Fn.kappa() == -1 || Fn.kappa() == -2 || Fn.kappa() == -3;
    if (Fn.n() == 3 && valence_kappa) {
      data.valence.push_back(Fn);
    }
  }
  data.i_orbs = qip::merge(data.core, data.valence);
  data.qk.fill(orbitals, data.yk, k_cut, false);
  return data;
}

//! Largest difference between the stored L^k_mnib (m,n excited; i in i_orbs;
//! b core) and Lkmnij() with the given rung; relative to the largest integral
inline double ladder_table_vs_scalar(const Coulomb::LkTable &lk,
                                     const LadderTestData &data,
                                     const std::vector<DiracSpinor> &i_orbs,
                                     const Coulomb::LkTable *rung,
                                     bool include_L4) {
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &i : i_orbs) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, i, b);
          for (int k = k0; k <= std::min(kI, lk.max_k()); ++k) {
            const auto expected =
              MBPT::Lkmnij(k, m, n, i, b, data.qk, data.core, data.excited,
                           include_L4, data.yk.SixJ(), rung);
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

//! Pointer to the orbital (n, kappa) of the list; nullptr if not present
inline const DiracSpinor *find_orbital(const std::vector<DiracSpinor> &orbitals,
                                       int n, int kappa) {
  for (const auto &Fn : orbitals) {
    if (Fn.n() == n && Fn.kappa() == kappa) {
      return &Fn;
    }
  }
  return nullptr;
}

} // namespace UnitTest

//==============================================================================
TEST_CASE("Ladder: iteration 6j array", "[Ladder][unit]") {
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
TEST_CASE("Ladder: iteration orbital pair index", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  const auto excited = MBPT::ladder_orbitals(data.excited);
  const auto core = MBPT::ladder_orbitals(data.core);
  const auto &Fn = data.excited.front();
  const auto &orbital = excited.front();
  const bool same_orbital =
    orbital.nk_index == Fn.nk_index() && orbital.label == Fn.nk_index() &&
    orbital.kappa == Fn.kappa() && orbital.twoj == Fn.twoj() &&
    orbital.l == Fn.l() && orbital.en == Fn.en();
  const auto num_excited = excited.size();
  const auto expected_num_excited = data.excited.size();
  const auto max_tj = MBPT::max_twoj(excited);
  const auto expected_max_tj = DiracSpinor::max_tj(data.excited);
  REQUIRE(num_excited == expected_num_excited);
  REQUIRE(same_orbital);
  REQUIRE(max_tj == expected_max_tj);

  // s_1/2 and p_3/2: odd parity; (-1)^(rank + 1/2 + 3/2) = (-1)^rank
  const MBPT::LadderOrbital s_half{0, 0, -1, 1, 0, 0.0};
  const MBPT::LadderOrbital p_three_halves{0, 0, -2, 3, 1, 0.0};
  const auto parity_sp = MBPT::pair_parity(s_half, p_three_halves);
  const auto parity_ss = MBPT::pair_parity(s_half, s_half);
  const auto sign_rank_1 = MBPT::pair_sign(1, s_half, p_three_halves);
  const auto sign_rank_2 = MBPT::pair_sign(2, s_half, p_three_halves);
  REQUIRE(parity_sp == 1);
  REQUIRE(parity_ss == 0);
  REQUIRE(sign_rank_1 == -1);
  REQUIRE(sign_rank_2 == 1);

  const int max_rank = 3;
  const MBPT::OrbitalPairIndex pairs(excited, core, max_rank);

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
          if ((x.l + y.l) % 2 != parity)
            continue;
          const bool triangle =
            Angular::triangle(x.twoj, y.twoj, 2 * rank) != 0;
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
  const auto pairs_max_rank = pairs.max_rank();
  const auto num_blocks = pairs.num_blocks();
  const auto block_2_odd = MBPT::OrbitalPairIndex::block(2, 1);
  const auto num_pairs_rank_1_even = pairs.size(1, 0);
  REQUIRE(pairs_max_rank == max_rank);
  REQUIRE(num_blocks == 2 * (std::size_t(max_rank) + 1));
  REQUIRE(block_2_odd == 5);
  REQUIRE(num_pairs_rank_1_even > 0);
  REQUIRE(wrong_membership == 0);
  REQUIRE(wrong_numbering == 0);

  // With a mask: only the allowed pairs are included
  std::vector<char> allowed(excited.size() * core.size(), 0);
  allowed[0] = 1;
  const MBPT::OrbitalPairIndex one_pair(excited, core, max_rank, allowed);
  std::size_t total = 0;
  for (int rank = 0; rank <= max_rank; ++rank) {
    total += one_pair.size(rank, 0) + one_pair.size(rank, 1);
  }
  const auto [rank_min, rank_max] =
    Coulomb::k_minmax_tj(excited[0].twoj, core[0].twoj);
  const auto expected_total =
    std::size_t(std::min(rank_max, max_rank) - rank_min + 1);
  REQUIRE(total == expected_total);
}

//==============================================================================
TEST_CASE("Ladder: iteration multipoles", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs, true,
                     -1, false);

  // An integral with Coulomb-allowed k, and with more k in the ladder range
  // (both parities) than in the Coulomb range: 3d+ 3d+ 3p+ 2p+
  const auto *pm = UnitTest::find_orbital(data.excited, 3, -3);
  const auto *pi = UnitTest::find_orbital(data.valence, 3, -2);
  const auto *pb = UnitTest::find_orbital(data.core, 2, -2);
  REQUIRE(pm != nullptr);
  REQUIRE(pi != nullptr);
  REQUIRE(pb != nullptr);
  const auto &m = *pm;
  const auto &n = *pm;
  const auto &i = *pi;
  const auto &b = *pb;
  const auto orbitals = MBPT::ladder_orbitals({m, n, i, b});
  const auto Q = MBPT::coulomb_multipoles(orbitals[0], orbitals[1], orbitals[2],
                                          orbitals[3], data.qk);
  const auto T = MBPT::rung_multipoles(orbitals[0], orbitals[1], orbitals[2],
                                       orbitals[3], data.qk, &lk);
  const auto [q_min, q_max] = Coulomb::k_minmax_Q(m, n, i, b);
  const auto [l_min, l_max] = MBPT::k_minmax_L(m, n, i, b);
  REQUIRE(q_max >= q_min);
  REQUIRE(l_max - l_min > (q_max - q_min) / 2);
  REQUIRE(Q.k_min == q_min);
  REQUIRE(Q.k_max == q_max);
  REQUIRE(T.k_min == l_min);
  REQUIRE(T.k_max == l_max);

  double max_difference_Q = 0.0;
  for (int k = q_min; k <= q_max; ++k) {
    const auto expected = data.qk.Q(k, m, n, i, b);
    max_difference_Q =
      std::max(max_difference_Q, std::abs(expected - Q.value[std::size_t(k)]));
  }
  double max_difference_T = 0.0;
  double largest_L = 0.0;
  for (int k = l_min; k <= l_max; ++k) {
    const auto expected = data.qk.Q(k, m, n, i, b) + lk.Q(k, m, n, i, b);
    max_difference_T =
      std::max(max_difference_T, std::abs(expected - T.value[std::size_t(k)]));
    largest_L = std::max(largest_L, std::abs(lk.Q(k, m, n, i, b)));
  }
  REQUIRE(largest_L > 0.0);
  REQUIRE(max_difference_Q == 0.0);
  REQUIRE(max_difference_T == 0.0);
}

//==============================================================================
TEST_CASE("Ladder: iteration table entries", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  const int max_k = 3;

  // Existing entries must be left unchanged
  Coulomb::LkTable lk;
  const auto &m = data.excited.at(0);
  const auto &b = data.core.at(0);
  const auto [k_first, k_last] = MBPT::k_minmax_L(m, m, b, b);
  REQUIRE(k_last >= k_first);
  lk.add(k_first, m, m, b, b, 1.5);

  MBPT::add_Lk_mnib_entries(&lk, data.excited, data.core, data.i_orbs, max_k);

  // Every required integral is present, and nothing else: count the required
  // integrals, each once (L^k_mnib = L^k_nmbi)
  std::size_t num_required = 0;
  std::size_t num_missing = 0;
  double sum_of_values = 0.0;
  for (const auto &m1 : data.excited) {
    for (const auto &n1 : data.excited) {
      for (const auto &i1 : data.i_orbs) {
        for (const auto &b1 : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m1, n1, i1, b1);
          for (int k = k0; k <= std::min(kI, max_k); ++k) {
            if (!lk.contains(k, m1, n1, i1, b1)) {
              ++num_missing;
            }
            const bool i_is_core =
              std::find(data.core.cbegin(), data.core.cend(), i1) !=
              data.core.cend();
            // (m,n,i,b) and (n,m,b,i) are both in this loop if i is core
            if (!i_is_core || lk.is_NormalOrdered(m1, n1, i1, b1)) {
              ++num_required;
              sum_of_values += lk.Q(k, m1, n1, i1, b1);
            }
          }
        }
      }
    }
  }
  const auto table_max_k = lk.max_k();
  const auto num_stored = lk.count();
  REQUIRE(table_max_k == max_k);
  REQUIRE(num_required > 0);
  REQUIRE(num_missing == 0);
  REQUIRE(num_stored == num_required);
  REQUIRE(sum_of_values == 1.5);

  // Calculate: only_new leaves the existing (non-zero) entry alone, and
  // fills all the others
  MBPT::calculate_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                          true, nullptr, max_k, true, 0.0);
  const auto kept_value = lk.Q(k_first, m, m, b, b);
  lk.update(k_first, m, m, b, b,
            MBPT::Lkmnij(k_first, m, m, b, b, data.qk, data.core, data.excited,
                         true, data.yk.SixJ(), nullptr));
  const auto eps =
    UnitTest::ladder_table_vs_scalar(lk, data, data.i_orbs, nullptr, true);
  REQUIRE(kept_value == 1.5);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: iteration fill matches scalar", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  for (const bool include_L4 : {false, true}) {
    Coulomb::LkTable lk;
    MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                       include_L4, -1, false);
    const auto eps = UnitTest::ladder_table_vs_scalar(lk, data, data.i_orbs,
                                                      nullptr, include_L4);
    const auto num_stored = lk.count();
    REQUIRE(num_stored > 0);
    REQUIRE(eps < 1.0e-12);
  }
}

//==============================================================================
TEST_CASE("Ladder: iteration update matches scalar", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  const bool include_L4 = true;
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                     include_L4, -1, false);

  // Core integrals only: valence integrals must be untouched
  auto lk_core = lk;
  MBPT::update_Lk_mnib(&lk_core, data.qk, data.excited, data.core, data.core,
                       include_L4, data.yk.SixJ(), &lk, 0.0, false);
  const auto eps_core =
    UnitTest::ladder_table_vs_scalar(lk_core, data, data.core, &lk, include_L4);
  // the valence integrals are exactly as before
  double max_change_valence = 0.0;
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &v : data.valence) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, v, b);
          for (int k = k0; k <= kI; ++k) {
            const auto change =
              std::abs(lk_core.Q(k, m, n, v, b) - lk.Q(k, m, n, v, b));
            max_change_valence = std::max(max_change_valence, change);
          }
        }
      }
    }
  }
  REQUIRE(eps_core < 1.0e-12);
  REQUIRE(max_change_valence == 0.0);

  // Then the valence integrals, with the rung from the previous table
  auto lk_all = lk_core;
  MBPT::update_Lk_mnib(&lk_all, data.qk, data.excited, data.core, data.valence,
                       include_L4, data.yk.SixJ(), &lk_core, 0.0, false);
  const auto eps_valence = UnitTest::ladder_table_vs_scalar(
    lk_all, data, data.valence, &lk_core, include_L4);
  REQUIRE(eps_valence < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: iteration damping", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  const bool include_L4 = false;
  const double damp = 0.4;
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                     include_L4, -1, false);
  auto lk_damped = lk;
  MBPT::update_Lk_mnib(&lk_damped, data.qk, data.excited, data.core, data.core,
                       include_L4, data.yk.SixJ(), &lk, damp, false);

  // L -> damp*L_old + (1-damp)*L_new
  double max_difference = 0.0;
  double largest = 0.0;
  for (const auto &m : data.excited) {
    for (const auto &n : data.excited) {
      for (const auto &a : data.core) {
        for (const auto &b : data.core) {
          const auto [k0, kI] = MBPT::k_minmax_L(m, n, a, b);
          for (int k = k0; k <= kI; ++k) {
            const auto L_new =
              MBPT::Lkmnij(k, m, n, a, b, data.qk, data.core, data.excited,
                           include_L4, data.yk.SixJ(), &lk);
            const auto expected =
              damp * lk.Q(k, m, n, a, b) + (1.0 - damp) * L_new;
            const auto found = lk_damped.Q(k, m, n, a, b);
            max_difference =
              std::max(max_difference, std::abs(expected - found));
            largest = std::max(largest, std::abs(expected));
          }
        }
      }
    }
  }
  REQUIRE(largest > 0.0);
  REQUIRE(max_difference < 1.0e-12 * largest);
}

//==============================================================================
TEST_CASE("Ladder: iteration multipole cut-off", "[Ladder][unit]") {
  const int max_k = 2;
  const auto data = UnitTest::ladder_test_data(max_k);
  const bool include_L4 = true;
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                     include_L4, max_k, false);
  const auto eps_fill = UnitTest::ladder_table_vs_scalar(lk, data, data.i_orbs,
                                                         nullptr, include_L4);
  auto lk_next = lk;
  MBPT::update_Lk_mnib(&lk_next, data.qk, data.excited, data.core, data.core,
                       include_L4, data.yk.SixJ(), &lk, 0.0, false);
  const auto eps_update =
    UnitTest::ladder_table_vs_scalar(lk_next, data, data.core, &lk, include_L4);
  const auto qk_max_k = data.qk.max_k();
  const auto lk_max_k = lk.max_k();
  REQUIRE(qk_max_k == max_k);
  REQUIRE(lk_max_k == max_k);
  REQUIRE(eps_fill < 1.0e-12);
  REQUIRE(eps_update < 1.0e-12);
}

//==============================================================================
namespace UnitTest {

//! Ladder table for the test data: filled, then iterated num_iterations times
//! (core integrals, then valence integrals)
inline Coulomb::LkTable ladder_test_table(const LadderTestData &data,
                                          bool include_L4, int num_iterations) {
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.i_orbs,
                     include_L4, -1, false);
  for (const auto *update_i : {&data.core, &data.valence}) {
    for (int it = 0; it < num_iterations; ++it) {
      const auto lk_prev = lk;
      MBPT::update_Lk_mnib(&lk, data.qk, data.excited, data.core, *update_i,
                           include_L4, data.yk.SixJ(), &lk_prev, 0.0, false);
    }
  }
  return lk;
}

} // namespace UnitTest

//==============================================================================
TEST_CASE("Ladder: external line first iteration", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data();
  const bool include_L4 = true;
  const auto lk = UnitTest::ladder_test_table(data, include_L4, 1);
  const auto projection = qip::merge(data.core, data.excited);

  // One iteration, starting from the stored integrals: must be Lkmnij with
  // the rung from the table and the external energy e_v
  const auto integrals = MBPT::external_line_integrals(
    data.valence, data.core, data.excited, projection, data.qk, lk, include_L4,
    1, 0.0, 0.0, false);

  double max_difference_i = 0.0, largest_i = 0.0;
  double max_difference_m = 0.0, largest_m = 0.0;
  std::size_t num_external = 0;
  for (std::size_t iv = 0; iv < data.valence.size(); ++iv) {
    const auto en_v = data.valence[iv].en();
    for (const auto &[x, label] : integrals.lines.of_valence[iv]) {
      ++num_external;
      for (const auto &n : data.excited) {
        for (const auto &a : data.core) {
          // external line in the i-slot: L^k_mnxa
          for (const auto &m : data.excited) {
            const auto [k0, kI] = MBPT::k_minmax_L(m, n, *x, a);
            for (int k = k0; k <= kI; ++k) {
              const auto expected =
                MBPT::Lkmnij(k, m, n, *x, a, data.qk, data.core, data.excited,
                             include_L4, data.yk.SixJ(), &lk, en_v);
              const auto found = integrals.i_slot.Q(
                k, m.nk_index(), n.nk_index(), label, a.nk_index());
              max_difference_i =
                std::max(max_difference_i, std::abs(expected - found));
              largest_i = std::max(largest_i, std::abs(expected));
            }
          }
          // external line in the m-slot: L^k_xnab
          for (const auto &b : data.core) {
            const auto [k0, kI] = MBPT::k_minmax_L(*x, n, a, b);
            for (int k = k0; k <= kI; ++k) {
              const auto expected =
                MBPT::Lkmnij(k, *x, n, a, b, data.qk, data.core, data.excited,
                             include_L4, data.yk.SixJ(), &lk, {}, en_v);
              const auto found = integrals.m_slot.Q(k, label, n.nk_index(),
                                                    a.nk_index(), b.nk_index());
              max_difference_m =
                std::max(max_difference_m, std::abs(expected - found));
              largest_m = std::max(largest_m, std::abs(expected));
            }
          }
        }
      }
    }
  }
  REQUIRE(num_external > data.valence.size());
  REQUIRE(largest_i > 0.0);
  REQUIRE(largest_m > 0.0);
  REQUIRE(max_difference_i < 1.0e-12 * largest_i);
  REQUIRE(max_difference_m < 1.0e-12 * largest_m);
}

//==============================================================================
TEST_CASE("Ladder: external line converged", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data(-1, 8.0);
  const bool include_L4 = true;
  // converged table (iteration converges geometrically)
  const auto lk = UnitTest::ladder_test_table(data, include_L4, 40);
  const auto projection = qip::merge(data.core, data.excited);
  const auto integrals = MBPT::external_line_integrals(
    data.valence, data.core, data.excited, projection, data.qk, lk, include_L4,
    60, 1.0e-14, 0.0, false);

  // For x = v, the integrals with the external line are the stored ones
  double max_difference = 0.0, largest = 0.0;
  for (std::size_t iv = 0; iv < data.valence.size(); ++iv) {
    const auto &v = data.valence[iv];
    for (const auto &[x, label] : integrals.lines.of_valence[iv]) {
      if (!(*x == v))
        continue;
      for (const auto &n : data.excited) {
        for (const auto &a : data.core) {
          for (const auto &m : data.excited) {
            const auto [k0, kI] = MBPT::k_minmax_L(m, n, v, a);
            for (int k = k0; k <= kI; ++k) {
              const auto expected = lk.Q(k, m, n, v, a);
              const auto found = integrals.i_slot.Q(
                k, m.nk_index(), n.nk_index(), label, a.nk_index());
              max_difference =
                std::max(max_difference, std::abs(expected - found));
              largest = std::max(largest, std::abs(expected));
            }
          }
          for (const auto &b : data.core) {
            const auto [k0, kI] = MBPT::k_minmax_L(v, n, a, b);
            for (int k = k0; k <= kI; ++k) {
              const auto expected = lk.Q(k, v, n, a, b);
              const auto found = integrals.m_slot.Q(k, label, n.nk_index(),
                                                    a.nk_index(), b.nk_index());
              max_difference =
                std::max(max_difference, std::abs(expected - found));
              largest = std::max(largest, std::abs(expected));
            }
          }
        }
      }
    }
  }
  REQUIRE(largest > 0.0);
  REQUIRE(max_difference < 1.0e-9 * largest);
}

//==============================================================================
TEST_CASE("Ladder: iterated potential diagonal", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data(-1, 8.0);
  const bool include_L4 = true;
  const auto lk = UnitTest::ladder_test_table(data, include_L4, 40);
  const auto projection = qip::merge(data.core, data.excited);

  const auto Sigmas = MBPT::Sigma_ladder_iterated(
    data.valence, data.core, data.excited, projection, data.qk, lk, include_L4,
    60, 1.0e-12, 0.0, 1.0e-4, 60.0, 1, true, false);
  const auto num_Sigmas = Sigmas.size();
  REQUIRE(num_Sigmas == data.valence.size());

  // <v|Sigma_L|v> is the ladder energy of v (with the lower components
  // included in Sigma_L, as they are in the energy)
  double worst_eps = 0.0;
  for (std::size_t iv = 0; iv < data.valence.size(); ++iv) {
    const auto &v = data.valence[iv];
    const auto expected = MBPT::de_valence_w(v, data.qk, lk, data.core,
                                             data.excited, &data.yk.SixJ());
    const auto found = v * (Sigmas[iv] * v);
    worst_eps = std::max(worst_eps, std::abs(found / expected - 1.0));
  }
  REQUIRE(worst_eps < 1.0e-4);
}

//==============================================================================
namespace UnitTest {

//! Sum of the six lowest-order ladder terms whose mirror images the ladder
//! equation does not generate, between bra and ket: closings (a) of L2;
//! (a),(b) of L4; (c),(d) of L1; (c) of L3; at external energy en.
//! bra = ket = v gives their contribution to the ladder energy of v.
inline double unpaired_ladder_terms(const DiracSpinor &bra,
                                    const DiracSpinor &ket, double en,
                                    const LadderTestData &data) {
  const auto &sjt = data.yk.SixJ();
  const auto tjp1 = bra.twojp1();
  double sum = 0.0;
  for (const auto &n : data.excited) {
    for (const auto &a : data.core) {
      for (const auto &m : data.excited) {
        const auto de = en + a.en() - m.en() - n.en();
        const auto [k0, kI] = MBPT::k_minmax_L(m, n, ket, a);
        for (int k = k0; k <= kI; ++k) {
          const auto f = (2 * k + 1) * tjp1 * de;
          // (a) of L2: direct closing
          sum += data.qk.Q(k, bra, a, m, n) *
                 MBPT::L2(k, m, n, ket, a, data.qk, data.core, data.excited,
                          sjt, nullptr) /
                 f;
          // (a),(b) of L4: direct and exchange closings
          sum += data.qk.W(k, bra, a, m, n, &sjt) *
                 MBPT::L4(k, m, n, ket, a, data.qk, data.core, sjt, nullptr) /
                 f;
        }
      }
      for (const auto &b : data.core) {
        const auto de = en + n.en() - a.en() - b.en();
        const auto [k0, kI] = MBPT::k_minmax_L(ket, n, a, b);
        for (int k = k0; k <= kI; ++k) {
          const auto f = (2 * k + 1) * tjp1 * de;
          // (c),(d) of L1
          sum +=
            data.qk.W(k, bra, n, a, b, &sjt) *
            MBPT::L1(k, ket, n, a, b, data.qk, data.excited, sjt, nullptr) / f;
          // (c) of L3
          sum += data.qk.Q(k, bra, n, a, b) *
                 MBPT::L3(k, ket, n, a, b, data.qk, data.core, data.excited,
                          sjt, nullptr) /
                 f;
        }
      }
    }
  }
  return sum;
}

} // namespace UnitTest

//==============================================================================
TEST_CASE("Ladder: mirror terms", "[Ladder][unit]") {
  const auto data = UnitTest::ladder_test_data(-1, 8.0);
  const auto projection = qip::merge(data.core, data.excited);
  const auto max_k = data.qk.max_k();

  const auto Sigmas =
    MBPT::Sigma_ladder_mirror(data.valence, data.core, data.excited, projection,
                              data.qk, true, max_k, 1.0e-4, 60.0, 1, true);
  const auto num_Sigmas = Sigmas.size();
  REQUIRE(num_Sigmas == data.valence.size());

  double worst_eps_diagonal = 0.0;
  double worst_eps_transpose = 0.0;
  std::size_t num_off_diagonal = 0;
  for (std::size_t iv = 0; iv < data.valence.size(); ++iv) {
    const auto &v = data.valence[iv];

    // Diagonal: the energy of the mirror images is that of the terms
    const auto expected = UnitTest::unpaired_ladder_terms(v, v, v.en(), data);
    const auto found = v * (Sigmas[iv] * v);
    worst_eps_diagonal =
      std::max(worst_eps_diagonal, std::abs(found / expected - 1.0));

    // Off-diagonal: the mirror image is the transpose, so
    // <y|Sigma|v> is the sum of the terms with v on the Coulomb vertex and
    // y in the ladder integral
    for (const auto &y : data.excited) {
      if (y.kappa() != v.kappa() || y == v)
        continue;
      const auto expected_yv =
        UnitTest::unpaired_ladder_terms(v, y, v.en(), data);
      const auto found_yv = y * (Sigmas[iv] * v);
      worst_eps_transpose =
        std::max(worst_eps_transpose, std::abs(found_yv / expected_yv - 1.0));
      ++num_off_diagonal;
    }
  }
  REQUIRE(num_off_diagonal > 0);
  REQUIRE(worst_eps_diagonal < 1.0e-4);
  REQUIRE(worst_eps_transpose < 1.0e-4);
}
