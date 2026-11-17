#include "MBPT/Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "MBPT/Ladder.testdata.hpp"
#include "MBPT/LadderIteration.hpp"
#include "MBPT/LadderReference.hpp"
#include "MBPT/Sigma2.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace UnitTest::LadderTest {

//! Pieces (1), (2), (3) of P_mnij, one side, in the m-scheme, for one set of
//! projections (twice m): brute-force sums over r in excited, c in the core
//! and their projections, of +g_cnrj X_mric, -g_cnjr X_mric and
//! -g_cnrj X_mrci, each over e_i + e_c - e_m - e_r. g from the Coulomb table;
//! X = g + x, with x from lk (nullptr: X = g)
inline std::array<double, 3>
P_pieces_mscheme(const DiracSpinor &m, const DiracSpinor &n,
                 const DiracSpinor &i, const DiracSpinor &j, int tmm, int tmn,
                 int tmi, int tmj, const std::vector<DiracSpinor> &excited,
                 const LadderTestData &data, const Coulomb::LkTable *lk) {
  std::array<double, 3> pieces{};
  for (const auto &r : excited) {
    for (const auto &c : data.core) {
      const auto de = i.en() + c.en() - m.en() - r.en();
      for (int tmr = -r.twoj(); tmr <= r.twoj(); tmr += 2) {
        for (int tmc = -c.twoj(); tmc <= c.twoj(); tmc += 2) {
          const auto g_cnrj = data.qk.g(c, n, r, j, tmc, tmn, tmr, tmj);
          const auto g_cnjr = data.qk.g(c, n, j, r, tmc, tmn, tmj, tmr);
          auto X_mric = data.qk.g(m, r, i, c, tmm, tmr, tmi, tmc);
          auto X_mrci = data.qk.g(m, r, c, i, tmm, tmr, tmc, tmi);
          if (lk != nullptr) {
            X_mric += lk->g(m, r, i, c, tmm, tmr, tmi, tmc);
            X_mrci += lk->g(m, r, c, i, tmm, tmr, tmc, tmi);
          }
          pieces[0] += g_cnrj * X_mric / de;
          pieces[1] -= g_cnjr * X_mric / de;
          pieces[2] -= g_cnrj * X_mrci / de;
        }
      }
    }
  }
  return pieces;
}

//! sum_k A^k_mnij T^k for one set of projections (twice m), with the angular
//! factor A^k of CoulombTable::g(); T_k[k - k_min] for k = k_min, k_min + 1,
//! ... (both parities)
inline double multipole_sum(const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &j, int tmm,
                            int tmn, int tmi, int tmj, int k_min,
                            const std::vector<double> &T_k) {
  if (tmi - tmm != tmn - tmj)
    return 0.0;
  const int twoq = tmi - tmm;
  double sum = 0.0;
  for (std::size_t index = 0; index < T_k.size(); ++index) {
    const auto k = k_min + int(index);
    if (std::abs(twoq) > 2 * k)
      continue;
    const auto tjs1 =
      Angular::threej_2(m.twoj(), 2 * k, i.twoj(), -tmm, -twoq, tmi);
    const auto tjs2 =
      Angular::threej_2(n.twoj(), 2 * k, j.twoj(), -tmn, twoq, tmj);
    const auto s = Angular::neg1pow_2(2 * k + tmn - tmm + twoq);
    sum += s * tjs1 * tjs2 * T_k[index];
  }
  return sum;
}

//! Hole-particle interaction of V^(a) on the pair (n,j), one side, in the
//! m-scheme, for one set of projections (twice m): the brute-force sum over r
//! in excited and its projection of -y_{nr} X_mrij / (e_i + e_j - e_m - e_r),
//! with y_{nr} the average of g_njrj over the projections of the hole j (the
//! hole keeps its own projection in X). Zero unless j is a core orbital. X =
//! g + x, with x from lk (nullptr: X = g)
inline double Vhp_mscheme(const DiracSpinor &m, const DiracSpinor &n,
                          const DiracSpinor &i, const DiracSpinor &j, int tmm,
                          int tmn, int tmi, int tmj,
                          const std::vector<DiracSpinor> &excited,
                          const LadderTestData &data,
                          const Coulomb::LkTable *lk) {
  if (std::find(data.core.cbegin(), data.core.cend(), j) == data.core.cend())
    return 0.0;
  double h = 0.0;
  for (const auto &r : excited) {
    const auto de = i.en() + j.en() - m.en() - r.en();
    for (int tmr = -r.twoj(); tmr <= r.twoj(); tmr += 2) {
      double y_nr = 0.0;
      for (int tmc = -j.twoj(); tmc <= j.twoj(); tmc += 2) {
        y_nr += data.qk.g(n, j, r, j, tmn, tmc, tmr, tmc);
      }
      y_nr /= double(j.twojp1());
      auto X_mrij = data.qk.g(m, r, i, j, tmm, tmr, tmi, tmj);
      if (lk != nullptr) {
        X_mrij += lk->g(m, r, i, j, tmm, tmr, tmi, tmj);
      }
      h -= y_nr * X_mrij / de;
    }
  }
  return h;
}

} // namespace UnitTest::LadderTest

//==============================================================================
// Independent (m-scheme) check of the reduced ladder integrals.
// Builds the lowest-order ladder integral l_mnij by brute-force summation
// over magnetic substates, using m-scheme Coulomb integrals g (from the Qk
// table), and compares against the angular-reduced form
//   l_mnij = sum_k A^k_mnij L^k_mnij ,
// with L^k from Lkmnij(). This independently validates the angular reduction:
// phases, 6j structure, denominators, all four diagrams (L1-L4), and the k
// range (including the Coulomb-forbidden-parity k; see k_minmax_L).
TEST_CASE("Ladder: m-scheme check", "[Ladder][unit]") {

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 2;
  const int num_ns = 3;

  // H-like orbitals: pure angular/algebra test - no HF required
  const auto orbs = DiracSpinor::HlikeBasis(lmax, num_ns, radial_grid, zeff);

  // "core" = 1s; "excited" = all others
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fn : orbs) {
    (Fn.n() == 1 ? core : excited).push_back(Fn);
  }
  REQUIRE(core.size() == 1);
  REQUIRE(excited.size() >= 5);

  const Coulomb::YkTable yk(orbs);
  const auto &sjt = yk.SixJ();
  Coulomb::QkTable qk;
  qk.fill(orbs, yk, -1, false);

  const auto find = [&orbs](int n, int kappa) {
    const auto it =
      std::find_if(orbs.cbegin(), orbs.cend(), [n, kappa](const auto &Fn) {
        return Fn.n() == n && Fn.kappa() == kappa;
      });
    REQUIRE(it != orbs.cend());
    return *it;
  };

  // Brute-force m-scheme ladder (lowest order), all four diagrams;
  // internal magnetic substates summed explicitly. cf. Eq. (l3) of the
  // ladder notes: pp + ph + ph + hh
  const auto l_brute = [&](const DiracSpinor &m, const DiracSpinor &n,
                           const DiracSpinor &i, const DiracSpinor &j, int tmm,
                           int tmn, int tmi, int tmj) {
    double l = 0.0;
    // pp ladder:
    for (const auto &r : excited) {
      for (const auto &s : excited) {
        const auto de = i.en() + j.en() - r.en() - s.en();
        for (int tmr = -r.twoj(); tmr <= r.twoj(); tmr += 2) {
          for (int tms = -s.twoj(); tms <= s.twoj(); tms += 2) {
            l += qk.g(m, n, r, s, tmm, tmn, tmr, tms) *
                 qk.g(r, s, i, j, tmr, tms, tmi, tmj) / de;
          }
        }
      }
    }
    // ph (ring) terms:
    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto de1 = c.en() + j.en() - m.en() - r.en();
        const auto de2 = c.en() + i.en() - n.en() - r.en();
        for (int tmr = -r.twoj(); tmr <= r.twoj(); tmr += 2) {
          for (int tmc = -c.twoj(); tmc <= c.twoj(); tmc += 2) {
            l -= qk.g(c, n, i, r, tmc, tmn, tmi, tmr) *
                 qk.g(m, r, c, j, tmm, tmr, tmc, tmj) / de1;
            l -= qk.g(m, c, r, j, tmm, tmc, tmr, tmj) *
                 qk.g(r, n, i, c, tmr, tmn, tmi, tmc) / de2;
          }
        }
      }
    }
    // hh ladder (L4):
    for (const auto &c : core) {
      for (const auto &d : core) {
        const auto de = c.en() + d.en() - m.en() - n.en();
        for (int tmc = -c.twoj(); tmc <= c.twoj(); tmc += 2) {
          for (int tmd = -d.twoj(); tmd <= d.twoj(); tmd += 2) {
            l += qk.g(c, d, i, j, tmc, tmd, tmi, tmj) *
                 qk.g(m, n, c, d, tmm, tmn, tmc, tmd) / de;
          }
        }
      }
    }
    return l;
  };

  // Reduced form: l_mnij = sum_k A^k_mnij L^k, with the same A^k convention
  // as CoulombTable::g() (3j symbols and phase)
  const auto l_reduced = [](const DiracSpinor &m, const DiracSpinor &n,
                            const DiracSpinor &i, const DiracSpinor &j,
                            const std::vector<double> &Lks, int k0, int tmm,
                            int tmn, int tmi, int tmj) {
    if (tmi - tmm != tmn - tmj)
      return 0.0;
    const int twoq = tmi - tmm;
    double l = 0.0;
    for (std::size_t ik = 0; ik < Lks.size(); ++ik) {
      const auto k = k0 + int(ik);
      if (std::abs(twoq) > 2 * k)
        continue;
      const auto tjs1 =
        Angular::threej_2(m.twoj(), 2 * k, i.twoj(), -tmm, -twoq, tmi);
      const auto tjs2 =
        Angular::threej_2(n.twoj(), 2 * k, j.twoj(), -tmn, twoq, tmj);
      const auto s = Angular::neg1pow_2(2 * k + tmn - tmm + twoq);
      l += s * tjs1 * tjs2 * Lks[ik];
    }
    return l;
  };

  // Test tuples {m, n, i, j}: chosen to cover both parities of k (including
  // the Coulomb-forbidden k=1 in all-s channels), ph-dominated channels, and
  // core-core (i,j) pairs
  const auto s1 = find(1, -1), s2 = find(2, -1), s3 = find(3, -1);
  const auto p2 = find(2, 1), P2 = find(2, -2), p3 = find(3, 1);
  const auto d3 = find(3, 2);
  using Four = std::array<DiracSpinor, 4>;
  const std::vector<Four> tuples{
    {s3, s2, s2, s1}, {p2, P2, s2, s1}, {d3, d3, s1, s1}, {p3, P2, s2, s1}};

  fmt::print("{:16s} {:>13s} {:>13s} {:>9s}\n", "l_mnij", "expected", "found",
             "eps");
  for (const auto &[m, n, i, j] : tuples) {

    // Reduced L^k over the full ladder k range (both parities)
    const auto [k0, kI] = MBPT::k_minmax_L(m, n, i, j);
    REQUIRE(kI >= k0);
    std::vector<double> Lks;
    for (int k = k0; k <= kI; ++k) {
      Lks.push_back(
        MBPT::Lkmnij(k, m, n, i, j, qk, core, excited, true, sjt, nullptr));
    }
    // Coulomb-forbidden-parity components must be present (parity fix)
    const auto [q0, qI] = Coulomb::k_minmax_Q(m, n, i, j);
    for (int k = k0; k <= kI; ++k) {
      const bool Q_allowed = k >= q0 && k <= qI && (k - q0) % 2 == 0;
      if (!Q_allowed) {
        REQUIRE(Lks[std::size_t(k - k0)] != 0.0);
      }
    }

    // Compare brute-force and reduced, for every magnetic substate combo
    double worst_del = 0.0, scale = 0.0;
    double x_expect = 0.0, x_found = 0.0;
    for (int tmm = -m.twoj(); tmm <= m.twoj(); tmm += 2) {
      for (int tmn = -n.twoj(); tmn <= n.twoj(); tmn += 2) {
        for (int tmi = -i.twoj(); tmi <= i.twoj(); tmi += 2) {
          for (int tmj = -j.twoj(); tmj <= j.twoj(); tmj += 2) {
            const auto lb = l_brute(m, n, i, j, tmm, tmn, tmi, tmj);
            const auto lr = l_reduced(m, n, i, j, Lks, k0, tmm, tmn, tmi, tmj);
            worst_del = std::max(worst_del, std::abs(lb - lr));
            if (std::abs(lb) > scale) {
              scale = std::abs(lb);
              x_expect = lb;
              x_found = lr;
            }
          }
        }
      }
    }
    REQUIRE(scale > 0.0);
    const auto eps = worst_del / scale;
    const auto label = "L_" + m.shortSymbol() + n.shortSymbol() +
                       i.shortSymbol() + j.shortSymbol();
    fmt::print("{:16s} {:13.6e} {:13.6e} {:9.1e}\n", label, x_expect, x_found,
               eps);
    REQUIRE(eps < 1.0e-10);
  }
}

//==============================================================================
// Independent (m-scheme) check of the reduced pieces of P [P1(), P2(), P3(),
// Pkmnij()]: each one-sided piece, of (m,n,i,j) and of its partner
// (n,m,j,i), and the whole P_mnij, summed by brute force over the projections
// of r and c, against
//   P_mnij = sum_k A^k_mnij P^k_mnij
// for every set of projections of m, n, i, j. At lowest order (X = g), and
// with X = g + x, x random at both parities of k.
TEST_CASE("Ladder: P m-scheme check", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  const auto x_random = random_ladder_table(data);

  // r runs over one excited orbital of each kappa (s to f): keeps the
  // brute-force sums quick. The reduced pieces use the same set
  std::vector<DiracSpinor> excited;
  for (const auto &group : MBPT::group_by_kappa(data.excited)) {
    excited.push_back(data.excited[group.front()]);
  }

  // {m, n, i, j}: i core (2p+), and i valence (3p+); k of both parities
  const auto *d5_3 = DiracSpinor::find(3, -3, data.excited);
  const auto *d3_3 = DiracSpinor::find(3, 2, data.excited);
  const auto *d5_4 = DiracSpinor::find(4, -3, data.excited);
  const auto *d3_4 = DiracSpinor::find(4, 2, data.excited);
  const auto *p3_2 = DiracSpinor::find(2, -2, data.core);
  const auto *p1_2 = DiracSpinor::find(2, 1, data.core);
  const auto *p3_3 = DiracSpinor::find(3, -2, data.valence);
  const std::vector<const DiracSpinor *> found{d5_3, d3_3, d5_4, d3_4,
                                               p3_2, p1_2, p3_3};
  const auto num_missing = std::count(found.cbegin(), found.cend(), nullptr);
  REQUIRE(num_missing == 0);
  using Four = std::array<DiracSpinor, 4>;
  const std::vector<Four> tuples{{*d5_3, *d3_4, *p3_2, *p1_2},
                                 {*d3_3, *d5_4, *p3_3, *p1_2}};

  // Worst relative difference over each quantity (3 pieces of each side, and
  // the whole P) of each tuple; and the largest P^k at l_m + l_i + k odd
  double worst_eps = 0.0;
  std::size_t num_zero_quantities = 0;
  double largest_P_odd_lowest_order = 0.0;
  double largest_P2_odd_dressed = 0.0;
  const std::vector<const Coulomb::LkTable *> dressing{nullptr, &x_random};
  for (const auto *x : dressing) {
    for (const auto &[m, n, i, j] : tuples) {
      // Reduced: pieces of (m,n,i,j) and of (n,m,j,i), and the whole P
      const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, i, j);
      std::array<std::vector<double>, 3> pieces, partner_pieces;
      std::vector<double> P_k;
      for (int k = k_min; k <= k_max; ++k) {
        pieces[0].push_back(
          MBPT::P1(k, m, n, i, j, data.qk, data.core, excited, x));
        pieces[1].push_back(
          MBPT::P2(k, m, n, i, j, data.qk, data.core, excited, sjt, x));
        pieces[2].push_back(
          MBPT::P3(k, m, n, i, j, data.qk, data.core, excited, sjt, x));
        partner_pieces[0].push_back(
          MBPT::P1(k, n, m, j, i, data.qk, data.core, excited, x));
        partner_pieces[1].push_back(
          MBPT::P2(k, n, m, j, i, data.qk, data.core, excited, sjt, x));
        partner_pieces[2].push_back(
          MBPT::P3(k, n, m, j, i, data.qk, data.core, excited, sjt, x));
        P_k.push_back(
          MBPT::Pkmnij(k, m, n, i, j, data.qk, data.core, excited, sjt, x));
        if ((m.l() + i.l() + k) % 2 == 0)
          continue;
        if (x == nullptr) {
          largest_P_odd_lowest_order =
            std::max(largest_P_odd_lowest_order, std::abs(P_k.back()));
        } else {
          largest_P2_odd_dressed =
            std::max(largest_P2_odd_dressed, std::abs(pieces[1].back()));
        }
      }

      // Brute force, for every set of projections
      std::array<double, 7> max_difference{}, largest{};
      for (int tmm = -m.twoj(); tmm <= m.twoj(); tmm += 2) {
        for (int tmn = -n.twoj(); tmn <= n.twoj(); tmn += 2) {
          for (int tmi = -i.twoj(); tmi <= i.twoj(); tmi += 2) {
            for (int tmj = -j.twoj(); tmj <= j.twoj(); tmj += 2) {
              const auto brute = P_pieces_mscheme(m, n, i, j, tmm, tmn, tmi,
                                                  tmj, excited, data, x);
              const auto brute_partner = P_pieces_mscheme(
                n, m, j, i, tmn, tmm, tmj, tmi, excited, data, x);
              std::array<double, 7> expected{}, reduced{};
              for (std::size_t p = 0; p < 3; ++p) {
                expected[p] = brute[p];
                reduced[p] = multipole_sum(m, n, i, j, tmm, tmn, tmi, tmj,
                                           k_min, pieces[p]);
                expected[3 + p] = brute_partner[p];
                reduced[3 + p] = multipole_sum(n, m, j, i, tmn, tmm, tmj, tmi,
                                               k_min, partner_pieces[p]);
                expected[6] += brute[p] + brute_partner[p];
              }
              reduced[6] =
                multipole_sum(m, n, i, j, tmm, tmn, tmi, tmj, k_min, P_k);
              for (std::size_t q = 0; q < expected.size(); ++q) {
                max_difference[q] = std::max(
                  max_difference[q], std::abs(expected[q] - reduced[q]));
                largest[q] = std::max(largest[q], std::abs(expected[q]));
              }
            }
          }
        }
      }
      for (std::size_t q = 0; q < largest.size(); ++q) {
        if (largest[q] == 0.0) {
          ++num_zero_quantities;
        } else {
          worst_eps = std::max(worst_eps, max_difference[q] / largest[q]);
        }
      }
    }
  }
  REQUIRE(num_zero_quantities == 0);
  REQUIRE(worst_eps < 1.0e-10);
  // At X = g all of P has the Coulomb parity of k; with x, piece (2) also
  // has the other one (through x^k_mric), so the check covers it
  REQUIRE(largest_P_odd_lowest_order == 0.0);
  REQUIRE(largest_P2_odd_dressed > 1.0e-6);
}

//==============================================================================
// Independent (m-scheme) check of the hole-particle interaction of V^(a)
// [Vhp()]: -y^0_jj on n (one side) and its partner -y^0_ii on m (i core),
// with y the average of g over the projections of the hole, summed by brute
// force over r and its projection, against
//   H_mnij = sum_k A^k_mnij H^k_mnij
// for every set of projections of m, n, i, j. At lowest order (X = g), and
// with X = g + x, x random at both parities of k.
TEST_CASE("Ladder: Vhp m-scheme check", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto x_random = random_ladder_table(data);

  // {m, n, i, j}: i core; i valence (no partner); m and n of one kappa
  const auto *d5_3 = DiracSpinor::find(3, -3, data.excited);
  const auto *d3_3 = DiracSpinor::find(3, 2, data.excited);
  const auto *d5_4 = DiracSpinor::find(4, -3, data.excited);
  const auto *d3_4 = DiracSpinor::find(4, 2, data.excited);
  const auto *p3_3 = DiracSpinor::find(3, -2, data.excited);
  const auto *p3_4 = DiracSpinor::find(4, -2, data.excited);
  const auto *s_1 = DiracSpinor::find(1, -1, data.core);
  const auto *s_2 = DiracSpinor::find(2, -1, data.core);
  const auto *p3_2 = DiracSpinor::find(2, -2, data.core);
  const auto *p1_2 = DiracSpinor::find(2, 1, data.core);
  const auto *p3_3_valence = DiracSpinor::find(3, -2, data.valence);
  const std::vector<const DiracSpinor *> found{
    d5_3, d3_3, d5_4, d3_4, p3_3, p3_4, s_1, s_2, p3_2, p1_2, p3_3_valence};
  const auto num_missing = std::count(found.cbegin(), found.cend(), nullptr);
  REQUIRE(num_missing == 0);
  using Four = std::array<DiracSpinor, 4>;
  const std::vector<Four> tuples{{*d5_3, *d3_4, *p3_2, *p1_2},
                                 {*d3_3, *d5_4, *p3_3_valence, *p1_2},
                                 {*p3_3, *p3_4, *s_1, *s_2}};

  // Worst relative difference over each side of each tuple; the largest
  // partner with i valence; the largest H^k at l_m + l_i + k odd
  double worst_eps = 0.0;
  std::size_t num_zero_quantities = 0;
  double largest_partner_i_valence = 0.0;
  double largest_odd_lowest_order = 0.0;
  double largest_odd_dressed = 0.0;
  const std::vector<const Coulomb::LkTable *> dressing{nullptr, &x_random};
  for (const auto *x : dressing) {
    for (const auto &[m, n, i, j] : tuples) {
      const bool i_core =
        std::find(data.core.cbegin(), data.core.cend(), i) != data.core.cend();
      // Reduced: H^k_mnij (on (n,j)) and the partner H^k_nmji (on (m,i))
      const auto [k_min, k_max] = MBPT::k_minmax_L(m, n, i, j);
      std::vector<double> side_k, partner_k;
      for (int k = k_min; k <= k_max; ++k) {
        side_k.push_back(
          MBPT::Vhp(k, m, n, i, j, data.qk, data.core, data.excited, x));
        partner_k.push_back(
          MBPT::Vhp(k, n, m, j, i, data.qk, data.core, data.excited, x));
        if (!i_core) {
          largest_partner_i_valence =
            std::max(largest_partner_i_valence, std::abs(partner_k.back()));
        }
        if ((m.l() + i.l() + k) % 2 == 0)
          continue;
        const auto largest_odd =
          std::max(std::abs(side_k.back()), std::abs(partner_k.back()));
        if (x == nullptr) {
          largest_odd_lowest_order =
            std::max(largest_odd_lowest_order, largest_odd);
        } else {
          largest_odd_dressed = std::max(largest_odd_dressed, largest_odd);
        }
      }

      // Brute force, for every set of projections
      std::array<double, 2> max_difference{}, largest{};
      for (int tmm = -m.twoj(); tmm <= m.twoj(); tmm += 2) {
        for (int tmn = -n.twoj(); tmn <= n.twoj(); tmn += 2) {
          for (int tmi = -i.twoj(); tmi <= i.twoj(); tmi += 2) {
            for (int tmj = -j.twoj(); tmj <= j.twoj(); tmj += 2) {
              const std::array<double, 2> expected{
                Vhp_mscheme(m, n, i, j, tmm, tmn, tmi, tmj, data.excited, data,
                            x),
                Vhp_mscheme(n, m, j, i, tmn, tmm, tmj, tmi, data.excited, data,
                            x)};
              const std::array<double, 2> reduced{
                multipole_sum(m, n, i, j, tmm, tmn, tmi, tmj, k_min, side_k),
                multipole_sum(n, m, j, i, tmn, tmm, tmj, tmi, k_min,
                              partner_k)};
              for (std::size_t q = 0; q < 2; ++q) {
                max_difference[q] = std::max(
                  max_difference[q], std::abs(expected[q] - reduced[q]));
                largest[q] = std::max(largest[q], std::abs(expected[q]));
              }
            }
          }
        }
      }
      // The partner exists only with i core
      const std::size_t num_quantities = i_core ? 2 : 1;
      for (std::size_t q = 0; q < num_quantities; ++q) {
        if (largest[q] == 0.0) {
          ++num_zero_quantities;
        } else {
          worst_eps = std::max(worst_eps, max_difference[q] / largest[q]);
        }
      }
    }
  }
  REQUIRE(num_zero_quantities == 0);
  REQUIRE(worst_eps < 1.0e-10);
  REQUIRE(largest_partner_i_valence == 0.0);
  // At X = g, H has the Coulomb parity of k; with x, also the other one
  REQUIRE(largest_odd_lowest_order == 0.0);
  REQUIRE(largest_odd_dressed > 1.0e-6);
}

//==============================================================================
TEST_CASE("Ladder: core energy with W", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto &sjt = data.yk.SixJ();
  // Ladder table with entries at both parities of k: lowest order, one sweep
  Coulomb::LkTable lk;
  MBPT::fill_Lk_mnib(&lk, data.qk, data.excited, data.core, data.core,
                     {true, true, true, true}, -1, false);
  const auto lk_lowest = lk;
  MBPT::update_Lk_mnib(&lk, data.qk, data.excited, data.core, data.core,
                       {true, true, true, true}, &lk_lowest, 0.0, false);

  const auto expected = MBPT::de_core(data.qk, lk, data.core, data.excited);
  const auto found =
    MBPT::de_core_w(data.qk, lk, data.core, data.excited, &sjt);
  const auto eps = std::abs(found / expected - 1.0);
  REQUIRE(eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: direct closing at second order", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();

  // With the Coulomb integrals in place of x_F and of the closing amplitude:
  // diagrams (a) and (c), the direct part of the second-order energy
  double max_eps = 0.0;
  for (const auto &v : data.valence) {
    const auto de2_direct =
      MBPT::Sigma_vw_direct_exchange(v, v, data.qk, data.core, data.excited)
        .first;
    const auto de2_closing =
      MBPT::de_valence_F(v, data.qk, data.qk, data.qk, data.core, data.excited);
    max_eps = std::max(max_eps, std::abs(de2_closing / de2_direct - 1.0));
  }
  REQUIRE(max_eps < 1.0e-12);
}

//==============================================================================
TEST_CASE("Ladder: direct closing of x_F", "[Ladder][unit]") {
  using namespace UnitTest::LadderTest;
  const auto data = ladder_test_data();
  const auto max_k = stored_max_k(data);

  // Lowest order of x_F and of its closing amplitude: the polarisation loop
  // alone, and with the hole-particle interaction
  MBPT::LadderDiagrams P1_only{false, false, false, false};
  P1_only.P1 = true;
  auto P1_Vhp = P1_only;
  P1_Vhp.hole_particle = true;
  double max_eps_loop = 0.0;
  double smallest_eps_Vhp = 1.0;
  for (const auto *diagrams : {&P1_only, &P1_Vhp}) {
    Coulomb::LkTable lk_F;
    MBPT::add_Lk_mnib_entries(&lk_F, data.excited, data.core, data.i_orbs,
                              max_k);
    MBPT::calculate_Lk_mnib(&lk_F, data.qk, data.excited, data.core,
                            data.i_orbs, *diagrams, nullptr, max_k, false, 0.0);
    Coulomb::NkTable L_vnab;
    MBPT::add_Lk_vnab_F_entries(&L_vnab, data.excited, data.core, data.valence,
                                max_k);
    MBPT::calculate_Lk_vnab_F(&L_vnab, data.qk, data.excited, data.core,
                              data.valence, *diagrams, nullptr, nullptr, 0.0);
    for (const auto &v : data.valence) {
      // Closed with the core amplitude x_F_vnab in place of the closing
      // amplitude: the same without the hole-particle interaction; with it,
      // they differ by its action on (v,a)
      const auto de_F =
        MBPT::de_valence_F(v, data.qk, lk_F, L_vnab, data.core, data.excited);
      const auto de_F_core_table =
        MBPT::de_valence_F(v, data.qk, lk_F, lk_F, data.core, data.excited);
      const auto eps = std::abs(de_F_core_table / de_F - 1.0);
      if (diagrams->hole_particle) {
        smallest_eps_Vhp = std::min(smallest_eps_Vhp, eps);
      } else {
        max_eps_loop = std::max(max_eps_loop, eps);
      }
    }
  }
  REQUIRE(max_eps_loop < 1.0e-12);
  REQUIRE(smallest_eps_Vhp > 1.0e-4);
}
