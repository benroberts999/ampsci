#include "MBPT/Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <vector>

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
