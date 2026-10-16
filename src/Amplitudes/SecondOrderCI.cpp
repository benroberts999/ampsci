#include "Amplitudes/SecondOrderCI.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "CI/SecondOrder.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include <cmath>
#include <string>

namespace Amplitudes {

//==============================================================================
double sos_ci(int K, const CI::PsiJPi &Psi_b, std::size_t ib,
              const CI::PsiJPi &Psi_a, std::size_t ia,
              const DiracOperator::TensorOperator *t,
              const Coulomb::meTable<double> &t_me,
              const DiracOperator::TensorOperator *s,
              const Coulomb::meTable<double> &s_me, double omega,
              double omega_s, const std::vector<CI::PsiJPi> &ciwfs,
              const std::vector<CI::Level> &levels_to_remove,
              std::ostream &outstream) {

  double A{0.0};

  const auto kt = t->rank();
  const auto ks = s->rank();
  const auto twoJa = Psi_a.twoJ();
  const auto twoJb = Psi_b.twoJ();
  const auto Ea = Psi_a.energy(ia);

  fmt::print(outstream, "\nSumming over the solved CI levels of each required "
                        "J^pi..\n\n");
  fmt::print(outstream, "{:>5} {:>7}  {:>16} {:>16}\n", "J^pi", "levels",
             "dA(ts)", "dA(st)");

  // Intermediate states are connected to a by one operator, and to b by the
  // other. The first ('ts') term has s acting on a, so those states have
  // parity pi_a*pi_s; the second ('st') term has t acting on a: pi_a*pi_t
  const auto twoJn_max = std::max(twoJa, twoJb) + 2 * std::max(kt, ks);
  for (int twoJn = twoJa % 2; twoJn <= twoJn_max; twoJn += 2) {
    for (const auto pi_n : {-1, 1}) {

      const bool do_ts = pi_n == Psi_a.parity() * s->parity() &&
                         Angular::triangle(twoJb, 2 * kt, twoJn) != 0 &&
                         Angular::triangle(twoJn, 2 * ks, twoJa) != 0;

      const bool do_st = pi_n == Psi_a.parity() * t->parity() &&
                         Angular::triangle(twoJb, 2 * ks, twoJn) != 0 &&
                         Angular::triangle(twoJn, 2 * kt, twoJa) != 0;

      if (!do_ts && !do_st) {
        continue;
      }

      // The solved levels of this (J, parity): only what the CI{} block
      // solved is available, so the sum is truncated to those
      const auto Jpi = CI::J_string(twoJn) + (pi_n == 1 ? "+" : "-");
      const CI::PsiJPi *Psi_n = nullptr;
      for (const auto &psi : ciwfs) {
        if (psi.twoJ() == twoJn && psi.parity() == pi_n) {
          Psi_n = &psi;
        }
      }
      if (Psi_n == nullptr || Psi_n->num_solutions() == 0) {
        fmt::print(outstream, "{:>5}  no CI solutions: contribution missing\n",
                   Jpi);
        continue;
      }

      const auto [c1, c2] = CI::A_K_coefs(K, kt, ks, twoJb, twoJn, twoJa);

      // Sums over the levels of this (J, parity)
      double A_ts{0.0};
      double A_st{0.0};
      std::size_t n_levels = 0;

      for (std::size_t n = 0; n < Psi_n->num_solutions(); ++n) {

        bool removed = false;
        for (const auto &level : levels_to_remove) {
          if (level.twoJ == twoJn && level.parity == pi_n && level.index == n) {
            removed = true;
          }
        }
        if (removed) {
          continue;
        }
        ++n_levels;

        const auto En = Psi_n->energy(n);

        // 'ts' term: c1 <b||t||n><n||s||a> / (E_a + omega_s - E_n)
        if (do_ts && c1 != 0.0) {
          const auto denom = Ea + omega_s - En;
          A_ts +=
            c1 * CI::ReducedME(Psi_b, ib, *Psi_n, n, t_me, kt, t->parity()) *
            CI::ReducedME(*Psi_n, n, Psi_a, ia, s_me, ks, s->parity()) / denom;
        }

        // 'st' term: c2 <b||s||n><n||t||a> / (E_a + omega - E_n)
        if (do_st && c2 != 0.0) {
          const auto denom = Ea + omega - En;
          A_st +=
            c2 * CI::ReducedME(Psi_b, ib, *Psi_n, n, s_me, ks, s->parity()) *
            CI::ReducedME(*Psi_n, n, Psi_a, ia, t_me, kt, t->parity()) / denom;
        }
      }

      fmt::print(outstream, "{:>5} {:>7}  {:16.6e} {:16.6e}\n", Jpi, n_levels,
                 A_ts, A_st);
      outstream << std::flush;
      A += A_ts + A_st;
    }
  }

  fmt::print(outstream, "{:>5} {:>7}  {:>16} {:16.6e}\n", "", "",
             fmt::format("A^{}", K), A);

  return A;
}

} // namespace Amplitudes
