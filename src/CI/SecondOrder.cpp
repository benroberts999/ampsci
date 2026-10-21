#include "SecondOrder.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include "DiracOperator/include.hpp"
#include "MixedStates.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include <algorithm>
#include <cmath>

namespace CI {

//==============================================================================
std::pair<double, double> A_K_coefs(int K, int kt, int ks, int twoJb, int twoJn,
                                    int twoJa) {

  const auto f = std::sqrt(2.0 * K + 1.0) * Angular::neg1pow_2(twoJb + twoJa);

  const auto c1 = Angular::neg1pow(K) * f *
                  Angular::sixj_2(2 * K, 2 * ks, 2 * kt, twoJn, twoJb, twoJa);

  const auto c2 = Angular::neg1pow(kt + ks) * f *
                  Angular::sixj_2(2 * K, 2 * kt, 2 * ks, twoJn, twoJb, twoJa);

  return {c1, c2};
}

//==============================================================================
double z_component(int K, int kt, int ks, int twoJb, int twoJa, int two_m) {
  return Angular::cg_2(2 * kt, 0, 2 * ks, 0, 2 * K, 0) *
         Angular::neg1pow_2(twoJb - two_m) *
         Angular::threej_2(twoJb, 2 * K, twoJa, -two_m, 0, two_m);
}

//==============================================================================
int symm_sign(const DiracOperator::TensorOperator *h, int twoJA, int twoJB) {
  const auto s_imaginary = h->imaginaryQ() ? -1 : 1;
  return s_imaginary * Angular::neg1pow_2(twoJA - twoJB);
}

//==============================================================================
double sigma_rme(const PsiJPi &Psi_b, std::size_t ib, const PsiJPi &Psi_a,
                 std::size_t ia, const std::vector<DiracSpinor> &ci_basis) {

  // sigma = 2S, a rank 1, even parity, operator
  const DiracOperator::s spin{};
  const auto spin_me = Amplitudes::me_table(ci_basis, &spin);

  return 2.0 * ReducedME(Psi_b, ib, Psi_a, ia, spin_me, 1, 1);
}

//==============================================================================
std::pair<double, double>
A_K(int K, const PsiJPi &Psi_b, std::size_t ib, const PsiJPi &Psi_a,
    std::size_t ia, const DiracOperator::TensorOperator *t,
    const Coulomb::meTable<double> &t_me,
    const DiracOperator::TensorOperator *s,
    const Coulomb::meTable<double> &s_me, double omega, double omega_s,
    const Integrals &ints, const std::vector<Level> &levels_to_remove,
    const std::vector<PsiJPi> &solved_Jpi, std::ostream &outstream) {

  // The amplitude, evaluated with the mixed states of s, and of t
  double A_s{0.0};
  double A_t{0.0};

  const auto kt = t->rank();
  const auto ks = s->rank();
  const auto twoJa = Psi_a.twoJ();
  const auto twoJb = Psi_b.twoJ();
  const auto Ea = Psi_a.energy(ia);
  const auto Eb = Psi_b.energy(ib);

  // Each row below requires the CI Hamiltonian of that (J, parity), and a
  // mixed state for each term: both are slow, so say so before starting
  fmt::print(outstream,
             "\nSolving the mixed states for each required J^pi..\n\n");
  fmt::print(outstream, "{:>5} {:>7} {:>4}  {:>16} {:>16} {:>8}\n", "J^pi",
             "CSFs", "term", "<B_s||t||A_s>", "<B_t||s||A_t>", "eps");
  outstream << std::flush;

  // Intermediate states are connected to a by one operator, and to b by the
  // other. The first ('ts') term has s acting on a, so those states have parity
  // pi_a*pi_s; the second ('st') term has t acting on a: pi_a*pi_t
  const auto twoJn_max = std::max(twoJa, twoJb) + 2 * std::max(kt, ks);
  for (int twoJn = twoJa % 2; twoJn <= twoJn_max; twoJn += 2) {
    for (const auto pi_n : {-1, 1}) {

      const bool do_ts = pi_n == Psi_a.parity() * s->parity() &&
                         Angular::triangle(twoJb, 2 * kt, twoJn) != 0 &&
                         Angular::triangle(twoJn, 2 * ks, twoJa) != 0;

      const bool do_st = pi_n == Psi_a.parity() * t->parity() &&
                         Angular::triangle(twoJb, 2 * ks, twoJn) != 0 &&
                         Angular::triangle(twoJn, 2 * kt, twoJa) != 0;

      if (!do_ts && !do_st)
        continue;

      // If this J/Pi was solved in the CI run, start from that solution: the
      // energy-dependent corrections to the Hamiltonian are then at its own
      // reference energy, as in the run (and as in the levels a
      // sum-over-states would use), rather than at the ground-state energy
      PsiJPi target(twoJn, pi_n, ints.ci_basis, ints.num_electrons,
                    ints.config_space);
      for (const auto &psi : solved_Jpi) {
        if (psi.twoJ() == twoJn && psi.parity() == pi_n) {
          target = psi;
          break;
        }
      }
      if (target.num_CSFs() == 0)
        continue;
      const auto Hci = construct_Hci(target, ints);

      // Levels of this (J, parity) that are to be removed from the mixed
      // states. The CI problem must be solved (as far as those levels) first
      std::vector<std::size_t> indices;
      for (const auto &level : levels_to_remove) {
        if (level.twoJ == twoJn && level.parity == pi_n) {
          indices.push_back(level.index);
        }
      }
      if (!indices.empty()) {
        const auto num_solutions =
          *std::max_element(indices.cbegin(), indices.cend()) + 1;
        target.solve(Hci, int(num_solutions));
      }

      // The two sums over the intermediate states of this (J, parity), each
      // formed with the mixed states of s, and with those of t
      double ts_s{0.0}, ts_t{0.0}, st_s{0.0}, st_t{0.0};

      if (do_ts) {
        // sum_n <b||t||n><n||s||a>/(E_a + omega_s - E_n)
        const auto da = project_out(
          solve_mixed_state(Psi_a, ia, target, Hci, s_me, ks, omega_s), target,
          indices);
        ts_s = ReducedME(Psi_b, ib, da, 0, t_me, kt, t->parity());

        const auto Db =
          project_out(solve_mixed_state(Psi_b, ib, target, Hci, t_me, kt,
                                        Ea + omega_s - Eb),
                      target, indices);
        ts_t = symm_sign(t, twoJb, twoJn) *
               ReducedME(Db, 0, Psi_a, ia, s_me, ks, s->parity());
      }

      if (do_st) {
        // sum_n <b||s||n><n||t||a>/(E_a + omega - E_n)
        const auto Da = project_out(
          solve_mixed_state(Psi_a, ia, target, Hci, t_me, kt, omega), target,
          indices);
        st_t = ReducedME(Psi_b, ib, Da, 0, s_me, ks, s->parity());

        const auto db = project_out(
          solve_mixed_state(Psi_b, ib, target, Hci, s_me, ks, Ea + omega - Eb),
          target, indices);
        st_s = symm_sign(s, twoJb, twoJn) *
               ReducedME(db, 0, Psi_a, ia, t_me, kt, t->parity());
      }

      // Contributions to A^K. With the mixed states of s, the 'ts' term is
      // formed from a and the 'st' term from b; with those of t, the other
      // way around
      const auto [c1, c2] = A_K_coefs(K, kt, ks, twoJb, twoJn, twoJa);
      A_s += c1 * ts_s + c2 * st_s;
      A_t += c1 * ts_t + c2 * st_t;

      const auto Jpi = J_string(twoJn) + (pi_n == 1 ? "+" : "-");
      if (do_ts) {
        fmt::print(outstream, "{:>5} {:>7} {:>4}  {:16.6e} {:16.6e} {:8.1e}\n",
                   Jpi, target.num_CSFs(), "ts", c1 * ts_s, c1 * ts_t,
                   std::abs(ts_s - ts_t) / std::max(std::abs(ts_s), 1.0e-30));
      }
      if (do_st) {
        fmt::print(outstream, "{:>5} {:>7} {:>4}  {:16.6e} {:16.6e} {:8.1e}\n",
                   Jpi, target.num_CSFs(), "st", c2 * st_s, c2 * st_t,
                   std::abs(st_s - st_t) / std::max(std::abs(st_s), 1.0e-30));
      }
      outstream << std::flush;
    }
  }

  fmt::print(outstream, "{:>5} {:>7} {:>4}  {:16.6e} {:16.6e} {:8.1e}\n", "",
             "", fmt::format("A^{}", K), A_s, A_t,
             std::abs(A_s - A_t) / std::max(std::abs(A_s), 1.0e-30));

  return {A_s, A_t};
}

//==============================================================================
double A_K_core(int K, int twoJ, const DiracOperator::TensorOperator *t,
                const DiracOperator::TensorOperator *s, double omega,
                double omega_s, const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &excited,
                const ExternalField::CorePolarisation *dVt,
                const ExternalField::CorePolarisation *dVs) {

  // The core is a closed shell: only a scalar, even-parity, amplitude survives
  if (K != 0 || t->parity() * s->parity() != 1)
    return 0.0;

  if (core.empty() || excited.empty())
    return 0.0;

  const auto kt = t->rank();
  const auto ks = s->rank();

  // Of each pair of matrix elements, only the one acting on the core orbital
  // carries RPA (see the note on A_K_core), at the frequency of its operator
  const auto t_0 = Amplitudes::me_table(core, excited, t);
  const auto t_v = Amplitudes::me_table(core, excited, t, dVt, nullptr, omega);
  const auto s_0 = Amplitudes::me_table(core, excited, s);
  const auto s_v =
    Amplitudes::me_table(core, excited, s, dVs, nullptr, omega_s);

  double A_core = 0.0;
  for (const auto &c : core) {
    double Uc = 0.0;
    for (const auto &m : excited) {
      const auto [c1, c2] = A_K_coefs(0, kt, ks, c.twoj(), m.twoj(), c.twoj());
      // 'ts': s excites the core electron, t returns it
      if (c1 != 0.0) {
        Uc +=
          c1 * t_0.getv(c, m) * s_v.getv(m, c) / (c.en() + omega_s - m.en());
      }
      // 'st': t excites the core electron, s returns it
      if (c2 != 0.0) {
        Uc += c2 * s_0.getv(c, m) * t_v.getv(m, c) / (c.en() + omega - m.en());
      }
    }
    A_core += std::sqrt(double(c.twojp1())) * Uc;
  }

  return std::sqrt(double(twoJ + 1)) * A_core;
}

//==============================================================================
double A_K_cv(int K, const PsiJPi &Psi_b, std::size_t ib, const PsiJPi &Psi_a,
              std::size_t ia, const DiracOperator::TensorOperator *t,
              const DiracOperator::TensorOperator *s, double omega,
              double omega_s, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &ci_basis,
              const ExternalField::CorePolarisation *dVt,
              const ExternalField::CorePolarisation *dVs) {

  const auto pi_ts = t->parity() * s->parity();

  if (core.empty() || ci_basis.empty())
    return 0.0;
  if (Psi_b.parity() * Psi_a.parity() != pi_ts)
    return 0.0;
  if (Angular::triangle(Psi_b.twoJ(), 2 * K, Psi_a.twoJ()) == 0)
    return 0.0;

  const auto kt = t->rank();
  const auto ks = s->rank();

  // Matrix elements between the CI basis and the core. RPA on both,
  // each at the frequency of its operator
  const auto t_vc =
    Amplitudes::me_table(ci_basis, core, t, dVt, nullptr, omega);
  const auto s_vc =
    Amplitudes::me_table(ci_basis, core, s, dVs, nullptr, omega_s);

  // The effective one-body operator of rank K, in the valence space
  Coulomb::meTable<double> W;
  for (const auto &vp : ci_basis) {
    for (const auto &v : ci_basis) {

      if (vp.parity() * v.parity() != pi_ts)
        continue;
      if (Angular::triangle(vp.twoj(), 2 * K, v.twoj()) == 0)
        continue;

      double Wvv = 0.0;
      for (const auto &c : core) {
        const auto [c1, c2] =
          A_K_coefs(K, kt, ks, vp.twoj(), c.twoj(), v.twoj());
        // 'ts': s excites the core electron to v', t drops v into the hole
        if (c2 != 0.0) {
          Wvv += c2 * s_vc.getv(vp, c) * t_vc.getv(c, v) /
                 (vp.en() - omega_s - c.en());
        }
        // 'st': t excites the core electron to v', s drops v into the hole
        if (c1 != 0.0) {
          Wvv += c1 * t_vc.getv(vp, c) * s_vc.getv(c, v) /
                 (vp.en() - omega - c.en());
        }
      }
      W.add(vp, v, Wvv);
    }
  }

  return ReducedME(Psi_b, ib, Psi_a, ia, W, K, pi_ts);
}

} // namespace CI
