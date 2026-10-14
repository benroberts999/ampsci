#include "MixedStates.hpp"
#include "CI_Integrals.hpp"
#include <cassert>

namespace CI {

//==============================================================================
// Forms T_I = sum_K <I||T||K> c_K, the action of the one-body operator T on
// the CI state Psi_0, in the CSF basis of the target J/Pi
LinAlg::Vector<double> TPsi_reduced(const std::vector<CSF2> &CSFs, int twoJ,
                                    const PsiJPi &Psi0, std::size_t i0,
                                    const Coulomb::meTable<double> &h,
                                    int K_rank) {

  const auto &CSF0s = Psi0.CSFs();
  const auto twoJ0 = Psi0.twoJ();
  const auto c0 = Psi0.coefs(i0);

  const auto N = CSFs.size();
  const auto N0 = CSF0s.size();

  LinAlg::Vector<double> Tc(N);

#pragma omp parallel for
  for (std::size_t i = 0; i < N; ++i) {
    double sum = 0.0;
    for (std::size_t k = 0; k < N0; ++k) {
      sum += c0[k] * RME_CSF2(CSFs.at(i), twoJ, CSF0s.at(k), twoJ0, h, K_rank);
    }
    Tc[i] = sum;
  }

  return Tc;
}

//==============================================================================
// Solves [<I|H|J> - (E_0 + w) d_IJ] c_J = -<I||T||K> c_K, for the mixed state
// c_I
PsiJPi solve_mixed_state(const PsiJPi &Psi0, std::size_t i0,
                         const PsiJPi &target,
                         const LinAlg::Matrix<double> &Hci,
                         const Coulomb::meTable<double> &h, int K_rank,
                         double omega) {

  assert(Hci.rows() == Hci.cols());
  assert(Hci.rows() == target.CSFs().size());

  const auto N = target.CSFs().size();
  const auto E0 = Psi0.energy(i0);

  // Right-hand side: -T|Psi_0>, in the CSF basis of the target J/Pi
  auto rhs =
    -1.0 * TPsi_reduced(target.CSFs(), target.twoJ(), Psi0, i0, h, K_rank);

  // Left-hand side: M_IJ = <I|H|J> - (E_0 + omega) * delta_IJ
  auto M = Hci;
  for (std::size_t i = 0; i < N; ++i) {
    M(i, i) -= E0 + omega;
  }

  // If the target J/Pi is that of Psi_0, and omega is zero, then M is singular
  // (M|Psi_0> = 0). Project Psi_0 out of the right-hand side, and add
  // |Psi_0><Psi_0| to M, which leaves M unchanged in the orthogonal subspace,
  // but makes it invertible; the solution is then orthogonal to Psi_0
  if (omega == 0.0 && target.twoJ() == Psi0.twoJ() &&
      target.parity() == Psi0.parity()) {
    assert(N == Psi0.CSFs().size());
    LinAlg::Vector<double> c0(N);
    for (std::size_t i = 0; i < N; ++i) {
      c0[i] = Psi0.coef(i0, i);
    }
    rhs -= (c0 * rhs) * c0;
    M += outer_product(c0, c0);
  }

  auto dPsi = target;
  dPsi.set_solution(E0, LinAlg::solve_Axeqb(M, rhs));

  return dPsi;
}

//==============================================================================
// Removes the given CI solutions from the mixed state (which is passed by
// value, and returned)
PsiJPi project_out(PsiJPi dPsi, const PsiJPi &levels,
                   const std::vector<std::size_t> &indices) {

  assert(dPsi.CSFs().size() == levels.CSFs().size());

  const auto N = dPsi.CSFs().size();

  LinAlg::Vector<double> c(N);
  for (std::size_t i = 0; i < N; ++i) {
    c[i] = dPsi.coef(0, i);
  }

  // Gram-Schmidt: since the levels are eigenstates of the CI Hamiltonian, and
  // the mixed state is the sum over that spectrum, this removes exactly the
  // terms of those levels
  for (const auto index : indices) {
    const auto c_p = levels.coefs(index);
    double overlap = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
      overlap += c_p[i] * c[i];
    }
    for (std::size_t i = 0; i < N; ++i) {
      c[i] -= overlap * c_p[i];
    }
  }

  dPsi.set_solution(dPsi.energy(0), c);

  return dPsi;
}

//==============================================================================
PsiJPi solve_mixed_state(const PsiJPi &Psi0, std::size_t i0, int twoJ,
                         int parity,
                         const std::vector<DiracSpinor> &ci_sp_basis,
                         const Coulomb::meTable<double> &h, int K_rank,
                         const Coulomb::meTable<double> &h1,
                         const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk,
                         const Coulomb::LkTable *Sk, double omega) {

  const PsiJPi target(twoJ, parity, ci_sp_basis);
  const auto Hci = construct_Hci(target, h1, qk, Bk, Sk);

  return solve_mixed_state(Psi0, i0, target, Hci, h, K_rank, omega);
}

} // namespace CI
