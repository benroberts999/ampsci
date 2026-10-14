#pragma once
#include "CI/CSF.hpp"
#include "Coulomb/meTable.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
namespace DiracOperator {
class TensorOperator;
}

namespace Amplitudes {

/*!
  @brief Second-order amplitude A^K between two CI states, by
  sum-over-states over the solved CI levels.
  @details
  Evaluates \f$ A^K \f$ (see CI::A_K_coefs) directly,

  \f[
    A^K = \sum_n \left[
      c_1(J_n)\,
      \frac{\redmatel{b}{t}{n}\redmatel{n}{s}{a}}{E_a + \omega_s - E_n}
    + c_2(J_n)\,
      \frac{\redmatel{b}{s}{n}\redmatel{n}{t}{a}}{E_a + \omega - E_n}
    \right],
  \f]

  the sum running over the solutions of each allowed intermediate
  (J, parity). This is the sum-over-states analogue of CI::A_K (the
  generalisation of the CI_Pol module): only the levels actually solved in
  the CI{} block are available, so the sum is truncated - both to the
  J/Pi that were solved (a note is printed for any that are
  missing) and to num_solutions levels within each. CI::A_K, by contrast,
  is complete; the comparison shows how much of the amplitude the low
  levels carry.

  The reduced matrix elements are the CI contractions of the
  single-particle tables (CI::ReducedME). Corrections to the matrix
  elements (RPA, structure radiation, normalisation of states) enter
  through those tables.

  @param K       Rank of the amplitude.
  @param Psi_b,ib  Final CI state (solution @p ib of @p Psi_b).
  @param Psi_a,ia  Initial CI state.
  @param t,t_me  The \f$ t \f$ operator, and its table of single-particle
                 reduced matrix elements (which may include RPA, structure
                 radiation, normalisation), formed at @p omega.
  @param s,s_me  The \f$ s \f$ operator, and its table (formed at
                 @p omega_s).
  @param omega,omega_s  Frequency of each operator; see CI::A_K.
  @param ciwfs   The solved CI J/Pi, e.g., Wavefunction::CIwfs(): the
                 intermediate states of the sum.
  @param levels_to_remove  CI levels skipped in the sum, so that they may be
                 treated separately - e.g., with experimental energies.
  @param outstream  Stream for the per-(J, parity) contributions.
  @return \f$ A^K \f$ (valence part).

  @note The polarisation of the closed core and the core-valence terms lie
        outside the CI space: see CI::A_K_core and CI::A_K_cv.

  @note A level degenerate with a denominator (\f$ E_a + \omega_s \f$ or
        \f$ E_a + \omega \f$) makes its term diverge, and nothing is skipped
        here: remove the offending level (@p levels_to_remove) and treat it
        separately.
*/
[[nodiscard]] double sos_ci(int K, const CI::PsiJPi &Psi_b, std::size_t ib,
                            const CI::PsiJPi &Psi_a, std::size_t ia,
                            const DiracOperator::TensorOperator *t,
                            const Coulomb::meTable<double> &t_me,
                            const DiracOperator::TensorOperator *s,
                            const Coulomb::meTable<double> &s_me, double omega,
                            double omega_s,
                            const std::vector<CI::PsiJPi> &ciwfs,
                            const std::vector<CI::Level> &levels_to_remove = {},
                            std::ostream &outstream = std::cout);

} // namespace Amplitudes
