#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cstddef>
#include <vector>

namespace CI {

/*!
  @brief Action of a one-body operator on a CI state, in the CSF basis.
  @details
  Forms the vector

  \f[
    T_I = \sum_K \redmatel{I}{T^{(K)}}{K} \, c^{(0)}_K,
  \f]

  where \f$ I \f$ runs over the CSFs in @p CSFs (which have total angular
  momentum @p twoJ /2), \f$ K \f$ runs over the CSFs of the reference state
  \f$ \Psi_0 \f$ (solution @p i0 of @p Psi0), and \f$ c^{(0)}_K \f$ are its CI
  expansion coefficients.

  The CSF matrix elements are reduced (see @ref RME_CSF), formed from the
  single-particle reduced matrix elements in @p h; so is the result.

  @param CSFs   CSFs spanning the J/Pi the operator maps into.
  @param twoJ   Twice the total angular momentum, 2J, of @p CSFs.
  @param Psi0   CI solutions containing the reference state.
  @param i0     Index of the reference solution within @p Psi0.
  @param h      Table of single-particle reduced matrix elements of T.
  @param K_rank Rank of the tensor operator T.
  @return Vector \f$ T_I \f$, of length @p CSFs .size().
*/
[[nodiscard]] LinAlg::Vector<double>
TPsi_reduced(const std::vector<CSF> &CSFs, int twoJ, const PsiJPi &Psi0,
             std::size_t i0, const Coulomb::meTable<double> &h, int K_rank);

/*!
  @brief Solves the CI mixed-states (Sternheimer) equation for a one-body
  operator.
  @details
  Finds the first-order correction to the CI state \f$ \Psi_0 \f$ (solution
  @p i0 of @p Psi0, with energy \f$ E_0 \f$) due to the one-body operator
  \f$ T^{(K)} \f$, expanded over the CSFs of a single J/Pi:

  \f[
    \ket{\delta\Psi} = \sum_I c_I \ket{I; J^\pi}.
  \f]

  The coefficients solve the linear system

  \f[
    \sum_J \left[ \matel{I}{H}{J} - (E_0 + \omega) \, \delta_{IJ} \right] c_J
      = - \sum_K \redmatel{I}{T^{(K)}}{K} \, c^{(0)}_K,
  \f]

  where \f$ H \f$ is the CI Hamiltonian in the J/Pi defined by @p target
  (given as the matrix @p Hci, e.g., from @ref construct_Hci), and the
  right-hand side is formed by @ref TPsi_reduced.

  Since the right-hand side is reduced (in \f$ T \f$), so is the solution: for
  any CI state \f$ A \f$ in the same J/Pi, the mixed state satisfies

  \f[
    \sum_I c^A_I \, c_I
      = \frac{\redmatel{A}{T^{(K)}}{\Psi_0}}{E_0 + \omega - E_A},
  \f]

  i.e., the sum over the entire spectrum of that J/Pi, without finding (or
  summing over) the individual CI solutions.

  Returned as a @ref PsiJPi holding a single "solution", the coefficients
  \f$ c_I \f$. Its stored energy is \f$ E_0 \f$, that of the reference state
  (not an eigenvalue).

  @param Psi0    CI solutions containing the reference state.
  @param i0      Index of the reference solution within @p Psi0.
  @param target  Defines the J/Pi the mixed state lives in (2J, parity, and
                 CSF list); its solutions, if any, are not used.
  @param Hci     CI Hamiltonian matrix in the CSF basis of @p target.
  @param h       Table of single-particle reduced matrix elements of T.
  @param K_rank  Rank of the tensor operator T.
  @param omega   Frequency: the mixed state due to a time-dependent operator,
                 \f$ T e^{-i\omega t} \f$, has denominators
                 \f$ E_0 + \omega - E_A \f$ [0].
  @return PsiJPi for the @p target J/Pi, holding the single mixed state.
  @see project_out, to remove individual levels from the mixed state.

  @note If @p target has the same J and parity as @p Psi0, and
        \f$ \omega = 0 \f$, the matrix on the left is singular, since
        \f$ \Psi_0 \f$ itself has zero eigenvalue. In that case,
        \f$ \Psi_0 \f$ is projected out (equivalent to subtracting
        \f$ \redmatel{\Psi_0}{T}{\Psi_0} \f$ from the right-hand side).

  @note Any other state degenerate with \f$ E_0 + \omega \f$ also makes the
        system singular. Its term in the sum over states is divergent, and must
        be dealt with separately, as in degenerate perturbation theory.

  @note If the operator cannot connect the two J/Pi (triangle rule or
        parity), the right-hand side vanishes, and the mixed state is zero.
*/
[[nodiscard]] PsiJPi solve_mixed_state(const PsiJPi &Psi0, std::size_t i0,
                                       const PsiJPi &target,
                                       const LinAlg::Matrix<double> &Hci,
                                       const Coulomb::meTable<double> &h,
                                       int K_rank, double omega = 0.0);

/*!
  @brief Removes CI levels from a mixed state, so that it is orthogonal to them.
  @details
  A mixed state is implicitly a sum over the entire spectrum of its (J, parity):

  \f[
    \ket{\delta\Psi} = \sum_A \ket{A}
      \frac{\redmatel{A}{T^{(K)}}{\Psi_0}}{E_0 + \omega - E_A}.
  \f]

  Subtracting the projection onto the listed levels removes exactly their terms
  from that sum, so they may be treated separately: e.g., with experimental
  energies or matrix elements.

  @param dPsi    Mixed state, from @ref solve_mixed_state (taken by value).
  @param levels  Solved CI levels of the same (J, parity): the eigenstates of
                 the CI Hamiltonian used for the mixed state.
  @param indices Which solutions of @p levels to remove.
  @return The mixed state, orthogonal to the listed levels.

  @note Removing a level degenerate with \f$ E_0 + \omega \f$ does not help:
        the mixed state itself does not exist in that case (see
        @ref solve_mixed_state).
*/
[[nodiscard]] PsiJPi project_out(PsiJPi dPsi, const PsiJPi &levels,
                                 const std::vector<std::size_t> &indices);

/*!
  @brief Solves the CI mixed-states equation; constructs the CI matrix
  internally.
  @details
  Convenience overload of @ref solve_mixed_state: forms the CSFs for the
  requested J/Pi (@p twoJ, @p parity) from @p ci_sp_basis, constructs the CI
  Hamiltonian matrix via @ref construct_Hci, then solves the mixed-states
  equation.

  Use the other overload if the CI matrix for the target J/Pi is already
  available (e.g., when several operators are considered).

  @param Psi0        CI solutions containing the reference state.
  @param i0          Index of the reference solution within @p Psi0.
  @param twoJ        Twice the total angular momentum, 2J, of the mixed state.
  @param parity      Parity of the mixed state: +1 or -1. This is the parity of
                     @p Psi0 times that of the operator.
  @param ci_sp_basis Single-particle basis used to construct the CSFs.
  @param h           Table of single-particle reduced matrix elements of T.
  @param K_rank      Rank of the tensor operator T.
  @param h1          One-body matrix element table (may include Sigma_1).
  @param qk          Coulomb Q^k integral table.
  @param Bk          Pointer to Breit W^k table; ignored if nullptr.
  @param Sk          Pointer to Sigma_2 L^k table; ignored if nullptr.
  @param omega       Frequency; see the other overload [0].
  @return PsiJPi for the J/Pi (@p twoJ, @p parity), holding the single mixed
          state.
*/
[[nodiscard]] PsiJPi
solve_mixed_state(const PsiJPi &Psi0, std::size_t i0, int twoJ, int parity,
                  const std::vector<DiracSpinor> &ci_sp_basis,
                  const Coulomb::meTable<double> &h, int K_rank,
                  const Coulomb::meTable<double> &h1,
                  const Coulomb::QkTable &qk,
                  const Coulomb::WkTable *Bk = nullptr,
                  const Coulomb::WkTable *Sk = nullptr, double omega = 0.0);

} // namespace CI
