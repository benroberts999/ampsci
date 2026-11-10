#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace CI {

/*!
  @brief Orbital data for the effective three-body correction, Sigma_3: the
  core orbitals of the hole line, and the orbital energies for the energy
  denominators.

  @details
  The energies are stored by DiracSpinor::nk_index, for the core orbitals and
  the CI orbitals. Formed by @ref sigma3_orbitals.
*/
struct Sigma3Orbitals {
  //! Core orbitals summed over (the hole line)
  std::vector<DiracSpinor::Index> core{};
  //! Orbital energies, indexed by nk_index
  std::vector<double> en{};
  //! Fermi-level energy of each CI orbital, indexed by nk_index: the energy of
  //! the lowest CI orbital with the same kappa
  std::vector<double> en_bar{};
  //! Energy of the lowest CI orbital
  double en_min{0.0};
};

/*!
  @brief Forms the orbital data for Sigma_3; see @ref Sigma3Orbitals.

  @param core      Core orbitals to sum over in the hole line (e.g., those
                   with n >= n_min_core).
  @param ci_basis  Single-particle basis of the CI expansion.
*/
Sigma3Orbitals sigma3_orbitals(const std::vector<DiracSpinor> &core,
                               const std::vector<DiracSpinor> &ci_basis);

/*!
  @brief Energy factor (inverse denominator) of the Sigma_3 diagram.

  @details
  For the diagram with ket orbitals (v, w, u) and bra orbitals (x, y, z), in
  which the core electron a and the valence electron u are excited to (y, z),
  and one of (v, w) then fills the hole, the other going to x. The
  intermediate state holds the valence orbitals v, w, y, z and the hole a:
  \f[
    E_{\rm int} = \varepsilon_v + \varepsilon_w + \varepsilon_y
                + \varepsilon_z - \varepsilon_a ,
  \f]
  always with the actual orbital energies. The denominator is
  \f$ E_{\rm target} - E_{\rm int} \f$, with the energy of the target state
  depending on the mode (\f$ \bar\varepsilon \f$ is the Fermi-level energy,
  see @ref Sigma3Orbitals):

  - RS     : \f$ \varepsilon_v + \varepsilon_w + \varepsilon_u \f$
             (the ket energy), so the denominator is
             \f$ \varepsilon_a + \varepsilon_u - \varepsilon_y -
             \varepsilon_z \f$.
  - DFK    : \f$ \bar\varepsilon_v + \bar\varepsilon_w + \bar\varepsilon_u \f$.
  - BW     : @p E0, the total valence energy of the target CI level.
  - Fermi  : every valence leg (target and intermediate) at its Fermi-level
             energy: \f$ \varepsilon_a + \bar\varepsilon_u -
             \bar\varepsilon_y - \bar\varepsilon_z \f$.
  - Fermi0 : every valence leg at the energy of the lowest CI orbital,
             \f$ \varepsilon_0 \f$: \f$ \varepsilon_a - \varepsilon_0 \f$.

  As for Sigma_2 (see @ref MBPT::Denominators), the result is averaged with
  that of the bra-ket partner, in which the bra is the target state,
  (v, w, u) <-> (y, z, x), so that the CI matrix is symmetric. For BW and
  Fermi0 the two coincide.

  @param denominators  Denominator mode.
  @param orbitals      Orbital energies; see @ref Sigma3Orbitals.
  @param a             Core orbital (hole line).
  @param v,w,u         Ket orbitals: the pair (v, w), and u.
  @param x,y,z         Bra orbitals: x, and the pair (y, z).
  @param E0            Total valence energy of the target CI level; used only
                       by BW.
  @return The averaged inverse denominator.
*/
double Sigma3_inverse_denominator(MBPT::Denominators denominators,
                                  const Sigma3Orbitals &orbitals,
                                  DiracSpinor::Index a, DiracSpinor::Index v,
                                  DiracSpinor::Index w, DiracSpinor::Index u,
                                  DiracSpinor::Index x, DiracSpinor::Index y,
                                  DiracSpinor::Index z, double E0 = 0.0);

/*!
  @brief Sigma_3 between two parent states of three-electron CSFs, for a
  single assignment of the electrons.

  @details
  The effective three-body interaction of second-order MBPT: two Coulomb
  interactions joined by a core (hole) line. As an operator on the valence
  electrons,
  \f[
    \Sigma_3 = -\sum_a \sum_{vwu,\,xyz}
    \frac{g_{xavw}\, g_{yzau}}{D}\,
    a^\dagger_x a^\dagger_y a^\dagger_z a_w a_v a_u ,
  \f]
  where a runs over the core, the other sums over the valence (CI) orbitals,
  and D is the energy denominator (see @ref Sigma3_inverse_denominator).

  For fully antisymmetric three-electron states this is \f$ 3! \f$ times the
  matrix element of a single assignment of the electrons: electrons 1 and 2
  interact, electron 2 going to the core state a, and then electrons 2 and 3
  interact,
  \f[
    G = -\sum_a g_{23}\,\frac{|a\rangle_2\langle a|_2}{D}\,g_{12} .
  \f]
  This function returns the matrix element of G between the parent states
  (see @ref CSF3::Parent) \f$ \langle x, [yz]K; J| \f$ (pair on electrons 2
  and 3) and \f$ |[vw]J_2, u; J\rangle \f$ (pair on electrons 1 and 2), each
  antisymmetric in its pair only. The first interaction takes the ket pair to
  (x, a), conserving \f$ J_2 \f$; recoupling to x with (a, u), the second
  takes (a, u) to the bra pair, conserving K:
  \f[
    \langle [yz]K, x; J|G|[vw]J_2, u; J\rangle =
    \frac{1}{2}\sum_a \frac{(-1)^{K + j_a + j_u}}{D}
    \sqrt{[J_2][K]}
    \begin{Bmatrix} j_x & j_a & J_2 \\ j_u & J & K \end{Bmatrix}
    \langle [yz]K | g | [au]K \rangle
    \langle [xa]J_2 | g | [vw]J_2 \rangle ,
  \f]
  where the two-electron matrix elements are those of CSF2_Coulomb() (the
  factor of 1/2 as the pairs (x, a) and (a, u) are not antisymmetrised). See
  @ref Sigma3_AB.

  @param bra,ket       Parent states (their coef members are not used).
  @param twoJ          Twice the total angular momentum 2J.
  @param qk            Coulomb \f$ Q^k \f$ table; must hold the integrals with
                       one core orbital and three CI orbitals.
  @param orbitals      Core orbitals and energies; see @ref Sigma3Orbitals.
  @param denominators  Denominator mode.
  @param E0            Total valence energy of the target CI level (BW only).
  @return The matrix element.
*/
double Sigma3_parents(const CSF3::Parent &bra, const CSF3::Parent &ket,
                      int twoJ, const Coulomb::QkTable &qk,
                      const Sigma3Orbitals &orbitals,
                      MBPT::Denominators denominators, double E0 = 0.0);

/*!
  @brief Sigma_3 (effective three-body) correction to Hab(), for
  three-electron CSFs.

  @details
  Each CSF is expanded over its parent states (see @ref CSF3). As both CSFs
  are fully antisymmetric,
  \f[
    \langle A|\Sigma_3|B\rangle = 3! \sum_{pq} C^A_p C^B_q \,
    \langle p | G | q \rangle ,
  \f]
  with the parent matrix element of @ref Sigma3_parents. Unlike the one- and
  two-body parts of the Hamiltonian, this does not vanish when all three
  orbitals of the CSFs differ.

  @param A,B           The two CSFs.
  @param twoJ          Twice the total angular momentum 2J.
  @param qk            Coulomb \f$ Q^k \f$ table; must hold the integrals with
                       one core orbital and three CI orbitals.
  @param orbitals      Core orbitals and energies; see @ref Sigma3Orbitals.
  @param denominators  Denominator mode.
  @param E0            Total valence energy of the target CI level (BW only).
  @return \f$ \Sigma_3 \f$ correction to \f$ H_{AB} \f$.
*/
double Sigma3_AB(const CSF3 &A, const CSF3 &B, int twoJ,
                 const Coulomb::QkTable &qk, const Sigma3Orbitals &orbitals,
                 MBPT::Denominators denominators, double E0 = 0.0);

} // namespace CI
