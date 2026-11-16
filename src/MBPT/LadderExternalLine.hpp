#pragma once
#include "Coulomb/include.hpp"
#include "MBPT/LadderIteration.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cstddef>
#include <vector>

namespace MBPT {

/*!
  @brief Ladder integrals of Sigma_L(e_v) for one valence state, one set per
         projection state x.

  @details
  Keys are the real DiracSpinor::nk_index() of every orbital, x included.
  L_mnxb is an LkTable: its symmetry L_mnxb = L_nmbx holds, since b is core
  and x is not. L_xnab is an NkTable: the symmetry L_xnab = L_nxba would map
  the integral with line x onto the one with line n, when n is also a
  projection state.
*/
struct ExternalLineLadder {
  //! Energy of the external line, e_v
  double en_v;
  //! The projection states x
  std::vector<DiracSpinor> x_orbitals;
  //! L^k_mnxb(e_v): m,n excited, x as the i orbital, b core
  Coulomb::LkTable L_mnxb{};
  //! L^k_xnab(e_v): x as the m orbital, n excited, a,b core
  Coulomb::NkTable L_xnab{};
};

/*!
  @brief Starting point of the ladder iteration with the external line: each
         L_x is L_v (from the ladder table) for x = v, and zero otherwise.

  @details
  All entries are added. The projection states x are the orbitals of
  @p excited with the kappa of @p v (v must be one of them).
*/
ExternalLineLadder external_line_start(const DiracSpinor &v,
                                       const std::vector<DiracSpinor> &core,
                                       const std::vector<DiracSpinor> &excited,
                                       const Coulomb::LkTable &lk);

/*!
  @brief One sweep of the ladder equation for L^k_mnxb(e_v), all x at once.

  @details
  calculate_Lk_mnib() with i_orbs = x_orbitals: the dressed integrals of L1
  and L3 that contain x come from @p L_mnxb_prev, at e_v; all others from
  @p lk (L2 and L4 hold x only in a bare Q). Only entries already in
  @p L_mnxb are written, damped by @p a_damp [see damped()].
*/
void calculate_Lk_mnxb(Coulomb::LkTable *L_mnxb, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &x_orbitals, double en_v,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable &lk,
                       const Coulomb::LkTable &L_mnxb_prev, double a_damp);

/*!
  @brief One sweep of the ladder equation for L^k_xnab(e_v), all x at once.

  @details
  The diagrams are those of calculate_Lk_mnib() with m = x and (i,j) = (a,b),
  for all n. The dressed integrals of L2 (row (x,b) of L2^K) and of L4
  [QL^J_{xn,cd}] contain x: they come from @p L_xnab_prev, at e_v. L1 and L3
  hold x only in a bare Q, and take theirs from @p lk. L3 is read from
  L2^K_{na,xb}, in which n is an internal line even when it is also a
  projection state, so it comes from a second L2^K matrix, with no external
  line. Only entries already in @p L_xnab are written, damped by @p a_damp
  [see damped()].
*/
void calculate_Lk_xnab(Coulomb::NkTable *L_xnab, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &x_orbitals, double en_v,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable &lk,
                       const Coulomb::NkTable &L_xnab_prev, double a_damp);

/*!
  @brief Writes L^k_xnab = (L1+L4) + L2 + L3 for one x and the n of one kappa
         group into the table.

  @details
  As store_Lk_mnij() with m = x and (i,j) = (a,b), x = excited[index_x],
  except that L2 is read from @p L2K_mb_na_line [element (x,b),(n,a); the line
  as m] and L3 from @p L2K_mb_na [element (n,a),(x,b); no external line].
  Each entry is damped by @p a_damp [see damped()].
*/
void store_Lk_xnab(Coulomb::NkTable *L_xnab, std::size_t index_x,
                   const std::vector<std::size_t> &n_group,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &core,
                   const OrbitalPairIndex &ab_pairs,
                   const OrbitalPairIndex &ph_pairs,
                   const std::vector<LinAlg::Matrix<double>> &L1_L4_J_xn_ab,
                   const std::vector<LinAlg::Matrix<double>> &L2K_mb_na_line,
                   const std::vector<LinAlg::Matrix<double>> &L2K_mb_na,
                   const SixJArray &sixj, double a_damp);

//! Largest change of any integral in either table of ladder from ladder_prev,
//! relative to the largest integral (same entries in both). ladder is
//! non-const only because CoulombTable gives no const access to its entries.
double external_line_change(ExternalLineLadder &ladder,
                            const ExternalLineLadder &ladder_prev);

/*!
  @brief The ladder integrals of Sigma_L(e_v) for valence state v, iterated
         with the ladder table fixed.

  @details
  Starts from external_line_start(), and sweeps both tables
  [calculate_Lk_mnxb(), calculate_Lk_xnab()] until the largest change of any
  integral, relative to the largest integral, is below @p eps_target, or for
  at most @p max_it sweeps; prints a warning if not converged. With
  @p max_it = 0 only the x = v integrals are non-zero, so Sigma_L is only the
  column Sigma_L|v>.

  @param v           Valence state (the basis version, one of @p excited)
  @param core        Core orbitals
  @param excited     Excited orbitals
  @param qk          Coulomb integral table
  @param lk          Converged ladder table (i orbitals: core and v)
  @param diagrams    Which diagrams to include
  @param max_it      Largest number of sweeps
  @param eps_target  Convergence target
  @param a_damp      Damping factor [see damped()]
  @param print       Print the convergence
*/
ExternalLineLadder
external_line_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                     const LadderDiagrams &diagrams, int max_it,
                     double eps_target, double a_damp, bool print);

/*!
  @brief Ladder correlation potential from ladder integrals with the
         projection states x in the third (L^k_mnxa) and first (L^k_xnab)
         index.

  @details
  \f[
    \Sigma = \sum_{x}\Big[\sum_{amn,k} |W^k_{\cdot amn}\rangle\,
      \frac{L^k_{mnxa}}
           {[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}
      + \sum_{nab,k} |W^k_{\cdot nab}\rangle\,
      \frac{L^k_{xnab}}
           {[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}
      \Big]\langle x| .
  \f]
  MnxbTable is Coulomb::LkTable [ExternalLineLadder::L_mnxb]; XnabTable is
  Coulomb::NkTable [ExternalLineLadder::L_xnab]. For x = v alone, both may be
  the ladder table itself (Coulomb::LkTable), which holds L^k_mnva and
  L^k_vnab.

  @param v               Valence state: kappa, j and energy
  @param x_orbitals      Projection states x
  @param L_mnxb          L^k_mnxa: m,n excited, x, a core
  @param L_xnab          L^k_xnab: x, n excited, a,b core
  @param core            Core orbitals
  @param excited         Excited orbitals
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G       Include the lower (g) component
*/
template <typename MnxbTable, typename XnabTable>
GMatrix Sigma_ladder(const DiracSpinor &v,
                     const std::vector<DiracSpinor> &x_orbitals,
                     const MnxbTable &L_mnxb, const XnabTable &L_xnab,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, double r0,
                     double rmax, std::size_t stride, bool include_G);

/*!
  @brief Ladder correlation potential Sigma_L(e_v) from the ladder integrals
         with the external line.

  @details
  Sigma_L(e_v) by projection onto the excited states x of kappa_v, with the
  ladder iterated for each x [external_line_ladder()]: the Sigma_ladder()
  template with the integrals L^k_mnxa(e_v) and L^k_xnab(e_v) of @p ladder.
  x is the external line: it carries the valence energy e_v in the energy
  denominators, not its own [see ExternalLine]. L^k_mnxb (x as the i orbital)
  and L^k_xnab (x as the m orbital) each obey the ladder equation by
  themselves: the dressed integrals that contain x come from the same table,
  all others from the converged ladder table [calculate_Lk_mnxb(),
  calculate_Lk_xnab()].

  For x = v these are the integrals of the ladder table, so the column
  Sigma_L|v> is exact, and <v|Sigma_L|v> is de_valence_w().

  @param v               Valence state: kappa, j and energy
  @param ladder          Ladder integrals of v [external_line_ladder()]
  @param core            Core orbitals
  @param excited         Excited orbitals
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G       Include the lower (g) component
*/
GMatrix Sigma_ladder(const DiracSpinor &v, const ExternalLineLadder &ladder,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, double r0,
                     double rmax, std::size_t stride, bool include_G);

} // namespace MBPT
