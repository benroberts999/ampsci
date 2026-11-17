#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace MBPT {

// Which diagrams to include [see LadderIteration.hpp]
struct LadderDiagrams;

/*!
  @brief Fills the ladder integral table for all _new_ index combinations.

  @details
  Iterates over all combinations of excited pairs \f$ (m,n) \f$ and
  orbitals in @p i_orbs, computing \f$ L^k_{mnib} \f$ (with \f$ b \f$ in
  @p core) and storing results in @p lk.
  Only calculates new integrals. Only lowest-order.
  All the integrals are calculated at once, see calculate_Lk_mnib().

  @param lk          Output ladder table (written in place)
  @param qk          Coulomb \f$ Q^k \f$ integral table
  @param excited     Excited orbitals
  @param core        Core orbitals
  @param i_orbs      Orbitals for the \f$ i \f$ index
  @param diagrams    Which diagrams to include [see calculate_Lk_mnib()]
  @param max_k       Maximum multipolarity; -1 uses qk.max_k()
  @param print       Print progress and timing to screen
*/
void fill_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs,
                  const LadderDiagrams &diagrams, int max_k = -1,
                  bool print = true);

/*!
  @brief Updates the ladder integral table with L(Q,Q) -> L(Q,Q+L)

  @details
  Re-calculates the existing entries \f$ L^k_{mnib} \f$ of @p lk (those with
  \f$ i \f$ in @p update_i), using @p lk_prev for the dressed integrals (Q+L),
  and stores results in @p lk. Designed for iterative refinement: pass the
  previous iteration's table as @p lk_prev.
  All the integrals are calculated at once, see calculate_Lk_mnib().

  @param lk          Output ladder table (written in place)
  @param qk          Coulomb \f$ Q^k \f$ integral table
  @param excited     Excited orbitals
  @param core        Core orbitals
  @param update_i    Restrict re-iteration to entries whose i index is in this
                     set (b is always core). Used to converge core
                     (update_i=core) before valence (update_i=valence).
  @param diagrams    Which diagrams to include [see calculate_Lk_mnib()]
  @param lk_prev     Ladder table from previous iteration
  @param a_damp      Damping factor [0,1) : 0 means no damping.
                     L -> a_damp*L + (1-a_damp)*L_new, where L is the
                     existing entry of @p lk (not @p lk_prev)
  @param print       Print progress and timing to screen
  @param rho_ma      Optional: core single excitations rho_ma (rows excited,
                     columns core) for the core relaxation R [see
                     calculate_Lk_mnib()]; nullptr: no R

  @note Does not calculate any new integrals - assumes all already present.
  Just updates them (based on iterative rule: L(Q,Q) -> L(Q,Q+L))
*/
void update_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                    const std::vector<DiracSpinor> &excited,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &update_i,
                    const LadderDiagrams &diagrams,
                    const Coulomb::LkTable *const lk_prev, double a_damp,
                    bool print, const LinAlg::Matrix<double> *rho_ma = nullptr);

/*!
  @brief Second-order (or ladder) correction to the valence energy.

  @details
  Computes the correlation energy shift
  \f[
    \delta\epsilon_v = \frac{1}{[j_v]} \sum_{amn,k}
      \frac{Q^k_{vamn}\,(L^k_{mnva} + P^k_{nmav})}
           {[k]\,(\epsilon_v + \epsilon_a - \epsilon_m - \epsilon_n)}
    + \frac{1}{[j_v]} \sum_{nab,k}
      \frac{Q^k_{vnab}\,(L^k_{vnab} + P^k_{vnab})}
           {[k]\,(\epsilon_v + \epsilon_n - \epsilon_a - \epsilon_b)},
  \f]
  where \f$ L^k \f$ and its exchange \f$ P^k \f$ (lk.P) both come from
  @p lk, \f$ a,b \f$ run over core, \f$ m,n \f$ over excited, and \f$ k \f$
  over the Coulomb-parity range only.
  When @p lk holds plain Coulomb integrals the result is the
  MBPT(2) correction; when @p lk holds ladder integrals it is the full
  ladder correction.

  @param v       Valence orbital
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param lk      \f$ Q^k \f$ or ladder \f$ L^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @return \f$ \delta\epsilon_v \f$
*/
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited);

/*!
  @brief Ladder (or MBPT2) valence energy, antisymmetrising the FIRST integral.

  @details
  Computes
  \f[
    \delta\epsilon_v = \sum_{amn,k}
      \frac{W^k_{vamn}\,L^k_{mnva}}{[k]\,[j_v]\,\Delta\epsilon}
      + \text{(c+d)},
  \f]
  with \f$ W = Q + P \f$ the antisymmetrised Coulomb integral (@p qk.W) and the
  ladder @p lk entering only through the direct integral @p lk.Q. Equivalent to
  de_valence() (which instead antisymmetrises the ladder via @p lk.P), since the
  exchange symmetry holds under the full sum. Unlike de_valence(), this works
  only with the Coulomb integrals in the first slot (the ladder lacks the
  required symmetry). No screening/eta is applied.

  @param v       Valence orbital
  @param qk      Coulomb integrals supplying \f$ W = Q+P \f$ (first integral)
  @param lk      Ladder \f$ L^k \f$ integrals (second, direct integral)
  @param core    Core orbitals
  @param excited Excited orbitals
  @param sj      Optional 6j table (speeds up the W exchange sum)
  @return \f$ \delta\epsilon_v \f$

  @warning @p sj must hold every 6j symbol of @p v, @p core and @p excited:
           a symbol missing from it is taken as 0, giving a wrong energy.
*/
template <typename Qintegrals, typename Lintegrals>
double de_valence_w(const DiracSpinor &v, const Qintegrals &qk,
                    const Lintegrals &lk, const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &excited,
                    const Angular::SixJTable *sj = nullptr);

/*!
  @brief de_valence_w() with the integrals L^k_vnab of the second term,
         diagrams (c+d), from a closing table of their own.

  @details
  With the core relaxation, the amplitude closed in the second term is not
  the entry of the ladder table: it has no single excitations of the core
  into the valence orbital [calculate_Lk_xnab()]. Their contribution is a
  chain through a core orbital, which the Brueckner equation generates.

  @param v       Valence orbital
  @param qk      Coulomb integrals supplying \f$ W = Q+P \f$
  @param lk      Ladder table: L^k_mnva
  @param L_vnab  Closing table: L^k_vnab, keyed (v,n,a,b)
  @param core    Core orbitals
  @param excited Excited orbitals
  @param sj      Optional 6j table [see de_valence_w()]
*/
template <typename Qintegrals, typename Lintegrals, typename VnabIntegrals>
double de_valence_w(const DiracSpinor &v, const Qintegrals &qk,
                    const Lintegrals &lk, const VnabIntegrals &L_vnab,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &excited,
                    const Angular::SixJTable *sj = nullptr);

/*!
  @brief Direct (Q) closing of x_F: the part of the valence energy that the
         Feynman direct diagram contains, in the basis.

  @details
  \f[
    \delta\epsilon^F_v = \frac{1}{[j_v]} \sum_{amn,k}
      \frac{Q^k_{vamn}\,L^{F,k}_{mnva}}
           {[k]\,(\epsilon_v + \epsilon_a - \epsilon_m - \epsilon_n)}
    + \frac{1}{[j_v]} \sum_{nab,k}
      \frac{Q^k_{vnab}\,L'^k_{vnab}}
           {[k]\,(\epsilon_v + \epsilon_n - \epsilon_a - \epsilon_b)},
  \f]
  diagrams (a) and (c) of de_valence(), k over the Coulomb range only. L^F is
  x_F [calculate_Lk_mnib() with P1 and hole_particle], and L' the closing
  amplitude:
  x_F_vnab without the hole-particle term on (v,a), which is part of the
  valence line there [calculate_Lk_vnab_F()]. The correction to the valence
  energy is de_valence_w() of x minus this.

  @param v        Valence orbital
  @param qk       Coulomb \f$ Q^k \f$ integral table
  @param lk_F     x_F table: L^k_mnva
  @param L_vnab   Closing table: L'^k_vnab
  @param core     Core orbitals
  @param excited  Excited orbitals
  @return \f$ \delta\epsilon^F_v \f$
*/
template <typename Qintegrals, typename Lintegrals, typename VnabIntegrals>
double de_valence_F(const DiracSpinor &v, const Qintegrals &qk,
                    const Lintegrals &lk_F, const VnabIntegrals &L_vnab,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &excited);

/*!
  @brief Second-order (or ladder) correction to the core energy.

  @details
  Sums the correlation energy shift over all core orbitals. When @p lk holds
  plain Coulomb integrals the result is the MBPT(2) core correction; when
  @p lk holds ladder integrals it is the ladder correction.

  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param lk      \f$ Q^k \f$ or ladder \f$ L^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @return Total core correlation energy shift
*/
template <typename Qintegrals, typename QorLintegrals>
double de_core(const Qintegrals &qk, const QorLintegrals &lk,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited);

/*!
  @brief Ladder (or MBPT2) core energy, antisymmetrising the FIRST integral.

  @details
  Computes
  \f[
    \delta E_{\rm core} = \frac{1}{2} \sum_{abmn,k}
      \frac{W^k_{abmn}\,L^k_{mnab}}
           {[k]\,(\epsilon_a + \epsilon_b - \epsilon_m - \epsilon_n)},
  \f]
  with \f$ W = Q + P \f$ the antisymmetrised Coulomb integral (@p qk.W) and the
  ladder @p lk entering only through the direct integral @p lk.Q; \f$ k \f$
  runs over both parities. Equivalent to de_core() (which instead
  antisymmetrises the ladder via @p lk.P), since the exchange symmetry holds
  under the full sum. Works only with the Coulomb integrals in the first slot
  (the ladder lacks the required symmetry).

  @param qk      Coulomb integrals supplying \f$ W = Q+P \f$ (first integral)
  @param lk      Ladder \f$ L^k \f$ integrals (second, direct integral)
  @param core    Core orbitals
  @param excited Excited orbitals
  @param sj      Optional 6j table (speeds up the W exchange sum)
  @return Total core correlation energy shift

  @warning @p sj must hold every 6j symbol of @p core and @p excited: a
           symbol missing from it is taken as 0, giving a wrong energy.
*/
template <typename Qintegrals, typename Lintegrals>
double de_core_w(const Qintegrals &qk, const Lintegrals &lk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const Angular::SixJTable *sj = nullptr);

//==============================================================================

//! Occupation blocks of the change of the core density at third order
//! [relaxation_density()]: m-scheme values, zero between orbitals of
//! different kappa
struct RelaxationDensity {
  //! gamma_cb, core orbitals: rows b, columns c
  LinAlg::Matrix<double> gamma_cb;
  //! gamma_ms, excited orbitals: rows m, columns s
  LinAlg::Matrix<double> gamma_ms;
};

/*!
  @brief Occupation blocks gamma_cb and gamma_ms of the change of the core
         density, at third order.

  @details
  Besides the core single excitations, gamma_am = gamma_ma = rho_ma
  [core_singles()], the second-order change of the core density matrix has
  the blocks
  \f[
    \gamma_{cb} = -\frac{1}{[j_b]} \sum_k \frac{1}{[k]} \sum_{amn}
      \frac{Q^k_{abmn}\,W^k_{mnac}}
           {\Delta\epsilon_{ab,mn}\,\Delta\epsilon_{ac,mn}},
  \f]
  \f[
    \gamma_{ms} = \frac{1}{[j_m]} \sum_k \frac{1}{[k]} \sum_{abn}
      \frac{Q^k_{abmn}\,W^k_{nsba}}
           {\Delta\epsilon_{ab,mn}\,\Delta\epsilon_{ab,ns}},
  \f]
  a, b, c core and m, n, s excited, kappa_c = kappa_b and kappa_s = kappa_m:
  the core orbitals are partly emptied, and the excited orbitals partly
  filled. In the m-scheme, with \f$ \rho_{mnab} = g_{mnab}/
  \Delta\epsilon_{ab,mn} \f$, they are
  \f$ \gamma_{cb} = -\sum_{amn} \rho_{mnab}\widetilde\rho_{mnac} \f$ and
  \f$ \gamma_{ms} = \sum_{abn} \rho_{mnab}\widetilde\rho_{snab} \f$, the
  same for every projection. Both are symmetric.

  @param qk       Coulomb integral table
  @param core     Core orbitals
  @param excited  Excited orbitals
*/
RelaxationDensity relaxation_density(const Coulomb::QkTable &qk,
                                     const std::vector<DiracSpinor> &core,
                                     const std::vector<DiracSpinor> &excited);

/*!
  @brief Core relaxation correction to the valence energy, <v|dV^R|v>: the
         interaction of the valence electron with the change of the core
         density.

  @details
  \f[
    \langle v|\delta V^{\rm R}|v\rangle = \sum_{pq}\widetilde g_{vpvq}
      \gamma_{qp}
      = 2\sum_{am} f_{va} W^0_{vavm} \rho_{ma}
      + \sum_{bc} f_{vb} W^0_{vbvc} \gamma_{cb}
      + \sum_{ms} f_{vm} W^0_{vmvs} \gamma_{ms},
  \f]
  with f the closed_shell_factor() and each pair of one kappa. The first
  term is the core-excited block, gamma_am = gamma_ma = rho_ma: the core
  single excitations [core_singles()], and their mirror image. The others
  are the occupation blocks [relaxation_density()]. With rho_ma at lowest
  order, this is the full third-order relaxation (all 24 time orderings),
  and the first term alone is twice the part the SD method contains.

  @param v        Valence orbital
  @param qk       Coulomb integral table
  @param rho_ma   Core single excitations, rows excited and columns core
  @param density  Occupation blocks [relaxation_density()]
  @param core     Core orbitals: columns of rho_ma
  @param excited  Excited orbitals: rows of rho_ma
  @return {rho_ma part, gamma_cb + gamma_ms part}; their sum is <v|dV^R|v>
*/
std::array<double, 2>
de_valence_relaxation(const DiracSpinor &v, const Coulomb::QkTable &qk,
                      const LinAlg::Matrix<double> &rho_ma,
                      const RelaxationDensity &density,
                      const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited);

/*!
  @brief Matrix elements <x|dV^R|y> of the core relaxation potential between
         orbitals of one kappa.

  @details
  dV^R is the interaction with the change of the core density, a one-body
  potential, the same for every valence state and independent of the energy:
  \f[
    \langle x|\delta V^{\rm R}|y\rangle = \sum_{pq}\widetilde g_{xpyq}
      \gamma_{qp}
      = \sum_{am} f_{xa} (W^0_{xaym} + W^0_{xmya}) \rho_{ma}
      + \sum_{bc} f_{xb} W^0_{xbyc} \gamma_{cb}
      + \sum_{ms} f_{xm} W^0_{xmys} \gamma_{ms},
  \f]
  with f the closed_shell_factor() and each pair of one kappa. The two terms
  of the first sum are the core single excitations and their mirror image
  (gamma_am = gamma_ma = rho_ma). Symmetric; the diagonal element of a
  valence state is the sum of the two parts of de_valence_relaxation().

  @param orbitals  Orbitals x, y: all of one kappa, and in @p qk with @p core
                   and @p excited
  @param qk        Coulomb integral table
  @param rho_ma    Core single excitations, rows excited and columns core
  @param density   Occupation blocks [relaxation_density()]
  @param core      Core orbitals: columns of rho_ma
  @param excited   Excited orbitals: rows of rho_ma
  @return <x|dV^R|y>, rows x and columns y by position in @p orbitals
*/
LinAlg::Matrix<double> relaxation_potential_matrix(
  const std::vector<DiracSpinor> &orbitals, const Coulomb::QkTable &qk,
  const LinAlg::Matrix<double> &rho_ma, const RelaxationDensity &density,
  const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &excited);

/*!
  @brief Core relaxation potential dV^R, projected onto a set of orbitals of
         one kappa.

  @details
  \f[
    \delta V^{\rm R} = \sum_{xy} |x\rangle\,
      \langle x|\delta V^{\rm R}|y\rangle\,\langle y|,
  \f]
  x, y in @p orbitals, with the matrix elements of
  relaxation_potential_matrix(): from the Coulomb integral table, so the
  near cancellation of the two occupation blocks is exact, not subject to
  the quadrature of the sub-grid. It acts only within the span of the
  orbitals. This is the part of the correction to the correlation potential
  that does not depend on the energy; the rest is from the ladder integrals
  with the external line [Sigma_ladder() template].

  @param orbitals        Projection states: orthonormal, all of one kappa
  @param qk              Coulomb integral table
  @param rho_ma          Core single excitations, rows excited, columns core
  @param density         Occupation blocks [relaxation_density()]
  @param core            Core orbitals
  @param excited         Excited orbitals
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G       Include the lower (g) component
*/
GMatrix relaxation_potential(const std::vector<DiracSpinor> &orbitals,
                             const Coulomb::QkTable &qk,
                             const LinAlg::Matrix<double> &rho_ma,
                             const RelaxationDensity &density,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<DiracSpinor> &excited, double r0,
                             double rmax, std::size_t stride, bool include_G);

/*!
  @brief Method used to construct the ladder correlation potential, Sigma_L.

  @details
  - basis : project onto the excited states of kappa_v, with the ladder
            iterated for each [external_line_ladder(), Sigma_ladder()]
  - ratio : no projection; rescale each Sigma(2) term by L/Q
            [Sigma_ladder()]
  - ratio2 : exact acting on the valence state, ratio otherwise
             [Sigma_ladder_ratio2()]
*/
enum class SigmaLMethod { basis, ratio, ratio2 };

/*!
  @brief Converts string (name) to SigmaLMethod enum (case-insensitive).

  @details
  Also accepts "full" (basis) and "Dzuba" (ratio). Warns and defaults to
  basis if unknown.
*/
SigmaLMethod parseSigmaLMethod(const std::string &method);
//! Converts SigmaLMethod enum to string (name)
std::string parseSigmaLMethod(SigmaLMethod method);

/*!
  @brief Ladder-diagram correction to the correlation potential, Sigma_L(e_v),
         by the L/Q ratio method.

  @details
  Uses the ratio method (following
  V. A. Dzuba, Phys. Rev. A 78, 042502 (2008)): no projection; each term of
  the regular second-order correlation potential (cf. Goldstone::Sigma_both)
  is rescaled by the scalar ratio of ladder to Coulomb integrals:
  \f[
    \Sigma_L = \sum_{amn,k} |Q^k_{\cdot amn}\rangle\,
      \frac{L^k_{mnva}/Q^k_{mnva}}
           {[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}\,
      \langle W^k_{\cdot amn}|
      + \sum_{nab,k} |Q^k_{\cdot nab}\rangle\,
      \frac{L^k_{vnab}/Q^k_{vnab}}
           {[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}\,
      \langle W^k_{\cdot nab}| .
  \f]
  By construction the diagonal reproduces the ladder energy exactly,
  <v|Sigma_L|v> = de_valence_w(v), except for the dropped terms (see note).
  The off-diagonal (radial) structure is approximate: each term keeps the
  shape of the corresponding second-order term, rescaled by a scalar. All
  integrals come straight from the stored tables (the L^k entries with i = v
  are already at the valence energy), so nothing is computed on-the-fly.

  The sub-grid (@p r0, @p rmax, @p stride) defaults match
  Wavefunction::formSigma.

  @param v               Valence state (basis version; must be in the qk/lk
                         tables)
  @param core            Core (hole) orbitals
  @param excited         Excited orbitals
  @param qk              Converged Coulomb \f$ Q^k \f$ table
  @param lk              Converged ladder \f$ L^k \f$ table: supplies the L^k
                         integrals (nullptr gives Sigma_L = 0)
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G       Include the lower (g) component of Sigma_L
  @return Sigma_L as a coordinate-space GMatrix

  @note Terms with \f$ Q^k_{mnva} = 0 \f$ (but \f$ L^k \ne 0 \f$) cannot be
        rescaled and are dropped; inherent to the method. Note that L^k exists
        at both parities of k (see k_minmax_L) while Q^k exists only at
        Coulomb parity, so the exchange-only (wrong-parity) k channels are
        always dropped here.
*/
GMatrix Sigma_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     double r0 = 1.0e-4, double rmax = 30.0,
                     std::size_t stride = 4, bool include_G = false);

/*!
  @brief Ladder correction to the correlation potential, Sigma_L(e_v): exact
         acting on the valence state, by the ratio method otherwise.

  @details
  In the Brueckner equation the change of the valence orbital,
  dF = F^Br - F_v, is small: Sigma_L acts exactly on F_v, and only its
  action on dF is from the ratio method. With S the ratio-method Sigma_L
  [Sigma_ladder()], made symmetric, and the exact column
  \f$ |c_v\rangle = \Sigma_L|v\rangle \f$ (the x = v term of the
  projection onto the basis, from the ladder integrals L^k_mnva and L^k_vnab
  of @p lk):
  \f[
    \Sigma_L = S + |c\rangle\langle v| + |v\rangle\langle c|
               - \langle v|c\rangle\,|v\rangle\langle v| ,
  \f]
  \f[
    |c\rangle = |c_v\rangle - S|v\rangle .
  \f]
  Sigma_L is symmetric, \f$ \Sigma_L|v\rangle = |c_v\rangle \f$, and
  <v|Sigma_L|v> = de_valence_w(v), up to the sub-grid quadrature. The column
  is that of the ladder diagrams as the ladder equation generates them: the
  projection onto the basis [sigma_method basis], made symmetric, acts on v
  with the average of this column and the matching row.

  @param v               Valence state (basis version; normalised; must be in
                         the qk/lk tables)
  @param core            Core (hole) orbitals
  @param excited         Excited orbitals
  @param qk              Converged Coulomb \f$ Q^k \f$ table
  @param lk              Converged ladder \f$ L^k \f$ table
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G       Include the lower (g) component of Sigma_L
  @return Sigma_L as a coordinate-space GMatrix
*/
GMatrix Sigma_ladder_ratio2(const DiracSpinor &v,
                            const std::vector<DiracSpinor> &core,
                            const std::vector<DiracSpinor> &excited,
                            const Coulomb::QkTable &qk,
                            const Coulomb::LkTable &lk, double r0 = 1.0e-4,
                            double rmax = 30.0, std::size_t stride = 4,
                            bool include_G = false);

//==============================================================================

/*!
  @brief Ladder correlation potential matrix, Sigma_L, for a single state.

  @details
  For a single (kappa, n) state, evaluated at energy en.
*/
struct SigmaLData {
  //! Dirac quantum number of the state Sigma_L was formed for
  int kappa;
  //! Principal quantum number of the state Sigma_L was formed for
  int n;
  //! Energy Sigma_L was evaluated at
  double en;
  //! Sigma_L matrix, on its own sub-grid
  GMatrix SL;
};

/*!
  @brief Writes Sigma_L (ladder) matrices to binary file.

  @details
  File contains the full-grid parameters (for checking on read), followed by
  each Sigma_L matrix with its own sub-grid parameters and include_G flag.

  @return True if written; false (and writes nothing) if @p fname is empty or
          'false', or @p SLs is empty.
*/
bool write_SigmaL(const std::string &fname, const std::vector<SigmaLData> &SLs,
                  const Grid &grid);

/*!
  @brief Reads Sigma_L (ladder) matrices from binary file.

  @details
  @p grid must match the full radial grid the matrices were calculated on.
  Each matrix carries its own sub-grid parameters and include_G flag (they
  need not match those of the base Sigma).

  @return The Sigma_L matrices. Empty vector if @p fname is empty or 'false',
          the file doesn't exist, or on grid mismatch.
*/
std::vector<SigmaLData> read_SigmaL(const std::string &fname,
                                    const std::shared_ptr<const Grid> &grid);

// template implementations:
#include "Ladder.ipp"

} // namespace MBPT
