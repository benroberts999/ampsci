#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/Sigma2.hpp" //temp - remove after refactor
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
class DiracDiracSpinor;
namespace MBPT {
class CorrelationPotential;
} // namespace MBPT
namespace HF {
class Breit;
}

namespace CI {

//==============================================================================
/*!
  @brief Resummed derivative (dSigma/dE) correction to a Sigma_1 matrix
  element (Kozlov formula).
  @details
  Returns the corrected matrix element,

  \f[
    \Sigma \to \Sigma \left[ 1 - \delta E \, (d\Sigma/dE)/\Sigma \right]^{-1},
  \f]

  which resums the linear expansion
  \f$ \Sigma(\epsilon_0 + \delta E) \approx \Sigma + \delta E \, d\Sigma/dE \f$.

  Guard: if the corrected value exceeds \f$ |\Sigma| \f$, then
  \f$ \delta E \, (d\Sigma/dE) \f$ is approaching \f$ \Sigma \f$ - the pole
  of the resummed (Pade) form, at
  \f$ \delta E = \Sigma / (d\Sigma/dE) \f$ - where the
  expression diverges; the correction is distrusted, and the uncorrected
  \f$ \Sigma \f$ is returned.

  @param Sigma   Matrix element \f$ \langle a|\Sigma_1|b\rangle \f$.
  @param dSigma  Energy derivative, \f$ \langle a|d\Sigma_1/dE|b\rangle \f$.
  @param dE      Energy shift \f$ \delta E \f$ from the energy Sigma_1 was
                 evaluated at.
  @return Corrected matrix element.
*/
double corrected_Sigma(double Sigma, double dSigma, double dE);

/*!
  @brief Resummed shift of a \f$ \Sigma_2 \f$ integral to a new reference
  energy E0 (Brillouin-Wigner denominators).
  @details
  \f$ S^k \to (S^k)^2 / (S^k - \delta E_0\, dS^k/dE_0) \f$, the same
  resummation as @ref corrected_Sigma. Exact for a single energy denominator
  (each term is \f$ N/(D + \delta E_0) \f$), so it holds up much better than
  the linear expansion when \f$ \delta E_0 \f$ is a sizeable fraction of the
  denominator.

  Unlike @ref corrected_Sigma there is no guard against the correction
  enhancing \f$ |S^k| \f$: for \f$ \Sigma_2 \f$ the shift legitimately goes
  either way. The only guard is on crossing the pole of the resummed form,
  at \f$ \delta E_0 = S^k / (dS^k/dE_0) \f$ (detected by the denominator
  changing sign), beyond which the expansion is meaningless: the unshifted
  \f$ S^k \f$ is returned there.

  @param Sk    Integral \f$ S^k \f$, evaluated at the reference E0.
  @param dSk   Derivative \f$ dS^k/dE_0 \f$, at the same reference.
  @param dE0   Shift from the reference, \f$ E_0 - E_0^{\rm ref} \f$.
  @return Shifted \f$ S^k \f$.
*/
double corrected_Sk(double Sk, double dSk, double dE0);

//==============================================================================
/*!
  @brief Derivative (dSigma/dE) correction data for the one-body Sigma_1
  matrix elements.
  @details
  Restores the state dependence of Sigma_1, which is otherwise evaluated at a
  fixed energy for each kappa. Each one-body matrix element entering a CI
  matrix element is corrected via @ref corrected_Sigma, with

  \f[
    \delta E = E_0 - E_\Sigma(\kappa) - \epsilon_{\rm spectator},
  \f]

  where \f$ E_0 \f$ is the reference total two-electron valence energy,
  \f$ E_\Sigma(\kappa) \f$ is the energy Sigma_1 was evaluated at, and
  \f$ \epsilon_{\rm spectator} \f$ is the orbital energy of the spectator
  electron in the determinant.

  Fill with @ref calculate_dSdE_correction; applied by Hab() (as a correction
  on top of an h1 table that already includes Sigma_1).

  @note The tables are agnostic to how Sigma_1 is calculated: to use the
        Feynman (all-orders) Sigma, fill S1 and dS1 from the
        CorrelationPotential (formSigma at two energies) instead of
        MBPT::Sigma_vw / MBPT::dSigma_dE_vw.
*/
struct Sigma1Correction {
  //! One-body Sigma_1 matrix elements (uncorrected), <a|Sigma_1|b>
  Coulomb::meTable<double> S1{};
  //! Energy derivative matrix elements, <a|dSigma_1/dE|b>
  Coulomb::meTable<double> dS1{};
  //! Energy Sigma_1 was evaluated at, for each kappa
  std::map<int, double> e_sigma{};
  //! Single-particle orbital energies, keyed by nk_index
  std::map<DiracSpinor::Index, double> en{};
  //! Reference total two-electron valence energy. 0.0 means "not set": it
  //! belongs to a single (J, parity), and is found by @ref iterate_E0
  double E0{0.0};

  [[nodiscard]] bool empty() const { return dS1.empty(); }

  /*!
    @brief Correction to the one-body matrix element <a|h1|b>, given the
    spectator orbital.
    @details
    Returns \f$ \Sigma_{\rm corrected} - \Sigma \f$, i.e., the amount to add
    to an h1 matrix element that already includes (uncorrected) Sigma_1.
    Returns zero for orbitals not in the tables.
  */
  [[nodiscard]] double delta_h1(DiracSpinor::Index a, DiracSpinor::Index b,
                                DiracSpinor::Index spectator) const;

  /*!
    @brief Reads or writes the correction tables to/from a binary file.
    @details
    Stores S1, dS1, e_sigma, and en. E0 is not stored: it belongs to a single
    (J, parity) and is set at solve time (see @ref iterate_E0). On read, the
    existing tables are replaced.

    @note No settings are stored in the file: the filename identifies the
          calculation (cf. @ref PsiJPi::read_write).

    @return True on success; false if the file does not exist or cannot be
            opened (read), or holds no tables.
  */
  bool read_write(const std::string &fname, IO::FRW::RoW rw);
};

/*!
  @brief Builds the Sigma_1 derivative-correction tables; see
  @ref Sigma1Correction.
  @details
  For each same-kappa pair in @p ci_basis, computes the (uncorrected) Sigma_1
  matrix element and its energy derivative (central finite difference of
  MBPT::Sigma_vw). Sigma_1 is evaluated at the energy of the first state of
  each kappa in @p ci_basis - the same convention as calculate_h1_table(), so
  the stored S1 matches the Sigma_1 included in the h1 table.

  @param ci_basis          Basis states for which table entries are needed.
  @param s1_basis_core     Core states used as internal lines for Sigma_1.
  @param s1_basis_excited  Excited states used as internal lines for Sigma_1.
  @param qk                Table of Coulomb \f$ Q^k \f$ integrals.
  @return Filled Sigma1Correction tables (E0 left unset; see @ref iterate_E0).
*/
[[nodiscard]] Sigma1Correction
calculate_dSdE_correction(const std::vector<DiracSpinor> &ci_basis,
                          const std::vector<DiracSpinor> &s1_basis_core,
                          const std::vector<DiracSpinor> &s1_basis_excited,
                          const Coulomb::QkTable &qk);

/*!
  @brief Builds the Sigma_1 derivative-correction tables directly from a
  correlation potential that holds dSigma/dE matrices.
  @details
  Fills S1 from the actual Sigma of the correlation potential
  (via CorrelationPotential::SigmaFv - so it matches the Sigma_1 in the h1
  table exactly, whatever the method: Goldstone, Feynman, all-orders), and
  dS1 from its stored dSigma/dE matrices (see the Correlations option
  `derivative`). Much faster than the qk-table overload (matrix
  applications only), and consistent with the all-orders Sigma.

  The reference energies e_sigma are taken from the stored Sigma data (the
  energy each Sigma matrix was formed at, first entry of each kappa).

  @param ci_basis  Basis states for which table entries are needed.
  @param Sigma     Correlation potential; must hold dSigma/dE matrices
                   (see CorrelationPotential::has_derivative()).
  @return Filled Sigma1Correction tables (E0 left unset; see @ref iterate_E0).

  @note The ladder part (Sigma_L) is included in S1 (via SigmaFv) but has no
        energy derivative, so it is absent from dS1.
*/
[[nodiscard]] Sigma1Correction
calculate_dSdE_correction(const std::vector<DiracSpinor> &ci_basis,
                          const MBPT::CorrelationPotential &Sigma);

/*!
  @brief Finds the reference energy E0 for the dSigma/dE correction, for a
  single (J, parity); returns the CI Hamiltonian built with it.
  @details
  E0 is the lowest energy of this (J, parity), which is only known once we
  have solved - so it is found self-consistently. If @p s1c has no E0 set
  (i.e., 0.0), the iteration starts from the lowest zeroth-order configuration
  energy. The first pass builds the
  CI Hamiltonian at the current E0 and diagonalises it for the lowest level.
  E0 enters only through the (small) dSigma/dE correction, so the state itself
  hardly changes from pass to pass: later passes rebuild the Hamiltonian with
  the updated E0 and take the new E0 as the expectation value of the new
  Hamiltonian in the (unchanged) state - first-order perturbation theory in
  the change to the Hamiltonian. Only the first pass is diagonalised.

  @param psi            Solved for the lowest level (updated in place).
  @param s1c            Correction tables. E0 belongs to a single (J, parity),
                        while these tables are shared between them, so a local
                        copy is made: @p s1c is left as it was.
  @param h1             One-body matrix element table (includes Sigma_1).
  @param qk             Coulomb \f$ Q^k \f$ table.
  @param Bk             Pointer to Breit table; ignored if nullptr.
  @param Sk             Pointer to \f$ \Sigma_2 \f$ table; ignored if nullptr.
  @param hk             Average S^k/Q^k ratios; see @ref MBPT::average_hk.
  @param dSk            Pointer to the dS^k/dE0 table (Brillouin-Wigner
                        Sigma_2); ignored if nullptr. Sigma_2 is shifted to
                        the current E0 each pass, so both corrections move
                        together as E0 converges.
  @param E0_sigma2      Reference E0 that @p Sk and @p dSk were tabulated at.
  @param outstream      Stream for the per-pass output.
  @return CI Hamiltonian matrix, built with the converged E0.
*/
[[nodiscard]] LinAlg::Matrix<double> iterate_E0(
  PsiJPi *psi, const Sigma1Correction &s1c, const Coulomb::meTable<double> &h1,
  const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk = nullptr,
  const Coulomb::LkTable *Sk = nullptr, const std::vector<double> &hk = {},
  const Coulomb::LkTable *dSk = nullptr, double E0_sigma2 = 0.0,
  std::ostream &outstream = std::cout);

//==============================================================================
/*!
  @brief The integral tables required to construct the CI Hamiltonian matrix.
  @details
  Everything needed to construct the CI Hamiltonian for any (J, parity), as it
  was constructed for the CI solutions: the single-particle basis, the one-body
  matrix elements (which may include \f$ \Sigma_1 \f$), and the two-body
  Coulomb, Breit and \f$ \Sigma_2 \f$ tables.

  Filled by @ref configuration_interaction and stored in the Wavefunction (see
  Wavefunction::CI_integrals), so that later calculations can construct CI
  Hamiltonians - e.g., for the mixed-states equation, @ref solve_mixed_state -
  without recalculating any integrals.

  @note The Breit and \f$ \Sigma_2 \f$ tables are empty if those corrections
        were not included; @ref construct_Hci then skips them.

  @note These tables are large (the Coulomb table especially): keeping them for
        the entire run costs memory.
*/
struct Integrals {
  //! Single-particle basis used for the CI expansion
  std::vector<DiracSpinor> ci_basis{};
  //! One-body matrix elements, <a|h1|b>; may include Sigma_1
  Coulomb::meTable<double> h1{};
  //! Two-body Coulomb integrals, Q^k
  Coulomb::QkTable qk{};
  //! Two-body Breit integrals, B^k; empty if not included
  Coulomb::WkTable Bk{};
  //! Two-body Sigma_2 integrals, S^k; empty if not included
  Coulomb::LkTable Sk{};
  //! Average S^k/Q^k ratios, indexed by k; empty if Sigma_2 is not being
  //! extrapolated beyond the cis2 basis. Diagrams with no stored S^k then
  //! use S^k = hk[k] * Q^k. See MBPT::average_hk
  std::vector<double> hk{};
  //! Energy derivatives dS^k/dE0 of the Sigma_2 integrals; empty unless the
  //! Brillouin-Wigner (E0-dependent) Sigma_2 correction is included.
  //! See CI::corrected_Sk
  Coulomb::LkTable dSk{};
  //! Reference E0 that Sk and dSk were tabulated at (Brillouin-Wigner)
  double E0_sigma2{0.0};
  //! Derivative (dSigma/dE) correction for Sigma_1; empty if not included
  Sigma1Correction s1_corr{};

  //! False if the tables were never calculated (e.g., CI was run 'read_only')
  [[nodiscard]] bool availableQ() const {
    return !ci_basis.empty() && !h1.empty();
  }
};

//==============================================================================
/*!
  @brief The result of a CI calculation: the solutions, and the integrals used
  to construct the CI Hamiltonian.
  @details
  Returned by @ref configuration_interaction, and stored in the Wavefunction
  (see Wavefunction::CIwfs and Wavefunction::CI_integrals). The integrals are
  kept so that CI Hamiltonians for other (J, parity) may be constructed later
  without recalculating them - e.g., for the mixed-states equation,
  @ref solve_mixed_state. See @ref Integrals.
*/
struct Solutions {
  //! One entry per {J, parity} requested
  std::vector<PsiJPi> levels{};
  //! Integral tables used to construct the CI Hamiltonians
  Integrals integrals{};
};

/*!
  @brief Antisymmetrised two-body Coulomb matrix element in the coupled CSF
  basis.
  @details
  Evaluates the angular-reduced, antisymmetrised Coulomb interaction between
  two two-electron CSFs \f$ |vw; J\rangle \f$ and \f$ |xy; J\rangle \f$:

  \f[
    \langle vw; J \| g \| xy; J \rangle
    = \eta_{vw}\eta_{xy}
      \sum_k (-1)^{j_v+j_x+k+J}
      \begin{Bmatrix} j_v & j_w & J \\ j_y & j_x & k \end{Bmatrix}
      Q^k_{vwxy} + \text{exchange},
  \f]

  where \f$ \eta_{ab} = 1/\sqrt{2} \f$ if \f$ a = b \f$ (identical-particle
  normalisation) and 1 otherwise, and \f$ Q^k \f$ are the Coulomb integrals stored in @p qk.

  @param qk    Table of Coulomb \f$ Q^k \f$ integrals.
  @param v,w   Indices of the bra single-particle states.
  @param x,y   Indices of the ket single-particle states.
  @param twoJ  Twice the total angular momentum 2J of the coupled pair.
  @return Antisymmetrised, angular-reduced two-body Coulomb matrix element.
*/
double CSF2_Coulomb(const Coulomb::QkTable &qk, DiracSpinor::Index v,
                    DiracSpinor::Index w, DiracSpinor::Index x,
                    DiracSpinor::Index y, int twoJ);

/*!
  @brief Two-body \f$ \Sigma_2 \f$ (MBPT) correction to CSF2_Coulomb().
  @details
  Evaluates the same angular reduction as CSF2_Coulomb(), but using the
  two-body \f$ \Sigma_2 \f$ integrals \f$ S^k \f$ stored in @p Sk in place of
  the Coulomb \f$ Q^k \f$ integrals.  Adds the second-order MBPT correction to
  the two-electron interaction.

  @param Sk    Table of two-body \f$ \Sigma_2 \f$ (\f$ L^k \f$) integrals.
  @param v,w   Indices of the bra single-particle states.
  @param x,y   Indices of the ket single-particle states.
  @param twoJ  Twice the total angular momentum 2J of the coupled pair.
  @return Antisymmetrised two-body \f$ \Sigma_2 \f$ matrix element.
*/
double CSF2_Sigma2(const Coulomb::LkTable &Sk, DiracSpinor::Index v,
                   DiracSpinor::Index w, DiracSpinor::Index x,
                   DiracSpinor::Index y, int twoJ,
                   const Coulomb::QkTable *qk = nullptr,
                   const std::vector<double> &hk = {},
                   const Coulomb::LkTable *dSk = nullptr, double dE0 = 0.0);

/*!
  @brief Antisymmetrised two-body Breit matrix element in the coupled CSF
  basis.
  @details
  Evaluates the same angular reduction as CSF2_Coulomb(), but using the Breit
  \f$ B^k \f$ integrals stored in @p Bk.

  @param Bk    Table of Breit \f$ W^k \f$ integrals.
  @param v,w   Indices of the bra single-particle states.
  @param x,y   Indices of the ket single-particle states.
  @param twoJ  Twice the total angular momentum 2J of the coupled pair.
  @return Antisymmetrised two-body Breit matrix element.
*/
double CSF2_Breit(const Coulomb::WkTable &Bk, DiracSpinor::Index v,
                  DiracSpinor::Index w, DiracSpinor::Index x,
                  DiracSpinor::Index y, int twoJ);

/*!
  @brief CI Hamiltonian matrix element between two two-electron CSFs.
  @details
  Computes \f$ H_{AB} = \langle A | \hat{H} | B \rangle \f$ using the
  Slater-Condon rules, including one-body terms from @p h1 (which may already
  incorporate \f$ \Sigma_1 \f$ corrections) and the two-body Coulomb
  interaction via CSF2_Coulomb().

  Does NOT include \f$ \Sigma_2 \f$ or Breit corrections; add those via
  Sigma2_AB() and Breit_AB() respectively.

  @param A,B   The two CSFs.
  @param twoJ  Twice the total angular momentum 2J.
  @param h1    Table of one-body matrix elements \f$ \langle a | h_1 | b \rangle \f$.
  @param qk    Table of Coulomb \f$ Q^k \f$ integrals.
  @param s1c   Optional derivative (dSigma/dE) correction to Sigma_1; applied
               to each one-body matrix element (with the spectator orbital
               energy) if given. See @ref Sigma1Correction.
  @return CI Hamiltonian matrix element \f$ H_{AB} \f$.
*/
double Hab(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *s1c = nullptr);

/*!
  @brief Two-body \f$ \Sigma_2 \f$ correction to Hab().
  @details
  Evaluates the MBPT \f$ \Sigma_2 \f$ contribution to the CI matrix element
  using CSF2_Sigma2().  Add to Hab() to form the full CI+MBPT Hamiltonian
  matrix element.

  @param A,B   The two CSFs.
  @param twoJ  Twice the total angular momentum 2J.
  @param Sk    Table of \f$ \Sigma_2 \f$ (\f$ L^k \f$) integrals.
  @return \f$ \Sigma_2 \f$ correction to \f$ H_{AB} \f$.
*/
double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::LkTable &Sk,
                 const Coulomb::QkTable *qk = nullptr,
                 const std::vector<double> &hk = {},
                 const Coulomb::LkTable *dSk = nullptr, double dE0 = 0.0);

/*!
  @brief Breit correction to Hab().
  @details
  Evaluates the two-body Breit contribution to the CI matrix element using
  CSF2_Breit().  Add to Hab() to include the Breit interaction.

  @param A,B   The two CSFs.
  @param twoJ  Twice the total angular momentum 2J.
  @param Bk    Table of Breit \f$ W^k \f$ integrals.
  @return Breit correction to \f$ H_{AB} \f$.
*/
double Breit_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                const Coulomb::WkTable &Bk);

/*!
  @brief Builds the one-body Hamiltonian matrix element table for the CI basis.
  @details
  Constructs a lookup table of single-particle matrix elements
  \f$ \langle a | h_1 | b \rangle \f$ for all pairs \f$ a, b \f$ in
  @p ci_basis.  The diagonal elements are the HF single-particle energies.

  If @p include_Sigma1 is true, the one-body MBPT \f$ \Sigma_1 \f$ correction
  is computed from the Coulomb integrals in @p qk using @p s1_basis_core and
  @p s1_basis_excited as the internal lines of the MBPT diagrams and added to
  the diagonal.

  @param ci_basis          Basis states for which table entries are needed.
  @param s1_basis_core     Core states used as internal lines for \f$ \Sigma_1 \f$.
  @param s1_basis_excited  Excited states used as internal lines for \f$ \Sigma_1 \f$.
  @param qk                Table of Coulomb \f$ Q^k \f$ integrals.
  @param include_Sigma1    If true, add one-body MBPT \f$ \Sigma_1 \f$ corrections.
  @return Table of \f$ \langle a | h_1 | b \rangle \f$ matrix elements.

  @warning Assumes @p ci_basis states are Hartree-Fock eigenstates, so
           off-diagonal HF terms vanish.
*/
[[nodiscard]] Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const std::vector<DiracSpinor> &s1_basis_core,
                   const std::vector<DiracSpinor> &s1_basis_excited,
                   const Coulomb::QkTable &qk, bool include_Sigma1);

/*!
  @brief Builds the one-body Hamiltonian table using a precomputed
  CorrelationPotential.
  @details
  Overload of calculate_h1_table() that uses a CorrelationPotential object
  (i.e., a precomputed \f$ \Sigma_1 \f$ operator) instead of computing MBPT
  diagrams on the fly.  Preferred when a CorrelationPotential is available, as
  it is generally faster and more complete.

  @param ci_basis       Basis states for which table entries are needed.
  @param Sigma          Precomputed one-body correlation potential \f$ \Sigma_1 \f$.
  @param include_Sigma1 If true, include \f$ \Sigma_1 \f$ corrections from @p Sigma.
  @return Table of \f$ \langle a | h_1 | b \rangle \f$ matrix elements.
*/
[[nodiscard]] Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const MBPT::CorrelationPotential &Sigma,
                   bool include_Sigma1);

/*!
  @brief Builds or loads the two-body Breit integral table.
  @details
  Computes Breit \f$ W^k \f$ integrals for all pairs in @p ci_basis using the
  Breit operator @p pBr.  Results are cached to/from @p bk_filename.

  If @p pBr is nullptr or @p no_new_integralsQ is true, no new integrals are
  computed; only cached values are loaded.

  @param bk_filename      Filename for caching the \f$ W^k \f$ table.
  @param pBr              Pointer to Breit operator; if nullptr, returns empty table.
  @param ci_basis         Basis for which Breit integrals are needed.
  @param max_k            Maximum multipolarity k to include.
  @param no_new_integralsQ If true, skip computing any new integrals.
  @return Table of Breit \f$ W^k \f$ integrals.
*/
[[nodiscard]] Coulomb::WkTable
calculate_Bk(const std::string &bk_filename, const HF::Breit *const pBr,
             const std::vector<DiracSpinor> &ci_basis, int max_k,
             bool no_new_integralsQ = false);

/*!
  @brief Returns the subset of @p basis matching @p include_str, excluding
  states in @p exclude_str.
  @details
  Filters @p basis to retain only states described by the ampsci basis-string
  notation (e.g., "20spdf") that are not part of the frozen core.

  @param basis               Full single-particle basis to filter.
  @param include_str         Basis-string specifying which states to keep;
                             if empty, all states in @p basis are kept (subject
                             to the frozen-core exclusion).
  @param exclude_str         Basis-string specifying core states to exclude.
  @return Filtered basis vector. Negative-energy states (if present) with
          kappa in @p include_str are kept (the n limit is not applied to them).
*/
[[nodiscard]] std::vector<DiracSpinor>
basis_subset(const std::vector<DiracSpinor> &basis,
             const std::string &include_str,
             const std::string &exclude_str = "");

/*!
  @brief Reduced matrix element between two CI states (low-level overload).
  @details
  Evaluates the reduced matrix element of a rank-@p K_rank tensor operator
  between two CI states:

  \f[
    \redmatel{A}{T^K}{B}
    = \sum_{ij} c_i^A \, c_j^B \, \redmatel{\text{CSF}_i}{T^K}{\text{CSF}_j},
  \f]

  where the single-particle reduced matrix elements are looked up from @p h.

  @param cA,cB      CI expansion coefficient vectors for states A and B.
  @param CSFAs,CSFBs CSF bases for states A and B respectively.
  @param twoJA,twoJB Twice the total angular momentum of states A and B.
  @param h          Lookup table of single-particle reduced matrix elements.
  @param K_rank     Rank of the tensor operator.
  @param Parity     Parity of the operator (+1 or -1).
  @return Reduced matrix element \f$ \redmatel{A}{T^K}{B} \f$.
*/
double ReducedME(const LinAlg::View<const double> &cA,
                 const std::vector<CI::CSF2> &CSFAs, int twoJA,
                 const LinAlg::View<const double> &cB,
                 const std::vector<CI::CSF2> &CSFBs, int twoJB,
                 const Coulomb::meTable<double> &h, int K_rank, int Parity);

/*!
  @brief Reduced matrix element between two CI states (PsiJPi overload).
  @details
  Convenience wrapper around the low-level ReducedME() overload.  Extracts
  expansion coefficients and CSF lists from @p As and @p Bs for the requested
  solution indices @p iA and @p iB.

  @param As,Bs   CI solution containers for the two states.
  @param iA,iB  Solution indices within @p As and @p Bs.
  @param h      Lookup table of single-particle reduced matrix elements.
  @param K_rank Rank of the tensor operator.
  @param Parity Parity of the operator (+1 or -1).
  @return Reduced matrix element \f$ \redmatel{A}{T^K}{B} \f$.
*/
inline double ReducedME(const PsiJPi &As, std::size_t iA, const PsiJPi &Bs,
                        std::size_t iB, const Coulomb::meTable<double> &h,
                        int K_rank, int Parity) {
  return ReducedME(As.coefs(iA), As.CSFs(), As.twoJ(), Bs.coefs(iB), Bs.CSFs(),
                   Bs.twoJ(), h, K_rank, Parity);
}

/*!
  @brief Reduced matrix element between two two-electron CSFs.
  @details
  Evaluates \f$ \redmatel{X; J_X}{T^K}{V; J_V} \f$ for a rank-@p K_rank
  one-body tensor operator using the standard 6j angular reduction, accounting
  for identical-particle normalisation factors.

  @warning This function may not handle all cases correctly; results should be
           verified for non-trivial configurations.
*/
double RME_CSF2(const CI::CSF2 &X, int twoJX, const CI::CSF2 &V, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank);

/*!
  @brief Leading non-relativistic configuration of a CI state.
  @details
  Sums |c|^2 over the relativistic CSFs belonging to each non-relativistic
  configuration, and returns the configuration with the largest total weight.

  @param coefs  CI expansion coefficients (one per CSF).
  @param csfs   The CSF basis (matching @p coefs).
  @return Pair {configuration label (non-rel notation), total |c|^2 weight}.
*/
std::pair<std::string, double>
leading_config(const LinAlg::View<const double> &coefs,
               const std::vector<CSF2> &csfs);

/*!
  @brief Determines the best-fit (S, L) term for a two-electron state by
  matching the g-factor.
  @details
  Iterates over all allowed (S, L) combinations for given orbital angular
  momenta @p l1, @p l2 and total @p twoJ /2, and returns the pair whose
  Lande g-factor is closest to @p gJ_target.

  @param l1,l2      Orbital angular momenta of the two electrons.
  @param twoJ       Twice the total angular momentum 2J.
  @param gJ_target  Target g-factor to match.
  @return Best-fit {2S, L} pair.
*/
std::pair<int, int> Term_S_L(int l1, int l2, int twoJ, double gJ_target);

/*!
  @brief Determines the (S, L) term for a two-electron state from the
  expectation values of L^2 and S^2.
  @details
  Returns the (S, L) pair, subject to the triangle condition with J, that
  minimises the combined distance |L(L+1) - @p L2| + |S(S+1) - @p S2|.
  For a mixed state, this is the nearest (dominant) term; the purity is
  indicated by @p L2, @p S2 themselves. See @ref expectation_L2S2.

  @param L2    Expectation value of L^2.
  @param S2    Expectation value of S^2.
  @param twoJ  Twice the total angular momentum 2J.
  @return Best-fit {S, L} pair.
*/
std::pair<int, int> Term_S_L_from_expectation(double L2, double S2, int twoJ);

//! Returns spectroscopic term symbol string, e.g. "3P_1"
std::string Term_Symbol(int two_J, int L, int two_S, int parity);
//! Returns term symbol without the J subscript, e.g. "3P"
std::string Term_Symbol(int L, int two_S, int parity);

/*!
  @brief Constructs the full CI Hamiltonian matrix in the CSF basis.
  @details
  Builds the symmetric matrix \f$ H_{AB} \f$ for all CSF pairs in @p psi,
  calling Hab() for each element and optionally adding Breit and
  \f$ \Sigma_2 \f$ corrections.

  @param psi   CI solution container holding the CSF basis and J/parity.
  @param h1    One-body matrix element table (may include \f$ \Sigma_1 \f$).
  @param qk    Coulomb \f$ Q^k \f$ integral table.
  @param Bk    Pointer to Breit \f$ W^k \f$ table; ignored if nullptr.
  @param Sk    Pointer to \f$ \Sigma_2 \f$ \f$ L^k \f$ table; ignored if nullptr.
  @param s1c   Pointer to derivative (dSigma/dE) correction for Sigma_1;
               ignored if nullptr. See @ref Sigma1Correction.
  @param hk    Average S^k/Q^k ratios; if non-empty, diagrams with no stored
               S^k use S^k = hk[k]*Q^k. See @ref MBPT::average_hk.
  @param dSk   Pointer to the dS^k/dE0 table; ignored if nullptr. If given,
               each stored S^k is shifted to the target E0. See
               @ref CI::corrected_Sk.
  @param dE0   Shift of the target E0 from the one Sk was tabulated at.
  @return Full CI Hamiltonian matrix in the CSF basis.
*/
LinAlg::Matrix<double>
construct_Hci(const PsiJPi &psi, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk = nullptr,
              const Coulomb::LkTable *Sk = nullptr,
              const Sigma1Correction *s1c = nullptr,
              const std::vector<double> &hk = {},
              const Coulomb::LkTable *dSk = nullptr, double dE0 = 0.0);

/*!
  @brief Constructs the CI Hamiltonian matrix from a set of integral tables.
  @details
  Overload of construct_Hci() taking the tables as @ref Integrals; the Breit
  and \f$ \Sigma_2 \f$ corrections are included if those tables are non-empty.

  If the iterative (dSigma/dE) correction is included, the reference energy
  E0 (which also sets the Sigma_2 BW shift, if present) belongs to this
  (J, parity) block: the lowest solved energy of @p psi is used if it has
  solutions; otherwise the stored fallback is used, which is the
  ground-state energy of the CI run. Blocks never solved directly - e.g.,
  the target block of a mixed state - thus have their corrections evaluated
  at the run's ground-state energy; the difference is higher order.

  @param psi   CI solution container holding the CSF basis and J/parity.
  @param ints  Integral tables, e.g., from Wavefunction::CI_integrals().
  @return Full CI Hamiltonian matrix in the CSF basis.
*/
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi, const Integrals &ints);

} // namespace CI
