#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/Sigma2.hpp" //temp - remove after refactor
#include <string>
#include <vector>
class DiracDiracSpinor;
namespace MBPT {
class CorrelationPotential;
class StructureRad;
} // namespace MBPT
namespace HF {
class Breit;
}

namespace CI {

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
                   DiracSpinor::Index y, int twoJ);

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
  @return CI Hamiltonian matrix element \f$ H_{AB} \f$.
*/
double Hab(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

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
                 const Coulomb::LkTable &Sk);

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
  @return Filtered basis vector.
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
  @brief Normalisation factor of a CI state,
  \f$ F_X = \matel{X}{\sum_i f(i)}{X} \f$.
  @details
  The states of the CI are normalised in the model space, while the true
  states have amplitude in the core-excited configurations that the CI does
  not span. The physical matrix element is

  \f[
    \redmatel{B}{h}{A}_{\rm phys} = (1 + F_B + F_A)\,\redmatel{B}{h}{A},
  \f]

  and the second-order amplitude, which has one such factor per vertex, is
  \f$ (1 + F_a + F_b + 2F_n) \f$; see @ref A_K.

  \f$ f \f$ is the one-body norm defect, \f$ f_v = -\tfrac12
  \braket{\chi_v}{\chi_v} \f$ (MBPT::StructureRad::f_norm). It is diagonal in
  the single-particle basis, so \f$ \sum_i f(i) \f$ is diagonal in the CSFs,
  with \f$ F_I = f_v + f_w \f$ for \f$ I = \ket{vw} \f$, and

  \f[
    F_X = \sum_I |c^X_I|^2 \, (f_v + f_w).
  \f]

  @param Psi,i  CI state (solution @p i of @p Psi).
  @param f      Table of the one-body norm defect; only the diagonal,
                @p f.getv(v,v), is used.
  @return \f$ F_X \f$.

  @note Every valence electron contributes, the spectators included. Adding
        \f$ f_{v'} + f_v \f$ to each single-particle matrix element instead
        keeps only the two orbitals the operator connects, and so drops the
        spectators: that is exact for one valence electron only.
*/
[[nodiscard]] double norm_factor(const PsiJPi &Psi, std::size_t i,
                                 const Coulomb::meTable<double> &f);

/*!
  @brief Table of the one-body norm defect \f$ f_v \f$, for @ref norm_factor.
  @details
  Fills the diagonal with MBPT::StructureRad::f_norm. States with
  \f$ n > \f$ @p n_max are set to zero: SR+N is meaningful only between
  physical states, and the high-n states of a CI basis are cavity states. Use
  the same @p n_max as for the structure radiation itself.

  @param sr     Structure radiation object; only f_norm() is used, so
                solve_core() need not have been called.
  @param basis  Orbitals to fill the table for; e.g., the CI basis.
  @param n_max  Maximum n for which \f$ f_v \f$ is non-zero.
*/
[[nodiscard]] Coulomb::meTable<double>
f_norm_table(const MBPT::StructureRad &sr,
             const std::vector<DiracSpinor> &basis, int n_max = 999);

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
  @return Full CI Hamiltonian matrix in the CSF basis.
*/
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi,
                                     const Coulomb::meTable<double> &h1,
                                     const Coulomb::QkTable &qk,
                                     const Coulomb::WkTable *Bk = nullptr,
                                     const Coulomb::LkTable *Sk = nullptr);

/*!
  @brief Constructs the CI Hamiltonian matrix from a set of integral tables.
  @details
  Overload of construct_Hci() taking the tables as @ref Integrals; the Breit
  and \f$ \Sigma_2 \f$ corrections are included if those tables are non-empty.

  @param psi   CI solution container holding the CSF basis and J/parity.
  @param ints  Integral tables, e.g., from Wavefunction::CI_integrals().
  @return Full CI Hamiltonian matrix in the CSF basis.
*/
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi, const Integrals &ints);

} // namespace CI
