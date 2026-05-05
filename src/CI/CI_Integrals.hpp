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
}
namespace HF {
class Breit;
}

namespace CI {

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
  @brief Builds or loads the two-body \f$ \Sigma_2 \f$ integral table.
  @details
  Computes the two-body MBPT \f$ \Sigma_2 \f$ matrix elements
  \f$ \langle vw | \Sigma_2 | xy \rangle \f$ for all pairs in @p cis2_basis,
  using @p s2_basis_core and @p s2_basis_excited as internal lines.  Results
  are cached to/from @p filename.

  If @p no_new_integrals is true, no new integrals are computed; only
  previously cached values are used.

  @param filename                Filename for caching the \f$ S^k \f$ table.
  @param cis2_basis              Basis for which \f$ \Sigma_2 \f$ matrix elements are needed.
  @param s2_basis_core           Core states used as internal lines.
  @param s2_basis_excited        Excited states used as internal lines.
  @param qk                      Coulomb integral table.
  @param max_k                   Maximum multipolarity k to include.
  @param exclude_wrong_parity_box If true, exclude box diagrams with "wrong" parity.
  @param denominators            Energy denominators to use in MBPT (RS, Fermi, Fermi0).
  @param no_new_integrals        If true, skip computing any new integrals.
  @return Table of two-body \f$ \Sigma_2 \f$ (\f$ L^k \f$) integrals.
*/
[[nodiscard]] Coulomb::LkTable calculate_Sk(
  const std::string &filename, const std::vector<DiracSpinor> &cis2_basis,
  const std::vector<DiracSpinor> &s2_basis_core,
  const std::vector<DiracSpinor> &s2_basis_excited, const Coulomb::QkTable &qk,
  int max_k, bool exclude_wrong_parity_box, MBPT::Denominators denominators,
  bool no_new_integrals = false);

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
  @brief Returns the subset of @p basis matching @p subset_string, excluding
  states in @p frozen_core_string.
  @details
  Filters @p basis to retain only states described by the ampsci basis-string
  notation (e.g., "20spdf") that are not part of the frozen core.

  @param basis               Full single-particle basis to filter.
  @param subset_string       Basis-string specifying which states to keep.
  @param frozen_core_string  Basis-string specifying core states to exclude.
  @return Filtered basis vector.
*/
[[nodiscard]] std::vector<DiracSpinor>
basis_subset(const std::vector<DiracSpinor> &basis,
             const std::string &subset_string,
             const std::string &frozen_core_string = "");

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

} // namespace CI
