#pragma once
#include "IO/FRW_fileReadWrite.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <iostream>
#include <optional>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace CI {

//------------------------------------------------------------------------------
/*!
  @brief Base class for configuration state functions (CSFs).

  @details
  A CSF is an antisymmetrised N-electron basis state with definite total
  angular momentum (\f$ J^2 \f$, \f$ J_z \f$) and parity \f$\pi\f$.
  Each CSF is built from a combination of \f$ N \f$-electron Slater determinants,
  which are themselves built from N single-particle relativistic orbitals.
  Each Slater determinant has a specific set of single-particle angular
  momentum projections \f$ j_{i,z} = m_i \f$, with \f$ \sum_i m_i = J_z \f$:

  \f[
    |I; J J_z\rangle = \sum_{\{m\}} d_{\{m\}} \,
    |n_1\kappa_1 m_1, n_2\kappa_2 m_2, \ldots, n_N\kappa_N m_N\rangle,
  \f]

  where the kets on the right are (normalised) Slater determinants, and

  \f[
    \sum_{\{m\}} d_{\{m\}}^2 = 1.
  \f]

  The parity is the product from the orbitals,

  \f[
    \pi = \prod_i (-1)^{l_i}.
  \f]

  The determinants are not formed explicitly: the coefficients
  \f$ d_{\{m\}} \f$ are fixed by angular momentum coupling, so matrix elements
  are evaluated directly in the coupled basis.

  Specialised by @ref CSF2 (two electrons) and @ref CSF3 (three electrons).
*/
class CSF {
protected:
  int m_parity{1};

  explicit CSF(int parity) : m_parity(parity) {}
  CSF(const CSF &) = default;
  CSF(CSF &&) = default;
  CSF &operator=(const CSF &) = default;
  CSF &operator=(CSF &&) = default;

public:
  virtual ~CSF() = default;

  //! Parity of the CSF, +/-1
  int parity() const { return m_parity; }

  //! Number of electrons in the CSF
  virtual std::size_t num_electrons() const = 0;

  //! Single-particle configuration as a string, in relativistic or non-rel form
  virtual std::string config(bool relativistic = false) const = 0;
};

//------------------------------------------------------------------------------
/*!
  @brief Two-electron configuration state function (CSF).

  @details
  A CSF is an antisymmetrised two-electron basis state with definite total
  angular momentum (\f$ J^2 \f$, \f$ J_z \f$) and parity,
  built from a pair of single-particle
  relativistic orbitals. See @ref CSF3 for three electrons.

  For orbitals v and w, the determinant coefficients of @ref CSF are
  simply Clebcsh-Gordan coefficients:
  \f[
    |vw; J J_z\rangle = \eta_{vw} \sum_{m_v m_w}
    C^{J J_z}_{j_v m_v, j_w m_w} \, |n_v\kappa_v m_v, n_w\kappa_w m_w\rangle,
  \f]
  where \f$ \eta_{vw} = 1/\sqrt{2} \f$ if v = w (both electrons in the same
  orbital), and 1 otherwise. 
  
  The orbitals are put in sorted order (v before w); this fixes the phase
  convention of all CSF2 matrix elements (e.g., CSF2_Coulomb()).

  Each CSF2 stores the indices (@ref DiracSpinor::Index) of its two constituent orbitals (always sorted to avoid double-counting) and the total parity, 
  which is the product of the parities of the two single-particle states.

  @note The orbital pair is stored as a sorted array of DiracSpinor::Index
        (uint16_t) rather than DiracSpinor references, so CSF2 objects are
        cheap to copy and store. There is a limit to maximum n<=256 - see @ref Angular::nk_to_index
*/
class CSF2 : public CSF {

public:
  // nb: array of states is always sorted
  std::array<DiracSpinor::Index, 2> states;

  CSF2(const DiracSpinor &a, const DiracSpinor &b);

  //! Index (nk_index) of the ith constituent orbital (i = 0 or 1)
  DiracSpinor::Index state(std::size_t i) const;

  friend bool operator==(const CSF2 &A, const CSF2 &B);
  friend bool operator!=(const CSF2 &A, const CSF2 &B);

  /*!
    @brief Returns the number of orbitals that differ between two CSFs (0, 1,
    or 2).
    @details
    Used to select the appropriate Slater-Condon rule when evaluating CI matrix
    elements: 0 -- diagonal; 1 -- single substitution; 2 -- double
    substitution; >2 -- zero by orthogonality.
  */
  static int num_different(const CSF2 &A, const CSF2 &B);

  /*!
    @brief For two CSFs differing by exactly one orbital, returns {n, a} where
    @p V contains orbital n and @p X contains orbital a.
    @details
    Identifies the "particle" index n (in @p V but not @p X) and the "hole"
    index a (in @p X but not @p V), as needed to apply the single-substitution
    Slater-Condon rule: \f$ \langle V | \hat{O} | X \rangle \f$ where
    \f$ |V\rangle = \hat{a}^\dag_n \hat{a}_a |X\rangle \f$.

    @warning Result is undefined if @p V and @p X do not differ by exactly one
             orbital; check with num_different() first.
  */
  static std::array<DiracSpinor::Index, 2> diff_1_na(const CSF2 &V,
                                                     const CSF2 &X);

  /*!
    @brief Returns the orbital index shared by two CSFs that differ by exactly
    one orbital.
    @details
    Extracts the common (spectator) orbital needed for single-substitution
    matrix elements.

    @warning Assumes @p A and @p B differ by exactly one orbital.
  */
  static DiracSpinor::Index same_1_j(const CSF2 &A, const CSF2 &B);

  std::size_t num_electrons() const override { return 2; }

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const override;
};

//------------------------------------------------------------------------------
/*!
  @brief Forms all two-electron CSFs with given total J and parity.
  @details
  Iterates over all pairs of single-particle states in @p cisp_basis and
  retains those whose angular momenta can be coupled to total \f$ J = \f$
  @p twoJ /2 and whose combined parity equals @p parity.  Duplicate pairs are
  excluded by construction.

  @param twoJ       Twice the total angular momentum 2J.
  @param parity     Total parity: +1 (even) or -1 (odd).
  @param cisp_basis Single-particle basis from which CSFs are constructed.
  @return Sorted list of all valid two-electron CSFs for the given J and parity.
*/
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis);

//------------------------------------------------------------------------------
/*!
  @brief Three-electron configuration state function (CSF).

  @details
  An antisymmetrised three-electron state with definite total angular momentum
  J, J_z, and parity, built from three single-particle orbitals (which may
  coincide).

  The determinant coefficients of @ref CSF are not stored directly. Instead,
  the CSF is stored as an expansion over "parent" states, which are
  antisymmetric in electrons 1 and 2 only: a two-electron CSF \f$ [xy] \f$
  (exactly as @ref CSF2) with angular momentum \f$ J_2 \f$, coupled to the
  orbital z of the third electron:
  \f[
    |[xy] J_2, z; J J_z\rangle = \sum_{M_2 m_z}
    C^{J J_z}_{J_2 M_2, j_z m_z} \, |xy; J_2 M_2\rangle_{12} \,
    |n_z\kappa_z m_z\rangle_3 ,
  \f]

  where \f$ |xy; J_2 M_2\rangle_{12} \f$ is the CSF2 state of 
  electrons 1 and 2. 
  
  The three-electron CSF is then:
  
  \f[
    |\Psi; J J_z\rangle = \sum_p C_p \, |[x_p y_p] J_{2p}, z_p; J J_z\rangle,
  \f]

  normalised such that

  \f[
    \sum_p C_p^2 = 1.
  \f]

  The parent coefficients \f$ C_p \f$ are chosen such that
  \f$ |\Psi\rangle \f$ is fully antisymmetric (see @ref form_CSF3s); together
  with the Clebsch-Gordan coefficients above, they fix the determinant
  coefficients \f$ d_{\{m\}} \f$. (In the literature, these are called
  "coefficients of fractional parentage": the three-electron state is built
  by adding one electron to each two-electron parent, and \f$ C_p \f$ is the
  amplitude of parent p.)

  Expanding over parents reduces three-electron matrix elements to sums of
  two-electron (or one-electron) ones. E.g., for a two-body operator,
  \f[
    \langle\Psi|\sum_{i<j} g_{ij}|\Psi'\rangle
    = 3\sum_{pq} C_p C'_q \, \delta_{z_p z_q}\delta_{J_{2p}J_{2q}}
    \langle [x_p y_p] J_{2p}|g|[x_q y_q] J_{2q}\rangle .
  \f]

  The expansion depends on J, so (as for @ref CSF2) each CSF belongs to a
  single J. Formed by @ref form_CSF3s.

  @note Orbitals are stored as DiracSpinor::Index (sorted), as in @ref CSF2.
*/
class CSF3 : public CSF {
public:
  //! One parent state of a CSF3, and its coefficient
  struct Parent {
    //! Orbitals of the parent two-electron CSF (sorted, as in CSF2)
    std::array<DiracSpinor::Index, 2> pair;
    //! Orbital of the third electron
    DiracSpinor::Index spectator;
    //! Twice the angular momentum of the parent pair, 2*J2
    int twoJ2;
    //! Coefficient of this parent, C_p
    double coef;
  };

  //! The three orbitals; always sorted
  std::array<DiracSpinor::Index, 3> states;

private:
  int m_twoJ12;
  std::vector<Parent> m_parents;

public:
  /*!
    @brief Constructs from the three orbitals and the parent expansion.
    @param orbitals  The three orbitals (any order; stored sorted).
    @param twoJ12    Label: 2*J12 of the coupling the CSF was formed from (see
                     @ref form_CSF3s).
    @param parents   Parent states and their coefficients.
  */
  CSF3(const std::array<DiracSpinor::Index, 3> &orbitals, int twoJ12,
       std::vector<Parent> parents);

  //! Index (nk_index) of the ith constituent orbital (i = 0, 1, or 2)
  DiracSpinor::Index state(std::size_t i) const;

  //! Label: 2*J12 of the two-electron coupling the CSF was formed from
  int twoJ12() const;

  //! Parent states and their coefficients
  const std::vector<Parent> &parents() const;

  /*!
    @brief Returns the number of orbitals that differ between two CSFs (0 to
    3), counting repeated orbitals.
    @details
    Matrix elements of one- and two-body operators vanish if more than two
    orbitals differ.
  */
  static int num_different(const CSF3 &A, const CSF3 &B);

  std::size_t num_electrons() const override { return 3; }

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const override;

  friend bool operator==(const CSF3 &A, const CSF3 &B);
  friend bool operator!=(const CSF3 &A, const CSF3 &B);
};

//------------------------------------------------------------------------------
/*!
  @brief Exchange of electrons 2 and 3, between two parent states of a
  CSF3.
  @details
  Matrix element
  \f$ \langle [x'y']J_2', z'; J | P_{23} | [xy]J_2, z; J \rangle \f$, where
  each state is antisymmetric in electrons 1 and 2 (see @ref CSF3::Parent). In
  terms of uncoupled-pair states, the exchange is the recoupling

  \f[
    P_{23}|(o_1 o_2)J_{12}, o_3; J\rangle
    = \sum_K (-1)^{j_3+j_2+K+J_{12}}\sqrt{[K][J_{12}]}
      \begin{Bmatrix} j_3 & j_1 & K \\ j_2 & J & J_{12} \end{Bmatrix}
      |(o_1 o_3)K, o_2; J\rangle .
  \f]

  @param bra,ket  Parent states (their coef members are not used).
  @param twoJ     Twice the total angular momentum 2J.
  @return The matrix element.
*/
double parent_exchange_23(const CSF3::Parent &bra, const CSF3::Parent &ket,
                          int twoJ);

/*!
  @brief Forms all three-electron CSFs of a single configuration, with given
  total J.
  @details
  Builds an orthonormal basis for the antisymmetric three-electron states of
  the configuration, each as its expansion over parent states (see
  @ref CSF3):

  1. Parent basis: all \f$ |[xy]J_2, z; J\rangle \f$, for each distinct choice
     of orbital z for the third electron, with [xy] the remaining pair. J2 is
     even if x = y, and is limited by the triangle rules. These are
     orthonormal, and span the states antisymmetric in electrons 1 and 2.
  2. On that space, the full antisymmetriser is
     \f[
       A = \frac{1}{3}(1 - 2P_{23})
     \f]
     (see @ref parent_exchange_23); each column \f$ A|p\rangle \f$ is fully
     antisymmetric.
  3. Gram-Schmidt orthonormalisation of the \f$ A|p\rangle \f$, starting with
     the standard couplings ((ab)J12, c)J (or ((bc)J12, a)J if only b = c);
     linearly dependent vectors are discarded.

  For three distinct orbitals this gives one CSF per allowed J12; for
  \f$ j^2 j' \f$, one per even J12; for \f$ j^3 \f$, the standard
  \f$ j^3 \f$ states, including repeated J (e.g., the two J = 9/2 states of
  \f$ (9/2)^3 \f$).

  @param orbitals  The three orbitals (nk_index; any order; may coincide).
  @param twoJ      Twice the total angular momentum 2J.
  @return The CSFs of this configuration and J; empty if there are none.
*/
std::vector<CSF3> form_CSF3s(std::array<DiracSpinor::Index, 3> orbitals,
                             int twoJ);

/*!
  @brief Forms all three-electron CSFs with given total J and parity.
  @details
  Iterates over all distinct triples (with repetition) of single-particle
  states in @p cisp_basis with total parity @p parity, and forms the CSFs of
  each; see @ref form_CSF3s(std::array<DiracSpinor::Index, 3>, int).

  @param twoJ       Twice the total angular momentum 2J (odd).
  @param parity     Total parity: +1 (even) or -1 (odd).
  @param cisp_basis Single-particle basis from which CSFs are constructed.
  @return All three-electron CSFs for the given J and parity.
*/
std::vector<CSF3> form_CSF3s(int twoJ, int parity,
                             const std::vector<DiracSpinor> &cisp_basis);

//==============================================================================
/*!
  @brief jj -> LS recoupling amplitude for an antisymmetrised two-electron CSF.
  @details
  Returns the amplitude of the antisymmetrised jj-coupled CSF
  \f$ |\{(n_1 l_1 j_1)(n_2 l_2 j_2)\}; J\rangle \f$ (orbitals in stored,
  i.e., sorted, order) onto the antisymmetrised LS-coupled state
  \f$ |\{(n_1 l_1)(n_2 l_2)\} L S; J\rangle \f$ of the same non-relativistic
  configuration:

  \f[
    A(L,S) = \eta \sqrt{[j_1][j_2][L][S]}
    \begin{Bmatrix} l_1 & l_2 & L \\ 1/2 & 1/2 & S \\ j_1 & j_2 & J \end{Bmatrix}
  \f]

  Taken in the non-relativistic limit: the radial orbitals of
  \f$ j = l \pm 1/2 \f$ are treated as identical (overlap = 1).

  For a common non-relativistic shell (\f$ n_1 = n_2 \f$, \f$ l_1 = l_2 \f$)
  only L+S even terms exist (Pauli). When additionally \f$ j_1 \neq j_2 \f$
  the L+S odd components cancel in the antisymmetrisation and the even ones
  carry \f$ \eta = \sqrt{2} \f$; otherwise \f$ \eta = 1 \f$.
  In all cases \f$ \sum_{LS} A^2 = 1 \f$.

  @param n1,l1,twoj1  Quantum numbers of the first stored orbital.
  @param n2,l2,twoj2  Quantum numbers of the second stored orbital.
  @param L,S          Total orbital and spin angular momenta of the LS term.
  @param twoJ         Twice the total angular momentum 2J.
  @return Recoupling amplitude A(L,S); zero if forbidden.

  @note The sign convention follows the stored (sorted) orbital order; since
        nk_index sorting keeps the (n, l) order identical for all CSFs of one
        non-relativistic configuration, relative signs between such CSFs are
        consistent.
*/
double LS_amplitude(int n1, int l1, int twoj1, int n2, int l2, int twoj2, int L,
                    int S, int twoJ);

/*!
  @brief Expectation values of L^2 and S^2 for a two-electron CI state.
  @details
  Recouples each CSF to LS coupling (see @ref LS_amplitude) and accumulates,
  per non-relativistic configuration g,
  \f[
    B_g(L,S) = \sum_{I \in g} c_I A_I(L,S),
  \f]
  giving
  \f[
    \langle L^2 \rangle = \sum_{g,L,S} B_g(L,S)^2 \, L(L+1),
  \f]
  \f[
    \langle S^2 \rangle = \sum_{g,L,S} B_g(L,S)^2 \, S(S+1).
  \f]

  These are expectation values of the CI state, not eigenvalues: deviation
  from L(L+1), S(S+1) measures the LS-purity of the state.

  @param coefs  CI expansion coefficients (one per CSF).
  @param csfs   The CSF basis (matching @p coefs).
  @param twoJ   Twice the total angular momentum 2J.
  @return Pair {<L^2>, <S^2>}.

  @note Non-relativistic limit: radial overlaps between j = l +- 1/2 orbitals
        are set to 1, so for a normalised state the total LS weight is exactly
        1 and no renormalisation is required.
*/
std::pair<double, double>
expectation_L2S2(const LinAlg::View<const double> &coefs,
                 const std::vector<CSF2> &csfs, int twoJ);

/*!
  @brief Matrix elements of L^2 and S^2 between two two-electron CSFs of the
  same J (non-relativistic limit).
  @details
  Both are scalars in the two-electron space, and (with the radial orbitals
  of \f$ j = l \pm 1/2 \f$ treated as identical, see @ref LS_amplitude)
  connect only CSFs of the same non-relativistic configuration, in the same
  stored order. Then
  \f[
    \langle x'y'; J | L^2 | xy; J \rangle
    = \sum_{LS} A_{x'y'}(L,S) A_{xy}(L,S) \, L(L+1),
  \f]
  and likewise for \f$ S^2 \f$ with \f$ S(S+1) \f$.

  @param x,y    Orbitals (nk_index, sorted) of the ket CSF.
  @param xp,yp  Orbitals of the bra CSF.
  @param twoJ   Twice the total angular momentum 2J of both CSFs.
  @return Pair {<L^2>, <S^2>}; zero if the configurations differ.
*/
std::pair<double, double> L2S2_pair(DiracSpinor::Index x, DiracSpinor::Index y,
                                    DiracSpinor::Index xp,
                                    DiracSpinor::Index yp, int twoJ);

/*!
  @brief Expectation values of L^2 and S^2 for a three-electron CI state.
  @details
  Written as sums of pair operators plus one-body terms,
  \f[
    L^2 = \sum_{i<j} L_{ij}^2 - \sum_i l_i(l_i+1),
  \f]
  \f[
    S^2 = \sum_{i<j} S_{ij}^2 - \tfrac{9}{4},
  \f]
  the pair terms reduce, by the expansion over parent states (see @ref CSF3),
  to the two-electron matrix elements of @ref L2S2_pair, between parents with
  the same third orbital. As for two electrons, only CSFs of the same
  non-relativistic configuration are connected.

  @param coefs  CI expansion coefficients (one per CSF).
  @param csfs   The CSF basis (matching @p coefs).
  @param twoJ   Twice the total angular momentum 2J.
  @return Pair {<L^2>, <S^2>}.
*/
std::pair<double, double>
expectation_L2S2(const LinAlg::View<const double> &coefs,
                 const std::vector<CSF3> &csfs, int twoJ);

//==============================================================================
/*!
  @brief Identifies one CI level: its (J, parity), and which solution.
  @details
  The standard text form is `J{+,-}:index`, e.g., `2+:3` is the fourth solution
  (index counts from zero) of the J=2, even parity, CI. 
  The index may be omitted, in which case it is zero:
  `0+` is the lowest even-parity J=0 solution. Half-integer J is written as a
  fraction: `3/2-:1`. See @ref parse_level and @ref to_string.
*/
struct Level {
  //! Twice the total angular momentum, 2J
  int twoJ{0};
  //! Parity: +1 or -1
  int parity{1};
  //! Which solution, counting from zero, in order of energy
  std::size_t index{0};
};

/*!
  @brief Parses the text form of a CI level reference; see @ref Level.
  @details
  Accepts `2+:3` (standard) and `e2:3`; the index is optional. Half-integer J
  is given as a fraction over 2, or a decimal: `3/2-:1`, `1.5-:1`, `o3/2:1`.
  Surrounding whitespace is ignored.

  @param str  Text form, e.g., `2+:3`, `e2:3`, `0+`, `1/2-`.
  @return The level; empty if @p str is not a valid level reference.
*/
[[nodiscard]] std::optional<Level> parse_level(std::string_view str);

//! Text form of a CI level reference, e.g., "2+:3", "1/2-:0"; see @ref Level
[[nodiscard]] std::string to_string(const Level &level);

//! J as text, given 2J: "2" for even 2J, "3/2" for odd
[[nodiscard]] std::string J_string(int twoJ);

/*!
  @brief Parses a J value given as text: an integer, a fraction over 2, or a
  decimal.
  @details
  Accepts `2`, `1/2`, `3/2`, `0.5`, `1.5`, `2.0`; a decimal must be an integer
  or half-integer. Surrounding whitespace is ignored.
  @return 2J; empty if @p str is not a valid J.
*/
[[nodiscard]] std::optional<int> parse_twoJ(std::string_view str);

//==============================================================================
/*!
  @brief Configuration metadata for a single CI level.
  @details
  Stores identifying information derived after solving the CI eigenvalue
  problem: the dominant non-relativistic configuration label, the squared CI
  coefficient of that configuration, and approximate good quantum numbers
  (g_J factor, L, S) where they can be assigned.

  Fields are left at their default (empty/negative) values if not yet computed;
  call @ref PsiJPi::update_config_info() to populate them.
*/
struct ConfigInfo {
  //! Dominant configuration label (typically non-relativistic notation)
  std::string config{};
  //! Squared CI coefficient of the dominant configuration (or sum over non-rel degenerates)
  double ci2{0.0};
  double gJ{0.0};
  //! Approximate orbital angular momentum L (-1 if not assigned)
  double L{-1.0};
  //! Twice the approximate spin S (-1 if not assigned)
  double twoS{-1.0};
  //! Expectation value of L^2 for the CI state (-1 if not computed)
  double L2{-1.0};
  //! Expectation value of S^2 for the CI state (-1 if not computed)
  double S2{-1.0};
};

//==============================================================================
/*!
  @brief Container for CI solutions in a single J/Pi.
  @details
  Holds the complete set of configuration state functions and the results of
  the CI diagonalisation for a fixed total angular momentum J and parity.

  Construction builds the CSF basis (form_CSFs() for two electrons,
  form_CSF3s() for three) but does not solve the eigenvalue problem; call
  solve() separately after constructing the CI Hamiltonian matrix.
  Configuration labels are not set automatically -- call update_config_info()
  for each solution after solving.

  A J/Pi holds CSFs of a single electron number: two (@ref CSF2) or three
  (@ref CSF3); see num_electrons(). The CSF list is accessed by type (CSFs(),
  CSF3s()); num_CSFs() and csf() work for either.
*/
class PsiJPi {

  int m_twoj{-1};
  int m_pi{0};

  // Number of solutions stored:
  std::size_t m_num_solutions{0};
  // List of CSFs: two- or three-electron
  std::variant<std::vector<CSF2>, std::vector<CSF3>> m_CSFs{};
  // Energy, and CI expansion coeficients
  std::pair<LinAlg::Vector<double>, LinAlg::Matrix<double>> m_Solution{};
  std::vector<ConfigInfo> m_Info{};

public:
  /*!
    @brief Constructs the CSF basis for the given J and parity; does not solve.
    @details
    Calls form_CSFs() (two electrons) or form_CSF3s() (three) to build the
    list of CSFs.  The eigenvalue problem is not solved until solve() is called
    with the CI Hamiltonian.

    @param twoJ           Twice the total angular momentum 2J.
    @param pi             Total parity: +1 or -1.
    @param cisp_basis     Single-particle basis used to construct the CSFs.
    @param num_electrons  Number of valence electrons: 2 or 3.
  */
  PsiJPi(int twoJ, int pi, const std::vector<DiracSpinor> &cisp_basis,
         int num_electrons = 2);

  PsiJPi() {}

  /*!
    @brief Solves the CI eigenvalue problem for the given Hamiltonian matrix.
    @details
    Diagonalises @p Hci and stores the resulting eigenvalues and eigenvectors.
    Does not populate ConfigInfo; call update_config_info() separately.

    - If @p num_solutions > 0, finds only the lowest @p num_solutions eigenpairs.
    - If @p all_below is set, finds all eigenpairs with energy below that value
      (in cm^-1); @p num_solutions is then ignored.
    - If both are unset (or @p num_solutions <= 0), all eigenpairs are computed.

    @param Hci           CI Hamiltonian matrix in the CSF basis.
    @param num_solutions Number of lowest solutions to find [0 = all].
    @param all_below     If set, find all solutions below this energy (cm^-1).
  */
  void solve(const LinAlg::Matrix<double> &Hci, int num_solutions = 0,
             std::optional<double> all_below = {});

  /*!
    @brief Stores a single solution directly, without diagonalising.
    @details
    Replaces any existing solutions with the single one given.  For states that
    are not eigenstates of the CI Hamiltonian: e.g., the mixed states of
    @ref solve_mixed_state, for which @p energy is that of the reference state.

    @param energy Energy to be associated with the solution (atomic units).
    @param coefs  CI expansion coefficients; one per CSF.
  */
  void set_solution(double energy, const LinAlg::Vector<double> &coefs);

  //! Set configuration info for the ith solution (must be called manually after solve())
  void update_config_info(std::size_t i, const ConfigInfo &info);

  //! Number of valence electrons in each CSF (2 or 3)
  int num_electrons() const;

  //! Number of CSFs spanning this J/Pi
  std::size_t num_CSFs() const;

  //! The ith CSF, as the common base (for parity, config string, ...)
  const CI::CSF &csf(std::size_t i) const;

  //! Full list of two-electron CSFs spanning this J/Pi; must have num_electrons() == 2
  const std::vector<CSF2> &CSFs() const;

  //! Full list of three-electron CSFs spanning this J/Pi; must have num_electrons() == 3
  const std::vector<CSF3> &CSF3s() const;

  //! Returns reference to the ith two-electron CSF; must have num_electrons() == 2
  const CSF2 &CSF(std::size_t i) const;

  //! Energy of the ith CI solution (atomic units)
  double energy(std::size_t i) const;

  //! CI expansion coefficients for the ith solution (one per CSF)
  LinAlg::View<const double> coefs(std::size_t i) const;

  //! CI coefficient for the ith solution corresponding to the jth CSF
  double coef(std::size_t i, std::size_t j) const;

  //! Parity of the J/Pi (+/-1)
  int parity() const;

  //! Twice the total angular momentum 2J for this J/Pi
  int twoJ() const;

  //! Number of CI solutions currently stored
  std::size_t num_solutions() const;

  //! Configuration info for the ith solution (must have been set via update_config_info())
  const ConfigInfo &info(std::size_t i) const;

  /*!
    @brief Reads or writes CI solutions (energies, eigenvectors) to/from a multi-block binary file.
    @details
    A single file holds multiple blocks, one per (twoJ, parity) pair.
    Each block is self-describing: (twoJ, pi, num_csfs, num_solutions, E[num_csfs], M[num_csfs x num_csfs]).
    Energies and eigenvectors are always stored at full num_csfs size (zero-padded
    if only a partial solve was done), so block size is determined from num_csfs
    alone, enabling O(N) scan and in-place overwrite.

    The CSF basis (m_CSFs) is not touched; it must be constructed via the normal
    constructor before calling this function. The basis (and hence num_csfs) must
    be consistent with the file on write; a mismatch causes failure.

    ConfigInfo is not stored -- call update_config_info() after reading if needed.

    On read: scans blocks until matching (twoJ, pi) is found, verifies num_csfs,
    reads num_solutions, E, and M; returns false if block not found or num_csfs
    mismatches.

    On write: if the block already exists and num_csfs matches, overwrites it
    in-place using @ref IO::FRW::update. If it is a new block, appends it to
    the end of the file -- no rewrite of existing data.

    @note No settings are stored in the file: the filename identifies the
    calculation. The default filenames encode the settings that change the
    solutions -- see @ref CI::configuration_interaction.

    @param fname      Path to the binary file.
    @param rw         @ref IO::FRW::read to read; @ref IO::FRW::write to write.
    @param outstream  Output stream for progress/notes.

    @return True on success; false if the file does not exist (read), the block
            is not found (read), or num_csfs mismatches.
  */
  bool read_write(const std::string &fname, IO::FRW::RoW rw,
                  std::ostream &outstream = std::cout);
};

} // namespace CI
